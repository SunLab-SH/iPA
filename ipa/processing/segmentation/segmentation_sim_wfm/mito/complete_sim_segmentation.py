#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIM图像线粒体分割 - 完整独立脚本
包含Frangi滤波、分割标记、网络分析等所有功能，不依赖外部模块
"""

import numpy as np
from pathlib import Path
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree, distance
from skimage import morphology as morph
from skimage import measure
from itertools import combinations_with_replacement
import warnings
warnings.filterwarnings('ignore')

# 全局配置
xp = np  # 使用numpy替代cupy
device_type = 'cpu'

# ============================================================================
# 阈值函数实现
# ============================================================================

def triangle_threshold(image):
    """三角阈值算法"""
    hist, bin_edges = np.histogram(image, bins=256)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 找到直方图的最大值点
    max_idx = np.argmax(hist)
    
    # 计算从最大值点到最右端的直线距离
    max_distance = 0
    threshold_idx = max_idx
    
    for i in range(max_idx + 1, len(hist)):
        # 计算点到直线的距离
        x1, y1 = max_idx, hist[max_idx]
        x2, y2 = len(hist) - 1, hist[-1]
        x0, y0 = i, hist[i]
        
        distance = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) / np.sqrt((y2-y1)**2 + (x2-x1)**2)
        
        if distance > max_distance:
            max_distance = distance
            threshold_idx = i
    
    return bin_centers[threshold_idx]

def otsu_threshold(image):
    """Otsu阈值算法"""
    hist, bin_edges = np.histogram(image, bins=256)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    total = hist.sum()
    current_max = 0
    threshold = 0
    sum_total = (bin_centers * hist).sum()
    sum_back = 0
    weight_back = 0
    
    for i in range(len(hist)):
        weight_back += hist[i]
        if weight_back == 0:
            continue
            
        weight_fore = total - weight_back
        if weight_fore == 0:
            break
            
        sum_back += bin_centers[i] * hist[i]
        mean_back = sum_back / weight_back
        mean_fore = (sum_total - sum_back) / weight_fore
        
        between_variance = weight_back * weight_fore * (mean_back - mean_fore) ** 2
        
        if between_variance > current_max:
            current_max = between_variance
            threshold = bin_centers[i]
    
    return threshold, current_max

# ============================================================================
# Frangi滤波器实现
# ============================================================================

class SimpleFilter:
    """简化的Frangi滤波器"""
    
    def __init__(self, im_info, min_radius_um=0.1, max_radius_um=0.8, 
                 alpha_sq=0.3, beta_sq=0.3, remove_edges=True):
        self.im_info = im_info
        self.min_radius_um = max(min_radius_um, im_info.dim_res['X'])
        self.max_radius_um = max_radius_um
        self.min_radius_px = self.min_radius_um / im_info.dim_res['X']
        self.max_radius_px = self.max_radius_um / im_info.dim_res['X']
        self.alpha_sq = alpha_sq
        self.beta_sq = beta_sq
        self.remove_edges = remove_edges
        
        if not im_info.no_z:
            self.z_ratio = im_info.dim_res['Z'] / im_info.dim_res['X']
        
        self._set_default_sigmas()
    
    def _set_default_sigmas(self):
        """设置默认sigma值（优化版）"""
        min_sigma_step_size = 0.3
        num_sigma = 3  # 减少sigma数量以加速
        
        sigma_1 = self.min_radius_px / 2
        sigma_2 = self.max_radius_px / 3
        sigma_min = min(sigma_1, sigma_2)
        sigma_max = max(sigma_1, sigma_2)
        
        sigma_step_size = max(min_sigma_step_size, (sigma_max - sigma_min) / num_sigma)
        self.sigmas = list(np.arange(sigma_min, sigma_max + sigma_step_size/2, sigma_step_size))
        print(f"   Sigma值: {[f'{s:.2f}' for s in self.sigmas]} (优化到{len(self.sigmas)}个)")
    
    def _get_sigma_vec(self, sigma):
        """计算sigma向量"""
        if self.im_info.no_z:
            return (sigma, sigma)
        else:
            return (sigma / self.z_ratio, sigma, sigma)
    
    def _calculate_gamma(self, gauss_volume):
        """计算gamma阈值"""
        valid_pixels = gauss_volume[gauss_volume > 0]
        if len(valid_pixels) == 0:
            return 1.0
        
        gamma_tri = triangle_threshold(valid_pixels)
        gamma_otsu, _ = otsu_threshold(valid_pixels)
        return min(gamma_tri, gamma_otsu)
    
    def _compute_hessian(self, image):
        """计算Hessian矩阵"""
        gradients = np.gradient(image)
        axes = range(image.ndim)
        h_elems = np.array([np.gradient(gradients[ax0], axis=ax1).astype('float32')
                           for ax0, ax1 in combinations_with_replacement(axes, 2)])
        
        # 简化的掩码
        mask = np.ones_like(image, dtype=bool)
        
        if self.im_info.no_z:
            hxx, hxy, hyy = [elem[..., np.newaxis, np.newaxis] for elem in h_elems[:, mask]]
            hessian_matrices = np.concatenate([
                np.concatenate([hxx, hxy], axis=-1),
                np.concatenate([hxy, hyy], axis=-1)
            ], axis=-2)
        else:
            hxx, hxy, hxz, hyy, hyz, hzz = [elem[..., np.newaxis, np.newaxis] for elem in h_elems[:, mask]]
            hessian_matrices = np.concatenate([
                np.concatenate([hxx, hxy, hxz], axis=-1),
                np.concatenate([hxy, hyy, hyz], axis=-1),
                np.concatenate([hxz, hyz, hzz], axis=-1)
            ], axis=-2)
        
        return mask, hessian_matrices
    
    def _compute_eigenvalues(self, hessian_matrices, chunk_size=100000):
        """分块计算特征值"""
        total_voxels = len(hessian_matrices)
        eigenvalues_list = []
        
        for start_idx in range(0, total_voxels, chunk_size):
            end_idx = min(start_idx + chunk_size, total_voxels)
            chunk = hessian_matrices[start_idx:end_idx]
            chunk_eigenvalues = np.linalg.eigvalsh(chunk)
            eigenvalues_list.append(chunk_eigenvalues)
        
        eigenvalues_flat = np.concatenate(eigenvalues_list, axis=0)
        sort_order = np.argsort(np.abs(eigenvalues_flat), axis=1)
        eigenvalues_flat = np.take_along_axis(eigenvalues_flat, sort_order, axis=1)
        
        return eigenvalues_flat
    
    def _filter_hessian(self, eigenvalues, gamma_sq):
        """应用Frangi滤波器"""
        if self.im_info.no_z:
            # 2D情况
            rb_sq = (np.abs(eigenvalues[:, 0]) / (np.abs(eigenvalues[:, 1]) + 1e-10)) ** 2
            s_sq = (eigenvalues[:, 0] ** 2) + (eigenvalues[:, 1] ** 2)
            filtered_im = (np.exp(-(rb_sq / self.beta_sq))) * (1 - np.exp(-(s_sq / gamma_sq)))
            filtered_im[eigenvalues[:, 1] > 0] = 0
        else:
            # 3D情况
            ra_sq = (np.abs(eigenvalues[:, 1]) / (np.abs(eigenvalues[:, 2]) + 1e-10)) ** 2
            rb_sq = (np.abs(eigenvalues[:, 1]) / (np.sqrt(np.abs(eigenvalues[:, 1] * eigenvalues[:, 2])) + 1e-10)) ** 2
            s_sq = np.sqrt((eigenvalues[:, 0] ** 2) + (eigenvalues[:, 1] ** 2) + (eigenvalues[:, 2] ** 2)) ** 2
            
            filtered_im = ((1 - np.exp(-(ra_sq / self.alpha_sq))) * 
                          np.exp(-(rb_sq / self.beta_sq)) * 
                          (1 - np.exp(-(s_sq / gamma_sq))))
            
            filtered_im[eigenvalues[:, 2] > 0] = 0
            filtered_im[eigenvalues[:, 1] > 0] = 0
        
        filtered_im = np.nan_to_num(filtered_im, nan=0.0, posinf=0.0, neginf=0.0)
        return filtered_im
    
    def run_frame(self, image):
        """运行单帧Frangi滤波"""
        print(f"   处理图像尺寸: {image.shape}")
        total_pixels = np.prod(image.shape)
        print(f"   总像素数: {total_pixels:,}")
        
        vesselness = np.zeros_like(image, dtype='float64')
        temp = np.zeros_like(image, dtype='float64')
        
        start_time = time.time()
        for i, sigma in enumerate(self.sigmas):
            step_start = time.time()
            print(f"   处理sigma {i+1}/{len(self.sigmas)}: {sigma:.2f}...", end=' ')
            
            # 高斯滤波
            sigma_vec = self._get_sigma_vec(sigma)
            gauss_volume = ndi.gaussian_filter(image.astype('float64'), sigma=sigma_vec, mode='reflect')
            
            # 计算gamma
            gamma = self._calculate_gamma(gauss_volume)
            gamma_sq = 2 * gamma ** 2
            
            # 计算Hessian矩阵
            h_mask, hessian_matrices = self._compute_hessian(gauss_volume)
            
            if len(hessian_matrices) == 0:
                print("跳过")
                continue
            
            # 计算特征值（分块处理大图像）
            chunk_size = min(50000, len(hessian_matrices))  # 减小块大小
            eigenvalues = self._compute_eigenvalues(hessian_matrices.astype('float32'), chunk_size=chunk_size)
            
            # 应用Frangi滤波器
            filtered_result = self._filter_hessian(eigenvalues, gamma_sq)
            
            temp[h_mask] = filtered_result
            
            # 取最大值
            max_indices = temp > vesselness
            vesselness[max_indices] = temp[max_indices]
            
            step_time = time.time() - step_start
            remaining_steps = len(self.sigmas) - i - 1
            eta = step_time * remaining_steps
            print(f"完成 ({step_time:.1f}s, 剩余约{eta:.1f}s)")
        
        total_time = time.time() - start_time
        print(f"   Frangi滤波总耗时: {total_time:.1f}秒")
        return vesselness

# ============================================================================
# 分割标记实现
# ============================================================================

class SimpleLabel:
    """简化的分割标记"""
    
    def __init__(self, im_info, snr_cleaning=True, otsu_thresh_intensity=False):
        self.im_info = im_info
        self.snr_cleaning = snr_cleaning
        self.otsu_thresh_intensity = otsu_thresh_intensity
        
        if not im_info.no_z:
            self.min_z_radius_um = min(im_info.dim_res['Z'], 0.2)
    
    def _get_labels(self, frame):
        """生成二值标签"""
        ndim = 2 if self.im_info.no_z else 3
        footprint = ndi.generate_binary_structure(ndim, 1)
        
        # 计算阈值
        valid_pixels = frame[frame > 0]
        if len(valid_pixels) == 0:
            return np.zeros_like(frame, dtype=bool), np.zeros_like(frame, dtype='int32')
        
        triangle = triangle_threshold(valid_pixels)
        otsu, _ = otsu_threshold(valid_pixels)
        min_thresh = min(triangle, otsu)
        
        print(f"   阈值 - Triangle: {triangle:.1f}, Otsu: {otsu:.1f}, 使用: {min_thresh:.1f}")
        
        mask = frame > min_thresh
        
        # 形态学操作
        if not self.im_info.no_z:
            mask = ndi.binary_fill_holes(mask)
            if hasattr(self, 'min_z_radius_um') and self.im_info.dim_res['Z'] >= self.min_z_radius_um:
                mask = ndi.binary_opening(mask, structure=np.ones((2, 2, 2)))
        else:
            mask = ndi.binary_opening(mask, structure=np.ones((2, 2)))
        
        # 连通组件标记
        labels, num_features = ndi.label(mask, structure=footprint)
        
        # 移除小对象
        areas = np.bincount(labels.ravel())[1:]
        if len(areas) > 0:
            keep_labels = np.where(areas >= 4)[0] + 1
            mask = np.isin(labels, keep_labels)
            labels, _ = ndi.label(mask, structure=footprint)
        
        print(f"   检测到 {np.max(labels)} 个连通组件")
        return mask, labels
    
    def _get_object_snrs(self, original_frame, labels_frame):
        """计算对象的信噪比"""
        if not self.snr_cleaning:
            return labels_frame
        
        print("   执行SNR清理...")
        
        subtraction_mask = original_frame.copy()
        subtraction_mask[labels_frame > 0] = 0
        
        unique_labels = np.unique(labels_frame)
        keep_labels = []
        
        for label in unique_labels:
            if label == 0:
                continue
            
            coords = np.where(labels_frame == label)
            
            if len(coords[0]) == 0:
                continue
            
            # 扩展边界框
            extend_bbox_by = 1
            if self.im_info.no_z:
                rmin, rmax = np.min(coords[0]) - extend_bbox_by, np.max(coords[0]) + extend_bbox_by
                cmin, cmax = np.min(coords[1]) - extend_bbox_by, np.max(coords[1]) + extend_bbox_by
                
                rmin = max(0, rmin)
                rmax = min(labels_frame.shape[0], rmax)
                cmin = max(0, cmin)  
                cmax = min(labels_frame.shape[1], cmax)
                
                local_intensity = subtraction_mask[rmin:rmax, cmin:cmax]
            else:
                zmin, zmax = np.min(coords[0]) - extend_bbox_by, np.max(coords[0]) + extend_bbox_by
                rmin, rmax = np.min(coords[1]) - extend_bbox_by, np.max(coords[1]) + extend_bbox_by
                cmin, cmax = np.min(coords[2]) - extend_bbox_by, np.max(coords[2]) + extend_bbox_by
                
                zmin = max(0, zmin)
                zmax = min(labels_frame.shape[0], zmax)
                rmin = max(0, rmin)
                rmax = min(labels_frame.shape[1], rmax)
                cmin = max(0, cmin)
                cmax = min(labels_frame.shape[2], cmax)
                
                local_intensity = subtraction_mask[zmin:zmax, rmin:rmax, cmin:cmax]
            
            local_valid = local_intensity[local_intensity > 0]
            if len(local_valid) == 0:
                continue
                
            local_mean = np.mean(local_valid)
            local_std = np.std(local_valid)
            label_mean = np.mean(original_frame[coords])
            
            if local_std > 0:
                intensity_cutoff = label_mean / (local_mean + local_std)
                if intensity_cutoff > 1:
                    keep_labels.append(label)
        
        # 保留符合SNR条件的标签
        keep_labels = np.array(keep_labels)
        labels_frame = np.where(np.isin(labels_frame, keep_labels), labels_frame, 0)
        
        print(f"   SNR清理后保留 {len(keep_labels)} 个对象")
        return labels_frame
    
    def run_frame(self, original_frame, frangi_frame):
        """运行单帧分割"""
        print(f"   分割图像尺寸: {frangi_frame.shape}")
        
        # 生成标签
        _, labels = self._get_labels(frangi_frame)
        
        # SNR清理
        if self.snr_cleaning and np.max(labels) > 0:
            labels = self._get_object_snrs(original_frame, labels)
        
        return labels

# ============================================================================
# 简化的ImInfo类
# ============================================================================

class SimpleImInfo:
    """简化的图像信息类"""
    
    def __init__(self, im_path, processed_image, dim_res):
        self.im_path = im_path
        self.processed_image = processed_image
        self.dim_res = dim_res
        self.shape = processed_image.shape
        
        # 根据维度设置属性
        if len(self.shape) == 2:
            self.no_z = True
            self.no_t = True
            self.axes = ['Y', 'X']
        elif len(self.shape) == 3:
            self.no_z = False
            self.no_t = True
            self.axes = ['Z', 'Y', 'X']
        else:
            self.no_z = False
            self.no_t = False
        
        # 设置输出路径
        base_dir = Path(im_path).parent / "temp_processing"
        base_dir.mkdir(exist_ok=True)
        
        self.pipeline_paths = {
            'im_preprocessed': str(base_dir / "frangi_filtered.tif"),
            'im_instance_label': str(base_dir / "instance_labels.tif")
        }
        
        # 存储处理结果
        self.frangi_result = None
        self.label_result = None

# ============================================================================
# 主处理函数
# ============================================================================

def prepare_sim_image(im_path, channel=0, downsample_factor=2, crop_center=True):
    """准备SIM图像数据（优化版）"""
    print(f"🔬 加载SIM图像: {Path(im_path).name}")
    
    with tifffile.TiffFile(im_path) as tif:
        image_array = tif.asarray()
    
    print(f"   原始尺寸: {image_array.shape}")
    print(f"   数据类型: {image_array.dtype}")
    print(f"   图像大小: {image_array.nbytes / 1024**3:.2f} GB")
    
    # 处理多通道图像
    if len(image_array.shape) == 4 and image_array.shape[0] <= 10:
        selected_image = image_array[channel]
        print(f"   选择通道 {channel}: {selected_image.shape}")
    else:
        selected_image = image_array
    
    # 优化：裁剪中心区域（减少边缘噪声）
    if crop_center and len(selected_image.shape) >= 2:
        h, w = selected_image.shape[-2:]
        if h > 1024 or w > 1024:
            # 裁剪到中心1024x1024区域
            center_h, center_w = h // 2, w // 2
            crop_h, crop_w = min(1024, h), min(1024, w)
            start_h = max(0, center_h - crop_h // 2)
            end_h = start_h + crop_h
            start_w = max(0, center_w - crop_w // 2)
            end_w = start_w + crop_w
            
            if len(selected_image.shape) == 3:
                selected_image = selected_image[:, start_h:end_h, start_w:end_w]
            else:
                selected_image = selected_image[start_h:end_h, start_w:end_w]
            
            print(f"   裁剪后尺寸: {selected_image.shape}")
    
    # 优化：下采样以加速处理
    if downsample_factor > 1 and len(selected_image.shape) >= 2:
        print(f"   应用 {downsample_factor}x 下采样...")
        
        if len(selected_image.shape) == 3:
            # 3D下采样
            new_shape = (selected_image.shape[0], 
                        selected_image.shape[1] // downsample_factor,
                        selected_image.shape[2] // downsample_factor)
            selected_image = transform.resize(selected_image, new_shape, 
                                            preserve_range=True, 
                                            anti_aliasing=True).astype(selected_image.dtype)
        else:
            # 2D下采样  
            new_shape = (selected_image.shape[0] // downsample_factor,
                        selected_image.shape[1] // downsample_factor)
            selected_image = transform.resize(selected_image, new_shape,
                                            preserve_range=True,
                                            anti_aliasing=True).astype(selected_image.dtype)
        
        print(f"   下采样后尺寸: {selected_image.shape}")
        print(f"   处理后大小: {selected_image.nbytes / 1024**2:.1f} MB")
    
    return selected_image

def optimize_for_large_image(image):
    """针对大图像的优化处理"""
    print(f"🚀 图像优化处理...")
    
    # 计算图像统计信息来决定是否需要进一步优化
    total_pixels = np.prod(image.shape)
    print(f"   总像素数: {total_pixels:,}")
    
    if total_pixels > 50_000_000:  # 5000万像素以上
        print("   ⚠️  检测到大图像，建议进一步优化")
        
        # 选择最有信号的Z层（如果是3D）
        if len(image.shape) == 3 and image.shape[0] > 10:
            print("   选择信号最强的Z层子集...")
            z_means = np.mean(image, axis=(1,2))
            z_std = np.std(image, axis=(1,2))
            z_scores = z_means + z_std  # 结合均值和标准差
            
            # 选择top 30%的Z层
            num_slices = max(10, image.shape[0] // 3)
            top_indices = np.argsort(z_scores)[-num_slices:]
            top_indices = np.sort(top_indices)
            
            image = image[top_indices]
            print(f"   保留 {len(top_indices)} 个Z层: {image.shape}")
    
    return image

def run_sim_mitochondria_segmentation(im_path, channel=0, dim_res=None, output_dir=None):
    """运行SIM图像线粒体分割"""
    
    print("🔬 SIM图像线粒体分割开始")
    print("=" * 50)
    
    # 设置默认参数
    if dim_res is None:
        dim_res = {'Z': 0.125, 'Y': 0.063, 'X': 0.063}
    
    if output_dir is None:
        output_dir = Path(im_path).parent / "mitochondria_results"
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # 1. 准备图像
    print("\n步骤1: 准备图像数据")
    # 更激进的下采样以加速处理
    processed_image = prepare_sim_image(im_path, channel, downsample_factor=4, crop_center=True)
    
    # 进一步优化大图像
    processed_image = optimize_for_large_image(processed_image)
    
    # 保存处理后的图像
    processed_path = Path(output_dir) / "processed_image.tif"
    tifffile.imwrite(processed_path, processed_image)
    print(f"   已保存处理后图像: {processed_path}")
    
    # 2. 创建ImInfo对象
    print("\n步骤2: 创建图像信息对象")
    im_info = SimpleImInfo(str(processed_path), processed_image, dim_res)
    print("   ✅ ImInfo对象创建成功")
    
    # 3. Frangi滤波
    print("\n步骤3: Frangi滤波预处理")
    try:
        filter_obj = SimpleFilter(
            im_info=im_info,
            min_radius_um=0.1,
            max_radius_um=0.8,
            alpha_sq=0.3,
            beta_sq=0.3,
            remove_edges=True
        )
        
        frangi_result = filter_obj.run_frame(processed_image)
        print("   ✅ Frangi滤波完成")
        
        # 保存滤波结果
        frangi_result_path = Path(output_dir) / "frangi_filtered.tif"
        tifffile.imwrite(frangi_result_path, frangi_result.astype('float32'))
        print(f"   已保存Frangi滤波结果: {frangi_result_path}")
        
        # 保存结果到im_info
        im_info.frangi_result = frangi_result
        
    except Exception as e:
        print(f"   ❌ Frangi滤波失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 4. 分割标记
    print("\n步骤4: 分割标记")
    try:
        label_obj = SimpleLabel(
            im_info=im_info,
            snr_cleaning=True,
            otsu_thresh_intensity=False
        )
        
        label_result = label_obj.run_frame(processed_image, frangi_result)
        print("   ✅ 分割标记完成")
        
        # 保存分割结果
        label_result_path = Path(output_dir) / "instance_labels.tif"
        tifffile.imwrite(label_result_path, label_result.astype('uint16'))
        print(f"   已保存分割结果: {label_result_path}")
        
        # 统计分割结果
        unique_labels = np.unique(label_result)
        num_objects = len(unique_labels) - 1  # 减去背景
        print(f"   最终检测到 {num_objects} 个线粒体对象")
        
        # 保存结果到im_info
        im_info.label_result = label_result
        
    except Exception as e:
        print(f"   ❌ 分割标记失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 5. 结果统计和可视化
    print("\n步骤5: 结果分析")
    try:
        analyze_results(processed_image, frangi_result, label_result, output_dir)
    except Exception as e:
        print(f"   ⚠️  结果分析失败: {e}")
    
    print(f"\n🎉 处理完成！结果保存在: {output_dir}")
    return True

def analyze_results(original_image, frangi_result, label_result, output_dir):
    """分析处理结果"""
    
    # 基本统计
    unique_labels = np.unique(label_result)
    num_objects = len(unique_labels) - 1
    
    print(f"   检测对象数量: {num_objects}")
    print(f"   原图像强度范围: {original_image.min():.0f} - {original_image.max():.0f}")
    print(f"   Frangi响应范围: {frangi_result.min():.3f} - {frangi_result.max():.3f}")
    
    if num_objects > 0:
        # 对象尺寸统计
        sizes = []
        for label in unique_labels[1:]:  # 跳过背景
            size = np.sum(label_result == label)
            sizes.append(size)
        
        sizes = np.array(sizes)
        print(f"   对象尺寸统计:")
        print(f"     平均像素数: {np.mean(sizes):.1f}")
        print(f"     最小像素数: {np.min(sizes)}")
        print(f"     最大像素数: {np.max(sizes)}")
    
    # 创建可视化图像 (如果是3D，取中间切片)
    if len(original_image.shape) == 3:
        mid_slice = original_image.shape[0] // 2
        vis_original = original_image[mid_slice]
        vis_frangi = frangi_result[mid_slice]
        vis_labels = label_result[mid_slice]
    else:
        vis_original = original_image
        vis_frangi = frangi_result  
        vis_labels = label_result
    
    # 保存可视化结果
    try:
        # 归一化到0-255
        vis_original_norm = ((vis_original - vis_original.min()) / (vis_original.max() - vis_original.min()) * 255).astype('uint8')
        vis_frangi_norm = ((vis_frangi - vis_frangi.min()) / (vis_frangi.max() - vis_frangi.min()) * 255).astype('uint8')
        
        tifffile.imwrite(Path(output_dir) / "visualization_original.tif", vis_original_norm)
        tifffile.imwrite(Path(output_dir) / "visualization_frangi.tif", vis_frangi_norm)
        tifffile.imwrite(Path(output_dir) / "visualization_labels.tif", vis_labels.astype('uint16'))
        
        print(f"   已保存可视化结果")
        
    except Exception as e:
        print(f"   保存可视化失败: {e}")

def main():
    """主函数"""
    print("🚀 SIM线粒体分割 - 完整独立版本")
    print("✅ 使用集成模式，无需外部文件")
    
    # SIM图像路径
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ 图像文件不存在: {sim_path}")
        return
    
    # 运行分割
    success = run_sim_mitochondria_segmentation(
        im_path=sim_path,
        channel=0,  # 选择第一个通道
        dim_res={'Z': 0.125, 'Y': 0.063, 'X': 0.063}
    )
    
    if success:
        print("\n✅ SIM图像线粒体分割成功完成！")
        print("\n📁 输出文件说明:")
        print("   - processed_image.tif: 预处理后的图像")
        print("   - frangi_filtered.tif: Frangi滤波结果") 
        print("   - instance_labels.tif: 分割标签结果")
        print("   - visualization_*.tif: 可视化图像")
    else:
        print("\n❌ 分割过程中出现错误")

if __name__ == "__main__":
    main()
    input("\n按Enter键退出...")