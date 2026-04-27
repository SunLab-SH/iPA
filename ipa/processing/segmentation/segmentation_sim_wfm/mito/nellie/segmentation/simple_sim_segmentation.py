#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIM图像线粒体分割 - 直接调用脚本
直接使用现有的filtering.py, labelling.py等模块处理SIM图像
"""

import sys
import os
import numpy as np
from pathlib import Path
import tifffile

def setup_environment():
    """设置环境和路径"""
    # 添加当前目录到Python路径
    current_dir = Path(__file__).parent
    sys.path.insert(0, str(current_dir))
    
    # 检查必需文件
    required_files = ['filtering.py', 'labelling.py', 'networking.py', 'mocap_marking.py']
    missing_files = []
    
    for file in required_files:
        if not (current_dir / file).exists():
            missing_files.append(file)
    
    if missing_files:
        print(f"❌ 缺少必需文件: {missing_files}")
        return False
    
    return True

def prepare_sim_image(im_path, channel=0):
    """准备SIM图像数据"""
    print(f"🔬 加载SIM图像: {Path(im_path).name}")
    
    # 加载图像
    with tifffile.TiffFile(im_path) as tif:
        image_array = tif.asarray()
    
    print(f"   原始尺寸: {image_array.shape}")
    print(f"   数据类型: {image_array.dtype}")
    
    # 如果是多通道图像，选择指定通道
    if len(image_array.shape) == 4:  # (C, Z, Y, X)
        if image_array.shape[0] <= 10:  # 第一维是通道
            selected_image = image_array[channel]
            print(f"   选择通道 {channel}: {selected_image.shape}")
        else:
            selected_image = image_array
    else:
        selected_image = image_array
    
    return selected_image

def create_mock_im_info(im_path, processed_image, dim_res):
    """创建模拟的ImInfo对象"""
    class MockImInfo:
        def __init__(self):
            self.im_path = im_path
            self.shape = processed_image.shape
            self.dim_res = dim_res
            
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
                
            # 设置管道路径
            base_dir = Path(im_path).parent / "temp_processing"
            base_dir.mkdir(exist_ok=True)
            
            self.pipeline_paths = {
                'im_preprocessed': str(base_dir / "frangi_filtered.tif"),
                'im_instance_label': str(base_dir / "instance_labels.tif"),
                'im_skel': str(base_dir / "skeleton.tif"),
                'im_pixel_class': str(base_dir / "pixel_class.tif"),
                'im_skel_relabelled': str(base_dir / "skeleton_relabelled.tif"),
                'im_marker': str(base_dir / "markers.tif"),
                'im_distance': str(base_dir / "distance.tif"),
                'im_border': str(base_dir / "border.tif")
            }
        
        def get_memmap(self, path):
            if path == self.im_path:
                return processed_image
            elif Path(path).exists():
                return tifffile.imread(path)
            return processed_image
        
        def allocate_memory(self, path, dtype='float32', description='', return_memmap=True):
            # 简单保存为TIFF文件
            shape = processed_image.shape
            zeros_array = np.zeros(shape, dtype=dtype)
            tifffile.imwrite(path, zeros_array)
            return zeros_array
    
    return MockImInfo()

def run_sim_mitochondria_segmentation(im_path, 
                                     channel=0,
                                     dim_res=None,
                                     output_dir=None):
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
    print("步骤1: 准备图像数据")
    processed_image = prepare_sim_image(im_path, channel)
    
    # 保存处理后的图像
    processed_path = Path(output_dir) / "processed_image.tif"
    tifffile.imwrite(processed_path, processed_image)
    print(f"   已保存处理后图像: {processed_path}")
    
    # 2. 创建ImInfo对象
    print("\n步骤2: 创建图像信息对象")
    im_info = create_mock_im_info(str(processed_path), processed_image, dim_res)
    print("   ✅ ImInfo对象创建成功")
    
    # 3. Frangi滤波
    print("\n步骤3: Frangi滤波预处理")
    try:
        from filtering import Filter
        
        filter_obj = Filter(
            im_info=im_info,
            min_radius_um=0.1,
            max_radius_um=0.8,
            alpha_sq=0.3,
            beta_sq=0.3,
            remove_edges=True
        )
        
        filter_obj.run(mask=True)
        print("   ✅ Frangi滤波完成")
        
        # 保存滤波结果
        frangi_result_path = Path(output_dir) / "frangi_filtered.tif"
        frangi_data = filter_obj.frangi_memmap[:]
        tifffile.imwrite(frangi_result_path, frangi_data)
        print(f"   已保存Frangi滤波结果: {frangi_result_path}")
        
    except Exception as e:
        print(f"   ❌ Frangi滤波失败: {e}")
        return False
    
    # 4. 分割标记
    print("\n步骤4: 分割标记")
    try:
        from labelling import Label
        
        label_obj = Label(
            im_info=im_info,
            snr_cleaning=True,
            otsu_thresh_intensity=False
        )
        
        label_obj.run()
        print("   ✅ 分割标记完成")
        
        # 保存分割结果
        label_result_path = Path(output_dir) / "instance_labels.tif"
        label_data = label_obj.instance_label_memmap[:]
        tifffile.imwrite(label_result_path, label_data.astype('uint16'))
        print(f"   已保存分割结果: {label_result_path}")
        
        # 统计分割结果
        unique_labels = np.unique(label_data)
        num_objects = len(unique_labels) - 1  # 减去背景
        print(f"   检测到 {num_objects} 个对象")
        
    except Exception as e:
        print(f"   ❌ 分割标记失败: {e}")
        return False
    
    # 5. 网络分析（可选）
    print("\n步骤5: 网络分析")
    try:
        from networking import Network
        
        network_obj = Network(
            im_info=im_info,
            min_radius_um=0.1,
            max_radius_um=0.8,
            clean_skel=True
        )
        
        network_obj.run()
        print("   ✅ 网络分析完成")
        
        # 保存骨架结果
        skel_result_path = Path(output_dir) / "skeleton.tif"
        skel_data = network_obj.skel_memmap[:]
        tifffile.imwrite(skel_result_path, skel_data.astype('uint16'))
        print(f"   已保存骨架结果: {skel_result_path}")
        
    except Exception as e:
        print(f"   ⚠️  网络分析失败: {e}")
        print("   继续处理...")
    
    # 6. 标记点检测（可选）
    print("\n步骤6: 标记点检测")
    try:
        from mocap_marking import Markers
        
        markers_obj = Markers(
            im_info=im_info,
            min_radius_um=0.1,
            max_radius_um=0.8,
            use_im='distance',
            num_sigma=5
        )
        
        markers_obj.run()
        print("   ✅ 标记点检测完成")
        
        # 保存标记点结果
        markers_result_path = Path(output_dir) / "markers.tif"
        markers_data = markers_obj.im_marker_memmap[:]
        tifffile.imwrite(markers_result_path, markers_data.astype('uint8'))
        print(f"   已保存标记点结果: {markers_result_path}")
        
    except Exception as e:
        print(f"   ⚠️  标记点检测失败: {e}")
        print("   继续处理...")
    
    print(f"\n🎉 处理完成！结果保存在: {output_dir}")
    return True

def main():
    """主函数"""
    # 检查环境
    if not setup_environment():
        return
    
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
    else:
        print("\n❌ 分割过程中出现错误")

if __name__ == "__main__":
    main()
    input("\n按Enter键退出...")