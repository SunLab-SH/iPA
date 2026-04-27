#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
完整的线粒体分割管道
用于处理静态3D SIM显微镜图像的线粒体分割

功能包括:
1. Frangi血管滤波器预处理
2. 阈值分割和标记
3. 网络骨架化分析
4. 运动捕获标记点检测

作者: GitHub Copilot

使用说明:
1. 确保所有依赖文件 (filtering.py, labelling.py, networking.py, mocap_marking.py) 在同一目录
2. 调用 run_mitochondria_segmentation(image_path) 来执行完整分割
3. 结果将保存在指定的输出目录中

示例:
    pipeline = run_mitochondria_segmentation(
        im_path="path/to/your/image.tif",
        dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1}
    )
"""

import os
import sys
import numpy as np
from pathlib import Path
import logging
import json
import time
from typing import Optional, Dict, Tuple, Union

# 设置日志
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MitochondriaSegmentationPipeline:
    """
    线粒体分割完整管道类
    
    用于处理静态3D SIM显微镜图像，执行完整的线粒体分割流程：
    Frangi滤波 → 分割标记 → 网络骨架化 → 关键点检测
    
    Attributes
    ----------
    im_path : str
        输入图像文件路径
    output_dir : str
        输出目录路径
    im_info : object
        图像信息对象
    config : dict
        处理配置参数
    results : dict
        处理结果存储
    """
    
    def __init__(self, 
                 im_path: str,
                 output_dir: Optional[str] = None,
                 dim_res: Optional[Dict[str, float]] = None,
                 dimension_order: str = 'ZYX',
                 channel: Optional[int] = None,
                 **kwargs):
        """
        初始化线粒体分割管道
        
        Parameters
        ----------
        im_path : str
            输入图像文件路径
        output_dir : str, optional
            输出目录路径，默认为输入文件同目录
        dim_res : dict, optional
            像素分辨率，格式如 {'Z': 0.2, 'Y': 0.1, 'X': 0.1}
        dimension_order : str, optional
            维度顺序，默认 'ZYX'，多通道图像可能是 'CZYX' 或 'ZCYX'
        channel : int, optional
            要处理的通道索引，None表示处理所有通道或单通道图像
        **kwargs : dict
            其他配置参数
        """
        self.im_path = str(Path(im_path).resolve())
        
        # 设置输出目录
        if output_dir is None:
            self.output_dir = str(Path(im_path).parent / "mitochondria_segmentation_results")
        else:
            self.output_dir = str(Path(output_dir).resolve())
        
        # 创建输出目录
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        # 设置默认分辨率
        if dim_res is None:
            dim_res = {'Z': 0.2, 'Y': 0.1, 'X': 0.1}  # 默认SIM显微镜分辨率
        
        self.dim_res = dim_res
        self.dimension_order = dimension_order
        self.channel = channel
        
        # 动态导入模块
        self._import_modules()
        
        # 预加载和分析图像
        self._analyze_image_structure()
        
        # 初始化图像信息
        try:
            if self.ImInfo is not None and hasattr(self, 'processed_image') and self.processed_image is not None:
                # 尝试使用Nellie的ImInfo（如果可用）
                try:
                    # 保存处理后的图像到临时文件
                    temp_path = Path(self.output_dir) / "temp_processed_image.tif"
                    import tifffile
                    tifffile.imwrite(temp_path, self.processed_image)
                    
                    self.im_info = self.ImInfo(str(temp_path), 
                                             dim_res=dim_res, 
                                             dimension_order=dimension_order)
                    logger.info(f"使用Nellie ImInfo，图像尺寸: {self.im_info.shape}")
                except Exception as e:
                    logger.warning(f"Nellie ImInfo初始化失败: {e}，使用简化模式")
                    self.im_info = self._create_simple_im_info()
            else:
                logger.info("使用简化ImInfo模式")
                self.im_info = self._create_simple_im_info()
                
        except Exception as e:
            logger.error(f"图像信息初始化失败: {e}")
            self.im_info = None
        
        # 设置处理配置
        self.config = self._setup_default_config(**kwargs)
        
        # 初始化结果字典
        self.results = {}
        
        # 设置日志
        self._setup_logging()
        
        logger.info(f"初始化线粒体分割管道")
        logger.info(f"输入图像: {self.im_path}")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info(f"像素分辨率: {dim_res}")
    
    def _import_modules(self):
        """动态导入所需的模块"""
        self.ImInfo = None
        self.Filter = None
        self.Label = None
        self.Network = None
        self.Markers = None
        
        # 使用相对导入
        try:
            # 尝试导入ImInfo
            try:
                from ..im_info.verifier import ImInfo
                self.ImInfo = ImInfo
                logger.info("成功导入 ImInfo")
            except ImportError:
                logger.warning("无法导入 ImInfo，将使用简化模式")
            
            # 导入处理模块
            try:
                from .filtering import Filter
                self.Filter = Filter
                logger.info("成功导入 Filter")
            except ImportError as e:
                logger.error(f"无法导入 Filter: {e}")
            
            try:
                from .labelling import Label  
                self.Label = Label
                logger.info("成功导入 Label")
            except ImportError as e:
                logger.error(f"无法导入 Label: {e}")
            
            try:
                from .networking import Network
                self.Network = Network
                logger.info("成功导入 Network")
            except ImportError as e:
                logger.error(f"无法导入 Network: {e}")
            
            try:
                from .mocap_marking import Markers
                self.Markers = Markers  
                logger.info("成功导入 Markers")
            except ImportError as e:
                logger.error(f"无法导入 Markers: {e}")
                
        except Exception as e:
            logger.error(f"导入模块时发生错误: {e}")
    
    def _setup_default_config(self, **kwargs) -> Dict:
        """设置默认配置参数"""
        default_config = {
            # Frangi滤波器参数
            'frangi': {
                'min_radius_um': 0.15,      # 最小线粒体半径 (微米)
                'max_radius_um': 1.0,       # 最大线粒体半径 (微米)  
                'alpha_sq': 0.5,            # Frangi α参数
                'beta_sq': 0.5,             # Frangi β参数
                'remove_edges': True,       # 移除边缘
            },
            # 分割参数
            'labelling': {
                'threshold': None,          # 强度阈值 (None为自动)
                'snr_cleaning': True,       # SNR清理
                'otsu_thresh_intensity': False,  # 使用Otsu阈值
            },
            # 网络分析参数
            'networking': {
                'min_radius_um': 0.15,
                'max_radius_um': 1.0,
                'clean_skel': True,         # 清理骨架
            },
            # 标记点检测参数
            'markers': {
                'min_radius_um': 0.15,
                'max_radius_um': 1.0,
                'use_im': 'distance',       # 使用距离图像检测峰值
                'num_sigma': 5,             # 多尺度σ数量
            },
            # 通用参数
            'num_t': None,                  # 时间点数量 (静态图像为None)
            'save_intermediate': True,      # 保存中间结果
        }
        
        # 更新用户提供的参数
        for section, params in kwargs.items():
            if section in default_config:
                default_config[section].update(params)
            else:
                default_config[section] = params
        
        return default_config
    
    def _setup_logging(self):
        """设置日志输出"""
        log_file = Path(self.output_dir) / "segmentation.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, encoding='utf-8'),
                logging.StreamHandler()
            ]
        )
    
    def _analyze_image_structure(self):
        """分析图像结构，检测通道、尺寸等信息"""
        try:
            import tifffile
            import numpy as np
            
            # 读取图像元数据
            with tifffile.TiffFile(self.im_path) as tif:
                image_array = tif.asarray()
                
            self.original_shape = image_array.shape
            self.image_array = image_array
            
            logger.info(f"原始图像尺寸: {self.original_shape}")
            
            # 分析维度结构
            self._detect_image_dimensions()
            
            # 如果是多通道图像，选择或处理通道
            if self.has_channels:
                self._process_channels()
            else:
                self.processed_image = self.image_array
                
            if hasattr(self, 'processed_image') and self.processed_image is not None:
                logger.info(f"处理后图像尺寸: {self.processed_image.shape}")
            
        except ImportError:
            logger.warning("未安装tifffile，尝试使用PIL")
            self._load_with_pil()
        except Exception as e:
            logger.error(f"图像分析失败: {e}")
            self.image_array = None
            self.processed_image = None
    
    def _detect_image_dimensions(self):
        """检测图像的维度结构"""
        shape = self.original_shape
        ndim = len(shape)
        
        # 根据维度数量和尺寸推断结构
        if ndim == 2:
            # 2D图像: (H, W)
            self.has_channels = False
            self.has_z = False
            self.actual_dimension_order = 'YX'
            
        elif ndim == 3:
            # 3D图像: 可能是 (Z, H, W) 或 (C, H, W) 或 (H, W, C)
            if shape[0] <= 10:  # 假设通道数不超过10
                self.has_channels = True
                self.has_z = False
                self.actual_dimension_order = 'CYX'
            elif shape[-1] <= 10:  # 最后一个维度是通道
                self.has_channels = True
                self.has_z = False
                self.actual_dimension_order = 'YXC'
            else:  # Z堆叠
                self.has_channels = False
                self.has_z = True
                self.actual_dimension_order = 'ZYX'
                
        elif ndim == 4:
            # 4D图像: (C, Z, H, W) 或 (Z, C, H, W) 或其他
            if shape[0] <= 10:  # 第一个维度是通道
                self.has_channels = True
                self.has_z = True
                self.actual_dimension_order = 'CZYX'
            elif shape[1] <= 10:  # 第二个维度是通道
                self.has_channels = True
                self.has_z = True
                self.actual_dimension_order = 'ZCYX'
            else:
                # 假设是时间序列
                self.has_channels = False
                self.has_z = True
                self.actual_dimension_order = 'TZYX'
        else:
            logger.warning(f"不支持的图像维度: {ndim}")
            self.has_channels = False
            self.has_z = False
            self.actual_dimension_order = 'YX'
        
        logger.info(f"检测到的维度结构: {self.actual_dimension_order}")
        logger.info(f"是否有通道: {self.has_channels}, 是否有Z轴: {self.has_z}")
    
    def _process_channels(self):
        """处理多通道图像"""
        if not self.has_channels or self.image_array is None:
            self.processed_image = self.image_array
            return
        
        # 获取通道信息
        if 'C' in self.actual_dimension_order:
            channel_axis = self.actual_dimension_order.index('C')
            num_channels = self.original_shape[channel_axis]
            
            logger.info(f"检测到 {num_channels} 个通道")
            
            # 选择处理的通道
            if self.channel is None:
                if num_channels > 1:
                    logger.info("多通道图像，自动选择第一个通道进行处理")
                    self.channel = 0
                else:
                    self.channel = 0
            
            if self.channel >= num_channels:
                logger.warning(f"指定通道 {self.channel} 超出范围，使用通道 0")
                self.channel = 0
            
            # 提取指定通道
            if self.image_array is not None:
                if self.actual_dimension_order == 'CYX':
                    self.processed_image = self.image_array[self.channel]
                elif self.actual_dimension_order == 'YXC':
                    self.processed_image = self.image_array[..., self.channel]
                elif self.actual_dimension_order == 'CZYX':
                    self.processed_image = self.image_array[self.channel]
                elif self.actual_dimension_order == 'ZCYX':
                    self.processed_image = self.image_array[:, self.channel]
                else:
                    logger.error(f"不支持的通道提取: {self.actual_dimension_order}")
                    self.processed_image = self.image_array
                    
                if self.processed_image is not None:
                    logger.info(f"选择通道 {self.channel}，处理后尺寸: {self.processed_image.shape}")
            else:
                self.processed_image = None
        else:
            self.processed_image = self.image_array
    
    def _load_with_pil(self):
        """使用PIL加载图像作为备用方案"""
        try:
            from PIL import Image
            import numpy as np
            
            img = Image.open(self.im_path)
            self.image_array = np.array(img)
            self.original_shape = self.image_array.shape
            self.has_channels = len(self.original_shape) == 3 and self.original_shape[-1] <= 10
            self.has_z = False
            
            if self.has_channels:
                self.actual_dimension_order = 'YXC'
                self._process_channels()
            else:
                self.processed_image = self.image_array
                self.actual_dimension_order = 'YX' if len(self.original_shape) == 2 else 'ZYX'
                
        except Exception as e:
            logger.error(f"PIL加载图像失败: {e}")
            self.image_array = None
            self.processed_image = None
    
    def _create_simple_im_info(self):
        """创建简化的图像信息对象，当无法使用Nellie ImInfo时"""
        class SimpleImInfo:
            def __init__(self, im_path, processed_image, dim_res, dimension_order):
                self.im_path = im_path
                self.processed_image = processed_image
                self.dim_res = dim_res
                self.dimension_order = dimension_order
                
                if processed_image is not None:
                    self.shape = processed_image.shape
                    self.ndim = len(self.shape)
                    
                    # 根据形状判断维度
                    if self.ndim == 2:
                        self.no_z = True
                        self.no_t = True
                        self.axes = ['Y', 'X']
                    elif self.ndim == 3:
                        self.no_z = False
                        self.no_t = True  
                        self.axes = ['Z', 'Y', 'X']
                    else:
                        self.no_z = False
                        self.no_t = False
                        self.axes = ['T', 'Z', 'Y', 'X']
                else:
                    self.shape = None
                    self.ndim = 0
                    self.no_z = True
                    self.no_t = True
                    self.axes = []
                
                # 创建管道路径字典
                self.pipeline_paths = {}
                base_path = Path(im_path).parent / "temp_processing"
                base_path.mkdir(exist_ok=True)
                
                self.pipeline_paths['im_preprocessed'] = str(base_path / "frangi_filtered.npy")
                self.pipeline_paths['im_instance_label'] = str(base_path / "instance_labels.npy") 
                self.pipeline_paths['im_skel'] = str(base_path / "skeleton.npy")
                self.pipeline_paths['im_pixel_class'] = str(base_path / "pixel_class.npy")
                self.pipeline_paths['im_skel_relabelled'] = str(base_path / "skeleton_relabelled.npy")
                self.pipeline_paths['im_marker'] = str(base_path / "markers.npy")
                self.pipeline_paths['im_distance'] = str(base_path / "distance.npy")
                self.pipeline_paths['im_border'] = str(base_path / "border.npy")
            
            def get_memmap(self, path):
                """获取内存映射数组"""
                if path == self.im_path:
                    return self.processed_image
                elif Path(path).exists():
                    return np.load(path, mmap_mode='r')
                else:
                    return None
            
            def allocate_memory(self, path, dtype='float32', description='', return_memmap=True):
                """分配内存映射数组"""
                if self.processed_image is not None:
                    shape = self.processed_image.shape
                    memmap_array = np.memmap(path, dtype=dtype, mode='w+', shape=shape)
                    if return_memmap:
                        return memmap_array
                return None
        
        return SimpleImInfo(self.im_path, self.processed_image, self.dim_res, self.dimension_order)
    
    def run_frangi_filtering(self) -> bool:
        """
        步骤1: 运行Frangi血管滤波器预处理
        
        Returns
        -------
        bool
            处理是否成功
        """
        try:
            logger.info("=" * 50)
            logger.info("步骤1: 开始Frangi血管滤波器预处理")
            
            if self.Filter is None:
                logger.error("Filter类未正确导入")
                return False
            
            if self.im_info is None:
                logger.error("图像信息未正确加载")
                return False
            
            filter_obj = self.Filter(
                im_info=self.im_info,
                num_t=self.config['num_t'],
                **self.config['frangi']
            )
            
            filter_obj.run(mask=True)
            
            self.results['frangi_filter'] = filter_obj
            logger.info("Frangi滤波完成")
            return True
            
        except Exception as e:
            logger.error(f"Frangi滤波失败: {e}")
            return False
    
    def run_segmentation_labelling(self) -> bool:
        """
        步骤2: 运行分割标记
        
        Returns
        -------
        bool
            处理是否成功
        """
        try:
            logger.info("=" * 50)
            logger.info("步骤2: 开始分割标记")
            
            if self.Label is None:
                logger.error("Label类未正确导入")
                return False
            
            if self.im_info is None:
                logger.error("图像信息未正确加载")
                return False
            
            label_obj = self.Label(
                im_info=self.im_info,
                num_t=self.config['num_t'],
                **self.config['labelling']
            )
            
            label_obj.run()
            
            self.results['labelling'] = label_obj
            logger.info("分割标记完成")
            return True
            
        except Exception as e:
            logger.error(f"分割标记失败: {e}")
            return False
    
    def run_network_analysis(self) -> bool:
        """
        步骤3: 运行网络骨架化分析
        
        Returns
        -------
        bool
            处理是否成功
        """
        try:
            logger.info("=" * 50) 
            logger.info("步骤3: 开始网络骨架化分析")
            
            if self.Network is None:
                logger.error("Network类未正确导入")
                return False
            
            if self.im_info is None:
                logger.error("图像信息未正确加载")
                return False
            
            network_obj = self.Network(
                im_info=self.im_info,
                num_t=self.config['num_t'],
                **self.config['networking']
            )
            
            network_obj.run()
            
            self.results['networking'] = network_obj
            logger.info("网络骨架化分析完成")
            return True
            
        except Exception as e:
            logger.error(f"网络骨架化分析失败: {e}")
            return False
    
    def run_mocap_marking(self) -> bool:
        """
        步骤4: 运行运动捕获标记点检测
        
        Returns
        -------
        bool
            处理是否成功
        """
        try:
            logger.info("=" * 50)
            logger.info("步骤4: 开始运动捕获标记点检测")
            
            if self.Markers is None:
                logger.error("Markers类未正确导入")
                return False
            
            if self.im_info is None:
                logger.error("图像信息未正确加载")
                return False
            
            markers_obj = self.Markers(
                im_info=self.im_info,
                num_t=self.config['num_t'],
                **self.config['markers']
            )
            
            markers_obj.run()
            
            self.results['markers'] = markers_obj
            logger.info("运动捕获标记点检测完成")
            return True
            
        except Exception as e:
            logger.error(f"运动捕获标记点检测失败: {e}")
            return False
    
    def save_results(self) -> bool:
        """
        保存处理结果
        
        Returns
        -------
        bool
            保存是否成功
        """
        try:
            logger.info("=" * 50)
            logger.info("开始保存处理结果")
            
            # 保存配置文件
            import json
            config_file = Path(self.output_dir) / "config.json"
            with open(config_file, 'w', encoding='utf-8') as f:
                json.dump(self.config, f, indent=2, ensure_ascii=False)
            
            # 保存结果摘要
            results_summary = {
                'input_image': self.im_path,
                'output_directory': self.output_dir,
                'image_shape': self.im_info.shape if self.im_info else "未知",
                'processing_steps': list(self.results.keys()),
                'success': len(self.results) == 4  # 所有4个步骤都成功
            }
            
            summary_file = Path(self.output_dir) / "results_summary.json"
            with open(summary_file, 'w', encoding='utf-8') as f:
                json.dump(results_summary, f, indent=2, ensure_ascii=False)
            
            logger.info(f"结果已保存到: {self.output_dir}")
            return True
            
        except Exception as e:
            logger.error(f"保存结果失败: {e}")
            return False
    
    def run_full_pipeline(self) -> bool:
        """
        运行完整的线粒体分割管道
        
        Returns
        -------
        bool
            整个管道是否成功完成
        """
        logger.info("🔬 开始线粒体分割完整管道")
        
        success_steps = []
        
        # 步骤1: Frangi滤波
        if self.run_frangi_filtering():
            success_steps.append("frangi_filtering")
        else:
            logger.error("Frangi滤波失败，停止管道")
            return False
        
        # 步骤2: 分割标记 
        if self.run_segmentation_labelling():
            success_steps.append("segmentation_labelling")
        else:
            logger.error("分割标记失败，停止管道")
            return False
        
        # 步骤3: 网络分析
        if self.run_network_analysis():
            success_steps.append("network_analysis")
        else:
            logger.warning("网络分析失败，继续其他步骤")
        
        # 步骤4: 标记点检测
        if self.run_mocap_marking():
            success_steps.append("mocap_marking")
        else:
            logger.warning("标记点检测失败，继续保存结果")
        
        # 保存结果
        if self.save_results():
            success_steps.append("save_results")
        
        logger.info("=" * 50)
        logger.info(f"✅ 管道完成！成功步骤: {', '.join(success_steps)}")
        logger.info(f"📁 结果保存在: {self.output_dir}")
        
        return len(success_steps) >= 3  # 至少完成前3个关键步骤
    
    def get_segmentation_stats(self) -> Dict:
        """
        获取分割统计信息
        
        Returns
        -------
        dict
            分割统计结果
        """
        stats = {}
        
        if 'labelling' in self.results:
            try:
                label_memmap = self.results['labelling'].instance_label_memmap
                unique_labels = np.unique(label_memmap)
                stats['num_mitochondria'] = len(unique_labels) - 1  # 减去背景标签0
                stats['label_range'] = [int(unique_labels.min()), int(unique_labels.max())]
            except Exception as e:
                logger.warning(f"无法获取分割统计: {e}")
        
        if 'networking' in self.results:
            try:
                skel_memmap = self.results['networking'].skel_memmap
                pixel_class_memmap = self.results['networking'].pixel_class_memmap
                stats['skeleton_pixels'] = int(np.sum(skel_memmap > 0))
                stats['junction_pixels'] = int(np.sum(pixel_class_memmap == 4))
                stats['branch_pixels'] = int(np.sum(pixel_class_memmap == 3))
                stats['endpoint_pixels'] = int(np.sum(pixel_class_memmap == 2))
            except Exception as e:
                logger.warning(f"无法获取网络统计: {e}")
        
        return stats


def run_mitochondria_segmentation(im_path: str,
                                output_dir: Optional[str] = None,
                                dim_res: Optional[Dict[str, float]] = None,
                                channel: Optional[int] = None,
                                **kwargs) -> MitochondriaSegmentationPipeline:
    """
    便捷函数：运行完整的线粒体分割管道
    
    Parameters
    ----------
    im_path : str
        输入图像文件路径
    output_dir : str, optional
        输出目录路径
    dim_res : dict, optional
        像素分辨率
    channel : int, optional
        要处理的通道索引（用于多通道图像）
    **kwargs : dict
        其他配置参数
        
    Returns
    -------
    MitochondriaSegmentationPipeline
        分割管道对象
        
    Example
    -------
    >>> # 处理单通道图像
    >>> pipeline = run_mitochondria_segmentation(
    ...     im_path="path/to/mitochondria.tif",
    ...     dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1},
    ...     frangi={'min_radius_um': 0.1, 'max_radius_um': 0.8}
    ... )
    >>> 
    >>> # 处理多通道SIM图像
    >>> pipeline = run_mitochondria_segmentation(
    ...     im_path="path/to/sim_image.tif",
    ...     dim_res={'Z': 0.125, 'Y': 0.063, 'X': 0.063},
    ...     channel=0,  # 选择第一个通道
    ...     frangi={'min_radius_um': 0.1, 'max_radius_um': 0.8}
    ... )
    >>> stats = pipeline.get_segmentation_stats()
    >>> print(f"检测到 {stats.get('num_mitochondria', 0)} 个线粒体")
    """
    
    # 创建管道对象
    pipeline = MitochondriaSegmentationPipeline(
        im_path=im_path,
        output_dir=output_dir,
        dim_res=dim_res,
        channel=channel,
        **kwargs
    )
    
    # 运行完整管道
    success = pipeline.run_full_pipeline()
    
    if success:
        print("🎉 线粒体分割完成！")
        stats = pipeline.get_segmentation_stats()
        if stats:
            print(f"📊 分割统计:")
            for key, value in stats.items():
                print(f"   {key}: {value}")
    else:
        print("❌ 线粒体分割失败，请检查日志文件")
    
    return pipeline


if __name__ == "__main__":
    # 示例用法
    
    # 基本用法
    im_path = r"F:\2024_06_26_SD_ExM_nhs_u2OS_488+578_cropped.tif"
    
    # 方式1: 使用便捷函数
    pipeline = run_mitochondria_segmentation(
        im_path=im_path,
        dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1},
        frangi={'min_radius_um': 0.15, 'max_radius_um': 1.0},
        labelling={'snr_cleaning': True},
        networking={'clean_skel': True}
    )
    
    # 方式2: 手动创建和运行
    # pipeline = MitochondriaSegmentationPipeline(
    #     im_path=im_path,
    #     dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1}
    # )
    # pipeline.run_full_pipeline()
    
    # 获取结果统计
    stats = pipeline.get_segmentation_stats()
    print("最终统计结果:", stats)

# ============================================================================
# 使用说明和示例
# ============================================================================

"""
## 线粒体分割管道使用指南

### 1. 文件准备
确保以下文件在同一目录下：
- mitochondria_segmentation_pipeline.py (本文件)
- filtering.py (Frangi滤波器)
- labelling.py (分割标记)
- networking.py (网络骨架化)
- mocap_marking.py (运动捕获标记)

### 2. 基本使用

```python
# 简单使用 - 使用默认参数
pipeline = run_mitochondria_segmentation(
    im_path="path/to/your/mitochondria_image.tif"
)

# 高级使用 - 自定义参数
pipeline = run_mitochondria_segmentation(
    im_path="path/to/your/image.tif",
    output_dir="path/to/output/directory",
    dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1},  # 像素分辨率(微米)
    frangi={
        'min_radius_um': 0.15,    # 最小线粒体半径
        'max_radius_um': 1.0,     # 最大线粒体半径
        'alpha_sq': 0.5,          # Frangi参数α²
        'beta_sq': 0.5,           # Frangi参数β²
        'remove_edges': True      # 移除边缘
    },
    labelling={
        'snr_cleaning': True,          # 启用SNR清理
        'otsu_thresh_intensity': False # 不使用Otsu阈值
    },
    networking={
        'clean_skel': True        # 清理骨架
    }
)

# 获取分割统计信息
stats = pipeline.get_segmentation_stats()
print(f"检测到 {stats.get('num_mitochondria', 0)} 个线粒体")
```

### 3. 输出文件说明

处理完成后，输出目录将包含：
- config.json: 处理配置参数
- results_summary.json: 结果摘要
- segmentation.log: 处理日志
- 各种中间处理结果的内存映射文件

### 4. 常见问题

Q: 如何调整参数以获得更好的分割效果？
A: 
- 如果线粒体太小/太大：调整 min_radius_um 和 max_radius_um
- 如果有太多噪声：启用 snr_cleaning
- 如果分割不够精细：降低 alpha_sq 和 beta_sq 值

Q: 处理很慢怎么办？
A: 
- 确保安装了GPU版本的依赖包
- 减少图像尺寸或裁剪感兴趣区域
- 调整 sigma 参数减少多尺度处理步骤

### 5. 支持的图像格式

- TIFF (.tif, .tiff)
- OME-TIFF (.ome.tif)
- 其他PIL支持的格式

### 6. 系统要求

- Python 3.7+
- NumPy, SciPy
- scikit-image
- 可选: CuPy (GPU加速)
"""