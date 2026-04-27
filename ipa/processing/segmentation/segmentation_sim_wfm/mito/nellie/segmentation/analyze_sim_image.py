#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIM图像结构分析测试脚本
用于检查多通道SIM图像的结构和选择合适的处理通道
"""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

def analyze_sim_image(im_path):
    """分析SIM图像结构"""
    print(f"🔍 分析图像: {im_path}")
    print("=" * 60)
    
    try:
        # 尝试使用tifffile加载
        import tifffile
        
        with tifffile.TiffFile(im_path) as tif:
            image_array = tif.asarray()
            
            print(f"📊 图像基本信息:")
            print(f"   形状: {image_array.shape}")
            print(f"   数据类型: {image_array.dtype}")
            print(f"   维度数: {image_array.ndim}")
            print(f"   数据范围: {image_array.min()} - {image_array.max()}")
            print()
            
            # 分析维度结构
            analyze_dimensions(image_array)
            
            # 如果是多通道，分析每个通道
            if analyze_channels(image_array):
                plot_channels(image_array, im_path)
            
            return image_array
            
    except ImportError:
        print("❌ 需要安装tifffile: pip install tifffile")
        return None
    except Exception as e:
        print(f"❌ 读取图像失败: {e}")
        return None

def analyze_dimensions(image_array):
    """分析图像维度结构"""
    shape = image_array.shape
    ndim = len(shape)
    
    print(f"🔧 维度分析:")
    
    if ndim == 2:
        print(f"   2D图像: (H={shape[0]}, W={shape[1]})")
        return 'YX', False, False
        
    elif ndim == 3:
        print(f"   3D图像: {shape}")
        if shape[0] <= 10:
            print(f"   推测结构: (C={shape[0]}, H={shape[1]}, W={shape[2]})")
            return 'CYX', True, False
        elif shape[-1] <= 10:
            print(f"   推测结构: (H={shape[0]}, W={shape[1]}, C={shape[2]})")
            return 'YXC', True, False
        else:
            print(f"   推测结构: (Z={shape[0]}, H={shape[1]}, W={shape[2]})")
            return 'ZYX', False, True
            
    elif ndim == 4:
        print(f"   4D图像: {shape}")
        if shape[0] <= 10:
            print(f"   推测结构: (C={shape[0]}, Z={shape[1]}, H={shape[2]}, W={shape[3]})")
            return 'CZYX', True, True
        elif shape[1] <= 10:
            print(f"   推测结构: (Z={shape[0]}, C={shape[1]}, H={shape[2]}, W={shape[3]})")
            return 'ZCYX', True, True
        else:
            print(f"   推测结构: (T={shape[0]}, Z={shape[1]}, H={shape[2]}, W={shape[3]})")
            return 'TZYX', False, True
    else:
        print(f"   不支持的维度数: {ndim}")
        return 'Unknown', False, False

def analyze_channels(image_array):
    """分析通道信息"""
    shape = image_array.shape
    ndim = len(shape)
    
    # 判断是否有通道
    has_channels = False
    channel_axis = -1
    num_channels = 1
    
    if ndim >= 3:
        if shape[0] <= 10:  # 第一个维度是通道
            has_channels = True
            channel_axis = 0
            num_channels = shape[0]
        elif ndim >= 3 and shape[-1] <= 10:  # 最后一个维度是通道
            has_channels = True  
            channel_axis = -1
            num_channels = shape[-1]
        elif ndim >= 4 and shape[1] <= 10:  # 第二个维度是通道
            has_channels = True
            channel_axis = 1
            num_channels = shape[1]
    
    if has_channels:
        print(f"\n📺 通道信息:")
        print(f"   通道数: {num_channels}")
        print(f"   通道轴: {channel_axis}")
        
        # 分析每个通道的统计信息
        for i in range(num_channels):
            if channel_axis == 0:
                channel_data = image_array[i]
            elif channel_axis == -1:
                channel_data = image_array[..., i]
            elif channel_axis == 1:
                channel_data = image_array[:, i]
            else:
                continue
            
            mean_val = np.mean(channel_data)
            std_val = np.std(channel_data)
            max_val = np.max(channel_data)
            
            print(f"   通道 {i}: 均值={mean_val:.1f}, 标准差={std_val:.1f}, 最大值={max_val}")
        
        return True
    else:
        print(f"\n📺 单通道图像")
        return False

def plot_channels(image_array, im_path):
    """可视化不同通道"""
    shape = image_array.shape
    ndim = len(shape)
    
    # 确定通道结构
    if ndim >= 3 and shape[0] <= 10:
        channel_axis = 0
        num_channels = shape[0]
    elif ndim >= 3 and shape[-1] <= 10:
        channel_axis = -1
        num_channels = shape[-1]
    elif ndim >= 4 and shape[1] <= 10:
        channel_axis = 1
        num_channels = shape[1]
    else:
        return
    
    try:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(1, min(num_channels, 4), figsize=(15, 4))
        if num_channels == 1:
            axes = [axes]
        
        for i in range(min(num_channels, 4)):
            # 提取通道数据
            if channel_axis == 0:
                channel_data = image_array[i]
            elif channel_axis == -1:
                channel_data = image_array[..., i]
            elif channel_axis == 1:
                channel_data = image_array[:, i]
            
            # 如果是3D，取中间切片
            if channel_data.ndim == 3:
                mid_slice = channel_data.shape[0] // 2
                channel_data = channel_data[mid_slice]
            
            # 显示
            ax = axes[i] if num_channels > 1 else axes[0]
            im = ax.imshow(channel_data, cmap='gray')
            ax.set_title(f'通道 {i}')
            ax.axis('off')
            plt.colorbar(im, ax=ax)
        
        plt.tight_layout()
        plt.suptitle(f'SIM图像通道预览: {Path(im_path).name}')
        
        # 保存预览图
        output_path = Path(im_path).parent / f"{Path(im_path).stem}_channels_preview.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n💾 通道预览已保存: {output_path}")
        
        plt.show()
        
    except ImportError:
        print("ℹ️  安装matplotlib可查看通道预览: pip install matplotlib")

def recommend_processing_params(image_array):
    """推荐处理参数"""
    shape = image_array.shape
    
    print(f"\n🎯 推荐处理参数:")
    
    # 推荐分辨率
    if len(shape) >= 3 and shape[0] > 50:  # 可能有Z轴
        print(f"   dim_res: {{'Z': 0.125, 'Y': 0.063, 'X': 0.063}}  # SIM典型分辨率")
    else:
        print(f"   dim_res: {{'Y': 0.063, 'X': 0.063}}  # SIM 2D分辨率")
    
    # 推荐通道
    if len(shape) >= 3 and shape[0] <= 10:
        print(f"   channel: 0  # 选择第一个通道，可尝试其他通道")
    
    # 推荐Frangi参数
    print(f"   frangi参数:")
    print(f"     min_radius_um: 0.1    # 适合高分辨率SIM")
    print(f"     max_radius_um: 0.8    # 线粒体典型尺寸")
    print(f"     alpha_sq: 0.3         # 对细长结构敏感")
    print(f"     beta_sq: 0.3")

def main():
    """主函数"""
    # SIM图像路径
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ 图像文件不存在: {sim_path}")
        return
    
    # 分析图像
    image_array = analyze_sim_image(sim_path)
    
    if image_array is not None:
        recommend_processing_params(image_array)
        
        print(f"\n✅ 分析完成！")
        print(f"📝 建议的处理代码:")
        print(f"""
pipeline = run_mitochondria_segmentation(
    im_path=r"{sim_path}",
    dim_res={{'Z': 0.125, 'Y': 0.063, 'X': 0.063}},
    channel=0,  # 根据预览选择合适通道
    frangi={{
        'min_radius_um': 0.1,
        'max_radius_um': 0.8,
        'alpha_sq': 0.3,
        'beta_sq': 0.3,
        'remove_edges': True
    }},
    labelling={{
        'snr_cleaning': True,
        'otsu_thresh_intensity': False
    }}
)
""")

if __name__ == "__main__":
    main()