#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIM图像线粒体分割测试脚本
用于测试多通道SIM图像的线粒体分割功能
"""

import sys
import numpy as np
from pathlib import Path

def test_sim_image_loading():
    """测试SIM图像加载和分析"""
    print("🔬 测试SIM图像加载...")
    
    # SIM图像路径
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ 图像文件不存在: {sim_path}")
        return False
    
    try:
        # 尝试加载图像
        import tifffile
        
        print(f"📂 加载图像: {Path(sim_path).name}")
        with tifffile.TiffFile(sim_path) as tif:
            image_array = tif.asarray()
        
        print(f"✅ 图像加载成功!")
        print(f"   形状: {image_array.shape}")
        print(f"   数据类型: {image_array.dtype}")
        print(f"   数据范围: {image_array.min()} - {image_array.max()}")
        
        return True, image_array
        
    except ImportError:
        print("❌ 需要安装tifffile: pip install tifffile")
        return False, None
    except Exception as e:
        print(f"❌ 图像加载失败: {e}")
        return False, None

def test_mitochondria_pipeline():
    """测试线粒体分割管道"""
    print("\n🔧 测试线粒体分割管道...")
    
    try:
        # 尝试导入管道
        sys.path.append(str(Path(__file__).parent))
        from mitochondria_segmentation_pipeline import MitochondriaSegmentationPipeline
        
        print("✅ 成功导入分割管道")
        
        # 创建管道实例
        sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
        
        pipeline = MitochondriaSegmentationPipeline(
            im_path=sim_path,
            dim_res={'Z': 0.125, 'Y': 0.063, 'X': 0.063},
            channel=0  # 尝试第一个通道
        )
        
        print("✅ 成功创建分割管道实例")
        
        # 测试图像分析
        if hasattr(pipeline, 'processed_image') and pipeline.processed_image is not None:
            print(f"✅ 图像预处理成功，处理后尺寸: {pipeline.processed_image.shape}")
        else:
            print("⚠️  图像预处理可能有问题")
        
        return True, pipeline
        
    except ImportError as e:
        print(f"❌ 导入失败: {e}")
        return False, None
    except Exception as e:
        print(f"❌ 管道创建失败: {e}")
        import traceback
        traceback.print_exc()
        return False, None

def test_simple_processing():
    """测试简单处理流程"""
    print("\n🚀 测试简单处理流程...")
    
    success, pipeline = test_mitochondria_pipeline()
    
    if not success:
        return False
    
    try:
        # 测试单个步骤
        print("📍 测试步骤1: Frangi滤波...")
        if pipeline.run_frangi_filtering():
            print("  ✅ Frangi滤波成功")
        else:
            print("  ❌ Frangi滤波失败")
            return False
        
        print("📍 测试步骤2: 分割标记...")
        if pipeline.run_segmentation_labelling():
            print("  ✅ 分割标记成功")
        else:
            print("  ❌ 分割标记失败")
            return False
        
        # 获取统计信息
        stats = pipeline.get_segmentation_stats()
        print(f"\n📊 分割结果统计:")
        for key, value in stats.items():
            print(f"   {key}: {value}")
        
        return True
        
    except Exception as e:
        print(f"❌ 处理失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主测试函数"""
    print("🧪 SIM图像线粒体分割测试")
    print("=" * 50)
    
    # 测试1: 图像加载
    load_success, image_array = test_sim_image_loading()
    
    if not load_success:
        print("\n❌ 基础测试失败，请检查图像文件和依赖")
        return
    
    # 测试2: 管道创建
    print("\n" + "="*50)
    pipeline_success = test_simple_processing()
    
    if pipeline_success:
        print("\n🎉 所有测试通过！")
        print("💡 现在可以使用完整的分割功能了")
    else:
        print("\n⚠️  部分测试失败，请检查依赖和配置")
    
    print(f"\n📝 下一步:")
    print(f"1. 运行 analyze_sim_image.py 查看图像详细信息") 
    print(f"2. 运行 example_usage.py 进行完整分割")
    print(f"3. 根据结果调整参数")

if __name__ == "__main__":
    main()
    input("\n按Enter键退出...")