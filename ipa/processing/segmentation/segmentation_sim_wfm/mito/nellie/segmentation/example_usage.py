#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
线粒体分割管道简单示例
演示如何使用线粒体分割功能
"""

import os
from pathlib import Path

# 导入分割管道
from mitochondria_segmentation_pipeline import run_mitochondria_segmentation, MitochondriaSegmentationPipeline


def example_basic_usage():
    """基本使用示例"""
    print("=== 基本使用示例 ===")
    
    # 设置图像路径（SIM多通道图像）
    im_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    # 检查文件是否存在
    if not Path(im_path).exists():
        print(f"错误：图像文件不存在: {im_path}")
        print("请修改 im_path 为你的实际图像文件路径")
        return
    
    try:
        # 运行分割（SIM图像专用配置）
        pipeline = run_mitochondria_segmentation(
            im_path=im_path,
            dim_res={'Z': 0.125, 'Y': 0.063, 'X': 0.063},  # 典型SIM分辨率
            channel=0,  # 选择第一个通道进行处理
            dimension_order='ZYX',  # 假设已经选择了单通道
            
            # SIM图像优化参数
            frangi={
                'min_radius_um': 0.1,     # 较小的最小半径适合高分辨率
                'max_radius_um': 0.8,     # 适合线粒体尺寸
                'alpha_sq': 0.3,          # 对细长结构更敏感
                'beta_sq': 0.3,           
                'remove_edges': True
            },
            
            labelling={
                'snr_cleaning': True,     # SIM图像噪声较多，启用SNR清理
                'otsu_thresh_intensity': False
            }
        )
        
        # 获取分割统计
        stats = pipeline.get_segmentation_stats()
        print("\n分割结果统计:")
        for key, value in stats.items():
            print(f"  {key}: {value}")
            
    except Exception as e:
        print(f"分割失败: {e}")
        import traceback
        traceback.print_exc()


def example_advanced_usage():
    """高级使用示例 - 自定义参数"""
    print("\n=== 高级使用示例 ===")
    
    # 图像路径
    im_path = r"F:\2024_06_26_SD_ExM_nhs_u2OS_488+578_cropped.tif"
    
    if not Path(im_path).exists():
        print(f"错误：图像文件不存在: {im_path}")
        return
    
    # 自定义输出目录
    output_dir = Path(im_path).parent / "custom_mitochondria_results"
    
    try:
        # 手动创建管道对象
        pipeline = MitochondriaSegmentationPipeline(
            im_path=im_path,
            output_dir=str(output_dir),
            dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1},  # SIM显微镜分辨率
            
            # Frangi滤波参数
            frangi={
                'min_radius_um': 0.1,     # 较小的最小半径
                'max_radius_um': 0.8,     # 较小的最大半径  
                'alpha_sq': 0.3,          # 更敏感的α参数
                'beta_sq': 0.3,           # 更敏感的β参数
                'remove_edges': True
            },
            
            # 分割参数
            labelling={
                'snr_cleaning': True,          # 启用SNR清理
                'otsu_thresh_intensity': False # 不使用Otsu阈值
            },
            
            # 网络分析参数
            networking={
                'clean_skel': True,       # 清理骨架
                'min_radius_um': 0.1,
                'max_radius_um': 0.8
            },
            
            # 标记点参数
            markers={
                'use_im': 'distance',     # 使用距离变换检测峰值
                'num_sigma': 5            # 多尺度σ步数
            }
        )
        
        # 逐步运行管道
        print("步骤1: Frangi滤波...")
        if pipeline.run_frangi_filtering():
            print("  ✓ Frangi滤波完成")
        else:
            print("  ✗ Frangi滤波失败")
            return
        
        print("步骤2: 分割标记...")
        if pipeline.run_segmentation_labelling():
            print("  ✓ 分割标记完成")
        else:
            print("  ✗ 分割标记失败")
            return
        
        print("步骤3: 网络骨架化...")
        if pipeline.run_network_analysis():
            print("  ✓ 网络骨架化完成")
        else:
            print("  ⚠ 网络骨架化失败，继续下一步")
        
        print("步骤4: 标记点检测...")
        if pipeline.run_mocap_marking():
            print("  ✓ 标记点检测完成")
        else:
            print("  ⚠ 标记点检测失败，继续下一步")
        
        # 保存结果
        pipeline.save_results()
        
        # 显示统计信息
        stats = pipeline.get_segmentation_stats()
        print(f"\n🎉 处理完成！结果保存在: {output_dir}")
        print("\n分割结果统计:")
        for key, value in stats.items():
            print(f"  {key}: {value}")
            
    except Exception as e:
        print(f"处理失败: {e}")


def example_batch_processing():
    """批量处理示例"""
    print("\n=== 批量处理示例 ===")
    
    # 图像文件夹路径（请修改为你的实际路径）
    image_folder = Path(r"F:\your_images_folder")
    
    if not image_folder.exists():
        print(f"错误：图像文件夹不存在: {image_folder}")
        print("请修改 image_folder 为你的实际图像文件夹路径")
        return
    
    # 查找所有TIFF文件
    image_files = list(image_folder.glob("*.tif")) + list(image_folder.glob("*.tiff"))
    
    if not image_files:
        print("未找到TIFF图像文件")
        return
    
    print(f"找到 {len(image_files)} 个图像文件")
    
    # 批量处理
    results = []
    for i, im_path in enumerate(image_files[:3]):  # 只处理前3个文件作为示例
        print(f"\n处理文件 {i+1}/{len(image_files)}: {im_path.name}")
        
        try:
            # 为每个文件创建单独的输出目录
            output_dir = image_folder / f"results_{im_path.stem}"
            
            pipeline = run_mitochondria_segmentation(
                im_path=str(im_path),
                output_dir=str(output_dir),
                dim_res={'Z': 0.2, 'Y': 0.1, 'X': 0.1}
            )
            
            stats = pipeline.get_segmentation_stats()
            results.append({
                'file': im_path.name,
                'stats': stats,
                'success': len(pipeline.results) >= 2  # 至少完成前2步
            })
            
        except Exception as e:
            print(f"  ❌ 处理失败: {e}")
            results.append({
                'file': im_path.name,
                'stats': {},
                'success': False
            })
    
    # 显示批量处理结果
    print("\n" + "="*50)
    print("批量处理结果摘要:")
    successful = sum(1 for r in results if r['success'])
    print(f"成功处理: {successful}/{len(results)} 个文件")
    
    for result in results:
        status = "✓" if result['success'] else "✗"
        mito_count = result['stats'].get('num_mitochondria', 0)
        print(f"  {status} {result['file']}: {mito_count} 个线粒体")


if __name__ == "__main__":
    print("🔬 线粒体分割管道示例")
    print("请确保以下文件在同一目录:")
    print("- filtering.py")  
    print("- labelling.py")
    print("- networking.py")
    print("- mocap_marking.py")
    print()
    
    # 运行示例
    try:
        # 基本使用
        example_basic_usage()
        
        # 高级使用 
        example_advanced_usage()
        
        # 批量处理（如果需要）
        # example_batch_processing()
        
    except KeyboardInterrupt:
        print("\n用户中断处理")
    except Exception as e:
        print(f"\n示例运行失败: {e}")
        print("请检查图像路径和依赖文件是否正确")
    
    input("\n按Enter键退出...")