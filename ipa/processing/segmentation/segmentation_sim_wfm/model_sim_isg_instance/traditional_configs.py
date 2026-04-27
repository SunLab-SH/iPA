# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sim_wfm\model_sim_isg_instance\traditional_configs.py
"""
传统ISG分割方法的配置参数
针对不同的ISG形态特征优化参数
"""

# 基础配置
BASE_CONFIG = {
    # 预处理参数
    'gaussian_blur_kernel': 3,
    'median_filter_size': 3,
    
    # 阈值分割参数
    'threshold_method': 'otsu',
    'manual_threshold': 128,
    'adaptive_block_size': 11,
    'adaptive_c': 2,
    
    # Watershed参数
    'watershed_connectivity': 2,
    'min_distance': 10,
    'watershed_threshold': 0.3,
    
    # 形态学操作参数
    'opening_kernel_size': 3,
    'closing_kernel_size': 5,
    'erosion_iterations': 1,
    'dilation_iterations': 2,
    
    # 后处理参数
    'min_area': 50,
    'max_area': 10000,
    'min_solidity': 0.3,
    'aspect_ratio_range': (0.3, 3.0),
    'circularity_threshold': 0.4,
    'use_shape_filter': True,
}

# 针对小ISG优化的配置
SMALL_ISG_CONFIG = {
    **BASE_CONFIG,
    'gaussian_blur_kernel': 1,  # 减少模糊保留细节
    'median_filter_size': 1,
    'min_distance': 5,  # 降低最小距离检测更小目标
    'watershed_threshold': 0.2,  # 降低阈值
    'opening_kernel_size': 2,  # 减小核大小
    'closing_kernel_size': 3,
    'min_area': 25,  # 降低最小面积
    'max_area': 2000,
    'circularity_threshold': 0.3,
}

# 针对大ISG优化的配置
LARGE_ISG_CONFIG = {
    **BASE_CONFIG,
    'gaussian_blur_kernel': 5,  # 增加模糊去除噪声
    'median_filter_size': 5,
    'min_distance': 20,  # 增加最小距离
    'watershed_threshold': 0.4,
    'opening_kernel_size': 5,  # 增大核大小
    'closing_kernel_size': 7,
    'min_area': 200,  # 提高最小面积
    'max_area': 50000,
    'circularity_threshold': 0.5,
}

# 高精度配置（更严格的形状过滤）
HIGH_PRECISION_CONFIG = {
    **BASE_CONFIG,
    'gaussian_blur_kernel': 3,
    'watershed_threshold': 0.35,
    'min_area': 75,
    'max_area': 5000,
    'min_solidity': 0.5,  # 提高凸度要求
    'aspect_ratio_range': (0.5, 2.0),  # 更严格的长宽比
    'circularity_threshold': 0.5,  # 提高圆形度要求
    'use_shape_filter': True,
}

# 高召回配置（宽松的过滤条件）
HIGH_RECALL_CONFIG = {
    **BASE_CONFIG,
    'gaussian_blur_kernel': 2,
    'watershed_threshold': 0.25,  # 降低阈值提高召回
    'min_area': 30,  # 降低最小面积
    'max_area': 20000,  # 提高最大面积
    'min_solidity': 0.2,  # 降低凸度要求
    'aspect_ratio_range': (0.2, 4.0),  # 更宽松的长宽比
    'circularity_threshold': 0.25,  # 降低圆形度要求
}

# 自适应阈值专用配置
ADAPTIVE_THRESHOLD_CONFIG = {
    **BASE_CONFIG,
    'threshold_method': 'adaptive',
    'adaptive_block_size': 15,  # 适应局部变化
    'adaptive_c': 3,
    'gaussian_blur_kernel': 2,
    'median_filter_size': 3,
}

# Watershed专用配置
WATERSHED_CONFIG = {
    **BASE_CONFIG,
    'min_distance': 12,
    'watershed_threshold': 0.3,
    'opening_kernel_size': 4,
    'closing_kernel_size': 6,
    'erosion_iterations': 2,
    'dilation_iterations': 3,
}

# 组合方法配置
COMBINED_CONFIG = {
    **BASE_CONFIG,
    'gaussian_blur_kernel': 3,
    'median_filter_size': 3,
    'threshold_method': 'otsu',
    'min_distance': 8,
    'watershed_threshold': 0.3,
    'opening_kernel_size': 3,
    'closing_kernel_size': 5,
    'min_area': 50,
    'max_area': 8000,
    'min_solidity': 0.4,
    'aspect_ratio_range': (0.4, 2.5),
    'circularity_threshold': 0.4,
}

# 配置字典
CONFIGS = {
    'base': BASE_CONFIG,
    'small_isg': SMALL_ISG_CONFIG,
    'large_isg': LARGE_ISG_CONFIG,
    'high_precision': HIGH_PRECISION_CONFIG,
    'high_recall': HIGH_RECALL_CONFIG,
    'adaptive': ADAPTIVE_THRESHOLD_CONFIG,
    'watershed': WATERSHED_CONFIG,
    'combined': COMBINED_CONFIG,
}

def get_config(config_name='base'):
    """获取指定配置"""
    if config_name in CONFIGS:
        return CONFIGS[config_name].copy()
    else:
        print(f"Warning: Unknown config '{config_name}', using base config")
        return BASE_CONFIG.copy()

def list_configs():
    """列出所有可用配置"""
    print("Available configurations:")
    for name, config in CONFIGS.items():
        print(f"  {name}: {config.get('description', 'No description')}")
    return list(CONFIGS.keys())