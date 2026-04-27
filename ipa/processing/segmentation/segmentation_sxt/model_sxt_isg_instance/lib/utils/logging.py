import os
import sys
import logging
import time


# 创建默认logger
logger = logging.getLogger('upsnet')
logger.setLevel(logging.INFO)

# 创建控制台handler（如果还没有的话）
if not logger.handlers:
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


def create_logger(output_dir, cfg_name, image_set):
    """
    创建日志记录器
    
    Args:
        output_dir: 输出目录
        cfg_name: 配置文件名
        image_set: 图像集名称
    
    Returns:
        logger: 日志记录器
        final_output_path: 最终输出路径
    """
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 创建时间戳
    time_str = time.strftime('%Y-%m-%d-%H-%M')
    
    # 创建最终输出路径
    cfg_name = cfg_name.split('.')[0] if cfg_name else 'default'
    final_output_path = os.path.join(output_dir, f'{cfg_name}_{image_set}_{time_str}')
    
    if not os.path.exists(final_output_path):
        os.makedirs(final_output_path)
    
    # 创建日志文件路径
    log_file = os.path.join(final_output_path, f'log_{time_str}.txt')
    
    # 创建logger
    train_logger = logging.getLogger('upsnet')
    train_logger.setLevel(logging.INFO)
    
    # 清除已有的handler
    for handler in train_logger.handlers[:]:
        train_logger.removeHandler(handler)
    
    # 创建文件handler
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    
    # 创建控制台handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    
    # 创建formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # 添加handler到logger
    train_logger.addHandler(file_handler)
    train_logger.addHandler(console_handler)
    
    train_logger.info(f'Created logger with output path: {final_output_path}')
    train_logger.info(f'Log file: {log_file}')
    
    return train_logger, final_output_path