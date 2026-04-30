import logging
import functools
from datetime import datetime
from pathlib import Path
import inspect
import sys

def setup_simple_logger(func_name, log_dir="logs", console_output=False):
    """
    Quick setup for logger with timestamped file output.
    
    Creates a logger that writes to a timestamped file in the specified directory.
    The logger is configured with INFO level and UTF-8 encoding for international
    character support.
    
    Args:
        func_name (str): Name of the function or analysis, used in log filename
        log_dir (str or Path): Directory to store log files. Defaults to "logs"
        console_output (bool): Whether to also print logs to console. Defaults to False.
    
    Returns:
        logging.Logger: Configured logger instance with file handler
        
    Example:
        >>> logger = setup_simple_logger("data_analysis", "analysis_logs")
        >>> logger.info("Starting analysis")
        # Creates: analysis_logs/data_analysis_20240115_143025.log
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{func_name}_{timestamp}.log"
    
    logger = logging.getLogger(f"{func_name}_{timestamp}")
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers to avoid duplicates if called multiple times
    logger.handlers.clear()
    
    # File Handler
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    
    # Console Handler (Optional)
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = logging.Formatter('%(levelname)s: %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
    
    return logger

def log_analysis(log_dir="logs"):
    """
    Decorator: Automatically log function execution process.
    
    This decorator provides automatic logging for data analysis functions.
    It logs function start/completion, automatically detects and logs file paths
    from function parameters, handles exceptions, and logs output files.
    
    Args:
        log_dir (str or Path): Directory to store log files. Defaults to "logs"
    
    Returns:
        function: Decorated function with automatic logging capabilities
        
        
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Setup logger
            logger = setup_simple_logger(func.__name__, log_dir)
            
            # Start logging
            logger.info(f"Starting execution: {func.__name__}")
            
            # Log parameters (if file paths exist)
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            
            for param_name, param_value in bound_args.arguments.items():
                if isinstance(param_value, (str, Path)) and ('path' in param_name.lower() or 'file' in param_name.lower()):
                    logger.info(f"Input path: {param_value}")
                elif param_name in ['output_dir', 'save_dir', 'result_path']:
                    logger.info(f"Output path: {param_value}")
            
            try:
                # Execute function
                result = func(*args, **kwargs)
                logger.info(f"Completed execution: {func.__name__}")
                
                # If result contains file paths, log them
                if isinstance(result, (str, Path)):
                    logger.info(f"Generated file: {result}")
                elif isinstance(result, (list, tuple)) and result and isinstance(result[0], (str, Path)):
                    for file_path in result:
                        logger.info(f"Generated file: {file_path}")
                
                return result
                
            except Exception as e:
                logger.error(f"Execution failed: {func.__name__} - {str(e)}")
                raise
                
        return wrapper
    return decorator

# Simpler version - manual logging
class QuickLogger:
    """
    Simple interface for step-by-step logging in data analysis workflows.
    
    Provides predefined methods for common logging scenarios including
    analysis steps, input/output files, and error handling. Creates
    timestamped log files automatically.
    
    Args:
        name (str): Logger name, used in log filename. Defaults to "analysis"
        log_dir (str or Path): Directory to store log files. Defaults to "logs"
        console_output (bool): Whether to print logs to console. Defaults to True.
        
    Attributes:
        logger (logging.Logger): Internal logger instance
        

    """
    def __init__(self, name="analysis", log_dir="logs", console_output=True):
        self.logger = setup_simple_logger(name, log_dir, console_output=console_output)
    
    def step(self, message):
        """
        Log analysis step or general information.
        
        Args:
            message (str): Step description or information message
            
        Example:
            >>> logger.step("Loading mask data...")
            >>> logger.step(f"Processing dataset: {dataid}")
        """
        self.logger.info(message)
    
    def file_in(self, file_path):
        """
        Log input file path.
        
        Args:
            file_path (str or Path): Path to input file
            
        Example:
            >>> logger.file_in("/data/images/sample.tif")
            >>> logger.file_in(mask_pm_path)
        """
        self.logger.info(f"Input file: {file_path}")
    
    def file_out(self, file_path):
        """
        Log output file path.
        
        Args:
            file_path (str or Path): Path to generated output file
            
        Example:
            >>> logger.file_out("/results/processed.tif")
            >>> logger.file_out(xvg_path)
        """
        self.logger.info(f"Output file: {file_path}")
    
    def error(self, message):
        """
        Log error message.
        
        Args:
            message (str): Error description
            
        Example:
            >>> logger.error("Failed to load input file")
            >>> logger.error(f"Processing failed: {str(e)}")
        """
        self.logger.error(message)
