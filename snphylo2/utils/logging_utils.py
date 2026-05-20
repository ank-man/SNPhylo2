"""
Logging utilities for SNPhylo2 using Loguru.
"""

import sys
from pathlib import Path
from typing import Optional

from loguru import logger


def setup_logging(
    verbose: bool = False,
    quiet: bool = False,
    log_file: Optional[Path] = None,
    rotation: str = "10 MB"
) -> None:
    """
    Configure logging for SNPhylo2.
    
    Args:
        verbose: Enable DEBUG level logging
        quiet: Only show ERROR level and above
        log_file: Optional file path for log output
        rotation: Log file rotation size
    """
    # Remove default handler
    logger.remove()
    
    # Determine log level
    if quiet:
        level = "ERROR"
        console_format = "<red>{message}</red>"
    elif verbose:
        level = "DEBUG"
        console_format = (
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
            "<level>{level: <8}</level> | "
            "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
            "<level>{message}</level>"
        )
    else:
        level = "INFO"
        console_format = "<level>{message}</level>"
    
    # Add console handler
    logger.add(
        sys.stderr,
        level=level,
        format=console_format,
        colorize=True,
        enqueue=True,
    )
    
    # Add file handler if requested
    if log_file:
        logger.add(
            log_file,
            level="DEBUG",
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}",
            rotation=rotation,
            enqueue=True,
        )


def get_logger():
    """Get the configured logger instance."""
    return logger
