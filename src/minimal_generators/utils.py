"""Utility functions and classes for the minimal_generators package."""

import sys


class Colors:
    """ANSI color codes for terminal output."""
    RESET = '\033[0m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    
    @staticmethod
    def is_terminal():
        """Check if output is going to a terminal (supports colors)."""
        return sys.stdout.isatty()


def log_error(msg: str):
    """Print error message to stderr in red.
    
    Args:
        msg: Error message to print
    """
    if Colors.is_terminal():
        print(f"{Colors.RED}{Colors.BOLD}[ERROR]{Colors.RESET} {msg}", file=sys.stderr)
    else:
        print(f"[ERROR] {msg}", file=sys.stderr)


def log_info(msg: str):
    """Print info message to stdout in purple/magenta.
    
    Args:
        msg: Info message to print
    """
    if Colors.is_terminal():
        print(f"{Colors.MAGENTA}[INFO]{Colors.RESET} {msg}")
    else:
        print(f"[INFO] {msg}")


def log_warning(msg: str):
    """Print warning message to stdout in yellow.
    
    Args:
        msg: Warning message to print
    """
    if Colors.is_terminal():
        print(f"{Colors.YELLOW}[WARNING]{Colors.RESET} {msg}")
    else:
        print(f"[WARNING] {msg}")


def log_success(msg: str):
    """Print success message to stdout in green with bold.
    
    Args:
        msg: Success message to print
    """
    if Colors.is_terminal():
        print(f"{Colors.GREEN}{Colors.BOLD}[SUCCESS]{Colors.RESET} {msg}")
    else:
        print(f"[SUCCESS] {msg}")


def log_debug(msg: str, enabled: bool = False):
    """Print debug message to stdout in cyan (only if enabled).
    
    Args:
        msg: Debug message to print
        enabled: Whether debug logging is enabled
    """
    if not enabled:
        return
    if Colors.is_terminal():
        print(f"{Colors.CYAN}[DEBUG]{Colors.RESET} {msg}")
    else:
        print(f"[DEBUG] {msg}")


def log_progress(bar: str, processed: int, total: int, pct: float, elapsed: float, stage: str):
    """Print progress bar with colored INFO tag.
    
    Args:
        bar: The progress bar string (e.g., '████░░░')
        processed: Number of steps processed
        total: Total number of steps
        pct: Percentage complete (0.0 to 1.0)
        elapsed: Elapsed time in seconds
        stage: Current stage description
    """
    # Always use sys.stdout.write for better control over line endings
    if Colors.is_terminal():
        # Colored output: magenta [INFO], grey bar, green percentage
        msg = f"\r{Colors.MAGENTA}[INFO]{Colors.RESET} [{bar}] {processed}/{total} ({Colors.GREEN}{pct*100:5.1f}%{Colors.RESET}) [{elapsed:6.2f}s] {stage:<20}"
    else:
        msg = f"\r[INFO] [{bar}] {processed}/{total} ({pct*100:5.1f}%) [{elapsed:6.2f}s] {stage:<20}"
    
    sys.stdout.write(msg)
    sys.stdout.flush()


