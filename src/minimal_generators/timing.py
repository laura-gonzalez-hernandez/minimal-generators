"""
Timing storage and reporting for minimal-generators benchmarks.
"""
import json
import os
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

# Store timing data in user's home directory
TIMING_FILE = '.minimal_generators_timing.json'

def load_timing_data() -> Dict:
    """Load timing data from file."""
    if not os.path.exists(TIMING_FILE):
        return {"runs": {}}
    
    try:
        with open(TIMING_FILE, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return {"runs": {}}

def save_timing_data(data: Dict):
    """Save timing data to file."""
    try:
        with open(TIMING_FILE, 'w') as f:
            json.dump(data, f, indent=2)
    except IOError as e:
        print(f"Warning: Could not save timing data: {e}")

def record_timing(n: int, elapsed_time: float, algorithm: str):
    """
    Record a timing measurement.
    
    Args:
        n: The n value used
        elapsed_time: Execution time in seconds
        algorithm: Either 'solve' or 'gauss_jordan_solve'
    """
    data = load_timing_data()
    
    n_key = str(n)
    if n_key not in data["runs"]:
        data["runs"][n_key] = {}
    
    data["runs"][n_key][algorithm] = {
        "time": elapsed_time,
        "timestamp": datetime.now().isoformat()
    }
    
    save_timing_data(data)

def get_timing_summary() -> Dict[int, Dict[str, float]]:
    """
    Get timing summary for all recorded runs.
    
    Returns:
        Dictionary mapping n values to algorithm timings
        Example: {10: {'solve': 0.5, 'gauss_jordan_solve': 0.3}}
    """
    data = load_timing_data()
    summary = {}
    
    for n_str, algorithms in data["runs"].items():
        n = int(n_str)
        summary[n] = {}
        for algo, info in algorithms.items():
            summary[n][algo] = info["time"]
    
    return summary

def clear_timing_data():
    """Clear all stored timing data."""
    if TIMING_FILE.exists():
        TIMING_FILE.unlink()
