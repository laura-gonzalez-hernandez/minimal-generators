import argparse
import time
import sys
import os
import subprocess
from typing import List
from .generators import get_minimal_generators, get_minimal_generators_new
from .utils import log_error, log_info, log_warning, log_success, log_debug
from .timing import record_timing, get_timing_summary, clear_timing_data

def compare_polynomial_files(n: int):
    """Compare OLD and NEW polynomial files for a given n."""
    out_dir = 'output_min_generators'
    current_file = os.path.join(out_dir, f'minimal_generators_n{n}.txt')
    new_file = os.path.join(out_dir, f'minimal_generators_n{n}_NEW.txt')
    
    # Check if files exist
    if not os.path.exists(current_file):
        log_error(f"File not found: {current_file}")
        return False
    if not os.path.exists(new_file):
        log_error(f"File not found: {new_file}")
        return False
    
    log_info(f"Comparing files for n={n}:")
    log_info(f"  CURRENT: {current_file}")
    log_info(f"  NEW: {new_file}")
    
    # Read both files
    with open(current_file, 'r', encoding='utf-8') as f:
        current_lines = f.readlines()
    with open(new_file, 'r', encoding='utf-8') as f:
        new_lines = f.readlines()
    
    # Compare number of lines
    if len(current_lines) != len(new_lines):
        log_error(f"Different number of polynomials: CURRENT has {len(current_lines)}, NEW has {len(new_lines)}")
        return False
    
    # Compare line by line
    differences = []
    for i, (current_file_line, new_line) in enumerate(zip(current_lines, new_lines), start=1):
        if current_file_line.strip() != new_line.strip():
            differences.append(i)
    
    if differences:
        log_error(f"Found {len(differences)} difference(s) at line(s): {differences}")
        # Show first few differences
        for line_num in differences[:5]:  # Show max 5 differences
            log_info(f"\nLine {line_num}:")
            log_info(f"  CURRENT: {current_lines[line_num-1].strip()}")
            log_info(f"  NEW: {new_lines[line_num-1].strip()}")
        if len(differences) > 5:
            log_info(f"\n... and {len(differences) - 5} more difference(s)")
        return False
    else:
        log_success(f"✓ Files are identical! ({len(current_lines)} polynomials match)")
        return True

def generate_report(output_file: str = 'benchmark_results.png'):
    """Generate performance plot from stored timing data."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        log_error("matplotlib is required for generating reports. Install with: uv add matplotlib")
        return False
    
    summary = get_timing_summary()
    
    if not summary:
        log_warning("No timing data available. Run some benchmarks first!")
        return False
    
    # Organize data
    n_values = sorted(summary.keys())
    times_solve = []
    times_gauss = []
    
    for n in n_values:
        times_solve.append(summary[n].get('solve', None))
        times_gauss.append(summary[n].get('gauss_jordan_solve', None))
    
    # Check if we have any data
    has_solve = any(t is not None for t in times_solve)
    has_gauss = any(t is not None for t in times_gauss)
    
    if not has_solve and not has_gauss:
        log_warning("No timing data available for any algorithm.")
        return False
    
    # Create plot
    plt.figure(figsize=(12, 7))
    
    if has_solve:
        n_solve = [n for n, t in zip(n_values, times_solve) if t is not None]
        t_solve = [t for t in times_solve if t is not None]
        plt.plot(n_solve, t_solve, 'ro-', linewidth=2, markersize=8, label='solve')
        
        # Add value labels
        for n, t in zip(n_solve, t_solve):
            plt.annotate(f'{t:.1f}s', 
                        xy=(n, t), 
                        xytext=(0, 10),
                        textcoords='offset points',
                        ha='center',
                        fontsize=9,
                        color='red')
    
    if has_gauss:
        n_gauss = [n for n, t in zip(n_values, times_gauss) if t is not None]
        t_gauss = [t for t in times_gauss if t is not None]
        plt.plot(n_gauss, t_gauss, 'bo-', linewidth=2, markersize=8, label='gauss_jordan_solve')
        
        # Add value labels
        for n, t in zip(n_gauss, t_gauss):
            plt.annotate(f'{t:.1f}s', 
                        xy=(n, t), 
                        xytext=(0, -15),
                        textcoords='offset points',
                        ha='center',
                        fontsize=9,
                        color='blue')
    
    plt.xlabel('n', fontsize=12)
    plt.ylabel('Execution Time (seconds)', fontsize=12)
    plt.title('minimal-generators Performance Comparison', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11, loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot - automatically detect format from file extension
    # SVG format doesn't use dpi parameter, PNG does
    if output_file.lower().endswith('.svg'):
        plt.savefig(output_file, format='svg')
    else:
        plt.savefig(output_file, dpi=300)
    log_success(f"Plot saved to: {output_file}")
    
    # Print summary table
    print("\n" + "=" * 60)
    print("Timing Summary:")
    print("=" * 60)
    print(f"{'n':>5} {'solve':>12} {'gauss_jordan_solve':>20} {'Speedup':>10}")
    print("-" * 60)
    for n in n_values:
        t_solve = summary[n].get('solve', None)
        t_gauss = summary[n].get('gauss_jordan_solve', None)
        
        solve_str = f"{t_solve:10.2f}s" if t_solve is not None else "     -"
        gauss_str = f"{t_gauss:18.2f}s" if t_gauss is not None else "          -"
        
        if t_solve is not None and t_gauss is not None:
            speedup = t_solve / t_gauss
            speedup_str = f"{speedup:9.2f}x"
        else:
            speedup_str = "     -"
        
        print(f"{n:5d} {solve_str} {gauss_str} {speedup_str}")
    
    return True

def run_benchmark(n_values: List[int], algorithms: List[str]):
    """Run benchmarks for multiple n values."""
    log_info(f"Running benchmarks for n = {n_values}")
    log_info(f"Algorithms: {', '.join(algorithms)}")
    print("=" * 60)
    
    for n in n_values:
        if n < 6:
            log_warning(f"Skipping n={n} (minimum is 6)")
            continue
        
        for algo in algorithms:
            use_new_algo = (algo == 'gauss_jordan_solve')
            algo_display = 'gauss_jordan_solve' if use_new_algo else 'solve'
            
            log_info(f"Running n={n} with {algo_display}...")
            
            # Build command
            cmd = ['minimal-generators', str(n)]
            if use_new_algo:
                cmd.append('--new-algo')
            
            try:
                start = time.perf_counter()
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                elapsed = time.perf_counter() - start
                
                # Record timing
                record_timing(n, elapsed, algo_display)
                
                log_success(f"  Completed in {elapsed:.2f}s")
                
            except subprocess.CalledProcessError as e:
                log_error(f"  Failed: {e.stderr}")
                continue
            except FileNotFoundError:
                log_error("  minimal-generators command not found. Make sure it's installed.")
                return False
        
        print()
    
    log_success("Benchmark complete! Use 'minimal-generators report' to generate the plot.")
    return True

def parse_args():
    # Check for backward compatibility: if first arg is a number, treat as old-style
    import sys as sys_module
    if len(sys_module.argv) > 1 and sys_module.argv[1].lstrip('-').isdigit():
        # Old-style command: minimal-generators <n> [options]
        p = argparse.ArgumentParser(description="Compute minimal generators polynomials symbolically.")
        p.add_argument('n', type=int, help='Positive integer n (minimum: 5)')
        p.add_argument('--singular', action='store_true', default=False, 
                       help='Print using singular-like format (default: False)')
        p.add_argument('--show-polys', action='store_true', default=False,
                       help='Print polynomials via terminal apart from storing them (default: False)')
        p.add_argument('--debug', action='store_true', default=False,
                       help='Show debug messages (default: False)')
        p.add_argument('--new-algo', action='store_true', default=False,
                       help='Use optimized symbolic solver variant (default: False)')
        args = p.parse_args()
        args.command = 'generate'
        return args
    
    # New-style with subcommands
    p = argparse.ArgumentParser(description="Compute minimal generators polynomials symbolically.")
    subparsers = p.add_subparsers(dest='command', help='Available commands')
    
    # Generate command (default)
    gen_parser = subparsers.add_parser('generate', help='Generate minimal generators')
    gen_parser.add_argument('n', type=int, help='Positive integer n (minimum: 5)')
    gen_parser.add_argument('--singular', action='store_true', default=False,
                            help='Print using singular-like format (default: False)')
    gen_parser.add_argument('--show-polys', action='store_true', default=False,
                            help='Print polynomials via terminal apart from storing them (default: False)')
    gen_parser.add_argument('--debug', action='store_true', default=False,
                            help='Show debug messages (default: False)')
    gen_parser.add_argument('--new-algo', action='store_true', default=False,
                            help='Use optimized symbolic solver variant (default: False)')
    
    # Compare command
    compare_parser = subparsers.add_parser('compare', help='Compare OLD and NEW polynomial files')
    compare_parser.add_argument('n', type=int, help='Positive integer n for which files were generated')
    
    # Report command
    report_parser = subparsers.add_parser('report', help='Generate performance plot from stored timing data')
    report_parser.add_argument('--output', type=str, default='benchmark_results.svg',
                              help='Output filename for the plot (default: benchmark_results.svg)')
    report_parser.add_argument('--clear', action='store_true',
                              help='Clear all stored timing data')
    
    # Benchmark command
    benchmark_parser = subparsers.add_parser('benchmark', help='Run benchmarks for multiple n values')
    benchmark_parser.add_argument('n_values', type=int, nargs='+',
                                 help='List of n values to benchmark (e.g., 10 20 40 60)')
    benchmark_parser.add_argument('--algorithms', nargs='+', default=['solve', 'gauss_jordan_solve'],
                                 choices=['solve', 'gauss_jordan_solve'],
                                 help='Algorithms to benchmark (default: both)')
    
    return p.parse_args()

def main():
    args = parse_args()

    if args.command == 'compare':
        success = compare_polynomial_files(args.n)
        sys.exit(0 if success else 1)
    
    if args.command == 'report':
        if args.clear:
            clear_timing_data()
            log_success("Timing data cleared.")
            sys.exit(0)
        success = generate_report(args.output)
        sys.exit(0 if success else 1)
    
    if args.command == 'benchmark':
        success = run_benchmark(args.n_values, args.algorithms)
        sys.exit(0 if success else 1)

    # Validate that n is at least 6
    if args.n < 6:
        log_error(f"Input n must be at least 6, but got {args.n}")
        return
    
    start = time.perf_counter()

    if args.new_algo:
        gen_func = get_minimal_generators_new
        algo_name = 'gauss_jordan_solve'
    else:
        gen_func = get_minimal_generators
        algo_name = 'solve'
    
    result = gen_func(args.n, singular=args.singular, emit=args.show_polys, progress=True, debug=args.debug)
    end = time.perf_counter()

    #--------- Write output file ---------#
    import os
    out_dir = 'output_min_generators'
    os.makedirs(out_dir, exist_ok=True)
    suffix = '_NEW' if args.new_algo else ''
    filename = os.path.join(out_dir, f'minimal_generators_n{args.n}{suffix}.txt')
    with open(filename, 'w', encoding='utf-8') as fh:
        for idx, poly in enumerate(result, start=1):
            line = poly.write(singular=args.singular)
            fh.write(f'f_{idx} = {line.strip()}\n')
    log_debug(f'Wrote {len(result)} polynomials to {filename}', enabled=args.debug)
    #------------------------------------#

    if start is not None and end is not None:
        elapsed = end - start
        log_success(f"Execution completed in {elapsed:.4f} seconds")
        
        # Record timing
        record_timing(args.n, elapsed, algo_name)
        log_debug(f"Timing recorded: n={args.n}, algorithm={algo_name}, time={elapsed:.4f}s", enabled=args.debug)

if __name__=='__main__':
    main()
