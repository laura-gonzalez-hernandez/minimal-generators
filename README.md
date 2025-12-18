# minimal-generators

(Here goes the explanation of what it does and reference to the article)

## Requirements

- **Python**: ≥ 3.11
- **Dependencies**: SymPy ≥ 1.12, NumPy ≥ 1.23

## Installation

### Option 1: Using `pip` (Standard Installation)

```bash
# Install directly from GitHub
pip install git+https://github.com/oriol-alm/minimal-generators-v0.git

# Run the installed command
minimal-generators 10
```

or alternatively, 
1. Download the repository
2. Open a terminal and navigate to the downloader folder
3. Run ```pip install .```

### Option 2: Using `uv` (Recommended for Development)

This project uses [uv](https://github.com/astral-sh/uv) for fast dependency management.

#### Install `uv`

**macOS** (via Homebrew):
```bash
brew install uv
```

**Linux/Ubuntu** (via curl):
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**macOS** (via curl):
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows** (via PowerShell):
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

#### Clone and Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/minimal-generators.git
cd minimal-generators

# Sync dependencies
uv sync

# Run the tool
uv run minimal-generators 10
```

### Option 3: Using `pipx` (Isolated CLI Tool)

```bash
# Install as an isolated command-line tool
pipx install git+https://github.com/yourusername/minimal-generators.git

# Run from anywhere
minimal-generators 10
```

### Option 4: Development Mode

```bash
# Clone and install in editable mode
git clone https://github.com/yourusername/minimal-generators.git
cd minimal-generators
pip install -e .

# Run the command
minimal-generators 10
```

## Usage

### Basic Command

```bash
minimal-generators <n> --show-polys --singular
```

or with `uv`:

```bash
uv run minimal-generators <n>
```

where `n` is an integer ≥ 5.

### Command-Line Options

| Option | Description |
|--------|-------------|
| `n` | **Required.** Positive integer (minimum: 5) |
| `--singular` | Output polynomials in [Singular](https://www.singular.uni-kl.de/) format for direct use in the Singular computer algebra system (default: False) |
| `--show-polys` | Print polynomials to terminal in addition to saving to file (default: False) |
| `--debug` | Show debug messages during computation (default: False) |
| `--new-algo` | Use optimized Gauss-Jordan solver (~2x faster) (default: False) |

### Subcommands

#### Generate Polynomials

```bash
minimal-generators generate <n> [options]
```

This is the default command when you run `minimal-generators <n>`.

#### Compare Output Files

```bash
minimal-generators compare <n>
```

Compares the standard and optimized algorithm outputs for a given `n` value.

#### Benchmark Performance

```bash
minimal-generators benchmark <n1> <n2> <n3> ... [options]
```

Run benchmarks for multiple n values and automatically record timing data.

**Options:**
- `--algorithms {solve,gauss_jordan_solve}`: Specify which algorithms to benchmark (default: both)

**Examples:**
```bash
# Benchmark both algorithms for n=10,20,40,60
minimal-generators benchmark 10 20 40 60

# Benchmark only the optimized algorithm
minimal-generators benchmark 10 20 40 --algorithms gauss_jordan_solve
```

#### Generate Performance Report

```bash
minimal-generators report [options]
```

Generate a performance comparison plot from stored timing data. Every time you run `minimal-generators` with a specific `n` value, the execution time is automatically recorded. Use this command to visualize all recorded timings.

**Options:**
- `--output <filename>`: Specify output filename for the plot (default: `benchmark_results.png`)
- `--clear`: Clear all stored timing data

**Examples:**
```bash
# Generate report with default filename
minimal-generators report

# Generate report with custom filename
minimal-generators report --output my_results.png

# Clear all stored timing data
minimal-generators report --clear
```

### Examples

**Basic usage (generates and saves to file):**
```bash
minimal-generators 10
```
Output: `output_min_generators/minimal_generators_n10.txt`

**Show polynomials in terminal:**
```bash
minimal-generators 18 --show-polys
```

**Use optimized algorithm with Singular format:**
```bash
minimal-generators 50 --new-algo --singular
```

**Compare algorithm outputs:**
```bash
# First generate with both algorithms
minimal-generators 20
minimal-generators 20 --new-algo

# Then compare
minimal-generators compare 20
```

**Run benchmarks and generate report:**
```bash
# Run benchmarks for multiple n values
minimal-generators benchmark 10 20 40 60

# Generate performance plot
minimal-generators report
```

**Automatic timing tracking:**
```bash
# Each run automatically records timing
minimal-generators 10          # Records timing for 'solve' algorithm
minimal-generators 10 --new-algo  # Records timing for 'gauss_jordan_solve'
minimal-generators 20
minimal-generators 20 --new-algo

# View all recorded timings
minimal-generators report
```

**Debug mode:**
```bash
minimal-generators 25 --debug
```

### Output Files

All polynomial files are automatically saved to:
```
output_min_generators/minimal_generators_n<n>.txt       # Standard algorithm
output_min_generators/minimal_generators_n<n>_NEW.txt   # Optimized algorithm (--new-algo)
```

**Example output format (with `--singular` flag):**
```
f_1 = x^5 + 2*y^3*z^2;
f_2 = -3*x^2*y + z^5;
...
```

### Timing Data Storage

Execution times are automatically stored in `~/.minimal_generators_timing.json` for performance tracking. This allows you to:
- Track performance improvements over time
- Compare different algorithm versions
- Generate performance reports without re-running computations

The timing data includes:
- The `n` value used
- Algorithm type (`solve` or `gauss_jordan_solve`)
- Execution time in seconds
- Timestamp of the run

## Performance Notes

To-do

## License

To-do

## Contributing

To-do

## Authors

- To-do
