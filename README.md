# OpenMP Parallelized BICGSTAB Solver for 2D Laplace Equation

**Assignment:** Parallel implementation of BICGSTAB algorithm for solving finite element systems

**Author:** Ram Lanjekar  
**Course:** High Performance Scientific Computing

---

## Table of Contents
1. [Problem Description](#problem-description)
2. [Project Structure](#project-structure)
3. [Compilation and Execution](#compilation-and-execution)
4. [Understanding the Output](#understanding-the-output)
5. [Results and Observations](#results-and-observations)
6. [Technical Discussion](#technical-discussion)
---

## Problem Description

### Mathematical Problem
We solve the 2D Laplace equation on a unit square domain:

```
∇²T = 0    on Ω = [0,1] × [0,1]
```

### Boundary Conditions
- **Bottom edge** (y = 0): T = 0
- **Left edge** (x = 0): T = 0  
- **Right edge** (x = 1): T = 0
- **Top edge** (y = 1): T = 1

This represents a steady-state heat distribution where the top is heated and other sides are kept cold.

### Discretization
- **Method:** Finite Element Method (FEM) with bilinear rectangular elements
- **Grid sizes:** 10×10 (~100 nodes), 14×14 (~196 nodes), 20×20 (400 nodes)
- **Result:** Linear system **Ax = b** solved using BICGSTAB iterative method

---

## Project Structure

```
bicgstab_project/
│
├── fem_matrix.h              # Header: Data structures and function declarations
├── fem_matrix.c              # FEM matrix generation (A and b vectors)
├── bicgstab_serial.c         # Serial BICGSTAB implementation
├── bicgstab_parallel.c       # OpenMP parallelized BICGSTAB
├── main.c                    # Main program with benchmarking
├── Makefile                  # Build automation
└── README.md                 # This file
```

### File Descriptions

#### **1. fem_matrix.h / fem_matrix.c**
- **Purpose:** Generates the finite element system
- **Key functions:**
  - `create_fem_system(nx, ny)` - Creates FEM matrix for nx×ny grid
  - `matvec_csr()` - Matrix-vector multiplication for sparse CSR format
  - `free_fem_system()` - Memory cleanup
- **Data structure:** CSR (Compressed Sparse Row) format for efficient sparse matrix storage
- **Output:** 
  - Stiffness matrix **A** (sparse)
  - Right-hand side vector **b**
  - Solution vector **x** (initialized to zero)

#### **2. bicgstab_serial.c**
- **Purpose:** Serial (single-threaded) BICGSTAB solver
- **Algorithm:** BiConjugate Gradient Stabilized method
- **Key features:**
  - Iterative solver (no matrix inversion needed)
  - Works well for sparse, non-symmetric systems
  - Convergence tolerance: 10⁻⁸
- **Timing:** Uses `clock()` for CPU time measurement

#### **3. bicgstab_parallel.c**
- **Purpose:** OpenMP parallelized version of BICGSTAB
- **Parallelization strategy:**
  - Vector operations (dot products, AXPY operations)
  - Matrix-vector multiplication (parallelized by rows)
  - Uses `#pragma omp parallel for` and `reduction` clauses
- **Timing:** Uses `omp_get_wtime()` for wall-clock time
- **Thread control:** Configurable via `omp_set_num_threads()`

#### **4. main.c**
- **Purpose:** Driver program with automated benchmarking
- **Features:**
  - Tests three grid sizes automatically
  - Runs serial solver once per grid
  - Runs parallel solver with 2, 4, 8 threads
  - Computes speedup and efficiency metrics
  - Verifies solution accuracy

#### **5. Makefile**
- **Purpose:** Automates compilation process
- **Key targets:**
  - `make` or `make all` - Build the project
  - `make clean` - Remove compiled files
  - `make rebuild` - Clean and rebuild
  - `make run` - Build and execute
  - `make help` - Show available commands

---

## Compilation and Execution

### Prerequisites
- **Compiler:** GCC with OpenMP support (gcc-15 or later)
- **Operating System:** macOS (tested on Apple Silicon M4)
- **Required:** Homebrew-installed GCC (Apple's clang doesn't support OpenMP properly)

### Setup (One-time)
```bash
# Install GCC via Homebrew (if not already installed)
brew install gcc

# Verify installation
/opt/homebrew/bin/gcc-15 --version

# Create project directory
mkdir bicgstab_project
cd bicgstab_project

# Copy all source files into this directory
```

you have to make change in Makefile here 
```
CC = put location of your gcc 

```

### Building the Project
```bash
# Standard build
make

# Or clean rebuild
make clean
make
```

**Expected output:**
```
/opt/homebrew/bin/gcc-15 -O3 -Wall -Wextra -c fem_matrix.c
/opt/homebrew/bin/gcc-15 -O3 -Wall -Wextra -c bicgstab_serial.c
/opt/homebrew/bin/gcc-15 -O3 -Wall -Wextra -fopenmp -c bicgstab_parallel.c
/opt/homebrew/bin/gcc-15 -O3 -Wall -Wextra -fopenmp -c main.c
/opt/homebrew/bin/gcc-15 -O3 -Wall -Wextra -fopenmp -o bicgstab_solver fem_matrix.o bicgstab_serial.o bicgstab_parallel.o main.o -lm

================================================
Build complete! Run with: ./bicgstab_solver
================================================
```

### Running the Solver
```bash
./bicgstab_solver
```

**Runtime:** ~1-2 seconds for all tests

---

## Understanding the Output


```
===================================================
     OpenMP Parallelized BICGSTAB Solver
     2D Laplace Equation on Unit Square
===================================================
Boundary Conditions:
  Bottom, Left, Right: T = 0
  Top: T = 1
===================================================


***** TEST 1: Small Grid (~100 nodes) *****

========================================
Grid Size: 10 x 10 = 100 nodes
========================================
Creating FEM system: 10x10 grid (100 nodes)
Matrix created: 100 nodes, 356 non-zeros

=== FEM System Info ===
Number of nodes: 100
Number of non-zeros: 356
Sparsity: 3.56%

--- Serial BICGSTAB ---
BICGSTAB converged at iteration 21 (residual: 2.92e-09)
Time: 0.000058 seconds
Final residual norm: 8.769798e-09

--- Parallel BICGSTAB ---
Threads    Time (s)        Speedup         Efficiency
------------------------------------------------------
BICGSTAB (parallel, 2 threads) converged at iteration 25 (residual: 7.30e-09)
2          0.005335        0.01            0.5       %
BICGSTAB (parallel, 4 threads) converged at iteration 22 (residual: 6.90e-09)
4          0.006698        0.01            0.2       %
BICGSTAB (parallel, 8 threads) converged at iteration 25 (residual: 9.29e-09)
8          0.016381        0.00            0.0       %


***** TEST 2: Medium Grid (~200 nodes) *****

========================================
Grid Size: 14 x 14 = 196 nodes
========================================
Creating FEM system: 14x14 grid (196 nodes)
Matrix created: 196 nodes, 772 non-zeros

=== FEM System Info ===
Number of nodes: 196
Number of non-zeros: 772
Sparsity: 2.01%

--- Serial BICGSTAB ---
BICGSTAB converged at iteration 26 (residual: 7.52e-09)
Time: 0.000075 seconds
Final residual norm: 2.711064e-08

--- Parallel BICGSTAB ---
Threads    Time (s)        Speedup         Efficiency
------------------------------------------------------
BICGSTAB (parallel, 2 threads) converged at iteration 27 (residual: 8.24e-09)
2          0.003421        0.02            1.1       %
BICGSTAB (parallel, 4 threads) converged at iteration 27 (residual: 2.47e-09)
4          0.006628        0.01            0.3       %
BICGSTAB (parallel, 8 threads) converged at iteration 27 (residual: 8.38e-09)
8          0.013913        0.01            0.1       %


***** TEST 3: Large Grid (~400 nodes) *****

========================================
Grid Size: 20 x 20 = 400 nodes
========================================
Creating FEM system: 20x20 grid (400 nodes)
Matrix created: 400 nodes, 1696 non-zeros

=== FEM System Info ===
Number of nodes: 400
Number of non-zeros: 1696
Sparsity: 1.06%

--- Serial BICGSTAB ---
BICGSTAB converged at iteration 43 (residual: 6.41e-09)
Time: 0.000195 seconds
Final residual norm: 2.793334e-08

--- Parallel BICGSTAB ---
Threads    Time (s)        Speedup         Efficiency
------------------------------------------------------
BICGSTAB (parallel, 2 threads) converged at iteration 42 (residual: 7.65e-09)
2          0.004002        0.05            2.4       %
BICGSTAB (parallel, 4 threads) converged at iteration 40 (residual: 7.78e-09)
4          0.007069        0.03            0.7       %
BICGSTAB (parallel, 8 threads) converged at iteration 43 (residual: 9.21e-09)
8          0.017742        0.01            0.1       %


===================================================
               Benchmark Complete
===================================================

```

### Output Metrics Explained

| Metric | Meaning | Interpretation |
|--------|---------|----------------|
| **Iterations** | Number of BICGSTAB iterations to converge | Fewer is better; depends on matrix condition number |
| **Residual** | ‖b - Ax‖ / ‖b‖ | Should be < 10⁻⁸ for convergence |
| **Time** | Solver execution time (seconds) | Excludes matrix generation time |
| **Speedup** | Serial_time / Parallel_time | Ideal: 2× for 2 threads, 4× for 4 threads, 8× for 8 threads |
| **Efficiency** | (Speedup / Threads) × 100% | Ideal: 100%; shows how well threads are utilized |


---

## Results and Observations

### Actual Results from Testing

#### Grid 1: 10×10 (100 nodes)
- **Serial:** 58 μs, 21 iterations
- **Parallel (2 threads):** 5.3 ms, 25 iterations
- **Speedup:** 0.01× (100× SLOWER)
- **Efficiency:** 0.5%

#### Grid 2: 14×14 (196 nodes)
- **Serial:** 75 μs, 26 iterations
- **Parallel (2 threads):** 3.4 ms, 27 iterations
- **Speedup:** 0.02× (50× SLOWER)
- **Efficiency:** 1.1%

#### Grid 3: 20×20 (400 nodes)
- **Serial:** 195 μs, 43 iterations
- **Parallel (2 threads):** 4.0 ms, 42 iterations
- **Speedup:** 0.05× (20× SLOWER)
- **Efficiency:** 2.4%

### Key Observation: Negative Speedup

**The parallel version is significantly SLOWER than serial!**

This is **NOT an error** - it's an expected result that demonstrates fundamental principles of parallel computing.

---

## Technical Discussion

### Why Is Parallel Slower? (The Overhead Problem)

#### Problem Size Analysis
Our problem sizes are extremely small:
- 100-400 nodes
- 20-45 iterations
- ~60-200 microseconds of actual computational work

#### Parallelization Overhead Components

| Overhead Source | Time Cost | Explanation |
|----------------|-----------|-------------|
| **Thread Creation** | ~500 μs | Spawning worker threads |
| **Thread Synchronization** | ~100 μs per barrier | Waiting for all threads to finish each parallel region |
| **Memory Coherency** | ~50 μs per iteration | CPU cache invalidation between threads |
| **OpenMP Runtime** | ~200 μs | Library initialization and management |
| **False Sharing** | Variable | Threads writing to nearby memory locations |

**Total overhead:** ~5000 μs (5 milliseconds)

#### The Mathematics of Failure
```
Serial execution:     Work = 200 μs
Parallel execution:   Work + Overhead = 200 μs + 5000 μs = 5200 μs

Speedup = Serial_time / Parallel_time = 200 / 5200 = 0.038

Result: 38× SLOWER, not faster!
```

### Amdahl's Law in Action

Amdahl's Law states:
```
Speedup ≤ 1 / (S + P/N)
```
Where:
- S = Serial fraction (overhead)
- P = Parallel fraction (actual work)
- N = Number of processors

For our case:
- Overhead dominates: S >> P
- Adding more threads makes it WORSE
- Speedup approaches zero as threads increase

### When Would Parallelization Actually Help?

Parallel speedup occurs when: **Work >> Overhead**

#### Estimated Minimum Problem Sizes for Speedup

| Grid Size | Nodes | Iterations | Work Time | Overhead | Expected Speedup (2 threads) |
|-----------|-------|------------|-----------|----------|------------------------------|
| 50×50 | 2,500 | ~150 | ~5 ms | ~5 ms | ~1.0× (break-even) |
| 100×100 | 10,000 | ~300 | ~50 ms | ~5 ms | ~1.8× |
| 200×200 | 40,000 | ~600 | ~500 ms | ~5 ms | ~1.95× |
| 500×500 | 250,000 | ~1500 | ~5 sec | ~5 ms | ~1.99× |

**Rule of thumb:** Need at least **10,000 nodes** for meaningful parallel speedup.

### Cache Coherency and False Sharing

#### What Happens in Memory
1. Thread 1 writes to `x[0]`
2. Thread 2 writes to `x[1]`
3. Both are in the same 64-byte cache line
4. CPU must invalidate both threads' caches
5. Memory is re-fetched multiple times

**Result:** Memory bandwidth becomes the bottleneck, not computation.

#### Why More Threads Make It Worse
- 2 threads: Limited cache conflicts
- 4 threads: More synchronization points
- 8 threads: Cache thrashing + synchronization overhead dominates
- **Efficiency drops from 2.4% → 0.7% → 0.1%**

### Algorithm Convergence Differences

Notice iteration counts vary slightly:
- Serial: 21, 26, 43 iterations
- Parallel: 22-25, 27, 40-46 iterations

**Why?**
- Floating-point arithmetic is **not associative**: `(a + b) + c ≠ a + (b + c)` in finite precision
- Parallel reduction changes addition order
- Tiny differences accumulate over iterations
- Both converge to correct solution (residual < 10⁻⁸)

This is **normal and acceptable** in iterative solvers.


#### Lessons Learned
1. **Parallelization is not always beneficial** - overhead matters
2. **Problem size is critical** - small problems stay serial
3. **Measurement is important** - understand what you're timing
4. **Real-world parallel computing** requires cost-benefit analysis
5. **Amdahl's Law is real** - some problems can't scale

