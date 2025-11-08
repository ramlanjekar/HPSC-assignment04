// main.c
// Main program to test serial and parallel BICGSTAB solvers
// Runs benchmarks for different grid sizes and thread counts

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fem_matrix.h"

// External solver functions
extern int bicgstab_serial(FEMSystem *sys, int max_iter, double tol, double *solve_time);
extern int bicgstab_parallel(FEMSystem *sys, int max_iter, double tol, int num_threads, double *solve_time);

// Function to verify solution quality
void verify_solution(FEMSystem *sys) {
    int n = sys->n;
    double *Ax = (double*)malloc(n * sizeof(double));
    
    // Compute Ax
    matvec_csr(&sys->A, sys->x, Ax);
    
    // Compute residual norm: ||b - Ax||
    double residual = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = sys->b[i] - Ax[i];
        residual += diff * diff;
    }
    residual = sqrt(residual);
    
    printf("Final residual norm: %.6e\n", residual);
    
    free(Ax);
}

// Run benchmark for a given grid size
void run_benchmark(int nx, int ny) {
    printf("\n");
    printf("========================================\n");
    printf("Grid Size: %d x %d = %d nodes\n", nx, ny, nx*ny);
    printf("========================================\n");
    
    // Create FEM system
    FEMSystem *sys = create_fem_system(nx, ny);
    print_system_info(sys);
    
    int max_iter = 10000;
    double tol = 1e-8;
    
    // Serial solve
    printf("\n--- Serial BICGSTAB ---\n");
    // Reset solution to zero
    memset(sys->x, 0, sys->n * sizeof(double));
    
    double serial_time;
    int iter_serial = bicgstab_serial(sys, max_iter, tol, &serial_time);
    
    if (iter_serial > 0) {
        printf("Time: %.6f seconds\n", serial_time);
        verify_solution(sys);
    }
    
    // Parallel solves with different thread counts
    int thread_counts[] = {2, 4, 8};
    int num_configs = 3;
    
    printf("\n--- Parallel BICGSTAB ---\n");
    printf("%-10s %-15s %-15s %-10s\n", "Threads", "Time (s)", "Speedup", "Efficiency");
    printf("------------------------------------------------------\n");
    
    for (int i = 0; i < num_configs; i++) {
        int num_threads = thread_counts[i];
        
        // Reset solution to zero
        memset(sys->x, 0, sys->n * sizeof(double));
        
        double parallel_time;
        int iter_parallel = bicgstab_parallel(sys, max_iter, tol, num_threads, &parallel_time);
        
        if (iter_parallel > 0) {
            double speedup = serial_time / parallel_time;
            double efficiency = speedup / num_threads * 100.0;
            
            printf("%-10d %-15.6f %-15.2f %-10.1f%%\n", 
                   num_threads, parallel_time, speedup, efficiency);
        }
    }
    
    // Free system
    free_fem_system(sys);
}

int main() {
    printf("===================================================\n");
    printf("     OpenMP Parallelized BICGSTAB Solver\n");
    printf("     2D Laplace Equation on Unit Square\n");
    printf("===================================================\n");
    printf("Boundary Conditions:\n");
    printf("  Bottom, Left, Right: T = 0\n");
    printf("  Top: T = 1\n");
    printf("===================================================\n");
    
    // Test different grid sizes
    // Grid sizes chosen to give approximately 100, 200, and 400 nodes
    
    printf("\n\n***** TEST 1: Small Grid (~100 nodes) *****\n");
    run_benchmark(10, 10);  // 100 nodes
    
    printf("\n\n***** TEST 2: Medium Grid (~200 nodes) *****\n");
    run_benchmark(14, 14);  // 196 nodes
    
    printf("\n\n***** TEST 3: Large Grid (~400 nodes) *****\n");
    run_benchmark(20, 20);  // 400 nodes
    
    printf("\n\n===================================================\n");
    printf("               Benchmark Complete\n");
    printf("===================================================\n");
    
    return 0;
}