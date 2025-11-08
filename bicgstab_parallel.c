// bicgstab_parallel.c
// OpenMP parallelized version of BICGSTAB solver
// Parallelizes vector operations and matrix-vector multiplication

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "fem_matrix.h"

// Parallel vector dot product
static double dot_product_parallel(double *a, double *b, int n) {
    double sum = 0.0;
    // OpenMP reduction: each thread computes partial sum, then combines
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

// Parallel vector copy
static void vector_copy_parallel(double *src, double *dst, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}

// Parallel y = a*x + y
static void vector_axpy_parallel(double a, double *x, double *y, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

// Parallel z = a*x + b*y
static void vector_axpby_parallel(double a, double *x, double b, double *y, double *z, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        z[i] = a * x[i] + b * y[i];
    }
}

// Parallel vector norm
static double vector_norm_parallel(double *x, int n) {
    return sqrt(dot_product_parallel(x, x, n));
}

// Parallel matrix-vector multiplication for CSR format
static void matvec_csr_parallel(CSRMatrix *A, double *x, double *y) {
    int n = A->n;
    // Each row is independent, so we can parallelize over rows
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++) {
            sum += A->values[j] * x[A->col_idx[j]];
        }
        y[i] = sum;
    }
}

// Parallel BICGSTAB solver
// num_threads: number of OpenMP threads to use
int bicgstab_parallel(FEMSystem *sys, int max_iter, double tol, int num_threads, double *solve_time) {
    int n = sys->n;
    CSRMatrix *A = &sys->A;
    double *b = sys->b;
    double *x = sys->x;
    
    // Set number of threads
    omp_set_num_threads(num_threads);
    
    // Allocate working vectors
    double *r = (double*)malloc(n * sizeof(double));
    double *r0 = (double*)malloc(n * sizeof(double));
    double *p = (double*)malloc(n * sizeof(double));
    double *v = (double*)malloc(n * sizeof(double));
    double *s = (double*)malloc(n * sizeof(double));
    double *t = (double*)malloc(n * sizeof(double));
    
    // Start timing (use omp_get_wtime for better precision)
    double start = omp_get_wtime();
    
    // Initial guess: x = 0
    // Initial residual: r = b - A*x = b
    vector_copy_parallel(b, r, n);
    
    // FIX: Use a constant r0 to avoid breakdown when boundary conditions
    // cause r to change its support (non-zero pattern)
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        r0[i] = 1.0;  // Uniform vector has support everywhere
    }
    
    vector_copy_parallel(r, p, n);
    
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    double rho_prev, beta;
    
    double bnorm = vector_norm_parallel(b, n);
    if (bnorm == 0.0) bnorm = 1.0;
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        rho_prev = rho;
        rho = dot_product_parallel(r0, r, n);
        
        if (fabs(rho) < 1e-30) {
            printf("BICGSTAB (parallel): rho breakdown at iteration %d\n", iter);
            break;
        }
        
        if (iter == 0) {
            vector_copy_parallel(r, p, n);
        } else {
            beta = (rho / rho_prev) * (alpha / omega);
            // p = r + beta*(p - omega*v)
            #pragma omp parallel for
            for (int i = 0; i < n; i++) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }
        
        // v = A*p
        matvec_csr_parallel(A, p, v);
        
        alpha = rho / dot_product_parallel(r0, v, n);
        
        // s = r - alpha*v
        vector_axpby_parallel(1.0, r, -alpha, v, s, n);
        
        // Check convergence
        double s_norm = vector_norm_parallel(s, n);
        if (s_norm / bnorm < tol) {
            vector_axpy_parallel(alpha, p, x, n);
            printf("BICGSTAB (parallel, %d threads) converged at iteration %d (residual: %.2e)\n", 
                   num_threads, iter+1, s_norm/bnorm);
            break;
        }
        
        // t = A*s
        matvec_csr_parallel(A, s, t);
        
        omega = dot_product_parallel(t, s, n) / dot_product_parallel(t, t, n);
        
        // x = x + alpha*p + omega*s
        vector_axpy_parallel(alpha, p, x, n);
        vector_axpy_parallel(omega, s, x, n);
        
        // r = s - omega*t
        vector_axpby_parallel(1.0, s, -omega, t, r, n);
        
        // Check convergence
        double r_norm = vector_norm_parallel(r, n);
        if (r_norm / bnorm < tol) {
            printf("BICGSTAB (parallel, %d threads) converged at iteration %d (residual: %.2e)\n", 
                   num_threads, iter+1, r_norm/bnorm);
            iter++;
            break;
        }
        
        if (fabs(omega) < 1e-30) {
            printf("BICGSTAB (parallel): omega breakdown at iteration %d\n", iter);
            break;
        }
    }
    
    // End timing
    double end = omp_get_wtime();
    *solve_time = end - start;
    
    // Free working vectors
    free(r);
    free(r0);
    free(p);
    free(v);
    free(s);
    free(t);
    
    if (iter >= max_iter) {
        printf("BICGSTAB (parallel) did not converge within %d iterations\n", max_iter);
        return -1;
    }
    
    return iter;
}