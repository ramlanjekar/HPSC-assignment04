// bicgstab_serial.c
// Serial implementation of BICGSTAB iterative solver
// BICGSTAB = Biconjugate Gradient Stabilized method

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fem_matrix.h"

// Vector operations
static double dot_product(double *a, double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

static void vector_copy(double *src, double *dst, int n) {
    for (int i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}

static void vector_axpy(double a, double *x, double *y, int n) {
    // y = a*x + y
    for (int i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

static void vector_axpby(double a, double *x, double b, double *y, double *z, int n) {
    // z = a*x + b*y
    for (int i = 0; i < n; i++) {
        z[i] = a * x[i] + b * y[i];
    }
}

static double vector_norm(double *x, int n) {
    return sqrt(dot_product(x, x, n));
}

// BICGSTAB solver
// Solves Ax = b using BICGSTAB method
// Returns: number of iterations, or -1 if failed
int bicgstab_serial(FEMSystem *sys, int max_iter, double tol, double *solve_time) {
    int n = sys->n;
    CSRMatrix *A = &sys->A;
    double *b = sys->b;
    double *x = sys->x;
    
    // Allocate working vectors
    double *r = (double*)malloc(n * sizeof(double));      // residual
    double *r0 = (double*)malloc(n * sizeof(double));     // shadow residual
    double *p = (double*)malloc(n * sizeof(double));      // search direction
    double *v = (double*)malloc(n * sizeof(double));      // A*p
    double *s = (double*)malloc(n * sizeof(double));      // intermediate residual
    double *t = (double*)malloc(n * sizeof(double));      // A*s
    
    // Start timing
    clock_t start = clock();
    
    // Initial guess: x = 0 (already initialized)
    // Initial residual: r = b - A*x = b
    vector_copy(b, r, n);
    
    // CRITICAL FIX: Use constant r0 instead of r0 = r
    // This prevents breakdown when boundary conditions cause
    // the residual's non-zero pattern to shift
    for (int i = 0; i < n; i++) {
        r0[i] = 1.0;  // Constant vector with support everywhere
    }
    
    vector_copy(r, p, n);
    
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    double rho_prev, beta;
    
    double bnorm = vector_norm(b, n);
    if (bnorm == 0.0) bnorm = 1.0;
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        rho_prev = rho;
        rho = dot_product(r0, r, n);
        
        if (fabs(rho) < 1e-30) {
            printf("BICGSTAB: rho breakdown at iteration %d\n", iter);
            break;
        }
        
        if (iter == 0) {
            vector_copy(r, p, n);
        } else {
            beta = (rho / rho_prev) * (alpha / omega);
            // p = r + beta*(p - omega*v)
            for (int i = 0; i < n; i++) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }
        
        // v = A*p
        matvec_csr(A, p, v);
        
        alpha = rho / dot_product(r0, v, n);
        
        // s = r - alpha*v
        vector_axpby(1.0, r, -alpha, v, s, n);
        
        // Check convergence
        double s_norm = vector_norm(s, n);
        if (s_norm / bnorm < tol) {
            // x = x + alpha*p
            vector_axpy(alpha, p, x, n);
            printf("BICGSTAB converged at iteration %d (residual: %.2e)\n", 
                   iter+1, s_norm/bnorm);
            break;
        }
        
        // t = A*s
        matvec_csr(A, s, t);
        
        omega = dot_product(t, s, n) / dot_product(t, t, n);
        
        // x = x + alpha*p + omega*s
        vector_axpy(alpha, p, x, n);
        vector_axpy(omega, s, x, n);
        
        // r = s - omega*t
        vector_axpby(1.0, s, -omega, t, r, n);
        
        // Check convergence
        double r_norm = vector_norm(r, n);
        if (r_norm / bnorm < tol) {
            printf("BICGSTAB converged at iteration %d (residual: %.2e)\n", 
                   iter+1, r_norm/bnorm);
            iter++;
            break;
        }
        
        if (fabs(omega) < 1e-30) {
            printf("BICGSTAB: omega breakdown at iteration %d\n", iter);
            break;
        }
    }
    
    // End timing
    clock_t end = clock();
    *solve_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    // Free working vectors
    free(r);
    free(r0);
    free(p);
    free(v);
    free(s);
    free(t);
    
    if (iter >= max_iter) {
        printf("BICGSTAB did not converge within %d iterations\n", max_iter);
        return -1;
    }
    
    return iter;
}