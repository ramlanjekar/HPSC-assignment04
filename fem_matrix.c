// fem_matrix.c
// Generates the FEM matrix for 2D Laplace equation on unit square
// Using 4-node rectangular elements (bilinear basis functions)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fem_matrix.h"

// Creates node numbering for structured grid
// Returns global node number for grid position (i,j)
static int get_node_number(int i, int j, int nx) {
    return i * nx + j;
}

// Checks if a node is on boundary
// Returns: 0=interior, 1=bottom, 2=right, 3=top, 4=left
static int get_boundary_type(int i, int j, int nx, int ny) {
    if (j == 0) return 1;      // Bottom: T=0
    if (i == ny-1) return 2;   // Right: T=0
    if (j == nx-1) return 3;   // Top: T=1
    if (i == 0) return 4;      // Left: T=0
    return 0;                  // Interior
}

// Main function to create the FEM system
FEMSystem* create_fem_system(int nx, int ny) {
    printf("Creating FEM system: %dx%d grid (%d nodes)\n", nx, ny, nx*ny);
    
    // Allocate system structure
    FEMSystem *sys = (FEMSystem*)malloc(sizeof(FEMSystem));
    sys->n = nx * ny;
    
    // For a structured grid, each interior node connects to ~5 nodes (itself + 4 neighbors)
    // Boundary nodes have fewer connections
    // Estimate: roughly 5*n non-zeros (we'll count exactly later)
    int max_nnz = 5 * sys->n;
    
    // Temporary storage for building matrix (we'll compress later)
    double *temp_values = (double*)calloc(max_nnz, sizeof(double));
    int *temp_col = (int*)calloc(max_nnz, sizeof(int));
    int *temp_row = (int*)calloc(sys->n + 1, sizeof(int));
    
    // Allocate solution vectors
    sys->b = (double*)calloc(sys->n, sizeof(double));
    sys->x = (double*)calloc(sys->n, sizeof(double));
    
    // Grid spacing
    double hx = 1.0 / (nx - 1);  // x-direction spacing
    double hy = 1.0 / (ny - 1);  // y-direction spacing
    
    // Element stiffness matrix entries for Laplace equation
    // For rectangular element with bilinear basis functions
    double ke = (hy/hx + hx/hy) / 3.0;  // Diagonal contribution
    double kn = -(hy/hx) / 6.0;         // North-south neighbor
    double kw = -(hx/hy) / 6.0;         // East-west neighbor
    
    int nnz_count = 0;
    temp_row[0] = 0;
    
    // Loop over all nodes to build matrix
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            int node = get_node_number(i, j, nx);
            int boundary = get_boundary_type(i, j, nx, ny);
            
            if (boundary > 0) {
                // Boundary node: apply Dirichlet condition
                // Row equation becomes: T[node] = boundary_value
                // Matrix: 1.0 on diagonal, 0 elsewhere
                temp_values[nnz_count] = 1.0;
                temp_col[nnz_count] = node;
                nnz_count++;
                
                // Right-hand side
                if (boundary == 3) {
                    sys->b[node] = 1.0;  // Top boundary: T=1
                } else {
                    sys->b[node] = 0.0;  // Other boundaries: T=0
                }
            } else {
                // Interior node: assemble stiffness matrix
                // 5-point stencil for 2D Laplace operator
                
                // West neighbor (j-1)
                if (j > 0) {
                    int neighbor = get_node_number(i, j-1, nx);
                    temp_values[nnz_count] = kw;
                    temp_col[nnz_count] = neighbor;
                    nnz_count++;
                }
                
                // South neighbor (i-1)
                if (i > 0) {
                    int neighbor = get_node_number(i-1, j, nx);
                    temp_values[nnz_count] = kn;
                    temp_col[nnz_count] = neighbor;
                    nnz_count++;
                }
                
                // Diagonal (self)
                temp_values[nnz_count] = ke;
                temp_col[nnz_count] = node;
                nnz_count++;
                
                // North neighbor (i+1)
                if (i < ny-1) {
                    int neighbor = get_node_number(i+1, j, nx);
                    temp_values[nnz_count] = kn;
                    temp_col[nnz_count] = neighbor;
                    nnz_count++;
                }
                
                // East neighbor (j+1)
                if (j < nx-1) {
                    int neighbor = get_node_number(i, j+1, nx);
                    temp_values[nnz_count] = kw;
                    temp_col[nnz_count] = neighbor;
                    nnz_count++;
                }
                
                // Right-hand side (source term = 0 for Laplace)
                sys->b[node] = 0.0;
            }
            
            temp_row[node + 1] = nnz_count;
        }
    }
    
    // Now copy to final CSR structure with exact size
    sys->A.n = sys->n;
    sys->A.nnz = nnz_count;
    sys->A.values = (double*)malloc(nnz_count * sizeof(double));
    sys->A.col_idx = (int*)malloc(nnz_count * sizeof(int));
    sys->A.row_ptr = (int*)malloc((sys->n + 1) * sizeof(int));
    
    for (int i = 0; i < nnz_count; i++) {
        sys->A.values[i] = temp_values[i];
        sys->A.col_idx[i] = temp_col[i];
    }
    
    for (int i = 0; i <= sys->n; i++) {
        sys->A.row_ptr[i] = temp_row[i];
    }
    
    // Free temporary storage
    free(temp_values);
    free(temp_col);
    free(temp_row);
    
    printf("Matrix created: %d nodes, %d non-zeros\n", sys->n, sys->A.nnz);
    return sys;
}

// Free all memory
void free_fem_system(FEMSystem *sys) {
    if (sys) {
        free(sys->A.values);
        free(sys->A.col_idx);
        free(sys->A.row_ptr);
        free(sys->b);
        free(sys->x);
        free(sys);
    }
}

// Print system information
void print_system_info(FEMSystem *sys) {
    printf("\n=== FEM System Info ===\n");
    printf("Number of nodes: %d\n", sys->n);
    printf("Number of non-zeros: %d\n", sys->A.nnz);
    printf("Sparsity: %.2f%%\n", 100.0 * sys->A.nnz / (sys->n * sys->n));
}

// Matrix-vector multiplication for CSR format: y = A*x
void matvec_csr(CSRMatrix *A, double *x, double *y) {
    for (int i = 0; i < A->n; i++) {
        y[i] = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++) {
            y[i] += A->values[j] * x[A->col_idx[j]];
        }
    }
}