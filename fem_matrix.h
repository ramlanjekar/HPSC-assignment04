// fem_matrix.h
// Header file for FEM matrix generation
// This defines the structures and functions we'll use

#ifndef FEM_MATRIX_H
#define FEM_MATRIX_H

// Structure to hold our sparse matrix in CSR format
// CSR = Compressed Sparse Row (efficient for sparse matrices)
typedef struct {
    int n;              // Number of nodes (equations)
    int nnz;            // Number of non-zero entries
    double *values;     // Non-zero values
    int *col_idx;       // Column indices for each non-zero
    int *row_ptr;       // Pointers to start of each row
} CSRMatrix;

// Structure to hold the linear system Ax=b
typedef struct {
    CSRMatrix A;        // The stiffness matrix
    double *b;          // Right-hand side vector
    double *x;          // Solution vector (initialized to zeros)
    int n;              // Problem size
} FEMSystem;

// Function declarations
// Creates the FEM system for given grid size (nx x ny)
FEMSystem* create_fem_system(int nx, int ny);

// Frees all allocated memory
void free_fem_system(FEMSystem *sys);

// Prints system info (for debugging)
void print_system_info(FEMSystem *sys);

// Matrix-vector multiplication: y = A*x
void matvec_csr(CSRMatrix *A, double *x, double *y);

#endif // FEM_MATRIX_H