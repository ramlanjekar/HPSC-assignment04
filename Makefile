# Makefile for BICGSTAB Project
# This automates compilation of your project

# Compiler settings
CC = /opt/homebrew/bin/gcc-15                  # Use your alias for real GCC
CFLAGS = -O3 -Wall -Wextra   # Optimization level 3, enable warnings
LDFLAGS = -lm                 # Link math library

# OpenMP flag
OMPFLAG = -fopenmp

# Source files
FEM_SRC = fem_matrix.c
SERIAL_SRC = bicgstab_serial.c
PARALLEL_SRC = bicgstab_parallel.c
MAIN_SRC = main.c

# Object files
FEM_OBJ = fem_matrix.o
SERIAL_OBJ = bicgstab_serial.o
PARALLEL_OBJ = bicgstab_parallel.o
MAIN_OBJ = main.o

# Executable name
TARGET = bicgstab_solver

# Default target: build everything
all: $(TARGET)
	@echo ""
	@echo "================================================"
	@echo "Build complete! Run with: ./$(TARGET)"
	@echo "================================================"

# Link all object files into final executable
$(TARGET): $(FEM_OBJ) $(SERIAL_OBJ) $(PARALLEL_OBJ) $(MAIN_OBJ)
	$(CC) $(CFLAGS) $(OMPFLAG) -o $@ $^ $(LDFLAGS)

# Compile FEM matrix generation (no OpenMP needed)
$(FEM_OBJ): $(FEM_SRC) fem_matrix.h
	$(CC) $(CFLAGS) -c $(FEM_SRC)

# Compile serial solver (no OpenMP needed)
$(SERIAL_OBJ): $(SERIAL_SRC) fem_matrix.h
	$(CC) $(CFLAGS) -c $(SERIAL_SRC)

# Compile parallel solver (needs OpenMP)
$(PARALLEL_OBJ): $(PARALLEL_SRC) fem_matrix.h
	$(CC) $(CFLAGS) $(OMPFLAG) -c $(PARALLEL_SRC)

# Compile main program (needs OpenMP for linking)
$(MAIN_OBJ): $(MAIN_SRC) fem_matrix.h
	$(CC) $(CFLAGS) $(OMPFLAG) -c $(MAIN_SRC)

# Clean up compiled files
clean:
	rm -f *.o $(TARGET)
	@echo "Cleaned up all build files"

# Rebuild from scratch
rebuild: clean all

# Run the program after building
run: $(TARGET)
	./$(TARGET)

# Help message
help:
	@echo "Available targets:"
	@echo "  make          - Build the project"
	@echo "  make clean    - Remove compiled files"
	@echo "  make rebuild  - Clean and rebuild"
	@echo "  make run      - Build and run the program"
	@echo "  make help     - Show this message"

.PHONY: all clean rebuild run help