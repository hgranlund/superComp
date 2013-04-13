#pragma once
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#ifndef M_PI
#define M_PI (4.0*atan(1))
#endif

typedef double Real;


Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
Real maxPointwiseError(Real** Matrix, int matrixSize, int* len, int* disp, int root);
void transpose (Real **b, int size, int *len, int *disp, int rank, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void init_app(int argc, char** argv, int* rank, int* size);
void splitVector(int globLen, int size, int** len, int** displ);
void matrixToVector(Real **matrix, Real *vector, int *len, int *disp, int size, int rank );
void vectorToMatrix(Real **matrix, Real *vector, int *len, int *disp, int size, int rank );
void gatherMatrix(Real** matrix, int matrixSize, Real* gatherRecvBuf, int* len, int* disp, int root);
