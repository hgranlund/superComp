#include "common.h"

Real funcf(Real x, Real y)
{
  return (5*M_PI*M_PI*sin(M_PI*x)*sin(2*M_PI*y));
}
Real funcfSolution(Real x, Real y)
{
  return (sin(x*M_PI)*sin(2*y*M_PI));
}

Real maxPointwiseError(Real** Matrix, int matrixSize, int* len, int* disp, int root){
  int rank;
  Real maxError, localeMaxError, x, y, n;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  localeMaxError = 0.0;
  n=matrixSize+1;
  for (int j=0; j < len[rank]; j++) {
    for (int i=0; i < matrixSize; i++) {
      x=(Real)(j+1+disp[rank])/n;
      y=(Real) (i+1)/n;
      if (Matrix[j][i]-funcfSolution(x, y) > localeMaxError) localeMaxError = Matrix[j][i]-funcfSolution(x, y);
    }
  }
  MPI_Reduce(&localeMaxError, &maxError, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  return maxError;
}

void runPoisson(int rank, int size, int n){
  double time=MPI_Wtime();
  Real **b, *diag, *RecvBuf, **z, h, maxError, x, y;
  int i, j, m, nn, *len, *disp;

  m  = n-1;
  nn = 4*n;
  splitVector(m, size, &len, &disp);
  diag = createRealArray (m);
  z    = createReal2DArray (omp_get_max_threads(), nn);
  b    = createReal2DArray (len[rank],m);
  h    = 1./(Real)n;

  for (i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*M_PI/(Real)n));
  }

  for (j=0; j < len[rank]; j++) {
    for (i=0; i < m; i++) {
      x=(Real)(j+1+disp[rank])/n;
      y=(Real) (i+1)/n;
      b[j][i] = h*h * funcf(x,y);
    }
  }

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    fst_(b[j], &n, z[omp_get_thread_num()], &nn);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (i=0; i < len[rank]; i++) {
    fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
  }

  for (j=0; j < len[rank]; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = b[j][i]/(diag[i]+diag[j+disp[rank]]);
    }
  }

  #pragma omp parallel for schedule(static)
  for (i=0; i < len[rank]; i++) {
    fst_(b[i], &n, z[omp_get_thread_num()], &nn);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
  }

  maxError = maxPointwiseError(b, m, len, disp, 0);

  if (rank==0)
  {
    printf ("Maximum pointwise error = %e, with n = %d \n",maxError, n);
    printf("Elapsed time: %f\n", MPI_Wtime()-time);
    printf ("================================================== \n");
  }
}

int main(int argc, char **argv )
{
  int rank, size, n, nrange;
  init_app (argc, argv, &rank, &size);
  if (rank == 0){
    printf("Number of threads per node = %d \n", omp_get_max_threads());
    printf("Number of MPI proc = %d \n", size);
    printf ("================================================== \n");
  }
  if( argc < 2 ) {
    if (rank == 0){
     printf("need a problem size\n");
   }
   MPI_Finalize();
   return 0;
 }

 if( argc == 2 ) {
  n  = atoi(argv[1]);
  runPoisson(rank, size, n);
}

else{
  n  = atoi(argv[1]);
  nrange  = atoi(argv[2]);
  if (rank==0){
    printf("Running with gridpionts n=2^k, with k ∈ [%d,%d]\n", n, nrange);
    printf ("================================================== \n");
  }
  for (int i = n; i <= nrange ; ++i){
    runPoisson(rank, size, pow(2,i));
  }
}

MPI_Finalize();
return 0;
}

