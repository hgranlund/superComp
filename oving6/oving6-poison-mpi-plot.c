#include "common.h"

void runPoisson(int rank, int size, int n){
  double time=MPI_Wtime();
  Real **b, *diag, *RecvBuf,*z, h, umax;
  int i, j, m, nn, *len, *disp;

  m  = n-1;
  nn = 4*n;
  splitVector(m, size, &len, &disp);
  diag = createRealArray (m);
  b    = createReal2DArray (len[rank],m);
  z    = createRealArray (nn);
  h    = 1./(Real)n;

  #pragma omp parallel for schedule(static)
  for (i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*M_PI/(Real)n));
  }

  #pragma omp for
  for (j=0; j < len[rank]; j++) {
  #pragma omp parallel for schedule(static)
    for (i=0; i < m; i++) {
      b[j][i] = h*h;
    }
  }

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    Real* zt = createRealArray (nn);
    fst_(b[j], &n, zt, &nn);
    free(zt);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (i=0; i < len[rank]; i++) {
    Real* zt  = createRealArray (nn);
    fstinv_(b[i], &n, zt, &nn);
    free(zt);
  }

  #pragma omp for
  for (j=0; j < len[rank]; j++) {
  #pragma omp parallel for schedule(static)
    for (i=0; i < m; i++) {
      b[j][i] = b[j][i]/(diag[i]+diag[j+disp[rank]]);
    }
  }

  #pragma omp parallel for schedule(static)
  for (i=0; i < len[rank]; i++) {
    Real* zt  = createRealArray (nn);
    fst_(b[i], &n, zt, &nn);
    free(zt);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    Real* zt  = createRealArray (nn);
    fstinv_(b[j], &n, zt, &nn);
    free(zt);
  }


  if (rank==0)
  {
    RecvBuf = createRealArray (m*m);
  }
  gatherMatrix(b, m, RecvBuf, len, disp,0);

  if (rank==0)
  {
    for (int j=0; j < m; j++) {
      for (int i=0; i < m; i++) {
        printf("%e %e %e \n",(Real)i/m,(Real)j/m,RecvBuf[j*m+i] );
      }
    }
  }
}

int main(int argc, char **argv )
{
  int rank, size, n, nrange;
  init_app (argc, argv, &rank, &size);
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
  for (int i = n; i <= nrange ; ++i){
    runPoisson(rank, size, pow(2,i));
  }
}

MPI_Finalize();
return 0;
}

