#include "common.h"

int main(int argc, char **argv )
{
  double time=MPI_Wtime();
  Real **b, *diag, *gatherRecvBuf,*z, h, umax;
  int i, j, n, m, nn, rank, size , *len, *disp;

  if( argc < 2 ) {
    printf("need a problem size\n");
    return 0;
  }

  init_app (argc, argv, &rank, &size);
  n  = atoi(argv[1]);
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

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = h*h;
    }
  }

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    fst_(b[j], &n, z, &nn);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (i=0; i < len[rank]; i++) {
    fstinv_(b[i], &n, z, &nn);
  }

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = b[j][i]/(diag[i]+diag[j+disp[rank]]);
    }
  }

  #pragma omp parallel for schedule(static)
  for (i=0; i < len[rank]; i++) {
    fst_(b[i], &n, z, &nn);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (j=0; j < len[rank]; j++) {
    fstinv_(b[j], &n, z, &nn);
  }

  if (rank==0)
  {
    gatherRecvBuf = createRealArray (m*m);
  }

  gatherMatrix(b, m, gatherRecvBuf, len, disp,0);


  if (rank==0)
  {
    umax = 0.0;
    for (i=0; i < m*m; i++) {
      if (gatherRecvBuf[i] > umax) umax = gatherRecvBuf[i];
    }
    printf("elapsed: %f\n", MPI_Wtime()-time);
    printf (" umax = %e \n",umax);
  }

  MPI_Finalize();
  return 0;
}
