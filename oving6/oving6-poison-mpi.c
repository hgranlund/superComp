#include "common.h"

int main(int argc, char **argv )
{
  double time=MPI_Wtime();
  Real *diag, **b, *gatherSendBuf, *gatherRecvBuf, *z;
  Real pi, h, umax;
  int i, j, n, m, nn, rank, size , *len, *disp, ml, il;

  if( argc < 2 ) {
    printf("need a problem size\n");
    return 0;
  }
  init_app (argc, argv, &rank, &size);
  n  = atoi(argv[1]);
  m  = n-1;
  nn = 4*n;
  ml = len[rank];
  il = disp[rank];
  diag = createRealArray (m);
  b    = createReal2DArray (ml,m);
  z    = createRealArray (nn);
  h    = 1./(Real)n;
  pi   = 4.*atan(1.);
  splitVector(m, size, &len, &disp);

  for (i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }
  for (j=0; j < ml; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = h*h;
    }
  }

  #pragma omp parallel for schedule(static)
  for (j=0; j < ml; j++) {
    fst_(b[j], &n, z, &nn);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (i=0; i < ml; i++) {
    fstinv_(b[i], &n, z, &nn);
  }

  for (j=0; j < ml; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = b[j][i]/(diag[i]+diag[j+il]);
    }
  }

  #pragma omp parallel for schedule(static)
  for (i=0; i < ml; i++) {
    fst_(b[i], &n, z, &nn);
  }

  transpose(b, size, len, disp, rank, m);

  #pragma omp parallel for schedule(static)
  for (j=0; j < ml; j++) {
    fstinv_(b[j], &n, z, &nn);
  }

  if (rank==0)
  {
    gatherRecvBuf = createRealArray (m*m);
  }

  gatherSendBuf = createRealArray (m*ml);
  matrixToVector(b, gatherSendBuf, len, disp, size, rank);
  int *sendcounts, *rdispls, index;
  sendcounts = calloc(size,sizeof(int));
  rdispls = calloc(size,sizeof(int));
  index=0;
  for (int i = 0; i < size; ++i)
  {
    sendcounts[i]= len[i]*m;
    rdispls[i]=index;
    index=index+sendcounts[i];
  }
  MPI_Gatherv(gatherSendBuf, m*ml, MPI_DOUBLE, gatherRecvBuf, sendcounts, rdispls, MPI_DOUBLE, 0,MPI_COMM_WORLD );

  umax = 0.0;

  if (rank==0)
  {
    for (i=0; i < m*m; i++) {
      if (gatherRecvBuf[i] > umax) umax = gatherRecvBuf[i];
    }
    printf("elapsed: %f\n", MPI_Wtime()-time);
    printf (" umax = %e \n",umax);
  }

  MPI_Finalize();
  return 0;
}
