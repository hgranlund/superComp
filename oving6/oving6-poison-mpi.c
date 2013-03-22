#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>

  typedef double Real;

/* function prototypes */
  Real *createRealArray (int n);
  Real **createReal2DArray (int m, int n);
  void transpose (Real **b, int size, int *len, int *disp, int rank, int m);
  void fst_(Real *v, int *n, Real *w, int *nn);
  void fstinv_(Real *v, int *n, Real *w, int *nn);
  void init_app(int argc, char** argv, int* rank, int* size);
  void splitVector(int globLen, int size, int** len, int** displ);




  int main(int argc, char **argv )
  {
    Real *diag, **b, **bt, *z;
    Real pi, h, umax;
    int i, j, n, m, nn, rank, size , *len, *disp, ml, il;


  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

    if( argc < 2 ) {
      printf("need a problem size\n");
      return 0;
    }
    init_app (argc, argv, &rank, &size);

    n  = atoi(argv[1]);
    m  = n-1;
    nn = 4*n;

    splitVector(m, size, &len, &disp);
    ml = len[rank];
    il = disp[rank];


    diag = createRealArray (m);
    b    = createReal2DArray (ml,m);
    z    = createRealArray (nn);


    h    = 1./(Real)n;
    pi   = 4.*atan(1.);



    for (i=0; i < m; i++) {
      diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }
    for (j=0; j < ml; j++) {
      for (i=0; i < m; i++) {
        b[j][i] = h*h;
      }
    }
    for (j=0; j < ml; j++) {
      fst_(b[j], &n, z, &nn);
    }


    transpose(b, size, len, disp, rank, m);

    for (i=0; i < ml; i++) {
      fstinv_(b[i], &n, z, &nn);
    }

    for (j=0; j < ml; j++) {
      for (i=0; i < m; i++) {
        b[j][i] = b[j][i]/(diag[i]+diag[j+il]);
      }
    }

    for (i=0; i < ml; i++) {
      fst_(b[i], &n, z, &nn);
    }

    transpose(b, size, len, disp, rank, m);

    for (j=0; j < ml; j++) {
      fstinv_(b[j], &n, z, &nn);
    }



    umax = 0.0;
    for (j=0; j < ml; j++) {
      for (i=0; i < m; i++) {
        if (b[j][i] > umax) umax = b[j][i];
      }
    }

    printf (" umax = %e \n",umax);
    MPI_Finalize();
    return 0;
  }

  void transpose (Real **b, int size, int *len, int *disp, int rank, int m){
    int i, *sendcounts, *rdispls;
    Real  *sendbuf, *recvbuf;
    int index = 0;
    sendbuf = createRealArray (m * len[rank]);
    recvbuf = createRealArray (m * len[rank]);
    sendcounts = calloc(size,sizeof(int));
    rdispls = calloc(size,sizeof(int));
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < len[rank]; ++j)
      {
        for (int k = 0; k < len[i]; ++k)
        {

          sendbuf[index] = b[j][disp[i]+k];
          index++;
        }
      }
    }
    index=0;
    for (int i = 0; i < size; ++i)
    {
      sendcounts[i]= len[rank]*len[i];
      rdispls[i]=index;
      index=index+sendcounts[i];
    }
    MPI_Alltoallv(sendbuf, sendcounts, rdispls, MPI_DOUBLE, recvbuf, sendcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    index = 0;
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < len[i]; ++j)
      {
        for (int k = 0; k < len[rank]; ++k)
        {

          b[k][disp[i]+j]=recvbuf[index];
          index++;
        }
      }
    }
  }

  Real *createRealArray (int n)
  {
    Real *a;
    int i;
    a = (Real *)malloc(n*sizeof(Real));
    for (i=0; i < n; i++) {
      a[i] = 0.0;
    }
    return (a);
  }

  Real **createReal2DArray (int n1, int n2)
  {
    int i, n;
    Real **a;
    a    = (Real **)malloc(n1   *sizeof(Real *));
    a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
    for (i=1; i < n1; i++) {
      a[i] = a[i-1] + n2;
    }
    n = n1*n2;
    memset(a[0],0,n*sizeof(Real));
    return (a);
  }

  void splitVector(int globLen, int size, int** len, int** displ)
  {
    *len = calloc(size,sizeof(int));
    *displ = calloc(size,sizeof(int));
    for (int i=0;i<size;++i) {
      (*len)[i] = globLen/size;
      if (globLen % size && i >= (size - globLen % size))
        (*len)[i]++;
      if (i < size-1)
        (*displ)[i+1] = (*displ)[i]+(*len)[i];
    }
  }


  void init_app(int argc, char** argv, int* rank, int* size)
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
  }

