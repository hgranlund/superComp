#include "common.h"



void vectorToMatrix(Real **matrix, Real *vector, int *len, int *disp, int size, int rank ){
  int index = 0;
  for (int i = 0; i < size; ++i)
  {
    for (int j = 0; j < len[i]; ++j)
    {
      for (int k = 0; k < len[rank]; ++k)
      {

        matrix[k][disp[i]+j]=vector[index];
        index++;
      }
    }
  }
}
void matrixToVector(Real **matrix, Real *vector, int *len, int *disp, int size, int rank ){
  int index = 0;
  for (int i = 0; i < size; ++i)
  {
    for (int j = 0; j < len[rank]; ++j)
    {
      for (int k = 0; k < len[i]; ++k)
      {

        vector[index] = matrix[j][disp[i]+k];
        index++;
      }
    }
  }
}

void transpose (Real **b, int size, int *len, int *disp, int rank, int m){
  int i, *sendcounts, *rdispls;
  Real  *sendbuf, *recvbuf;
  sendbuf = createRealArray (m * len[rank]);
  recvbuf = createRealArray (m * len[rank]);
  sendcounts = calloc(size,sizeof(int));
  rdispls = calloc(size,sizeof(int));
  matrixToVector(b,sendbuf,len,disp, size, rank);

  int index = 0;
  for (int i = 0; i < size; ++i)
  {
    sendcounts[i]= len[rank]*len[i];
    rdispls[i]=index;
    index=index+sendcounts[i];
  }
  MPI_Alltoallv(sendbuf, sendcounts, rdispls, MPI_DOUBLE, recvbuf, sendcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
  vectorToMatrix(b,recvbuf,len,disp, size, rank);
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



