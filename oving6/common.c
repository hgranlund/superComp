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

void gatherMatrix(Real** Matrix, int matrixSize, Real* gatherRecvBuf, int* len, int* disp, int root){
  int size, rank, *sendcounts, *rdispls, index;
  Real *gatherSendBuf;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  gatherSendBuf = createRealArray (matrixSize * len[rank]);
  for (int i = 0; i < len[rank]; ++i)
  {
    for (int j = 0; j < matrixSize; ++j)
    {
      gatherSendBuf[i*matrixSize+j]=Matrix[i][j];
    }
  }
  sendcounts = calloc(size,sizeof(int));
  rdispls = calloc(size,sizeof(int));
  index=0;
  for (int i = 0; i < size; ++i)
  {
    sendcounts[i]= len[i]*matrixSize;
    rdispls[i]=i*matrixSize;
  }
  MPI_Gatherv(gatherSendBuf, matrixSize * len[rank], MPI_DOUBLE, gatherRecvBuf, sendcounts, rdispls, MPI_DOUBLE, 0,MPI_COMM_WORLD );
}

Real maxMatrix(Real** Matrix, int matrixSize, int* len, int root){
  int rank;
  Real uMax, localuMax;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  localuMax = 0.0;
  for (int j=0; j < len[rank]; j++) {
    for (int i=0; i < matrixSize; i++) {
      if (Matrix[j][i] > localuMax) localuMax = Matrix[j][i];
    }
  }
  MPI_Reduce(&localuMax, &uMax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  return uMax;
}

void printRoot(const char *string, int rank)
{
  if (rank==0)
  {
    printf("%s",string );
  }
}


void init_app(int argc, char** argv, int* rank, int* size)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, size);
  MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

  // for (int i = 0; i < matrixSize * len[rank]; ++i)
  // {
  //   printf("%d:%e \n ", i, Matrix[i]);
  // }



