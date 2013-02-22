#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include <mpi.h>
#define M_PI 3.14159265358979323846


Vector genVector(int start, int stop)
{
  int len = stop-start;
  int sumIter = start+1;
  Vector vec = createVector(len);
  for (int i = 0; i < len; ++i)
  {
    vec->data[i]=1./((sumIter)*(sumIter));
    sumIter +=1;
  }
  return vec;
}

double doSum(Vector vec){
  double sum=0;
  #pragma omp parallel for schedule(static) reduction(+:sum)
  for (int i = 0; i < vec->len; ++i)
  {
    sum += vec->data[i];
  }
  freeVector(vec);
  return sum;
}

int main(int argc, char** argv)
{
  int size, rank;
  #ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  if (!(size & (size-1))==0) {
    printf("Number of processes must be power of two");
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return 1;
  }

  double time = WallTime();
  double Sn=(M_PI*M_PI)/6;
  double sum=0;

  for (int i = 4; i <15 ; ++i)
  {
    int n= pow(2, i);
    int *startIndex, *len;
    splitVector(n, size, &len, &startIndex);
    Vector vec = genVector(startIndex[rank],startIndex[rank]+len[rank]);
    sum = doSum(vec);

    #ifdef HAVE_MPI
    double s2=sum;
    MPI_Reduce(&s2, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    #endif

    if (rank == 0)
    {
      printf("Diff (n=%d) = %f",n, sum-Sn);
      printf(", Elapsed: %fs\n", WallTime()-time);
    }
  }



  #ifdef HAVE_MPI
  MPI_Finalize();
  #endif
  return 0;
}
