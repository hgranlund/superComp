#include <stdlib.h>
#include <stdio.h>

#include "common.h"



// calculates \sum_K v_i'*A*v_i
Vector genVec(int n)
{
    Vector vec = createVector(n);
    for (int i = 0; i < n; ++i)
    {
      vec->data[i]=1./((i+1)*(i+1));
    }
    return vec;
}

double doSum(Vector vec){
  double sum=0;

  for (int i = 0; i < vec->len; ++i)
  {
      sum += vec->data[i];
  }
  return sum;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    return 1;
  }
  int N=atoi(argv[1]);

  double time = WallTime();
  Vector vec = genVec(N);
  double sum = doSum(vec);

  printf("sum: %f\n", sum);
  printf("elapsed: %f\n", WallTime()-time);

  freeVector(vec);
  return 0;
}
