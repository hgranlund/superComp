#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

#define M_PI 3.14159265358979323846

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
  freeVector(vec);
  return sum;
}

void printDiff(int n){
  double Sn=(M_PI*M_PI)/6;
  double sum=0;
  for (int i = 4; i < 15; ++i)
  {
    sum = doSum(genVec(pow(2, i)));
    printf("Diff n=%f: %f\n",pow(2, i), sum-Sn);
  }
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
  printDiff(N);
  printf("sum: %f\n", sum);
  printf("elapsed: %f\n", WallTime()-time);

  return 0;
}
