#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

#define M_PI 3.14159265358979323846



Vector genVector(int length)
{
    Vector vec = createVector(length);
    for (int i = 0; i < length; ++i)
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
    sum = doSum(genVector(pow(2, i)));
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
  // Vector vec = genVector(N);
  // double sum = doSum(vec);
  // // printf("sum: %f\n", sum);
  printDiff(N);
  printf("elapsed: %f\n", WallTime()-time);

  return 0;
}
