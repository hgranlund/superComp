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

Vector createPointerVector(int len)
{
  Vector result = (Vector)malloc(sizeof(vector_t));
  result->data = calloc(len, sizeof(double*));
  result->len = len;
  result->stride = 1;
  return result;
}

double doSum(Vector vec){
  double sum=0;
  double one=1.;
  //Generating a vector of one's, to take advantace of blas-ddot
  Vector oneVec = createPointerVector(vec->len);
  for (int i = 0; i < vec->len; ++i)
  {
    oneVec->data[i]=*&one;
  }
  for (int i = 0; i < vec->len; ++i)
  {
      sum = innerproduct(vec, oneVec);
  }
  freeVector(vec);
  return sum;
}

void printDiff(){
  double Sn=(M_PI*M_PI)/6;
  double sum=0;
  for (int i = 4; i < 15; ++i)
  {
    sum = doSum(genVector(pow(2, i)));
    printf("Diff (n=%f) = %f\n",pow(2, i), sum-Sn);
  }
}

int main(int argc, char** argv)
{
  double time = WallTime();
  printDiff();
  printf("elapsed: %f\n", WallTime()-time);
  return 0;
}
