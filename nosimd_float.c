#include <stdio.h>
#include <stdlib.h>
// no simd

void vec_add(const size_t n, double *z, const double *x, const double *y){
  for(size_t i=0; i<n; ++i)
    z[i] = x[i] + y[i];
}

int main(void){
  const size_t n = 1024000;
  double *x, *y, *z;

  x = (double *)malloc(sizeof(double) * n);
  y = (double *)malloc(sizeof(double) * n);
  z = (double *)malloc(sizeof(double) * n);

  for(size_t i=0; i<n; ++i) x[i] = i;
  for(size_t i=0; i<n; ++i) y[i] = i+1;
  for(size_t i=0; i<n; ++i) z[i] = 0.0;

  vec_add(n, z, x, y);
  
  // for(size_t i=0; i<n; ++i) printf("%g\n", z[i]);

  free(x);
  free(y);
  free(z);

  return 0;
}
