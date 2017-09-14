/*--------------------
  simd_double.c
  --------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>//AVX: -mavx

void vec_add(const size_t n, double *z, const double *x, const double *y){
  static const size_t double_size = 4;
  const size_t end = n / double_size;

  __m256d *vz = (__m256d *)z;
  __m256d *vx = (__m256d *)x;
  __m256d *vy = (__m256d *)y;
  
  for(size_t i=0; i<end; ++i)
    vz[i] = _mm256_add_pd(vx[i], vy[i]);
}

int main(void){
  const size_t n = 1024;
  double *x, *y, *z;

  x = (double *)_mm_malloc(sizeof(double) * n, 32);
  y = (double *)_mm_malloc(sizeof(double) * n, 32);
  z = (double *)_mm_malloc(sizeof(double) * n, 32);

  for(size_t i=0; i<n; ++i) x[i] = i;
  for(size_t i=0; i<n; ++i) y[i] = i+1;
  for(size_t i=0; i<n; ++i) z[i] = 0.0;

  vec_add(n, z, x, y);
  
  for(size_t i=0; i<n; ++i) printf("%g\n", z[i]);

  _mm_free(x);
  _mm_free(y);
  _mm_free(z);

  return 0;
}
