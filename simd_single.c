#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>//AVX: -mavx

void vec_add(const size_t n, float *z, const float *x, const float *y){
  static const size_t single_size = 8;
  const size_t end = n / single_size; 

  __m256 *vz = (__m256 *)z;
  __m256 *vx = (__m256 *)x;
  __m256 *vy = (__m256 *)y;
  
  for(size_t i=0; i<end; ++i)
    vz[i] = _mm256_add_ps(vx[i], vy[i]);
}

int main(void){
  const size_t n = 1024000;
  float *x, *y, *z;

  x = (float *)_mm_malloc(sizeof(float) * n, 32);
  y = (float *)_mm_malloc(sizeof(float) * n, 32);
  z = (float *)_mm_malloc(sizeof(float) * n, 32);

  for(size_t i=0; i<n; ++i) x[i] = i;
  for(size_t i=0; i<n; ++i) y[i] = i+1;
  for(size_t i=0; i<n; ++i) z[i] = 0.0;

  vec_add(n, z, x, y);
  
  // for(size_t i=0; i<n; ++i) printf("%g\n", z[i]);

  _mm_free(x);
  _mm_free(y);
  _mm_free(z);

  return 0;
}
