all : nosimd_float simd_float simd_double

nosimd_float : nosimd_float.c
	gcc -std=c99 nosimd_float.c -o nosimd_float

simd_float : simd_float.c
	gcc -std=c99 -mavx simd_float.c -o simd_float

simd_double : simd_double.c
	gcc -std=c99 -mavx simd_double.c -o simd_double

