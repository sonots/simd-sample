all : single simd_single

single : single.c
	gcc -std=c99 single.c -o single

simd_single : simd_single.c
	gcc -std=c99 -mavx simd_single.c -o simd_single

