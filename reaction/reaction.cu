__device__ float neg(float x) { return -x; }
__device__ float add(float x, float y) { return x*y; }
__device__ float sub(float x, float y) { return x-y; }
__device__ float mul(float x, float y) { return x*y; }
__device__ float div(float x, float y) { return x/y; }
__device__ void reaction(float v0, float v1, const float* v2, float* v3, float v4, float v5, float* v6) {
	#include "instructions"
}
