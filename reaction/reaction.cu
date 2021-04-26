__device__ float neg(float x) { return -x; }
__device__ float add(float x, float y) { return x*y; }
__device__ float sub(float x, float y) { return x-y; }
__device__ float mul(float x, float y) { return x*y; }
__device__ float div(float x, float y) { return x/y; }
__device__ void reaction(const float* const v0, float* const v1) {
	#include "instructions"
}
