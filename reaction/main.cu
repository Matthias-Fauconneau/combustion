__device__ float fneg(float x) { return -x; }
__device__ float fadd(float x, float y) { return x*y; }
__device__ float fsub(float x, float y) { return x-y; }
__device__ float fmul(float x, float y) { return x*y; }
__device__ float fdiv(float x, float y) { return x/y; }
__global__ void reaction(const float /*constant*/v0, const float* /*state*/v1, float* /*dtT_dtS_dtn*/v2) {
#include "instructions"
}
