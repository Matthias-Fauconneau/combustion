__device__ float fneg(float x) { return -x; }
__device__ float fadd(float x, float y) { return x*y; }
__device__ float fsub(float x, float y) { return x-y; }
__device__ float fmul(float x, float y) { return x*y; }
__device__ float fdiv(float x, float y) { return x/y; }
__global__ void reaction(/*const float constant,*/ const float* state, float* rates) {
	const uint i = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
	/*const float v0 = constant;
	const float* const v1 = state+i;
	float* const v2 = rates+i;*/
	const float* const v0 = state+i;
	float* const v1 = rates+i;
	#include "instructions"
}
