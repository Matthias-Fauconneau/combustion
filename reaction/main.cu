#include "reaction.cu"
__global__ void kernel(const float* state, float* rates) {
	const uint i = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
	reaction(state+i, rates+i);
}
