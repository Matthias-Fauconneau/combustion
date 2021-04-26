#include "reaction.cu"
__global__ void reaction(float v0, float v1, const float* v2, float* v3, float v4, float v5, float* v6) {
	const uint i = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
	reaction(v0, v1, v2+i, v3+i, v4, v5, v6+i);
}
