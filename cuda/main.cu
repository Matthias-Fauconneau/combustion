__device__ float sq(float x) { return x*x; }
__device__ float dot5(float a[5], float b[5]) { float sum = 0.; for(uint k=0;k<5;k++) { sum += a[k]*b[k]; } return sum; }
__device__ float eval_poly(float P[5], float x) {
	float generate_k[5]; for(uint k=0;k<5;k++) generate_k[k] = pow(x, k);
	return dot5(P, generate_k);
}

__device__ float dot(float a[SPECIES], float b[SPECIES]) { float sum = 0.; for(uint k=0;k<SPECIES;k++) { sum += a[k]*b[k]; } return sum; }

__constant__ 	float molar_mass[SPECIES];
__constant__	float sqrt_viscosity_T14[SPECIES][5];
__constant__ float thermal_conductivity_T12[SPECIES][5];
__constant__ float binary_thermal_diffusion_coefficients_T32[SPECIES][SPECIES][5];

__device__ float sqrt_viscosity(uint a, float T) { return sqrt(sqrt(T)) * eval_poly(sqrt_viscosity_T14[a], log(T)); }
__device__ float thermal_conductivity(uint a, float T) { return sqrt(T) * eval_poly(thermal_conductivity_T12[a], log(T)); }
__device__ float binary_thermal_diffusion_coefficient(uint a, uint b, float T) { return T*sqrt(T) * eval_poly(binary_thermal_diffusion_coefficients_T32[a>b?a:b][a>b?b:a], log(T)); }

extern "C" __global__ void rates_transport(
		const uint len,
		const float pressure_R,
		float* temperature, float* _amounts,
		float* viscosity, float* _thermal_conductivity, float* mixture_molar_averaged_thermal_diffusion_coefficients) {
	const float volume = 1.;
	const uint i = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
	float T = temperature[i];
	float amounts[SPECIES]; for(int k=0;k<SPECIES;k++) amounts[k] = _amounts[k*len+i];
	float generate_k[SPECIES]; for(uint k=0;k<SPECIES;k++) {
		float generate_j[SPECIES];
		for(uint j=0;j<SPECIES;j++) {
			generate_j[j] = sq(1. + (sqrt_viscosity(k, T)/sqrt_viscosity(j, T)) * sqrt(sqrt(molar_mass[j]/molar_mass[k]))) / (sqrt(8.) * sqrt(1. + molar_mass[k]/molar_mass[j]));
		}
		generate_k[k] = sq(sqrt_viscosity(k, T)) / dot(amounts, generate_j);
	}
	viscosity[i] = dot(amounts, generate_k);

	float amount = pressure_R * volume / T;
	float dot_amounts_thermal_conductivities = 0;
	float dot_amounts_rcp_thermal_conductivities = 0;
	for(uint k=0;k<SPECIES;k++) {
		float thermal_conductivity_k_T = thermal_conductivity(k, T);
		dot_amounts_thermal_conductivities += amounts[k] * thermal_conductivity_k_T;
		dot_amounts_rcp_thermal_conductivities += amounts[k] / thermal_conductivity_k_T;
	}
	_thermal_conductivity[i] = 1./2. * (dot_amounts_thermal_conductivities / amount + amount / dot_amounts_rcp_thermal_conductivities);
	for(uint k=0;k<SPECIES;k++) {
		float generate_j[SPECIES];
		for(uint j=0;j<SPECIES;j++) {
			if (j != k) {
				generate_j[j] = 1. / binary_thermal_diffusion_coefficient(k, j, T);
			} else {
				generate_j[j] = 0.;
			}
		}
		mixture_molar_averaged_thermal_diffusion_coefficients[k*len+i] = (1. - amounts[k]/amount) / dot(amounts, generate_j);
	}
}
