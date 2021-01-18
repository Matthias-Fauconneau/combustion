__device__ double sq(double x) { return x*x; }
__device__ double dot4(double a[4], double b[4]) { double sum = 0.; for(uint k=0;k<4;k++) { sum += a[k]*b[k]; } return sum; }
__device__ double eval_poly(double P[4], double x) {
	double generate_k[4]; for(uint k=0;k<4;k++) generate_k[k] = pow(x, k);
	return dot4(P, generate_k);
}

//const double ideal_gas_constant = 8.31446261815324;

//const double reference_pressure_R = 101325. / ideal_gas_constant;
const double T_split = 1000.;
__device__ /*&*/const double/*[7]*/* A(/*&'static*/const double s[2][7], double T) { if (T < T_split) { return s[0]; } else { return s[1]; } }
__device__ double specific_heat_capacity(const double s[2][7], double T) { const double/*[7]*/* a = A(s, T); return a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T; } // /R
__device__ double specific_enthalpy(const double s[2][7], double T) { const double/*[7]*/* a = A(s, T); return a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T; } // /R
__device__ double specific_entropy(const double s[2][7], double T) { const double/*[7]*/* a = A(s, T); return a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T; } // /R

__device__ double dot(double a[SPECIES], double b[SPECIES]) { double sum = 0.; for(uint k=0;k<SPECIES;k++) { sum += a[k]*b[k]; } return sum; }
__device__ double dot1(double a[SPECIES-1], double b[SPECIES-1]) { double sum = 0.; for(uint k=0;k<SPECIES-1;k++) { sum += a[k]*b[k]; } return sum; }

struct RateConstant {
	double log_preexponential_factor;
	double temperature_exponent;
	double activation_temperature;
};

__device__ double log_arrhenius(RateConstant r, double T) { return r.log_preexponential_factor + r.temperature_exponent*log(T) - r.activation_temperature*(1./T); }

struct Reaction {
	double reactants[SPECIES];
	double products[SPECIES];
	double net[SPECIES-1];
	double sum_reactants;
	double sum_products;
	double sum_net;
	RateConstant rate_constant;
};

struct Elementary {
	Reaction reaction;
};

struct ThreeBody {
	Reaction reaction;
	double efficiencies[SPECIES];
};

struct PressureModification {
	Reaction reaction;
	double efficiencies[SPECIES];
	RateConstant k0;
};

struct Falloff {
	Reaction reaction;
	uint efficiency;
	RateConstant k0;
	double A;
	double T3;
	double T1;
	double T2;
};

__device__ double pressure_modification(PressureModification r, double T, double concentrations[SPECIES], double log_k_inf) {
	double Pr = dot(r.efficiencies, concentrations) * exp(log_arrhenius(r.k0, T) - log_k_inf);
	return Pr / (1.+Pr);
}

__device__ double falloff(Falloff r, double T, double efficiency, double log_k_inf) {
	double Pr = efficiency * exp(log_arrhenius(r.k0, T) - log_k_inf);
	double Fcent = (1.-r.A)*exp(-T/r.T3)+r.A*exp(-T/r.T1)+exp(-r.T2/T);
	double log10Fcent = log10(Fcent);
	double C = -0.4-0.67*log10Fcent;
	double N = 0.75-1.27*log10Fcent;
	double log10PrC = log10(Pr) + C;
	double f1 = log10PrC/(N-0.14*log10PrC);
	double F = exp10(log10Fcent/(1.+f1*f1));
	return Pr / (1.+Pr) * F;
}

__device__ double dot_m(double a[SPECIES], double b[SPECIES]) { double sum = 0.; for(uint k=0;k<SPECIES;k++) { if(a[k]!=0.) sum += a[k]*b[k]; } return sum; } // 0*inf=0

__device__ void reaction(double (&net_rates)[SPECIES-1], Reaction r, double logP0_RT, double G[SPECIES-1], double log_kf, double log_concentrations[SPECIES], double c) {
	double Rf = exp(dot_m(r.reactants, log_concentrations) + log_kf);
	double log_equilibrium_constant = -dot1(r.net, G) + r.sum_net*logP0_RT;
	double Rr = exp(dot_m(r.products, log_concentrations) + log_kf - log_equilibrium_constant);
	double R = Rf - Rr;
	double cR = c * R;
	for(uint k=0;k<SPECIES-1;k++) {
		double net = r.net[k];
		net_rates[k] += net * cR;
	}
}

struct TransportPolynomials {
	double sqrt_viscosity_T14[SPECIES][4];
	double thermal_conductivity_T12[SPECIES][4];
	double binary_thermal_diffusion_coefficients_T32[SPECIES][SPECIES][4];
};

struct System {
	double molar_mass[SPECIES];
	double thermodynamics[SPECIES][2][7];
	Elementary elementary[ELEMENTARY];
	ThreeBody three_body[THREE_BODY];
	PressureModification pressure_modification[PRESSURE_MODIFICATION];
	double efficiencies[EFFICIENCIES][SPECIES];
	Falloff falloff[FALLOFF];
	TransportPolynomials transport_polynomials;
};

__device__ double sqrt_viscosity(System system, uint a, double T) { return sqrt(sqrt(T)) * eval_poly(system.transport_polynomials.sqrt_viscosity_T14[a], log(T)); }
__device__ double thermal_conductivity(System system, uint a, double T) { return sqrt(T) * eval_poly(system.transport_polynomials.thermal_conductivity_T12[a], log(T)); }
__device__ double binary_thermal_diffusion_coefficient(System system, uint a, uint b, double T) { return pow(T,3./2.) * eval_poly(system.transport_polynomials.binary_thermal_diffusion_coefficients_T32[a>b?a:b][a>b?b:a], log(T)); }

extern "C" __global__ void rates_transport(const uint len, const double pressure_R, double* temperature, double* _amounts, double* d_temperature, double* d_amounts, double* viscosity, double* _thermal_conductivity, double* mixture_averaged_thermal_diffusion_coefficients) {
const System system = System{
#include "model.h"
};
const double volume = 1.;
for (uint i = blockIdx.x * /*workgroup size*/blockDim.x + /*SIMD lane*/threadIdx.x; i < len; i += /*workgroup size*/blockDim.x * /*SIMD width*/gridDim.x) {
	double T = temperature[i];
	/*double C = pressure_R / T;
	double logP0_RT = log(reference_pressure_R) - log(T);
	double H[SPECIES-1];
	for(uint k=0;k<SPECIES-1;k++) H[k] = specific_enthalpy(system.thermodynamics[k], T);
	double H_T[SPECIES-1];
	for(uint k=0;k<SPECIES-1;k++) H_T[k] = H[k] / T;
	double G[SPECIES-1];
	for(uint k=0;k<SPECIES-1;k++) G[k] = H_T[k] - specific_entropy(system.thermodynamics[k], T);
	double active_concentrations[SPECIES-1];
	for(int k=0;k<SPECIES-1;k++) active_concentrations[k] = amounts[k*len+i];
	double sumC = 0.;
	for(int k=0;k<SPECIES-1;k++) sumC += active_concentrations[k];
	double Ca = C - sumC;
	double concentrations[SPECIES];
	for(int k=0;k<SPECIES;k++) concentrations[k] = active_concentrations[k];
	concentrations[SPECIES-1] = Ca;
	double log_concentrations[SPECIES];
	for(uint k=0;k<SPECIES;k++) log_concentrations[k] = log(concentrations[k]);
	double net_rates[SPECIES-1];
	for(uint k=0;k<SPECIES-1;k++) net_rates[k] = 0.;
	for(uint j=0;j<ELEMENTARY;j++) {
		Elementary r = system.elementary[j];
		double log_kf = log_arrhenius(r.reaction.rate_constant, T);
		reaction(net_rates, r.reaction, logP0_RT, G, log_kf, log_concentrations, 1.);
	}
	for(uint j=0;j<THREE_BODY;j++) {
		ThreeBody r = system.three_body[j];
		double log_kf = log_arrhenius(r.reaction.rate_constant, T);
		reaction(net_rates, r.reaction, logP0_RT, G, log_kf, log_concentrations, dot(r.efficiencies, concentrations));
	}
	for(uint j=0;j<PRESSURE_MODIFICATION;j++) {
		PressureModification r = system.pressure_modification[j];
		double log_kf = log_arrhenius(r.reaction.rate_constant, T);
		reaction(net_rates, r.reaction, logP0_RT, G, log_kf, log_concentrations, dot(r.efficiencies, concentrations));
	}
	double efficiencies[EFFICIENCIES];
	for(uint j=0;j<EFFICIENCIES;j++) {
		efficiencies[j] = dot(system.efficiencies[j], concentrations);
	}
	for(uint j=0;j<FALLOFF;j++) {
		Falloff r = system.falloff[j];
		double log_kf = log_arrhenius(r.reaction.rate_constant, T);
		reaction(net_rates, r.reaction, logP0_RT, G, log_kf, log_concentrations, falloff(r, T, efficiencies[r.efficiency], log_kf));
	}
	double Cp[SPECIES];
	for(uint k=0;k<SPECIES;k++) Cp[k] = specific_heat_capacity(system.thermodynamics[k], T);
	double rcp_sumCCp = 1. / dot(concentrations, Cp);
	double dtT_T = - rcp_sumCCp * dot1(H_T, net_rates);
	d_temperature[i] = dtT_T;
	for(uint k=0;k<SPECIES-1;k++) d_amounts[k*len+i] = net_rates[k];*/
	double amounts[SPECIES]; for(int k=0;k<SPECIES;k++) amounts[k] = _amounts[k*len+i];
	double generate_k[SPECIES]; for(uint k=0;k<SPECIES;k++) {
		double generate_j[SPECIES];
		for(uint j=0;j<SPECIES;j++) {
			generate_j[j] = sq(1. + (sqrt_viscosity(system, k, T)/sqrt_viscosity(system, j, T)) * sqrt(sqrt(system.molar_mass[j]/system.molar_mass[k]))) /
				(sqrt(8.) * sqrt(1. + system.molar_mass[k]/system.molar_mass[j]));
		}
		generate_k[k] = sq(sqrt_viscosity(system, k, T)) / dot(amounts, generate_j);
	}
#if 0
	viscosity[i] = dot(amounts, generate_k);
	double amount = pressure_R / T * volume;
	double thermal_conductivities[SPECIES]; for(uint k=0;k<SPECIES;k++) { thermal_conductivities[k] = thermal_conductivity(system, k, T); }
	double rcp_thermal_conductivities[SPECIES]; for(uint k=0;k<SPECIES;k++) { rcp_thermal_conductivities[k] = 1. / thermal_conductivity(system, k, T); }
	_thermal_conductivity[i] = 1./2. * (dot(amounts, thermal_conductivities) / amount + amount / dot(amounts, rcp_thermal_conductivities));
	for(uint k=0;k<SPECIES;k++) {
		double generate_j[SPECIES];
		for(uint j=0;j<SPECIES;j++) {
			if (j != k) {
				generate_j[j] = 1. / binary_thermal_diffusion_coefficient(system, k, j, T);
			} else {
				generate_j[j] = 0.;
			}
		}
		mixture_averaged_thermal_diffusion_coefficients[k*len+i] = (1. - amounts[k]/amount) / dot(amounts, generate_j);
	}
#endif
}
}
