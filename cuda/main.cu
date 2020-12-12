double log(double x) { return log(float(x)); }
double exp(double x) { return exp(float(x)); }

const double ideal_gas_constant = 8.31446261815324;

const double reference_pressure_R = 101325. / ideal_gas_constant;
const double T_split = 1000.;
__device__ /*&*/const double/*[7]*/* A(/*&'static*/const double s[2][7], double T) { if (T < T_split) { return s[0]; } else { return s[1]; } }
__device__ double specific_heat_capacity(const double s[2][7], double T) { const double/*[7]*/* a = A(s, T); return a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T; } // /R
__device__ double specific_enthalpy(const double s[2][7], double T) { const double/*[7]*/* a = A(s, T); return a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T; } // /R
__device__ double specific_entropy(const double s[2][7], double T) { const double/*[7]*/* a = A(s, T); return a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T; } // /R

const uint S = SPECIES;

__device__ double dot(double a[S], double b[S]) { double sum = 0.; for(uint k=0;k<S;k++) { sum += a[k]*b[k]; } return sum; }
__device__ double dot1(double a[S-1], double b[S-1]) { double sum = 0.; for(uint k=0;k<S-1;k++) { sum += a[k]*b[k]; } return sum; }

struct RateConstant {
	double log_preexponential_factor;
	double temperature_exponent;
	double activation_temperature;
};

__device__ double log_arrhenius(RateConstant r, double T) { return r.log_preexponential_factor + r.temperature_exponent*log(T) - r.activation_temperature*(1./T); }

struct Reaction {
	double reactants[S];
	double products[S];
	double net[S-1];
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
	double efficiencies[S];
};

struct Falloff {
	Reaction reaction;
	double efficiencies[S];
	RateConstant k0;
	double A;
	double T3;
	double T1;
	double T2;
};

double __device__ falloff_efficiency(Falloff f, double T, double concentrations[S], double log_k_inf) {
	double Pr = dot(f.efficiencies, concentrations) * exp(log_arrhenius(f.k0, T) - log_k_inf);
	double Fcent = (1.-f.A)*exp(-T/f.T3)+f.A*exp(-T/f.T1)+exp(-f.T2/T);
	double log10Fcent = log10(Fcent);
	double C = -0.4-0.67*log10Fcent;
	double N = 0.75-1.27*log10Fcent;
	double log10PrC = log10(Pr) + C;
	double f1 = log10PrC/(N-0.14*log10PrC);
	double F = exp10(log10Fcent/(1.+f1*f1));
	return Pr / (1.+Pr) * F;
}

__device__ void reaction(double (&net_rates)[S-1], Reaction r, double logP0_RT, double G[S-1], double log_kf, double log_concentrations[S], double c) {
	double Rf = exp(dot(r.reactants, log_concentrations) + log_kf);
	double log_equilibrium_constant = -dot(r.net, G) + r.sum_net*logP0_RT;
	double Rr = exp(dot(r.products, log_concentrations) + log_kf - log_equilibrium_constant);
	double R = Rf - Rr;
	double cR = c * R;
	for(uint k=0;k<S-1;k++) {
		double net = r.net[k];
		net_rates[k] += net * cR;
	}
}

struct System {
	double molar_masses[S];
	double thermodynamics[S][2][7];
	Elementary elementary[ELEMENTARY];
	ThreeBody three_body[THREE_BODY];
	Falloff falloff[FALLOFF];
};

extern "C" __global__ void dt(const size_t len, const double volume, const double pressure_R, double* temperature, double* amounts) {
const System system = System{
#include "system.h"
};
for (uint i = blockIdx.x * blockDim.x + threadIdx.x; i < len; i += blockDim.x * gridDim.x) {
	double V = volume;
	double T = temperature[i];
	double C = pressure_R / T;
	double logP0_RT = log(reference_pressure_R) - log(T);
	double H[S-1];
	for(uint k=0;k<S-1;k++) H[k] = specific_enthalpy(system.thermodynamics[k], T);
	double H_T[S-1];
	for(uint k=0;k<S-1;k++) H_T[k] = H[k] / T;
	double G[S-1];
	for(uint k=0;k<S-1;k++) G[k] = H_T[k] - specific_entropy(system.thermodynamics[k], T);
	double active_concentrations[S-1];
	for(int k=0;k<S-1;k++) active_concentrations[k] = amounts[k*len+i] / V;
	double sumC = 0.;
	for(int k=0;k<S-1;k++) sumC += active_concentrations[k];
	double Ca = C - sumC;
	double concentrations[S];
	for(int k=0;k<S;k++) concentrations[k] = active_concentrations[k];
	concentrations[S-1] = Ca;
	double log_concentrations[S];
	for(uint k=0;k<S;k++) log_concentrations[k] = log(concentrations[k]);
	double net_rates[S-1];
	for(uint k=0;k<S-1;k++) net_rates[k] = 0.;
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
	for(uint j=0;j<FALLOFF;j++) {
		Falloff r = system.falloff[j];
		double log_kf = log_arrhenius(r.reaction.rate_constant, T);
		reaction(net_rates, r.reaction, logP0_RT, G, log_kf, log_concentrations, falloff_efficiency(r, T, concentrations, log_kf));
	}
	double Cp[S];
	for(uint k=0;k<S;k++) Cp[k] = specific_heat_capacity(system.thermodynamics[k], T);
	double rcp_sumCCp = 1. / dot(concentrations, Cp);
	double dtT_T = - rcp_sumCCp * dot1(H_T, net_rates);
	temperature[i] = dtT_T;
	double dtn[S-1];
	for(uint k=0;k<S-1;k++) dtn[k] = V*net_rates[k];
	for(uint k=0;k<S-1;k++) amounts[k*len+i] = dtn[k];
}
}
