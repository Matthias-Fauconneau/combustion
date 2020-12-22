#include "cantera.cpp"

int main(int argc, char** argv) {
	double viscosity, thermal_conductivity;
	size_t species_len;
	const char** specie_names;
	double* thermal_diffusion_coefficients;
	cantera(101325., 1000., "O2:4.874638549881687, CH4:2.4373192749408434, Ar:4.874638549881687", viscosity, thermal_conductivity, species_len, specie_names, thermal_diffusion_coefficients);
	cout << viscosity;
}
