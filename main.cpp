#include "table.h"
//static double dot5(const double x[5], const double y[5]) { return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3] + x[4]*y[4]; }
static double poly5(double x, const double c[5]) { return c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x + c[4]*x*x*x*x + c[5]*x*x*x*x*x; }

typedef unsigned int uint;
#include <Eigen/Dense>
static void polyfit(uint n, uint deg, const double* xp, const double* yp, const double* wp, double* pp) {
	using namespace Eigen;
	Map<const VectorXd> x(xp, n);
	VectorXd y = Map<const VectorXd>(yp, n);
	Map<VectorXd> p(pp, deg+1);

	// Construct A such that each row i of A has the elements
	// 1, x[i], x[i]^2, x[i]^3 ... + x[i]^deg
	MatrixXd A(n, deg+1);
	A.col(0).setConstant(1.0);

	if (deg > 0) {
			A.col(1) = x;
	}
	for (uint i = 1; i < deg; i++) {
			A.col(i+1) = A.col(i).array() * x.array();
	}

	if (wp != nullptr && wp[0] > 0) {
			// For compatibility with old Fortran dpolft, input weights are the
			// squares of the weight vector used in this algorithm
			VectorXd w = Map<const VectorXd>(wp, n).cwiseSqrt().eval();

			// Multiply by the weights on both sides
			A = w.asDiagonal() * A;
			y.array() *= w.array();
	}

	// Solve W*A*p = W*y to find the polynomial coefficients
	p = A.colPivHouseholderQr().solve(y);
}

static void fitDelta(int table, int ntstar, int degree, double* c) {
	double* begin = 0;
	switch (table) {
	case 0:
			begin = omega22_table + 8*ntstar;
			break;
	case 1:
			begin = astar_table + 8*(ntstar + 1);
			break;
	case 2:
			begin = bstar_table + 8*(ntstar + 1);
			break;
	case 3:
			begin = cstar_table + 8*(ntstar + 1);
			break;
	default:
			return;
	}
	double delta[8] = {0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5};
	polyfit(8, degree, delta, begin, nullptr, c);
}

static double quadInterp(double x0, const double* x, const double* y) {
    double dx21 = x[1] - x[0];
    double dx32 = x[2] - x[1];
    double dx31 = dx21 + dx32;
    double dy32 = y[2] - y[1];
    double dy21 = y[1] - y[0];
    double a = (dx21*dy32 - dy21*dx32)/(dx21*dx31*dx32);
    return a*(x0 - x[0])*(x0 - x[1]) + (dy21/dx21)*(x0 - x[1]) + y[1];
}

static double omega22(const double (logTemp)[37], const double (o22poly)[37][5], double ts, double deltastar) {
    int i;
    for (i = 0; i < 37; i++) {
        if (ts < tstar22[i]) {
            break;
        }
    }
    int i1 = std::max(i - 1, 0);
    int i2 = i1+3;
    if (i2 > 36) {
        i2 = 36;
        i1 = i2 - 3;
    }
    double values[3];
    for (i = i1; i < i2; i++) {
        if (deltastar == 0.0) {
            values[i-i1] = omega22_table[8*i];
        } else {
            values[i-i1] = poly5(deltastar, o22poly[i]);
        }
    }
    return quadInterp(log(ts), &logTemp[i1], values);
}

static double astar(const double (logTemp)[37], const double (apoly)[37][5], double ts, double deltastar) {
    int i;
    for (i = 0; i < 37; i++) if (ts < tstar22[i]) {
            break;
        }
    int i1 = std::max(i - 1, 0);
    int i2 = i1+3;
    if (i2 > 36) {
        i2 = 36;
        i1 = i2 - 3;
    }
    double values[3];
    for (i = i1; i < i2; i++) {
        if (deltastar == 0.0) {
            values[i-i1] = astar_table[8*(i + 1)];
        } else {
            values[i-i1] = poly5(deltastar, apoly[i]);
        }
    }
    return quadInterp(log(ts), &logTemp[i1], values);
}

#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/transport/MultiTransport.h>
#include <cantera/transport/TransportData.h>
using namespace Cantera;
static double vismix(std::shared_ptr<ThermoPhase> thermo) {
	auto S = thermo->nSpecies();
	std::cout << S << std::endl;
	auto minTemp = thermo->minTemp();
	auto maxTemp = thermo->maxTemp();
	std::cout << minTemp << " " << maxTemp << std::endl;
	const auto& W = thermo->molecularWeights();
	//std::cout << W << std::endl;
	double molefracs[S];
	thermo->getMoleFractions(z);
	//std::cout << molefracs << std::endl;
	const uint np = 50;
	const double dt = (maxTemp - minTemp)/(np-1);
	double cp_R_all[S][np];
	double crot[S], sigma[S], eps[S], dipole[S][S], polar[S], alpha[S], zrot[S];
	for(uint k = 0; k < S; k++) {
		shared_ptr<Species> s = thermo->species(k);
		const GasTransportData* sptran = dynamic_cast<GasTransportData*>(s->transport.get());
		if (sptran->geometry == "atom") { crot[k] = 0.0; }
		else if (sptran->geometry == "linear") { crot[k] = 1.0; }
		else if (sptran->geometry == "nonlinear") { crot[k] = 1.5; }
		sigma[k] = sptran->diameter;
		eps[k] = sptran->well_depth;
		dipole[k][k] = sptran->dipole;
		polar[k] = (sptran->dipole > 0);
		alpha[k] = sptran->polarizability;
		zrot[k] = sptran->rotational_relaxation;
		for(uint n = 0; n < np; n++) {
			double cp_R_all_tmp[S];
			double t = thermo->minTemp() + dt*n;
			thermo->setTemperature(t);
			thermo->getCp_R_ref(cp_R_all_tmp);
			thermo->setTemperature(1000.);
			cp_R_all[k][n] = cp_R_all_tmp[k];
		}
	}
	/*std::cout << crot << std::endl;
	std::cout << eps << std::endl;
	std::cout << dipole << std::endl;
	std::cout << polar << std::endl;
	std::cout << alpha << std::endl;
	std::cout << zrot << std::endl;*/

	double epsilon[S][S], delta[S][S], reducedMass[S][S], diam[S][S];
	double f_eps, f_sigma;

	for (uint i = 0; i < S; i++) {
		for (uint j = i; j < S; j++) {
			reducedMass[i][j] = W[i] * W[j] / (Avogadro * (W[i] + W[j]));

			// hard-sphere diameter for [i][j] collisions
			diam[i][j] = 0.5*(sigma[i] + sigma[j]);

			// the effective well depth for [i][j] collisions
			epsilon[i][j] = sqrt(eps[i]*eps[j]);

			// the effective dipole moment for [i][j] collisions
			dipole[i][j] = sqrt(dipole[i][i]*dipole[j][j]);

			// reduced dipole moment delta* (nondimensional)
			double d = diam[i][j];
			const double epsilon_0 = 1.0 / (lightSpeed * lightSpeed * permeability_0);
			delta[i][j] = 0.5 * dipole[i][j]*dipole[i][j] / (4 * Pi * epsilon_0 * epsilon[i][j] * d * d * d);
			// no correction if both are nonpolar, or both are polar
			if (polar[i] == polar[j]) {
					f_eps = 1.0;
					f_sigma = 1.0;
			} else {
				// corrections to the effective diameter and well depth
				// if one is polar and one is non-polar
				uint kp = (polar[i] ? i : j); // the polar one
				uint knp = (i == kp ? j : i); // the nonpolar one
				double d3np, d3p, alpha_star, mu_p_star, xi;
				d3np = pow(sigma[knp],3);
				d3p = pow(sigma[kp],3);
				alpha_star = alpha[knp]/d3np;
				mu_p_star = dipole[kp][kp]/sqrt(4 * Pi * epsilon_0 * d3p * eps[kp]);
				xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star * sqrt(eps[kp]/eps[knp]);
				f_sigma = pow(xi, -1.0/6.0);
				f_eps = xi*xi;
			}

			diam[i][j] *= f_sigma;
			epsilon[i][j] *= f_eps;

			// properties are symmetric
			reducedMass[j][i] = reducedMass[i][j];
			diam[j][i] = diam[i][j];
			epsilon[j][i] = epsilon[i][j];
			dipole[j][i] = dipole[i][j];
			delta[j][i] = delta[i][j];
		}
	}

	double tstar_min = 1.e8, tstar_max = 0.0;
	for (uint i = 0; i < S; i++) {
			for (uint j = i; j < S; j++) {
					// The polynomial fits of collision integrals vs. T*
					// will be done for the T* from tstar_min to tstar_max
					tstar_min = std::min(tstar_min, Boltzmann * minTemp/epsilon[i][j]);
					tstar_max = std::max(tstar_max, Boltzmann * maxTemp/epsilon[i][j]);
			}
	}

	auto nmin = -1;
	auto nmax = -1;

	double tstar[39] = {
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
    18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0, 500.0
	};

	for (int n = 0; n < 37; n++) {
			if (tstar_min > tstar[n+1]) {
					nmin = n;
			}
			if (tstar_max > tstar[n+1]) {
					nmax = n+1;
			}
	}
	if (nmin < 0 || nmin >= 36 || nmax < 0 || nmax > 36) {
			nmin = 0;
			nmax = 36;
	}
	double logTemp[37];
	double o22poly[37][5];
	double apoly[37][5];
	double bpoly[37][5];
	double cpoly[37][5];
	for (int i = 0; i < 37; i++) {
		logTemp[i] = log(tstar[i+1]);
		fitDelta(0, i, 5, o22poly[i]);
		fitDelta(1, i, 5, apoly[i]);
		fitDelta(2, i, 5, bpoly[i]);
		fitDelta(3, i, 5, cpoly[i]);
	}

	// number of points to use in generating fit data
	double tlog[np], spvisc[np], spcond[np];
	double w[np], w2[np];

	// generate array of log(t) values
	for (uint n = 0; n < np; n++) {
			double t = minTemp + dt*n;
			tlog[n] = log(t);
	}

	double diff[np + 1];
	double visccoeffs[S][5], condcoeffs[S][5], diffcoeffs[S][S][5];

	for (uint k = 0; k < S; k++) {
		const double tstar = Boltzmann * 298.0 / eps[k];
		// Scaling factor for temperature dependence of z_rot. [Kee2003] Eq.
		// 12.112 or [Kee2017] Eq. 11.115
		double fz_298 = 1.0 + pow(Pi, 1.5) / sqrt(tstar) * (0.5 + 1.0 / tstar) +
				(0.25 * Pi * Pi + 2) / tstar;

		for (uint n = 0; n < np; n++) {
			const double t = thermo->minTemp() + dt*n;
			const double cp_R = cp_R_all[k][n];
			const double tstar = Boltzmann * t / eps[k];
			double sqrt_T = sqrt(t);
			double om22 = omega22(logTemp, o22poly, tstar, delta[k][k]);
			double om11 = omega22(logTemp, o22poly, tstar, delta[k][k])/astar(logTemp, apoly, tstar, delta[k][k]);

			// self-diffusion coefficient, without polar corrections
			double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/reducedMass[k][k]) * pow((Boltzmann * t), 1.5)/ (Pi * sigma[k] * sigma[k] * om11);

			// viscosity
			double visc = 5.0/16.0 * sqrt(Pi * W[k] * Boltzmann * t / Avogadro) / (om22 * Pi * sigma[k]*sigma[k]);

			// thermal conductivity
			double f_int = W[k]/(GasConstant * t) * diffcoeff/visc;
			double cv_rot = crot[k];
			double A_factor = 2.5 - f_int;
			double fz_tstar = 1.0 + pow(Pi, 1.5) / sqrt(tstar) * (0.5 + 1.0 / tstar) + (0.25 * Pi * Pi + 2) / tstar;
			double B_factor = zrot[k] * fz_298 / fz_tstar + 2.0/Pi * (5.0/3.0 * cv_rot + f_int);
			double c1 = 2.0/Pi * A_factor/B_factor;
			double cv_int = cp_R - 2.5 - cv_rot;
			double f_rot = f_int * (1.0 + c1);
			double f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);
			double cond = (visc/W[k])*GasConstant*(f_trans * 1.5 + f_rot * cv_rot + f_int * cv_int);

				// the viscosity should be proportional approximately to
				// sqrt(T); therefore, visc/sqrt(T) should have only a weak
				// temperature dependence. And since the mixture rule requires
				// the square root of the pure-species viscosity, fit the square
				// root of (visc/sqrt(T)) to avoid having to compute square
				// roots in the mixture rule.
				spvisc[n] = sqrt(visc/sqrt_T);

				// the pure-species conductivity scales approximately with
				// sqrt(T). Unlike the viscosity, there is no reason here to fit
				// the square root, since a different mixture rule is used.
				spcond[n] = cond/sqrt_T;
				w[n] = 1.0/(spvisc[n]*spvisc[n]);
				w2[n] = 1.0/(spcond[n]*spcond[n]);
		}
		polyfit(np, 4, tlog, spvisc, w, visccoeffs[k]);
		polyfit(np, 4, tlog, spcond, w, condcoeffs[k]);

		for (uint j = k; j < S; j++) {
			for (uint n = 0; n < np; n++) {
				double t = minTemp + dt*n;
				double eps = epsilon[j][k];
				double tstar = Boltzmann * t/eps;
				double sigma = diam[j][k];
				//double om11 = integrals.omega11(tstar, delta[j][k]);
				//double omega11(double ts, double deltastar) { return omega22(ts, deltastar)/astar(ts, deltastar); }
				double om11 = omega22(logTemp, o22poly, tstar, delta[j][k])/astar(logTemp, apoly, tstar, delta[j][k]);
				double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/reducedMass[k][j]) * pow(Boltzmann * t, 1.5) / (Pi * sigma * sigma * om11);
				diff[n] = diffcoeff/pow(t, 1.5);
				w[n] = 1.0/(diff[n]*diff[n]);
			}
			polyfit(np, 4, tlog, diff, w, diffcoeffs[k][j]);
		}
	}

	double phi[S][S];
	double wratjk[S][S], wratkj1[S][S];
	for (uint j = 0; j < S; j++) {
			for (uint k = j; k < S; k++) {
					wratjk[j][k] = sqrt(W[j]/W[k]);
					wratjk[k][j] = sqrt(wratjk[j][k]);
					wratkj1[j][k] = sqrt(1.0 + W[k]/W[j]);
			}
	}

	double T = 1000.;
	double polytempvec[5];
	polytempvec[0] = 1.0;
	polytempvec[1] = log(T);
	polytempvec[2] = log(T)*log(T);
	polytempvec[3] = log(T)*log(T)*log(T);
	polytempvec[4] = log(T)*log(T)*log(T)*log(T);

	double sqvisc[S], visc[S];
	for (uint k = 0; k < S; k++) {
		// the polynomial fit is done for sqrt(visc/sqrt(T))
		auto sum = 0.;
		for(uint i=0; i<5; i++) sum += polytempvec[i]*visccoeffs[k][i];
		sqvisc[k] = sqrt(sqrt(T)) * sum;
		visc[k] = (sqvisc[k] * sqvisc[k]);
	}
	// see Eq. (9-5.15) of Reid, Prausnitz, and Poling
	for (uint j = 0; j < S; j++) {
		for (uint k = j; k < S; k++) {
			double vratiokj = visc[k]/visc[j];
			double wratiojk = W[j]/W[k];

			// Note that wratjk(k,j) holds the square root of wratjk[j][k]!
			double factor1 = 1.0 + (sqvisc[k]/sqvisc[j]) * wratjk[k][j];
			phi[k][j] = factor1*factor1 / (sqrt(8.0) * wratkj1[j][k]);
			phi[j][k] = phi[k][j]/(vratiokj * wratiojk);
		}
	}

	double spwork[S];
	for(uint i=0; i<S; i++) for(uint j=0; i<S; i++) spwork[i] = phi[i][j] * molefracs[j];

	double vismix = 0.0;
	for (uint k = 0; k < S; k++) {
			vismix += molefracs[k] * visc[k]/spwork[k];
	}
	return vismix;
}

int main(int argc, char** argv) {
	auto mechanism = newSolution("gri30.yaml","gri30","None");
	auto thermo = mechanism->thermo();
	thermo->setState_TPX(101325., 1000., "O2:4.874638549881687, CH4:2.4373192749408434, Ar:4.874638549881687");
	std::cout << vismix(thermo) << std::endl;
}
