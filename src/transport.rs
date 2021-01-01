fn quadratic_interpolation(x: [f64; 3], y: [f64; 3], x0: f64) -> f64 { ((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1] }

fn eval_poly<const N: usize>(P: &[f64; N], x: f64) -> f64 { P.dot(generate(|k| x.powi(k as i32))) }

trait Vector<const N: usize> = vec::Vector<N>+iter::IntoIterator<Item=f64>;
 fn weighted_polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>, w: impl Vector<N>) -> [f64; D] {
	use nalgebra::{DMatrix, DVector, SVD};
	let w = DVector::from_iterator(N, w.into_iter());
	let A = DMatrix::from_iterator(N, D, x.into_iter().zip(w.iter()).map(|(x, w)| (0..D).map(move |k| w*x.powi(k as i32))).flatten());
	let b = DVector::from_iterator(N, y.into_iter().zip(w.iter()).map(|(x, w)| w*x));
	use std::convert::TryInto;
	trace::timeout(1, || SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap())
}
// Regression with 1/y² weights (towards relative vertical error)
fn polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Copy) -> [f64; D] { weighted_polynomial_regression(x, y, map(y, |y| 1./sq(y))) }
fn polynomial_fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { polynomial_regression(map(t, x), map(t, y)) }

mod collision_integrals; // Reduced collision integrals table computed from Chapman-Enskog theory with Stockmayer potential by L. Monchick and E.A. Mason. Transport properties of polar gases. J. Chem. Phys.
use super::{kB, NA};
const light_speed : f64 = 299_792_458.;
const μ0 : f64 = 1.2566370621e-6;
const ε0 : f64 = 1./(light_speed*light_speed*μ0);
use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
use {iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}}};

//struct Transport<const S: usize> { D: [f64; S], η: f64, λ: f64 }
impl<const S: usize> super::System<S> where [(); S-1]: {
	pub fn transport(&self, _pressure: f64, temperature: f64, amounts: &[f64; S]) -> f64 {//Transport<S> {
		let (header_log_T⃰, [Ω⃰22,_A⃰ ,_B⃰,_C⃰]) = {use collision_integrals::*; (eval(header_T⃰.copied(), log),
			// Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
			 [&Ω⃰22,&A⃰ ,&B⃰,&C⃰].map(|table| table.map(|T⃰_row| polynomial_regression(header_δ⃰.copied(), T⃰_row)))
		)};
		let Self{molar_mass, thermodynamics: _, diameter, well_depth, polarizability, dipole, /*rotational_relaxation, internal_degrees_of_freedom,*/ ..} = self;
		let χ = |a, b| { // Corrections to the effective diameter and well depth to account for interaction between a polar and a non-polar molecule
			if dipole[a] == dipole[b] { 1. } else {
				let (polar, non_polar) = if dipole[a] != 0. { (a,b) } else { (b,a) };
				1. + polarizability[non_polar]/cb(diameter[non_polar]) * sq(dipole[polar]/sqrt(4.*π*ε0*cb(diameter[polar]) * well_depth[polar]))
						 * sqrt(well_depth[polar]/well_depth[non_polar]) / 4.
			}
		};
		let reduced_well_depth = |a, b| sqrt(well_depth[a]*well_depth[b]) * sq(χ(a, b));
		let T⃰ = |a, b, T| kB * T / reduced_well_depth(a, b);
		let reduced_dipole_moment = |a, b| dipole[a]*dipole[b] / (8. * π * ε0 * sqrt(well_depth[a]*well_depth[b]) * cb((diameter[a] + diameter[b])/2.)); // ̃δ⃰
		let collision_integral = |table : &[[f64; 8]; 39], a, b, T| {
			let log_T⃰ = log(T⃰ (a, a, T));
			assert!(*header_log_T⃰ .first().unwrap() <= log_T⃰  && log_T⃰  <= *header_log_T⃰ .last().unwrap());
			let i = header_log_T⃰ [..header_log_T⃰.len()-4].iter().position(|&header_log_T⃰ | header_log_T⃰  < log_T⃰).unwrap();
			let δ⃰ = reduced_dipole_moment(a, b);
			use std::convert::TryInto;
			let polynomials : &[_; 3] = &table[i..i+3].try_into().unwrap();
			quadratic_interpolation(header_log_T⃰[i..][..3].try_into().unwrap(), polynomials.map(|P| eval_poly(&P, δ⃰ )), log_T⃰)
		};
		let Ω⃰22 = |a, b, T| collision_integral(&Ω⃰22, a, b, T);
		let viscosity = |a, T| 5./16. * sqrt(π * molar_mass[a] * kB * T / NA) / (Ω⃰22(a, a, T) * π * sq(diameter[a]));
		/*let reduced_mass = |a, b| molar_mass[a] * molar_mass[b] / (NA * (molar_mass[a] + molar_mass[b])); // reduced_mass (a,a) = W/2NA
		let _Ω⃰11 = |a, b, T| Ω⃰22(a, b, T)/collision_integral(&A⃰, a, b, T);
		let conductivity = |a, T| {
			let self_diffusion_coefficient = 3./16. * sqrt(2.*π/reduced_mass(a,a)) * pow(kB*T, 3./2.) / (π * sq(diameter[a]) * Ω⃰11(a, a, T));
			let f_internal = molar_mass[a]/(kB*NA * T) * self_diffusion_coefficient / viscosity(a, T);
			let T⃰ = T⃰ (a, a, T);
			let fz_T⃰ = 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
			// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
			let fz_298 = (|T⃰| 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰)(kB * 298. / well_depth[a]);
			let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[a] * fz_298 / fz_T⃰ + 2./π * (5./3. * internal_degrees_of_freedom[a] + f_internal));
			let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[a]/(3./2.));
			let f_rotation = f_internal * (1. + c1);
			let Cv_internal = thermodynamics[a].specific_heat_capacity(T) - 5./2. - internal_degrees_of_freedom[a];
			(viscosity(a, T)/molar_mass[a])*kB*NA*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[a] + f_internal * Cv_internal)
		};
		let diffusion_coefficient = |a, b, T| {
			let reduced_diameter = (diameter[a] + diameter[b])/2. * pow(χ(a, b), -1./6.);
			3./16. * sqrt(2.*π/reduced_mass(a,b)) * pow(kB*T, 3./2.) / (π*sq(reduced_diameter)*Ω⃰11(a, b, T))
		};*/

		// polynomial fit in the temperature range for every specie (pairs)
		/*const*/let [temperature_min, temperature_max] : [f64; 2] = [300., 3000.];
		const D : usize = 4; // FIXME: Remez
		const N : usize = D+2; // FIXME: Remez
		let T = generate(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min));
		let sqrt_viscosity_sqrt_T_polynomials : [[_; D]; S] = generate(|a| polynomial_fit::<_,_,_,D,N>(T, log, |T| sqrt(viscosity(a,T)/sqrt(T)))).collect();
		//let conductivity_polynomials = generate(|a| polyfit(|T| conductivity(a, T)/sqrt(T)));
		//let diffusion_coefficient_polynomials = generate(|a| iter::box_collect((0..=a).map(|b| polyfit(|T| diffusion_coefficient(a,b,T)/pow(T, 3./2.)))));

		let T = temperature;
		let sqrt_viscosity = |a| sq(sqrt(T) * eval_poly(&sqrt_viscosity_sqrt_T_polynomials[a], log(T)));
		amounts.dot(generate(|k|
			sq(sqrt_viscosity(k)) /
			amounts.dot(generate(|j|
				sq(1. + (sqrt_viscosity(k)/sqrt_viscosity(j)) * sqrt(sqrt(molar_mass[j]/molar_mass[k]))) /
				(sqrt(8.) * sqrt(1. + molar_mass[k]/molar_mass[j])))
			)
		))
	}}
