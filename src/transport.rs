//fn linear_interpolation(x: &[f64; 2], y: &[f64; 2], x0: f64) -> f64 { assert!(x[0] != x[1]); y[0] + (y[1]-y[0])/(x[1]-x[0])*(x0-x[0]) }
fn quadratic_interpolation(x: &[f64; 3], y: &[f64; 3], x0: f64) -> f64 {
	assert!(x[0] != x[1]); assert!(x[1] != x[2]); assert!(x[0] != x[2]);
	((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1]
}

fn eval_poly<const N: usize>(P: &[f64; N], x: f64) -> f64 { P.dot(generate(|k| x.powi(k as i32))) }

trait Vector<const N: usize> = vec::Vector<N>+iter::IntoIterator<Item=f64>;
 fn weighted_polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>, w: impl Vector<N>) -> [f64; D] {
	use nalgebra::{DMatrix, DVector, SVD};
	let ref w = DVector::from_iterator(N, w.into_iter());
	//let A = DMatrix::from_iterator(N, D, x.into_iter().zip(w.iter()).map(|(x, w)| (0..D).map(move |k| w*x.powi(k as i32))).flatten());
	let ref x = x.collect();
	let A = DMatrix::from_iterator(N, D, (0..D).map(|k| x.iter().zip(w.iter()).map(move |(x, w)| w*x.powi(k as i32))).flatten());
	let b = DVector::from_iterator(N, y.into_iter().zip(w.iter()).map(|(y, w)| w*y));
	use std::convert::TryInto;
	assert!(!A.iter().any(|a| a.is_nan()), "{}", A);
	SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap()
}
// Regression with 1/y² weights (towards relative vertical error)
fn polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Copy) -> [f64; D] { weighted_polynomial_regression(x, y, map(y, |y| 1./sq(y))) }
fn polynomial_fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { polynomial_regression(map(t, x), map(t, y)) }

mod collision_integrals; // Reduced collision integrals table computed from Chapman-Enskog theory with Stockmayer potential by L. Monchick and E.A. Mason. Transport properties of polar gases. J. Chem. Phys.
use super::{kB, NA};
/*const light_speed : f64 = 299_792_458.;
const μ0 : f64 = 1.2566370621e-6; //  H/m (Henry=kg⋅m²/(s²A²))
const ε0 : f64 = 1./(light_speed*light_speed*μ0); // F/m (Farad=s⁴A²/(m²kg)*/
use {std::{cmp::min, f64::consts::PI as π}, num::{sq, cb, sqrt, log, pow, powi}};
use {iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}}};
use super::System;

#[derive(Debug, Clone, Copy)] pub struct Transport<const S: usize> {pub viscosity: f64, pub thermal_conductivity: f64, pub mixture_averaged_thermal_diffusion_coefficients: [f64; S] }
impl<const S: usize> System<S> where [(); S-1]: {
	pub fn transport(&self, pressure_R: f64, temperature: f64, amounts: &[f64; S]) -> Transport<S> {
		let header_log_T⃰ = eval(collision_integrals::header_T⃰.copied(), log);
		// Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
		let [Ω⃰22,A⃰ ,_B⃰,_C⃰] = {use collision_integrals::*; [&Ω⃰22,&A⃰ ,&B⃰,&C⃰].map(|table| table.map(|T⃰_row| polynomial_regression(header_δ⃰.copied(), T⃰_row)))};
		let Self{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, ..} = self;
		let χ = |a, b| { // Corrections to the effective diameter and well depth to account for interaction between a polar and a non-polar molecule
			if (permanent_dipole_moment[a]>0.) == (permanent_dipole_moment[b]>0.) { 1. } else {
				let (polar, non_polar) = if permanent_dipole_moment[a] != 0. { (a,b) } else { (b,a) };
				1. + 1./4. * polarizability[non_polar]/cb(diameter[non_polar]) * permanent_dipole_moment[polar]/sqrt(well_depth_J[polar]*cb(diameter[polar])) * sqrt(well_depth_J[polar]/well_depth_J[non_polar])
			}
		};
		let interaction_well_depth = |a, b| sqrt(well_depth_J[a]*well_depth_J[b]) * sq(χ(a, b));
		let T⃰ = |a, b, T| kB * T / interaction_well_depth(a, b);
		//let reduced_dipole_moment = |a, b| permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * π * ε0 * sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.)); // ̃δ⃰
		let reduced_diameter = |a,b| (diameter[a] + diameter[b])/2. * pow(χ(a, b), -1./6.);
		let reduced_dipole_moment = |a, b| 1./4. * permanent_dipole_moment[a]*permanent_dipole_moment[b] / (interaction_well_depth(a,b) * cb(reduced_diameter(a,b))) * χ(a,b); // ̃δ⃰
		let collision_integral = |table : &[[f64; /*8 FIXME*/1]; 39], a, b, T| {
			let log_T⃰ = log(T⃰ (a, a, T));
			//assert!(*header_log_T⃰ .first().unwrap() <= log_T⃰  && log_T⃰  <= *header_log_T⃰ .last().unwrap(), "{} {} {} {} {}", header_log_T⃰ .first().unwrap(), log_T⃰ ,  *header_log_T⃰ .last().unwrap(), T, T⃰ (a, a, T));
			use std::convert::TryInto;
			let δ⃰ = reduced_dipole_moment(a, b);
			/*let i = min((1+header_log_T⃰ [1..header_log_T⃰.len()].iter().position(|&header_log_T⃰ | log_T⃰ < header_log_T⃰ ).unwrap())-1, header_log_T⃰.len()-/*3*/2);
			let header_log_T⃰ : &[_; /*3*/2] = header_log_T⃰[i..][../*3*/2].try_into().unwrap();
			let polynomials: &[_; /*3*/2] = &table[i..][../*3*/2].try_into().unwrap();
			let image = linear_interpolation(header_log_T⃰, &polynomials.map(|P| eval_poly(&P, δ⃰ )), log_T⃰);*/
			let i = min((1+header_log_T⃰ [1..header_log_T⃰.len()].iter().position(|&header_log_T⃰ | log_T⃰ < header_log_T⃰ ).unwrap())-1, header_log_T⃰.len()-3);
			let header_log_T⃰ : &[_; 3] = header_log_T⃰[i..][..3].try_into().unwrap();
			let polynomials: &[_; 3] = &table[i..][..3].try_into().unwrap();
			let image = quadratic_interpolation(header_log_T⃰, &polynomials.map(|P| eval_poly(&P, δ⃰ )), log_T⃰);
			assert!(*header_log_T⃰ .first().unwrap() <= log_T⃰  && log_T⃰  <= *header_log_T⃰ .last().unwrap(), "{} {} {} {} {}", header_log_T⃰ .first().unwrap(), log_T⃰ ,  header_log_T⃰ .last().unwrap(), T, T⃰ (a, a, T));
			assert!(*collision_integrals::header_δ⃰ .first().unwrap() <= δ⃰  && δ⃰  <= *collision_integrals::header_δ⃰ .last().unwrap(), "{} {} {}", collision_integrals::header_δ⃰ .first().unwrap(), δ⃰ , collision_integrals::header_δ⃰ .last().unwrap());
			assert!(image > 0., "{:?}", (image, log_T⃰, header_log_T⃰, polynomials, δ⃰));
			image
		};
		let Ω⃰22 = |a, b, T| collision_integral(&Ω⃰22, a, b, T);
		let viscosity = |a, T| 5./16. * sqrt(π * molar_mass[a]/NA * kB * T) / (Ω⃰22(a, a, T) * π * sq(diameter[a]));
		//amounts.enumerate().filter(|(_, &n)| n > 0.).for_each(|(a, _)| { dbg!(viscosity(a,temperature)); } );
		let reduced_mass = |a, b| molar_mass[a]/NA * molar_mass[b]/NA / (molar_mass[a]/NA + molar_mass[b]/NA); // reduced_mass (a,a) = W/2NA
		let Ω⃰11 = |a, b, T| Ω⃰22(a, b, T)/collision_integral(&A⃰, a, b, T);
		let thermal_conductivity = |a: usize, T| {
			let self_diffusion_coefficient = 3./16. * sqrt(2.*π/reduced_mass(a,a)) * pow(kB*T, 3./2.) / (π * sq(diameter[a]) * Ω⃰11(a, a, T));
			let f_internal = molar_mass[a]/NA/(kB*T) * self_diffusion_coefficient / viscosity(a, T);
			let T⃰ = T⃰ (a, a, T);
			let fz_T⃰ = 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
			// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
			let fz_298 = (|T⃰| 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰)(kB * 298. / well_depth_J[a]);
			let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[a] * fz_298 / fz_T⃰ + 2./π * (5./3. * internal_degrees_of_freedom[a] + f_internal));
			let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[a]/(3./2.));
			let f_rotation = f_internal * (1. + c1);
			let Cv_internal = thermodynamics[a].specific_heat_capacity(T) - 5./2. - internal_degrees_of_freedom[a];
			(viscosity(a, T)/(molar_mass[a]/NA))*kB*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[a] + f_internal * Cv_internal)
		};
		//amounts.enumerate().filter(|(_, &n)| n > 0.).for_each(|(a, _)| { dbg!(thermal_conductivity(a,temperature)); } );
		let binary_thermal_diffusion_coefficient = |a, b, T| {
			3./16. * sqrt(2.*π/reduced_mass(a,b)) * pow(kB*T, 3./2.) / (π*sq(reduced_diameter(a,b))*Ω⃰11(a, b, T))
		};

		// polynomial fit in the temperature range for every specie (pairs)
		//for a in 0..S {use collision_integrals::*; dbg!(header_T⃰[header_T⃰.len()-2] *  well_depth_J[a] / kB); }
		/*const*/let [temperature_min, temperature_max] : [f64; 2] = [300., /*3000.*/1000.];
		const D : usize = 4;
		const N : usize = /*D+2 FIXME: Remez*/50;
		let T = generate(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min));
		let sqrt_viscosity_T14_polynomials : [[_; D]; S] = generate(|a| polynomial_fit::<_,_,_,D,N>(T, log, |T| sqrt(viscosity(a,T))/sqrt(sqrt(T)))).collect();
		let thermal_conductivity_sqrt_T_polynomials : [[_; D]; S] = generate(|a| polynomial_fit::<_,_,_,D,N>(T, log, |T| thermal_conductivity(a,T)/sqrt(T))).collect();
		let binary_thermal_diffusion_coefficient_polynomials : [Box<[[f64; D]]>; S] = generate(|a| iter::box_collect((0..=a).map(|b| polynomial_fit::<_,_,_,D,N>(T, log, |T| binary_thermal_diffusion_coefficient(a,b,T)/pow(T,3./2.))))).collect();
		let sqrt_viscosity = |a, T| sqrt(sqrt(T)) * eval_poly(&sqrt_viscosity_T14_polynomials[a], log(T));
		//let sqrt_viscosity = |a| sqrt(viscosity(a, T));
		let thermal_conductivity = |a, T| sqrt(T) * eval_poly(&thermal_conductivity_sqrt_T_polynomials[a], log(T));
		let binary_thermal_diffusion_coefficient = |a:usize, b:usize, T:f64| pow(T,3./2.) * eval_poly(&binary_thermal_diffusion_coefficient_polynomials[if a>b {a} else { b }][if a>b {b} else {a}], log(T));
		let T = temperature;
		let viscosity = amounts.dot(generate(|k|
			sq(sqrt_viscosity(k, T)) /
			amounts.dot(generate(|j|
				sq(1. + (sqrt_viscosity(k, T)/sqrt_viscosity(j, T)) * sqrt(sqrt(molar_mass[j]/molar_mass[k]))) /
				(sqrt(8.) * sqrt(1. + molar_mass[k]/molar_mass[j])))
			))
		);
		let amount = pressure_R / T * System::<S>::volume;
		assert_eq!(amounts.iter().sum::<f64>(), amount);
		let thermal_conductivity = 1./2. * (amounts.dot(generate(|k| thermal_conductivity(k, T))) / amount + amount / amounts.dot(generate(|k| 1. / thermal_conductivity(k, T))));
		let mixture_averaged_thermal_diffusion_coefficients = generate(|k| (1. - amounts[k]/amount) / amounts.dot(generate(|j| if j != k { 1. / binary_thermal_diffusion_coefficient(k, j, T) } else { 0. }))).collect();
		Transport{viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients}
	}}
