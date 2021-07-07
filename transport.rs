#![feature(trait_alias,once_cell,array_methods,array_map,in_band_lifetimes,bindings_after_at)]#![allow(uncommon_codepoints,non_upper_case_globals,non_snake_case)]
fn quadratic_interpolation(x: &[f64; 3], y: &[f64; 3], x0: f64) -> f64 {
	assert!(x[0] != x[1]); assert!(x[1] != x[2]); assert!(x[0] != x[2]);
	((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1]
}
mod polynomial {
//#![feature(trait_alias)]#![allow(non_snake_case)]
use {num::sq, iter::{IntoIterator, ConstRange, IntoConstSizeIterator, Dot}};
pub fn evaluate<const N: usize>(P: &[f64; N], x: f64) -> f64 { Dot::dot(P, ConstRange.map(|k| x.powi(k as i32))) }

pub trait Vector<const N: usize> = IntoConstSizeIterator<N>+IntoIterator<Item=f64>;
pub fn weighted_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>, w: impl Vector<N>) -> [f64; D] {
	use nalgebra::{DMatrix, DVector, SVD};
	let ref w = DVector::from_iterator(N, w.into_iter());
	let ref x = x.collect();
	let A = DMatrix::from_iterator(N, D, (0..D).map(|k| x.iter().zip(w.iter()).map(move |(x, w)| w*x.powi(k as i32))).flatten());
	let b = DVector::from_iterator(N, y.into_iter().zip(w.iter()).map(|(y, w)| w*y));
	SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap()
}

// Regression with 1/y² weights (towards relative vertical error)
pub fn regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Clone) -> [f64; D] { weighted_regression(x, y.clone(), y.map(|y| 1./sq(y))) }
pub fn fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { regression(t.map(x), t.map(y)) }
}
use {std::{cmp::min, f64::consts::PI as π}, num::{sq, cb, sqrt, log2, ln, pow}, iter::{IntoConstSizeIterator, ConstRange, into::IntoCopied, list, map}};

const light_speed : f64 = 299_792_458.;
const μ0 : f64 = 1.2566370621e-6; //  H/m (Henry=kg⋅m²/(s²A²))
const ε0 : f64 = 1./(light_speed*light_speed*μ0); // F/m (Farad=s⁴A²/(m²kg)
mod collision_integrals; // Reduced collision integrals table computed from Chapman-Enskog theory with Stockmayer potential by L. Monchick and E.A. Mason. Transport properties of polar gases. J. Chem. Phys.
pub use collision_integrals::{header_T⃰, header_δ⃰};
// Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
fn polynomial_regression_δ⃰(table: &[[f64; 8]; 39]) -> [[f64; 7]; 39] { table.each_ref().map(|T⃰_row:&[_; 8]| polynomial::regression(header_δ⃰.copied(), T⃰_row.copied())) }
use std::lazy::SyncLazy;
/*const*/static Ω⃰22: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::Ω⃰22));
/*const*/static A⃰: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::A⃰));
//*const*/static B⃰: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::B⃰));
//*const*/static C⃰: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::C⃰));

use super::{kB, NA, Species};
struct Interactions<'t>(&'t Species);
impl Interactions<'t> {
	fn	reduced_mass(&self, a: usize, b: usize) -> f64 {
		let Species{molar_mass, ..} = self.0;
		molar_mass[a]/NA * molar_mass[b]/NA / (molar_mass[a]/NA + molar_mass[b]/NA) // reduced_mass (a,a) = W/2NA
	}
	fn χ(&self, a: usize, b: usize) -> f64 { // Corrections to the effective diameter and well depth to account for interaction between a polar and a non-polar molecule
		let Species{diameter, well_depth_J, polarizability, permanent_dipole_moment, ..} = self.0;
		if (permanent_dipole_moment[a]>0.) == (permanent_dipole_moment[b]>0.) { 1. } else {
			let (polar, non_polar) = if permanent_dipole_moment[a] != 0. { (a,b) } else { (b,a) };
			1. + 1./4. * polarizability[non_polar]/cb(diameter[non_polar]) * permanent_dipole_moment[polar]/sqrt(well_depth_J[polar]*cb(diameter[polar])) * sqrt(well_depth_J[polar]/well_depth_J[non_polar])
		}
	}
	fn interaction_well_depth(&self, a: usize, b: usize) -> f64 {
		let Species{well_depth_J, ..} = &self.0;
		sqrt(well_depth_J[a]*well_depth_J[b]) * sq(self.χ(a, b))
	}
	pub fn T⃰(&self, a: usize, b: usize, T: f64) -> f64 { T * kB / self.interaction_well_depth(a, b) }
	pub fn reduced_dipole_moment(&self, a: usize, b: usize) -> f64 {
		let Species{well_depth_J, permanent_dipole_moment, diameter, ..} = self.0;
		permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * π * ε0 * sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.))
	}
	fn reduced_diameter(&self, a: usize, b: usize) -> f64 {
		let Species{diameter, ..} = self.0;
		(diameter[a] + diameter[b])/2. * pow(self.χ(a, b), -1./6.)
	}
	fn collision_integral(&self, table: &[[f64; 7]; 39], a: usize, b: usize, T: f64) -> f64 {
		let ln_T⃰ = ln(self.T⃰ (a, b, T));
		/*const*/let header_ln_T⃰ = header_T⃰.each_ref().map(|&T| ln(T));
		let interpolation_start_index = min((1+header_ln_T⃰ [1..header_ln_T⃰.len()].iter().position(|&header_ln_T⃰ | ln_T⃰ < header_ln_T⃰ ).unwrap())-1, header_ln_T⃰.len()-3);
		let header_ln_T⃰ : &[_; 3] = header_ln_T⃰[interpolation_start_index..][..3].try_into().unwrap();
		let polynomials: &[_; 3] = &table[interpolation_start_index..][..3].try_into().unwrap();
		let δ⃰ = self.reduced_dipole_moment(a, b);
		let image = quadratic_interpolation(header_ln_T⃰, &polynomials.each_ref().map(|P| polynomial::evaluate(P, δ⃰ )), ln_T⃰);
		assert!(*header_ln_T⃰ .first().unwrap() <= ln_T⃰  && ln_T⃰  <= *header_ln_T⃰ .last().unwrap());
		assert!(*header_δ⃰ .first().unwrap() <= δ⃰  && δ⃰  <= *header_δ⃰ .last().unwrap());
		assert!(image > 0.);
		image
	}
	pub fn Ω⃰22(&self, a: usize, b: usize, T: f64) -> f64 { self.collision_integral(&Ω⃰22, a, b, T) }
	pub fn viscosity(&self, k: usize, T: f64) -> f64 {
		let Species{molar_mass, diameter, ..} = &self.0;
		5./16. * sqrt(π * molar_mass[k]/NA * kB*T) / (self.Ω⃰22(k, k, T) * π * sq(diameter[k]))
	}
	fn Ω⃰11(&self, a: usize, b: usize, T: f64) -> f64 { self.Ω⃰22(a, b, T)/self.collision_integral(&A⃰, a, b, T) }
	fn binary_thermal_diffusion_coefficient(&self, a: usize, b: usize, T: f64) -> f64 {
		3./16. * sqrt(2.*π/self.reduced_mass(a,b)) * pow(kB*T, 3./2.) / (π*sq(self.reduced_diameter(a,b))*self.Ω⃰11(a, b, T))
	}
	fn thermal_conductivity(&self, k: usize, T: f64) -> f64 {
		let Species{molar_mass, thermodynamics, well_depth_J, rotational_relaxation, internal_degrees_of_freedom, ..} = &self.0;
		let f_internal = molar_mass[k]/NA/(kB * T) * self.binary_thermal_diffusion_coefficient(k,k,T) / self.viscosity(k, T);
		let T⃰ = self.T⃰ (k, k, T);
		let fz = |T⃰| 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
		// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
		let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[k] * fz(298.*kB / well_depth_J[k]) / fz(T⃰) + 2./π * (5./3. * internal_degrees_of_freedom[k] + f_internal));
		let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[k]/(3./2.));
		let f_rotation = f_internal * (1. + c1);
		let Cv_internal = thermodynamics[k].molar_heat_capacity_at_constant_pressure_R(T) - 5./2. - internal_degrees_of_freedom[k];
		(self.viscosity(k, T)/(molar_mass[k]/NA))*kB*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[k] + f_internal * Cv_internal)
	}
}

pub struct Polynomials<const D: usize> {
	pub sqrt_viscosity_T_14: Box<[[f64; D]]>,
	pub thermal_conductivity_T_12: Box<[[f64; D]]>,
	pub binary_thermal_diffusion_coefficients_T32: Box<[[f64; D]]>
}
impl<const D: usize> Polynomials<D> {
pub fn new(species: &Species) -> Self {
	let K = species.len();
	let [temperature_min, temperature_max] : [f64; 2] = [300., 3500.];
	const N : usize = /*D+2 FIXME: Remez*/50;
	let T : [_; N] = ConstRange.map(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min)).collect();
	let ref s = Interactions(species);
	Self{
		sqrt_viscosity_T_14: map(0..K, |k| polynomial::fit(T, log2, |T| sqrt(s.viscosity(k, T))/sqrt(sqrt(T)))),
		thermal_conductivity_T_12: map(0..K, |k| polynomial::fit(T, log2, |T| s.thermal_conductivity(k,T)/sqrt(T))),
		binary_thermal_diffusion_coefficients_T32: map(0..K, map(|k| (0..K).map(move |j|
			polynomial::fit(T, log2, |T| s.binary_thermal_diffusion_coefficient(k, j, T) / (T*sqrt(T)))
		)).flatten())
	}
}
}

use ast::*;

pub fn thermal_conductivity_T_12_2<const D: usize>(thermal_conductivity_T_12: &[[f64; D]], [log_T, log_T_2, log_T_3]: &[Value; 3], mole_fractions: &[Value], f: &mut Block) -> Expression  {
	let K = thermal_conductivity_T_12.len();
	let ref thermal_conductivity_T_12 = map(thermal_conductivity_T_12, |P| f.def(P[0] + P[1]*log_T + P[2]*log_T_2 + P[3]*log_T_3));
	      (0..K).map(|k| &mole_fractions[k] * &thermal_conductivity_T_12[k]).sum::<Expression>() +
	1. / (0..K).map(|k| &mole_fractions[k] / &thermal_conductivity_T_12[k]).sum::<Expression>()
}

pub fn viscosity_T_12<const D: usize>(molar_mass: &[f64], sqrt_viscosity_T_14: &[[f64; D]], [log_T, log_T_2, log_T_3]: &[Value; 3], mole_fractions: &[Value], f: &mut Block) -> Expression {
	let K = sqrt_viscosity_T_14.len();
	let ref sqrt_viscosity_T_14 = map(sqrt_viscosity_T_14, |P| f.def(P[0] + P[1]*log_T + P[2]*log_T_2 + P[3]*log_T_3));
	sum((0..K).map(|k|
		&mole_fractions[k] * sq(&sqrt_viscosity_T_14[k]) / sum((0..K).map(|j| {
			let sqrt_a = sqrt(1./sqrt(8.) * 1./sqrt(1. + molar_mass[k]/molar_mass[j]));
			let mut sq = |x| { let x=f.def(x); &x*&x };
			&mole_fractions[j] * sq(sqrt_a + (sqrt_a*sqrt(sqrt(molar_mass[j]/molar_mass[k]))) * &sqrt_viscosity_T_14[k]/&sqrt_viscosity_T_14[j])
		}))
	))
}

pub fn P_T_32_mixture_diffusion_coefficients<'t, const D: usize>(binary_thermal_diffusion_coefficients_T32: &[[f64; D]], [log_T, log_T_2, log_T_3]: &[Value; 3], mole_fractions: &'t [Value], mass_fractions: impl 't+IntoIterator<Item=Expression>, f: &mut Block) -> impl 't+Iterator<Item=Expression> {
	let K = mole_fractions.len();
	let binary_thermal_diffusion_coefficients_T32 = map(0..K, |k| map(0..K, |j|
		if k>j { Some({let P = binary_thermal_diffusion_coefficients_T32[k*K+j]; f.def(P[0] + P[1]*log_T + P[2]*log_T_2 + P[3]*log_T_3)}) } else { None }
	));
	mass_fractions.into_iter().enumerate().map(move |(k, mass_fraction)| (1. - mass_fraction) / (0..K).filter(|&j| j != k).map(|j|
		&mole_fractions[j] / binary_thermal_diffusion_coefficients_T32[if k>j {k} else {j}][if k>j {j} else{k}].as_ref().unwrap()
	).sum::<Expression>())
}

pub fn properties_<const D: usize>(molar_mass: &[f64], polynomials: &Polynomials<D>) -> Function {
	let K = molar_mass.len();
	let_!{ input@[ref pressure_R, ref total_amount, ref T, ref active_amounts @ ..] = &*map(0..(3+K-1), Value) => {
	let mut function = FunctionBuilder::new(&input);
	let mut f = Block::new(&mut function);
	let log_T = f.def(log2(r#use(T)));
	let log_T_2 = f.def(sq(&log_T));
	let log_T_3 = f.def(&log_T_2*&log_T);
	let ref log_T = [log_T, log_T_2, log_T_3];
	let ref rcp_amount = f.def(1./total_amount);
	let active_fractions= map(0..K-1, |k| f.def(rcp_amount*max(0., &active_amounts[k])));
	let inert_fraction= f.def(1. - active_fractions.iter().sum::<Expression>());
	let ref mole_fractions = list(active_fractions.into_vec().into_iter().chain(std::iter::once(inert_fraction)));
	let ref rcp_mean_molar_mass = f.def(1./dot(&molar_mass, &mole_fractions));
	let mass_fractions =	mole_fractions.iter().zip(&*molar_mass).map(|(x,&m)| m * rcp_mean_molar_mass * x);
	let ref T_12 = f.def(sqrt(r#use(T)));
	let push= |v,f: &mut Block| f.push(v);
	push(output(0, T_12*viscosity_T_12(molar_mass, &polynomials.sqrt_viscosity_T_14, log_T, mole_fractions, &mut f)), &mut f);
	push(output(1, (T_12/2.)*thermal_conductivity_T_12_2(&polynomials.thermal_conductivity_T_12, log_T, mole_fractions, &mut f)), &mut f);
	let ref T_32_P = f.def(T*T_12/(pressure_R*NA*kB));
	for (k, P_T_32_D) in P_T_32_mixture_diffusion_coefficients(&polynomials.binary_thermal_diffusion_coefficients_T32, log_T, mole_fractions, mass_fractions, &mut f).enumerate() {
		f.push(output(2+k, T_32_P*P_T_32_D))
	}
	Function::new(2+K, f.into(), function)
}}}

pub fn properties<const D: usize>(species: &Species) -> Function { properties_(&species.molar_mass, &Polynomials::<D>::new(&species)) }
