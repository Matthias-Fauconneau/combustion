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
	assert!(!A.iter().any(|a| a.is_nan()), "{}", A);
	SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap()
}

use iter::into::IntoMap;
// Regression with 1/y² weights (towards relative vertical error)
fn polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Clone) -> [f64; D] { weighted_polynomial_regression(x, y.clone(), y.map(|y| 1./sq(y))) }
fn polynomial_fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { polynomial_regression(t.map(x), t.map(y)) }

use {std::{convert::TryInto, cmp::min, f64::consts::PI as π}, num::{sq, cb, sqrt, log, pow}};
use {iter::{from_iter, eval, zip, dot, into::IntoCopied, vec::{self, Dot, generate}}};
const light_speed : f64 = 299_792_458.;
const μ0 : f64 = 1.2566370621e-6; //  H/m (Henry=kg⋅m²/(s²A²))
const ε0 : f64 = 1./(light_speed*light_speed*μ0); // F/m (Farad=s⁴A²/(m²kg)
use super::{K, NA, Species, State};
mod collision_integrals; // Reduced collision integrals table computed from Chapman-Enskog theory with Stockmayer potential by L. Monchick and E.A. Mason. Transport properties of polar gases. J. Chem. Phys.
pub use collision_integrals::{header_T⃰, header_δ⃰};
// Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
fn polynomial_regression_δ⃰(table: &[[f64; 8]; 39]) -> [[f64; 7]; 39] { table.each_ref().map(|T⃰_row:&[_; 8]| polynomial_regression(header_δ⃰.copied(), T⃰_row.copied())) }
use std::lazy::SyncLazy;
/*const*/pub static Ω⃰22: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::Ω⃰22));
/*const*/pub static A⃰: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::A⃰));
/*const*/pub static B⃰: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::B⃰));
/*const*/pub static C⃰: SyncLazy<[[f64; 7]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::C⃰));

const D : usize = 5;
#[derive(Debug)] pub struct TransportPolynomials {
	pub sqrt_viscosity_T14: Box<[[f64; D]]>,
	pub thermal_conductivity_T12: Box<[[f64; D]]>,
	pub binary_thermal_diffusion_coefficients_T32: Box<[Box<[[f64; D]]>]>
}
impl Species {
	fn	reduced_mass(&self, a: usize, b: usize) -> f64 {
		let Self{molar_mass, ..} = self;
		molar_mass[a]/NA * molar_mass[b]/NA / (molar_mass[a]/NA + molar_mass[b]/NA) // reduced_mass (a,a) = W/2NA
	}

	fn χ(&self, a: usize, b: usize) -> f64 { // Corrections to the effective diameter and well depth to account for interaction between a polar and a non-polar molecule
		let Self{diameter, well_depth_J, polarizability, permanent_dipole_moment, ..} = self;
		if (permanent_dipole_moment[a]>0.) == (permanent_dipole_moment[b]>0.) { 1. } else {
			let (polar, non_polar) = if permanent_dipole_moment[a] != 0. { (a,b) } else { (b,a) };
			1. + 1./4. * polarizability[non_polar]/cb(diameter[non_polar]) * permanent_dipole_moment[polar]/sqrt(well_depth_J[polar]*cb(diameter[polar])) * sqrt(well_depth_J[polar]/well_depth_J[non_polar])
		}
	}
	fn interaction_well_depth(&self, a: usize, b: usize) -> f64 {
		let Self{well_depth_J, ..} = self;
		sqrt(well_depth_J[a]*well_depth_J[b]) * sq(self.χ(a, b))
	}
	fn T⃰(&self, a: usize, b: usize, T: f64) -> f64 { T * K / self.interaction_well_depth(a, b) }
	fn reduced_dipole_moment(&self, a: usize, b: usize) -> f64 { // ̃δ⃰ Cantera
		let Self{well_depth_J, permanent_dipole_moment, diameter, ..} = self;
		permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * π * ε0 * sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.))
	}
	fn reduced_diameter(&self, a: usize, b: usize) -> f64 {
		let Self{diameter, ..} = self;
		(diameter[a] + diameter[b])/2. * pow(self.χ(a, b), -1./6.)
	}
	fn collision_integral(&self, table: &[[f64; 7]; 39], a: usize, b: usize, T: f64) -> f64 {
		let log_T⃰ = log(self.T⃰ (a, a, T));
		let δ⃰ = self.reduced_dipole_moment(a, b);
		/*const*/let header_log_T⃰ = vec::eval(header_T⃰.copied(), log);
		let interpolation_start_index = min((1+header_log_T⃰ [1..header_log_T⃰.len()].iter().position(|&header_log_T⃰ | log_T⃰ < header_log_T⃰ ).unwrap())-1, header_log_T⃰.len()-3);
		let header_log_T⃰ : &[_; 3] = header_log_T⃰[interpolation_start_index..][..3].try_into().unwrap();
		let polynomials: &[_; 3] = &table[interpolation_start_index..][..3].try_into().unwrap();
		let image = quadratic_interpolation(header_log_T⃰, &from_iter(polynomials.map(|P| eval_poly(P, δ⃰ ))), log_T⃰);
		assert!(*header_log_T⃰ .first().unwrap() <= log_T⃰  && log_T⃰  <= *header_log_T⃰ .last().unwrap());
		assert!(*header_δ⃰ .first().unwrap() <= δ⃰  && δ⃰  <= *header_δ⃰ .last().unwrap());
		assert!(image > 0.);
		image
	}
	fn Ω⃰22(&self, a: usize, b: usize, T: f64) -> f64 { self.collision_integral(&Ω⃰22, a, b, T) }
	fn viscosity(&self, a: usize, T: f64) -> f64 {
		let Self{molar_mass, diameter, ..} = self;
		5./16. * sqrt(π * molar_mass[a]/NA * K*T) / (self.Ω⃰22(a, a, T) * π * sq(diameter[a]))
	}
	fn Ω⃰11(&self, a: usize, b: usize, T: f64) -> f64 { self.Ω⃰22(a, b, T)/self.collision_integral(&A⃰, a, b, T) }
	fn thermal_conductivity(&self, a: usize, T: f64) -> f64 {
		let Self{molar_mass, thermodynamics, diameter, well_depth_J, rotational_relaxation, internal_degrees_of_freedom, ..} = self;
		let self_diffusion_coefficient = 3./16. * sqrt(2.*π/self.reduced_mass(a,a)) * pow(T, 3./2.) / (π * sq(diameter[a]) * self.Ω⃰11(a, a, T));
		let f_internal = molar_mass[a]/NA/T * self_diffusion_coefficient / self.viscosity(a, T);
		let T⃰ = self.T⃰ (a, a, T);
		let fz = |T⃰| 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
		// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
		let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[a] * fz(298.*K / well_depth_J[a]) / fz(T⃰) + 2./π * (5./3. * internal_degrees_of_freedom[a] + f_internal));
		let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[a]/(3./2.));
		let f_rotation = f_internal * (1. + c1);
		let Cv_internal = thermodynamics[a].specific_heat_capacity(T) - 5./2. - internal_degrees_of_freedom[a];
		(self.viscosity(a, T)/(molar_mass[a]/NA))*K*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[a] + f_internal * Cv_internal)
	}
	fn binary_thermal_diffusion_coefficient(&self, a: usize, b: usize, T: f64) -> f64 {
		3./16. * sqrt(2.*π/self.reduced_mass(a,b)) * pow(T, 3./2.) / (π*sq(self.reduced_diameter(a,b))*self.Ω⃰11(a, b, T))
	}
	pub fn transport_polynomials(&self) -> TransportPolynomials {
		/*const*/let [temperature_min, temperature_max] : [f64; 2] = [300., /*3000.*/1000.]; //|*K
		const N : usize = /*D+2 FIXME: Remez*/50;
		let T : [_; N] = generate(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min)).collect();
		TransportPolynomials{
			sqrt_viscosity_T14: eval(self.len(), |a| polynomial_fit::<_,_,_,D,N>(T, log, |T| sqrt(self.viscosity(a, T))/sqrt(sqrt(T)))),
			thermal_conductivity_T12: eval(self.len(), |a| polynomial_fit::<_,_,_,D,N>(T, log, |T| self.thermal_conductivity(a,T)/sqrt(T))),
			binary_thermal_diffusion_coefficients_T32: eval(self.len(), |a| eval(self.len(), |b| polynomial_fit::<_,_,_,D,N>(T, log, |T| self.binary_thermal_diffusion_coefficient(a,b,T)/pow(T,3./2.))))
		}
	}
}

impl TransportPolynomials {
	fn sqrt_viscosity(&self, a: usize, T: f64) -> f64 { sqrt(sqrt(T)) * eval_poly(&self.sqrt_viscosity_T14[a], log(T)) }
	fn thermal_conductivity(&self, a: usize, T: f64) -> f64 {  sqrt(T) * eval_poly(&self.thermal_conductivity_T12[a], log(T)) }
	fn binary_thermal_diffusion_coefficient(&self, a: usize, b: usize, T: f64) -> f64 { pow(T,3./2.) * eval_poly(&self.binary_thermal_diffusion_coefficients_T32[if a>b {a} else {b}][if a>b {b} else {a}], log(T)) }
}

#[derive(Debug)] pub struct Transport { pub viscosity: f64, pub thermal_conductivity: f64, pub mixture_averaged_thermal_diffusion_coefficients: Box<[f64]> }
pub fn transport(molar_mass: &[f64], transport_polynomials: &TransportPolynomials, State{temperature, pressure, volume, amounts}: &State) -> Transport {
	let T = *temperature;
	let viscosity = dot(zip(amounts.copied(), |k|
		sq(transport_polynomials.sqrt_viscosity(k, T)) /
		dot(zip(amounts.copied(), |j|
			sq(1. + transport_polynomials.sqrt_viscosity(k, T)/transport_polynomials.sqrt_viscosity(j, T) * sqrt(sqrt(molar_mass[j]/molar_mass[k]))) /
			(sqrt(8.) * sqrt(1. + molar_mass[k]/molar_mass[j]))
		))
	));
	let amount = pressure * volume / (K * T);
	{let e = f64::abs(amounts.iter().sum::<f64>()-amount)/amount; assert!(e < 2e-16, "{} {} {:e}", amounts.iter().sum::<f64>(), amount, e);}
	let thermal_conductivity = 1./2. * (
		dot(zip(amounts.copied(), |k| transport_polynomials.thermal_conductivity(k, T))) / amount +
		amount / dot(zip(amounts.copied(), |k| 1. / transport_polynomials.thermal_conductivity(k, T)))
	);
	let mixture_averaged_thermal_diffusion_coefficients = eval(amounts.len(), |k|
		(1. - amounts[k]/amount) /
		dot(zip(amounts.copied(), |j| if j != k { 1. / transport_polynomials.binary_thermal_diffusion_coefficient(k, j, T) } else { 0. }))
	);
	Transport{viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients}
}

#[cfg(test)] mod test;

/*pub trait AbsError {
	fn error(&self, o: &Self) -> f64;
}

impl AbsError for f64 {
	fn error(&self, o: &Self) -> f64 {
		f64::abs(self-o)
	}
}

impl<const N: usize> AbsError for [f64; N] {
	fn error(&self, o: &Self) -> f64 {
		self.iter().zip(o).map(|(s,o)| s.error(o)).reduce(f64::max).unwrap()
	}
}

impl AbsError for Box<[f64]> {
	fn error(&self, o: &Self) -> f64 {
		self.iter().zip(o.iter()).map(|(s,o)| s.error(o)).reduce(f64::max).unwrap()
	}
}

impl AbsError for Transport {
	fn error(&self, o: &Self) -> f64 {
		[
			self.viscosity.error(&o.viscosity),
			self.thermal_conductivity.error(&o.thermal_conductivity),
			self.mixture_averaged_thermal_diffusion_coefficients.error(&o.mixture_averaged_thermal_diffusion_coefficients)
		].iter().copied().reduce(f64::max).unwrap()
	}
}*/

pub trait RelError {
	fn error(&self, o: &Self) -> f64;
}

impl RelError for f64 {
	fn error(&self, o: &Self) -> f64 {
		num::relative_error(*self, *o)
	}
}

impl<const N: usize> RelError for [f64; N] {
	fn error(&self, o: &Self) -> f64 {
		self.iter().zip(o).map(|(s,o)| s.error(o)).reduce(f64::max).unwrap()
	}
}

impl RelError for Box<[f64]> {
	fn error(&self, o: &Self) -> f64 {
		self.iter().zip(o.iter()).map(|(s,o)| s.error(o)).reduce(f64::max).unwrap()
	}
}

impl RelError for Transport {
	fn error(&self, o: &Self) -> f64 {
		[
			self.viscosity.error(&o.viscosity),
			self.thermal_conductivity.error(&o.thermal_conductivity),
			self.mixture_averaged_thermal_diffusion_coefficients.error(&o.mixture_averaged_thermal_diffusion_coefficients)
		].iter().copied().reduce(f64::max).unwrap()
	}
}
