//#![feature(trait_alias,once_cell,array_methods,array_map,in_band_lifetimes,bindings_after_at,iter_zip)]#![allow(uncommon_codepoints,non_upper_case_globals,non_snake_case)]
fn quadratic_interpolation(x: &[f64; 3], y: &[f64; 3], x0: f64) -> f64 {
	assert!(x[0] != x[1]); assert!(x[1] != x[2]); assert!(x[0] != x[2]);
	((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1]
}
mod polynomial {
//#![feature(trait_alias)]#![allow(non_snake_case)]
use {num::sq, iter::{IntoIterator, generate, IntoConstSizeIterator, DotN}};
pub fn evaluate<const N: usize>(P: &[f64; N], x: f64) -> f64 { P.dot(generate(|k| x.powi(k as i32))) }

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
use {std::{iter::zip, cmp::min, f64::consts::PI as π}, num::{sq, cb, pow}, iter::{IntoConstSizeIterator, Copied, list, map, DotN, Enumerate, Cloned, eval}};

use super::{light_speed, kB, NA};
const fine_structure : f64 = 7.2973525693e-3;
const Planck : f64 = 6.62607015e-34;
const electron_charge : f64 = 1.602176634e-19;
const μ0 : f64 = 2. * fine_structure * Planck / (electron_charge * electron_charge * light_speed); // H/m (Henry=kg⋅m²/(s²A²))
const ε0 : f64 = 1./(light_speed*light_speed*μ0); // F/m (Farad=s⁴A²/(m²kg)
mod collision_integrals; // Reduced collision integrals table computed from Chapman-Enskog theory with Stockmayer potential by L. Monchick and E.A. Mason. Transport properties of polar gases. J. Chem. Phys.
pub use collision_integrals::{header_T⃰, header_δ⃰};
// Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
fn polynomial_regression_δ⃰<const N: usize>(table: &[[f64; 8]; N]) -> [[f64; 6]; N] {
	// Cantera fits degree 6 but evaluates only up to 5 :/
	table.each_ref().map(|T⃰_row:&[_; 8]| { let poly6 : [_; 7] = polynomial::weighted_regression(header_δ⃰.copied(), T⃰_row.copied(), [1.; 8]); poly6[0..6].try_into().unwrap() })
}
use std::lazy::SyncLazy;
/*const*/static Ω⃰22: SyncLazy<[[f64; 6]; 37]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::Ω⃰22));
/*const*/static A⃰: SyncLazy<[[f64; 6]; 39]> = SyncLazy::new(|| polynomial_regression_δ⃰(&collision_integrals::A⃰));

use super::Species;
impl Species {
	fn	reduced_mass(&self, a: usize, b: usize) -> f64 {
		let Self{molar_mass, ..} = self;
		molar_mass[a]/NA * molar_mass[b]/NA / (molar_mass[a]/NA + molar_mass[b]/NA) // reduced_mass (a,a) = W/2NA
	}
	fn χ(&self, a: usize, b: usize) -> f64 { // Corrections to the effective diameter and well depth to account for interaction between a polar and a non-polar molecule
		let Self{diameter, well_depth_J, polarizability, permanent_dipole_moment, ..} = self;
		if (permanent_dipole_moment[a]>0.) == (permanent_dipole_moment[b]>0.) { 1. } else {
			let (polar, non_polar) = if permanent_dipole_moment[a] != 0. { (a,b) } else { (b,a) };
			1. + 1./4. * polarizability[non_polar]/cb(diameter[non_polar]) * sq(permanent_dipole_moment[polar]/f64::sqrt(4.*π*ε0*well_depth_J[polar]*cb(diameter[polar]))) * f64::sqrt(well_depth_J[polar]/well_depth_J[non_polar])
		}
	}
	fn interaction_well_depth(&self, a: usize, b: usize) -> f64 {
		let Self{well_depth_J, ..} = self;
		f64::sqrt(well_depth_J[a]*well_depth_J[b]) * sq(self.χ(a, b))
	}
	pub fn T⃰(&self, a: usize, b: usize, T: f64) -> f64 { T * kB / self.interaction_well_depth(a, b) }
	pub fn reduced_dipole_moment(&self, a: usize, b: usize) -> f64 {
		let Self{well_depth_J, permanent_dipole_moment, diameter, ..} = self;
		permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * π * ε0 * f64::sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.))
	}
	fn reduced_diameter(&self, a: usize, b: usize) -> f64 {
		let Self{diameter, ..} = self;
		(diameter[a] + diameter[b])/2. * pow(self.χ(a, b), -1./6.)
	}
	fn collision_integral<const I0: usize, const N: usize>(&self, table: &[[f64; 8]], fit: &[[f64; 6]; N], a: usize, b: usize, T: f64) -> f64 {
		let ln_T⃰ = f64::ln(self.T⃰ (a, b, T));
		/*const*/let header_ln_T⃰ = header_T⃰.each_ref().map(|&T| f64::ln(T));
		let interpolation_start_index = min((1+header_ln_T⃰ [1..header_ln_T⃰.len()].iter().position(|&header_ln_T⃰ | ln_T⃰ < header_ln_T⃰ ).unwrap())-1, I0+table.len()-3);
		let header_ln_T⃰ : &[_; 3] = header_ln_T⃰[interpolation_start_index..][..3].try_into().unwrap();
		assert!(*header_ln_T⃰ .first().unwrap() <= ln_T⃰  && ln_T⃰  <= *header_ln_T⃰ .last().unwrap());
		let polynomials: &[_; 3] = &fit[interpolation_start_index-I0..][..3].try_into().unwrap();
		let δ⃰ = self.reduced_dipole_moment(a, b);
		assert!(*header_δ⃰ .first().unwrap() <= δ⃰  && δ⃰  <= *header_δ⃰ .last().unwrap(),"{a} {b} {δ⃰}");
		let table : [_; 3] = table[interpolation_start_index-I0..][..3].try_into().unwrap();
		let y = if δ⃰ == 0. { table.map(|row| row[0]) } else { polynomials.each_ref().map(|P| polynomial::evaluate(P, δ⃰ )) };
		let image = quadratic_interpolation(header_ln_T⃰, &y, ln_T⃰);
		assert!(image > 0.);
		image
	}
	pub fn Ω⃰22(&self, a: usize, b: usize, T: f64) -> f64 { self.collision_integral::<1, 37>(&collision_integrals::Ω⃰22, &Ω⃰22, a, b, T) }
	pub fn viscosity(&self, k: usize, T: f64) -> f64 {
		let Self{molar_mass, diameter, ..} = self;
		5./16. * f64::sqrt(π * molar_mass[k]/NA * kB*T) / (self.Ω⃰22(k, k, T) * π * sq(diameter[k]))
	}
	fn Ω⃰11(&self, a: usize, b: usize, T: f64) -> f64 { self.Ω⃰22(a, b, T)/self.collision_integral::<0, 39>(&collision_integrals::A⃰, &A⃰, a, b, T) }
	fn binary_thermal_diffusion(&self, a: usize, b: usize, T: f64) -> f64 {
		3./16. * f64::sqrt(2.*π/self.reduced_mass(a,b)) * pow(kB*T, 3./2.) / (π*sq(self.reduced_diameter(a,b))*self.Ω⃰11(a, b, T))
	}
	fn thermal_conductivity(&self, k: usize, T: f64) -> f64 {
		let Self{molar_mass, thermodynamics, rotational_relaxation, internal_degrees_of_freedom, ..} = self;
		let f_internal = molar_mass[k]/NA/(kB * T) * self.binary_thermal_diffusion(k,k,T) / self.viscosity(k, T);
		let fz = |T⃰| 1. + pow(π, 3./2.) / f64::sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
		// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
		let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[k] * fz(self.T⃰(k,k, 298.)) / fz(self.T⃰(k,k, T)) + 2./π * (5./3. * internal_degrees_of_freedom[k] + f_internal));
		let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[k]/(3./2.));
		let f_rotation = f_internal * (1. + c1);
		let Cv_internal = thermodynamics[k].molar_heat_capacity_at_constant_pressure_R(T) - 5./2. - internal_degrees_of_freedom[k];
		(self.viscosity(k, T)/(molar_mass[k]/NA))*kB*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[k] + f_internal * Cv_internal)
	}
}

pub struct Polynomials<const D: usize> {
	pub thermal_conductivityIVT: Box<[[f64; D]]>,
	pub VviscosityIVVT: Box<[[f64; D]]>,
	pub binary_thermal_diffusionITVT: Box<[[f64; D]]>
}
impl<const D: usize> Polynomials<D> {
pub fn new(species: &Species) -> Self {
	let K = species.len();
	let [temperature_min, temperature_max] : [f64; 2] = [300., 3000.];
	const N : usize = /*D+2 FIXME: Remez*/50;
	let T : [_; N] = eval(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min));
	for (n,&T) in T.iter().enumerate() { if T < 1900. { for k in 0..K { assert!(species.T⃰(k,k, T) <= 50., "{k} {n} {T}"); } } }
	Self{
		thermal_conductivityIVT: map(0..K, |k| polynomial::fit(T, f64::ln, |T| species.thermal_conductivity(k,T)/f64::sqrt(T))),
		VviscosityIVVT: map(0..K, |k| polynomial::fit(T, f64::ln, |T| f64::sqrt(species.viscosity(k, T))/f64::sqrt(f64::sqrt(T)))),
		binary_thermal_diffusionITVT: list((0..K).map(|k| (0..K).map(move |j|
			polynomial::fit(T, f64::ln, |T| species.binary_thermal_diffusion(k, j, T) / (T*f64::sqrt(T))) )).flatten())
	}
}
}

use ast::*;

pub fn thermal_conductivityIVT<const D: usize>(thermal_conductivityIVT: &[[f64; D]], lnT: &[Expr; D], mole_fractions: &[Value], f: &mut Block) -> Expression  {
	// zip(mole_fractions, thermal_conductivityIVT).map(|(X, P)| { let y=l!(f P.dot(lnT.cloned()):Expression); (X*y, X/y) }).reduce(|(A,B),(a,b)| (l!(f;A+a),l!(f;B+b))).unwrap();
	let [mut A, mut B]:[Option<Expression>;2] = [None,None];
	for (X, P) in zip(mole_fractions, thermal_conductivityIVT) {
		let y = l!(f P.dot(lnT.cloned()):Expression);
		let [a,b] = [X*y, X/y];
		A = Some(if let Some(A) = A { l!(f; A+a) } else { a });
		B = Some(if let Some(B) = B { l!(f; B+b) } else { b });
	}
	A.unwrap() + 1./B.unwrap()
}

pub fn viscosityIVT<const D: usize>(molar_mass: &[f64], VviscosityIVVT: &[[f64; D]], lnT: &[Expr; D], mole_fractions: &[Value], f: &mut Block) -> Expression {
	let K = VviscosityIVVT.len();
	let VviscosityIVVT = map(VviscosityIVVT, |P| l!(f P.dot(lnT.cloned()):Expression));
	let ref rcp_VviscosityIVVT = map(&*VviscosityIVVT, |x| l!(f 1./x));
	sum((0..K).map(|k|
		&mole_fractions[k] * sq(&VviscosityIVVT[k]) / sum((0..K).map(|j| {
			let Va = f64::sqrt(1./f64::sqrt(8.) * 1./f64::sqrt(1. + molar_mass[k]/molar_mass[j]));
			let mut sq = |x| { let ref x=l!(f x); x*x };
			&mole_fractions[j] * sq(Va + (Va*f64::sqrt(f64::sqrt(molar_mass[j]/molar_mass[k]))) * VviscosityIVVT[k] * rcp_VviscosityIVVT[j])
		}))
	))
}

pub fn PITVT_mixture_diffusion<'t, const D: usize>(binary_thermal_diffusionITVT: &[[f64; D]], lnT: &[Expr; D], mole_fractions: &'t [Value], mass_fractions: impl 't+IntoIterator<Item=Expression>, f: &mut Block) -> impl 't+Iterator<Item=Expression> {
	let K = mole_fractions.len();
	let binary_thermal_diffusionITVT = map(0..K, |k| map(0..k, |j|
		{let P = binary_thermal_diffusionITVT[k*K+j]; l!(f P.dot(lnT.cloned()):Expression)}
	));
	mass_fractions.into_iter().enumerate().map(move |(k, mass_fraction)| (1. - mass_fraction) / (0..K).filter(|&j| j != k).map(|j|
		&mole_fractions[j] / &binary_thermal_diffusionITVT[std::cmp::max(k,j)][min(k,j)]
	).sum::<Expression>())
}

pub fn properties_<const D: usize>(molar_mass: &[f64], Polynomials{thermal_conductivityIVT, VviscosityIVVT, binary_thermal_diffusionITVT} : &Polynomials<D>,
	temperature0: f64, Vviscosity: f64, thermal_conductivity: f64, density_mixture_diffusion: f64) -> Function {
	// (x+y)^n = (0..=n).sum(|k| binomial(n,k)*x^k*y^(n-k)) = (0..=n).sum(|k| binomial(n,k)*x^(n-k)*y^k)
	// (lnT')^n = ln(T0*T)^n = (lnT0+lnT)^n = (0..=n).sum(|k| binomial(n,k)*lnT0^(n-k)*lnT^k)
	let scale = |u, P:&[_; D]| {
		fn factorial(n: usize) -> usize { assert!(n<5); (1..=n).product() }
		fn binomial(n: usize, k: usize) -> usize { assert!(k<=n); factorial(n) / (factorial(n - k) * factorial(k)) }
		let lnT0 = |n| f64::powi(f64::ln(temperature0), n as i32);
		P.enumerate().map(|(k,P)| (k..D).map(|n| binomial(n,k) as f64*lnT0(n-k)).sum::<f64>()*P/u).collect()
	};
	let VviscosityIVVT = map(&**VviscosityIVVT, |P| scale(Vviscosity, P));
	let thermal_conductivityIVT = map(&**thermal_conductivityIVT, |P| scale(2.*thermal_conductivity, P));
	let binary_thermal_diffusionITVT = map(&**binary_thermal_diffusionITVT, |P| scale(density_mixture_diffusion/(kB*NA), P));

	let K = molar_mass.len();
	let_!{ input@[ref total_amount, ref temperature, ref nonbulk_amounts @ ..] = &*map(0..(3+K-1), Value) => {
	let mut values = ["pressure_","total_amount","temperature"].iter().map(|s| s.to_string()).chain((0..K-1).map(|i| format!("nonbulk_amounts[{i}]"))).collect();
	let mut function = Block::new(&mut values);
	let ref mut f = function;
	let T = temperature;
	let lnT = l!(f ln(1024., T, f));
	fn replace_with<T, F: FnOnce(T) -> T>(x: &mut T, f: F) {unsafe{std::ptr::write(x, f(std::ptr::read(x)))}}
	//let ref lnT = scan((1.).into(), |x| { let y = x.clone(); replace_with(x, |x| (x * lnT).expr()); l!(f; y) }); // Would need a scan(||->T) i.e without early return and thus with impl ExactSize
	let ref lnT = {let mut x:Expr=(1.).into(); eval(|_| { let y = x.clone(); replace_with(&mut x, |x| (x * lnT).expr()); l!(f; y) })};
	let ref rcp_total_amount = l!(f 1./total_amount);
	let nonbulk_fractions= map(0..K-1, |k| l!(f rcp_total_amount*max(0., &nonbulk_amounts[k])));
	let bulk_fraction= l!(f 1. - nonbulk_fractions.iter().sum::<Expression>());
	let ref mole_fractions = list(nonbulk_fractions.into_vec().into_iter().chain(std::iter::once(bulk_fraction)));
	use iter::Dot;
	let ref mean_molar_mass = l!(f molar_mass.copied().dot(mole_fractions):Expression);
	let ref rcp_mean_molar_mass = l!(f 1./mean_molar_mass);
	let mass_fractions =	mole_fractions.iter().zip(&*molar_mass).map(|(x,&m)| m * rcp_mean_molar_mass * x);
	let ref VT = l!(f sqrt(T));
	let ref density_TVTIP = l!(f mean_molar_mass * VT);
	Function{
		output: list([
			VT*self::thermal_conductivityIVT(&thermal_conductivityIVT, lnT, mole_fractions, f),
			VT*viscosityIVT(molar_mass, &VviscosityIVVT, lnT, mole_fractions, f),
		].into_iter().chain(
			PITVT_mixture_diffusion(&binary_thermal_diffusionITVT, lnT, mole_fractions, mass_fractions, f).map(|PITVTID| density_TVTIP*PITVTID)
		)),
		statements: function.statements.into(),
		input: vec![Type::F64; input.len()].into(),
		values: values.into()
	}
}}}

pub fn properties<const D: usize>(species: &Species, temperature: f64, sqrt_viscosity: f64, thermal_conductivity: f64, density_mixture_diffusion_coefficient: f64) -> Function {
	properties_(&species.molar_mass, &Polynomials::<D>::new(&species), temperature, sqrt_viscosity, thermal_conductivity, density_mixture_diffusion_coefficient)
}
