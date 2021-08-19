// -> iter
// In map(|..| f(&mut context)).reduce(|..| g(&mut context)): borrow checker cannot detect context is borrowed exclusively.
// => map_reduce(|context,..| f(context), |context,..| g(context))
trait MapReduce : Iterator {
fn map_reduce<C, B, M: FnMut(&mut C, Self::Item) -> B, R: FnMut(&mut C, B, B) -> B>(self, context: &mut C, map: M, reduce: R) -> Option<B>;
}
impl<I:Iterator> MapReduce for I {
fn map_reduce<C, B, M: FnMut(&mut C, Self::Item) -> B, R: FnMut(&mut C, B, B) -> B>(mut self, context: &mut C, mut map: M, mut reduce: R) -> Option<B>  {
	let first = map(context, self.next()?);
	Some(self.fold(first, |acc,x| { let x = map(context, x); reduce(context, acc, x) }))
}
}

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
use {std::{iter::zip, cmp::min, f64::consts::PI as π}, num::{sq, cb, sqrt, ln, pow}, iter::{Copied, list, map, DotN, Cloned, eval}};

use super::{light_speed, kB, NA, R};
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
			1. + 1./4. * polarizability[non_polar]/cb(diameter[non_polar]) * sq(permanent_dipole_moment[polar]/sqrt(4.*π*ε0*well_depth_J[polar]*cb(diameter[polar]))) * sqrt(well_depth_J[polar]/well_depth_J[non_polar])
		}
	}
	fn interaction_well_depth(&self, a: usize, b: usize) -> f64 {
		let Self{well_depth_J, ..} = self;
		sqrt(well_depth_J[a]*well_depth_J[b]) * sq(self.χ(a, b))
	}
	pub fn T⃰(&self, a: usize, b: usize, T: f64) -> f64 { T * kB / self.interaction_well_depth(a, b) }
	pub fn reduced_dipole_moment(&self, a: usize, b: usize) -> f64 {
		let Self{well_depth_J, permanent_dipole_moment, diameter, ..} = self;
		permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * π * ε0 * sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.))
	}
	fn reduced_diameter(&self, a: usize, b: usize) -> f64 {
		let Self{diameter, ..} = self;
		(diameter[a] + diameter[b])/2. * pow(self.χ(a, b), -1./6.)
	}
	fn collision_integral<const I0: usize, const N: usize>(&self, table: &[[f64; 8]], fit: &[[f64; 6]; N], a: usize, b: usize, T: f64) -> f64 {
		let ln_T⃰ = ln(self.T⃰ (a, b, T));
		/*const*/let header_ln_T⃰ = header_T⃰.each_ref().map(|&T| ln(T));
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
		5./16. * sqrt(π * molar_mass[k]/NA * kB*T) / (self.Ω⃰22(k, k, T) * π * sq(diameter[k]))
	}
	fn Ω⃰11(&self, a: usize, b: usize, T: f64) -> f64 { self.Ω⃰22(a, b, T)/self.collision_integral::<0, 39>(&collision_integrals::A⃰, &A⃰, a, b, T) }
	fn diffusivity(&self, a: usize, b: usize, T: f64) -> f64 {
		3./16. * sqrt(2.*π/self.reduced_mass(a,b)) * pow(kB*T, 3./2.) / (π*sq(self.reduced_diameter(a,b))*self.Ω⃰11(a, b, T))
	}
	fn conductivity(&self, k: usize, T: f64) -> f64 {
		let Self{molar_mass, thermodynamics, rotational_relaxation, internal_degrees_of_freedom, ..} = self;
		let f_internal = molar_mass[k]/NA/(kB * T) * self.diffusivity(k,k,T) / self.viscosity(k, T);
		let fz = |T⃰| 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
		// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
		let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[k] * fz(self.T⃰(k,k, 298.)) / fz(self.T⃰(k,k, T)) + 2./π * (5./3. * internal_degrees_of_freedom[k] + f_internal));
		let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[k]/(3./2.));
		let f_rotation = f_internal * (1. + c1);
		let Cv_internal = thermodynamics[k].molar_heat_capacity_at_constant_pressure_R(T) - 5./2. - internal_degrees_of_freedom[k];
		(self.viscosity(k, T)/(molar_mass[k]/NA))*kB*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[k] + f_internal * Cv_internal)
	}
}

pub struct Polynomials<const D: usize> {
	pub conductivityIVT: Box<[[f64; D]]>,
	pub VviscosityIVVT: Box<[[f64; D]]>,
	pub diffusivityITVT: Box<[[f64; D]]>
}
impl<const D: usize> Polynomials<D> {
pub fn new(species: &Species, T0: f64) -> Self {
	let K = species.len();
	let [temperature_min, temperature_max] : [f64; 2] = [300., 3000.];
	const N : usize = /*D+2 FIXME: Remez*/50;
	let T : [_; N] = eval(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min));
	//use itertools::Itertools; println!("[{}]", (0..K).format_with(",\n",|i, f| f(&format_args!("[{}]", (0..K).format_with(", ",|j, f| f(&format_args!("[{:e}]", T.iter().map(|&T| species.diffusivity(i,j, T)).format(", "))))))));
	Self{
		conductivityIVT: map(0..K, |k| polynomial::fit(T, |T| ln(T/T0), |T| species.conductivity(k,T)/sqrt(T))),
		VviscosityIVVT: map(0..K, |k| polynomial::fit(T, |T| ln(T/T0), |T| sqrt(species.viscosity(k, T)/sqrt(T)))),
		diffusivityITVT: list((0..K).map(|k| (0..K).map(move |j| polynomial::fit(T, |T| ln(T/T0), |T| species.diffusivity(k, j, T) / (T*sqrt(T))) )).flatten())
	}
}
}

use ast::*;

pub fn conductivityNIVT<const D: usize>(conductivityIVT: &[[f64; D]], total_amount: &Value, lnT: &[Expr; D], mole_proportions: &[Value], f: &mut Block) -> Expression  {
	let (_,[A,B]) = zip(mole_proportions, conductivityIVT).enumerate().map_reduce(f,
		|f,(k,(X, P))| { let c = f.def(P.dot(lnT.cloned()):Expression, format!("c{k}")); (k,[X*c, X/c]) },
		|f,(_,[A,B]),(k,[a,b])| (k,[f.def(A+a, format!("a{k}")).into(), f.def(B+b, format!("b{k}")).into()]) ).unwrap();
	A/total_amount + total_amount/B
}

pub fn viscosityIVT<const D: usize>(molar_mass: &[f64], VviscosityIVVT: &[[f64; D]], lnT: &[Expr; D], mole_fractions: &[Value], f: &mut Block) -> Expression {
	let K = VviscosityIVVT.len();
	let VviscosityIVVT = map(VviscosityIVVT.iter().enumerate(), |(k,P)| f.def(P.dot(lnT.cloned()):Expression, format!("VviscosityIVVT{k}")));
	let rcp_VviscosityIVVT = map(VviscosityIVVT.iter().enumerate(), |(k,x)| f.def(1./x, format!("rcp_VviscosityIVVT{k}")));
	sum((0..K).map(|k|
		&mole_fractions[k] * ast::sq(&VviscosityIVVT[k]) / sum((0..K).map(|j| {
			let Va = sqrt(1./sqrt(8.) * 1./sqrt(1. + molar_mass[k]/molar_mass[j]));
			&mole_fractions[j] * ast::sq(Va + (Va*sqrt(sqrt(molar_mass[j]/molar_mass[k]))) * VviscosityIVVT[k] * rcp_VviscosityIVVT[j])
		}))
	))
}

fn replace_with<T, F: FnOnce(T) -> T>(x: &mut T, f: F) {unsafe{std::ptr::write(x, f(std::ptr::read(x)))}}

pub fn density_diffusivity<'t, const D: usize>(molar_mass: &'t [f64], diffusivityITVT: &[[f64; D]], mean_molar_mass_VTN: &'t Value, VT: &'t Value, lnT: &[Expr; D], mole_proportions: &'t [Value], f: &mut Block) -> impl 't+Iterator<Item=Expression> {
	let K = mole_proportions.len();
	let rcp_diffusivityITVT = |f:&mut Block, k,j| {
		assert!(j<k);
		let P:[f64;D] = diffusivityITVT[k*K+j];
		f.def(1./(P.dot(lnT.cloned()):Expression), format!("R{k}_{j}"))
	};
	let mut S: Box<[Option<Value>]> = vec![None; K].into();
	for k in 0..K { for j in 0..k {
		let rcp_diffusivityITVT = rcp_diffusivityITVT(f, k,j);
		replace_with(&mut S[k], |S| {let t = mole_proportions[j]*rcp_diffusivityITVT; Some(f.def(if let Some(S) = S { S+t } else { t }, format!("S{k}_{j}")))});
		replace_with(&mut S[j], |S| {let t = mole_proportions[k]*rcp_diffusivityITVT; Some(f.def(if let Some(S) = S { S+t } else { t }, format!("S{j}_{k}")))});
	}}
	//let rcp_diffusivityITVT = map(0..K, |k| map(0..k, |j| rcp_diffusivityITVT(k,j)));
	S.into_vec().into_iter().enumerate().map(move |(k, S)|
		(mean_molar_mass_VTN - VT * 	molar_mass[k] * mole_proportions[k])
	/ S.unwrap()//(0..K).filter(|&j| j != k).map(|j| mole_proportions[j] * rcp_diffusivityITVT[std::cmp::max(k,j)][std::cmp::min(k,j)]).sum::<Expression>()
	)
}

pub fn properties_<const D: usize>(molar_mass: &[f64], Polynomials{conductivityIVT, VviscosityIVVT, diffusivityITVT} : &Polynomials<D>, temperature0: f64, viscosity: f64, conductivity: f64) -> Function {
	let VviscosityIVVT = map(&**VviscosityIVVT, |P| P.map(|p| (sqrt(sqrt(temperature0))/sqrt(viscosity))*p));
	let conductivityIVT = map(&**conductivityIVT, |P| P.map(|p| (sqrt(temperature0)/(2.*conductivity))*p));
	let diffusivityITVT = map(&**diffusivityITVT, |P| P.map(|p| (sqrt(temperature0)/(R*viscosity))*p));

	let K = molar_mass.len();
	let_!{ input@[ref sum_mole_proportions, ref temperature, ref nonbulk_mole_proportions @ ..] = &*map(0..(2+K-1), Value) => {
	let mut values = ["sum_mole_proportions","temperature"].iter().map(|s| s.to_string()).chain((0..K-1).map(|i| format!("nonbulk_mole_proportions[{i}]"))).collect::<Vec<_>>();
	assert!(input.len() == values.len());
	let mut function = Block::new(&mut values);
	let ref mut f = function;
	let T = temperature;
	let lnT = l!(f ast::ln(1024., T, f));
	//let ref lnT = scan((1.).into(), |x| { let y = x.clone(); replace_with(x, |x| (x * lnT).expr()); l!(f, y) }); // Would need a scan(||->T) i.e no early return i.e impl ExactSize
	let ref lnT = {let mut x:Expr=(1.).into(); eval(|_| { let y = x.clone(); replace_with(&mut x, |x| l!(f, x * lnT).expr()); y })};
	let ref mole_proportions = list(nonbulk_mole_proportions.iter().copied().chain([f.def(sum_mole_proportions-nonbulk_mole_proportions.iter().sum::<Expression>(), format!("mole_proportions{}",K-1))]));
	use iter::Dot;
	let ref VT = l!(f ast::sqrt(T));
	let ref mean_molar_massN = l!(f molar_mass.copied().dot(mole_proportions):Expression);
	let ref mean_molar_mass_VTN = f.def(mean_molar_massN * VT, "mean_molar_mass_VTN");
	Function{
		output: list([
			self::conductivityNIVT(&conductivityIVT, sum_mole_proportions, lnT, mole_proportions, f)*VT,
			VT*viscosityIVT(molar_mass, &VviscosityIVVT, lnT, mole_proportions, f),
		].into_iter().chain(
			density_diffusivity(molar_mass, &diffusivityITVT, mean_molar_mass_VTN, VT, lnT, mole_proportions, f)
		)),
		statements: function.statements.into(),
		input: vec![Type::F64; input.len()].into(),
		values: values.into()
	}
}}}

pub fn properties<const D: usize>(species: &Species, temperature: f64, viscosity: f64, conductivity: f64) -> Function {
	properties_(&species.molar_mass, &Polynomials::<D>::new(&species, temperature), temperature, viscosity, conductivity)
}
