//#![feature(trait_alias,once_cell,array_methods,array_map,in_band_lifetimes,bindings_after_at)]#![allow(uncommon_codepoints,non_upper_case_globals,non_snake_case)]
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
	if false {
		let A = DMatrix::from_iterator(N, D, (0..D).map(|k| x.iter().zip(w.iter()).map(move |(x, w)| w*x.powi(k as i32))).flatten());
		let b = DVector::from_iterator(N, y.into_iter().zip(w.iter()).map(|(y, w)| w*y));
		SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap()
	} else {
		let ref y = y.collect();
		#[link(name = "cantera")] extern "C" { fn cantera_polyfit(n: usize, deg: usize, x: *const f64 , y: *const f64, w: *const f64 , p: *mut f64 ) -> f64; }
		let mut p = [0.; D];
		unsafe{cantera_polyfit(N, D-1, x.as_ptr(), y.as_ptr(), w.as_ptr(), p.as_mut_ptr())};
		p
	}
}

// Regression with 1/y² weights (towards relative vertical error)
pub fn regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Clone) -> [f64; D] { weighted_regression(x, y.clone(), y.map(|y| 1./sq(y))) }
pub fn fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { regression(t.map(x), t.map(y)) }
}
use {std::{cmp::min, f64::consts::PI as π}, num::{sq, cb, pow}, iter::{IntoConstSizeIterator, ConstRange, into::IntoCopied, list, map}};

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

use std::cell::RefCell;
thread_local! {
	static debug_data: RefCell<f64> = Default::default();
	static debug_info: RefCell<String> = Default::default();
}

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
		debug_data.with(|d| d.replace((interpolation_start_index-I0) as f64));
		let header_ln_T⃰ : &[_; 3] = header_ln_T⃰[interpolation_start_index..][..3].try_into().unwrap();
		assert!(*header_ln_T⃰ .first().unwrap() <= ln_T⃰  && ln_T⃰  <= *header_ln_T⃰ .last().unwrap());
		let polynomials: &[_; 3] = &fit[interpolation_start_index-I0..][..3].try_into().unwrap();
		let δ⃰ = self.reduced_dipole_moment(a, b);
		assert!(*header_δ⃰ .first().unwrap() <= δ⃰  && δ⃰  <= *header_δ⃰ .last().unwrap(),"{a} {b} {δ⃰}");
		let table : [_; 3] = table[interpolation_start_index-I0..][..3].try_into().unwrap();
		let y = if δ⃰ == 0. { table.map(|row| row[0]) } else { polynomials.each_ref().map(|P| polynomial::evaluate(P, δ⃰ )) };
		debug_info.with(|d| d.replace(format!("{:?} {header_ln_T⃰:?} {y:?} {a} {b} {T} {} {ln_T⃰}", header_ln_T⃰.each_ref().map(|&x| f64::exp(x)), self.T⃰ (a, b, T))));
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
	let T : [_; N] = ConstRange.map(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min)).collect();
	for (n,&T) in T.iter().enumerate() { if T < 1900. { for k in 0..K { assert!(species.T⃰(k,k, T) <= 50., "{k} {n} {T}"); } } }
	if true {
		#[link(name = "cantera")] extern "C" { fn get_data(i: usize) -> f64; fn get_info(i: usize) -> *const i8; }
		let k = (0..K).map(|k| (0..1).map(move |_| (0..30).map(move |n| (({ // Cantera incorrectly clips its omega22 indexing (breaks for T⃰*>50(~temperature>1900))
			let T = T[n];
			//species.thermodynamics[k].molar_heat_capacity_at_constant_pressure_R(T)
			//species.T⃰ (k, j, T)
			//species.T⃰ (k, k, T)
			//species.reduced_dipole_moment(k, k)
			species.Ω⃰22(k, k, T);
			debug_data.with(|d| d.take())
			//species.Ω⃰11(k, j, T)
			//species.binary_thermal_diffusion(k, j, T)
			//species.binary_thermal_diffusion(k, j, T) / (T*f64::sqrt(T))
			//species.thermal_conductivity(k,T)
		}, debug_info.with(|d| d.take())),
			(unsafe{get_data(k*50+n)}, unsafe{std::ffi::CStr::from_ptr(get_info(k*50+n))})
			//unsafe{*get_global().offset(((k*53+j)*50+n) as isize)}
		))));
		let k = k.map(|j| j.map(|n| n.map(|(A@(a,_),B@(b,_))| (num::relative_error(a,b), (A,B)))));
		let k = k.enumerate().flat_map(|(k,j)| std::iter::repeat(k).zip(j.enumerate().flat_map(|(j,n)| std::iter::repeat(j).zip(n.enumerate()))));
		let (k,(j,(n,(e,((a,A),(b,B)))))) = k.reduce(|A@(_,(_,(_, (a, _)))), B@(_,(_,(_, (b, _))))| if a>b { A } else { B }).unwrap();
		assert!(e<1e-14, "{}: {k} {j} {n} {a} {b} {e:.0e} {}\n{A}\n{B:?}", T[n], species.T⃰(k,k,T[n]),

			 );
	}
	Self{
		thermal_conductivityIVT: map(0..K, |k| polynomial::fit(T, f64::ln, |T| species.thermal_conductivity(k,T)/f64::sqrt(T))),
		VviscosityIVVT: map(0..K, |k| polynomial::fit(T, f64::ln, |T| f64::sqrt(species.viscosity(k, T))/f64::sqrt(f64::sqrt(T)))),
		binary_thermal_diffusionITVT: list((0..K).map(|k| (0..K).map(move |j|
			polynomial::fit(T, f64::ln, |T| species.binary_thermal_diffusion(k, j, T) / (T*f64::sqrt(T))) )).flatten())
	}
}
}

use ast::*;

pub fn thermal_conductivityIVTI2<const D: usize>(thermal_conductivityIVT: &[[f64; D]], [lnT, lnT2, lnT3, lnT4]: &[Value; 4], mole_fractions: &[Value], f: &mut Block) -> Expression  {
	assert!(D == 5);
	let K = thermal_conductivityIVT.len();
	let ref thermal_conductivityIVT = map(thermal_conductivityIVT, |P| l!(f P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3 + P[4]*lnT4));
	      (0..K).map(|k| &mole_fractions[k] * &thermal_conductivityIVT[k]).sum::<Expression>() +
	1. / (0..K).map(|k| &mole_fractions[k] / &thermal_conductivityIVT[k]).sum::<Expression>()
}

pub fn viscosityIVT<const D: usize>(molar_mass: &[f64], VviscosityIVVT: &[[f64; D]], [lnT, lnT2, lnT3, lnT4]: &[Value; 4], mole_fractions: &[Value], f: &mut Block) -> Expression {
	assert!(D == 5);
	let K = VviscosityIVVT.len();
	let ref VviscosityIVVT = map(VviscosityIVVT, |P| l!(f P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3 + P[4]*lnT4));
	sum((0..K).map(|k|
		&mole_fractions[k] * sq(&VviscosityIVVT[k]) / sum((0..K).map(|j| {
			let Va = f64::sqrt(1./f64::sqrt(8.) * 1./f64::sqrt(1. + molar_mass[k]/molar_mass[j]));
			let mut sq = |x| { let ref x=l!(f x); x*x };
			&mole_fractions[j] * sq(Va + (Va*f64::sqrt(f64::sqrt(molar_mass[j]/molar_mass[k]))) * &VviscosityIVVT[k]/&VviscosityIVVT[j])
		}))
	))
}

pub fn PITVT_mixture_diffusion_coefficients<'t, const D: usize>(binary_thermal_diffusionITVT: &[[f64; D]], [lnT, lnT2, lnT3, lnT4]: &[Value; 4], mole_fractions: &'t [Value], mass_fractions: impl 't+IntoIterator<Item=Expression>, f: &mut Block) -> impl 't+Iterator<Item=Expression> {
	assert!(D == 5);
	let K = mole_fractions.len();
	/*{
		let T = 1000.;
		let lnT = f64::ln(T);
		let lnT2 = sq(lnT);
		let lnT3 = lnT2*lnT;
		let lnT4 = lnT2*lnT2;
		let ref binary_thermal_diffusion_coefficients_TVT = map(0..K, |k| map(0..k, |j|
			{let P = binary_thermal_diffusion_coefficients_TVT[k*K+j]; P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3+ P[4]*lnT4}
		));
		#[link(name = "cantera")] extern "C" { fn get_global() -> *const f64; }
		let ((k,j,n), (a,b,e)) = (0..K).map(move |k| (0..k).map(move |j| (0..1).map(move |n| ((k,j,n), {
			let (a, b) = (T*f64::sqrt(T)*binary_thermal_diffusion_coefficients_TVT[k][j], unsafe{*get_global().offset((k*9+j) as isize)});
			(a,b,num::relative_error(a,b))
		})))).flatten().flatten().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
		assert!(e<=1e-12, "{k} {j} {n} {a} {b} {e:.0e}");
	}*/
	let binary_thermal_diffusionITVT = map(0..K, |k| map(0..k, |j|
		{let P = binary_thermal_diffusionITVT[k*K+j]; l!(f P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3+ P[4]*lnT4)}
	));
	mass_fractions.into_iter().enumerate().map(move |(k, mass_fraction)| (1. - mass_fraction) / (0..K).filter(|&j| j != k).map(|j|
		&mole_fractions[j] / &binary_thermal_diffusionITVT[std::cmp::max(k,j)][min(k,j)]
	).sum::<Expression>())
}

pub fn properties_<const D: usize>(molar_mass: &[f64], polynomials: &Polynomials<D>) -> Function {
	let K = molar_mass.len();
	/*{
		let pressure = 101325.;
		let pressure_R = pressure/(kB*NA);
		let temperature : f64 = 1000.;
		let volume = 1.;
		let amount = pressure_R * volume / temperature;
		let amount_proportions = map(0..K, |_| 1.);
		let amounts = map(&*amount_proportions, |amount_proportion| amount * amount_proportion/amount_proportions.iter().sum::<f64>());
		let total_amount = amounts.iter().sum::<f64>();
		let rcp_amount = 1./total_amount;
		let nonbulk_amounts = &amounts[0..amounts.len()-1];
		let nonbulk_fractions= map(0..K-1, |k| rcp_amount*f64::max(0., nonbulk_amounts[k]));
		let bulk_fraction= 1. - nonbulk_fractions.iter().sum::<f64>();
		let mole_fractions = list(nonbulk_fractions.into_vec().into_iter().chain(std::iter::once(bulk_fraction)));
		std::dbg!(&mole_fractions);
		//let mean_molar_mass : f64 = iter::dot(molar_mass.iter().copied().zip(mole_fractions.iter().copied()));
		let T = temperature;
		let lnT = f64::ln(T);
		let lnT2 = sq(lnT);
		let lnT3 = lnT2*lnT;
		let lnT4 = lnT2*lnT2;
		let ref binary_thermal_diffusionITVT : Box<[Box<[f64]>]> = map(0..K, |k| map(0..k, |j|
			{let P = polynomials.binary_thermal_diffusionITVT[k*K+j]; P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3+ P[4]*lnT4}
		));
		#[link(name = "cantera")] extern "C" { fn get_global() -> *const f64; }
		//let k = (0..K).map(|k| (0..K).map(move |j| (0..5).map(move |n| (polynomials.binary_thermal_diffusionITVT[k*K+j][n], unsafe{*get_global().offset(((k*9+j)*5+n) as isize)}))));
		//let p = pressure;
		//let mmw = mean_molar_mass;
		//let ref m_molefracs = mole_fractions;
		//let ref m_mw = molar_mass;
		//let k = (0..K).map(|k| (0..1).map(move |_| (0..1).map(move |_| ((mmw - m_molefracs[k] * m_mw[k])*1000., unsafe{*get_global().offset(k as isize)}))));
		/*let k = (0..K).map(|k| (0..1).map(move |_| (0..1).map(move |_| ({
			let m_bdiff = |j:usize,k:usize| -> f64 { T*f64::sqrt(T)*binary_thermal_diffusionITVT[std::cmp::max(k,j)][min(k,j)] };
			let sum2 : f64 = (0..K).filter(|&j| j != k).map(|j| m_molefracs[j] / m_bdiff(j,k)).sum();
			sum2
			//p * mmw * sum2 * 1000.
			//(mmw - m_molefracs[k] * m_mw[k]) / (p * mmw * sum2)
		}, unsafe{*get_global().offset(k as isize)}))));*/
		let k = (0..K).map(|k| (0..K).map(move |j| (0..1).map(move |_| ({
			let m_bdiff = |j:usize,k:usize| -> f64 { T*f64::sqrt(T)*binary_thermal_diffusionITVT[std::cmp::max(k,j)][min(k,j)] };
			//if j != k { m_molefracs[j] / m_bdiff(j,k) } else { 0. }
			if j != k { m_bdiff(j,k) } else { 0. }
		}, unsafe{*get_global().offset((k*9+j) as isize)}))));
		let k = k.map(|j| j.map(|n| n.map(|(a,b)| (num::relative_error(a,b), (a,b)))));
		let k = k.enumerate().flat_map(|(k,j)| std::iter::repeat(k).zip(j.enumerate().flat_map(|(j,n)| std::iter::repeat(j).zip(n.enumerate()))));
		let (k,(j,(n,(e,(a,b))))) = k.reduce(|A@(_,(_,(_, (a, _)))), B@(_,(_,(_, (b, _))))| if a>b { A } else { B }).unwrap();
		assert!(e<=1e-6, "{k} {j} {n} {a} {b} {e:.0e}");
	}*/
	let_!{ input@[ref pressure_R, ref total_amount, ref T, ref nonbulk_amounts @ ..] = &*map(0..(3+K-1), Value) => {
	let mut values = ["pressure_","total_amount","T"].iter().map(|s| s.to_string()).chain((0..K-1).map(|i| format!("nonbulk_amounts[{i}]"))).collect();
	let mut function = Block::new(&mut values);
	let ref mut f = function;
	let lnT = l!(f ln(1024., T, f));
	let lnT2 = l!(f sq(&lnT));
	let lnT3 = l!(f &lnT2*&lnT);
	let lnT4 = l!(f &lnT2*&lnT2);
	let ref lnT = [lnT, lnT2, lnT3, lnT4];
	let ref rcp_amount = l!(f 1./total_amount);
	let nonbulk_fractions= map(0..K-1, |k| l!(f rcp_amount*max(0., &nonbulk_amounts[k])));
	let bulk_fraction= l!(f 1. - nonbulk_fractions.iter().sum::<Expression>());
	let ref mole_fractions = list(nonbulk_fractions.into_vec().into_iter().chain(std::iter::once(bulk_fraction)));
	let ref rcp_mean_molar_mass = l!(f 1./dot(&molar_mass, &mole_fractions));
	let mass_fractions =	mole_fractions.iter().zip(&*molar_mass).map(|(x,&m)| m * rcp_mean_molar_mass * x);
	//d[k] = (mmw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
	//d[k] = (1 - (m_molefracs[k] * m_mw[k])/mmw)/(p * sum2);
	//d[k] = [T√T/p] (1 - (m_molefracs[k] * m_mw[k])/mmw)/([T√T*]sum2);
	// sum2 = sum m_molefracs[j] / m_bdiff(j,k) [/T√T];
	// sum2 += [T√T*] m_molefracs[j] / m_bdiff(j,k);
	let ref VT = l!(f sqrt(T));
	let ref TVTIP = l!(f T*VT/(pressure_R*NA*kB));
	Function{output: list([
		(VT/2.)*thermal_conductivityIVTI2(&polynomials.thermal_conductivityIVT, lnT, mole_fractions, f),
		VT*viscosityIVT(molar_mass, &polynomials.VviscosityIVVT, lnT, mole_fractions, f),
	].into_iter().chain(
		PITVT_mixture_diffusion_coefficients(&polynomials.binary_thermal_diffusionITVT, lnT, mole_fractions, mass_fractions, f).map(|PITVTID| TVTIP*PITVTID)
	)), statements: function.statements.into(), input: input.len(), values: values.into()}
}}}

pub fn properties<const D: usize>(species: &Species) -> Function { properties_(&species.molar_mass, &Polynomials::<D>::new(&species)) }

pub fn properties_rust<const D: usize>(species@Species{molar_mass,..}: &Species) -> impl Fn(f64, f64, f64, &[f64]) -> (f64, f64, Box<[f64]>) +'_ {
	assert!(D == 5);
	let Polynomials{thermal_conductivityIVT, VviscosityIVVT, binary_thermal_diffusionITVT} = Polynomials::<D>::new(&species);
	let K = molar_mass.len();
	move |pressure_R: f64, total_amount: f64, T: f64, nonbulk_amounts: &[f64]| -> (f64, f64, Box<[f64]>) {
	let lnT = f64::ln(T);
	let lnT2 = sq(lnT);
	let lnT3 = lnT2*lnT;
	let lnT4 = lnT2*lnT2;
	let ref rcp_amount = 1./total_amount;
	let nonbulk_fractions= map(0..K-1, |k| rcp_amount*f64::max(0., nonbulk_amounts[k]));
	let bulk_fraction= 1. - nonbulk_fractions.iter().sum::<f64>();
	let ref mole_fractions = list(nonbulk_fractions.into_vec().into_iter().chain(std::iter::once(bulk_fraction)));
	let mean_molar_mass : f64 = iter::dot(molar_mass.iter().copied().zip(mole_fractions.iter().copied()));
	let ref rcp_mean_molar_mass = 1./mean_molar_mass;
	let mass_fractions =	mole_fractions.iter().zip(&**molar_mass).map(|(x,&m)| m * rcp_mean_molar_mass * x);
	(
		(f64::sqrt(T)/2.)*{
			let thermal_conductivityIVT = map(&*thermal_conductivityIVT, |P| P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3 + P[4]*lnT4);
	      (0..K).map(|k| mole_fractions[k] * thermal_conductivityIVT[k]).sum::<f64>() + 1. / (0..K).map(|k| mole_fractions[k] / thermal_conductivityIVT[k]).sum::<f64>()
		},
		f64::sqrt(T)*{
			let VviscosityIVVT = map(&*VviscosityIVVT, |P| P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3 + P[4]*lnT4);
			(0..K).map(|k|
				mole_fractions[k] * sq(VviscosityIVVT[k]) / (0..K).map(|j| {
					let Va = f64::sqrt(1./f64::sqrt(8.) * 1./f64::sqrt(1. + molar_mass[k]/molar_mass[j]));
					mole_fractions[j] * sq(Va + (Va*f64::sqrt(f64::sqrt(molar_mass[j]/molar_mass[k]))) * VviscosityIVVT[k]/VviscosityIVVT[j])
				}).sum::<f64>()
			).sum::<f64>()
		},
		{
			let binary_thermal_diffusionITVT = map(0..K, |k| map(0..k, |j| {let P = binary_thermal_diffusionITVT[k*K+j]; P[0] + P[1]*lnT + P[2]*lnT2 + P[3]*lnT3+ P[4]*lnT4}));
			map(mass_fractions.into_iter().enumerate(), move |(k, mass_fraction)| T*f64::sqrt(T)/(pressure_R*NA*kB) * (1. - mass_fraction) / (0..K).filter(|&j| j != k).map(|j|
				mole_fractions[j] / binary_thermal_diffusionITVT[std::cmp::max(k,j)][min(k,j)]
			).sum::<f64>())
		}
	)
}}
