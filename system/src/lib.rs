//#![feature(const_generics, const_generics_defaults, const_evaluatable_checked, type_ascription, array_methods, in_band_lifetimes, once_cell, array_map, map_into_keys_values, bindings_after_at, destructuring_assignment, trait_alias, non_ascii_idents)]
#![feature(non_ascii_idents)]
//#![allow(incomplete_features, non_snake_case, confusable_idents, mixed_script_confusables, non_upper_case_globals, unused_imports, uncommon_codepoints)]
#![allow(non_upper_case_globals, non_snake_case)]
//use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
use num::log;
//use iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}};
//use model::{/*Map, Element,*/ Troe};
//mod transport; pub use transport::{TransportPolynomials, Transport};
//pub mod reaction; pub use reaction::{Reaction, Model, RateConstant};*/

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;

#[derive(PartialEq, Debug, /*Eq*/)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325. / NA; // 1 atm
	pub const T_split : f64 = 1000.*K;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T } // /RT
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R
	//fn dT_specific_heat_capacity(&self, T: f64) -> f64 { 	let a = self.a(T); kB*NA * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T) } // /R
	//fn dT_Gibbs_free_energy(&self, T: f64) -> f64 { let a = self.a(T); (1.-a[0])/T - a[1]/2. - a[2]/12.*T - a[3]/36.*T*T - a[4]/80.*T*T*T - a[5]/(T*T) } // dT((H-TS)/RT)
}

#[derive(PartialEq, /*Eq,*/ Debug, Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

pub use model;

impl From<model::RateConstant> for RateConstant {
	fn from(model::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: model::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{log_preexponential_factor: log(preexponential_factor)-temperature_exponent*log(K), temperature_exponent, activation_temperature: activation_energy*J_per_cal/NA}
	}
}

