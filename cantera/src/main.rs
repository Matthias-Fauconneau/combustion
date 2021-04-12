#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents, uncommon_codepoints)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked, destructuring_assignment, test)]
#![feature(unboxed_closures, fn_traits)] // CVODE shim

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
//include!(concat!(env!("OUT_DIR"), "/cantera.rs"));
use std::os::raw::c_char;
#[link(name = "cantera")]
extern "C" {
fn thermo_newFromFile(file_name: *const c_char, phase_name: *const c_char) -> i32;
fn thermo_nSpecies(n: i32) -> usize;
//fn thermo_getMolecularWeights(n: i32, lenm: usize, mw: *mut f64) -> i32;
fn thermo_setTemperature(n: i32, t: f64) -> i32;
fn thermo_setMoleFractions(n: i32, len: usize, x: *const f64, norm: i32) -> i32;
fn thermo_getSpeciesName(n: i32, m: usize, len: usize, buffer: *mut c_char) -> i32;
fn thermo_setPressure(n: i32, p: f64) -> i32;
//fn thermo_getEnthalpies_RT(n: i32, len: usize , H_RT: *mut f64) -> i32;
//fn thermo_getCp_R(n: i32, len: usize, Cp_R: *mut f64) -> i32;

/*fn kin_newFromFile(file_name: *const c_char, phase_name: *const c_char, reactingPhase: i32, neighbor0: i32, neighbor1: i32, neighbor2: i32, neighbor3: i32) -> i32;
fn kin_getEquilibriumConstants(n: i32, len: usize, kc: *mut f64) -> i32;
fn kin_getNetProductionRates(n: i32, len: usize, w_dot: *mut f64) -> i32;
fn kin_getFwdRatesOfProgress(n: i32, len: usize, rate: *mut f64) -> i32;
fn kin_getRevRatesOfProgress(n: i32, len: usize, rate: *mut f64) -> i32;
fn kin_getNetRatesOfProgress(n: i32, len: usize, rate: *mut f64) -> i32;
fn kin_getReactionString(n: i32, i: usize, len: usize, buffer: *mut c_char) -> i32;*/
fn trans_newDefault(th: i32, loglevel: i32) -> i32;
fn trans_viscosity(n: i32) -> f64;
fn trans_thermalConductivity(n: i32) -> f64;
fn trans_getThermalDiffCoeffs(n: i32, ldt: i32, dt: *mut f64) -> i32; // Mass-averaged
//fn trans_getBinDiffCoeffs(n: i32, ld: i32, d: *mut f64) -> i32;
}

use {fehler::throws, error::Error};
pub use combustion::*;
//mod reaction;
mod transport;

#[throws] fn main() {
	let model = std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let ref simulation = Simulation::new(&model)?;
	//reaction::check(model, &simulation)?;
	transport::check(model, &simulation);
}
