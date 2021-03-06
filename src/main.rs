#![feature(type_ascription)]#![feature(non_ascii_idents)]#![allow(confusable_idents,non_snake_case,unused_variables,unused_mut)]
use {fehler::throws, error::Error, combustion::{*, Property::*}};

fn exp2(x: f32) -> f32 {
	let x = f32::max(x, -126.99999);
	let ipart = (x-1./2.) as i32;
	let fpart = x - ipart as f32;
	let expipart = f32::from_bits(((ipart+127)<<23) as u32);
	let exp = [0.99999994, 0.69315308, 0.24015361, 0.055826318, 0.0089893397, 0.0018775767];
	let expfpart = (((((exp[5]*fpart + exp[4])*fpart + exp[3])*fpart + exp[2])*fpart + exp[1]) * fpart) + exp[0];
	//let exp: [f32; 3] = [1.0017247, 0.65763628, 0.33718944];
	//let expfpart = ((exp[2]*fpart + exp[1]) * fpart) + exp[0];
	expipart * expfpart
}
//fn exp(x: f32) -> f32 { exp2(x/std::f32::consts::LN_2) }

fn log2(x: f32) -> f32 {
	let i = x.to_bits();
	let exponent = 0x7F800000;
	let e = (((i&exponent)>>23) as i32-127) as f32;
	let mantissa = 0x007FFFFF;
	let m = f32::from_bits( (i&mantissa) | (1f32).to_bits());
	let log = [3.1157899, -3.3241990, 2.5988452, -1.2315303,  0.31821337, -0.034436006];
	let p = (((((log[5]*m + log[4])*m + log[3])*m + log[2])*m + log[1]) * m) + log[0];
	//let log = [2.28330284476918490682, -1.04913055217340124191, 0.204446009836232697516];
	//let p = (log[2] * m + log[1]) * m + log[0];
	let p = p * (m - 1.); //?
	e + p
}
//fn log(x: f32) -> f32 { std::f32::consts::LN_2*log2(x) }

use std::f64::consts::LN_2;
pub fn log_arrhenius_(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: f32) -> f32 {
	log_preexponential_factor as f32 + temperature_exponent as f32*log2(T) - ((activation_temperature/LN_2) as f32)*(1./T)
}

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let ref state = Simulation::new(&model)?.state;
	let model = Model::new(model);
	let (_traps, (_function, _size), rate) = model.rate::<{Volume}>();
	let mut derivative = /*Derivative*/StateVector::<{Volume}>(std::iter::repeat(0.).take(/*model.len()*/2+model.reactions.len()).collect());
	/*{
		let function = unsafe{std::slice::from_raw_parts(function as *const u8, size)}; // refcheck: leaked from dropped JITModule
		let (constant, ref state, derivative) = (state.constant::<{Volume}>(), state.into(): StateVector::<{Volume}>, &mut derivative);
		let constant = constant.0 as f32;
		use x86emu::*;
		let mut guest = State::new();
		allocate_stack(&mut guest);
		load(&mut guest, function);
		let mut heap = Heap::new(&mut guest);
		let state = heap.push_slice(&mut guest, &state.0);
		let derivative = heap.push_slice(&mut guest, &derivative.0);
		call(&mut guest, &[state as i64, derivative as i64], &[constant]);
		pretty_env_logger::init();
		guest.print_instructions = true;
		guest.execute()
	}*/
	//use itertools::Itertools;
	//println!("{:7.0e}", {rate(state.constant(), &state.into(), &mut derivative); derivative.0.iter().format(" ")});
	use std::ops::Deref;
	let vector = {rate(state.constant(), &state.into(), &mut derivative); &derivative.0.deref()[2..]};
	/*for (i, &v) in {rate(state.constant(), &state.into(), &mut derivative); &derivative.0.deref()[2..]}.iter().enumerate() {
		if v.abs() > 1e-29 { print!("{}:{:.0e} ", i, v); }
		//{let v = v-1.; if v.abs() > 1e-29 { print!("{:8}", format!("{}:{:.1}", i, v)); }}
	}
	println!("");*/
	// Test
	let (T, amounts) :(f32,[f32;53])= (1000.,
		[0.,0.,0.,4.874638549881687,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.4373192749408434,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		4.874638549881688]
	);
	//println!("{} {:#?}", T, amounts);
	let logP0_RT = f64::log2(NASA7::reference_pressure) as f32 - log2(T);
	//println!("{}", logP0_RT/LN_2);
	let len = model.len();
	let Model{species: Species{thermodynamics, ..}, reactions} = model;
	let H = thermodynamics.iter().map(|s| s.specific_enthalpy(T as f64) as f32).collect():Box<[f32]>;
	let H_T = H[..len-1].iter().map(|H| H/T).collect():Box<[f32]>;
	use std::f32::consts::LN_2;
	let G = thermodynamics[..len-1].iter().zip(H_T.iter()).map(|(s, h_T)| h_T - s.specific_entropy(T as f64) as f32).map(|v| v/LN_2).collect():Box<[f32]>;
	//use itertools::Itertools;
	//println!("{}", [1./T, log(T)/LN_2].iter().format(" "));
	//println!("{}", [T, T*T, T*T*T, T*T*T*T, 1./T, log(T)/LN_2].iter().format(" "));
	//println!("{}", G.iter().format(" "));
	let concentrations = amounts;
	//let mut sum = 0.; for x in concentrations.iter() { sum += x; }
	//let pressure_R = 101325. / (K*NA);
	//let C = pressure_R / T;
	//let Ca = C - sum;
	//println!("{}", Ca);
	let log_concentrations = concentrations.iter().map(|&x| log2(x)).collect():Box<[f32]>;
	//println!("{}", log_concentrations.iter().map(|&x| x/LN_2).format(" "));
	let mut dtω = vec![0.; len-1];
	//for Reaction{reactants, products, net, Σnet, rate_constant, model, ..} in reactions.iter() {
	for (i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..}) in reactions.iter().enumerate() {
		let log_k_inf = log_arrhenius_(*rate_constant, T);
		//print!("{} ", log_k_inf/LN_2);
		use model::Troe;
		let c = match model {
			ReactionModel::Elementary => 1.,
			ReactionModel::ThreeBody{efficiencies} => {
				let mut sum = 0.; for (&a, b) in efficiencies.iter().zip(concentrations.iter()) { sum += a as f32 * b; }
				sum
			}
			ReactionModel::PressureModification{efficiencies, k0} => {
				let mut sum = 0.; for (&a, b) in efficiencies.iter().zip(concentrations.iter()) { sum += a as f32 * b; }
				let Pr = sum * exp2(log_arrhenius_(*k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
				Pr / (1.+Pr)
			}
			ReactionModel::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
				let mut sum = 0.; for (&a, b) in efficiencies.iter().zip(concentrations.iter()) { sum += a as f32* b; }
				let Pr = sum * exp2(log_arrhenius_(*k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
				fn rcp(x: f64) -> f64 { 1./x }
				use std::f64::consts::LN_2;
				let T = T as f32;
				let Pr = Pr as f32;

				let Fcent = (1.-*A as f32)*f32::exp2((-T)*(rcp(LN_2*T3) as f32)) +
													(*A as f32)*f32::exp2((-T)*(rcp(LN_2*T1) as f32)) +
													exp2((-1./T)*((T2/LN_2) as f32));
				let logFcent = log2(Fcent);
				let c = -0.67 * logFcent + (-0.4*f64::log2(10.) as f32);
				let N = -1.27 * logFcent + (0.75*f64::log2(10.) as f32);
				let logPrc = log2(Pr) + c;
				let f1 = logPrc / (-0.14 * logPrc + N);
				let F = exp2(logFcent / (f1*f1 + 1.));
				/*{
					let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
					let log10Fcent = f64::log10(Fcent);
					let C = -0.4-0.67*log10Fcent;
					let N = 0.75-1.27*log10Fcent;
					let log10PrC = f64::log10(Pr) + C;
					let f1 = log10PrC/(N-0.14*log10PrC);
					let F = num::exp10(log10Fcent/(1.+f1*f1));
				}*/
				(Pr / (1. + Pr) * F) as f32
				//Fcent as f64
			}
		};
		//{let v=c-1.; if v != 0. { print!("{}:{:.0e} ", i, v); }}
		//{let v=c-1.; if v != 0. { print!("{:8}", format!("{}:{:.1}", i, v)); }}
		let mut sum = 0.; for (&a, b) in reactants.iter().zip(log_concentrations.iter()) { if a != 0 { sum += a as f32 * b; } }
		//{let v=sum + log_k_inf; if v != 0. { print!("{}:{:.0e} ", i, v); }}
		let Rf = exp2(sum + log_k_inf);
		//{let v=Rf; if v != 0. { print!("{}:{:.0e} ", i, v); }}
		let mut sum = 0.; for (&a, b) in net.iter().zip(G.iter()) { if a != 0 { sum += a as f32 * b; } }
		let log_equilibrium_constant = -sum + (*Σnet as f32)*logP0_RT;
		//{let v=-log_equilibrium_constant; if v != 0. { print!("{}:{:.0e} ", i, v); }}
		// exp(x) = exp2(x/LN_2)
		// exp2(x) = exp(x*LN_2)
		//use std::f32::consts::LN_2;
		//{let v=log_equilibrium_constant; print!("{:8}",format!("{:.1}|{:.1} ", v, -vector[i]*LN_2)); }
		let mut sum = 0.; for (&a, b) in products.iter().zip(log_concentrations.iter()) { if a != 0 { sum += a as f32 * b; } }
		let Rr = exp2(sum + log_k_inf - log_equilibrium_constant);
		let R = Rf - Rr;
		//{let v=c-1.; if v != 0. { print!("{}:{}|{} ", i, v, vector[i]-1.); }}
		//{let v=R; if v != 0. { print!("{}:{}|{} ", i, v, vector[i]); }}
		//{let v=R; if v != 0. { print!("{}:{} ", i, num::relative_error(v, vector[i] as f64)); }}
		let cR = c * R;
		//{let v=cR; if v != 0. { print!("{}:{:.0e}|{:.0e} ", i, v, vector[i]); }}
		//{let v=cR; if v != 0. { print!("{}:{} ", i, num::relative_error(v, vector[i] as f64)); }}
		for (specie, &ν) in net.iter().enumerate() { dtω[specie] += (ν as f32) * cR; }
	}
	for (i,(&a, &b)) in dtω.iter().zip(vector).enumerate() { if a.abs() > 1e-29 || b.abs() > 1e-29 { print!("{}:{:.3e}|{:.3e} ", i, a, b); } }
	//for (i,(&a, &b)) in dtω.iter().zip(vector).enumerate() { if a != 0. || b.abs() > 1e-29 { print!("{}:{} ", i, num::relative_error(a as f64, b as f64)); } }
	println!("");
	//use itertools::Itertools;
	//println!("{}", dtω.iter().zip(vector).format(" "));
	/*let Cp = thermodynamics.iter().map(|s| s.specific_heat_capacity(T)).collect():Box<[f64]>;
	let mut sum = 0.; for (a, b) in concentrations.iter().zip(Cp.iter()) { sum += a * b; }
	let rcp_ΣCCp = 1./sum;
	let mut sum = 0.; for (a, b) in dtω.iter().zip(H_T.iter()) { sum += a * b; }
	let dtT_T = - rcp_ΣCCp * sum; // R/RT
	let dtn = dtω;*/
	//dbg!(dtT_T*T, dtn);
}
