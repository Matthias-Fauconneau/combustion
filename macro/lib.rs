#![feature(non_ascii_idents, type_ascription, once_cell, array_map, proc_macro_quote, in_band_lifetimes)]
#![allow(confusable_idents, non_upper_case_globals, non_snake_case, unused_variables)]
use std::ops::Deref;
use system::{Model, Reaction, System};

use proc_macro::TokenStream;
fn literal(value: impl std::fmt::Debug) -> TokenStream { format!("{:?}", value).parse().unwrap() }
//use syn::{ItemFn, Expr, parse_quote as quote};
use proc_macro::quote; type Expr = TokenStream; type ItemFn = TokenStream;

fn dot<T: num::IsZero + num::IsOne + num::IsMinusOne + std::fmt::Debug/*, const N: usize*/>(coefficients: &[T/*; N*/]) -> ItemFn {
	let adds = coefficients.into_iter().enumerate().filter_map(|(index, coefficient)| { let index = literal(index);
		//use proc_macro::quote;
		if coefficient.is_zero() { None }
    else if coefficient.is_one() { Some(quote!(+ vector[$index])) }
    else if coefficient.is_minus_one() { Some(quote!(+ vector[$index])) }
    else { let coefficient = literal(coefficient); Some(quote!(+ ($coefficient as f64) * vector[$index])) }
	}).collect(): TokenStream;
	let literal_N = literal(/*N*/coefficients.len());
	quote!({fn dot<const N: usize>(vector: &[f64; N]) -> f64 { 0. $adds } dot::<{$literal_N}>})
}

fn efficiency(model: &Model) -> Expr/*fn<const S: usize>(T: f64, concentrations: &[f64; S], log_k_inf: f64) -> f64*/ {
	use Model::*;
	match model {
		Elementary => quote!(1.),
		ThreeBody{efficiencies} => { let dot_efficiencies = dot(efficiencies.deref()); quote!($dot_efficiencies(concentrations)) },
		PressureModification{efficiencies, k0} => {
			let dot_efficiencies = dot(efficiencies.deref());
			let k0 = literal(k0);
			quote!({
				let Pr = $dot_efficiencies(concentrations) * f64::exp(log_arrhenius($k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
				Pr / (1.+Pr)
			})
		}
		Falloff{efficiencies, k0, troe} => {
			let dot_efficiencies = dot(efficiencies.deref());
			let k0 = literal(k0);
			let troe = literal(troe);
			quote!({
				let Pr = $dot_efficiencies(concentrations) * f64::exp(log_arrhenius($k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
				let model::Troe{A, T3, T1, T2} = $troe;
				let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
				let log10Fcent = f64::log10(Fcent);
				let C = -0.4-0.67*log10Fcent;
				let N = 0.75-1.27*log10Fcent;
				let log10PrC = f64::log10(Pr) + C;
				let f1 = log10PrC/(N-0.14*log10PrC);
				let F = num::exp10(log10Fcent/(1.+f1*f1));
				Pr / (1.+Pr) * F
			})
		}
	}
}

//#[proc_macro_attribute]
//pub fn r#macro(_attribute: TokenStream, item: /*ItemFn*/TokenStream) -> /*ItemFn*/TokenStream {
#[proc_macro]
pub fn r#macro(ident: /*Ident*/TokenStream) -> /*ItemFn*/TokenStream {
	let model = std::fs::read([&std::env::var("MODEL").unwrap(),".ron"].concat()).unwrap();
	let system = System::new(::ron::de::from_bytes(&model).unwrap());
	let len = literal(system.len());
	let molar_mass = literal(system.species.molar_mass);
	let thermodynamics = literal(system.species.thermodynamics);
	let heat_capacity_ratio = literal(system.species.heat_capacity_ratio);
	let reactions = system.reactions.iter().map(|Reaction{reactants, products, net, Σnet, rate_constant, model, ..}| {
		let rate_constant = literal(rate_constant);
		let Σnet = literal(Σnet);
		let efficiency_model = efficiency(&model); // todo: CSE
		let dot_reactants = dot(reactants.deref());
		let dot_net = dot(net.deref());
		let dot_products = dot(products.deref());
		let accumulate = net.iter().enumerate().map(|(index, coefficient)| { let index = literal(index); let compound_assignment_expression = match coefficient {
				0 => {quote!()},
				1 => quote!(dtω[$index] += cR),
				-1 => quote!(dtω[$index] -= cR),
				coefficient => { let coefficient = literal(coefficient); quote!(dtω[$index] += ($coefficient as f64) * cR) },
		};
		quote!( $compound_assignment_expression; ) } ).collect(): TokenStream;
		quote!{
			let log_k_inf = log_arrhenius($rate_constant, T);
			let c = $efficiency_model/*(T, concentrations, log_k_inf)*/;
			let Rf = f64::exp($dot_reactants(log_concentrations) + log_k_inf);
			let log_equilibrium_constant = -$dot_net(G_RT) + ($Σnet as f64)*logP0_RT;
			let Rr = f64::exp($dot_products(log_concentrations) + log_k_inf - log_equilibrium_constant);
			let R = Rf - Rr;
			let cR = c * R;
			$accumulate/*(&mut dtω, cR)*/
		}
	}).collect(): TokenStream;
	//use ::quote::ToTokens;
	/*let mut iter = ident.into_iter();
	let ident_len = iter.next().unwrap();
	let ident_dtω = iter.next().unwrap();
	quote!{
		pub const $ident_len: usize = $len;
		pub fn $ident_dtω(T: f64, logP0_RT: f64, G_RT: &[f64; $len-1], concentrations: &[f64; $len], log_concentrations: &[f64; $len]) -> [f64; $len-1] {
			use {model::Troe, system::RateConstant};
			fn log_arrhenius(::system::RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: ::system::RateConstant, T: f64) -> f64 {
				log_preexponential_factor + temperature_exponent*num::log(T) - activation_temperature*(1./T)
			}
			let mut dtω = [0.; $len-1];
			$reactions
			dtω
		}
	}
	//syn::parse_macro_input!(q as syn::ItemFn).into_token_stream().into()*/
quote!{
use {num::log, iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoChain, map, Sum}, vec::{eval, Dot}, eval}};
use system::{NASA7, RateConstant, Property, State, Derivative};
#[fehler::throws(as Option)]
pub fn $ident<const CONSTANT: Property>(state_constant: f64, u: &State<CONSTANT, $len>) -> (Derivative<CONSTANT, $len>, /*[[f64; 2+$len-1]; 2+$len-1]*/) {
	const molar_mass: [f64; $len] = $molar_mass;
	const thermodynamics: [NASA7; $len] = $thermodynamics;
	const heat_capacity_ratio: [f64; $len] = $heat_capacity_ratio;
	let (T, thermodynamic_state_variable, amounts) = (u[0], u[1]/*volume|pressure*/, u.suffix());
	let ref H_T = eval(thermodynamics.prefix(), |s| s.specific_enthalpy_T(T));
	let (ref E_T/*H|U*/, Cc/*p|v*/, pressure, volume) = {
		let Cp = eval(&thermodynamics, |s:&NASA7| s.specific_heat_capacity(T));
		{use Property::*; match CONSTANT {
			Pressure => (*H_T, Cp, state_constant, thermodynamic_state_variable),
			Volume => (eval!(H_T, heat_capacity_ratio.prefix(); |h_T, γ| h_T / γ), eval!(Cp, heat_capacity_ratio; |cp, γ| cp / γ), thermodynamic_state_variable, state_constant),
		}}
	};
	let C = pressure / T; // n/V = P/kT
	let logP0_RT = log(NASA7::reference_pressure) - log(T);
	let ref G_RT = eval!(thermodynamics.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
	let concentrations: [_; $len-1] = eval(amounts, |&n| n/*.max(0.)*/ / volume); // Skips most abundant specie (last index) (will be deduced from conservation)
	let Ca = C - Sum::<f64>::sum(concentrations);
	//if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
	let ref concentrations = from_iter(concentrations.chain([Ca]));
	let ref log_concentrations = eval(concentrations, |&x| log(x));
	//let dtω = dtω(T, logP0_RT, G_RT, concentrations, log_concentrations);
	fn log_arrhenius(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: f64) -> f64 {
		log_preexponential_factor + temperature_exponent*num::log(T) - activation_temperature*(1./T)
	}
	let mut dtω = [0.; $len-1];
	$reactions
	let rcp_ΣCCc = 1./concentrations.dot(Cc);
	let dtT_T = - rcp_ΣCCc * dtω.dot(E_T);
	let R_S_Tdtn = T / pressure * map(molar_mass.prefix(), |w| 1. - w/molar_mass[$len-1]).dot(dtω); // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
	let dtS_S = R_S_Tdtn + dtT_T;
	let dtn = eval(dtω, |dtω| dtω * volume);
	(State/*Derivative*/(from_iter([dtT_T*T, dtS_S*thermodynamic_state_variable].chain(dtn))), /*J*/)
}
}
}
