/*use std::convert::TryInto;
use fehler::{throws, throw};
use num::log;
use iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoChain, map}, eval, vec::{eval, Dot}};
use model::Troe;
use system::{NASA7, RateConstant, Property};
use super::Species;*/

/*#[derive(PartialEq, /*Eq,*/ Debug, Clone, Copy)] pub enum Model<const S: usize> {
	Elementary,
	ThreeBody { efficiencies: [f64; S] },
	PressureModification { efficiencies: [f64; S], k0: RateConstant },
	Falloff { efficiencies: [f64; S], k0: RateConstant, troe: Troe },
}

#[derive(PartialEq, /*Eq,*/ Clone, Copy, Debug)] pub struct Reaction<const S: usize> where [(); S-1]: {
	pub model: Model<S>,
	pub rate_constant: RateConstant,
	pub reactants: [f64; S],
	pub products: [f64; S],
	pub net: [f64; S-1],
	pub Σreactants: f64,
	pub Σproducts: f64,
	pub Σnet: f64,
}*/

/*#[r#macro::r#macro(MODEL)]
fn dtω/*(/*T: f64, logP0_RT: f64, G_RT: &[f64; S-1], concentrations: &[f64; S-1], log_concentrations: &[f64; S-1]*/)/* -> [f64; S-1] where [(); S-1]:*/ {
	/*use {model::Troe, system::RateConstant};
	let mut dtω = [0.; S-1];
	#fn dot(coefficients)<T, const N: usize>(vector: &[T; N]) -> T {
		(0 as T) #for (index, coefficient) in $coefficients { #match $coefficient {
			0 => ,
			1 => + vector[$index],
			-1 => - vector[$index],
			coefficient => + $coefficient * $vector[$index],
		}}
	}
	#for (_, {reactants, products, net, Σnet, rate_constant, model, ..}) in $reactions {
		let log_k_inf = log_arrhenius($rate_constant, T);
		let c = #match $model { // todo: CSE
			Elementary => 1.,
			ThreeBody{efficiencies} => #dotf($efficiencies)(concentrations),
			PressureModification{efficiencies, k0} => {
					let Pr = #dotf($efficiencies)(concentrations) * f64::exp(log_arrhenius($k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
					Pr / (1.+Pr)
			},
			Falloff{efficiencies, k0, troe} => {
					let Pr = #dotf($efficiencies, concentrations) * f64::exp(log_arrhenius($k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
					let Troe{A, T3, T1, T2} = $troe;
					let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
					let log10Fcent = f64::log10(Fcent);
					let C = -0.4-0.67*log10Fcent;
					let N = 0.75-1.27*log10Fcent;
					let log10PrC = f64::log10(Pr) + C;
					let f1 = log10PrC/(N-0.14*log10PrC);
					let F = num::exp10(log10Fcent/(1.+f1*f1));
					Pr / (1.+Pr) * F
			}
		};
		let Rf = f64::exp(#dot(reactants)(log_concentrations) + log_k_inf);
		let log_equilibrium_constant = -#doti(net)(G_RT) + $Σnet*logP0_RT;
		let Rr = f64::exp(#dot(products)(log_concentrations) + log_k_inf - log_equilibrium_constant);
		let R = Rf - Rr;
		let cR = c * R;
		#for (index, coefficient) in net { #match coefficient: {
			0 => ,
				1 => dtω[$index] += cR,
			-1 => dtω[$index] -= cR,
			coefficient => dtω[$index] += ($coefficient as f64) * cR,
		}; }
	}
	dtω*/
}*/*/
r#macro::r#macro!{rate_and_jacobian}

//impl<const S: usize, const REACTIONS_LEN: usize> System<S, REACTIONS_LEN> where [(); S-1]:, [(); 2+S-1]: {
/*impl<const S: usize, const R: usize> System<S, R> where [(); S-1]:, [(); 2+S-1]: {
	const LEN: usize = S;
}*/
//use CH4::LEN;
//impl System<LEN,325> {
/*#[throws(as Option)]
pub fn rate_and_jacobian<const CONSTANT: Property>(/*const &self,*/ state_constant/*pressure|volume*/: f64, u: &State<CONSTANT, {LEN}>) -> (Derivative<CONSTANT, {Self::LEN}>, /*[[f64; 2+S-1]; 2+S-1]*/) {
		use iter::into::Sum;
		let a = Self::LEN-1;
		let Self{species: super::Species{molar_mass, ref thermodynamics, heat_capacity_ratio,..}, reactions, ..} = self;
		let (T, thermodynamic_state_variable, amounts) = (u[0], u[1]/*volume|pressure*/, u.suffix());
		let ref H_T = eval(thermodynamics.prefix(), |s| s.specific_enthalpy_T(T));
		let (ref E_T/*H|U*/, Cc/*p|v*/, pressure, volume) = {
			let Cp = eval(thermodynamics, |s:&NASA7| s.specific_heat_capacity(T));
			{use Property::*; match CONSTANT {
				Pressure => (*H_T, Cp, state_constant, thermodynamic_state_variable),
				Volume => (eval!(H_T, heat_capacity_ratio.prefix(); |h_T, γ| h_T / γ), eval!(Cp, heat_capacity_ratio; |cp, γ| cp / γ), thermodynamic_state_variable, state_constant),
			}}
		};
		let C = pressure / T; // n/V = P/kT
		let logP0_RT = log(NASA7::reference_pressure) - log(T);
		let ref G_RT = eval!(thermodynamics.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
		let concentrations /*: [_; S-1]*/ = eval(amounts, |&n| n/*.max(0.)*/ / volume); // Skips most abundant specie (last index) (will be deduced from conservation)
		let Ca = C - Sum::<f64>::sum(concentrations);
		//if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
		let ref concentrations = from_iter(concentrations.chain([Ca]));
		let ref log_concentrations = eval(concentrations, |&x| log(x));
		let dtω = dtω(T, logP0_RT, G_RT, concentrations, log_concentrations);
		let rcp_ΣCCc = 1./concentrations.dot(Cc);
		let dtT_T = - rcp_ΣCCc * dtω.dot(E_T);
		let R_S_Tdtn = T / pressure * map(molar_mass.prefix(), |w| 1. - w/molar_mass[a]).dot(dtω); // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
		let dtS_S = R_S_Tdtn + dtT_T;
		let dtn = eval(dtω, |dtω| dtω * volume);
		(State/*Derivative*/(from_iter([dtT_T*T, dtS_S*thermodynamic_state_variable].chain(dtn))), /*J*/)
	}
//}*/
