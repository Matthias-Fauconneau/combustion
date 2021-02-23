use std::convert::TryInto;
use fehler::{throws, throw};
use num::log;
use iter::{Prefix, Suffix, array_from_iter as from_iter, into::{Enumerate, IntoChain, map}, eval, vec::{eval, Dot}};
use model::Troe;
use system::{NASA7, RateConstant};
use super::System;

pub fn log_arrhenius(&RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: &RateConstant, T: f64) -> f64 {
	log_preexponential_factor + temperature_exponent*log(T) - activation_temperature*(1./T)
}

#[derive(Debug, Clone, Copy)] pub enum Model<const S: usize> {
	Elementary,
	ThreeBody { efficiencies: [f64; S] },
	PressureModification { efficiencies: [f64; S], k0: RateConstant },
	Falloff { efficiencies: [f64; S], k0: RateConstant, troe: Troe },
}

impl<const S: usize> Model<S> {
pub fn efficiency(&self, T: f64, concentrations: &[f64; S], log_k_inf: f64) -> f64 {
	match self {
		Self::Elementary => 1.,
		Self::ThreeBody{efficiencies} => efficiencies.dot(concentrations),
		Self::PressureModification{efficiencies, k0} => {
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
			Pr / (1.+Pr)
		}
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
			let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
			let log10Fcent = f64::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			let N = 0.75-1.27*log10Fcent;
			let log10PrC = f64::log10(Pr) + C;
			let f1 = log10PrC/(N-0.14*log10PrC);
			let F = num::exp10(log10Fcent/(1.+f1*f1));
			Pr / (1.+Pr) * F
		}
	}
}
}

#[derive(Clone, Copy, Debug)] pub struct Reaction<const S: usize> where [(); S-1]: {
	pub reactants: [f64; S],
	pub products: [f64; S],
	pub net: [f64; S-1],
	pub Σreactants: f64,
	pub Σproducts: f64,
	pub Σnet: f64,
	pub rate_constant: RateConstant,
	pub model: Model<S>,
}

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }

#[derive(Debug)] pub struct State<const CONSTANT: Property, const S: usize>(pub [f64; 2+S-1]) where [(); 2+S-1]:;
pub type Derivative<const CONSTANT: Property, const S: usize> = State<CONSTANT, S>;

impl<const CONSTANT: Property, const S: usize> std::ops::Deref for State<CONSTANT, S> where [(); 2+S-1]: {
	type Target = [f64; 2+S-1]; fn deref(&self) -> &Self::Target { &self.0 }
}

impl<const CONSTANT: Property, const S: usize> From<&super::State<S>> for State<CONSTANT, S> where [(); S-1]:, [(); 2+S-1]: {
	fn from(super::State{temperature, pressure, volume, amounts}: &super::State<S>) -> Self {
		Self(from_iter([*temperature, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure}}].chain(amounts[..S-1].try_into().unwrap():[_;S-1])))
	}
}

impl<const S: usize, const REACTIONS_LEN: usize> super::System<S, REACTIONS_LEN> where [(); S-1]:, [(); 2+S-1]: {
#[throws(as Option)]
pub fn rate_and_jacobian<const CONSTANT: Property>(/*const*/ &self, state_constant/*pressure|volume*/: f64, u: &State<CONSTANT, S>) -> (Derivative<CONSTANT, S>, /*[[f64; 2+S-1]; 2+S-1]*/) {
		use iter::into::{IntoMap, Sum};
		let a = S-1;
		let System{species: super::Species{molar_mass: W, thermodynamics, heat_capacity_ratio: γ,..}, reactions, ..} = self;
		//let rcpV = 1. / V;
		let (T, thermodynamic_state_variable, amounts) = (u[0], u[1]/*volume|pressure*/, u.suffix());
		//let ref H = eval(thermodynamics, |s| s.specific_enthalpy(T));
		//let ref H_T = eval(H.prefix(), |H| H/T);
		let ref H_T = eval(thermodynamics.prefix(), |s| s.specific_enthalpy_T(T));
		let (ref E_T/*H|U*/, Cc/*p|v*/, pressure, volume) = {
			let Cp = eval(thermodynamics, |s:&NASA7| s.specific_heat_capacity(T));
			{use Property::*; match CONSTANT {
				Pressure => (*H_T, Cp, state_constant, thermodynamic_state_variable),
				Volume => (eval!(H_T, γ.prefix(); |h_T, γ| h_T / γ), eval!(Cp, γ; |cp, γ| cp / γ), thermodynamic_state_variable, state_constant),
			}}
		};
		assert!(T>0.);
		let C = pressure / T; // n/V = P/kT
		//let rcp_C = 1. / C;
		//let rcp_amount = rcpV * rcp_C;
		let logP0_RT = log(NASA7::reference_pressure) - log(T);
		let ref G_RT = eval!(thermodynamics.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
		//let ref dT_G = eval!(thermodynamics.prefix(); |s| s.dT_Gibbs_free_energy(T));
		//let ref amounts = eval(amounts, |n| n.max(0.));
		for &n in amounts { assert!(n>=0.); }
		let concentrations : [_; S-1] = eval(amounts, |&n| n/*.max(0.)*/ / volume); // Skips most abundant specie (last index) (will be deduced from conservation)
		let Ca = C - Sum::<f64>::sum(concentrations);
		if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
		assert!(Ca>0.,"{:?}", (C, Sum::<f64>::sum(concentrations), concentrations));
		let ref concentrations = from_iter(concentrations.chain([Ca]));
		let ref log_concentrations = eval(concentrations, |&x| if x==0. { -f64::INFINITY } else { assert!(x>0.); log(x) }); // Explicit to avoid signal
		let ref mut dtω = [0.; S-1];
		/*let mut dTω = [0.; S-1];
		let mut dVω = [0.; S-1];
		let mut dnω = [[0.; S-1]; S-1];*/
		for Reaction{reactants, products, net, /*Σreactants, Σproducts,*/ Σnet, rate_constant/*: rate_constant@RateConstant{temperature_exponent, activation_temperature, ..}*/, model, ..} in reactions.iter() {
			let log_kf = log_arrhenius(rate_constant, T);
			assert!(log_kf.is_finite(), "{:?}", (rate_constant, T));
			let c = model.efficiency(T, concentrations, log_kf);
			let mask = |mask, v| iter::zip!(mask, v).map(|(&mask, v):(_,&_)| if mask != 0. { *v } else { 0. });
			let Rf = f64::exp(reactants.dot(mask(reactants, log_concentrations)) + log_kf);
			let log_equilibrium_constant = -net.dot(G_RT) + Σnet*logP0_RT;
			let Rr = f64::exp(products.dot(mask(products, log_concentrations)) + log_kf - log_equilibrium_constant);
			assert!(Rf.is_finite() && Rr.is_finite(), "{:?}", (Rf, Rr, reactants, products, log_concentrations));
			let R = Rf - Rr;
			let cR = c * R;

			/*let efficiencies = match model {Model::Elementary => &[0.; S], Model::ThreeBody{efficiencies}|Model::Falloff{efficiencies,..} => efficiencies};
			let has = eval(efficiencies, |&e| if e != 0. { 1. } else { 0. });
			let νfRfνrRr = eval!(reactants, products; |νf,νr| νf*Rf - νr*Rr);

			let dTc = has[a] * -c/T;
			// dT(R) = (β+Ea/(T.Rideal))/T.R + Rr. Σ ν.dT(G) - νfRfνrRr[a]/T
			let dTR = (temperature_exponent+activation_temperature/T)/T*R + Rr*net.dot(dT_G) - νfRfνrRr[a] / T;
			let RdTccdTR = R * dTc + c * dTR;

			let dVc = rcpV * (has[a]*c  - efficiencies.sum():f64);
			// dV(R) = 1/V . ( (kf.Sf-kr.Sr) - (Σνf.Rf - Σνr.Rr) )
			let dVR = rcpV * ( νfRfνrRr[a] - (Σreactants*Rf - Σproducts*Rr));
			let RdVccdVR = R * dVc + c * dVR;

			let dnc = map(has.prefix(), |has_k| rcpV * (has_k - has[a]));
			// dn(R) = 1/n . ( kf.(Sf-Sfa) - kr.(Sr-Sra) )
			let dnR = map(νfRfνrRr.prefix(), |νfRfνrRrj| rcp_amount * (νfRfνrRrj - νfRfνrRr[a]));
			let RdnccdnR : [_; S-1] = eval!(dnc, dnR; |dnc,dnR| R*dnc + c*dnR);*/

			for (specie, ν) in net.enumerate() {
				// dtω = Σ ν c R
				dtω[specie] += ν * cR;
				/*// dT(ω̇̇̇̇̇̇̇̇̇̇) = Σ ν.(R.dT(c)+c.dT(R))
				dTω[specie] += ν * RdTccdTR;
				// dV(ω) = Σ ν.(R.dV(c)+c.dV(R))
				dVω[specie] += ν * RdVccdVR;
				// dn(ω) = Σ ν.(R.dn(c)+c.dn(R))
				for (dnω, RdnccdnR) in zip!(&mut dnω[specie], RdnccdnR) { *dnω += ν * RdnccdnR; }*/
			}
		}
		let dtω = *dtω;

		let rcp_ΣCCc = 1./concentrations.dot(Cc);
		let dtT_T = - rcp_ΣCCc * dtω.dot(E_T);
		assert!(dtT_T.is_finite(), "{:?}", (dtT_T, dtω, E_T));
		let R_S_Tdtn = T / pressure * map(W.prefix(), |w| 1. - w/W[a]).dot(dtω); // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
		assert!(R_S_Tdtn.is_finite(), "{}", R_S_Tdtn);
		let dtS_S = R_S_Tdtn + dtT_T;
		let dtn = eval(dtω, |dtω| dtω * volume);

		/*let mut J = [[f64::NAN; 2+S-1]; 2+S-1];
		let dtT = - rcp_ΣCCp * H.prefix().dot(dtω);
		let Cpa = Cp[a];
		let HaWa = H[a]/W[a];
		let HaWaWH = eval!(W.prefix(), H.prefix(); |W, H| HaWa*W - H);
		let Cpa_Wa = Cpa/W[a];
		let concentrations = concentrations.prefix();
		let Cp = Cp.prefix();
		// dT(dtT) = 1 / (Σ C.Cp) . [ dtT . Σ C.(Cpa/T - dT(Cp)) + Σ_ ( (Ha/Wa*W - H).dT(ω) + (Cpa/Wa.W - Cp).ω ) ]
		let dTdtT = rcp_ΣCCp * (dtT * concentrations.dot(thermodynamics.map(|s:&NASA7| Cpa/T - s.dT_specific_heat_capacity(T))) + HaWaWH.dot(dTω) + map!(W.prefix(), Cp; |W, Cp| Cpa_Wa*W - Cp).dot(dtω));
		J[0][0] = dTdtT;
		// dn(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dn(ω) + dtT/V . (Cpa-Cp) ]
		let ref dndtT = eval!(dnω, Cp; |dnω, Cp| rcp_ΣCCp * (HaWaWH.dot(dnω) + rcpV * dtT * (Cpa-Cp) ));
		J[0][2..2+S-1].copy_from_slice(dndtT);
		let dTdtn = dTω.scale(V);
		for (k, dTdtn) in dTdtn.enumerate() { J[2+k][0] = dTdtn; }
		let dndtn = generate(|k| generate(|l| V*dnω[l][k])):[[_;S-1];S-1]; // Transpose [l][k] -> [k][l]
		for l in 0..S-1 { for (k, dndtn) in dndtn[l].enumerate() { J[2+k][l] = dndtn; } } // Transpose back*/
		assert!(dtS_S.is_finite(), "{:?}",(dtS_S, rcp_ΣCCc, dtω, E_T));
		for dtn in &dtn { assert!(dtn.is_finite()); }
		(State/*Derivative*/(from_iter([dtT_T*T, dtS_S*thermodynamic_state_variable].chain(dtn))), /*J*/)
	}
}
