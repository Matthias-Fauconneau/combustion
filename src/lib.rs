#![allow(incomplete_features)]#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals,unused_imports)]
pub mod ron;
use {iter::{Prefix, Suffix, array_from_iter as from_iter, into::{Enumerate, IntoChain, IntoMap, map}, zip, map, eval, vec::{eval, Dot, generate, Scale, Sub}}, self::ron::{Map, Element, Troe}};

pub const ideal_gas_constant : f64 = 8.31446261815324; // J⋅K−1⋅mol−1

#[derive(Debug)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure_R : f64 = 101325. / ideal_gas_constant; // 1 atm
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R
	//fn dT_specific_heat_capacity(&self, T: f64) -> f64 { 	let a = self.a(T); ideal_gas_constant * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T) } // /R
	//fn dT_Gibbs_free_energy(&self, T: f64) -> f64 { let a = self.a(T); (1.-a[0])/T - a[1]/2. - a[2]/12.*T - a[3]/36.*T*T - a[4]/80.*T*T*T - a[5]/(T*T) } // dT((H-TS)/RT)
}

#[derive(Debug, Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

impl From<ron::RateConstant> for RateConstant { fn from(ron::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: ron::RateConstant) -> Self {
	Self{log_preexponential_factor: f64::ln(preexponential_factor), temperature_exponent, activation_temperature: activation_energy*4.184/ideal_gas_constant}
}}

pub fn log_arrhenius(&RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: &RateConstant, T: f64) -> f64 {
	log_preexponential_factor + temperature_exponent*f64::ln(T) - activation_temperature*(1./T)
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
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(k0, T) - log_k_inf); // [k0/kinf] = [C] (m3/mol)
			Pr / (1.+Pr)
		}
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(k0, T) - log_k_inf); // [k0/kinf] = [C] (m3/mol)
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

#[derive(Clone, Copy)] pub struct Reaction<const S: usize> where [(); S-1]: {
	pub reactants: [f64; S],
	pub products: [f64; S],
	pub net: [f64; S-1],
	pub Σreactants: f64,
	pub Σproducts: f64,
	pub Σnet: f64,
	pub rate_constant: RateConstant,
	pub model: Model<S>,
}

pub struct System<const S: usize> where [(); S-1]: {
	pub molar_masses: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S>]>,
}

impl<const S: usize> System<S> where [(); S-1]:, [(); 1+S-1]: {
	pub fn dt_J(&self, pressure_R: f64, y: &[f64; 1+S-1]) -> ([f64; 1+S-1], /*[[f64; 1+S-1]; 1+S-1]*/) {
		use iter::into::Sum;
		//let a = S-1;
		let Self{thermodynamics: species, reactions/*, molar_masses: W*/, ..} = self;
		//let rcpV = 1. / V;
		let (T, amounts) = (y[0], y.suffix());
		let C = pressure_R / T;
		//let rcp_C = 1. / C;
		//let rcp_amount = rcpV * rcp_C;
		let logP0_RT = f64::ln(NASA7::reference_pressure_R) - f64::ln(T);
		let ref H = eval(species, |s| s.specific_enthalpy(T));
		let ref H_T = eval(H.prefix(), |H| H/T);
		let ref G = eval!(species.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
		//let ref dT_G = eval!(species.prefix(); |s| s.dT_Gibbs_free_energy(T));
		let concentrations : [_; S-1] = eval(amounts, |&n| n/*.max(0.)*/); // Skips most abundant specie (last index) (will be deduced from conservation)
		let Ca = C - concentrations.sum():f64;
		let ref concentrations = from_iter(concentrations.chain([Ca]));
		let ref log_concentrations = eval(concentrations, |&x| f64::ln(x));
		let ref mut dtω = [0.; S-1];
		/*let mut dTω = [0.; S-1];
		let mut dVω = [0.; S-1];
		let mut dnω = [[0.; S-1]; S-1];*/
		for Reaction{reactants, products, net, /*Σreactants, Σproducts,*/ Σnet, rate_constant/*: rate_constant@RateConstant{temperature_exponent, activation_temperature, ..}*/, model, ..} in reactions.iter() {
			let log_kf = log_arrhenius(rate_constant, T);
			let c = model.efficiency(T, concentrations, log_kf);
			let mask = |mask, v| iter::zip!(mask, v).map(|(&mask, v):(_,&_)| if mask != 0. { *v } else { 0. });
			let Rf = f64::exp(reactants.dot(mask(reactants, log_concentrations)) + log_kf);
			let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
			let Rr = f64::exp(products.dot(mask(products, log_concentrations)) + log_kf - log_equilibrium_constant);
			//assert!(Rf.is_finite() && Rr.is_finite(), Rf, Rr, reactants, products, log_concentrations);
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

		let Cp = eval(species, |s:&NASA7| s.specific_heat_capacity(T));
		let rcp_ΣCCp = 1./concentrations.dot(Cp); // All?
		let dtT_T = - rcp_ΣCCp * dtω.dot(H_T); // R/RT
		let dtn = dtω;

		/*let mut J = [[f64::NAN; 1+S-1]; 1+S-1];
		let dtT = - rcp_ΣCCp * H.prefix().dot(dtω);
		let Cpa = Cp[a];
		let HaWa = H[a]/W[a];
		let HaWaWH = eval!(W.prefix(), H.prefix(); |W, H| HaWa*W - H);
		let Cpa_Wa = Cpa/W[a];
		let concentrations = concentrations.prefix();
		let Cp = Cp.prefix();
		// dT(dtT) = 1 / (Σ C.Cp) . [ dtT . Σ C.(Cpa/T - dT(Cp)) + Σ_ ( (Ha/Wa*W - H).dT(ω) + (Cpa/Wa.W - Cp).ω ) ]
		let dTdtT = rcp_ΣCCp * (dtT * concentrations.dot(species.map(|s:&NASA7| Cpa/T - s.dT_specific_heat_capacity(T))) + HaWaWH.dot(dTω) + map!(W.prefix(), Cp; |W, Cp| Cpa_Wa*W - Cp).dot(dtω));
		J[0][0] = dTdtT;
		// dn(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dn(ω) + dtT/V . (Cpa-Cp) ]
		let ref dndtT = eval!(dnω, Cp; |dnω, Cp| rcp_ΣCCp * (HaWaWH.dot(dnω) + rcpV * dtT * (Cpa-Cp) ));
		J[0][1..1+S-1].copy_from_slice(dndtT);
		let dTdtn = dTω.scale(V);
		for (k, dTdtn) in dTdtn.enumerate() { J[1+k][0] = dTdtn; }
		let dndtn = generate(|k| generate(|l| V*dnω[l][k])):[[_;S-1];S-1]; // Transpose [l][k] -> [k][l]
		for l in 0..S-1 { for (k, dndtn) in dndtn[l].enumerate() { J[1+k][l] = dndtn; } } // Transpose back*/
		(from_iter([dtT_T].chain(dtn)), /*J*/)
	}
}

#[derive(Clone)] pub struct State<const S: usize> where [(); S-1]: {
	pub temperature: f64,
	pub amounts: [f64; S-1]
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species: [&'t str; S],
	pub system: System<S>,
	pub time_step: f64,
	pub pressure_r: f64,
	pub state: State<S>
}

use std::lazy::SyncLazy;
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl<const S: usize> Simulation<'t, S> where [(); S-1]: {
	pub fn new(system: &'b [u8]) -> ::ron::Result<Self> where 'b: 't, [(); S]: {
		let ron::System{species: species_data, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
		let ron::Phase::IdealGas{species, state, ..} = phases.into_vec().into_iter().next().unwrap();
		use std::convert::TryInto;
		let species : [_; S] = species.as_ref().try_into().unwrap_or_else(|_| panic!("Compiled for {} species, got {}", S, species.len()));
		let molar_masses = eval(species, |s| species_data[s].composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = eval(species, |s| { let ron::Specie{thermodynamic: ron::NASA7{temperature_ranges, pieces},..} = &species_data[s]; match temperature_ranges[..] {
			[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
			[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
			ref ranges => panic!("{:?}", ranges),
		}});
		#[allow(unused_variables)]
		let reactions = iter::into::Collect::collect(reactions.map(|self::ron::Reaction{ref equation, rate_constant, model}| {
			let [reactants, products] = eval(equation, |e| eval(species, |s| *e.get(s).unwrap_or(&0) as f64));
			let net = iter::vec::Sub::sub(products.prefix(), reactants.prefix());
			let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
			let Σnet = Σproducts-Σreactants;
			let from = |efficiencies:Map<_,_>| eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.));
			Reaction{
				reactants, products, net, Σreactants, Σproducts, Σnet,
				rate_constant: rate_constant.into(),
				model: {use self::ron::Model::*; match model {
					Elementary => Model::Elementary,
					ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: from(efficiencies)},
					PressureModification{efficiencies, k0} => Model::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
					Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe},
				}},
			}
		}));

		let ron::InitialState{temperature, pressure, mole_proportions} = state;
		let pressure_r = pressure / ideal_gas_constant;
		let amount = pressure_r / temperature;
		let mole_proportions = eval(&species, |specie| *mole_proportions.get(specie).unwrap_or(&0.));
		let amounts = eval(mole_proportions.prefix(), |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);

		Ok(Self{
			species,
			system: System{molar_masses, thermodynamics, reactions},
			time_step, pressure_r,
			state: State{temperature, amounts}
		})
	}
}
