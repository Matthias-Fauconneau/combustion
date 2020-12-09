#![allow(incomplete_features)]#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals,unused_imports)]
pub mod ron;
use {iter::{Prefix, Suffix, array_from_iter as from_iter, into::{Enumerate, IntoChain, IntoMap, map}, zip, map, eval, vec::{eval, Dot, generate, Scale, Sub}}, self::ron::{Map, Element, Troe}};

#[derive(Debug)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325.; // 1 atm
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn dimensionless_specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T }
	pub fn dimensionless_specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T }
	//pub fn dimensionless_specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T }
	pub fn dimensionless_specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T }
	fn dT_Cp(&self, T: f64) -> f64 { 	let a = self.a(T); ideal_gas_constant * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T) }
	fn dT_G(&self, T: f64) -> f64 { let a = self.a(T); (1.-a[0])/T - a[1]/2. - a[2]/12.*T - a[3]/36.*T*T - a[4]/80.*T*T*T - a[5]/(T*T) } // dT((H-TS)/RT)
}

pub const ideal_gas_constant : f64 = 8.31446261815324; // J⋅K−1⋅mol−1

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
	Falloff { efficiencies: [f64; S], k0: RateConstant, troe: Troe },
}

impl<const S: usize> Model<S> {
pub fn efficiency(&self, T: f64, concentrations: &[f64; S], log_k_inf: f64) -> f64 {
	match self {
		Self::Elementary => 1.,
		Self::ThreeBody{efficiencies} => efficiencies.dot(concentrations),
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
	Σreactants: f64,
	Σproducts: f64,
	pub Σnet: f64,
	pub rate_constant: RateConstant,
	pub model: Model<S>,
}

pub struct System<const S: usize> where [(); S-1]: {
	pub molar_masses: [f64; S],
	pub reduced_molar_masses: [f64; S-1],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S>]>,
}

impl<const S: usize> System<S> where [(); S-1]:, [(); 2+S-1]: {
	pub fn dt(&self, pressure: f64, y: &[f64; 2+S-1]) -> ([f64; 2+S-1], [[f64; 2+S-1]; 2+S-1]) {
		use iter::into::Sum;
		let a = S-1;
		let Self{thermodynamics: species, reactions, molar_masses: W, reduced_molar_masses, ..} = self;
		let (T, V, amounts) = (y[0], y[1], y.suffix());
		let C = pressure / (ideal_gas_constant * T);
		let amount = V * C;
		let rcp_amount = 1. / amount;
		let concentrations : [_; S-1] = eval(amounts, |n| (n / V)/*.max(0.)*/); // Skips most abundant specie (last index) (will be deduced from conservation)
		let Ca = C - concentrations.sum():f64;
		let ref concentrations = from_iter(concentrations.chain([Ca]));
		let rcpC = 1. / C;
		let rcpV = 1. / V;
		let logP0_RT = f64::ln(NASA7::reference_pressure/ideal_gas_constant) - f64::ln(T);
		let ref H = eval(species, |s| s.dimensionless_specific_enthalpy(T));
		//let ref H_T = eval(species.prefix(), |s| s.dimensionless_specific_enthalpy_T(T));
		let ref H_T = eval(H.prefix(), |H| H/T);
		let ref G = eval!(species.prefix(), H_T; |s, h_T| h_T - s.dimensionless_specific_entropy(T)); // (H-TS)/RT
		let ref dT_G = eval!(species.prefix(); |s| s.dT_G(T));
		let ref log_concentrations = eval(concentrations, |&x| f64::ln(x));
		let ref mut dtω = [0.; S-1];
		let mut dTω = [0.; S-1];
		let mut dVω = [0.; S-1];
		let mut dnω = [[0.; S-1]; S-1];
		for Reaction{reactants, products, net, Σreactants, Σproducts, Σnet, rate_constant: rate_constant@RateConstant{temperature_exponent, activation_temperature, ..}, model, ..} in reactions.iter() {
			let log_kf = log_arrhenius(rate_constant, T);
			let c = model.efficiency(T, concentrations, log_kf);
			let mask = |mask, v| iter::zip!(mask, v).map(|(&mask, v):(_,&_)| if mask != 0. { *v } else { 0. });
			let Rf = f64::exp(reactants.dot(mask(reactants, log_concentrations)) + log_kf);
			let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
			let Rr = f64::exp(products.dot(mask(products, log_concentrations)) + log_kf - log_equilibrium_constant);
			//assert!(Rf.is_finite() && Rr.is_finite(), Rf, Rr, reactants, products, log_concentrations);
			let R = Rf - Rr;
			let cR = c * R;
			let νfRfνrRr = reactants.scale(Rf).sub(&products.scale(Rr));
			// dT(R) = (β+Ea/(T.Rideal))/T.R + Rr. Σ ν.dT(G) - νfRfνrRr[a]/T
			let dTR = (temperature_exponent+activation_temperature/T)/T*R + Rr*net.dot(dT_G) - νfRfνrRr[a] / T;
			// dV(R) = 1/V . ( (kf.Sf-kr.Sr) - (Σνf.Rf - Σνr.Rr) )
			let dVR = rcpV * ( νfRfνrRr[a] - (Σreactants*Rf - Σproducts*Rr));
			// dn(R) = 1/n . ( kf.(Sf-Sfa) - kr.(Sr-Sra) )
			let dnR = map(νfRfνrRr.prefix(), |νfRfνrRrj| rcp_amount * (νfRfνrRrj - νfRfνrRr[a]));
			let (dTc, dVc, dnc) = match model {
				Model::Elementary => (0., 0., [0.; S-1]),
				Model::ThreeBody{efficiencies}|Model::Falloff{efficiencies,..} => {
					let has = eval(efficiencies, |&e| if e != 0. { 1. } else { 0. });
					(
						// dT(c) = has(a) . -c/T
						has[a] * -c/T,
						// dV(c) = 1/V . ( has(a).c - Σefficiencies)
						rcpV * (has[a]*c  - efficiencies.sum():f64),
						// dn(c) = 1/V . ( has(n) - has(a) )
						eval(has.prefix(), |has_n| rcpV * (has_n - has[a]))
					)
				}
			};
			let RdTccdTR = R * dTc + c * dTR;
			let RdVccdVR = R * dVc + c * dVR;
			let RdnccdnR = eval!(dnc, dnR; |dnc,dnR| R*dnc + c*dnR);
			for (specie, ν) in net.enumerate() {
				// dtω = Σ ν c R
				dtω[specie] += ν * cR;
				// dT(ω̇̇̇̇̇̇̇̇̇̇) = Σ ν.(R.dT(c)+c.dT(R))
				dTω[specie] += ν * RdTccdTR;
				// dV(ω) = Σ ν.(R.dV(c)+c.dV(R))
				dVω[specie] += ν * RdVccdVR;
				// dn(ω) = Σ ν.(R.dn(c)+c.dn(R))
				for (dnω, RdnccdnR) in zip!(&mut dnω[specie], RdnccdnR) { *dnω += ν * RdnccdnR; }
			}
		}
		let ref dtω = *dtω;

		let Cp = eval(species, |s:&NASA7| s.dimensionless_specific_heat_capacity(T));
		let rcp_ΣCCp = 1./concentrations.dot(Cp);
		let dtT_T = - rcp_ΣCCp * dtω.dot(H_T); // R/RT
		let dtE = reduced_molar_masses.dot(dtω);
		let dtV = V * (dtT_T + T * ideal_gas_constant / pressure * dtE);
		let dtn = eval(dtω, |dtω| V*dtω);
		let mut J = [[f64::NAN; 2+S-1]; 2+S-1];
		// 1 / (Σ C.Cp)
		let rcp_ΣCCp = 1./concentrations.dot(Cp);
		let HaWa = H[a]/W[a];
		// Ha/Wa*W - H
		let HaWaWH = eval!(W.prefix(), H.prefix(); |W, H| HaWa*W - H);
		// dtT = - 1 / (Σ C.Cp) . Σ H.dtω
		let dtT = - rcp_ΣCCp * H.prefix().dot(dtω);

		let Cpa = Cp[a];
		let Cpa_Wa = Cpa/W[a];
		let concentrations = concentrations.prefix(); //::<{S-1}>();
		let Cp = Cp.prefix();
		// dT(dtT) = 1 / (Σ C.Cp) . [ dtT . Σ C.(Cpa/T - dT(Cp)) + Σ_ ( (Ha/Wa*W - H).dT(ω) + (Cpa/Wa.W - Cp).ω ) ]
		let dTdtT = rcp_ΣCCp * (dtT * concentrations.dot(species.prefix().map(|s:&NASA7| Cpa/T - s.dT_Cp(T))) + HaWaWH.dot(dTω) + map!(W.prefix(), Cp; |W, Cp| Cpa_Wa*W - Cp).dot(dtω));
		J[0][0] = dTdtT;
		// dV(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dV(ω) + dtT/V . Σ_ C.(Cp-Cpa) ]
		let dVdtT = rcp_ΣCCp * (HaWaWH.dot(dVω) + rcpV * dtT * concentrations.dot(Cp.map(|Cp| Cp - Cpa)));
		J[0][1] = dVdtT;
		// dn(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dn(ω) + dtT/V . (Cpa-Cp) ]
		let ref dndtT = eval!(dnω, Cp; |dnω, Cp| rcp_ΣCCp * (HaWaWH.dot(dnω) + rcpV * dtT * (Cpa-Cp) ));
		J[0][2..S+1].copy_from_slice(dndtT);

		// dT(dtV) = V/C . Σ_ (1-W/Wa).(dT(ω)+ω/T) + V/T.(dT(dtT) - dtT/T)
		let V_T = V / T;
		let dTdtV = V*rcpC* map!(reduced_molar_masses, dTω, dtω; |reduced_molar_mass, dTω, dtω| reduced_molar_mass*(dTω+dtω/T)).sum():f64 + V_T * (dTdtT - dtT/T);
		J[1][0] = dTdtV;
		// dV(dtn) = VdV(ω)+ω
		let dVdtn = from_iter(map!(dVω,dtω; |dVω,dtω| V*dVω+dtω));
		// dV(dtV) = 1/C . Σ_ (1-W/Wa).dV(dtn) + 1/T.(V.dV(dtT)+dtT)
		let dVdtV = rcpC * reduced_molar_masses.dot(dVdtn) + 1./T*(V*dVdtT+dtT);
		J[1][1] = dVdtV;
		// dn(dtn) = Vdn(ω)
		let dndtn = generate(|j| generate(|k| V*dnω[k][j])); // Transpose [k][j] -> [j][k]
		// dn(dtV) = 1/C . Σ_ (1-W/Wa).dn(dtn) + V/T.dn(dtT))
		let ref dndtV = eval!(dndtn, dndtT; |dndtn, dndtT| rcpC * reduced_molar_masses.dot(dndtn)+ V_T*dndtT);
		J[1][2..S+1].copy_from_slice(dndtV);

		// dT(dtn) = VdT(ω)
		let dTdtn = dTω.scale(V);
		for (i, dTdtn) in dTdtn.enumerate() { J[2+i][0] = dTdtn; }
		// dV(dtn)
		for (i, dVdtn) in dVdtn.enumerate() { J[2+i][1] = dVdtn; }
		// dn(dtn)
		for j in 0..S-1 { for (i, dndtn) in dndtn[j].enumerate() { J[2+i][j] = dndtn; } }
		(from_iter([dtT_T*T, dtV].chain(dtn)), J)
	}
}

#[derive(Clone)] pub struct State<const S: usize> where [(); S-1]: {
	pub temperature: f64,
	//pub volume: f64,
	pub amounts: [f64; S-1]
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species: [&'t str; S],
	pub system: System<S>,
	pub time_step: f64,
	pub pressure: f64,
	pub volume: f64,
	pub state: State<S>
}

use std::lazy::SyncLazy;
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl<const S: usize> Simulation<'t, S> where [(); S-1]: {
	pub const species_len : usize = S;

	pub fn new(system: &'b [u8]) -> ::ron::Result<Self> where 'b: 't, [(); S]: {
		let ron::System{species: species_data, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
		let ron::Phase::IdealGas{species, state, ..} = phases.into_vec().into_iter().next().unwrap();
		use std::convert::TryInto;
		let species : [_; S] = species.as_ref().try_into().unwrap();
		let molar_masses = eval(species, |s| species_data[s].composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let reduced_molar_masses = eval(molar_masses.prefix(), |w:&f64| 1.-w/molar_masses[S-1]);
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
					Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe},
				}},
			}
		}));

		let ron::InitialState{temperature, pressure, mole_proportions, volume} = state;
		let amount = pressure * volume / (ideal_gas_constant * temperature);
		let mole_proportions = eval(&species, |specie| *mole_proportions.get(specie).unwrap_or(&0.));
		let amounts = eval(mole_proportions.prefix(), |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);

		Ok(Self{
			species,
			system: System{molar_masses, reduced_molar_masses, thermodynamics, reactions},
			time_step, pressure, volume,
			state: State{temperature, /*volume,*/ amounts}
		})
	}
}
