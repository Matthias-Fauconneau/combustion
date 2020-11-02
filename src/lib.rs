#![feature(in_band_lifetimes)]
#![feature(trait_alias)]
#![feature(box_syntax)]
#![feature(map_into_keys_values)]
#![allow(non_snake_case)]
#![feature(associated_type_bounds)]
#![feature(bindings_after_at)]

pub fn scale(s: f64, v: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { v.into_iter().map(move |v| s*v) }
pub fn recip(x: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { x.into_iter().map(|x| f64::recip(x)) }
pub fn mul(a: impl IntoIterator<Item=f64,IntoIter:'t>, b: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { a.into_iter().zip(b.into_iter()).map(|(a,b)| a*b) }
pub fn dot(a: &[f64], b: &[f64]) { mul(a.iter().copied(), b.iter().copied()).sum() }
pub fn dotia(a: impl IntoIterator<Item=f64>, b: &[f64]) -> f64 { mul(a, b.iter().copied()).sum() }
pub fn dotai(a: &[f64], b: impl IntoIterator<Item=f64>) -> f64 { mul(a.iter().copied(), b).sum() }

mod ron {
use serde::Deserialize;
pub use std::collections::BTreeMap as Map;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] pub enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] pub struct InitialState<'t> { pub temperature: f64, pub pressure: f64, #[serde(borrow)] pub mole_proportions: Map<&'t str, f64>, pub volume: f64 }
#[derive(Deserialize, Debug)] pub enum Phase<'t> {
	IdealGas {
		elements: Box<[Element]>,
		species: Box<[&'t str]>,
		#[serde(borrow)] state: InitialState<'t>,
	}
}
#[derive(Deserialize, Debug)] pub struct NASA7 {
	pub temperature_ranges: Box<[f64]>,
	pub coefficients: Box<[Box<[f64]>]>,
}

#[derive(Deserialize, Debug)] enum Transport {
	Atom { well_depth: f64, diameter: f64},
	Linear { well_depth: f64, diameter: f64, polarizability: f64, rotational_relaxation: f64},
	Nonlinear { well_depth: f64, diameter: f64, rotational_relaxation: f64},
}
#[derive(Deserialize, Debug)] pub struct Specie {
	pub composition: Map<Element, u8>,
	pub thermodynamic: NASA7,
	transport: Transport
}

#[derive(Deserialize, Debug)] pub struct RateConstant {
	#[serde(rename="A")] pub preexponential_factor: f64,
	#[serde(rename="β")] pub temperature_exponent: f64,
	#[serde(rename="Ea")] pub activation_energy: f64
}

#[derive(Deserialize, Debug)] pub struct Troe { pub A: f64, pub T3: f64, pub T1: f64, pub T2: f64 }

#[derive(Deserialize, Debug)] pub enum Model<'t> {
	Elementary,
	ThreeBody { #[serde(borrow)] efficiencies: Map<&'t str, f64> },
	Falloff { #[serde(borrow)] efficiencies: Map<&'t str, f64>, k0: RateConstant, troe: Troe },
}

#[derive(Deserialize, Debug)] pub struct Reaction<'t> {
	#[serde(borrow)] pub equation: [Map<&'t str, u8>; 2],
	pub rate_constant: RateConstant,
	pub model: Model<'t>,
}

#[derive(Deserialize, Debug)] pub struct System<'t> {
	pub time_step: f64,
	#[serde(borrow)] pub phases: Box<[Phase<'t>]>,
	#[serde(borrow)] pub species: Map<&'t str, Specie>,
	#[serde(borrow)] pub reactions: Box<[Reaction<'t>]>,
}
}

#[allow(non_upper_case_globals)] pub const ideal_gas_constant : f64 = 8.31446261815324;

use self::ron::*;

impl NASA7 {
	pub fn specific_heat_capacity(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T)
	}
	pub fn specific_enthalpy(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T)
	}
	pub fn specific_entropy(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T)
	}
	pub fn b(&self, T: f64) -> f64 { // S/R - H/RT
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		a[6] - a[0] + (a[0]-1.)*f64::ln(T) + a[1]/2.*T + a[2]/6.*T*T + a[3]/12.*T*T*T + a[4]/20.*T*T*T*T - a[5]/T
	}
	pub fn dT_Cp(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[1]+2*a[2]*T+3*a[3]*T*T+4*a[4]*T*T*T)
	}
	pub fn dT_b(&self, T: f64) -> f64 { // dT(S/R - H/RT)
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		(a[0]-1.)/T + a[1]/2. + a[2]/12.*T + a[3]/36.*T*T + a[4]/80.*T*T*T + a[5]/(T*T)
	}
}

pub fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &RateConstant, temperature: f64) -> f64 {
	preexponential_factor*temperature.powf(temperature_exponent)*f64::exp(-activation_energy/(ideal_gas_constant/4.184*temperature))
}

#[derive(Debug)] pub<N: usize> enum Model<N> {
	Elementary,
	ThreeBody { efficiencies: [f64; N] },
	Falloff { efficiencies: [f64; N], k0: RateConstant, troe: Troe },
}

impl<N: usize> Model<N> {
fn efficiency(&self, T: f64, concentrations: &[f64; N], k_inf: f64) -> f64 {
	match self {
		Self::Elementary => 1.,
		Self::ThreeBody{efficiencies} => dot(efficiencies, concentrations),
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = dot(efficiencies, concentrations) * arrhenius(k0, T) / kf;
			let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
			let log10Fcent = f64::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			/*let N = 0.75-1.27*log10Fcent;
			let f1 = (f64::log10(Pr) + C)/(N-0.14*(f64::log10(Pr) + C));
			let F = num::exp10(log10Fcent/(1.+f1*f1));*/
			let ATroe = f64::log10(Pr) + C;
			let BTroe = 0.806 - 1.1762*log10Fcent - 0.14*f64::log10(Pr);
			let F = f64::powf(Fcent, f64::recip(1.+num::sq(ATroe/BTroe)));
			Pr * F / (1.+Pr); // Chemically activated bimolecular reaction
		}
	}
}
}

pub<N: usize> struct Reaction<N> {
	pub equation: [(Box<[usize]>, Box<[u8]>); 2],
	rate_constant: RateConstant,
	pub model: Model<N>,
	specie_net_coefficients: [f64; N-1]
	Σνf: f64,
	Σνr: f64,
	PRν: f64,
}

pub<N: usize> struct System<N> {
	pub molar_masses: [f64; N],
	pub thermodynamics: [NASA7; N],
	pub reactions: Box<[Reaction<N>]>,
	pub amount: f64,
	time_step: f64,
	pub pressure: f64,
	pub volume: f64,
}

#[derive(Clone)] pub<N: usize> struct State<N> {
	pub temperature: f64,
	pub amounts: [f64; N]
}

pub<N: usize> struct Simulation<'t, N> {
	pub species: Box<[&'t str]>,
	pub system: System<N>,
	pub state: State<N>
}

impl<N: usize> Simulation<'t, N> {
#[fehler::throws(anyhow::Error)] pub fn new(system: &'b [u8]) -> Self where 'b: 't {
	let ron::System{species, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;

	let standard_atomic_weights : Map<Element, f64> = ::ron::de::from_str("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}")?;
	let standard_atomic_weights : Map<_, f64> = standard_atomic_weights.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect();

	use iter::from_iter;
	let specie_names = from_iter(species.keys().copied());
	let molar_masses = from_iter(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum()));
	let thermodynamics = from_iter(species.into_values().map(|Specie{thermodynamic,..}| thermodynamic));
	let species = specie_names;

	let reactions = from_iter(Vec::from(reactions).into_iter().map(|self::ron::Reaction{equation, rate_constant, model}| {
		let atmospheric_pressure = 101325.;
		let equation = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| (e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key)).collect(), e.values().copied().collect()))),
		let specie_net_coefficients =
			from_iter(species.iter().map(|specie| {let [reactant, product] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| *e.get(specie).unwrap_or(&0) as i8)); (product-reactant) as f64})),
		let [(_, νf), (_, νr)] = equation;
		let [Σνf, Σνr] = [νf.iter().sum(), νr.iter().sum();
		let PRν = f64::powf(atmospheric_pressure / ideal_gas_constant, specie_net_coefficients.iter().sum());
		Reaction{
			equation: iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| (e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key)).collect(), e.values().copied().collect()))),
			rate_constant: rate_constant
			model: {use self::ron::Model::*; match model {
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)).collect()},
				Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)).collect(), k0, troe},
			}},
			specie_net_coefficients,
			Σνf, Σνr, PRν
	}));

	let Phase::IdealGas{state, ..} = Vec::from(phases).into_iter().next().unwrap();
	let InitialState{temperature, pressure, mole_proportions, volume} = state;
	let mole_proportions = from_iter(species.iter().map(|specie| *mole_proportions.get(specie).unwrap_or(&0.)));
	let amount = pressure * volume / (ideal_gas_constant * temperature);
	let amounts = from_iter(scale(amount/mole_proportions.iter().sum::<f64>(), mole_proportions.iter().copied()));

	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions, amount, time_step, pressure, volume},
		state: State{temperature, amounts}
	}
}
}

impl<N: usize> State<N> {
	pub fn step(&mut self, System{thermodynamics: species, reactions, amount, volume: V}: &System<N>) {
		use iter::array::{from_iter, map, generate};
		let scale = |s, v| from_iter(scale(s, v.iter().copied()));
		macro_rules! map2 { ($a:ident, $b:ident, $expr:expr) => { $a.iter().zip($b.iter()).map(|($a,$b)| $expr) } }
		macro_rules! map3 { ($a:ident, $b:ident, $c:ident, $expr:expr) => { $a.iter().zip($b.iter()).zip($c.iter()).map(|(($a,$b),$c)| $expr) } }

		let rcpn = 1. / amount;
		let T = self.temperature;
		let B = map(species, |s| s.b(T));
		let dTB = map(species, |s| s.dT_b(T));
		let rcpV = 1. / V;
		// C = n / V
		let C = scale(rcpV, self.amounts);
		let a = N-1;
		let Ca = C[a];
		let mut ω = vec!(0.; N-1); // Skips most abundant specie (last index) (will be deduced from conservation)
		let mut ∂Tω = vec!(0.; N-1);
		let mut ∂Vω = vec!(0.; N-1);
		let mut ∂nω = vec!(vec!(0.; N-1); N-1); //[k][j]
		for Reaction{equation, rate_constant: rate_constant@RateConstant{temperature_exponent: β, activation_energy: Ea}, specie_net_coefficients: ν, Σνf, Σνr, PRν} in reactions.iter() {
			let equilibrium_constant = PRν * f64::exp(dot(ν, B));
			let kf = arrhenius(rate_constant, T);
			let kr = kf / equilibrium_constant;
			// ΠC^ν
			let [ΠCνf, ΠCνr] = from_iter(equation.iter().map(|(species, ν)| species.iter().zip(ν.iter()).map(|(&specie, &ν)| C[specie].powi(ν as i32)).product::<f64>() ));
			let Rf = kf * ΠCνf;
			let Rr = kr * ΠCνr;
			let νfRfνrRr = vec!(0.; N);
			let [forward, reverse] = equation;
			for (specie, ν) in forward { νfRfνrRr[specie] += ν * Rf; }
			for (specie, ν) in reverse { νfRfνrRr[specie] -= ν * Rr; }
			let c = model.efficiency(T, &C, kf);
			let R = Rf - Rr;
			// ∂T(R) = (β+Ea/(T.Rideal))/T.R + Rr. Σ ν.dT(B) - νfRfνrRr[a]/T
			let ∂TR = (β+Ea/(T*ideal_gas_constant))/T*R + Rr*dot(ν,dTB) - νfRfνrRr[a] / T
			// ∂V(R) = 1/V . ( (kf.Sf-kr.Sr) - (Σνf.Rf - Σνr.Rr) )
			let ∂VR = rcpV * ( νfRfνrRr[a] - (Σνf*Rf - Σνr*Rr));
			// ∂n(R) = 1/n . ( kf.(Sf-Sfa) - kr.(Sr-Sra) )
			let ∂nR = map(νfRfνrRr, |νfRfνrRrj| rcpn * (νfRfνrRrj - νfRfνrRr[a]));
			let (∂Tc, ∂Vc, ∂nc) = match model {
				Model::Elementary => (0., 0., vec!(0.; N-1)),
				Model::ThirdBody{efficiencies}|Model::Falloff{efficiencies} => {
					let has = map(efficiencies, |e| e != 0. { 1. } else { 0. });
					(
						// ∂T(c) = has(a) . (-C/T)
						has[a] * -C/T
						// ∂V(c) = 1/V . ( has(a). C - Σ C )
						rcpV * (has[a] * C  - dot(has, C)),
						// ∂n(c) = 1/V . ( has(n) - has(a) )
						has[..a].map(|has_n| rcpV * (has_n - has[a]));
					)
				}
			};
			let cR = c * R;
			let R∂Tcc∂TR = R * ∂Tc + c * ∂TR;
			let R∂Vcc∂VR = R * ∂Vc + c * ∂VR;
			let R∂ncc∂nR = from_iter(map2!(∂nc,∂nR, R*∂nc + c*∂nR));
			for (specie, ν) in ν.iter().enumerate() {
				//ω̇̇̇̇̇̇̇̇̇̇ = Σ ν c R
				ω[specie] += ν * cR;
				// ∂T(ω̇̇̇̇̇̇̇̇̇̇) = Σ ν.(R.∂T(c)+c.∂T(R))
				∂Tω[specie] += ν * R∂Tcc∂TR;
				// ∂V(ω) = Σ ν.(R.∂V(c)+c.∂V(R))
				∂Vω[specie] += ν * R∂Vcc∂VR;
				// ∂n(ω) = Σ ν.(R.∂n(c)+c.∂n(R))
				for ∂nω in ∂nω[specie] { ∂nω += ν * R∂ncc∂nR; }
			}
		}
		let mut J = unsafe{nalgebra::MatrixMN::<f64, 2+N-1, 2+N-1>::new_uninitialized()}; // fixme
		let Cp = map(species, |s| s.specific_heat_capacity(T));
		// 1 / (Σ C.Cp)
		let rcp_ΣCCp = 1./dot(C, Cp);
		let H = species.map(|s| s.specific_enthalpy(T));
		let W = system.molar_masses;
		let Wa = W[a];
		let HaWa = H[a]/Wa;
		// Ha/Wa*W - H
		let HaWaWH = from_iter(map2!(W,H, HaWa*W - H));
		// dtT = - 1 / (Σ C.Cp) . Σ H.ω̇
		let dtT = - rcp_ΣCCp * dot(H, ω̇);
		let CpaWa = Cp[a]/Wa;
		// ∂T(dtT) = 1 / (Σ C.Cp) . [ dtT . Σ C.(Cpa/T - ∂T(Cp)) + Σ_ ( (Ha/Wa*W - H).∂T(ω) + (Cpa/Wa.W - Cp).ω ) ]
		let ∂TdtT = rcp_ΣCCp * ( dtT * dot(C, species.iter().map(|s| Cp[a]/T - s.dt_Cp())) + dot(HaWaWH, ∂Tω̇̇̇̇̇̇̇̇̇̇) + dotia(map2!(W,Cp, CpaWa*W - Cp), ω̇))
		J[(0,0)] = ∂TdtT;
		// ∂V(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).∂V(ω) + dtT/V . Σ_ C.(Cp-Cpa) ]
		let ∂VdtT = rcp_ΣCCp * ( dot(HaWaWH, ∂Vω̇̇̇̇̇̇̇̇̇̇) + rcpV * dtT * dotai(C, Cp.iter().map(|Cp| Cp - Cpa)))
		J[(0,1)] = ∂VdtT;
		// ∂n(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).∂n(ω) + dtT/V . (Cpa-Cp) ]
		let ∂ndtT = from_iter(map2!(∂nω, Cp, rcp_ΣCCp * ( dot(HaWaWH, ∂nω) + rcpV * dtT . (Cpa-Cp) )));
		J.row_part_mut(0,2,N+1).copy_from_slice(∂ndtT);

		// ∂T(dtV) = V/C . Σ_ (1-W/Wa).(∂T(ω)+ω/T) + V/T.(∂T(dtT) - dtT/T)
		let rcpC = V*rcpn;
		let WWa = map(W, |W| 1.-W/Wa);
		let VT = V/T;
		let ∂TdtV = V*rcpC * map3!(WWa,∂Tω̇,ω̇, WWa*(∂Tω+ω/T)) + VT * (∂TdtT - dtT/T);
		J[(1,0)] = ∂TdtV;
		// ∂V(dtn) = V∂V(ω)+ω
		let ∂Vdtn = from_iter(map2!(∂Vω,ω, V*∂Vω+ω));
		// ∂V(dtV) = 1/C . Σ_ (1-W/Wa).∂V(dtn) + 1/T.(V.∂V(dtT)+dtT)
		let ∂VdtV = rcpC * dot(WWa, ∂Vdtn) + 1./T*(V*∂VdtT+dtT);
		J[(1,1)] = ∂VdtV;
		// ∂n(dtn) = V∂n(ω)
		let ∂ndtn = generate(N-1, |j| generate(N-1, |k| V*∂nω[k][j])); // Transpose [k][j] -> [j][k]
		// ∂n(dtV) = 1/C . Σ_ (1-W/Wa).∂n(dtn) + V/T.∂n(dtT))
		let ∂ndtV = map2!(∂ndtn, ∂ndtT, |∂ndtn,∂ndtT| rcpC * dot(WWa, ∂ndtn)+ VT*∂ndtT);
		J.row_part_mut(1,2,N+1) = Matrix::from_iterator(∂ndtV);

		// ∂T(dtn) = V∂T(ω)
		let ∂Tdtn = scale(V, ∂Tω);
		J.column_part_mut(0,2,N+1) = Matrix::from_iterator(∂Tdtn);
		// ∂V(dtn)
		J.column_part_mut(1,2,N+1) = Matrix::from_iterator(∂Vdtn);
		// ∂n(dtn)
		for j in 2..N+1 { J.column_part_mut(j,2,N+1).copy_from_slice(∂ndtn[j]);

		// dt(V) = V.(TR/P.Σ{(1-Wk/Wa).ω}+1/T.dt(T))
		// dt(T) = - Σ_ ((H-(W/Wa).Ha).w) / ( Ct.Cpa + Σ_ (Cp_Cpa)C )
		// dt(n) = Vω
		for (state, rate) in self.amounts[..specie_count-1].iter_mut().zip(production_rates.iter()) { *state = 0f64.max(*state + system.time_step * system.volume * rate); }
		let total_amount = system.pressure * system.volume / (ideal_gas_constant * self.temperature);
		self.amounts[specie_count-1] = 0f64.max(total_amount - self.amounts[..specie_count-1].iter().sum::<f64>()); // Enforces mass conservation constraint by rescaling most abundant specie (last index)
		for &amount in self.amounts.iter() { assert!(amount>= 0. && amount< total_amount,"{} {:?} {}", amount, &self.amounts, total_amount); }
		self.temperature += system.time_step * dtT;
	}
}

impl From<State> for Box<[Box<[f64]>]> { fn from(s: State) -> Self { box [box [s.temperature] as Box<[_]>, s.amounts] as Box<[_]> } }
