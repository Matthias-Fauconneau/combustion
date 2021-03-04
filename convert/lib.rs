#![feature(non_ascii_idents, type_ascription, once_cell, array_map)]
#![allow(confusable_idents, non_upper_case_globals, non_snake_case)]
use std::ops::Deref;
use linear_map::LinearMap as Map;
use system::K;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)
use model::{/*Map,*/ Element, Troe};
use system::{RateConstant, NASA7};

use std::{convert::TryInto, lazy::SyncLazy};
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, N: 14.0067, O: 15.999, Ar: 39.95}").unwrap()
	.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

#[derive(Debug)] struct Species {
	pub molar_mass: Box<[f64]>,
	pub thermodynamics: Box<[NASA7]>,
	diameter: Box<[f64]>,
	well_depth_J: Box<[f64]>,
	polarizability: Box<[f64]>,
	permanent_dipole_moment: Box<[f64]>,
	rotational_relaxation: Box<[f64]>,
	internal_degrees_of_freedom: Box<[f64]>,
	heat_capacity_ratio: Box<[f64]>,
}

#[derive(Debug)] enum Model {
	Elementary,
	ThreeBody { efficiencies: Box<[f64]> },
	PressureModification { efficiencies: Box<[f64]>, k0: RateConstant },
	Falloff { efficiencies: Box<[f64]>, k0: RateConstant, troe: Troe },
}

impl Model {
	pub fn efficiency(&self) -> Tok fn(T: f64, concentrations: &[f64; S], log_k_inf: f64) -> f64 {
	match self {
		Self::Elementary => 1.,
		Self::ThreeBody{efficiencies} => efficiencies.dot(concentrations),
		Self::PressureModification{efficiencies, k0} => {
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(/*const*/ *k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
			Pr / (1.+Pr)
		}
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(/*const*/ *k0, T) - log_k_inf); // [k0/kinf] = [1/C] (m3/mol)
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

#[derive(Debug)] struct Reaction {
	pub reactants: Box<[f64]>,
	pub products: Box<[f64]>,
	pub net: Box<[f64/*; S-1*/]>,
	pub Σreactants: f64,
	pub Σproducts: f64,
	pub Σnet: f64,
	pub rate_constant: RateConstant,
	pub model: Model,
}

#[derive(Debug)] struct System {
	pub species: Species,
	pub reactions: Box<[Reaction]>,
	//pub transport_polynomials: TransportPolynomials_,
}

impl System {
fn new(model::Model{species, reactions, ..}: &model::Model) {
	let species: Box<[_]> = (species.into():Vec<_>).into();
	let species = species.deref();
	pub fn eval<T, U>(v: impl IntoIterator<Item=T>, f: impl Fn(T)->U) -> Box<[U]> { v.into_iter().map(f).collect() }
	let species = eval(species, |(k,specie):&(_,model::Specie)| {
		let mut specie = specie.clone();
		for T in specie.thermodynamic.temperature_ranges.iter_mut() { *T *= K; } // K->J
		for piece in specie.thermodynamic.pieces.iter_mut() { for (order, a) in piece[0..6].iter_mut().enumerate() { *a /= f64::powi(K, (1+order) as i32); } } // /K^n->/J^n
		(k, specie)
	});
	let species = species.deref();
	let molar_mass = eval(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
	let thermodynamics = eval(species, |(_, model::Specie{thermodynamic: model::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
		[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
		[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
		ref ranges => panic!("{:?}", ranges),
	});
	let diameter = eval(species, |(_,s)| s.transport.diameter_Å*1e-10);
	let well_depth_J = eval(species, |(_,s)| s.transport.well_depth_K * K);
	use model::Geometry::*;
	let polarizability = eval(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
	let permanent_dipole_moment = eval(species, |(_,s)|
		if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
	let rotational_relaxation = eval(species, |(_,s)| if let Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
	let internal_degrees_of_freedom = eval(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
	let heat_capacity_ratio = eval(species, |(_,s)| {
		let f = match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. };
		1. + 2./f
	});
	let species_names = eval(species, |(name,_)| *name);
	let species_names = species_names.deref();
	let species_composition = eval(species, |(_,s)| &s.composition);
	let species = Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, heat_capacity_ratio};
	//let transport_polynomials = species.transport_polynomials();
	let reactions = eval(reactions.into():Vec<_>, |model::Reaction{ref equation, rate_constant, model}| {
		for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
		let [reactants, products] = iter::vec::eval(equation, |e| eval(species_names, |&s| *e.get(s).unwrap_or(&0) as f64));
		let net = eval(products.into_iter().zip(reactants.into_iter()).take(reactants.len()-1), |(a, b)| a-b);
		{let mut net_composition = Map::new();
			for (s, &ν) in net.into_iter().enumerate() {
				for (element, &count) in species_composition[s] {
					if !net_composition.contains_key(&element) { net_composition.insert(element, 0); }
					*net_composition.get_mut(&element).unwrap() += ν as i8 * count as i8;
				}
			}
			for (_, &ν) in &net_composition { assert!(ν == 0, "{:?} {:?}", net_composition, equation); }}
		let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
		let Σnet = Σproducts-Σreactants;
		let from = |efficiencies:Map<_,_>| eval(species_names, |&specie| *efficiencies.get(specie).unwrap_or(&1.));
		Reaction{
			reactants, products, net, Σreactants, Σproducts, Σnet,
			rate_constant: rate_constant.into(),
			model: {use model::ReactionModel::*; match model {
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: from(efficiencies)},
				PressureModification{efficiencies, k0} => Model::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
				Falloff{efficiencies, k0, troe: Troe{A,T3,T1,T2}} => Model::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: Troe{A,T3:T3*K,T1:T1*K,T2:T2*K}},
			}},
		}
	});
	System{species, /*transport_polynomials,*/ reactions}
}

#[proc_macro]
pub fn system(model: proc_macro::TokenStream) -> proc_macro::TokenStream {
	let model = syn::parse_macro_input!(model as syn::LitStr).value();
	let model = std::fs::read(&model).expect(&format!("{:?}", (std::env::current_dir(), model)));
	let system = System::new(&::ron::de::from_bytes(&model).unwrap());
	let system: proc_macro2::TokenStream = format!("{:?}", system).parse().unwrap();
	quote::quote!({use {model::Troe, system::{NASA7, RateConstant}, crate::{Species, reaction::Model::*}}; #system}).into()
}

#[proc_macro]
pub fn dtω(model: proc_macro::TokenStream) -> proc_macro::TokenStream {
	let model = syn::parse_macro_input!(model as syn::LitStr).value();
	let model = std::fs::read(&model).expect(&format!("{:?}", (std::env::current_dir(), model)));
	let system = System::new(&::ron::de::from_bytes(&model).unwrap());
	let reactions = proc_macro::TokenStream::new();
	for Reaction{reactants, products, net, /*Σreactants, Σproducts,*/ Σnet, rate_constant/*: rate_constant@RateConstant{temperature_exponent, activation_temperature, ..}*/, model, ..} in system.reactions.iter() {
		//let rate_constant: proc_macro2::TokenStream = format!("{:?}", rate_constant).parse().unwrap();
		let literal = |value| -> proc_macro2::TokenStream { format!("{:?}", value).parse().unwrap() }
		reactions.append(quote!{
			let log_kf = log_arrhenius(/*const*/ #literal(rate_constant), T);
			let c = model.efficiency(T, concentrations, log_kf);
		//fn dot<const N: usize, const a: [f64; N]>(b: [f64; N]) { let sum = 0.; unroll! { for index in 0..N { if a[index] != 0. { sum += a[index]*b[index] } } sum } }
		macro_rules! dot { ($a:ident, $b:ident) => ({ let sum = 0.; debug_assert_eq!(S, 53); unroll! { for index in 0..53 { if $a[index] != 0. { sum += $a[index]*$b[index] } } } sum }) }
		let Rf = f64::exp(dot!(reactants, log_concentrations) + log_kf);
		let log_equilibrium_constant = -net.dot(G_RT) + Σnet*logP0_RT;
		let Rr = f64::exp(dot!(products, log_concentrations) + log_kf - log_equilibrium_constant);
		//assert!(Rf.is_finite() && Rr.is_finite(), "{:?}", (Rf, Rr, reactants, products, log_concentrations));
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

		/*for (specie, ν) in net.enumerate() {
			// dtω = Σ ν c R
			dtω[specie] += ν * cR;
			/*// dT(ω̇̇̇̇̇̇̇̇̇̇) = Σ ν.(R.dT(c)+c.dT(R))
			dTω[specie] += ν * RdTccdTR;
			// dV(ω) = Σ ν.(R.dV(c)+c.dV(R))
			dVω[specie] += ν * RdVccdVR;
			// dn(ω) = Σ ν.(R.dn(c)+c.dn(R))
			for (dnω, RdnccdnR) in zip!(&mut dnω[specie], RdnccdnR) { *dnω += ν * RdnccdnR; }*/
		}*/
		debug_assert_eq!(S-1, 52); unroll! { for specie in 0..52 { if net[specie] != 0. { dtω[specie] += net[specie] * cR; } } }
	}}
	quote!{
		fn dtω<const S: usize>(T: f64, log_concentrations: [f64; S-1]) {
			let ref mut dtω = [0.; S-1];
			/*let mut dTω = [0.; S-1];
			let mut dVω = [0.; S-1];
			let mut dnω = [[0.; S-1]; S-1];*/
			#reactions
			*dtω
		}
	}
}
