#![allow(incomplete_features)]
#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![feature(trait_alias)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals,unused_imports,uncommon_codepoints)]
pub mod ron;
mod collision_integrals; // Reduced collision integrals table computed from Chapman-Enskog theory with Stockmayer potential by L. Monchick and E.A. Mason. Transport properties of polar gases. J. Chem. Phys.
use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
use {iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}}, self::ron::{Map, Element, Troe}};

fn quadratic_interpolation(x: [f64; 3], y: [f64; 3], x0: f64) -> f64 { ((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1] }

fn eval_poly<const N: usize>(P: &[f64; N], x: f64) -> f64 { P.dot(generate(|k| x.powi(k as i32))) }

trait Vector<const N: usize> = vec::Vector<N>+iter::IntoIterator<Item=f64>;
 fn weighted_polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>, w: impl Vector<N>) -> [f64; D] {
	use nalgebra::{DMatrix, DVector, SVD};
	let w = DVector::from_iterator(N, w.into_iter());
	let A = DMatrix::from_iterator(N, D, x.into_iter().zip(w.iter()).map(|(x, w)| (0..D).map(move |k| w*x.powi(k as i32))).flatten());
	let b = DVector::from_iterator(N, y.into_iter().zip(w.iter()).map(|(x, w)| w*x));
	use std::convert::TryInto;
	SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap()
}
// Regression with 1/y² weights (towards relative vertical error)
fn polynomial_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Copy) -> [f64; D] { weighted_polynomial_regression(x, y, map(y, |y| 1./sq(y))) }
fn polynomial_fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { polynomial_regression(map(t, x), map(t, y)) }

pub const kB : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const light_speed : f64 = 299_792_458.;
const μ0 : f64 = 1.2566370621e-6;
const ε0 : f64 = 1./(light_speed*light_speed*μ0);

#[derive(Debug)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure_R : f64 = 101325. / (kB*NA); // 1 atm
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R
	//fn dT_specific_heat_capacity(&self, T: f64) -> f64 { 	let a = self.a(T); kB*NA * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T) } // /R
	//fn dT_Gibbs_free_energy(&self, T: f64) -> f64 { let a = self.a(T); (1.-a[0])/T - a[1]/2. - a[2]/12.*T - a[3]/36.*T*T - a[4]/80.*T*T*T - a[5]/(T*T) } // dT((H-TS)/RT)
}

#[derive(Debug, Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

impl From<ron::RateConstant> for RateConstant { fn from(ron::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: ron::RateConstant) -> Self {
	const J_per_cal: f64 = 4.184;
	Self{log_preexponential_factor: log(preexponential_factor), temperature_exponent, activation_temperature: activation_energy*J_per_cal/(kB*NA)}
}}

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
	//mass: f64,
	pub molar_mass: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S>]>, // .net[S-1]
	diameter: [f64; S],
	well_depth: [f64; S],
	polarizability: [f64; S],
	dipole: [f64; S],
	//rotational_relaxation: [f64; S],
	//internal_degrees_of_freedom: [f64; S],
}
impl<const S: usize> System<S> where [(); S-1]: {
	const volume : f64 = 1.;
}

impl<const S: usize> System<S> where [(); S-1]:, [(); 1+S-1]: {
	pub fn dt_J(&self, pressure_R: f64, y: &[f64; 1+S-1]) -> ([f64; 1+S-1], /*[[f64; 1+S-1]; 1+S-1]*/) {
		use iter::into::{IntoMap, Sum};
		//let a = S-1;
		let Self{thermodynamics, reactions/*, molar_masses: W*/, ..} = self;
		//let rcpV = 1. / V;
		let (T, amounts) = (y[0], y.suffix());
		let C = pressure_R / T;
		//let rcp_C = 1. / C;
		//let rcp_amount = rcpV * rcp_C;
		let logP0_RT = log(NASA7::reference_pressure_R) - log(T);
		let ref H = eval(thermodynamics, |s| s.specific_enthalpy(T));
		let ref H_T = eval(H.prefix(), |H| H/T);
		let ref G = eval!(thermodynamics.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
		//let ref dT_G = eval!(thermodynamics.prefix(); |s| s.dT_Gibbs_free_energy(T));
		let concentrations : [_; S-1] = eval(amounts, |&n| n/*.max(0.)*/ / Self::volume); // Skips most abundant specie (last index) (will be deduced from conservation)
		let Ca = C - Sum::<f64>::sum(concentrations);
		let ref concentrations = from_iter(concentrations.chain([Ca]));
		let ref log_concentrations = eval(concentrations, |&x| log(x));
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

		let Cp = eval(thermodynamics, |s:&NASA7| s.specific_heat_capacity(T));
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
		let dTdtT = rcp_ΣCCp * (dtT * concentrations.dot(thermodynamics.map(|s:&NASA7| Cpa/T - s.dT_specific_heat_capacity(T))) + HaWaWH.dot(dTω) + map!(W.prefix(), Cp; |W, Cp| Cpa_Wa*W - Cp).dot(dtω));
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

//struct Transport<const S: usize> { D: [f64; S], η: f64, λ: f64 }
impl<const S: usize> System<S> where [(); S-1]: {
pub fn transport(&self, _pressure: f64, temperature: f64, amounts: &[f64; S]) -> f64 {//Transport<S> {
	let (header_log_T⃰, [Ω⃰22,_A⃰ ,_B⃰,_C⃰]) = {use collision_integrals::*; (eval(header_T⃰.copied(), log),
		// Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
		 [&Ω⃰22,&A⃰ ,&B⃰,&C⃰].map(|table| table.map(|T⃰_row| polynomial_regression(header_δ⃰.copied(), T⃰_row)))
	)};
	let Self{molar_mass, thermodynamics: _, diameter, well_depth, polarizability, dipole, /*rotational_relaxation, internal_degrees_of_freedom,*/ ..} = self;
	let χ = |a, b| { // Corrections to the effective diameter and well depth to account for interaction between a polar and a non-polar molecule
		if dipole[a] == dipole[b] { 1. } else {
			let (polar, non_polar) = if dipole[a] != 0. { (a,b) } else { (b,a) };
			1. + polarizability[non_polar]/cb(diameter[non_polar]) * sq(dipole[polar]/sqrt(4.*π*ε0*cb(diameter[polar]) * well_depth[polar]))
					 * sqrt(well_depth[polar]/well_depth[non_polar]) / 4.
		}
	};
	let reduced_well_depth = |a, b| sqrt(well_depth[a]*well_depth[b]) * sq(χ(a, b));
	let T⃰ = |a, b, T| kB * T / reduced_well_depth(a, b);
	let reduced_dipole_moment = |a, b| dipole[a]*dipole[b] / (8. * π * ε0 * sqrt(well_depth[a]*well_depth[b]) * cb((diameter[a] + diameter[b])/2.)); // ̃δ⃰
	let collision_integral = |table : &[[f64; 8]; 39], a, b, T| {
		let log_T⃰ = log(T⃰ (a, a, T));
		assert!(*header_log_T⃰ .first().unwrap() <= log_T⃰  && log_T⃰  <= *header_log_T⃰ .last().unwrap());
		let i = header_log_T⃰ [..header_log_T⃰.len()-4].iter().position(|&header_log_T⃰ | header_log_T⃰  < log_T⃰).unwrap();
		let δ⃰ = reduced_dipole_moment(a, b);
		let polynomials : &[_; 3] = &table[i..i+3].try_into().unwrap();
		quadratic_interpolation(header_log_T⃰[i..][..3].try_into().unwrap(), polynomials.map(|P| eval_poly(&P, δ⃰ )), log_T⃰)
	};
	let Ω⃰22 = |a, b, T| collision_integral(&Ω⃰22, a, b, T);
	let viscosity = |a, T| 5./16. * sqrt(π * molar_mass[a] * kB * T / NA) / (Ω⃰22(a, a, T) * π * sq(diameter[a]));
	/*let reduced_mass = |a, b| molar_mass[a] * molar_mass[b] / (NA * (molar_mass[a] + molar_mass[b])); // reduced_mass (a,a) = W/2NA
	let _Ω⃰11 = |a, b, T| Ω⃰22(a, b, T)/collision_integral(&A⃰, a, b, T);
	let conductivity = |a, T| {
		let self_diffusion_coefficient = 3./16. * sqrt(2.*π/reduced_mass(a,a)) * pow(kB*T, 3./2.) / (π * sq(diameter[a]) * Ω⃰11(a, a, T));
		let f_internal = molar_mass[a]/(kB*NA * T) * self_diffusion_coefficient / viscosity(a, T);
		let T⃰ = T⃰ (a, a, T);
		let fz_T⃰ = 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰;
		// Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
		let fz_298 = (|T⃰| 1. + pow(π, 3./2.) / sqrt(T⃰) * (1./2. + 1./T⃰) + (1./4. * sq(π) + 2.) / T⃰)(kB * 298. / well_depth[a]);
		let c1 = 2./π * (5./2. - f_internal)/(rotational_relaxation[a] * fz_298 / fz_T⃰ + 2./π * (5./3. * internal_degrees_of_freedom[a] + f_internal));
		let f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[a]/(3./2.));
		let f_rotation = f_internal * (1. + c1);
		let Cv_internal = thermodynamics[a].specific_heat_capacity(T) - 5./2. - internal_degrees_of_freedom[a];
		(viscosity(a, T)/molar_mass[a])*kB*NA*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[a] + f_internal * Cv_internal)
	};
	let diffusion_coefficient = |a, b, T| {
		let reduced_diameter = (diameter[a] + diameter[b])/2. * pow(χ(a, b), -1./6.);
		3./16. * sqrt(2.*π/reduced_mass(a,b)) * pow(kB*T, 3./2.) / (π*sq(reduced_diameter)*Ω⃰11(a, b, T))
	};*/

	// polynomial fit in the temperature range for every specie (pairs)
	/*const*/let [temperature_min, temperature_max] : [f64; 2] = [300., 3000.];
	const D : usize = 4; // FIXME: Remez
	const N : usize = D+2; // FIXME: Remez
	let T = generate(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min));
	let sqrt_viscosity_sqrt_T_polynomials : [[_; D]; S] = generate(|a| polynomial_fit::<_,_,_,D,N>(T, log, |T| sqrt(viscosity(a,T)/sqrt(T)))).collect();
	//let conductivity_polynomials = generate(|a| polyfit(|T| conductivity(a, T)/sqrt(T)));
	//let diffusion_coefficient_polynomials = generate(|a| iter::box_collect((0..=a).map(|b| polyfit(|T| diffusion_coefficient(a,b,T)/pow(T, 3./2.)))));

	let T = temperature;
	let sqrt_viscosity = |a| sq(sqrt(T) * eval_poly(&sqrt_viscosity_sqrt_T_polynomials[a], log(T)));
	amounts.dot(generate(|k|
		sq(sqrt_viscosity(k)) /
		amounts.dot(generate(|j|
			sq(1. + (sqrt_viscosity(k)/sqrt_viscosity(j)) * sqrt(sqrt(molar_mass[j]/molar_mass[k]))) /
			(sqrt(8.) * sqrt(1. + molar_mass[k]/molar_mass[j])))
		)
	))
}}

#[derive(Clone)] pub struct State<const S: usize> where [(); S-1]: {
	pub temperature: f64,
	pub amounts: [f64; S/*-1*/]
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species: [&'t str; S],
	pub system: System<S>,
	pub time_step: f64,
	pub pressure_r: f64,
	pub state: State<S>
}

use std::{convert::TryInto, lazy::SyncLazy};
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl<const S: usize> Simulation<'t, S> where [(); S-1]: {
	pub fn new(system: &'b [u8]) -> ::ron::Result<Self> where 'b: 't, [(); S]: {
		let ron::System{species: species_data, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
		let ron::Phase::IdealGas{species, state, ..} = phases.into_vec().into_iter().next().unwrap();
		use {std::convert::TryInto, iter::into::IntoMap};
		let species : [_; S] = species.as_ref().try_into().unwrap_or_else(|_| panic!("Compiled for {} species, got {}", S, species.len()));
		let molar_mass = eval(species, |s| species_data.get(s).unwrap_or_else(|| panic!("{}", s)).composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = eval(species, |s| { let ron::Specie{thermodynamic: ron::NASA7{temperature_ranges, pieces},..} = &species_data[s]; match temperature_ranges[..] {
			[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
			[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
			ref ranges => panic!("{:?}", ranges),
		}});
		let reactions = iter::into::Collect::collect(reactions.map(|self::ron::Reaction{ref equation, rate_constant, model}| {
			for side in equation { for (specie, _) in side { assert!(species.contains(specie), "{}", specie) } }
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
		let diameter = eval(species, |s| species_data[s].transport.diameter);
		let well_depth = eval(species, |s| species_data[s].transport.well_depth);
		use self::ron::Geometry::*;
		let polarizability = eval(species, |s| if let Linear{polarizability,..}|Nonlinear{polarizability,..} = species_data[s].transport.geometry { polarizability } else { 0. });
		let _rotational_relaxation = eval(species, |s| if let Nonlinear{rotational_relaxation,..} = species_data[s].transport.geometry { rotational_relaxation } else { 0. });
		let dipole = eval(species, |s| if let Nonlinear{dipole,..} = species_data[s].transport.geometry { dipole } else { 0. });
		let _internal_degrees_of_freedom = eval(species, |s| match species_data[s].transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
		let ron::InitialState{temperature, pressure, mole_proportions} = state;
		let pressure_r = pressure / (kB*NA);
		let amount = pressure_r / temperature * System::<S>::volume;
		let mole_proportions = eval(&species, |specie| *mole_proportions.get(specie).unwrap_or(&0.));
		let amounts = eval(mole_proportions/*.prefix()*/, |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);
		//let mass = {use iter::into::Sum; map!(amounts, molar_masses; |n, W| n*W).sum()};

		Ok(Self{
			species,
			system: System{molar_mass, thermodynamics, reactions, diameter, well_depth, polarizability, dipole, /*rotational_relaxation, internal_degrees_of_freedom*/},
			time_step, pressure_r,
			state: State{temperature, amounts}
		})
	}
}
