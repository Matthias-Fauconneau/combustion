#![feature(associated_type_bounds,bindings_after_at)]#![allow(uncommon_codepoints,confusable_idents,non_snake_case)]
fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use ast::*;

fn product_of_exponentiations(c: &[impl Copy+Into<i16>], v: &[Value]) -> Expression {
	let (num, div) : (Vec::<_>,Vec::<_>) = c.iter().map(|&c| c.into()).zip(v).filter(|&(c,_)| c!=0).partition(|&(c,_)| c>0);
	let num : Option<Expression> = num.into_iter().fold(None, |mut p, (c,e)|{ for _ in 0..c { p = Some(match p { Some(p) => p*e, None => e.into() }); } p });
	let div : Option<Expression> = div.into_iter().fold(None, |mut p, (c,e)|{ for _ in 0..-c { p = Some(match p { Some(p) => p*e, None => e.into() }); } p });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(1./div),
		(Some(num), Some(div)) => Some(num/div)
	}.unwrap()
}

use iter::map;
use chemistry::*;
#[derive(Clone, Copy)] struct T<'t> { log_T: &'t Value, T: &'t Value, T2: &'t Value, T3: &'t Value, T4: &'t Value, rcp_T: &'t Value, rcp_T2: &'t Value}
fn thermodynamics(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<'_>)->Expression, T: T<'_>, slots: Box<[Variable]>) -> (Box<[Statement]>, Box<[Variable]>) {
	let mut slots = map(slots.into_vec(), Some);
	let mut results = map(thermodynamics, |_| None);
	(map(bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())).into_iter(), |(temperature_split, ref species)| {
		let slots = map(species, |&specie| slots[specie].take().unwrap());
		let spline = Statement::Branch{
			condition: less(r#use(T.T), f64::from_bits(temperature_split)),
			consequent: map(species.iter().zip(&*slots), |(&specie, slot)| store(slot, expression(&thermodynamics[specie].pieces[0], T))),
			alternative: map(species.iter().zip(&*slots), |(&specie, slot)| store(slot, expression(&thermodynamics[specie].pieces[1], T)))
		};
		for (&specie, slot) in species.iter().zip(slots.into_vec().into_iter()) { assert!(results[specie].replace(slot).is_none()) }
		spline
	}), map(results.into_vec().into_iter(), Option::unwrap))
}
fn eval((algorithm, results): (Box<[Statement]>, Box<[Variable]>), f: &mut Block) -> Box<[Value]> {
	f.extend(algorithm.into_vec());
	map(results.into_vec(), |v| f.load(v))
}

//fn molar_heat_capacity_at_constant_pressure_R(a: &[f64], T{T,T2,T3,T4,..}: T) -> Expression { a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 }
fn enthalpy_RT(a: &[f64], T{T,T2,T3,T4,rcp_T,..}: T<'_>) -> Expression { a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T }
use std::f64::consts::LN_2;
fn exp_Gibbs_RT(a: &[f64], T{log_T,T,T2,T3,T4,rcp_T,..}: T<'_>) -> Expression {
	exp2((a[0]-a[6])/LN_2-a[0]*log_T-a[1]/2./LN_2*T+(1./3.-1./2.)*a[2]/LN_2*T2+(1./4.-1./3.)*a[3]/LN_2*T3+(1./5.-1./4.)*a[4]/LN_2*T4+a[5]/LN_2*rcp_T)
}

// A.T^β.exp(-Ea/kT)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{log_T,T,T2,T4,rcp_T,rcp_T2,..}: T<'_>) -> Expression {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		if temperature_exponent == 0. { A.into() }
		else if temperature_exponent == -1. { A * rcp_T }
		else if temperature_exponent == 1. { A * T }
		else if temperature_exponent == 2. { A * T2 }
		else if temperature_exponent == 4. { A * T4 }
		else if temperature_exponent == -2. { A * rcp_T2 }
		else { unreachable!() }
	} else {
		let βlogT𐊛logA = if temperature_exponent == 0. { f64::log2(A).into() } else { temperature_exponent * log_T + f64::log2(A) };
		let log_arrhenius = if activation_temperature == 0. { βlogT𐊛logA } else { -activation_temperature/LN_2 * rcp_T + βlogT𐊛logA };
		exp2(log_arrhenius)
	}
}

fn rcp_arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{log_T,T,T2,rcp_T,rcp_T2,..}: T<'_>) -> Expression {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		if temperature_exponent == 0. { (1./A).into() }
		else if temperature_exponent == -1. { T / A }
		else if temperature_exponent == 1. { rcp_T / A }
		else if temperature_exponent == 2. { rcp_T2 / A }
		else if temperature_exponent == -2. { T2 / A }
		else { unreachable!() }
	} else {
		let m_βlogT𐊛logA = if temperature_exponent == 0. { (-f64::log2(A)).into() } else { - temperature_exponent * log_T - f64::log2(A) };
		let m_log_arrhenius = if activation_temperature == 0. { m_βlogT𐊛logA } else { activation_temperature/LN_2 * rcp_T + m_βlogT𐊛logA };
		exp2(m_log_arrhenius)
	}
}

fn efficiency(model: &ReactionModel, k_inf: &RateConstant, T: T, concentrations: &[Value], f: &mut Block) -> Expression {
	use ReactionModel::*; match model {
		Elementary|Irreversible => arrhenius(k_inf, T),
		ThreeBody{efficiencies} => arrhenius(k_inf, T) * dot(efficiencies, concentrations),
		PressureModification{efficiencies, k0} => f.block(|def|{
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T));
			Pr / (rcp_arrhenius(k_inf, T) * Pr + 1.)
		}),
		Falloff{efficiencies, k0, troe} => {let Troe{A, T3, T1, T2} = *troe;/*ICE inside*/ f.block(|def|{
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T));
			let Fcent = {let T{T,rcp_T,..}=T; (1.-A) * exp2(r#use(T)/(-LN_2*T3)) + A * exp2(r#use(T)/(-LN_2*T1)) + exp2((-T2/LN_2)*r#use(rcp_T))};
			let ref logFcent = def(log2(Fcent));
			let c =-0.67*logFcent - 0.4*f64::log2(10.);
			let N = -1.27*logFcent + 0.75*f64::log2(10.);
			let ref logPr𐊛c = def(log2(Pr) + c);
			let ref f1 = def(logPr𐊛c / (-0.14*logPr𐊛c+N));
			let F = exp2(logFcent/(f1*f1+1.));
			Pr / (rcp_arrhenius(k_inf, T) * Pr + 1.) * F
		})}
	}
}

fn species_rates(reactions: &[Reaction], T: T, C0: &Value, rcp_C0: &Value, exp_Gibbs0_RT: &[Value], concentrations: &[Value], f: &mut Block) -> Box<[Expression]> {
	let rates = iter::map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let efficiency = efficiency(model, rate_constant, T, concentrations, f); // todo: CSE
		let forward_rate = product_of_exponentiations(reactants, concentrations);
		let rate = if let ReactionModel::Irreversible = model { forward_rate } else {
			let rcp_equilibrium_constant_0 = product_of_exponentiations(net, exp_Gibbs0_RT);
			let rcp_equilibrium_constant = match -Σnet {
				0 => rcp_equilibrium_constant_0,
				1 => C0 * rcp_equilibrium_constant_0,
				-1 => rcp_C0 * rcp_equilibrium_constant_0,
				_ => unreachable!()
			};
			let reverse_rate = rcp_equilibrium_constant * product_of_exponentiations(products, concentrations);
			forward_rate - reverse_rate
		};
		f.def(efficiency * rate)
	});
	map(0..exp_Gibbs0_RT.len(), |specie| idot(reactions.iter().map(|Reaction{net, ..}| net[specie] as f64).zip(&*rates)))
}

pub fn rates(active_species: &[NASA7], reactions: &[Reaction]) -> Function {
	let_!{ input@[ref pressure_R, ref total_amount, ref T, ref active_amounts @ ..] = &*map(0..(3+active_species.len()), Value) => {
	let mut function = FunctionBuilder::new(input);
	let mut f = Block::new(&mut function);
	let ref log_T = f.def(log2(T));
	let ref T2 = f.def(T*T);
	let ref T3 = f.def(T*T2);
	let ref T4 = f.def(T*T3);
	let ref rcp_T = f.def(1./T);
	let ref rcp_T2 = f.def(num::sq(rcp_T));
	let ref rcp_P0_RT = f.def((1./NASA7::reference_pressure) * T);
	let ref P0_RT = f.def(NASA7::reference_pressure * rcp_T);
	let ref total_concentration = f.def(pressure_R / T); // Constant pressure
	let T = T{log_T,T,T2,T3,T4,rcp_T,rcp_T2};
	let ref exp_Gibbs0_RT = eval(thermodynamics(active_species, exp_Gibbs_RT, T, map(active_species, |_| f.decl())), &mut f);
	let ref density = f.def(total_concentration / total_amount);
	let active_concentrations = map(0..active_species.len(), |k| f.def(density*max(0., &active_amounts[k])));
	let inert_concentration = f.def(total_concentration - active_concentrations.iter().sum::<Expression>());
	let ref concentrations = iter::box_(active_concentrations.into_vec().into_iter().chain(std::iter::once(inert_concentration)));
	let rates = map(species_rates(reactions, T, P0_RT, rcp_P0_RT, exp_Gibbs0_RT, concentrations, &mut f).into_vec(), |e| f.def(e));
	let enthalpy_RT = eval(thermodynamics(active_species, enthalpy_RT, T, map(active_species, |_| f.decl())), &mut f);
	let energy_rate_RT : Expression = iter::dot(rates.iter().zip(&*enthalpy_RT));
	f.push(output(0, energy_rate_RT));
	f.extend(rates.iter().enumerate().map(|(i, r)| output(1+i, r)));
	Function::new(1+rates.len(), f.into(), function)
	}}
}
