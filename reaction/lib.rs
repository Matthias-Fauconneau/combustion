#![feature(associated_type_bounds)]#![allow(uncommon_codepoints,confusable_idents,non_snake_case)]
fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use ast::*;
fn cdot(c: &[f64], v: &Value) -> Expression {
	c.iter().enumerate().fold(None, |sum, (i,&c)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + index(v, i), None => index(v, i)}) }
		else if c == -1. { Some(match sum { Some(sum) => sum - index(v, i), None => -index(v, i)}) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => c * index(v, i) + sum, None => c * index(v, i) }) }
	).unwrap()
}

fn product_of_exponentiations(c: &[impl Copy+Into<i16>], v: &Value) -> Expression {
	let (num, div) : (Vec::<_>,Vec::<_>) = c.iter().map(|&c| c.into()).enumerate().filter(|&(_,c)| c!=0).partition(|&(_,c)| c>0);
	let num = num.into_iter().fold(None, |mut a, (i,c)|{ for _ in 0..c { a = Some(match a { Some(a) => a*index(v,i), None => index(v,i) }); } a });
	let div = div.into_iter().fold(None, |mut a, (i,c)|{ for _ in 0..-c { a = Some(match a { Some(a) => a*index(v,i), None => index(v,i) }); } a });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(1./div),
		(Some(num), Some(div)) => Some(num/div)
	}.unwrap()
}

use chemistry::*;
struct T<'t> { log_T: &'t Value, T: &'t Value, T2: &'t Value, T3: &'t Value, T4: &'t Value, rcp_T: &'t Value, rcp_T2: &'t Value}
fn thermodynamics(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<'t>)->Expression, T: T<'t>, slots: Box<[Variable]>) -> (Box<[Statement]>, Box<[Variable]>) {
	let mut slots = map(slots, Some);
	let mut results = vec![None; slots.len()];
	(bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())).into_iter().map(|(temperature_split, species)| {
		let species = map(species, |&specie| (slots[specie].take(), &thermodynamics[specie]));
		let spline = Statement::ConditionalStatement{
			condition: less(r#use(T), f64::from_bits(temperature_split)),
			consequent: map(&*species, |(slot, spline)| assign(slot, expression(&spline.pieces[0], T))),
			alternative: map(&*species, |(slot, spline)| assign(slot, expression(&spline.pieces[1], T)))
		})
		for (slot, _) in species { assert!(results[specie].replace(slot).is_none()) }
		spline
	}, map(results, Option::unwrap))
}
fn eval(f: &mut Block, algorithm: Box<[Statement]>, results: Box<[Variable]>) -> Box<[Value]> {
	f.extend(algorithm);
	map(results, |v| f.load(v))
}

fn molar_heat_capacity_at_constant_pressure_R(a: &[f64], T{T,T2,T3,T4,..}: T) -> Expression { a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 }
fn enthalpy_RT(a: &[f64], T{T,T2,T3,T4,rcp_T,..}: T<'_>) -> Expression { a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T }
use std::f64::consts::LN_2;
fn exp_Gibbs_RT(a: &[f64], T{log_T,T,T2,T3,T4,rcp_T,..}: T<'_>) -> Expression {
	exp2((a[0]-a[6])/LN_2-a[0]*log_T-a[1]/2./LN_2*T+(1./3.-1./2.)*a[2]/LN_2*T2+(1./4.-1./3.)*a[3]/LN_2*T3+(1./5.-1./4.)*a[4]/LN_2*T4+a[5]/LN_2*rcp_T)
}

// A.T^β.exp(-Ea/kT)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant,
											T{log_T,T,T2,T4,rcp_T,rcp_T2,..}: T<Expression>) -> Expression {
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

fn rcp_arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant,
														T{log_T,T,T2,rcp_T,rcp_T2,..}: T<Expression>) -> Expression {
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

fn efficiency(model: &ReactionModel, k_inf: &RateConstant, T: &T, concentrations: &Value, f: &Block) -> Expression {
	use ReactionModel::*; match model {
		Elementary|Irreversible => arrhenius(k_inf, T.into()),
		ThreeBody{efficiencies} => arrhenius(k_inf, T.into()) * cdot(efficiencies, concentrations),
		PressureModification{efficiencies, k0} => f.block(|def|{
			let ref Pr = def(cdot(efficiencies, concentrations) * arrhenius(k0, T.into()));
			Pr / (rcp_arrhenius(k_inf, T.into()) * Pr + 1.)
		}),
		Falloff{efficiencies, k0, troe} => f.block(|def|{
			let ref Pr = def(cdot(efficiencies, concentrations) * arrhenius(k0, T.into()));
			let Troe{A, T3, T1, T2} = *troe;
			let Fcent = {let &T{T,rcp_T,..}=T; (1.-A) * exp2(r#use(T)/(-LN_2*T3)) + A * exp2(r#use(T)/(-LN_2*T1)) + exp2((-T2/LN_2)*r#use(rcp_T))};
			let ref logFcent = def(log2(Fcent));
			let c =-0.67*logFcent - 0.4*f64::log2(10.);
			let N = -1.27*logFcent + 0.75*f64::log2(10.);
			let ref logPr𐊛c = def(log2(Pr) + c);
			let ref f1 = def(logPr𐊛c / (-0.14*logPr𐊛c+N));
			let F = exp2(logFcent/(f1*f1+1.));
			Pr / (rcp_arrhenius(k_inf, T.into()) * Pr + 1.) * F
		})
	}
}

pub fn species_rates(f: &mut Block, reactions: &[Reaction], T: &T, P0_RT: &Value, rcp_P0_RT: &Value, exp_Gibbs0_RT: &[Value], concentrations: &[Value]) -> Box<Expression> {
	let species_len = reactions[0].reactants.len();
	let rates = iter::map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let efficiency = efficiency(model, rate_constant, T, concentrations, &f); // todo: CSE
		let forward_rate = product_of_exponentiations(reactants, concentrations);
		let rate = if let ReactionModel::Irreversible = model { forward_rate } else {
			let rcp_equilibrium_constant = product_of_exponentiations(net, exp_Gibbs0_RT)
																											* match -Σnet { 0=>(1.).into(), 1=>r#use(P0_RT), -1=>r#use(rcp_P0_RT), _ => unreachable!()};
			let reverse_rate = rcp_equilibrium_constant * product_of_exponentiations(products, concentrations);
			forward_rate - reverse_rate
		};
		f.def(efficiency * rate)
	});
	map(0..concentrations.len(), |specie| rates.iter().zip(reactions.iter()).fold(None, |dtω, (rate, Reaction{net, ..})| {
		let rate = r#use(rate);
		match net[specie] {
			0 => dtω,
			1 => match dtω { None => Some(rate), Some(dtω) => Some(dtω + rate) }
			-1 => match dtω { None => Some(-rate), Some(dtω) => Some(dtω - rate) }
			ν => match dtω { None => Some((ν as f64) * rate), Some(dtω) => Some((ν as f64) * rate + dtω) }
		}
	}).unwrap())
}

pub fn rates(species: &[NASA7], reactions: &[Reaction]) -> Subroutine<2, 1> {
	let (parameters, [ref T, ref pressure_R], [ref active_amounts])
		 = parameters!([ref T, ref pressure_R], [ref active_amounts]);
	let ref mut f = Block::new(&parameters);
	let ref log_T = f.def(f64::log2(*T));
	let ref T2 = f.def(T*T);
	let ref T3 = f.def(T*T2);
	let ref T4 = f.def(T*T3);
	let ref rcp_T = f.def(1./T);
	let ref rcp_T2 = f.def(rcp_T*rcp_T);
	let rcp_P0_RT = (1./NASA7::reference_pressure) * T;
	let P0_RT = NASA7::reference_pressure * rcp_T;
	let total_concentration = f.def(pressure_R / T); // Constant pressure
	let ref T = T{log_T,T,T2,T3,T4,rcp_T,rcp_T2};
	let exp_Gibbs0_RT = eval(f, thermodynamics(species, exp_Gibbs_RT, T, map(species, |_| Some(f.decl()))));
	let density = f.def(total_concentration / total_amount);
	let active_concentrations = map(active_amounts, |&n| f.def(density*max(0., n)));
	let inert_concentration = total_concentration - active_concentrations.iter().sum::<Expression>();
	let concentrations = [&active_concentrations as &[_],&[inert_concentration]].concat();
	let rates = map(species_rates(f, reactions, T, P0_RT, rcp_P0_RT, exp_Gibbs0_RT, concentrations), |e| f.def(e));
	let enthalpy_RT = eval(f, thermodynamics(species, enthalpy_RT, T, map(species, |_| Some(f.decl()))));
	let energy_rate_RT = dot(rates.iter(&*enthalpy_RT));
	f.push(output(0, energy_rate_RT));
	f.extend(map(rates.enumerate(), output));
	//f.extend(map(rates.enumerate(), |(slot, rate)| output(slot, rate)));
	Function{parameters: parameters.into(), output: 1+species_len-1, statements: f.into()}
}
