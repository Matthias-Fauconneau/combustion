#![feature(associated_type_bounds,bindings_after_at,iter_is_partitioned)]#![allow(uncommon_codepoints,confusable_idents,non_snake_case)]
fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use ast::{*, Type::*};

fn exp2(x: impl Into<Expression>, f: &mut Block) -> Expression {
	let ref x = f(max(fdemote(x), -126.99999f32));
	let ref ipart = f(fcvt_to_sint(x-1./2f32));
	let ref x = f(x - fcvt_from_sint(ipart));
	fpromote(cast(F32, ishl_imm(iadd(ipart, u32(127)), 23)) * fma(fma(fma(fma(fma(0.0018775767f32, x, 0.0089893397f32), x, 0.055826318f32), x, 0.24015361f32), x, 0.69315308f32), x, 0.99999994f32))
}
fn log2(x: impl Into<Expression>, f: &mut Block) -> Expression {
	let ref i = f(cast(I32, fdemote(x)));
	let ref m = f(cast(F32, or(and(i, u32(0x007FFFFF)), u32(1f32.to_bits()))));
	fpromote(fma(fma(fma(fma(fma(-0.034436006f32, m, 0.31821337f32), m, -1.2315303f32), m, 2.5988452f32), m, -3.3241990f32), m, 3.1157899f32) * (m - 1f32) + fcvt_from_sint(isub(ushr_imm(and(i, u32(0x7F800000)), 23), u32(127))))
}

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

use iter::{box_, map};
use chemistry::*;
#[derive(Clone, Copy)] struct T<'t> { log_T: &'t Value, T: &'t Value, T2: &'t Value, T3: &'t Value, T4: &'t Value, rcp_T: &'t Value, rcp_T2: &'t Value}
fn thermodynamics(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<'_>, &mut Block)->Expression, T: T<'_>, f: &mut Block) -> Box<[Value]> {
	let mut specie_results = map(thermodynamics, |_| None);
	for (temperature_split, ref species) in bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())) {
		let results = map(0..species.len(), |_| f.value());
		for (&specie, result) in species.iter().zip(&*results) { assert!(specie_results[specie].replace(Value(result.0)).is_none()) }
		push(Statement::Branch{
			condition: less_or_equal(T.T, f64::from_bits(temperature_split)),
			consequent: map(species, |&specie| expression(&thermodynamics[specie].pieces[0], T, f)),
			alternative: map(species, |&specie| expression(&thermodynamics[specie].pieces[1], T, f)),
			results
		}, f);
	}
	map(specie_results.into_vec().into_iter(), Option::unwrap)
}

fn molar_heat_capacity_at_constant_pressure_R(a: &[f64], T{T,T2,T3,T4,..}: T, _: &mut Block) -> Expression { a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 }
fn enthalpy_RT(a: &[f64], T{T,T2,T3,T4,rcp_T,..}: T<'_>, _: &mut Block) -> Expression { a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T }
use std::f64::consts::LN_2;
fn Gibbs_RT(a: &[f64], T{log_T,T,T2,T3,T4,rcp_T,..}: T<'_>, _: &mut Block) -> Expression {
	(a[0]-a[6])/LN_2-a[0]*log_T-a[1]/2./LN_2*T+(1./3.-1./2.)*a[2]/LN_2*T2+(1./4.-1./3.)*a[3]/LN_2*T3+(1./5.-1./4.)*a[4]/LN_2*T4+a[5]/LN_2*rcp_T
}
//fn exp_Gibbs_RT(a: &[f64], T: T<'_>, f: &mut Block) -> Expression { exp2(Gibbs_RT(a, T, f), f) }

// A.T^β.exp(-Ea/kT)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{log_T,T,T2,T4,rcp_T,rcp_T2,..}: T<'_>, f: &mut Block) -> Expression {
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
		exp2(log_arrhenius, f)
	}
}

fn forward_rate_constant(model: &ReactionModel, k_inf: &RateConstant, T: T, concentrations: &[Value], f: &mut Block) -> Expression {
	use ReactionModel::*; match model {
		Elementary|Irreversible => arrhenius(k_inf, T, f),
		ThreeBody{efficiencies} => arrhenius(k_inf, T, f) * dot(efficiencies, concentrations),
		PressureModification{efficiencies, k0} => f.block(|f|{
			let ref C_k0 = def(dot(efficiencies, concentrations) * arrhenius(k0, T, f), f);
			let ref k_inf = def(arrhenius(k_inf, T, f), f);
			(C_k0 * k_inf) / (C_k0 + k_inf)
		}),
		Falloff{efficiencies, k0, troe} => {let Troe{A, T3, T1, T2} = *troe;/*ICE inside*/ f.block(|f|{
			let ref k_inf = def(arrhenius(k_inf, T, f), f);
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T, f) / k_inf, f);
			let Fcent = {let T{T,rcp_T,..}=T; (1.-A) * exp2(T/(-LN_2*T3), f) + A * exp2(T/(-LN_2*T1), f) + exp2((-T2/LN_2)*rcp_T, f)};
			let ref logFcent = def(log2(Fcent, f), f);
			let c =-0.67*logFcent - 0.4*f64::log2(10.);
			let N = -1.27*logFcent + 0.75*f64::log2(10.);
			let ref logPr𐊛c = def(log2(Pr, f) + c, f);
			let ref f1 = f(logPr𐊛c / (-0.14*logPr𐊛c+N));
			let F = exp2(logFcent/(f1*f1+1.), f);
			k_inf * Pr / (Pr + 1.) * F
		})}
	}
}

fn reaction_rates(reactions: &[Reaction], T: T, C0: &Value, rcp_C0: &Value, /*exp_*/Gibbs0_RT: &[Value], concentrations: &[Value], f: &mut Block) -> Box<[Value]> {
	map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let forward_rate_constant = forward_rate_constant(model, rate_constant, T, concentrations, f); // todo: CSE
		let forward = product_of_exponentiations(reactants, concentrations);
		let coefficient = if let ReactionModel::Irreversible = model { forward } else {
			let rcp_equilibrium_constant_0 = exp2(idot(net.iter().map(|&net| net as f64).zip(Gibbs0_RT)), f); //product_of_exponentiations(net, exp_Gibbs0_RT);
			let rcp_equilibrium_constant = match -Σnet {
				0 => rcp_equilibrium_constant_0,
				1 => C0 * rcp_equilibrium_constant_0,
				-1 => rcp_C0 * rcp_equilibrium_constant_0,
				_ => unreachable!()
			};
			let reverse = rcp_equilibrium_constant * product_of_exponentiations(products, concentrations);
			forward - reverse
		};
		f(forward_rate_constant * coefficient)
	})
}

pub fn rates(species: &[NASA7], reactions: &[Reaction]) -> Function {
	let active = {
		let active = map(0..species.len()-1, |k| reactions.iter().any(|Reaction{net,..}| net[k] != 0));
		assert!(active.iter().is_partitioned(|&active| active));
		active.iter().position(|active| !active).unwrap_or(species.len()-1)
	};
	let_!{ input@[ref pressure_R, ref total_amount, ref T, ref active_amounts @ ..] = &*map(0..(3+species.len()-1), Value) => {
	let mut function = FunctionBuilder::new(input);
	let mut f = Block::new(&mut function);
	let ref log_T = def(log2(T, &mut f), &mut f);
	let ref T2 = f(T*T);
	let ref T3 = f(T*T2);
	let ref T4 = f(T*T3);
	let ref rcp_T = f(1./T);
	let ref rcp_T2 = f(num::sq(rcp_T));
	let ref rcp_C0 = f((1./NASA7::reference_pressure) * T);
	let ref C0 = f(NASA7::reference_pressure * rcp_T);
	let ref total_concentration = f(pressure_R / T); // Constant pressure
	let T = T{log_T,T,T2,T3,T4,rcp_T,rcp_T2};
	let ref Gibbs0_RT = thermodynamics(&species[0..active], Gibbs_RT, T, &mut f);
	//let ref exp_Gibbs0_RT = thermodynamics(&species[0..active], exp_Gibbs_RT, T, &mut f);
	let ref density = f(total_concentration / total_amount);
	let active_concentrations = map(0..active, |k| f(density*max(0., &active_amounts[k])));
	let inert_concentration = f(total_concentration - active_concentrations.iter().sum::<Expression>());
	let concentrations = box_(active_concentrations.into_vec().into_iter().chain([inert_concentration].into_iter()));
	let rates = reaction_rates(reactions, T, C0, rcp_C0, /*exp_*/Gibbs0_RT, &concentrations, &mut f);
	let rates = map(0..active, |specie| f(idot(reactions.iter().map(|Reaction{net, ..}| net[specie] as f64).zip(&*rates))));
	let enthalpy_RT = thermodynamics(&species[0..active], enthalpy_RT, T, &mut f);
	let dot = |a:&[Value], b:&[Value]| iter::dot(a.iter().zip(b.iter()));
	let energy_rate_RT : Expression = dot(&rates, &enthalpy_RT);
	let Cp : Expression = dot(&concentrations, &thermodynamics(species, molar_heat_capacity_at_constant_pressure_R, T, &mut f));
	let dtT_T = - energy_rate_RT / Cp;
	Function{output: box_([T.T * dtT_T].into_iter().chain(rates.iter().map(|v| v.into()))), statements: f.into(), input: input.len()}
}}}
