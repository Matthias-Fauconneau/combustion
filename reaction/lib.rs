#![feature(associated_type_bounds,bindings_after_at,iter_is_partitioned,array_map,format_args_capture,trait_alias)]#![allow(uncommon_codepoints,confusable_idents,non_snake_case)]
fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use ast::{*, Expression::Float as c};

fn product_of_exponentiations(b: &[Value], n: &[impl Copy+Into<i16>]) -> Expression {
	let (num, div) : (Vec::<_>,Vec::<_>) = n.iter().map(|&n| n.into()).zip(b).filter(|&(n,_)| n!=0).partition(|&(n,_)| n>0);
	let num : Option<Expression> = num.into_iter().fold(None, |mut p, (n,e)|{ for _ in 0..n { p = Some(match p { Some(p) => p*e, None => e.into() }); } p });
	let div : Option<Expression> = div.into_iter().fold(None, |mut p, (n,e)|{ for _ in 0..-n { p = Some(match p { Some(p) => p*e, None => e.into() }); } p });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(1./div),
		(Some(num), Some(div)) => Some(num/div)
	}.unwrap()
}

use iter::{box_, map};
use chemistry::*;
#[derive(Clone, Copy)] struct T<'t> { ln_T: &'t Value, T: &'t Value, T2: &'t Value, T3: &'t Value, T4: &'t Value, rcp_T: &'t Value, rcp_T2: &'t Value}
fn thermodynamics(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<'_>, &mut Block)->Expression, T: T<'_>, f: &mut Block) -> Box<[Value]> {
	let mut specie_results = map(thermodynamics, |_| None);
	for (temperature_split, ref species) in bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())) {
		let results = map(species, |_| f.value());
		for (&specie, result) in species.iter().zip(&*results) { assert!(specie_results[specie].replace(result.clone()).is_none()) }
		push(Statement::Select{
			condition: less_or_equal(T.T, c(f64::from_bits(temperature_split))),
			true_exprs: map(species, |&specie| expression(&thermodynamics[specie].pieces[0], T, f)),
			false_exprs: map(species, |&specie| expression(&thermodynamics[specie].pieces[1], T, f)),
			results
		}, f);
	}
	map(specie_results.into_vec().into_iter(), Option::unwrap)
}

fn molar_heat_capacity_at_constant_pressure_R(a: &[f64], T{T,T2,T3,T4,..}: T, _: &mut Block) -> Expression { a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 }
fn enthalpy_RT(a: &[f64], T{T,T2,T3,T4,rcp_T,..}: T<'_>, _: &mut Block) -> Expression { a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T }
fn Gibbs_RT(a: &[f64], T{ln_T,T,T2,T3,T4,rcp_T,..}: T<'_>, _: &mut Block) -> Expression {
	(a[0]-a[6])-a[0]*ln_T-a[1]/2.*T+(1./3.-1./2.)*a[2]*T2+(1./4.-1./3.)*a[3]*T3+(1./5.-1./4.)*a[4]*T4+a[5]*rcp_T
}
fn exp_Gibbs_RT(a: &[f64], T: T<'_>, f: &mut Block) -> Expression { exp(Gibbs_RT(a, T, f), f) }

// A.T^β.exp(-Ea/kT)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{ln_T,T,T2,T4,rcp_T,rcp_T2,..}: T<'_>, f: &mut Block) -> Expression {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		if temperature_exponent == 0. { c(A) }
		else if temperature_exponent == -1. { A * rcp_T }
		else if temperature_exponent == 1. { A * T }
		else if temperature_exponent == 2. { A * T2 }
		else if temperature_exponent == 4. { A * T4 }
		else if temperature_exponent == -2. { A * rcp_T2 }
		else { unreachable!() }
	} else {
		let βlnT𐊛lnA = if temperature_exponent == 0. { c(f64::ln(A)) } else { temperature_exponent * ln_T + c(f64::ln(A)) };
		let ln_arrhenius = if activation_temperature == 0. { βlnT𐊛lnA } else { -activation_temperature * rcp_T + βlnT𐊛lnA };
		exp(ln_arrhenius, f)
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
			let Fcent = {let T{T,rcp_T,..}=T; sum([
				(T3 > 1e-30).then(|| { let y = 1.-A; if T3<1e30 { y * exp(T/(-T3), f) } else { y.into() }}),
				(T1 > 1e-30).then(|| { let y = A; if T1<1e30 { y * exp(T/(-T1), f) } else { y.into() }}),
				(T2 != 0.).then(|| exp((-T2)*rcp_T, f))
			].into_iter().filter_map(|x| x))};
			let ref lnFcent = def(ln(1./2., Fcent, f), f); // 0.1-0.7 => e-3
			let C =-0.67*lnFcent - 0.4*f64::ln(10.);
			let N = -1.27*lnFcent + 0.75*f64::ln(10.);
			let ref lnPr𐊛C = def(ln(1., Pr, f) + C, f); // 2m - 2K
			let ref f1 = f(lnPr𐊛C / (-0.14*lnPr𐊛C+N));
			let F = exp(lnFcent/(f1*f1+1.), f);
			k_inf * Pr / (Pr + 1.) * F
		})}
	}
}

fn reaction_rates(reactions: &[Reaction], T: T, C0: &Value, rcp_C0: &Value, exp_Gibbs0_RT: &[Value], concentrations: &[Value], f: &mut Block) -> Box<[Value]> {
	map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let forward_rate_constant = forward_rate_constant(model, rate_constant, T, concentrations, f); // todo: CSE
		let forward = product_of_exponentiations(concentrations, reactants);
		let coefficient = if let ReactionModel::Irreversible = model { forward } else {
			let rcp_equilibrium_constant_0 = product_of_exponentiations(exp_Gibbs0_RT, net);
			//let rcp_equilibrium_constant_0 = exp2(idot(net.iter().map(|&net| net as f64).zip(Gibbs0_RT)), f);
			let rcp_equilibrium_constant = match -Σnet { // reverse_rate_constant / forward_rate_constant
				0 => rcp_equilibrium_constant_0,
				1 => C0 * rcp_equilibrium_constant_0,
				-1 => rcp_C0 * rcp_equilibrium_constant_0,
				_ => unreachable!()
			};
			let reverse = rcp_equilibrium_constant * product_of_exponentiations(concentrations, products);
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
	let ref ln_T = def(ln(1024., T, &mut f), &mut f);
	let ref T2 = f(T*T);
	let ref T3 = f(T*T2);
	let ref T4 = f(T*T3);
	let ref rcp_T = f(1./T);
	let ref rcp_T2 = f(num::sq(rcp_T));
	let ref rcp_C0 = f((1./NASA7::reference_pressure) * T);
	let ref C0 = f(NASA7::reference_pressure * rcp_T);
	let ref total_concentration = f(pressure_R / T); // Constant pressure
	let T = T{ln_T,T,T2,T3,T4,rcp_T,rcp_T2};
	let ref exp_Gibbs0_RT = thermodynamics(&species[0..active], exp_Gibbs_RT, T, &mut f);
	let ref density = f(total_concentration / total_amount);
	let active_concentrations = map(0..active, |k| f(density*max(0., &active_amounts[k])));
	let inert_concentration = f(total_concentration - active_concentrations.iter().sum::<Expression>());
	let concentrations = box_(active_concentrations.into_vec().into_iter().chain([inert_concentration].into_iter()));
	let rates = reaction_rates(reactions, T, C0, rcp_C0, exp_Gibbs0_RT, &concentrations, &mut f);
	let rates = map(0..active, |specie| f(idot(reactions.iter().map(|Reaction{net, ..}| net[specie] as f64).zip(&*rates))));
	let enthalpy_RT = thermodynamics(&species[0..active], enthalpy_RT, T, &mut f);
	let dot = |a:&[Value], b:&[Value]| iter::dot(a.iter().zip(b.iter()));
	let energy_rate_RT : Expression = dot(&rates, &enthalpy_RT);
	let Cp : Expression = dot(&concentrations, &thermodynamics(species, molar_heat_capacity_at_constant_pressure_R, T, &mut f));
	let dtT_T = - energy_rate_RT / Cp;
	Function{output: box_([T.T * dtT_T].into_iter().chain(rates.iter().map(|v| v.into()))), statements: f.into(), input: input.len()}
}}}
