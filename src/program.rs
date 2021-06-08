use std::default::default;
fn box_<T>(t: T) -> Box<T> { Box::new(t) }

fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(std::default::default()).push(index) }
	map
}

type Value = usize;
enum Expression {
	Literal(f64),
	Use(Value),
	Index { base: Value, index: usize },
	Neg(Box<Expression>),
	Add(Box<Expression>, Box<Expression>),
	Sub(Box<Expression>, Box<Expression>),
	Less(Box<Expression>, Box<Expression>),
	Mul(Box<Expression>, Box<Expression>),
	Div(Box<Expression>, Box<Expression>),
	Call { function: &'static str, arguments: Box<[Expression]> },
	Block { statements: Box<[Statement]>, result: Box<Expression> },
}

enum Statement {
	Definition { id: Value, value: Expression },
	Output { index: usize, value: Expression },
	ConditionalStatement { condition: Expression, consequent: Vec<Statement>, alternative: Vec<Statement> },
}

pub struct Subroutine {
	parameters: Vec<&'static str>,
	output: usize,
	statements: Vec<Statement>
}

fn r#use(index: Value) -> Expression { Expression::Use(index) }
fn index(base: Value, index: usize) -> Expression { Expression::Index { base, index } }

fn neg(x: impl Into<Expression>) -> Expression { Expression::Neg(box_(x.into())) }

fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Add(box_(a.into()), box_(b.into())) }
fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Sub(box_(a.into()), box_(b.into())) }
fn less(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Less(box_(a.into()), box_(b.into())) }
fn mul(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Mul(box_(a.into()), box_(b.into())) }
fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Div(box_(a.into()), box_(b.into())) }

fn exp2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "exp2", arguments: box_([x.into()]) } }
fn log2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "log2", arguments: box_([x.into()]) } }

fn define(id: Value, value: impl Into<Expression>) -> Statement { Statement::Definition{ id, value: value.into() } }
fn output(index: usize, value: impl Into<Expression>) -> Statement { Statement::Output{ index, value: value.into() } }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }

impl std::ops::Add<Expression> for Expression { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for Expression { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for Expression { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

struct Definition(Value);
impl From<&Definition> for Expression { fn from(d: &Definition) -> Expression { r#use(d.0) } }

impl<E:Into<Expression>> std::ops::Mul<E> for &Definition { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for &Definition { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

struct Block {
	base: usize,
	statements: Vec<Statement>,
}
impl FnOnce<(Expression,)> for Block {
	type Output = Definition;
	extern "rust-call" fn call_once(mut self, args: (Expression,)) -> Self::Output { self.call_mut(args) }
}
impl FnMut<(Expression,)> for Block {
	extern "rust-call" fn call_mut(&mut self, (value,): (Expression,)) -> Self::Output  {
		let id = self.base+self.statements.len();
		self.statements.push(define(id, value));
		Definition(id)
	}
}

struct Builder(usize);
impl Builder {
	fn new(parameters: &[&'static str]) -> Self { Self(parameters.len()) }
	fn block(&self, block: impl Fn(&mut Block)->Expression) -> Expression {
		let mut defines = Block{ base: self.0, statements: vec![] };
		let result = block(&mut defines);
		Expression::Block { statements: defines.statements.into(), result: box_(result) }
	}
}

impl From<f64> for Expression { fn from(v: f64) -> Self { Self::Literal(v) } }
impl std::ops::Add<f64> for Expression { type Output = Expression; fn add(self, b: f64) -> Self::Output { add(self, b) } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<f64> for Expression { type Output = Expression; fn sub(self, b: f64) -> Self::Output { sub(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Mul<&Definition> for f64 { type Output = Expression; fn mul(self, b: &Definition) -> Self::Output { mul(self, b) } }
impl std::ops::Div<f64> for Expression { type Output = Expression; fn div(self, b: f64) -> Self::Output { mul(1./b, self) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

impl FnOnce<(&[f64],&[&[f64]])> for Subroutine {
	type Output = Box<[f64]>; extern "rust-call" fn call_once(mut self, args: (&[f64],&[&[f64]])) -> Self::Output { self.call_mut(args) } }
impl FnMut<(&[f64],&[&[f64]])> for Subroutine { extern "rust-call" fn call_mut(&mut self, args: (&[f64],&[&[f64]])) -> Self::Output  { self.call(args) } }
impl Fn<(&[f64],&[&[f64]])> for Subroutine {
	extern "rust-call" fn call(&self, (arguments, arrays): (&[f64],&[&[f64]])) -> Self::Output  {
		struct State<'t> {
			arguments: &'t [f64],
			arrays: &'t [&'t [f64]],
			definitions: linear_map::LinearMap<Value, f64>,
			output: Box<[f64]>
		}
		impl State<'_> {
			fn eval(&self, expression: &Expression) -> f64 {
				use Expression::*;
				match expression {
					&Use(id) => {
						if id < self.arguments.len() { self.arguments[id] }
						else if id < self.arguments.len()+self.arrays.len() { unreachable!() }
						else { self.definitions[&id] }
					}
					Less(a, b) => if self.eval(a) < self.eval(b) { 1. } else { 0. },
					e => panic!("{:?}", e),
				}
			}
			fn run(&mut self, statements: &[Statement]) {
				for statement in statements {
					use Statement::*;
					match statement {
						ConditionalStatement { condition, consequent, alternative } => {
							if self.eval(condition) != 0. { self.run(consequent); } else { self.run(alternative); }
						}
						Output { index, value } => self.output[*index] = self.eval(value),
						s => panic!("{:?}", s),
					}
				}
			}
		}
		let mut state = State{arguments, arrays, definitions: default(), output: vec![0.; self.output].into_boxed_slice()};
		state.run(&self.statements);
		state.output
	}
}


use itertools::Itertools;

impl std::fmt::Debug for Expression {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		use Expression::*;
		match self {
			Literal(v) => v.fmt(f),
			Use(v) => write!(f, "#{}", v),
			Index { base, index } => write!(f, "{}[{}]", base, index),
			Neg(x) => write!(f, "-{:?}", x),
			Add(a, b) => write!(f, "({:?} + {:?})", a, b),
			Sub(a, b) => write!(f, "({:?} - {:?})", a, b),
			Less(a, b) => write!(f, "{:?} < {:?}", a, b),
			Mul(a, b) => write!(f, "({:?} * {:?})", a, b),
			Div(a, b) => write!(f, "({:?} / {:?})", a, b),
			Call { function, arguments } => write!(f, "{:?}({:?})", function, arguments.iter().format(", ")),
			Block { statements, result } => write!(f, "{{ {:?}; {:?} }}", statements.iter().format(";"), result),
		}
	}
}

impl std::fmt::Debug for Statement {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		use Statement::*;
		match self {
			Definition { id, value } => write!(f, "#{:?} = {:?}", id, value),
			Output { index, value } => write!(f, "@{:?} = {:?}", index, value),
			ConditionalStatement { condition, consequent, alternative } => write!(f, "if {:?} {{\n{:?}\n}} else {{\n{:?}\n}}", condition, consequent.iter().format("\n"), alternative.iter().format("\n")),
		}
	}
}

impl std::fmt::Debug for Subroutine {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		write!(f, "({:?}) -> {:?} {{\n{:?}\n}}", self.parameters.iter().format(", "), self.output, self.statements.iter().format("\n"))
	}
}

fn dot(c: &[f64], v: Value) -> Expression {
	c.iter().enumerate().fold(None, |sum, (i,&c)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + index(v, i), None => index(v, i)}) }
		else if c == -1. { Some(match sum { Some(sum) => sum - index(v, i), None => -index(v, i)}) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => c * index(v, i) + sum, None => c * index(v, i) }) }
	).unwrap()
}

fn product_of_exponentiations(c: &[impl Copy+Into<i16>], v: Value) -> Expression {
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

macro_rules! stringify{ [$($parameter:ident),*] => ([$(std::stringify!($parameter)),*]) }
fn parameters<const N: usize>(parameters: [&'static str; N]) -> ([&'static str; N], [Value; N]) { (parameters, iter::from_iter(iter::ConstRange::<N>)) }

use std::f64::consts::LN_2;
use super::*;

struct T<T> { log_T: T, T: T, T2: T, T3: T, T4: T, rcp_T: T, rcp_T2: T }
macro_rules! T{ {$($field:ident),*} => (T{$($field: r#use($field)),*}) }
impl From<&T<Value>> for T<Expression> { fn from(&T{log_T,T,T2,T3,T4,rcp_T,rcp_T2}:&T<Value>) -> Self { T!{log_T,T,T2,T3,T4,rcp_T,rcp_T2} } }

fn thermodynamic_function(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<Expression>)->Expression) -> Subroutine {
	let (parameters, [log_T,T,T2,T3,T4,rcp_T]) = parameters(stringify![log_T,T,T2,T3,T4,rcp_T]);
	let ref Ts = T{log_T,T,T2,T3,T4,rcp_T,rcp_T2:0};
	Subroutine{
		parameters: parameters.to_vec(),
		output: thermodynamics.len(),
		statements: bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())).into_iter().map(|(temperature_split, species)| Statement::ConditionalStatement{
			condition: less(r#use(T), f64::from_bits(temperature_split)),
			consequent: species.iter().map(|&specie| output(specie, expression(&thermodynamics[specie].pieces[0], Ts.into()))).collect(),
			alternative: species.iter().map(|&specie| output(specie, expression(&thermodynamics[specie].pieces[1], Ts.into()))).collect()
		}).collect()
	}
}

pub fn molar_heat_capacity_at_constant_pressure_R(thermodynamics: &[NASA7]) -> Subroutine {
		thermodynamic_function(&thermodynamics, |a, T{T,T2,T3,T4,..}| a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 )
}
pub fn enthalpy_RT(thermodynamics: &[NASA7]) -> Subroutine {
	thermodynamic_function(&thermodynamics, |a, T{T,T2,T3,T4,rcp_T,..}| a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T )
}
pub fn exp_Gibbs_RT(thermodynamics: &[NASA7]) -> Subroutine {
	thermodynamic_function(&thermodynamics, |a, T{log_T,T,T2,T3,T4,rcp_T,..}|
			exp2((a[0]-a[6])/LN_2-a[0]*log_T+a[1]/2./LN_2*T+(1./3.-1./2.)*a[2]/LN_2*T2+(1./4.-1./3.)*a[3]/LN_2*T3+(1./5.-1./4.)*a[4]/LN_2*T4+a[5]*rcp_T) )
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
														T{log_T,T,T2,T4,rcp_T,rcp_T2,..}: T<Expression>) -> Expression {
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

impl ReactionModel {
fn efficiency(&self, k_inf: &RateConstant, T: &T<Value>, concentrations: Value, f: &Builder) -> Expression {
	use ReactionModel::*; match self {
		Elementary|Irreversible => arrhenius(k_inf, T.into()),
		ThreeBody{efficiencies} => arrhenius(k_inf, T.into()) * dot(efficiencies, concentrations),
		PressureModification{efficiencies, k0} => f.block(|def|{
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T.into()));
			Pr / (rcp_arrhenius(k_inf, T.into()) * Pr + 1.)
		}),
		Falloff{efficiencies, k0, troe} => f.block(|def|{
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T.into()));
			let model::Troe{A, T3, T1, T2} = *troe;
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
}

pub fn rates(reactions: &[Reaction]) -> Subroutine {
	let species_len = reactions[0].reactants.len();
	let (parameters,          [log_T,T,T2,T4,rcp_T,rcp_T2, exp_Gibbs0_RT, P0_RT,rcp_P0_RT, concentrations])
	= parameters(stringify![log_T,T,T2,T4,rcp_T,rcp_T2, exp_Gibbs0_RT, P0_RT,rcp_P0_RT, concentrations]);
	let ref f = Builder::new(&parameters);
	let ref T = T{log_T,rcp_T,T,T2,T3:0,T4,rcp_T2};
	Subroutine {
		parameters: parameters.into(),
		output: species_len-1,
		statements: reactions.iter().enumerate().map(|(reaction, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
			let c = model.efficiency(rate_constant, T, concentrations, f); // todo: CSE
			let Rf = product_of_exponentiations(reactants, concentrations);
			let R = if let ReactionModel::Irreversible = model { Rf } else {
				let rcp_equilibrium_constant = product_of_exponentiations(net, exp_Gibbs0_RT)
																												- match Σnet { 0=>(0.).into(), 1=>r#use(P0_RT), -1=>r#use(rcp_P0_RT), _ => unreachable!()};
				let Rr = rcp_equilibrium_constant * product_of_exponentiations(products, concentrations);
				Rf - Rr
			};
			define(reaction, c * R)
		}).chain((0..species_len-1).map(|specie|
			output(specie, reactions.iter().enumerate().fold(None, |dtω, (reaction, Reaction{net, ..})| {
					let cR = r#use(reaction);
					match net[specie] {
						0 => dtω,
						1 => match dtω { None => Some(cR), Some(dtω) => Some(dtω + cR) }
						-1 => match dtω { None => Some(-cR), Some(dtω) => Some(dtω - cR) }
						ν => match dtω { None => Some((ν as f64) * cR), Some(dtω) => Some((ν as f64) * cR + dtω) }
					}
			}).unwrap())
		)).collect()
	}
}
