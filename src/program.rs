use std::default::default;
fn box_<T>(t: T) -> Box<T> { Box::new(t) }

#[derive(PartialEq, Eq, Clone, Copy)] pub struct Value(usize);
impl Value {
	pub const NONE: Value = Value(!0);
}

pub enum Expression {
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

pub enum Statement {
	Definition { id: Value, value: Expression },
	Output { index: usize, value: Expression },
	ConditionalStatement { condition: Expression, consequent: Vec<Statement>, alternative: Vec<Statement> },
}

pub struct Subroutine {
	pub parameters: Vec<&'static str>,
	pub output: usize,
	pub statements: Vec<Statement>
}

pub fn r#use(index: Value) -> Expression { Expression::Use(index) }
fn index(base: Value, index: usize) -> Expression { Expression::Index { base, index } }

fn neg(x: impl Into<Expression>) -> Expression { Expression::Neg(box_(x.into())) }

fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Add(box_(a.into()), box_(b.into())) }
fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Sub(box_(a.into()), box_(b.into())) }
pub fn less(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Less(box_(a.into()), box_(b.into())) }
fn mul(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Mul(box_(a.into()), box_(b.into())) }
fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Div(box_(a.into()), box_(b.into())) }

#[must_use] fn define(id: Value, value: impl Into<Expression>) -> Statement { Statement::Definition{ id, value: value.into() } }
#[must_use] pub fn output(index: usize, value: impl Into<Expression>) -> Statement { Statement::Output{ index, value: value.into() } }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }

impl std::ops::Add<Expression> for Expression { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for Expression { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for Expression { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

impl From<&Value> for Expression { fn from(value: &Value) -> Expression { r#use(*value) } }
impl<E:Into<Expression>> std::ops::Mul<E> for &Value { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for &Value { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

pub struct Block {
	base: usize,
	pub statements: Vec<Statement>,
}
impl Block {
	pub fn new(parameters: &[&'static str]) -> Self { Self { base: parameters.len(), statements: vec![] } }
	pub fn block(&self, build: impl Fn(&mut Block)->Expression) -> Expression {
		let mut block = Block{ base: self.base+self.statements.len(), statements: vec![] };
		let result = build(&mut block);
		Expression::Block { statements: block.statements.into(), result: box_(result) }
	}
	pub fn def(&mut self, value: Expression) -> Value {
		let id = Value(self.base+self.statements.len());
		self.statements.push(define(id, value));
		id
	}
	pub fn extend(&mut self, iter: impl IntoIterator<Item=Statement>) { self.statements.extend(iter) }
}
impl FnOnce<(Expression,)> for Block { type Output = Value; extern "rust-call" fn call_once(mut self, args: (Expression,)) -> Self::Output { self.call_mut(args) }}
impl FnMut<(Expression,)> for Block { extern "rust-call" fn call_mut(&mut self, (value,): (Expression,)) -> Self::Output  { self.def(value) } }

impl From<f64> for Expression { fn from(v: f64) -> Self { Self::Literal(v) } }
impl std::ops::Add<f64> for Expression { type Output = Expression; fn add(self, b: f64) -> Self::Output { add(self, b) } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<f64> for Expression { type Output = Expression; fn sub(self, b: f64) -> Self::Output { sub(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Mul<&Value> for f64 { type Output = Expression; fn mul(self, b: &Value) -> Self::Output { mul(self, b) } }
impl std::ops::Div<f64> for Expression { type Output = Expression; fn div(self, b: f64) -> Self::Output { mul(1./b, self) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

impl FnOnce<(&[f64],&[&[f64]], &mut [f64])> for Subroutine {
	type Output = (); extern "rust-call" fn call_once(mut self, args: (&[f64], &[&[f64]], &mut [f64])) -> Self::Output { self.call_mut(args) } }
impl FnMut<(&[f64],&[&[f64]])> for Subroutine { extern "rust-call" fn call_mut(&mut self, args: (&[f64], &[&[f64]], &mut [f64])) -> Self::Output  { self.call(args) } }
impl Fn<(&[f64], &[&[f64]], &mut [f64])> for Subroutine {
	extern "rust-call" fn call(&self, (arguments, arrays, output): (&[f64], &[&[f64]], &mut [f64])) -> Self::Output  {
		assert!(arguments.len()+arrays.len() == self.parameters.len());
		struct State<'t> {
			arguments: &'t [f64],
			arrays: &'t [&'t [f64]],
			definitions: linear_map::LinearMap<Value, f64>,
			output: &mut [f64]
		}
		impl State<'_> {
			fn eval(&mut self, expression: &Expression) -> f64 {
				use Expression::*;
				match expression {
					&Literal(value) => value,
					&Use(id) => {
						if id.0 < self.arguments.len() { self.arguments[id.0] }
						else if id.0 < self.arguments.len()+self.arrays.len() { unreachable!() }
						else { self.definitions[&id] }
					}
					&Index { base, index } => {
						assert!(base.0 >= self.arguments.len(), "{:?} {:?}", base.0, self.arguments);
						let base = base.0 - self.arguments.len();
						self.arrays[base][index]
					}
					Neg(x) => -self.eval(x),
					Add(a, b) => self.eval(a) + self.eval(b),
					Sub(a, b) => self.eval(a) - self.eval(b),
					Less(a, b) => if self.eval(a) < self.eval(b) { 1. } else { 0. },
					Mul(a, b) => self.eval(a) * self.eval(b),
					Div(a, b) => self.eval(a) / self.eval(b),
					Call { function: "exp2", arguments } => f64::exp2(self.eval(&arguments[0])),
					Call { function: "log2", arguments } => f64::log2(self.eval(&arguments[0])),
					Block { statements, result } => {
						self.run(statements);
						let result = self.eval(result);
						for statement in statements.iter() { if let Statement::Definition { id, .. } = statement { self.definitions.remove(id).unwrap(); } else { unreachable!() } }
						result
					}
					e => panic!("{:?}", e),
				}
			}
			fn run(&mut self, statements: &[Statement]) {
				for statement in statements {
					use Statement::*;
					match statement {
						ConditionalStatement { condition, consequent, alternative } => {
							if self.eval(condition) != 0. { self.run(consequent); } else { self.run(alternative); }
						},
						Definition { id, value } => {
							let value = self.eval(value);
							assert!(self.definitions.insert(*id, value).is_none());//, "{} {:?}", id.0, &self.definitions)
						},
						Output { index, value } => self.output[*index] = self.eval(value),
					}
				}
			}
		}
		let mut state = State{arguments, arrays, definitions: default(), output};
		state.run(&self.statements);
	}
}


use itertools::Itertools;

impl std::fmt::Debug for Expression {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		use Expression::*;
		match self {
			&Literal(v) => {
				if v == f64::floor(v) {
					if v >= 1e4 { write!(f, "{:e}", v) }
					else { (v as i64).fmt(f) }
				}
				else { use std::fmt::Display; float_pretty_print::PrettyPrintFloat(v).fmt(f) }
			}
			Use(v) => write!(f, "#{}", v.0),
			Index { base, index } => write!(f, "{}[{}]", base.0, index),
			Neg(x) => write!(f, "-{:?}", x),
			Add(a, b) => write!(f, "({:?} + {:?})", a, b),
			Sub(a, b) => write!(f, "({:?} - {:?})", a, b),
			Less(a, b) => write!(f, "{:?} < {:?}", a, b),
			Mul(a, b) => write!(f, "({:?} * {:?})", a, b),
			Div(a, b) => write!(f, "({:?} / {:?})", a, b),
			Call { function, arguments } => write!(f, "{}({:?})", function, arguments.iter().format(", ")),
			Block { statements, result } => write!(f, "{{ {:?}; {:?} }}", statements.iter().format(";\n"), result),
		}
	}
}

impl std::fmt::Debug for Statement {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		use Statement::*;
		match self {
			Definition { id, value } => write!(f, "#{} = {:?}", id.0, value),
			Output { index, value } => write!(f, "@{} = {:?}", index, value),
			ConditionalStatement { condition, consequent, alternative } => write!(f, "if {:?} {{\n{:?}\n}} else {{\n{:?}\n}}", condition, consequent.iter().format("\n"), alternative.iter().format("\n")),
		}
	}
}

impl std::fmt::Debug for Subroutine {
	fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
		write!(f, "({}) -> {:?} {{\n{:?}\n}}", self.parameters.iter().format(", "), self.output, self.statements.iter().format("\n"))
	}
}

#[macro_export] macro_rules! stringify{ [$($parameter:ident),*] => ([$(std::stringify!($parameter)),*]) }
pub fn parameters<const N: usize>(parameters: [&'static str; N]) -> ([&'static str; N], [Value; N]) { (parameters, iter::from_iter(iter::ConstRange::<N>).map(|id| Value(id))) }

pub fn exp2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "exp2", arguments: box_([x.into()]) } }
pub fn log2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "log2", arguments: box_([x.into()]) } }

pub fn dot(c: &[f64], v: Value) -> Expression {
	c.iter().enumerate().fold(None, |sum, (i,&c)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + index(v, i), None => index(v, i)}) }
		else if c == -1. { Some(match sum { Some(sum) => sum - index(v, i), None => -index(v, i)}) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => c * index(v, i) + sum, None => c * index(v, i) }) }
	).unwrap()
}

pub fn product_of_exponentiations(c: &[impl Copy+Into<i16>], v: Value) -> Expression {
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
