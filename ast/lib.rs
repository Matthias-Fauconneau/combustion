#![feature(unboxed_closures, default_free_fn, fn_traits, in_band_lifetimes, associated_type_bounds)]
fn box_<T>(t: T) -> Box<T> { Box::new(t) }
#[macro_export] macro_rules! let_ { { $p:pat = $e:expr => $b:block } => { if let $p = $e { $b } else { unreachable!() } } }
use std::default::default;

#[derive(PartialEq, Eq, Debug)] pub struct Value(pub usize);
#[derive(PartialEq, Eq, Debug)] pub struct Variable(pub usize);

#[derive(Debug)] pub enum Expression {
	Literal(f64),
	Use(Value),
	Load(Variable),
	Neg(Box<Expression>),
	Add(Box<Expression>, Box<Expression>),
	Sub(Box<Expression>, Box<Expression>),
	Less(Box<Expression>, Box<Expression>),
	Mul(Box<Expression>, Box<Expression>),
	Div(Box<Expression>, Box<Expression>),
	Call { function: &'static str, arguments: Box<[Expression]> },
	Block { statements: Box<[Statement]>, result: Box<Expression> },
}

#[derive(Debug)] pub enum Statement {
	Define { id: Value, value: Expression },
	Store { variable: Variable, value: Expression },
	Output { index: usize, value: Expression },
	Branch { condition: Expression, consequent: Box<[Statement]>, alternative: Box<[Statement]> },
}

pub fn r#use(value: &Value) -> Expression { Expression::Use(Value(value.0)) }

fn neg(x: impl Into<Expression>) -> Expression { Expression::Neg(box_(x.into())) }

fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Add(box_(a.into()), box_(b.into())) }
fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Sub(box_(a.into()), box_(b.into())) }
pub fn less(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Less(box_(a.into()), box_(b.into())) }
fn mul(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Mul(box_(a.into()), box_(b.into())) }
fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Div(box_(a.into()), box_(b.into())) }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }

impl<E:Into<Expression>> std::ops::Add<E> for Expression { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for Expression { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for Expression { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }

impl From<&Value> for Expression { fn from(value: &Value) -> Expression { r#use(value) } }
impl std::ops::Neg for &Value { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for &Value { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for &Value { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for &Value { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for &Value { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }

impl From<f64> for Expression { fn from(v: f64) -> Self { Self::Literal(v) } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }
impl std::ops::Add<&Value> for f64 { type Output = Expression; fn add(self, b: &Value) -> Self::Output { add(self, b) } }
impl std::ops::Sub<&Value> for f64 { type Output = Expression; fn sub(self, b: &Value) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<&Value> for f64 { type Output = Expression; fn mul(self, b: &Value) -> Self::Output { mul(self, b) } }
impl std::ops::Div<&Value> for f64 { type Output = Expression; fn div(self, b: &Value) -> Self::Output { div(self, b) } }

#[must_use] fn define(id: Value, value: impl Into<Expression>) -> Statement { Statement::Define{ id, value: value.into() } }
#[must_use] pub fn store(variable: &Variable, value: impl Into<Expression>) -> Statement { Statement::Store{ variable: Variable(variable.0), value: value.into() } }
#[must_use] pub fn output(index: usize, value: impl Into<Expression>) -> Statement { Statement::Output{ index, value: value.into() } }

pub struct FunctionBuilder {
	input: usize,
	variables: usize,
}
impl FunctionBuilder {
	pub fn new(input: &[Value]) -> Self { Self { input: input.len(), variables: 0 } }
}

#[derive(derive_more::Deref,derive_more::DerefMut)] pub struct Block<'t> {
	base: usize,
	#[deref]#[deref_mut] statements: Vec<Statement>,
	variables: &'t mut usize,
}
impl Block<'t> {
	pub fn new(f: &'t mut FunctionBuilder) -> Self { Self { base: f.input, statements: vec![], variables: &mut f.variables } }
	pub fn block(&mut self, build: impl Fn(&mut Block)->Expression) -> Expression {
		let mut block = Block{ base: self.base+self.statements.len(), statements: vec![], variables: self.variables };
		let result = build(&mut block);
		self.base += block.statements.len();
		Expression::Block { statements: block.statements.into(), result: box_(result) }
	}
	pub fn def(&mut self, value: Expression) -> Value {
		let id = Value(self.base+self.statements.len());
		self.statements.push(define(Value(id.0), value));
		id
	}
	pub fn decl(&mut self) -> Variable {
		let id = Variable(*self.variables);
		*self.variables += 1;
		id
	}
	pub fn load(&mut self, variable: Variable) -> Value { self.def(Expression::Load(variable)) }
}
impl FnOnce<(Expression,)> for Block<'_> { type Output = Value; extern "rust-call" fn call_once(mut self, args: (Expression,)) -> Self::Output { self.call_mut(args) }}
impl FnMut<(Expression,)> for Block<'_> { extern "rust-call" fn call_mut(&mut self, (value,): (Expression,)) -> Self::Output  { self.def(value) } }
impl From<Block<'_>> for Box<[Statement]> { fn from(b: Block) -> Self { b.statements.into() } }

pub struct Function {
	pub input: usize,
	pub output: usize,
	pub variables: usize,
	pub statements: Box<[Statement]>,
}
impl Function {
	pub fn new(output: usize, statements: Box<[Statement]>, FunctionBuilder{input, variables}: FunctionBuilder) -> Self { Function{input, output, variables, statements} }
}

impl FnOnce<(&[f64], &mut [f64])> for Function {
	type Output = ();
	extern "rust-call" fn call_once(mut self, args: (&[f64], &mut [f64])) -> Self::Output { self.call_mut(args) }
}
impl FnMut<(&[f64], &mut [f64])> for Function {
	extern "rust-call" fn call_mut(&mut self, args: (&[f64], &mut [f64])) -> Self::Output { self.call(args) }
}
impl Fn<(&[f64], &mut [f64])> for Function {
	extern "rust-call" fn call(&self, (input, output): (&[f64], &mut [f64])) -> Self::Output {
		struct State<'t> {
			input: &'t [f64],
			output: &'t mut [f64],
			definitions: linear_map::LinearMap<Value, f64>,
			variables: Box<[Option<f64>]>,
		}
		impl State<'_> {
			fn eval(&mut self, expression: &Expression) -> f64 {
				use Expression::*;
				match expression {
					&Literal(value) => value,
					Use(id) => {
						if id.0 < self.input.len() { self.input[id.0] }
						else { self.definitions[&id] }
					}
					Load(variable) => { self.variables[variable.0].unwrap() }
					Neg(x) => -self.eval(x),
					Add(a, b) => self.eval(a) + self.eval(b),
					Sub(a, b) => self.eval(a) - self.eval(b),
					Less(a, b) => if self.eval(a) < self.eval(b) { 1. } else { 0. },
					Mul(a, b) => self.eval(a) * self.eval(b),
					Div(a, b) => self.eval(a) / self.eval(b),
					Call { function, arguments } => match *function {
						"max" => f64::max(self.eval(&arguments[0]), self.eval(&arguments[1])),
						"sqrt" => f64::sqrt(self.eval(&arguments[0])),
						"exp2" => f64::exp2(self.eval(&arguments[0])),
						"log2" => f64::log2(self.eval(&arguments[0])),
						function => panic!("{}", function)
					},
					Block { statements, result } => {
						self.run(statements);
						let result = self.eval(result);
						for statement in statements.iter() { if let Statement::Define { id, .. } = statement { self.definitions.remove(id).unwrap(); } else { unreachable!() } }
						result
					},
					//e => panic!("{:?}", e),
				}
			}
			fn run(&mut self, statements: &[Statement]) {
				for statement in statements {
					use Statement::*;
					match statement {
						Define { id, value } => {
							let value = self.eval(value);
							assert!(self.definitions.insert(Value(id.0), value).is_none());
						},
						Store { variable, value } => {
							let value = self.eval(value);
							assert!(self.variables[variable.0].replace(value).is_none());
						}
						Output { index, value } => self.output[*index] = self.eval(value),
						Branch { condition, consequent, alternative } => {
							if self.eval(condition) != 0. { self.run(consequent); } else { self.run(alternative); }
						},
					}
				}
			}
		}
		assert!(input.len() == self.input);
		let mut state = State{input, output, definitions: default(), variables: vec![None; self.variables].into()};
		state.run(&self.statements);
	}
}

pub fn wrap(f: Function) -> impl Fn(&[f64]) -> Box<[f64]> { move |input| { let mut output = vec![0.; f.output].into_boxed_slice(); f(input, &mut output); output } }

impl<E:Into<Expression>> std::iter::Sum<E> for Expression { fn sum<I:Iterator<Item=E>>(iter: I) -> Self { iter.into_iter().map(|e| e.into()).reduce(add).unwrap() } }
pub fn sum(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().sum() }

pub fn max(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Call{ function: "max", arguments: box_([a.into(), b.into()]) } }
pub fn sqrt(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "sqrt", arguments: box_([x.into()]) } }
pub fn exp2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "exp2", arguments: box_([x.into()]) } }
pub fn log2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "log2", arguments: box_([x.into()]) } }

#[cfg(feature="num")] impl num::Sqrt for Expression { fn sqrt(self) -> Self { sqrt(self) } }
#[cfg(feature="num")] impl num::Log for Expression { fn log2(self) -> Self { log2(self) } }

pub fn idot<'t>(iter: impl IntoIterator<Item=(f64, &'t Value)>) -> Expression {
	iter.into_iter().fold(None, |sum, (c, e)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + e, None => e.into() }) }
		else if c == -1. { Some(match sum { Some(sum) => sum - e, None => -e }) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => c * e + sum, None => c * e }) }
	).unwrap()
}
pub fn dot(c: &[f64], v: &[Value]) -> Expression { idot(c.iter().copied().zip(v)) }
