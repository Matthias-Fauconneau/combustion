#![feature(unboxed_closures, default_free_fn, fn_traits)]
//, associated_type_bounds, once_cell, , , in_band_lifetimes, array_methods, array_map, trait_alias)]
fn box_<T>(t: T) -> Box<T> { Box::new(t) }
#[macro_export] macro_rules! stringify{ [$($id:ident),*] => ([$(std::stringify!($id)),*]) }
use std::default::default;

#[derive(PartialEq, Eq)] pub struct Value(pub usize);
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

pub fn r#use(index: &Value) -> Expression { Expression::Use(Value(index.0)) }

pub fn index(base: &Value, index: usize) -> Expression { Expression::Index { base: Value(base.0), index } }

fn neg(x: impl Into<Expression>) -> Expression { Expression::Neg(box_(x.into())) }

fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Add(box_(a.into()), box_(b.into())) }
fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Sub(box_(a.into()), box_(b.into())) }
pub fn less(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Less(box_(a.into()), box_(b.into())) }
fn mul(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Mul(box_(a.into()), box_(b.into())) }
fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Div(box_(a.into()), box_(b.into())) }

//pub fn sq(x: Value) -> Expression { Expression::Mul(box_(r#use(x)), box_(r#use(x))) }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }

impl std::ops::Add<Expression> for Expression { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for Expression { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for Expression { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

impl From<&Value> for Expression { fn from(value: &Value) -> Expression { r#use(value) } }
impl<E:Into<Expression>> std::ops::Mul<E> for &Value { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for &Value { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }
impl std::ops::Div<&Value> for Expression { type Output = Expression; fn div(self, b: &Value) -> Self::Output { div(self, b) } }

impl From<f64> for Expression { fn from(v: f64) -> Self { Self::Literal(v) } }
impl std::ops::Add<f64> for Expression { type Output = Expression; fn add(self, b: f64) -> Self::Output { add(self, b) } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<f64> for Expression { type Output = Expression; fn sub(self, b: f64) -> Self::Output { sub(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Mul<&Value> for f64 { type Output = Expression; fn mul(self, b: &Value) -> Self::Output { mul(self, b) } }
impl std::ops::Div<f64> for Expression { type Output = Expression; fn div(self, b: f64) -> Self::Output { mul(1./b, self) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }

#[must_use] fn define(id: Value, value: impl Into<Expression>) -> Statement { Statement::Definition{ id, value: value.into() } }
#[must_use] pub fn output(index: usize, value: impl Into<Expression>) -> Statement { Statement::Output{ index, value: value.into() } }

#[derive(derive_more::Deref,derive_more::DerefMut)] pub struct Block {
	base: usize,
	#[deref]#[deref_mut] statements: Vec<Statement>,
}
impl Block {
	pub fn new<const V: usize, const A: usize>(_: &Parameters<V,A>) -> Self { Self { base: V+A, statements: vec![] } }
	pub fn block(&self, build: impl Fn(&mut Block)->Expression) -> Expression {
		let mut block = Block{ base: self.base+self.statements.len(), statements: vec![] };
		let result = build(&mut block);
		Expression::Block { statements: block.statements.into(), result: box_(result) }
	}
	pub fn def(&mut self, value: Expression) -> Value {
		let id = Value(self.base+self.statements.len());
		self.statements.push(define(Value(id.0), value));
		id
	}
}
impl FnOnce<(Expression,)> for Block { type Output = Value; extern "rust-call" fn call_once(mut self, args: (Expression,)) -> Self::Output { self.call_mut(args) }}
impl FnMut<(Expression,)> for Block { extern "rust-call" fn call_mut(&mut self, (value,): (Expression,)) -> Self::Output  { self.def(value) } }
impl From<Block> for Box<[Statement]> { fn from(b: Block) -> Self { b.statements.into() } }

pub struct Parameters<const V: usize, const A: usize> {
	pub value: [&'static str; V],
	pub array: [&'static str; A],
}
pub struct Subroutine<const V: usize, const A: usize> {
	pub parameters: Parameters<V, A>,
	pub output: usize,
	pub statements: Box<[Statement]>
}

impl<const V: usize, const A: usize> FnOnce<([f64; V], [&[f64]; A], &mut [f64])> for Subroutine<V,A> {
	type Output = ();
	extern "rust-call" fn call_once(mut self, args: ([f64; V], [&[f64]; A], &mut [f64])) -> Self::Output { self.call_mut(args) }
}
impl<const V: usize, const A: usize> FnMut<([f64; V], [&[f64]; A], &mut [f64])> for Subroutine<V,A> {
	extern "rust-call" fn call_mut(&mut self, args: ([f64; V], [&[f64]; A], &mut [f64])) -> Self::Output { self.call(args) }
}
impl<const V: usize, const A: usize> Fn<([f64; V], [&[f64]; A], &mut [f64])> for Subroutine<V,A> {
	extern "rust-call" fn call(&self, (ref value_arguments, ref array_arguments, output): ([f64; V], [&[f64]; A], &mut [f64])) -> Self::Output {
		struct State<'t> {
			value_arguments: &'t [f64],
			array_arguments: &'t [&'t [f64]],
			definitions: linear_map::LinearMap<Value, f64>,
			output: &'t mut [f64]
		}
		impl State<'_> {
			fn eval(&mut self, expression: &Expression) -> f64 {
				use Expression::*;
				match expression {
					&Literal(value) => value,
					Use(id) => {
						if id.0 < self.value_arguments.len() { self.value_arguments[id.0] }
						else if id.0 < self.value_arguments.len()+self.array_arguments.len() { unreachable!() }
						else { self.definitions[&id] }
					}
					Index { base, index } => {
						assert!(base.0 >= self.value_arguments.len(), "{:?} {:?}", base.0, self.value_arguments);
						let base = base.0 - self.value_arguments.len();
						self.array_arguments[base][*index]
					}
					Neg(x) => -self.eval(x),
					Add(a, b) => self.eval(a) + self.eval(b),
					Sub(a, b) => self.eval(a) - self.eval(b),
					Less(a, b) => if self.eval(a) < self.eval(b) { 1. } else { 0. },
					Mul(a, b) => self.eval(a) * self.eval(b),
					Div(a, b) => self.eval(a) / self.eval(b),
					Call { function: "sq", arguments } => num::sq(self.eval(&arguments[0])),
					Call { function: "sqrt", arguments } => f64::sqrt(self.eval(&arguments[0])),
					Call { function: "exp2", arguments } => f64::exp2(self.eval(&arguments[0])),
					Call { function: "log2", arguments } => f64::log2(self.eval(&arguments[0])),
					Block { statements, result } => {
						self.run(statements);
						let result = self.eval(result);
						for statement in statements.iter() { if let Statement::Definition { id, .. } = statement { self.definitions.remove(id).unwrap(); } else { unreachable!() } }
						result
					}
					_ => unimplemented!(),
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
							assert!(self.definitions.insert(Value(id.0), value).is_none());
						},
						Output { index, value } => self.output[*index] = self.eval(value),
					}
				}
			}
		}
		let mut state = State{value_arguments, array_arguments, definitions: default(), output};
		state.run(&self.statements);
	}
}

pub fn wrap<const V: usize, const A: usize>(f: Subroutine<V, A>) -> impl Fn([f64; V], [&[f64]; A])->Vec<f64> {
	move |values, arrays| { let mut output = vec![0.; f.output]; f(values, arrays, &mut output); output }
}

use iter::{ConstRange, ConstSizeIterator};
pub fn parameters<const V: usize, const A: usize>(value: [&'static str; V], array: [&'static str; A]) -> (Parameters<V, A>, [Value; V], [Value; A]) {
	(Parameters{value, array}, ConstRange.map(|id| Value(id)).collect(), ConstRange.map(|id| Value(V+id)).collect())
}
#[macro_export] macro_rules! parameters { ($([$(ref $parameter:ident),*]),*) => ( $crate::parameters( $([$(std::stringify!($parameter)),*]),* ) ) }

pub fn sum(iter: impl Iterator<Item=Expression>) -> Expression { iter.reduce(add).unwrap() }

pub fn sq(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "sq", arguments: box_([x.into()]) } }
pub fn sqrt(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "sqrt", arguments: box_([x.into()]) } }
pub fn exp2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "exp2", arguments: box_([x.into()]) } }
pub fn log2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "log2", arguments: box_([x.into()]) } }
