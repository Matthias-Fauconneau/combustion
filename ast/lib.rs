#![allow(incomplete_features,non_snake_case)]
#![feature(unboxed_closures,fn_traits,in_band_lifetimes,associated_type_bounds,format_args_capture,array_map,if_let_guard)]
fn box_<T>(t: T) -> Box<T> { Box::new(t) }
#[macro_export] macro_rules! let_ { { $p:pat = $e:expr => $b:block } => { if let $p = $e { $b } else { unreachable!() } } }

#[derive(PartialEq, Eq, Debug, Clone)] pub struct Value(pub usize);

#[derive(Debug, PartialEq, Eq, Clone, Copy)] pub enum Type { I32, F32, F64 }

#[derive(Debug)] pub enum Expression {
	I32(u32),
	F32(f32),
	F64(f64),
	Float(f64),
	Value(Value),
	Cast(Type, Box<Expression>),
	And(Box<Expression>, Box<Expression>),
	Or(Box<Expression>, Box<Expression>),
	IShLImm(Box<Expression>, u8),
	UShRImm(Box<Expression>, u8),
	IAdd(Box<Expression>, Box<Expression>),
	ISub(Box<Expression>, Box<Expression>),
	//Load(Variable),
	Neg(Box<Expression>),
	Max(Box<Expression>, Box<Expression>),
	Add(Box<Expression>, Box<Expression>),
	Sub(Box<Expression>, Box<Expression>),
	LessOrEqual(Box<Expression>, Box<Expression>),
	Mul(Box<Expression>, Box<Expression>),
	MulAdd(Box<Expression>, Box<Expression>, Box<Expression>),
	Div(Box<Expression>, Box<Expression>),
	FPromote(Box<Expression>),
	FDemote(Box<Expression>),
	FCvtToSInt(Box<Expression>),
	FCvtFromSInt(Box<Expression>),
	Sqrt(Box<Expression>),
	Exp(Box<Expression>),
	Ln { x0: f64, x: Box<Expression> },
	Block { statements: Box<[Statement]>, result: Box<Expression> },
}

#[derive(Debug)] pub enum Statement {
	Value { id: Value, value: Expression },
	Select { condition: Expression, true_exprs: Box<[Expression]>, false_exprs: Box<[Expression]>, results: Box<[Value]> },
	//Display(Value)
}

//pub fn u32(integer: u32) -> Expression { Expression::I32(integer) }
#[track_caller] pub fn float(float: f64) -> Option<Expression> { /*assert!((float as f32).is_finite(), "{float:e}");*/ (float as f32).is_finite().then(|| Expression::Float(float)) }
/*pub fn cast(to: Type, x: impl Into<Expression>) -> Expression { Expression::Cast(to, box_(x.into())) }
pub fn and(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::And(box_(a.into()), box_(b.into())) }
pub fn or(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Or(box_(a.into()), box_(b.into())) }
pub fn ishl_imm(a: impl Into<Expression>, b: u8) -> Expression { Expression::IShLImm(box_(a.into()), b) }
pub fn ushr_imm(a: impl Into<Expression>, b: u8) -> Expression { Expression::UShRImm(box_(a.into()), b) }
pub fn iadd(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::IAdd(box_(a.into()), box_(b.into())) }
pub fn isub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::ISub(box_(a.into()), box_(b.into())) }*/
fn neg(x: impl Into<Expression>) -> Expression { Expression::Neg(box_(x.into())) }
pub fn max(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Max(box_(a.into()), box_(b.into())) }
fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	let [a,b] = [a.into(), b.into()];
	for x in [&a,&b] { if let Expression::Float(x) = x { let x = *x as f32; assert!(x.is_finite() && x != 0.); } }
	Expression::Add(box_(a), box_(b))
}
fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	let [a,b] = [a.into(), b.into()];
	for x in [&a,&b] { if let Expression::Float(x) = x { let x = *x as f32; assert!(x.is_finite() && x != 0.); } }
	Expression::Sub(box_(a), box_(b))
}
pub fn less_or_equal(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::LessOrEqual(box_(a.into()), box_(b.into())) }
#[track_caller] fn mul<A:Into<Expression>, B:Into<Expression>>(a: A, b: B) -> Expression {
	let [a,b] = [a.into(), b.into()];
	for x in [&a,&b] { if let Expression::Float(x) = x { let x = *x as f32; assert!(x.is_finite() && x != 0. && x != 1., "{x}"); } }
	if let [Expression::Float(a), Expression::Float(b)] = [&a,&b] { Expression::Float(a*b) }
	else { Expression::Mul(box_(a), box_(b)) }
}
fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	let [a,b] = [a.into(), b.into()];
	for x in [&a,&b] { if let Expression::Float(x) = x { let x = *x as f32; assert!(x.is_finite() && x != 0.); } }
	if let Expression::Float(x) = &b { let x = *x as f32; assert!(x != 1.); }
	Expression::Div(box_(a), box_(b))
}
/*pub fn fpromote(x: impl Into<Expression>) -> Expression { Expression::FPromote(box_(x.into())) }
pub fn fdemote(x: impl Into<Expression>) -> Expression { Expression::FDemote(box_(x.into())) }
pub fn fcvt_to_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtToSInt(box_(x.into())) }
pub fn fcvt_from_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtFromSInt(box_(x.into())) }
pub fn fma(a: impl Into<Expression>, b: impl Into<Expression>, c: impl Into<Expression>) -> Expression { Expression::MulAdd(box_(a.into()), box_(b.into()), box_(c.into())) }*/
pub fn sqrt(x: impl Into<Expression>) -> Expression { Expression::Sqrt(box_(x.into())) }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }

impl<E:Into<Expression>> std::ops::Add<E> for Expression { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for Expression { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for Expression { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }

impl From<&Value> for Expression { fn from(value: &Value) -> Expression { Expression::Value(value.clone()) } }
impl std::ops::Neg for &Value { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for &Value { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for &Value { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for &Value { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
//impl<E:Into<Expression>> std::ops::Div<E> for &Value { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }
impl std::ops::Div<Expression> for &Value { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }
impl std::ops::Div<&Value> for &Value { type Output = Expression; fn div(self, b: &Value) -> Self::Output { div(self, b) } }
//impl std::ops::Div<f32> for &Value { type Output = Expression; fn div(self, b: f32) -> Self::Output { mul(1./b, self) } }
impl std::ops::Div<f64> for &Value { type Output = Expression; #[track_caller] fn div(self, b: f64) -> Self::Output { mul(1./b, self) } }

//impl From<f32> for Expression { fn from(v: f32) -> Self { Expression::F32(v) } }
impl From<f64> for Expression { fn from(v: f64) -> Self { float(v).unwrap() } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }
impl std::ops::Add<&Value> for f64 { type Output = Expression; fn add(self, b: &Value) -> Self::Output { add(self, b) } }
impl std::ops::Sub<&Value> for f64 { type Output = Expression; fn sub(self, b: &Value) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<&Value> for f64 { type Output = Expression; #[track_caller] fn mul(self, b: &Value) -> Self::Output { assert!(self.is_finite() && self != 0. && self != 1., "{self}"); mul(self, b) } }
impl std::ops::Div<&Value> for f64 { type Output = Expression; fn div(self, b: &Value) -> Self::Output { div(self, b) } }

//type FunctionBuilder = Vec<String>;
//impl FunctionBuilder { pub fn new(input: &[&str]) -> Self { Self { values: input.iter().map(|s| s.to_string()).collect() } } }
pub struct Block<'t> {
	pub statements: Vec<Statement>,
	pub values: &'t mut Vec<String>,
}
pub fn push(s: Statement, block: &mut Block) { block.statements.push(s) }
impl Block<'t> {
	pub fn new(values: &'t mut Vec<String>) -> Self { Self{statements: vec![], values} }
	pub fn block(&mut self, build: impl Fn(&mut Block)->Expression) -> Expression {
		let mut block = Block{ statements: vec![], values: &mut self.values };
		let result = build(&mut block);
		Expression::Block { statements: block.statements.into(), result: box_(result) }
	}
	pub fn value(&mut self, debug: String) -> Value {
		let id = Value(self.values.len());
		self.values.push(debug);
		id
	}
}
pub fn def(value: impl Into<Expression>, block: &mut Block, debug: String) -> Value {
	let id = block.value(debug);
	push(Statement::Value{id: id.clone(), value: value.into()}, block);
	id
}
#[macro_export] macro_rules! l { ($f:ident $e:expr) => ( def($e, $f, format!("{}:{}: {}", file!(), line!(), stringify!($e))) ) }
/*pub fn display<const N: usize>(values: [Value; N], f: &mut Block) -> [Value; N] {
	f.statements.extend(values.iter().cloned().map(Statement::Display));
	values
}
#[macro_export] macro_rules! dbg { ($f:ident; $($id:ident),* ) => { let [$(ref $id),*] = display([$(l!($f $id)),*], $f); } }*/

pub struct Function {
	pub input: usize,
	pub statements: Box<[Statement]>,
	pub output: Box<[Expression]>,
	pub values: Box<[String]>,
}

use std::iter::{Sum, Product};

// cannot impl Sum<Option<Expression>> for Option<Expression> as std:iter_arith_traits_option defines impl Sum<Option> for Option<Sum>
// Forgetting .filter_map(|x| x) will compile with unexpected results
pub fn Σ(iter: impl IntoIterator<Item:Into<Expression>>) -> Option<Expression> {
	iter.into_iter().map(|e| e.into()).filter(|e| if let Expression::Float(e) = e {*e!=0.} else {true}).reduce(add)
}
impl Sum<Expression> for Option<Expression> { fn sum<I:Iterator<Item=Expression>>(iter: I) -> Self { Σ(iter) } }
impl<E:Into<Expression>> Sum<E> for Expression { fn sum<I:Iterator<Item=E>>(iter: I) -> Self { Σ(iter).unwrap() } }
pub fn sum(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().sum() }

pub fn Π(iter: impl IntoIterator<Item:Into<Expression>>) -> Option<Expression> {
	iter.into_iter().map(|e| e.into()).filter(|e| if let Expression::Float(e) = e {*e!=1.} else {true}).reduce(mul)
}
impl Product<Expression> for Option<Expression> { fn product<I:Iterator<Item=Expression>>(iter: I) -> Self { Π(iter) } }
impl<E:Into<Expression>> Product<E> for Expression { fn product<I:Iterator<Item=E>>(iter: I) -> Self { Π(iter).unwrap() } }
pub fn product(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().product() }

pub fn zdot(iter: impl IntoIterator<Item=(f64, impl Into<Expression>)>) -> Option<Expression> {
	iter.into_iter().fold(None, |sum, (c, e)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + e, None => e.into() }) }
		else if c == -1. { Some(match sum { Some(sum) => sum - e, None => neg(e) }) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => mul(float(c).unwrap(), e) + sum, None => mul(float(c).unwrap(), e) }) }
	)
}
#[track_caller] pub fn dot(c: &[f64], v: impl IntoIterator<Item:Into<Expression>>) -> Option<Expression> { zdot(c.iter().copied().zip(v)) }

pub fn exp_approx(x: impl Into<Expression>, f: &mut Block) -> Expression { //e-12 (19*,1/) (-9->-7)
	let ref x = l!(f (1./2048.)*x.into());
	let ref x2 = l!(f x*x);
	let ref x3 = l!(f x2*x);
	let ref a = l!(f 1.+(3./28.)*x2+(1./1680.)*x3*x);
	let ref b = l!(f (1./2.)*x+(1./84.)*x3);
	let sq = |x,f:&mut Block| { let ref x=l!(f x); x*x };
	sq(sq(sq(sq(sq(sq(sq(sq(sq(sq(sq((a+b) / (a-b),f),f),f),f),f),f),f),f),f),f),f)
}
pub fn ln_approx(x0: f64, x: impl Into<Expression>, f: &mut Block) -> Expression { // -5
	let x = (1./x0)*x.into();
	let ref x = l!(f sqrt(sqrt(sqrt(sqrt(x)))));
	let ref x = l!(f (x-1.)/(x+1.));
	let ref x2 = l!(f x*x);
	let ref x4 = l!(f x2*x2);
	let ref x6 = l!(f x4*x2);
	f64::ln(x0) + (16.*2.)*x * (1. + (1./3.)*x2 + (1./5.)*x4 + (1./7.)*x4*x2 + (1./9.)*x6*x2)
}

pub fn exp(x: impl Into<Expression>, f: &mut Block) -> Expression { if true { Expression::Exp(box_(x.into())) } else { exp_approx(x, f) } }
pub fn ln(x0: f64, x: impl Into<Expression>, f: &mut Block) -> Expression { if true { Expression::Ln{x0, x: box_(x.into())} } else { ln_approx(x0, x, f) } }

#[allow(non_camel_case_types)] pub type float = f32;

pub mod interpret;
