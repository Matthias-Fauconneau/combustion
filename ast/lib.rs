#![allow(incomplete_features,non_snake_case,mixed_script_confusables)]
#![feature(unboxed_closures,fn_traits,in_band_lifetimes,associated_type_bounds,format_args_capture,if_let_guard,type_ascription)]
//#![recursion_limit="5"]
fn box_<T>(t: T) -> Box<T> { Box::new(t) }
#[macro_export] macro_rules! let_ { { $p:pat = $e:expr => $b:block } => { if let $p = $e { $b } else { unreachable!() } } }

#[derive(PartialEq,Eq,Hash,Debug,Clone,Copy)] pub struct Value(pub usize);
impl From<&Value> for Value { fn from(x: &Value) -> Value { x.clone() } }
impl From<&mut Value> for Value { fn from(x: &mut Value) -> Value { x.clone() } }

#[derive(PartialEq,Eq,Debug,Clone,Copy)] pub enum Type { I32, F32, F64 }

use ordered_float::NotNan;
type R32 = NotNan<f32>;
#[derive(Eq,Hash,Debug,Clone,Copy)] pub struct R64(NotNan<f64>);
impl PartialEq for R64 { fn eq(&self, b: &R64) -> bool { R32::new(*self.0 as _) == R32::new(*b.0 as _) } }
//impl PartialEq for R64 { fn eq(&self, b: &R64) -> bool { self.0 == b.0 } }
impl R64 { fn new(val: f64) -> Self { Self(NotNan::new(val).unwrap()) } }
//impl std::fmt::Display for R64 { fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> { self.0.fmt(f) } }
impl From<R64> for f64 { fn from(o: R64) -> f64 { o.0.into() } }
use std::ops::{Deref,DerefMut};
//impl Deref for R64 { type Target=NotNan<f64>; fn deref(&self) -> &Self::Target { &self.0 } }
impl Deref for R64 { type Target=f64; fn deref(&self) -> &Self::Target { &self.0 } }

#[derive(PartialEq,Eq,Hash,Debug,Clone)] pub enum Expr/*ExpressionWithoutBlock*/ {
	//I32(u32),
	F32(R32),
	F64(R64),
	Value(Value),
	/*Cast(Type, Box<Expression>),
	And(Box<Expression>, Box<Expression>),
	Or(Box<Expression>, Box<Expression>),
	IShLImm(Box<Expression>, u8),
	UShRImm(Box<Expression>, u8),
	IAdd(Box<Expression>, Box<Expression>),).into_inner
	ISub(Box<Expression>, Box<Expression>),*/
	//Load(Variable),
	Neg(Box<Expression>),
	Max(Box<Expression>, Box<Expression>),
	Add(Box<Expression>, Box<Expression>),
	Sub(Box<Expression>, Box<Expression>),
	LessOrEqual(Box<Expression>, Box<Expression>),
	Mul(Box<Expression>, Box<Expression>),
	//MulAdd(Box<Expression>, Box<Expression>, Box<Expression>),
	Div(Box<Expression>, Box<Expression>),
	/*FPromote(Box<Expression>),
	FDemote(Box<Expression>),
	FCvtToSInt(Box<Expression>),
	FCvtFromSInt(Box<Expression>),*/
	Sqrt(Box<Expression>),
	Exp(Box<Expression>),
	Ln { x0: NotNan<f64>, x: Box<Expression> },
}
impl Expr {
	pub fn is_leaf(&self) -> bool { use Expr::*; matches!(self, F32(_)|F64(_)|Value(_)) }
	pub fn shallow(&self) -> Expr { assert!(self.is_leaf()); self.clone() }
	pub fn to_string(&self, names: &[String]) -> String {
		use Expr::*; match self {
			F64(x) => x.to_string(),
			Value(id) => names[id.0].clone(),
			Add(a, b) => format!("{} + {}", a.to_string(names), b.to_string(names)),
			Sub(a, b) => format!("{} - {}", a.to_string(names), b.to_string(names)),
			Mul(a, b) => format!("({}) * ({})", a.to_string(names), b.to_string(names)),
			Div(a, b) => format!("({}) / ({})", a.to_string(names), b.to_string(names)),
			Exp(x) => format!("exp({})", x.to_string(names)),
			_ => format!("{self:?}")
		}
	}
}
//impl From<&Expr> for Expr { fn from(e: &Expr) -> Self { e.shallow() } }

#[derive(Debug, PartialEq, Eq, Hash)] pub enum Expression {
	Expr(Expr),
	Block { statements: Box<[Statement]>, result: Box<Expression> },
}
impl Deref for Expression { type Target=Expr; #[track_caller] fn deref(&self) -> &Expr { if let Expression::Expr(e) = self { e } else { panic!("block") } } }
impl DerefMut for Expression { fn deref_mut(&mut self) -> &mut Expr { if let Expression::Expr(e) = self { e } else { panic!("block") } } }
impl<E: Into<Expr>> From<E> for Expression { fn from(e: E) -> Expression { Expression::Expr(e.into()) } }
//impl From<&Expression> for Expression { fn from(e: &Expression) -> Expression { Deref::deref(e).shallow().into() } }
impl Default for Expression { fn default() -> Expression { Expr::Value(Value(usize::MAX)).into() } }
impl Clone for Expression { fn clone(&self) -> Expression { Deref::deref(self).clone().into() } }
impl Expression {
	pub fn is_leaf(&self) -> bool { if let Expression::Expr(e) = self { e.is_leaf() } else { false } }
	pub fn to_string(&self, names: &[String]) -> String {
		use itertools::Itertools;
		match self {
			Self::Expr(e) => e.to_string(names),
			Self::Block{statements, result} => format!("{{{} {}}}", statements.iter().map(|x| x.to_string(names)).format(" "), result.to_string(names))
		}
	}
	pub fn visit<T>(&self, mut visitor: impl FnMut(&Self)->T) -> [Option<T>; 2] {
		match self {
			Self::Expr(e) => {use Expr::*; match e {
				F32(_)|F64(_)|Value(_) => [None, None],
				Neg(x)|Sqrt(x)|Exp(x)|Ln{x,..} => [Some(visitor(x)), None],
				Max(a, b)|Add(a, b)|Sub(a, b)|LessOrEqual(a, b)|Mul(a, b)|Div(a, b) => { [visitor(a), visitor(b)].map(Some) }
			}}
			Self::Block{..} => [None, None]//unimplemented!()
		}
	}
	pub fn visit_mut(&mut self, mut visitor: impl FnMut(&mut Self)) {
		match self {
			Self::Expr(e) => {use Expr::*; match e {
				F32(_)|F64(_)|Value(_) => {},
				Neg(x)|Sqrt(x)|Exp(x)|Ln{x,..} => visitor(x),
				Max(a, b)|Add(a, b)|Sub(a, b)|LessOrEqual(a, b)|Mul(a, b)|Div(a, b) => { visitor(a); visitor(b); }
			}}
			Self::Block{..} => unimplemented!()
		}
	}
	pub fn has_block(&self) -> bool { matches!(self, Self::Block{..}) || self.visit(Self::has_block).into_iter().filter_map(|x| x).any(|x| x) }
}

#[derive(Debug,PartialEq,Eq,Hash)] pub enum Statement {
	Value { id: Value, value: Expression },
	Select { condition: Expression, true_exprs: Box<[Expression]>, false_exprs: Box<[Expression]>, results: Box<[Value]> },
	//Display(Value)
}
impl Statement {
	fn to_string(&self, names: &[String]) -> String {
		use Statement::*;
		match self {
			Value{id, value} => format!("{} = {}", names[id.0], value.to_string(names)),
			Select{..} => unimplemented!()
		}
	}
}

pub fn f32(x: f32) -> Option<Expr> { x.is_finite().then(|| Expr::F32(R32::new(x).unwrap()) ) }
pub fn f64(x: f64) -> Option<Expr> { (x as f32).is_finite().then(|| Expr::F64(R64::new(x)) ) }
pub trait From_<F> { fn from(_: F) -> Self; }
//impl From_<f32> for Option<Expr> { fn from(x: f32) -> Option<Expr> { f32(x) } }
impl From_<f64> for Option<Expr> { fn from(x: f64) -> Option<Expr> { f64(x) } }
impl<F> From<F> for Expr where Option<Expr>: From_<F> { fn from(x: F) -> Expr { let x:Option<Expr> = From_::from(x); x.unwrap() } }
//impl From<f32> for Expr where Option<Expr>: From_<F> { fn from(x: F) -> Expr { let x:Option<Expr> = From_::from(x); x.unwrap() } }
impl From<Value> for Expr { fn from(x: Value) -> Expr { Expr::Value(x) } }
impl From<&Value> for Expr { fn from(x: &Value) -> Expr { x.clone().into() } }
impl From<&mut Value> for Expr { fn from(x: &mut Value) -> Expr { (&*x).into() } }

impl Expr {
	pub fn f32(&self) -> Option<f32> { use Expr::*; match self { F32(x) => Some(f32::from(*x)), F64(x) => Some(f64::from(*x) as _), _ => None } }
	pub fn f64(&self) -> Option<f64> { use Expr::*; match self { F32(x) => Some(f32::from(*x) as _), F64(x) => Some(f64::from(*x)), _ => None } }
}

//impl From<&Value> for Expression { fn from(x: &Value) -> Expression { Expr::from(x).into() } }
//impl From<&mut Value> for Expression { fn from(x: &mut Value) -> Expression { Expr::from(&*x).into() } }

/*pub fn cast(to: Type, x: impl Into<Expression>) -> Expression { Expression::Cast(to, box_(x.into())) }
pub fn and(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::And(box_(a.into()), box_(b.into())) }
pub fn or(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Or(box_(a.into()), box_(b.into())) }
pub fn ishl_imm(a: impl Into<Expression>, b: u8) -> Expression { Expression::IShLImm(box_(a.into()), b) }
pub fn ushr_imm(a: impl Into<Expression>, b: u8) -> Expression { Expression::UShRImm(box_(a.into()), b) }
pub fn iadd(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::IAdd(box_(a.into()), box_(b.into())) }
pub fn isub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::ISub(box_(a.into()), box_(b.into())) }*/

fn neg(x: impl Into<Expression>) -> Expression {
	Expr::Neg(box_(x.into())).into()
}

pub fn max(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	Expr::Max(box_(a.into()), box_(b.into())).into()
}

pub fn less_or_equal(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	Expr::LessOrEqual(box_(a.into()), box_(b.into())).into()
}

fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	let [a,b] = [a.into(), b.into()];
	if let [Some(a), Some(b)] = [a.f64(),b.f64()] { (a+b).into() } else {
		if let Some(a) = a.f32() { if a==0. { return b; } }
		if let Some(b) = b.f32() { if b==0. { return a; } }
		Expr::Add(box_(a), box_(b))
	}.into()
}

#[track_caller] fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	let [a,b] = [a.into(), b.into()];
	if let [Some(a), Some(b)] = [a.f64(),b.f64()] { (a-b).into() } else {
		if let Some(a) = a.f32() { if a==0. { return -b; } }
		if let Some(b) = b.f32() { if b==0. { return a; } }
		Expr::Sub(box_(a), box_(b))
	}.into()
}

#[track_caller] fn mul<A:Into<Expression>, B:Into<Expression>>(a: A, b: B) -> Expression {
	let [a,b] = [a.into(), b.into()];
	assert!(!matches!(&*a, Expr::Div(_,_)) && !matches!(&*b, Expr::Div(_,_)), "({a:?}) * ({b:?})");
	if let [Some(a), Some(b)] = [a.f64(),b.f64()] { (a*b).into() } else {
		for x in [&a,&b] { if let Some(x) = x.f32() { if x == 0. { return f64(0.).unwrap().into() } } }
		if let Some(a) = a.f32() { if a==1. { return b; } else if a==-1. { return -b } }
		if let Some(b) = b.f32() { if b==1. { return a; } else if b==-1. { return -a } }
		Expr::Mul(box_(a), box_(b))
	}.into()
}

fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression {
	let [a,b] = [a.into(), b.into()];
	if let Some(b) = b.f64() { return a*(1./b) } else {
		for x in [&a,&b] { if let Some(x) = x.f32() { assert!(x != 0.); } }
		if let Some(x) = b.f32() { assert!(x != 1.); }
		Expr::Div(box_(a), box_(b))
	}.into()
}

pub fn sqrt(x: impl Into<Expression>) -> Expression {
	Expr::Sqrt(box_(x.into())).into()
}

/*pub fn fpromote(x: impl Into<Expression>) -> Expression { Expression::FPromote(box_(x.into())) }
pub fn fdemote(x: impl Into<Expression>) -> Expression { Expression::FDemote(box_(x.into())) }
pub fn fcvt_to_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtToSInt(box_(x.into())) }
pub fn fcvt_from_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtFromSInt(box_(x.into())) }
pub fn fma(a: impl Into<Expression>, b: impl Into<Expression>, c: impl Into<Expression>) -> Expression { Expression::MulAdd(box_(a.into()), box_(b.into()), box_(c.into())) }*/

fn oadd(a: Option<impl Into<Expression>>, b: impl Into<Expression>) -> Expression { if let Some(a) = a { a.into()+b.into() } else { b.into() } }
//fn addo(a: impl Into<Expression>, b: Option<impl Into<Expression>>) -> Expression { b.map(|b| a+b).unwrap_or(a.into()) }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for Expression { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for Expression { type Output = Expression; #[track_caller] fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; #[track_caller] fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for Expression { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }
impl std::ops::Neg for Expr { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for Expr { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for Expr { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expr { type Output = Expression; #[track_caller] fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for Expr { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }
/*impl std::ops::Neg for &Expr { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for &Expr { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for &Expr { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
//impl<E:Into<Expression>> std::ops::Mul<E> for &Expr { type Output = Expression; #[track_caller] fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for &Expr { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }*/
impl std::ops::Neg for Value { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for Value { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for Value { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Value { type Output = Expression; #[track_caller] fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for Value { type Output = Expression; #[track_caller] fn div(self, b: E) -> Self::Output { div(self, b) } }
impl std::ops::Neg for &Value { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for &Value { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for &Value { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for &Value { type Output = Expression; #[track_caller] fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for &Value { type Output = Expression; #[track_caller] fn div(self, b: E) -> Self::Output { div(self, b) } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }
impl std::ops::Add<Expr> for f64 { type Output = Expression; fn add(self, b: Expr) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expr> for f64 { type Output = Expression; fn sub(self, b: Expr) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expr> for f64 { type Output = Expression; fn mul(self, b: Expr) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expr> for f64 { type Output = Expression; fn div(self, b: Expr) -> Self::Output { div(self, b) } }
/*impl std::ops::Add<&Expr> for f64 { type Output = Expression; fn add(self, b: &Expr) -> Self::Output { add(self, b) } }
impl std::ops::Sub<&Expr> for f64 { type Output = Expression; fn sub(self, b: &Expr) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<&Expr> for f64 { type Output = Expression; fn mul(self, b: &Expr) -> Self::Output { mul(self, b) } }
impl std::ops::Div<&Expr> for f64 { type Output = Expression; fn div(self, b: &Expr) -> Self::Output { div(self, b) } }*/
impl std::ops::Add<Value> for f64 { type Output = Expression; fn add(self, b: Value) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Value> for f64 { type Output = Expression; fn sub(self, b: Value) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Value> for f64 { type Output = Expression; #[track_caller] fn mul(self, b: Value) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Value> for f64 { type Output = Expression; fn div(self, b: Value) -> Self::Output { div(self, b) } }
impl std::ops::Add<&Value> for f64 { type Output = Expression; fn add(self, b: &Value) -> Self::Output { add(self, b) } }
impl std::ops::Sub<&Value> for f64 { type Output = Expression; fn sub(self, b: &Value) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<&Value> for f64 { type Output = Expression; fn mul(self, b: &Value) -> Self::Output { mul(self, b) } }
impl std::ops::Div<&Value> for f64 { type Output = Expression; fn div(self, b: &Value) -> Self::Output { div(self, b) } }

impl std::ops::Add<Expression> for Option<Expression> { type Output = Expression; fn add(self, b: Expression) -> Self::Output { oadd(self, b) } }
//impl<E:Into<Expression>> std::ops::Add<Option<E>> for Expression { type Output = Expression; fn add(self, b: E) -> Self::Output { addo(self, b) } }

pub struct Block<'t> {
	pub statements: Vec<Statement>,
	pub names: &'t mut Vec<String>,
}
pub fn push(s: Statement, block: &mut Block) { block.statements.push(s) }
impl Block<'t> {
	pub fn new(names: &'t mut Vec<String>) -> Self { Self{statements: vec![], names} }
	pub fn block(&mut self, build: impl FnOnce(&mut Block)->Expression) -> Expression {
		let mut block = Block{ statements: vec![], names: &mut self.names };
		let result = build(&mut block);
		if block.statements.is_empty() { result } else { Expression::Block { statements: block.statements.into(), result: box_(result) } }
	}
	pub fn value(&mut self, name: String) -> Value {
		let id = Value(self.names.len());
		self.names.push(name);
		id
	}
}
pub fn def(value: impl Into<Expression>, block: &mut Block, name: String) -> Result<Value, f64> {
	let value = value.into();
	if let Expression::Expr(ref e) = value { if let Some(x) = e.f64() { return Err(x) } }
	let id = block.value(name);
	push(Statement::Value{id: id.clone(), value}, block);
	Ok(id)
}
#[macro_export] macro_rules! l {
	($f:ident $e:expr) => ( def($e, $f, format!("{}:{}: {}", file!(), line!(), stringify!($e))).unwrap() );
	($f:ident; $e:expr) => ( l!($f $e).into() );
}
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

// cannot impl Sum<Option<Expression>> for Option<Expression> as std:iter_arith_traits_option defines impl Sum<Option> for Option<Sum>
// Forgetting .filter_map(|x| x) will compile with unexpected results
pub fn Σ(iter: impl IntoIterator<Item:Into<Expression>>) -> Option<Expression> {
	iter.into_iter().map(|e| e.into()).filter(|e| if let Some(x) = e.f32() {x!=0.} else {true}).reduce(add)
}
use std::iter::Sum;
impl Sum<Expression> for Option<Expression> { fn sum<I:Iterator<Item=Expression>>(iter: I) -> Self { Σ(iter) } }
impl<E:Into<Expression>> Sum<E> for Expression { fn sum<I:Iterator<Item=E>>(iter: I) -> Self { Σ(iter).unwrap() } }
pub fn sum(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().sum() }

pub fn Π(iter: impl IntoIterator<Item=Option<impl Into<Expression>>>) -> Option<Expression> {
	iter.into_iter().filter_map(|x| x).map(|e| e.into()).filter(|e| if let Some(x) = e.f32() {x!=1.} else {true}).reduce(mul)
}
//use std::iter::Product;
//impl Product<Expression> for Option<Expression> { fn product<I:Iterator<Item=Expression>>(iter: I) -> Self { Π(iter.map(Some)) } }
//impl<E:Into<Expression>> Product<Option<E>> for Expression { fn product<I:Iterator<Item=E>>(iter: I) -> Self { Π(iter.map(Some)).unwrap() } }
//impl<E:Into<Expression>> Product<E> for Expression { fn product<I:Iterator<Item=E>>(iter: I) -> Self { iter.map(Some).product() } }
//pub fn product(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().product() }

/*pub fn exp_approx(x: impl Into<Expression>, f: &mut Block) -> Expression { //e-12 (19*,1/) (-9->-7)
	let ref x = l!(f (1./2048.)*x.into());
	let ref x2 = l!(f x*x);
	let ref x3 = l!(f x2*x);
	let ref a = l!(f 1.+(3./28.)*x2+(1./1680.)*x3*x);
	let ref b = l!(f (1./2.)*x+(1./84.)*x3);
	let sq = |x,f:&mut Block| { let ref x=l!(f x); x*x };
	sq(sq(sq(sq(sq(sq(sq(sq(sq(sq(sq((a+b) / (a-b),f),f),f),f),f),f),f),f),f),f),f)
}
pub fn ln_approx(x0: f64, x: impl Into<Expression>, f: &mut Block) -> Expression { // -5
	/*let x = (1./x0)*x.into();
	let ref x = l!(f sqrt(sqrt(sqrt(sqrt(x)))));
	let ref x = l!(f (x-1.)/(x+1.));
	let ref x2 = l!(f x*x);
	let ref x4 = l!(f x2*x2);
	let ref x6 = l!(f x4*x2);
	f64::ln(x0) + (16.*2.)*x * (1. + (1./3.)*x2 + (1./5.)*x4 + (1./7.)*x4*x2 + (1./9.)*x6*x2)*/unimplemented!()
}*/

pub fn exp(x: impl Into<Expression>, _f: &mut Block) -> Expression {
	let x = x.into();
	assert!(x.f64().is_none());
	Expr::Exp(box_(x)).into()	//exp_approx(x, f)
}
#[track_caller] pub fn ln(x0: f64, x: impl Into<Expression>, _f: &mut Block) -> Expression {
	let x = x.into();
	if let Some(x) = x.f64() { f64::ln(x).into() }
	else { Expr::Ln{x0: NotNan::new(x0).unwrap(), x: box_(x.into())}.into() } // ln_approx(x0, x, f)
}

pub fn eliminate_common_subexpression(a: &mut Expression, b: &mut Expression, f: &mut Block) {
	if !a.is_leaf() && !b.is_leaf() && *a == *b {
		let common: Expr = l!(f; std::mem::take(a));
		*a = common.shallow().into();
		*b = common.into();
	} else {
		a.visit_mut(|a| eliminate_common_subexpression(a, b, f));
		b.visit_mut(|b| eliminate_common_subexpression(a, b, f));
	}
}

pub fn eliminate_common_subexpressions(a: &mut [Expression], b: &mut [Expression], f: &mut Block) {
	for a in &mut *a { for b in &mut *b { eliminate_common_subexpression(a, b, f); } }
}

//pub fn dot(iter: impl IntoIterator<Item=(f64, impl Into<Expression>)>) -> Option<Expression> { iter.into_iter().map(|(c,e)| c*e.into()).sum() }
//#[track_caller] pub fn dot(c: &[f64], v: impl IntoIterator<Item:Into<Expression>>) -> Option<Expression> { zdot(c.iter().copied().zip(v)) }

#[cfg(feature="interpret")] pub mod interpret;
