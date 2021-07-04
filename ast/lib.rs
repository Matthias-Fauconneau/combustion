#![allow(incomplete_features)]#![feature(unboxed_closures,fn_traits,in_band_lifetimes,associated_type_bounds,format_args_capture,array_map)]
fn box_<T>(t: T) -> Box<T> { Box::new(t) }
#[macro_export] macro_rules! let_ { { $p:pat = $e:expr => $b:block } => { if let $p = $e { $b } else { unreachable!() } } }

#[derive(PartialEq, Eq, Debug, Clone)] pub struct Value(pub usize);
//#[derive(PartialEq, Eq, Debug)] pub struct Variable(pub usize);

#[derive(Debug, PartialEq, Eq, Clone, Copy)] pub enum Type { I32, F32, F64 }

#[derive(Debug)] pub enum Expression {
	I32(u32),
	F32(f32),
	F64(f64),
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
	Log2(Box<Expression>),
	//Call { function: &'static str, arguments: Box<[Expression]> },
	Block { statements: Box<[Statement]>, result: Box<Expression> },
}

#[derive(Debug)] pub enum Statement {
	Value { id: Value, value: Expression },
	//Store { variable: Variable, value: Expression },
	Branch { condition: Expression, consequent: Box<[Expression]>, alternative: Box<[Expression]>, results: Box<[Value]> },
	Trace { id: Value },
}

pub fn u32(integer: u32) -> Expression { Expression::I32(integer) }

//pub fn r#use(value: &Value) -> Expression { Expression::Value(Value(value.0)) }

pub fn cast(to: Type, x: impl Into<Expression>) -> Expression { Expression::Cast(to, box_(x.into())) }
pub fn and(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::And(box_(a.into()), box_(b.into())) }
pub fn or(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Or(box_(a.into()), box_(b.into())) }
pub fn ishl_imm(a: impl Into<Expression>, b: u8) -> Expression { Expression::IShLImm(box_(a.into()), b) }
pub fn ushr_imm(a: impl Into<Expression>, b: u8) -> Expression { Expression::UShRImm(box_(a.into()), b) }
pub fn iadd(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::IAdd(box_(a.into()), box_(b.into())) }
pub fn isub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::ISub(box_(a.into()), box_(b.into())) }
fn neg(x: impl Into<Expression>) -> Expression { Expression::Neg(box_(x.into())) }
pub fn max(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Max(box_(a.into()), box_(b.into())) }
fn add(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Add(box_(a.into()), box_(b.into())) }
fn sub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Sub(box_(a.into()), box_(b.into())) }
pub fn less_or_equal(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::LessOrEqual(box_(a.into()), box_(b.into())) }
fn mul(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Mul(box_(a.into()), box_(b.into())) }
fn div(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Div(box_(a.into()), box_(b.into())) }
pub fn fpromote(x: impl Into<Expression>) -> Expression { Expression::FPromote(box_(x.into())) }
pub fn fdemote(x: impl Into<Expression>) -> Expression { Expression::FDemote(box_(x.into())) }
pub fn fcvt_to_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtToSInt(box_(x.into())) }
pub fn fcvt_from_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtFromSInt(box_(x.into())) }
pub fn fma(a: impl Into<Expression>, b: impl Into<Expression>, c: impl Into<Expression>) -> Expression { Expression::MulAdd(box_(a.into()), box_(b.into()), box_(c.into())) }
pub fn sqrt(x: impl Into<Expression>) -> Expression { Expression::Sqrt(box_(x.into())) }
//pub fn exp2(x: impl Into<Expression>) -> Expression { Expression::Exp2(box_(x.into())) }
//pub fn log2(x: impl Into<Expression>) -> Expression { Expression::Log2(box_(x.into())) }

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
impl<E:Into<Expression>> std::ops::Div<E> for &Value { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }

impl From<f32> for Expression { fn from(v: f32) -> Self { Self::F32(v) } }
impl From<f64> for Expression { fn from(v: f64) -> Self { Self::F64(v) } }
impl std::ops::Add<Expression> for f64 { type Output = Expression; fn add(self, b: Expression) -> Self::Output { add(self, b) } }
impl std::ops::Sub<Expression> for f64 { type Output = Expression; fn sub(self, b: Expression) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<Expression> for f64 { type Output = Expression; fn mul(self, b: Expression) -> Self::Output { mul(self, b) } }
impl std::ops::Div<Expression> for f64 { type Output = Expression; fn div(self, b: Expression) -> Self::Output { div(self, b) } }
impl std::ops::Add<&Value> for f64 { type Output = Expression; fn add(self, b: &Value) -> Self::Output { add(self, b) } }
impl std::ops::Sub<&Value> for f64 { type Output = Expression; fn sub(self, b: &Value) -> Self::Output { sub(self, b) } }
impl std::ops::Mul<&Value> for f64 { type Output = Expression; fn mul(self, b: &Value) -> Self::Output { mul(self, b) } }
impl std::ops::Div<&Value> for f64 { type Output = Expression; fn div(self, b: &Value) -> Self::Output { div(self, b) } }

//#[must_use] fn define(id: Value, value: impl Into<Expression>) -> Statement { Statement::Define{ id, value: value.into() } }
//#[must_use] pub fn store(variable: &Variable, value: impl Into<Expression>) -> Statement { Statement::Store{ variable: Variable(variable.0), value: value.into() } }
//#[must_use] pub fn output(index: usize, value: impl Into<Expression>) -> Statement { Statement::Output{ index, value: value.into() } }

pub struct FunctionBuilder {
	//input: usize,
	values: usize,
	//variables: usize,
}
impl FunctionBuilder {
	pub fn new(input: &[Value]) -> Self { Self { /*input: input.len(),*/ values: input.len()/*, variables: 0*/ } }
}

pub struct Block<'t> {
	statements: Vec<Statement>,
	values: &'t mut usize,
}
pub fn push(s: Statement, block: &mut Block) { block.statements.push(s) }
impl Block<'t> {
	pub fn new(f: &'t mut FunctionBuilder) -> Self { Self { /*base: f.input,*/ statements: vec![], values: &mut f.values/*variables: &mut f.variables*/ } }
	pub fn block(&mut self, build: impl Fn(&mut Block)->Expression) -> Expression {
		let mut block = Block{ statements: vec![], values: &mut self.values };
		let result = build(&mut block);
		Expression::Block { statements: block.statements.into(), result: box_(result) }
	}
	pub fn value(&mut self) -> Value {
		let id = Value(*self.values);
		*self.values += 1;
		id
	}
}
pub fn def(value: impl Into<Expression>, block: &mut Block) -> Value {
	let id = block.value();
	push(Statement::Value{id: id.clone(), value: value.into()}, block);
	id
}

impl<E:Into<Expression>> FnOnce<(E,)> for Block<'_> { type Output = Value; extern "rust-call" fn call_once(mut self, args: (E,)) -> Self::Output { self.call_mut(args) }}
impl<E:Into<Expression>> FnMut<(E,)> for Block<'_> { extern "rust-call" fn call_mut(&mut self, (value,): (E,)) -> Self::Output  { def(value, self) } }

impl From<Block<'_>> for Box<[Statement]> { fn from(b: Block) -> Self { b.statements.into() } }

impl Block<'_> {
	pub fn trace(&mut self, value: impl Into<Expression>) -> Value {
		let id = def(value, self);
		push(Statement::Trace{id: id.clone()}, self);
		id
	}
}
pub struct Function {
	pub input: usize,
	pub statements: Box<[Statement]>,
	pub output: Box<[Expression]>,
}

impl<E:Into<Expression>> std::iter::Sum<E> for Expression { fn sum<I:Iterator<Item=E>>(iter: I) -> Self { iter.into_iter().map(|e| e.into()).reduce(add).unwrap() } }
pub fn sum(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().sum() }

#[cfg(feature="num")] impl num::Sqrt for Expression { fn sqrt(self) -> Self { sqrt(self) } }
//#[cfg(feature="num")] impl num::Log for Expression { fn log2(self) -> Self { log2(self) } }

#[track_caller] pub fn idot<'t>(iter: impl IntoIterator<Item=(f64, &'t Value)>) -> Expression {
	iter.into_iter().fold(None, |sum, (c, e)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + e, None => e.into() }) }
		else if c == -1. { Some(match sum { Some(sum) => sum - e, None => -e }) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => c * e + sum, None => c * e }) }
	).unwrap()
}
pub fn dot(c: &[f64], v: &[Value]) -> Expression { idot(c.iter().copied().zip(v)) }

pub fn exp(x: impl Into<Expression>, f: &mut Block) -> Expression { //e-12 (19*,1/) (-9->-7)
	let ref x = f((1./2048.)*x.into());
	let ref x2 = f(x*x);
	let ref x3 = f(x2*x);
	let ref a = f(1.+3./28.*x2+1./1680.*x3*x);
	let ref b = f(1./2.*x+1./84.*x3);
	let sq = |x,f:&mut Block| { let ref x=f(x); x*x };
	sq(sq(sq(sq(sq(sq(sq(sq(sq(sq(sq((a+b) / (a-b),f),f),f),f),f),f),f),f),f),f),f)
	//Expression::Exp(box_(x.into()))
}
pub fn log2(x: impl Into<Expression>, _: &mut Block) -> Expression {
	/*let ref i = f(cast(I32, fdemote(x)));
	let ref m = f(cast(F32, or(and(i, u32(0x007FFFFF)), u32(1f32.to_bits()))));
	fpromote(fma(fma(fma(fma(fma(-0.034436006f32, m, 0.31821337f32), m, -1.2315303f32), m, 2.5988452f32), m, -3.3241990f32), m, 3.1157899f32) * (m - 1f32) + fcvt_from_sint(isub(ushr_imm(and(i, u32(0x7F800000)), 23), u32(127))))*/
	//ast::log2(x)
	Expression::Log2(box_(x.into()))
}

pub mod interpret;
