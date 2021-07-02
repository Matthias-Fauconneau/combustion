#![feature(unboxed_closures,fn_traits,in_band_lifetimes,associated_type_bounds,format_args_capture)]
fn box_<T>(t: T) -> Box<T> { Box::new(t) }
#[macro_export] macro_rules! let_ { { $p:pat = $e:expr => $b:block } => { if let $p = $e { $b } else { unreachable!() } } }

#[derive(PartialEq, Eq, Debug, Clone)] pub struct Value(pub usize);
//#[derive(PartialEq, Eq, Debug)] pub struct Variable(pub usize);

#[derive(Debug)] pub enum Type { I32, F32 }

#[derive(Debug)] pub enum Expression {
	Float(f64),
	Integer(u32),
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
	FCvtToSInt(Box<Expression>),
	FCvtFromSInt(Box<Expression>),
	Sqrt(Box<Expression>),
	//Call { function: &'static str, arguments: Box<[Expression]> },
	Block { statements: Box<[Statement]>, result: Box<Expression> },
}

#[derive(Debug)] pub enum Statement {
	Value { id: Value, value: Expression },
	//Store { variable: Variable, value: Expression },
	Branch { condition: Expression, consequent: Box<[Expression]>, alternative: Box<[Expression]>, results: Box<[Value]> },
}

pub fn u32(integer: u32) -> Expression { Expression::Integer(integer) }

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
pub fn fcvt_to_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtToSInt(box_(x.into())) }
pub fn fcvt_from_sint(x: impl Into<Expression>) -> Expression { Expression::FCvtFromSInt(box_(x.into())) }
pub fn sqrt(x: impl Into<Expression>) -> Expression { Expression::Sqrt(box_(x.into())) }
pub fn fma(a: impl Into<Expression>, b: impl Into<Expression>, c: impl Into<Expression>) -> Expression { Expression::MulAdd(box_(a.into()), box_(b.into()), box_(c.into())) }

impl std::ops::Neg for Expression { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }

impl<E:Into<Expression>> std::ops::Add<E> for Expression { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for Expression { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for Expression { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for Expression { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }

impl From<&Value> for Expression { fn from(value: &Value) -> Expression { Expression::Value(Value(value.0)) } }
impl std::ops::Neg for &Value { type Output = Expression; fn neg(self) -> Self::Output { neg(self) } }
impl<E:Into<Expression>> std::ops::Add<E> for &Value { type Output = Expression; fn add(self, b: E) -> Self::Output { add(self, b) } }
impl<E:Into<Expression>> std::ops::Sub<E> for &Value { type Output = Expression; fn sub(self, b: E) -> Self::Output { sub(self, b) } }
impl<E:Into<Expression>> std::ops::Mul<E> for &Value { type Output = Expression; fn mul(self, b: E) -> Self::Output { mul(self, b) } }
impl<E:Into<Expression>> std::ops::Div<E> for &Value { type Output = Expression; fn div(self, b: E) -> Self::Output { div(self, b) } }

impl From<f64> for Expression { fn from(v: f64) -> Self { Self::Float(v) } }
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

//#[derive(derive_more::Deref,derive_more::DerefMut)]
pub struct Block<'t> {
	//base: usize,
	//#[deref]#[deref_mut]
	statements: Vec<Statement>,
	values: &'t mut usize,
	//variables: &'t mut usize,
}
pub fn push(s: Statement, block: &mut Block) { block.statements.push(s) }
impl Block<'t> {
	pub fn new(f: &'t mut FunctionBuilder) -> Self { Self { /*base: f.input,*/ statements: vec![], values: &mut f.values/*variables: &mut f.variables*/ } }
	pub fn block(&mut self, build: impl Fn(&mut Block)->Expression) -> Expression {
		let mut block = Block{ /*base: self.base+self.statements.len(),*/ statements: vec![], values: &mut self.values/*variables: self.variables*/ };
		let result = build(&mut block);
		//self.base += block.statements.len();
		Expression::Block { statements: block.statements.into(), result: box_(result) }
	}
	pub fn value(&mut self) -> Value {
		let id = Value(*self.values);
		*self.values += 1;
		id
	}
	/*pub fn def(&mut self, value: Expression) -> Value {
		let id = Value(self.base+self.statements.len());
		self.statements.push(define(id.clone(), value));
		id
	}*/
	/*pub fn decl(&mut self) -> Variable {
		let id = Variable(*self.variables);
		*self.variables += 1;
		id
	}
	pub fn load(&mut self, variable: Variable) -> Value { self.def(Expression::Load(variable)) }*/
}
pub fn def(value: impl Into<Expression>, block: &mut Block) -> Value {
	let id = block.value();
	push(Statement::Value{id: id.clone(), value: value.into()}, block);
	id
}

impl<E:Into<Expression>> FnOnce<(E,)> for Block<'_> { type Output = Value; extern "rust-call" fn call_once(mut self, args: (E,)) -> Self::Output { self.call_mut(args) }}
impl<E:Into<Expression>> FnMut<(E,)> for Block<'_> { extern "rust-call" fn call_mut(&mut self, (value,): (E,)) -> Self::Output  { def(value, self) } }

impl From<Block<'_>> for Box<[Statement]> { fn from(b: Block) -> Self { b.statements.into() } }

pub struct Function {
	pub input: usize,
	//pub values: usize,
	//pub variables: usize,
	pub statements: Box<[Statement]>,
	pub output: Box<[Expression]>,
}
/*impl Function {
	pub fn new(output: Box<[Expression]>, statements: Box<[Statement]>, FunctionBuilder{input, ../*values, variables*/}: FunctionBuilder) -> Self { Function{input, output, /*values,*/ /*variables,*/ statements} }
}*/

impl FnOnce<(&[f64], &mut [f64])> for Function {
	type Output = ();
	extern "rust-call" fn call_once(mut self, args: (&[f64], &mut [f64])) -> Self::Output { self.call_mut(args) }
}
impl FnMut<(&[f64], &mut [f64])> for Function {
	extern "rust-call" fn call_mut(&mut self, args: (&[f64], &mut [f64])) -> Self::Output { self.call(args) }
}
impl Fn<(&[f64], &mut [f64])> for Function {
	extern "rust-call" fn call(&self, (input, output): (&[f64], &mut [f64])) -> Self::Output {
		/*struct State {
			//input: &'t [f64],
			//output: &'t mut [f64],
			values: linear_map::LinearMap<Value, f64>,
			//variables: Box<[Option<f64>]>,
		}
		impl State {*/
		/*fn debug(&mut self, expression: &Expression) -> String {
			use Expression::*;
			match expression {
				Value(id) => format!("{}", if id.0 < self.input.len() { self.input[id.0] } else { self.values[&id] }),
				Add(a, b) => format!("{} + {}", eval(state, a), eval(state, b)),
				Mul(a, b) => format!("{} * {}", eval(state, a), eval(state, b)),
				Call { function, arguments } => match *function {
					"exp2" => format!("exp2({})",self.debug(&arguments[0])),
					function => panic!("{}", function)
				}
				e => panic!("{:?}", e),
			}
		}*/
		type State = Vec<f64>;
		fn eval(state: &mut State, expression: &Expression) -> f64 {
			use Expression::*;
			let result = match expression {
				&Float(value) => value,
				&Integer(_)|Cast(_,_)|And(_,_)|Or(_,_)|IShLImm(_,_)|UShRImm(_,_)|IAdd(_,_)|ISub(_,_)|FCvtToSInt(_)|FCvtFromSInt(_) => panic!(),
				Value(id) => { assert!(state[id.0].is_finite()); state[id.0] }, //if id.0 < self.input.len() { self.input[id.0] } else { self.values[&id] },
				//Load(variable) => { self.variables[variable.0].unwrap() }
				Neg(x) => -eval(state, x),
				Max(a,b) => f64::max(eval(state, a), eval(state, b)),
				Add(a, b) => eval(state, a) + eval(state, b),
				Sub(a, b) => eval(state, a) - eval(state, b),
				LessOrEqual(a, b) => if eval(state, a) <= eval(state, b) { 1. } else { 0. },
				Mul(a, b) => eval(state, a) * eval(state, b),
				MulAdd(a, b, c) => eval(state, a).mul_add(eval(state, b), eval(state, c)),
				Div(a, b) => eval(state, a) / eval(state, b),
				Sqrt(x) => f64::sqrt(eval(state, x)),
				Block { statements, result } => {
					run(state, statements);
					let result = eval(state, result);
					for statement in statements.iter() { if let Statement::Value{ id, .. } = statement { state[id.0] = f64::NAN; } else { unreachable!() } }
					result
				},
				//e => panic!("{:?}", e),
			};
			assert!(result.is_finite());//, "{result}: {expression:?} {}", self.debug(expression));
			result
		}
		fn run(state: &mut State, statements: &[Statement]) {
			for statement in statements {
				use Statement::*;
				match statement {
					Value { id, value } => {
						let value = eval(state, value);
						assert!(state.len() < id.0);
						state.resize(id.0+1, f64::NAN);
						state[id.0] = value;
					},
					/*Store { variable, value } => {
						let value = eval(state, value);
						assert!(self.variables[variable.0].replace(value).is_none());
					}*/
					Branch { condition, consequent, alternative, results } => {
						let values = if eval(state, condition) != 0. { consequent } else { alternative };
						for (id, value) in results.iter().zip(&**values) {
							let value = eval(state, value);
							assert!(state.len() < id.0);
							state.resize(id.0+1, f64::NAN);
							state[id.0] = value;
						}
					},
				}
			}
		}
		assert!(input.len() == self.input);
		let ref mut state = input.to_vec();
		run(state, &self.statements);
		for (slot, e) in output.iter_mut().zip(&*self.output) { *slot = eval(state, e); }
	}
}

pub fn wrap(f: Function) -> impl Fn(&[f64]) -> Box<[f64]> { move |input| { let mut output = vec![0.; f.output.len()].into_boxed_slice(); f(input, &mut output); output } }

impl<E:Into<Expression>> std::iter::Sum<E> for Expression { fn sum<I:Iterator<Item=E>>(iter: I) -> Self { iter.into_iter().map(|e| e.into()).reduce(add).unwrap() } }
pub fn sum(iter: impl IntoIterator<Item:Into<Expression>>) -> Expression { iter.into_iter().sum() }

/*pub fn max(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Call{ function: "max", arguments: box_([a.into(), b.into()]) } }
pub fn sqrt(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "sqrt", arguments: box_([x.into()]) } }
pub fn band(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Call{ function: "band", arguments: box_([a.into(), b.into()]) } }
pub fn bor(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Call{ function: "bor", arguments: box_([a.into(), b.into()]) } }
pub fn fcvt_to_sint(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "fcvt_to_sint", arguments: box_([x.into()]) } } // fcvt_to_sint I32
pub fn fcvt_from_sint(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "fcvt_from_sint", arguments: box_([x.into()]) } } // fcvt_from_sint F32
pub fn iadd(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Call{ function: "iadd", arguments: box_([a.into(), b.into()]) } }
pub fn isub(a: impl Into<Expression>, b: impl Into<Expression>) -> Expression { Expression::Call{ function: "isub", arguments: box_([a.into(), b.into()]) } }*/
//pub fn exp2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "exp2", arguments: box_([x.into()]) } }
//pub fn log2(x: impl Into<Expression>) -> Expression { Expression::Call{ function: "log2", arguments: box_([x.into()]) } }

#[cfg(feature="num")] impl num::Sqrt for Expression { fn sqrt(self) -> Self { sqrt(self) } }
#[cfg(feature="num")] impl num::Log for Expression { fn log2(self) -> Self { log2(self) } }

#[track_caller] pub fn idot<'t>(iter: impl IntoIterator<Item=(f64, &'t Value)>) -> Expression {
	iter.into_iter().fold(None, |sum, (c, e)|
		if c == 0. { sum }
		else if c == 1. { Some(match sum { Some(sum) => sum + e, None => e.into() }) }
		else if c == -1. { Some(match sum { Some(sum) => sum - e, None => -e }) } // fixme: reorder -a+b -> b-a to avoid neg
		else { Some(match sum { Some(sum) => c * e + sum, None => c * e }) }
	).unwrap()
}
pub fn dot(c: &[f64], v: &[Value]) -> Expression { idot(c.iter().copied().zip(v)) }
