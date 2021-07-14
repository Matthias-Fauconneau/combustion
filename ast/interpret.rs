use super::*;

#[derive(PartialEq, Debug, Clone)] enum DataValue { None, Bool(bool), I32(u32), F32(f32), F64(f64) }
impl std::fmt::Display for DataValue { fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result { match self { Self::F32(v) => write!(f,"{v:e}"), _ => unimplemented!() } } }

impl From<&DataValue> for f32 { fn from(v: &DataValue) -> Self { if let DataValue::F32(v) = v { *v } else { f64::from(v) as _ } } }
impl From<&DataValue> for f64 { fn from(v: &DataValue) -> Self { if let DataValue::F64(v) = v { *v } else { panic!("{v:?}"); } } }
impl From<DataValue> for bool { fn from(v: DataValue) -> Self { if let DataValue::Bool(v) = v { v } else { panic!("{v:?}") } } }
impl From<DataValue> for f32 { fn from(v: DataValue) -> Self { (&v).into() } }
impl From<DataValue> for f64 { fn from(v: DataValue) -> Self { (&v).into() } }
impl From<f32> for DataValue { fn from(v: f32) -> Self { DataValue::F32(v) } }
impl From<f64> for DataValue { fn from(v: f64) -> Self { DataValue::F64(v) } }

#[derive(Clone)] struct State<'t> {
	values: Box<[DataValue]>,
	debug: &'t [String],
}
impl std::ops::Deref for State<'_> { type Target=[DataValue]; fn deref(&self) -> &Self::Target { &self.values} }
impl std::ops::DerefMut for State<'_> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.values} }

fn to_string(state: &State, expression: &Expression) -> String {
	use Expression::*;
	match expression {
		Float(x) => format!("{:e}", x),
		Value(id) => format!("{} = {}", state.debug[id.0], state.values[id.0]),
		Add(a, b) => format!("{} + {}", to_string(state, a), to_string(state, b)),
		Sub(a, b) => format!("{} - {}", to_string(state, a), to_string(state, b)),
		Mul(a, b) => format!("({}) * ({})", to_string(state, a), to_string(state, b)),
		Div(a, b) => format!("({}) / ({})", to_string(state, a), to_string(state, b)),
		Exp(x) => format!("exp({})", to_string(state, x)),
		e => panic!("{:?}", e),
	}
}

impl DataValue {
	fn is_valid(&self) -> bool { match self {
		Self::F32(x) => x.is_finite(),
		Self::F64(x) => x.is_finite(),
		_ => true,
	} }
}

fn eval(state: &State, expression: &Expression) -> DataValue {
	use {Expression::*, DataValue::{Bool, I32, F64, F32}};
	let result = match expression {
		&Expression::I32(value) => I32(value),
		&Expression::F32(value) => F32(value),
		&Expression::F64(value) => F64(value),
		&Expression::Float(value) => {assert!((value as f32).is_finite()); (value as float).into()},
		/*Cast(Type::I32, from) => I32(eval(state, from).f32().to_bits()),
		Cast(Type::F32, from) => F32(f32::from_bits(eval(state, from).i32())),
		And(a,b) => I32(eval(state, a).i32()&eval(state, b).i32()),
		Or(a,b) => I32(eval(state, a).i32()|eval(state, b).i32()),*/
		//IShLImm(_,_)|UShRImm(_,_)|IAdd(_,_)|ISub(_,_)|FPromote(_)|FCvtToSInt(_)|FCvtFromSInt(_) => panic!("{e:?}"),
		Value(id) => state[id.0].clone(),
		//Load(variable) => { self.variables[variable.0].unwrap() }
		Neg(x) if let F32(x) = eval(state, x) => F32(-x),
		Neg(x) if let F64(x) = eval(state, x) => F64(-x),
		Max(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(f32::max(a,b)),
		Max(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(f64::max(a,b)),
		Add(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a+b),
		Add(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a+b),
		Sub(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a-b),
		Sub(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a-b),
		LessOrEqual(a,b)  if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => Bool(a<=b),
		LessOrEqual(a,b)  if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => Bool(a<=b),
		Mul(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a*b),
		Mul(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a*b),
		//MulAdd(a,b,c) => { let [a,b,c] = [a,b,c].map(|x| eval(state, x).f64()); F64(f64::mul_add(a,b,c)) },
		Div(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a/b),
		Div(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a/b),
		//FDemote(x) => F32(eval(state, x).f64() as f32),
		//FPromote(x) => F64(eval(state, x).f32() as f64),
		Sqrt(x) if let F32(x) = eval(state, x) => F32(f32::sqrt(x)),
		Sqrt(x) if let F64(x) = eval(state, x) => F64(f64::sqrt(x)),
		//Exp(x) if let F32(x) = eval(state, x) => F32(f32::exp(x)),
		Exp(x) if let F32(x) = eval(state, x) => F32(f32::exp({assert!(x < 1e5, "{x} {}", to_string(state,expression)); x})),
		Exp(x) if let F64(x) = eval(state, x) => F64(f64::exp(x)),
		Ln{x,..} if let F32(x) = eval(state, x) => F32(f32::ln(x)),
		Ln{x,..} if let F64(x) = eval(state, x) => F64(f64::ln(x)),
		Block { statements, result } => {
			let ref mut state = (*state).clone();
			run(state, statements);
			let result = eval(state, result);
			/*for statement in statements.iter() {
				match statement {
					Statement::Value{ id, .. } => { state[id.0] = DataValue::None; }
					//Statement::Trace{..} => {},
					_ => { unreachable!() }
				}
			}*/
			result
		},
		expression => panic!("{expression:?}"),
	};
	assert!(result.is_valid(), "{result}: {expression:?} {}", to_string(state, expression));
	result
}

fn run(state: &mut State, statements: &[Statement]) {
	for statement in statements {
		use Statement::*;
		match statement {
			Value { id, value: expression } => {
				let result= eval(state, expression);
				assert!(state[id.0] == DataValue::None);
				assert!(result.is_valid() /*&& f32::abs(x)<1e28*/, "{} = {result}: {expression:?} {}", state.debug[id.0], to_string(state, expression));
				state[id.0] = result;
			},
			Select { condition, true_exprs, false_exprs, results } => {
				let expressions = if eval(state, condition).into() { true_exprs } else { false_exprs };
				for (id, expression) in results.iter().zip(&**expressions) {
					let result = eval(state, expression);
					assert!(result.is_valid() /*&& f32::abs(x)<1e28*/, "{} = {result}: {}", state.debug[id.0], to_string(state, expression));
					assert!(state[id.0] == DataValue::None);
					state[id.0] = result;
				}
			},
			//Display(id) => if false { println!("{} = {}", state.debug[id.0], state.values[id.0]); },
		}
	}
}

pub fn call(f: &Function, input: &[float], output: &mut [float]) {
	assert!(input.len() == f.input);
	let mut state = State{
		values: input.iter().map(|&v| v.into()).chain((0..f.values.len()).map(|_| DataValue::None)).collect(),
		debug: &f.values
	};
	run(&mut state, &f.statements);
	for (slot, e) in output.iter_mut().zip(&*f.output) { *slot = eval(&state, e).into(); }
}

impl FnOnce<(&[float], &mut [float])> for Function {
	type Output = (); extern "rust-call" fn call_once(mut self, args: (&[float], &mut [float])) -> Self::Output { self.call_mut(args) }
}
impl FnMut<(&[float], &mut [float])> for Function { extern "rust-call" fn call_mut(&mut self, args: (&[float], &mut [float])) -> Self::Output { self.call(args) } }
impl Fn<(&[float], &mut [float])> for Function { extern "rust-call" fn call(&self, (input, output): (&[float], &mut [float])) -> Self::Output { call(&self, input, output); } }
