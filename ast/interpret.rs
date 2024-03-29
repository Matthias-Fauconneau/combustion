use super::*;

#[derive(PartialEq, Debug, Clone)] pub enum DataValue { None, Bool(bool), /*I32(u32),*/ F32(f32), F64(f64) }
impl std::fmt::Display for DataValue { fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result { match self {
	Self::F32(v) => write!(f,"{v:e}"),
	Self::F64(v) => write!(f,"{v:e}"),
	_ => unimplemented!("{self:?}")
} } }

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

fn to_string(state: &State, expr: &Expr) -> String {
	use Expr::*;
	match expr{
		F32(x) => x.to_string(),
		F64(x) => x.to_string(),
		Value(id) => format!("{} = {}", state.debug[id.0], state.values[id.0]),
		Neg(x) => format!("-{}", to_string(state, x)),
		Add(a, b) => format!("({}) + ({})", to_string(state, a), to_string(state, b)),
		Sub(a, b) => format!("({}) - ({})", to_string(state, a), to_string(state, b)),
		Mul(a, b) => format!("({}) * ({})", to_string(state, a), to_string(state, b)),
		Div(a, b) => format!("({}) / ({})", to_string(state, a), to_string(state, b)),
		Exp(x) => format!("exp({})", to_string(state, x)),
		Ln{x,..} => format!("ln({})", to_string(state, x)),
		Min(a, b) => format!("min({}, {})", to_string(state, a), to_string(state, b)),
		Max(a, b) => format!("max({}, {})", to_string(state, a), to_string(state, b)),
		e => panic!("{:?}", e),
	}
}

impl DataValue {
	fn is_valid(&self) -> bool { match self {
		Self::F32(x) => !x.is_nan(), //x.is_finite(),
		Self::F64(x) => !x.is_nan(), //x.is_finite(),
		_ => true,
	} }
}

fn eval(state: &State, expression: &Expression) -> DataValue {
	use {Expr::*, DataValue::{Bool, /*I32,*/ F64, F32}};
	let result = match expression { Expression::Expr(expr) => match expr {
		//&I32(value) => I32(value),
		&Expr::F32(value) => F32(*value),
		&Expr::F64(value) => F64(*value),
		/*Cast(Type::I32, from) => I32(eval(state, from).f32().to_bits()),
		Cast(Type::F32, from) => F32(f32::from_bits(eval(state, from).i32())),
		And(a,b) => I32(eval(state, a).i32()&eval(state, b).i32()),
		Or(a,b) => I32(eval(state, a).i32()|eval(state, b).i32()),*/
		//IShLImm(_,_)|UShRImm(_,_)|IAdd(_,_)|ISub(_,_)|FCvtToSInt(_)|FCvtFromSInt(_) => panic!("{e:?}"),
		Value(id) => state[id.0].clone(),
		//Load(variable) => { self.variables[variable.0].unwrap() }
		Neg(x) if let F32(x) = eval(state, x) => F32(-x),
		Neg(x) if let F64(x) = eval(state, x) => F64(-x),
		Min(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(f32::min(a,b)),
		Min(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(f64::min(a,b)),
		Max(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(f32::max(a,b)),
		Max(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(f64::max(a,b)),
		Add(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a+b),
		Add(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a+b),
		Sub(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a-b),
		Sub(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a-b),
		LessOrEqual(a,b)  if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => Bool(a<=b),
		LessOrEqual(a,b)  if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => Bool(a<=b),
		Mul(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a*b),
		Mul(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a*b),
		Mul(a,b) => panic!("{:?}", [a,b].map(|x| eval(state, x))),
		//MulAdd(a,b,c) => { let [a,b,c] = [a,b,c].map(|x| eval(state, x).f64()); F64(f64::mul_add(a,b,c)) },
		Div(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a/b),
		Div(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a/b),
		FDemote(x) => F32(f64::from(eval(state, x)) as f32),
		FPromote(x) => F64(f32::from(eval(state, x)) as f64),
		Sqrt(x) if let F32(x) = eval(state, x) => F32(f32::sqrt(x)),
		Sqrt(x) if let F64(x) = eval(state, x) => F64(f64::sqrt(x)),
		Exp(x) if let F32(x) = eval(state, x) => F32(f32::exp(x)),
		//Exp(x) if let F32(x) = eval(state, x) => F32(f32::exp({assert!(x < 1e5, "{x} {}", to_string(state,expression)); x})),
		Exp(x) if let F64(x) = eval(state, x) => F64(f64::exp(x)),
		Ln{x,..} if let F32(x) = eval(state, x) => F32(f32::ln(x)),
		Ln{x,..} if let F64(x) = eval(state, x) => F64(f64::ln(x)),
		Sq(x) if let F32(x) = eval(state, x) => F32(x*x),
		Sq(x) if let F64(x) = eval(state, x) => F64(x*x),
		_ => panic!("{expr:?}"),
		},
		Expression::Block { statements, result } => {
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
	};
	assert!(result.is_valid(), "{result}: {expression:?} {}", to_string(state, expression));
	result
}

fn run(state: &mut State, statements: &[Statement]) {
	for (_program_counter, statement) in statements.iter().enumerate() {
		//eprintln!("{_program_counter}: {statement:?}");
		use Statement::*;
		match statement {
			Value { id, value: expression } => {
				let result= eval(state, expression);
				assert!(state[id.0] == DataValue::None);
				assert!(result.is_valid() /*&& f32::abs(x)<1e28*/, "{} = {result}: {expression:?} {}", state.debug[id.0], to_string(state, expression));
				//println!("{} = {} = {result}", state.debug[id.0], to_string(state, expression));
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

pub fn call<T:Copy+Into<DataValue>+From<DataValue>>(f: &Function, input: &[T], output: &mut [T]) {
	assert!(input.len() == f.input.len());
	let mut state = State{
		values: input.iter().map(|&v| v.into()).chain((0..f.values.len()).map(|_| DataValue::None)).collect(),
		debug: &f.values
	};
	run(&mut state, &f.statements);
	for (slot, e) in output.iter_mut().zip(&*f.output) { *slot = eval(&state, e).into(); }
}

impl<T:Copy+Into<DataValue>+From<DataValue>> FnOnce<(&[T], &mut [T])> for Function {
	type Output = (); extern "rust-call" fn call_once(mut self, args: (&[T], &mut [T])) -> Self::Output { self.call_mut(args) }
}
impl<T:Copy+Into<DataValue>+From<DataValue>> FnMut<(&[T], &mut [T])> for Function { extern "rust-call" fn call_mut(&mut self, args: (&[T], &mut [T])) -> Self::Output { self.call(args) } }
impl<T:Copy+Into<DataValue>+From<DataValue>> Fn<(&[T], &mut [T])> for Function { extern "rust-call" fn call(&self, (input, output): (&[T], &mut [T])) -> Self::Output { call(&self, input, output); } }
