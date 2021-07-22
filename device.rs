use {anyhow::Result, iter::map};
type Output<T> = Result<Box<[Box<[T]>]>>;
pub trait Convert { fn convert(f: ast::Function) -> ast::Function; }
impl Convert for f64 { fn convert(f: ast::Function) -> ast::Function { f } }
impl Convert for f32 { fn convert(mut f: ast::Function) -> ast::Function {
	use {ast::*, Expr::*};
	f.input = vec![Type::F32; f.input.len()].into();
	fn demote(e: &mut Expression) { if let Expression::Expr(F64(ref x)) = e { *e = f32(f64::from(*x) as _).unwrap().into() } else { e.visit_mut(demote) } }
	for s in &mut *f.statements { use Statement::*; match s {
		Value{value,..} => demote(value),
		Select { condition, true_exprs, false_exprs, .. } => {
			demote(condition);
			for e in &mut **true_exprs { demote(e); }
			for e in &mut **false_exprs { demote(e); }
		}
	}}
	for e in &mut *f.output{ demote(e); }
	f
}}

#[cfg(not(feature="vpu"))] mod device {
	use {std::default::default, iter::{list, map}, ast::*, super::*};
	#[cfg(not(feature="interpret"))] pub trait Convert = super::Convert;
	#[cfg(feature="interpret")] pub trait Convert = super::Convert+Into<interpret::DataValue>+From<interpret::DataValue>;
	pub fn assemble<T:'t+Copy+Default+Convert>(function: Function) -> impl 't+Fn(&[T], &[&[T]]) -> Output<T> {
		let input_len = function.input.len();
		let output_len = function.output.len();
		let function = T::convert(function);
		#[cfg(feature="ir")] let function = ir::assemble(ir::compile(&function));
		move |constants:&[T], inputs:&[&[T]]| {
			assert!(constants.len() == 1 && constants.len()+inputs.len() == input_len);
			let states_len = inputs[0].len();
			let mut outputs = map(0..output_len, |_| vec![default(); states_len].into_boxed_slice());
			let time = std::time::Instant::now();
			let mut input = list(constants.iter().copied().chain(inputs.iter().map(|_| default())));
			let mut output = vec![default(); output_len].into_boxed_slice();
			for state_id in 0..states_len {
				for (input, array) in input[constants.len()..].iter_mut().zip(inputs) { *input = array[state_id]; }
				#[cfg(any(feature="interpret",feature="ir"))] function(&input, &mut output);
				#[cfg(not(any(feature="interpret",feature="ir")))] compile_error!("Either feature \"interpret\" or \"ir\" must be enabled for this crate.");
				for (array, output) in outputs.iter_mut().zip(&*output) { array[state_id] = *output; }
			}
			let time= time.elapsed().as_secs_f64();
			if states_len > 1 { println!("{:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e9, (states_len as f64)/1e6/time); }
			Ok(outputs)
		}
	}
}
#[cfg(feature="vpu")] mod device {
use {std::default::default, std::borrow::Borrow, iter::{list, map}, vulkan::*, super::*};
pub struct Function<D: Borrow<Device>, T:Plain> {
	device: D,
	input_len: usize,
	output_len: usize,
	function: Box<[u32]>,
	_marker: std::marker::PhantomData<T>,
}
pub fn call<D: Borrow<Device>, T:Plain+Default>(Function{input_len, output_len, device, function,..}: &Function<D,T>, constants: &[T], input: &[&[T]]) -> Output<T> {
	assert!(std::mem::size_of::<T>() == 4 && constants.len() == 1 && constants.len()+input.len() == *input_len);
	let states_len = input[0].len();
	let local_size = std::cmp::min(512, states_len as u32);
	let device = device.borrow();
	let input = map(&*input, |array| Buffer::new(device, array).unwrap());
	let output = map(0..*output_len, |_| Buffer::new(device, &vec![default(); states_len]).unwrap());
	let buffers = list(input.iter().chain(&*output));
	let pipeline = device.pipeline(&function, local_size, as_u8(&constants), &buffers)?;
	device.bind(pipeline.descriptor_set, &buffers)?;
	let command_buffer = device.command_buffer(&pipeline, as_u8(&constants), (states_len as u32)/local_size)?;
	for _ in 0..1 { // First iteration loads from host visible to device local memory
		let time = device.submit_and_wait(command_buffer)?;
		if states_len > 1 { println!("{local_size}: {:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e9, (states_len as f32)/1e6/time); }
	}
	Ok(map(&*output, |array| (*array.map(device).unwrap()).into()))
}

/*#[throws] pub fn assemble<'t>(device: &'t Device, function: &ast::Function) -> Function<&'t Device> {
	Function{device, input_len: function.input, output_len: function.output.len(), function: spirv::compile(1, function)?}
}*/
pub fn assemble<T:Plain+Convert>(function: ast::Function) -> Function<Device, T> {
	let function = T::convert(function);
	Function{device: vulkan::Device::new().unwrap(), input_len: function.input.len(), output_len: function.output.len(), function: spirv::compile(1, &function).unwrap(), _marker: default()}
}

impl<D: Borrow<Device>, T:Plain+Default> FnOnce<(&[T], &[&[T]])> for Function<D,T> { type Output = self::Output<T>; extern "rust-call" fn call_once(mut self, args: (&[T], &[&[T]])) -> Self::Output { self.call_mut(args) } }
impl<D: Borrow<Device>, T:Plain+Default> FnMut<(&[T], &[&[T]])> for Function<D,T> { extern "rust-call" fn call_mut(&mut self, args: (&[T], &[&[T]])) -> Self::Output { self.call(args) } }
impl<D: Borrow<Device>, T:Plain+Default> Fn<(&[T], &[&[T]])> for Function<D,T> { extern "rust-call" fn call(&self, (constants, input): (&[T], &[&[T]])) -> Self::Output { call(&self, constants, input) } }
}
pub use {device::*, ast::let_};
pub fn all_same<T:PartialEq+Copy>(array:&[T], times: usize) -> T { assert!(array.len() == times); for &v in array { if v != array[0] {panic!()} } array[0] }
pub fn with_repetitive_input<T:Copy+PartialEq>(f: impl Fn(&[T],&[&[T]])->Output<T>, times: usize) -> impl Fn(&[T],&[T])->Result<Box<[T]>> {
	move |constants, inputs| Ok(map(&*f(&map(constants, |x| *x as _), &map(&*map(inputs, |x| vec![*x as _; times]), |x| &**x))?, |y| all_same(y, times)))
}
