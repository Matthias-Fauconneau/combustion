#![feature(array_map,format_args_capture)]
use {iter::map, ast::*, ::spirv::{*, Decoration::*, BuiltIn}, ::rspirv::dr::{self as rspirv, *}};

struct Builder {
	builder: rspirv::Builder,
	gl: Word,
	values: linear_map::LinearMap<ast::Value, Word>,
}
impl std::ops::Deref for Builder { type Target=rspirv::Builder; fn deref(&self) -> &Self::Target { &self.builder } }
impl std::ops::DerefMut for Builder { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.builder } }

impl Builder {
fn expr(&mut self, e: &Expression) -> Word {
	let [f64, gl] = [self.type_float(64), self.gl];
	use Expression::*;
	match e {
		&F64(value) => self.constant_f64(f64, value),
		Value(v) => self.values[v],
		Neg(x) => { let x = self.expr(x); self.f_negate(f64, None, x).unwrap() }
		Max(a, b) => { let operands = [a,b].map(|x| Operand::IdRef(self.expr(x))); self.ext_inst(f64, None, gl, GLOp::FMax as u32, operands).unwrap() }
		Add(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_add(f64, None, a, b).unwrap() }
		Sub(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_sub(f64, None, a, b).unwrap() }
		LessOrEqual(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_ord_less_than_equal(f64, None, a, b).unwrap() }
		Mul(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_mul(f64, None, a, b).unwrap() }
		Div(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_div(f64, None, a, b).unwrap() }
		Sqrt(x) => { let x = Operand::IdRef(self.expr(x)); self.ext_inst(f64, None, gl, GLOp::Sqrt as u32, [x]).unwrap() }
		Block { statements, result } => {
			for s in &**statements { self.push(s) }
			self.expr(result)
		}
		_ => panic!("{e:?}")
	}
}
fn push(&mut self, s: &Statement) {
	let f64 = self.type_float(64);
	use Statement::*;
	match s {
		Value { id, value } => { let value = self.expr(value); assert!(!self.values.contains_key(id)); self.values.insert(id.clone(), value); },
		Select { condition, true_exprs, false_exprs, results } => {
			let condition = self.expr(condition);
			let merge = self.id();
			self.selection_merge(merge, SelectionControl::DONT_FLATTEN).unwrap();
			let true_block = self.id();
			let false_block = self.id();
			self.branch_conditional(condition, true_block, false_block, []).unwrap();

			let scope = self.values.clone();
			self.begin_block(Some(true_block)).unwrap();
			let true_values = map(&**true_exprs, |e| self.expr(e));
			self.values = scope;
			assert!(results.len() == true_values.len());
			self.branch(merge).unwrap();

			let scope = self.values.clone();
			self.begin_block(Some(false_block)).unwrap();
			let false_values = map(&**false_exprs, |e| self.expr(e));
			self.values = scope;
			assert!(results.len() == false_values.len());
			self.branch(merge).unwrap();

			self.begin_block(Some(merge)).unwrap();
			for id in &**results { assert!(!self.values.contains_key(id)); }
			let phi = map(true_values.iter().zip(&*false_values),|(&true_value, &false_value)| self.phi(f64, None, [(true_value, true_block),(false_value,false_block)]).unwrap());
			self.values.extend(results.iter().zip(&*phi).map(|(id, &phi)| (id.clone(), phi)));
		}
	}
}
}

pub fn compile(uniform_len: usize, ast: &ast::Function) -> Result<Box<[u32]>, rspirv::Error> {
	assert!(uniform_len==1);
	let mut b = rspirv::Builder::new();
	b.set_version(1, 5);
	b.capability(Capability::Shader); b.capability(Capability::VulkanMemoryModel);
	b.memory_model(AddressingModel::Logical, MemoryModel::Vulkan);
	let u32 = b.type_int(32, 0);
	let uvec3 = b.type_vector(u32, 3);
	let id = b.variable(uvec3, None, StorageClass::UniformConstant, None);
	b.decorate(id, BuiltIn, [Operand::BuiltIn(BuiltIn::GlobalInvocationId)]);
	let f64 = b.type_float(64);
	let array = b.type_runtime_array(f64);
	let gl = b.ext_inst_import("GLSL.std.450");
	let uniform = b.variable(f64, None, StorageClass::PushConstant, None);
	let input = map(0..ast.input-uniform_len, |_| {
		let v = b.variable(array, None, StorageClass::StorageBuffer, None);
		b.decorate(v, NonWritable, []);
		v
	});
	let output = map(&*ast.output,|_| {
		let v = b.variable(array, None, StorageClass::StorageBuffer, None);
		b.decorate(v, NonReadable, []);
		v
	});
	for (binding, &variable) in input.iter().chain(&*output).enumerate() {
		b.decorate(variable, DescriptorSet, [Operand::LiteralInt32(0)]);
		b.decorate(variable, Binding, [Operand::LiteralInt32(binding as u32)]);
	}
	let return_type = b.type_void();
	let function_type = b.type_function(return_type, vec![]);
	let f = b.begin_function(return_type, /*id: */None, FunctionControl::DONT_INLINE | FunctionControl::CONST, function_type)?;
	b.begin_block(None)?;
	let id = b.composite_extract(u32, None, id, [0])?;
	let uniform_constant = b.type_pointer(None, StorageClass::UniformConstant, f64);
	let input = map(&*input, |&input| {
		let input = b.access_chain(uniform_constant, None, input, [id]).unwrap();
		b.load(f64, None, input, Some(MemoryAccess::NONTEMPORAL), []).unwrap()
	});
	let values = [uniform].iter().chain(&*input).enumerate().map(|(value, &input)| (Value(value), input)).collect();
	let mut b = Builder{builder: b, gl, values};
	for s in &*ast.statements { b.push(s); }
	let uniform = b.type_pointer(None, StorageClass::Uniform, f64);
	for (expr, &output) in ast.output.iter().zip(&*output) {
		let value = b.expr(expr);
		let output = b.access_chain(uniform, None, output, [id]).unwrap();
		b.store(output, value, Some(MemoryAccess::NONTEMPORAL), []).unwrap();
	}
	b.ret()?;
	b.end_function()?;
	b.entry_point(spirv::ExecutionModel::GLCompute, f, "main", vec![]);
	Ok(::rspirv::binary::Assemble::assemble(&b.builder.module()).into())
}
