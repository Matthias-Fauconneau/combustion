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
	let [bool, f32, gl] = [self.type_bool(), self.type_float(32), self.gl];
	use Expression::*;
	match e {
		&F32(value) => self.constant_f32(f32, value),
		//&F64(value) => self.constant_f32(f32, value as f32),
		&Float(value) => self.constant_f32(f32, value as f32),
		Value(v) => self.values[v],
		Neg(x) => { let x = self.expr(x); self.f_negate(f32, None, x).unwrap() }
		Max(a, b) => { let operands = [a,b].map(|x| Operand::IdRef(self.expr(x))); self.ext_inst(f32, None, gl, GLOp::FMax as u32, operands).unwrap() }
		Add(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_add(f32, None, a, b).unwrap() }
		Sub(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_sub(f32, None, a, b).unwrap() }
		LessOrEqual(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_ord_less_than_equal(bool, None, a, b).unwrap() }
		Mul(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_mul(f32, None, a, b).unwrap() }
		Div(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_div(f32, None, a, b).unwrap() }
		Sqrt(x) => { let x = Operand::IdRef(self.expr(x)); self.ext_inst(f32, None, gl, GLOp::Sqrt as u32, [x]).unwrap() }
		Exp(x) => { let x = Operand::IdRef(self.expr(x)); self.ext_inst(f32, None, gl, GLOp::Exp as u32, [x]).unwrap() }
		Ln(x) => { let x = Operand::IdRef(self.expr(x)); self.ext_inst(f32, None, gl, GLOp::Log as u32, [x]).unwrap() }
		Block { statements, result } => {
			for s in &**statements { self.push(s) }
			self.expr(result)
		}
		_ => panic!("{e:?}")
	}
}
fn push(&mut self, s: &Statement) {
	let f32 = self.type_float(32);
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
			let phi = map(true_values.iter().zip(&*false_values),|(&true_value, &false_value)| self.phi(f32, None, [(true_value, true_block),(false_value,false_block)]).unwrap());
			self.values.extend(results.iter().zip(&*phi).map(|(id, &phi)| (id.clone(), phi)));
		}
	}
}
}

pub fn compile(constants_len: usize, ast: &ast::Function) -> Result<Box<[u32]>, rspirv::Error> {
	assert!(constants_len==1);
	let mut b = rspirv::Builder::new();
	b.set_version(1, 5);
	b.capability(Capability::Shader); b.capability(Capability::VulkanMemoryModel);
	b.memory_model(AddressingModel::Logical, MemoryModel::Vulkan);
	let void = b.type_void();
	let u32 = b.type_int(32, 0);
	let function_void = b.type_function(void, vec![]);
	let f = b.begin_function(void, /*id: */None, FunctionControl::DONT_INLINE | FunctionControl::CONST, function_void)?;
	let local_size = b.spec_constant_u32(u32, 1);
	b.decorate(local_size, SpecId, [0u32.into()]);
	let u32_1 = b.constant_u32(u32, 1);
	b.execution_mode(f, ExecutionMode::LocalSizeId, [local_size, u32_1, u32_1]);
	let gl = b.ext_inst_import("GLSL.std.450");
	let v3u = b.type_vector(u32, 3);
	let workgroup_size = b.spec_constant_composite(v3u, [local_size, u32_1, u32_1]);
	b.decorate(workgroup_size, BuiltIn, [Operand::BuiltIn(BuiltIn::WorkgroupSize)]);
	let pv3u = b.type_pointer(None, StorageClass::Input, v3u);
	let global_invocation_id_ref = b.variable(pv3u, None, StorageClass::Input, None);
	b.decorate(global_invocation_id_ref, BuiltIn, [Operand::BuiltIn(BuiltIn::GlobalInvocationId)]);
	let f32 = b.type_float(32);
	let bsf = b.type_struct([f32]);
	b.decorate(bsf, Block, []);
	b.member_decorate(bsf, 0, Offset, [0u32.into()]);
	let pcf = b.type_pointer(None, StorageClass::PushConstant, f32);
	let pcsf = b.type_pointer(None, StorageClass::PushConstant, bsf);
	let push_constants = b.variable(pcsf, None, StorageClass::PushConstant, None);
	let sbf = b.type_pointer(None, StorageClass::StorageBuffer, f32);
	let af = b.type_runtime_array(f32);
	b.decorate(af, ArrayStride, [(std::mem::size_of::<f32>() as u32).into()]);
	let bsaf = b.type_struct([af]);
	b.decorate(bsaf, Block, []);
	b.member_decorate(bsaf, 0, Offset, [0u32.into()]);
	let sbsaf = b.type_pointer(None, StorageClass::StorageBuffer, bsaf);
	let input = map(0..ast.input-constants_len, |_| {
		let v = b.variable(sbsaf, None, StorageClass::StorageBuffer, None);
		b.decorate(v, NonWritable, []);
		v
	});
	let output = map(&*ast.output,|_| {
		let v = b.variable(sbsaf, None, StorageClass::StorageBuffer, None);
		b.decorate(v, NonReadable, []);
		v
	});
	for (binding, &variable) in input.iter().chain(&*output).enumerate() {
		b.decorate(variable, DescriptorSet, [0u32.into()]);
		b.decorate(variable, Binding, [(binding as u32).into()]);
	}
	b.begin_block(None)?;
	let global_invocation_id3 = b.load(v3u, None, global_invocation_id_ref, None, [])?;
	let id = b.composite_extract(u32, None, global_invocation_id3, [0])?;
	let index0 = b.constant_u32(u32, 0);
	let push_constant_0 = b.access_chain(pcf, None, push_constants, [index0]).unwrap();
	let push_constant_0 = b.load(f32, None, push_constant_0, None, [])?;
	let input_values = map(&*input, |&input| {
		let input = b.access_chain(sbf, None, input, [index0, id]).unwrap();
		b.load(f32, None, input, Some(MemoryAccess::NONTEMPORAL), []).unwrap()
	});
	let values = [push_constant_0].iter().chain(&*input_values).enumerate().map(|(value, &input)| (Value(value), input)).collect();
	let mut b = Builder{builder: b, gl, values};
	for s in &*ast.statements { b.push(s); }
	for (expr, &output) in ast.output.iter().zip(&*output) {
		let value = b.expr(expr);
		let output = b.access_chain(sbf, None, output, [index0, id]).unwrap();
		b.store(output, value, Some(MemoryAccess::NONTEMPORAL), []).unwrap();
	}
	b.ret()?;
	b.end_function()?;
	b.entry_point(spirv::ExecutionModel::GLCompute, f, "main", [&[global_invocation_id_ref, push_constants] as &[_],&input,&output].concat());
	Ok(::rspirv::binary::Assemble::assemble(&b.builder.module()).into())
}
