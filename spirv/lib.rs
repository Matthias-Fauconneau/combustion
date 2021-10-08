#![allow(incomplete_features)]#![feature(format_args_capture,default_free_fn,if_let_guard)]
use {std::default::default, iter::{list, map}, ast::*, ::spirv::{*, Decoration::*, BuiltIn}, ::rspirv::dr::{self as rspirv, *, Operand::IdRef}};

type Type= Word;
type Value = Word;
//use num_traits::cast::ToPrimitive;
type R32 = ordered_float::NotNan<f32>;
type R64 = ordered_float::NotNan<f64>;

fn stype(b: &mut rspirv::Builder, rtype: &ast::Type) -> Type { b.type_float(rtype.bits()) }
fn wrap_f32_op(b: &mut rspirv::Builder, rtype: &ast::Type, x: Value, y: impl Fn(&mut rspirv::Builder, Value)->Value) -> Value {
	use ast::Type::*;
	match rtype {
		F32 => y(b, x),
		F64 => {
			let f32 = b.type_float(32);
			let x = b.f_convert(f32, None, x).unwrap();
			let y = y(b, x);
			let f64 = stype(b, rtype);
			b.f_convert(f64, None, y).unwrap()
		}
	}
}

struct Builder<'t> {
	builder: rspirv::Builder,
	values: Box<[Option<(ast::Type, Value)>]>,
	gl: Word,
	//debug: Word,
	//constants_u32: std::collections::HashMap<u32, Value>,
	constants_f32: std::collections::HashMap<R32, Value>,
	constants_f64: std::collections::HashMap<R64, Value>,
	//expressions: std::collections::HashSet<Expr/*(Op, LeafValue, Option<LeafValue>)*/>,
	names: &'t [String],
}
impl std::ops::Deref for Builder<'_> { type Target=rspirv::Builder; fn deref(&self) -> &Self::Target { &self.builder } }
impl std::ops::DerefMut for Builder<'_> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.builder } }

impl Builder<'_> {
//fn u32(&mut self, value: u32) -> Value { *self.constants_u32.entry(value).or_insert_with(|| self.builder.constant_u32(u32, value)) }
fn f32(&mut self, value: f32) -> Value {
	let f32 = self.type_float(32);
	//assert!(value != 0. && value != -0.);
	let value = R32::new(if value == -0. {0. } else { value }).unwrap();
	*self.constants_f32.entry(value).or_insert_with(|| self.builder.constant_f32(f32, *value))
}
fn f64(&mut self, value: f64) -> Value {
	let f64 = self.type_float(64);
	//assert!(value != 0. && value != -0.);
	let value = R64::new(if value == -0. {0. } else { value }).unwrap();
	*self.constants_f64.entry(value).or_insert_with(|| self.builder.constant_f64(f64, *value))
}
fn rtype(&self, e: &Expression) -> ast::Type { e.rtype(&|value| self.values[value.0].unwrap().0) }
fn stype(&mut self, rtype: &ast::Type) -> Type { stype(self, rtype) }

fn expr(&mut self, expr: &Expression) -> Value {
	let rtype = self.rtype(expr);
	let [bool, f32, gl, stype] = [self.type_bool(), self.type_float(32), self.gl, self.stype(&rtype)];
	match expr {
		Expression::Expr(e) => {
			use Expr::*;
			//if !(e.is_leaf() || (!expr.has_block() && self.expressions.insert(e.clone()))) { panic!("{}", e.to_string(self.names)) }
			match e {
				&F32(value) => self.f32(*value),
				&F64(value) => self.f64(*value),
				Value(v) => self.values[v.0].unwrap().1,
				Neg(x) => { let x = self.expr(x); self.f_negate(stype, None, x).unwrap() }
				FPromote(x) => { let x = self.expr(x); self.f_convert(stype, None, x).unwrap() }
				FDemote(x) => { let x = self.expr(x); self.f_convert(stype, None, x).unwrap() }
				Min(a, b) => { let [a,b] = [a,b].map(|x| IdRef(self.expr(x))); self.ext_inst(stype, None, gl, GLOp::FMin as u32, [a,b]).unwrap() }
				Max(a, b) => { let [a,b] = [a,b].map(|x| IdRef(self.expr(x))); self.ext_inst(stype, None, gl, GLOp::FMax as u32, [a,b]).unwrap() }
				Add(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_add(stype, None, a, b).unwrap() }
				Sub(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_sub(stype, None, a, b).unwrap() }
				LessOrEqual(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_ord_less_than_equal(bool, None, a, b).unwrap() }
				Mul(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_mul(stype, None, a, b).unwrap() }
				Div(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); self.f_div(stype, None, a, b).unwrap() }
				Sqrt(x) => { let x = self.expr(x); self.ext_inst(stype, None, gl, GLOp::Sqrt as u32, [IdRef(x)]).unwrap() }
				Exp(x) => { let x = self.expr(x); wrap_f32_op(self, &rtype, x, |b,x| b.ext_inst(f32, None, gl, GLOp::Exp as u32, [IdRef(x)]).unwrap()) },
				Ln{x,..} => { let x = self.expr(x); wrap_f32_op(self, &rtype, x, |b,x| b.ext_inst(f32, None, gl, GLOp::Log as u32, [IdRef(x)]).unwrap()) },
				Sq(x) => { let x = self.expr(x); self.f_mul(stype, None, x, x).unwrap() }
			}
		}
		Expression::Block { statements, result } => {
			for s in &**statements { self.extend(s) }
			self.expr(result)
		}
	}
}

fn extend(&mut self, s: &Statement) {
	use Statement::*;
	match s {
		Value { id, value } => {
			let result = self.expr(value);
			self.names;/*let [void, debug] = [self.type_void(), self.debug];
			let format = self.builder.string(format!("{} %f", self.names[id.0]));
			self.ext_inst(void, None, debug, DebugPrintFOp::DebugPrintf as u32, [IdRef(format), IdRef(result)]).unwrap();*/
			let rtype = self.rtype(value);
			assert!(self.values[id.0].replace((rtype,result)).is_none());
		},
		Select { condition, true_exprs, false_exprs, results } => {
			let condition = self.expr(condition);
			let merge = self.id();
			self.selection_merge(merge, SelectionControl::DONT_FLATTEN).unwrap();
			let true_block = self.id();
			let false_block = self.id();
			self.branch_conditional(condition, true_block, false_block, []).unwrap();

			let scope = self.values.clone();
			self.begin_block(Some(true_block)).unwrap();
			let true_values = map(&**true_exprs, |e| (self.expr(e), self.rtype(e)));
			self.values = scope;
			assert!(results.len() == true_values.len());
			self.branch(merge).unwrap();

			let scope = self.values.clone();
			self.begin_block(Some(false_block)).unwrap();
			let false_values = map(&**false_exprs, |e| (self.expr(e), self.rtype(e)));
			self.values = scope;
			assert!(results.len() == false_values.len());
			self.branch(merge).unwrap();

			self.begin_block(Some(merge)).unwrap();
			for (id, (&(true_value,true_type), &(false_value,false_type))) in results.iter().zip(true_values.iter().zip(&*false_values)) {
				assert!(true_type == false_type); let rtype = true_type;
				let stype = self.stype(&rtype);
				assert!(self.values[id.0].replace((rtype, self.builder.phi(stype, None, [(true_value, true_block),(false_value,false_block)]).unwrap())).is_none());
			}
		}
	}
}
}

pub fn compile(constants_len: usize, ast: &ast::Function) -> Result<Box<[u32]>, rspirv::Error> {
	let mut b = rspirv::Builder::new();
	b.set_version(1, 5);
	b.capability(Capability::Shader); b.capability(Capability::VulkanMemoryModel);
	if ast.input.contains(&ast::Type::F64) { b.capability(Capability::Float64); }
	b.memory_model(AddressingModel::Logical, MemoryModel::Vulkan);
	b.extension("SPV_KHR_non_semantic_info");
	let void = b.type_void();
	let u32 = b.type_int(32, 0);
	let function_void = b.type_function(void, vec![]);
	let f = b.begin_function(void, /*id: */None, FunctionControl::DONT_INLINE | FunctionControl::CONST, function_void)?;
	let local_size = b.spec_constant_u32(u32, 1);
	b.decorate(local_size, SpecId, [0u32.into()]);
	let u32_1 = b.constant_u32(u32, 1);
	//b.execution_mode_id(f, ExecutionMode::LocalSizeId, [local_size, u32_1, u32_1]);
	b.execution_mode(f, ExecutionMode::LocalSize, [1, 1, 1]);
	let gl = b.ext_inst_import("GLSL.std.450");
	//let debug = b.ext_inst_import("NonSemantic.DebugPrintf");
	let v3u = b.type_vector(u32, 3);
	let workgroup_size = b.spec_constant_composite(v3u, [local_size, u32_1, u32_1]);
	b.decorate(workgroup_size, BuiltIn, [Operand::BuiltIn(BuiltIn::WorkgroupSize)]);
	let pv3u = b.type_pointer(None, StorageClass::Input, v3u);
	let global_invocation_id_ref = b.variable(pv3u, None, StorageClass::Input, None);
	b.decorate(global_invocation_id_ref, BuiltIn, [Operand::BuiltIn(BuiltIn::GlobalInvocationId)]);
	let output_types = {
		let mut types = Types(ast.input.iter().copied().map(Some).chain((ast.input.len()..ast.values.len()).map(|_| None)).collect());
		for s in &*ast.statements { types.push(s); }
		map(&*ast.output, |output| {
			types.expr(output);
			types.rtype(output)
		})
	};

	let constants_types = &ast.input[0..constants_len];
	let constants = {
		let constants = map(constants_types, |constant| stype(&mut b, constant));
		let block_struct_constants = b.type_struct(constants.into_vec());
		b.decorate(block_struct_constants, Block, []);
		for (member, offset) in constants_types.iter().scan(0, |offset, t| { let current = *offset; *offset+=t.bits()/8; Some(current)}).enumerate() { b.member_decorate(block_struct_constants, member as u32, Offset, [offset.into()]); }
		let constant_pointer_to_block_struct_constants = b.type_pointer(None, StorageClass::PushConstant, block_struct_constants);
		b.variable(constant_pointer_to_block_struct_constants, None, StorageClass::PushConstant, None)
	};

	let mut storage_buffer_pointer_to_array = {
		let mut cache: linear_map::LinearMap<ast::Type, Type> = default();
		move |b:&mut rspirv::Builder, rtype:&ast::Type| *cache.entry(*rtype).or_insert_with(|| {
			let stype = stype(b, &rtype);
			let array = b.type_runtime_array(stype);
			b.decorate(array, ArrayStride, [rtype.bytes().into()]);
			let block_struct_array = b.type_struct([array]);
			b.decorate(block_struct_array, Block, []);
			b.member_decorate(block_struct_array, 0, Offset, [0u32.into()]);
			b.type_pointer(None, StorageClass::StorageBuffer, block_struct_array)
		})
	};

	let input_pointers = map(&ast.input[constants_len..], |input| {
		let pointer_type = storage_buffer_pointer_to_array(&mut b, input);
		let pointer = b.variable(pointer_type, None, StorageClass::StorageBuffer, None);
		b.decorate(pointer, NonWritable, []);
		(*input, pointer)
	});
	let output = map(&*output_types,|output| {
		let pointer_type = storage_buffer_pointer_to_array(&mut b, &output);
		let pointer = b.variable(pointer_type, None, StorageClass::StorageBuffer, None);
		b.decorate(pointer, NonReadable, []);
		(*output, pointer)
	});
	for (binding, &(_,variable)) in input_pointers.iter().chain(&*output).enumerate() {
		b.decorate(variable, DescriptorSet, [0u32.into()]);
		b.decorate(variable, Binding, [(binding as u32).into()]);
	}
	b.begin_block(None)?;
	let global_invocation_id3 = b.load(v3u, None, global_invocation_id_ref, None, [])?;
	let id = b.composite_extract(u32, None, global_invocation_id3, [0])?;
	let constant_values = map(constants_types.into_iter().enumerate(), |(member, constant_type)| {
		let stype = stype(&mut b, constant_type);
		let constant_pointer_to_constant = b.type_pointer(None, StorageClass::PushConstant, stype);
		let member = b.constant_u32(u32, member as u32);
		let pointer = b.access_chain(constant_pointer_to_constant, None, constants, [member]).unwrap();
		(*constant_type, b.load(stype, None, pointer, None, []).unwrap())
	});
	let index0 = b.constant_u32(u32, 0);
	let input_values = map(&*input_pointers, |&(rtype, input)| {
		let stype = stype(&mut b, &rtype);
		let storage_buffer_pointer = b.type_pointer(None, StorageClass::StorageBuffer, stype);
		let input = b.access_chain(storage_buffer_pointer, None, input, [index0, id]).unwrap();
		(rtype, b.load(stype, None, input, Some(MemoryAccess::NONTEMPORAL), []).unwrap())
	});
	let mut b = Builder{
		builder: b,
		gl, //debug,
		values: list(constant_values.into_vec().into_iter().chain(input_values.into_vec()).map(Some).chain((ast.input.len()..ast.values.len()).map(|_| None))),
		constants_f32: default(), constants_f64: default(), /*expressions: default(),*/ names: &ast.values};
	for s in &*ast.statements { b.extend(s); }
	for (expr, &(rtype, output)) in ast.output.iter().zip(&*output) {
		let value = b.expr(expr);
		let stype = b.stype(&rtype);
		let storage_buffer_pointer = b.type_pointer(None, StorageClass::StorageBuffer, stype);
		let output = b.access_chain(storage_buffer_pointer, None, output, [index0, id]).unwrap();
		b.store(output, value, Some(MemoryAccess::NONTEMPORAL), []).unwrap();
	}
	b.ret()?;
	b.end_function()?;
	let interface = [global_invocation_id_ref, constants].into_iter().chain(input_pointers.into_vec().into_iter().map(|(_,v)|v)).chain(output.into_vec().into_iter().map(|(_,v)|v));
	b.entry_point(spirv::ExecutionModel::GLCompute, f, "main", list(interface));
	let code = ::rspirv::binary::Assemble::assemble(&b.builder.module());
	if false {
		Ok(code.into())
	} else {
		let path = std::env::var("XDG_RUNTIME_DIR").unwrap()+"/spirv.spv";
		let path = std::path::Path::new(&path);
		pub fn as_u8<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		let bytes = as_u8(&code);
		if !path.exists() || bytes != std::fs::read(path).unwrap() {
			std::fs::write(path, bytes).unwrap();
			/*{use spirv_tools::{*, opt::*, binary::*};
				let code = create(Some(TargetEnv::Vulkan_1_2)).optimize(code, &mut |_| {}, Some(Options{preserve_spec_constants: true, ..Default::default()})).map_err(|e| e.diagnostic.unwrap().message).unwrap();
				if let Binary::OwnedU32(code) = code { Ok(code.into()) } else { unreachable!() }
			}*/
		}
		if true {
			Ok(code.into())
		} else {
			let opt = std::env::var("XDG_RUNTIME_DIR").unwrap()+"/spirv-opt.spv";
			let opt = std::path::Path::new(&opt);
			if !opt.exists() || opt.metadata().unwrap().modified().unwrap() < path.metadata().unwrap().modified().unwrap() {
				let [path,opt] = [path,opt].map(|path| path.to_str().unwrap());
				assert!(std::process::Command::new("spirv-opt").args(["--skip-validation","--target-env=vulkan1.2",&format!("-Oconfig={}/.spirv-opt",std::env::var("HOME").unwrap()),path,"-o",opt]).status().unwrap().success());
			}
			pub fn as_u32(slice: &[u8]) -> &[u32] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u32, slice.len() / std::mem::size_of::<u8>())} }
			Ok(as_u32(&std::fs::read(std::env::var("XDG_RUNTIME_DIR").unwrap()+"/spirv-opt.spv").unwrap()).into())
		}
	}
}
