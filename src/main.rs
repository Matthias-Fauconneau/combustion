#![feature(type_ascription)]
use {fehler::throws, error::Error, combustion::{*, Property::*}};

#[throws] fn main() {
	let model = Model::new(model::Model::new(&std::fs::read("CH4+O2.ron")?)?);
	let (_traps, (function, size), _rate) = model.rate::<{Volume}>();
	let ref state = Simulation::new(&std::fs::read("CH4+O2.ron")?)?.state;
	let mut derivative = /*Derivative*/StateVector::<{Volume}>(std::iter::repeat(0.).take(model.len()).collect());
	/*{
		let handler = move |info:&libc::siginfo_t| {
			let offset = unsafe{(info.si_addr() as *const u8).offset_from(function as *const u8)};
			dbg!(traps.iter().find(|&x| dbg!(x.code_offset) == offset as u32).unwrap().source_location);
			std::process::abort()
		};
		unsafe{signal_hook_registry::register_unchecked(libc::SIGILL, handler)}?;
		rate(state.constant(), &state.into(), &mut derivative);
	}*/
	{
		let function = unsafe{std::slice::from_raw_parts(function as *const u8, size)}; // refcheck: leaked from dropped JITModule
		let (constant, ref state, derivative) = (state.constant::<{Volume}>(), state.into(): StateVector::<{Volume}>, &mut derivative);
		let constant = constant.0 as f32;
		let mut guest = x86emu::State::new();
		static HEAP_BASE : u64 = 0x2_0000_0000;
    static HEAP_SIZE : usize = 0x0_0010_0000;
    guest.memory.host_allocate_physical(HEAP_BASE, HEAP_SIZE);
    let mut heap_next = 0;
    fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		let mut new = |slice| { let offset = HEAP_BASE+heap_next; guest.memory.write_unaligned_bytes(offset, slice); heap_next += slice.len() as u64; offset };
		let state = new(as_bytes(&state.0));
		let derivative = new(as_bytes(&derivative.0));
		guest.call(function, &[state as i64, derivative as i64], &[constant]);
	}
	//println!("{:?}", {rate(state.constant(), &state.into(), &mut derivative); derivative});
}
