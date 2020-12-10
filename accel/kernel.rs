#[accel::kernel] pub unsafe fn dt(len: usize, pressure: *const f32, temperature: *mut f32, volume: *mut f32, amounts: *const *mut f32) {
	let i = accel_core::index();
	if (i as usize) < len {
			//*c.offset(i) = *a.offset(i) + *b.offset(i);
	}
}
