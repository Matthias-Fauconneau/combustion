pub struct Device(Option<(accel::Device, accel::Context)>);
impl std::ops::Deref for Device { type Target = accel::Context; fn deref(&self) -> &Self::Target { &self.0.as_ref().unwrap().1 } }
impl Device { #[fehler::throws(anyhow::Error)] pub fn new() -> Self {
	std::process::Command::new("gpu-on").spawn()?.wait()?;
	assert!(std::str::from_utf8(&std::fs::read("/proc/modules")?)?.lines().any(|line| line.starts_with("nvidia ")));
	let device = accel::Device::nth(0)?;
	let context = device.create_context();
	Self(Some((device, context)))
}}
impl Drop for Device { fn drop(&mut self) {
	self.0 = None; // Drop accel before trying to unload nvidia module and switching GPU off
	std::process::Command::new("gpu-off").spawn().unwrap().wait().unwrap();
	assert!(!std::str::from_utf8(&std::fs::read("/proc/modules").unwrap()).unwrap().lines().any(|line| line.starts_with("nvidia ")));
} }
