fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = combustion_yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read("/usr/share/cantera/data/LiDryer.yaml")?)?)?;
	let model = combustion_yaml::parse(&model)?;
	/*let Model{time_step, state, ..}  = ::ron::de::from_bytes(&std::fs::read("CH4.ron")?)?;
	model.time_step = time_step;
	model.state = state;*/
	println!("{}", combustion_yaml::to_string(&model)?);
	Ok(())
}
