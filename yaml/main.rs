fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args()[0])?)?)?;
	let model = yaml_model::parse(&model)?;
	println!("{}", yaml_model::to_string(&model)?);
	Ok(())
}
