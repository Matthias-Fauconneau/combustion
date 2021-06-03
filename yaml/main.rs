#[throws] fn main() {
	let model = parse(Loader::load_from_str(std::str::from_utf8(&std::fs::read("/usr/share/cantera/data/gri30.yaml")?)?).unwrap()).unwrap();
	/*let Model{time_step, state, ..}  = ::ron::de::from_bytes(&std::fs::read("CH4.ron")?)?;
	model.time_step = time_step;
	model.state = state;*/
	println!("{}", to_string(&model)?);
}
