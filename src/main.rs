#[fehler::throws(anyhow::Error)] fn main() { ui::app::run(combustion::Simulation::new(&std::fs::read("H2+O2.ron")?)?.plot())? }
