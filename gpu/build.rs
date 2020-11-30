#[fehler::throws(Box<dyn std::error::Error>)] fn main() { spirv_builder::SpirvBuilder::new("f").build()?; }
