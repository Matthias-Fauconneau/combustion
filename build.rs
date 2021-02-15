#[fehler::throws(error::Error)] fn main() { println!("cargo:rustc-env=SPECIES_LEN={}", system::System::new(&system::default()?)?.species.len()); }
