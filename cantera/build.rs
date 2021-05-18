fn main() -> std::result::Result<(),std::io::Error> {
	println!("cargo:rustc-link-lib=stdc++");
	println!("cargo:rustc-link-lib=yaml-cpp");
	println!("cargo:rustc-link-lib=sundials_cvodes");
	println!("cargo:rustc-link-lib=cantera");
	Ok(())
}
