#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias)]
#![allow(non_snake_case,non_upper_case_globals)]
mod yaml;
use {anyhow::Error, std::env::*, combustion::*};
#[fehler::throws] fn main() {
	for path in args().skip(1) {
		let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
		let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path)?)?)?;
		let model = yaml::parse(&model);
		let (ref species, _, _, reactions, _) = new(&model);
		//use itertools::Itertools; println!("{:?}", bucket(reactions.iter().map(|r| r.products.clone())).into_iter().map(|(_,r)| r.len()).format(" "));
		//use reaction::bucket; let reused_divisors = bucket(reactions.iter().map(|r| r.products.clone())).into_iter().map(|(_,r)| r.len()-1).sum::<usize>();
		//let diff = reused_divisors as i64-species.len() as i64;
		//println!("{path}: {reused_divisors}-{}={} {:.0}%{}", species.len(), diff, 100.*diff as f64/reactions.len() as f64, reactions.len());
		let diff = reactions.len()-2*species.len();
		println!("{path}: {}-{}={} {:.0}%{}", reactions.len(), 2*species.len(), diff, 100.*diff as f64/reactions.len() as f64, reactions.len());
	}
}
