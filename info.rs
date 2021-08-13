#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip)]
#![allow(non_snake_case,non_upper_case_globals)]
mod yaml;
use {std::{iter::zip, cmp::max}, anyhow::Error, std::env::*, combustion::*};
#[fehler::throws] fn main() {
	for path in args().skip(1) {
		let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
		let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path)?)?)?;
		let model = yaml::parse(&model);
		let (ref species, _, _, reactions, _) = new(&model);
		let mut live = vec![vec![false; species.len()]; 2];
		let mut max_live_count = 0;
		for (i, r) in reactions.iter().enumerate() {
			let f = |c| if c==1||c==2 { 0 } else if c==-1||c==-2 { 1 } else { panic!("{c}") };
			for (k,&c) in r.net.iter().enumerate() { if c!=0 { live[f(c)][k] = true; } }
			let live_count = live.iter().map(|l| l.iter()).flatten().filter(|&&x| x).count();
			max_live_count = max(max_live_count, live_count);
			//print!("{} ", live_count);
			let mut used = vec![vec![false; species.len()]; 2];
			for r in reactions[i+1..].iter() {
				for (k,&c) in r.net.iter().enumerate() { if c!=0 { used[f(c)][k] = true; } }
			}
			// Kill defs after last use
			for (live, used) in zip(live.iter_mut(), used) { for (live, used) in zip(live.iter_mut(), used) { if !used { *live = false; } } }
		}
		println!("{path} {} {max_live_count}", species.len());
		/*use reaction::bucket;
		use itertools::Itertools;
		let n = bucket(reactions.iter().map(|r| r.reactants.iter().enumerate().filter(|(_,&c)| c>0).map(|(s,_)| s)).flatten());
		let d = bucket(reactions.iter().map(|r| r.products.iter().enumerate().filter(|(_,&c)| c>0).map(|(s,_)| s)).flatten());
		println!("{:?}", n.into_iter().map(|(_,r)| r.len()).format(" "));
		println!("{:?}", d.into_iter().map(|(_,r)| r.len()).format(" "));*/
		/*let n = n.into_iter().filter(|(_,r)| !r.is_empty()).count();
		let d = d.into_iter().filter(|(_,r)| !r.is_empty()).count();
		println!("{}+{}-{}={}", n, d, species.len(), n+d-species.len());*/
		//println!("{:?}", bucket(reactions.iter().map(|r| r.products.clone())).into_iter().map(|(_,r)| r.len()).format(" "));
		//let reused_divisors = bucket(reactions.iter().map(|r| r.products.clone())).into_iter().map(|(_,r)| r.len()-1).sum::<usize>();
		//let diff = reused_divisors as i64-species.len() as i64;
		//println!("{path}: {reused_divisors}-{}={} {:.0}%{}", species.len(), diff, 100.*diff as f64/reactions.len() as f64, reactions.len());
		//let diff = reactions.len()-2*species.len();
		//println!("{path}: {}-{}={} {:.0}%{}", reactions.len(), 2*species.len(), diff, 100.*diff as f64/reactions.len() as f64, reactions.len());
	}
}
