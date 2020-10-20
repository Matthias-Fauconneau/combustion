use {::xy::{xy,vec2}, image::{Image, bgrf, bgra8}};

fn line(target: &mut Image<&mut[bgra8]>, p0: vec2, p1: vec2, color: bgrf, opacity: f32) {
	use num::{abs, fract};
  let d = p1 - p0;
  let (transpose, p0, p1, d) = if abs(d.x) < abs(d.y) { (true, p0.yx(), p1.yx(), d.yx()) } else { (false, p0, p1, d) };
  if d.x == 0. { return; } // p0==p1
  let (p0, p1) = if p0.x > p1.x { (p1, p0) } else { (p0, p1) };
  let gradient = d.y / d.x;
  fn blend(target: &mut Image<&mut[bgra8]>, x: u32, y: u32, color: bgrf, opacity: f32, transpose: bool) {
		let xy{x,y} = if transpose { xy{x: y, y: x} } else { xy{x,y} };
		if x < target.size.x && y < target.size.y { target[xy{x,y}].saturating_add_assign( (opacity*color).into() ); }
	}
	let (i0, intery) = {
		let xend = f32::round(p0.x);
		let yend = p0.y + gradient * (xend - p0.x);
    let xgap = p0.x - xend;
    let fract_yend = yend - f32::floor(yend);
    blend(target, xend as u32, yend as u32, color, (1.-fract_yend) * xgap * opacity, transpose);
    blend(target, xend as u32, yend as u32+1, color, fract_yend * xgap * opacity, transpose);
    (xend as i32, yend + gradient)
	};
	let i1 = {
		let xend = f32::round(p1.x);
    let yend = p1.y + gradient * (xend - p1.x);
    let xgap = p1.x - xend;
    let fract_yend = yend - f32::floor(yend);
    blend(target, xend as u32, yend as u32, color, (1.-fract_yend) * xgap * opacity, transpose);
    blend(target, xend as u32, yend as u32+1, color, fract_yend * xgap * opacity, transpose);
    xend as u32
	};
	let x = i0+1;
	let (mut intery, mut x) = if x < 0 { (intery+(0-x as i32) as f32 * gradient, 0) } else { (intery, x as u32) };
	while x < i1.min(if transpose { target.size.y } else { target.size.x }) {
		blend(target, x, intery as u32, color, (1.-fract(intery)) * opacity, transpose);
    blend(target, x, intery as u32+1, color, fract(intery) * opacity, transpose);
    intery += gradient;
    x += 1;
	}
}

use num::real;
type Frame = (real, Box<[Box<[real]>]>);
pub struct Plot<'t>{
	pub keys: Box<[&'t [&'t str]]>,
	pub values: Vec<Frame>,
}
impl ui::widget::Widget for Plot<'_> {
	#[fehler::throws(anyhow::Error)] fn paint(&mut self, mut target: &mut ui::widget::Target) {
		let Self{keys, values} = self;
		let _ = keys;
		let x_minmax = vector::minmax(values.iter().map(|&(x,_)| x)).unwrap();
		let mut serie_of_sets = values.iter().map(|(_,sets)| sets);
		let mut sets_minmax = serie_of_sets.next().unwrap().iter().map(|set| vector::minmax(set.iter().copied()).unwrap()).collect::<Box<[_]>>();
		for sets in serie_of_sets { for (minmax, set) in sets_minmax.iter_mut().zip(sets.iter()) { *minmax = minmax.minmax(vector::minmax(set.iter().copied()).unwrap()) } }
		fn map(vector::MinMax{min,max} : vector::MinMax<real>, v: real) -> Option<f32> { if min < max { Some(((v-min) / (max-min)).0) } else { None } }
		let size = target.size;
		use itertools::Itertools;
		values.iter().map(|(x, sets)| sets.iter().zip(sets_minmax.iter()).map(move |(set, &minmax)| set.iter().map(move |&y| Some(xy{x: map(x_minmax, *x)?, y: map(minmax, y)?}*size.into()))).flatten())
			.tuple_windows().for_each(|(s0, s1)| s0.zip(s1).for_each(|line| if let (Some(p0), Some(p1)) = line { self::line(&mut target, p0, p1, bgrf{b:1., g:1., r: 1.}, 1.) }));
	}
}
