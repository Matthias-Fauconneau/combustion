use {::xy::{xy,vec2}, image::{Image, bgrf, bgra8}};

fn line(target: &mut Image<&mut[bgra8]>, p0: vec2, p1: vec2, color: bgrf) {
	use num::{abs, fract};
  let d = p1 - p0;
  let (transpose, p0, p1, d) = if abs(d.x) < abs(d.y) { (true, p0.yx(), p1.yx(), d.yx()) } else { (false, p0, p1, d) };
  if d.x == 0. { return; } // p0==p1
  let (p0, p1) = if p0.x > p1.x { (p1, p0) } else { (p0, p1) };
  let gradient = d.y / d.x;
  fn blend(target: &mut Image<&mut[bgra8]>, x: u32, y: u32, color: bgrf, coverage: f32, transpose: bool) {
		let xy{x,y} = if transpose { xy{x: y, y: x} } else { xy{x,y} };
		if x < target.size.x && y < target.size.y { target[xy{x,y}].saturating_add_assign( (coverage*color).into() ); }
	}
	let (i0, intery) = {
		let xend = f32::round(p0.x);
		let yend = p0.y + gradient * (xend - p0.x);
    let xgap = p0.x - xend;
    let fract_yend = yend - f32::floor(yend);
    blend(target, xend as u32, yend as u32, color, (1.-fract_yend) * xgap, transpose);
    blend(target, xend as u32, yend as u32+1, color, fract_yend * xgap, transpose);
    (xend as i32, yend + gradient)
	};
	let i1 = {
		let xend = f32::round(p1.x);
    let yend = p1.y + gradient * (xend - p1.x);
    let xgap = p1.x - xend;
    let fract_yend = yend - f32::floor(yend);
    blend(target, xend as u32, yend as u32, color, (1.-fract_yend) * xgap, transpose);
    blend(target, xend as u32, yend as u32+1, color, fract_yend * xgap, transpose);
    xend as u32
	};
	let x = i0+1;
	let (mut intery, mut x) = if x < 0 { (intery+(0-x as i32) as f32 * gradient, 0) } else { (intery, x as u32) };
	while x < i1.min(if transpose { target.size.y } else { target.size.x }) {
		blend(target, x, intery as u32, color, 1.-fract(intery), transpose);
    blend(target, x, intery as u32+1, color, fract(intery), transpose);
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
	let Self{keys: _, values} = self;

	let ticks = |vector::MinMax{max,..}| {
		let log10 = real::log10(real::abs(max));
		let exp_fract_log10 = real::exp10(log10 - real::floor(log10));
		let (max, tick_count) = *[(1.,5),(1.2,6),(1.4,7),(2.,10),(2.5,5),(3.,3),(3.2,8),(4.,8),(5.,5),(6.,6),(8.,8),(10.,5)].iter().find(|(max,_)| exp_fract_log10-real::log2(real(-52.)) <= real(*max)).unwrap();
		let max = real(max)*real::exp10(real::floor(log10));
		let precision = (real::ceil(-real::log10(max/real(tick_count as f32))).0 as usize).max(1);
		iter::from_iter((0..=tick_count).map(|i| max*real(i as f32)/real(tick_count as f32)).map(|value| (value, format!("{:.1$e}", value, precision))))
	};

	let x_minmax = vector::minmax(values.iter().map(|&(x,_)| x)).unwrap();
	let x_labels = ticks(x_minmax);
	let mut x_ticks = iter::from_iter(x_labels.iter().map(|(_,label)| ui::text::View::new(ui::text::Borrowed::new(label))));
	let x_label_size = vector::minmax(x_ticks.iter_mut().map(|tick| tick.size())).unwrap().max;

	let mut serie_of_sets = values.iter().map(|(_,sets)| sets);
	let mut sets_minmax = serie_of_sets.next().unwrap().iter().map(|set| vector::minmax(set.iter().copied()).unwrap()).collect::<Box<[_]>>();
	for sets in serie_of_sets { for (minmax, set) in sets_minmax.iter_mut().zip(sets.iter()) { *minmax = minmax.minmax(vector::minmax(set.iter().copied()).unwrap()) } }
	let sets_labels = iter::from_iter(sets_minmax.iter().map(|&minmax| ticks(minmax)));
	let mut sets_ticks = iter::from_iter(sets_labels.iter().map(|set_labels| iter::from_iter(set_labels.iter().map(|(_,label)| ui::text::View::new(ui::text::Borrowed::new(label))))));
	let sets_label_size = iter::from_iter(sets_ticks.iter_mut().map(|set_ticks| vector::minmax(set_ticks.iter_mut().map(|tick| tick.size())).unwrap().max));

	let size = target.size;
	let top = 0;
	let bottom = ui::text::fit_width(size.x/(x_labels.len() as u32), x_label_size).y;
	let [left, right] = iter::array::generate(|i| if let (Some(labels), Some(&label_size)) = (sets_labels.get(i), sets_label_size.get(i)) { ui::text::fit_height(size.y/(labels.len() as u32), label_size).x } else { 0 });

	let fg = bgra8{b:0xFF, g:0xFF, r:0xFF, a:0xFF};
	//let rect = |top_le| target.slice_mut(top_left, vector::component_wise_min(bottom_right, target.size)-top_left).set(|_| fg);
	//line(&mut target, xy{x: left, y: size.y-bottom}, xy{x: right, y: size.y-bottom}, fg);
	target.slice_mut(xy{x: left, y: size.y-bottom}, xy{x: right-left, y: 1}).set(|_| fg);
	//line(&mut target, xy{x: left, y: top}, xy{x: left, y: size.y-bottom}, fg);
	target.slice_mut(xy{x: left, y: top}, xy{x: 1, y: size.y-bottom-top}).set(|_| fg);
	//line(&mut target, xy{x: size.x-right, y: top}, xy{x: size.x-right, y: size.y-bottom}, fg);
	target.slice_mut(xy{x: right, y: top}, xy{x: 1, y: size.y-bottom-top}).set(|_| fg);

	fn map(vector::MinMax{min,max} : vector::MinMax<real>, v: real) -> Option<f32> { if min < max { Some(((v-min) / (max-min)).0) } else { None } }
	let map_x = |minmax, v| Some(left as f32+map(minmax, v)?*(size.x-left-right-1) as f32);
	let map_y = |minmax, v| Some((size.y-bottom) as f32-map(minmax, v)?*(size.y-top-bottom-1) as f32);

	let tick_length = 16;
	for (&(value,_), tick) in x_labels.iter().zip(x_ticks.iter_mut()) {
		let p = xy{x: map_x(x_minmax, value).unwrap() as u32, y: size.y-bottom};
		target.slice_mut(p-xy{x:0, y:tick_length}, xy{x:1, y:tick_length}).set(|_| fg);
		let scale = num::Ratio{num: size.x/(x_labels.len() as u32)-1, div: x_label_size.x-1};
		let p = p - xy{x: scale*tick.size().x/2, y: 0};
		tick.paint(&mut target, scale, p);
	}

	for i in 0..2 {
		if let (Some(&minmax), Some(labels), Some(ticks), Some(label_size)) = (sets_minmax.get(i), sets_labels.get(i), sets_ticks.get_mut(i), sets_label_size.get(i)) {
			for (&(value,_), tick) in labels.iter().zip(ticks.iter_mut()) {
				let p = xy{x: left, y: map_y(minmax, value).unwrap() as u32};
				target.slice_mut(p-xy{x:0,y:[0,tick_length][i]}, xy{x:tick_length, y:1}).set(|_| fg);
				let scale = num::Ratio{num: size.y/(labels.len() as u32)-1, div: label_size.y-1};
				let p = p - xy{x:0, y: scale*tick.size().y/2};
				tick.paint(&mut target, scale, p);
			}
		}
	}

	let fg = bgrf{b:1., g:1., r: 1.};
	use itertools::Itertools;
	values.iter().map(|(x, sets)| sets.iter().zip(sets_minmax.iter()).map(move |(set, &minmax)| set.iter().map(move |&y| Some(xy{x: map_x(x_minmax, *x)?, y: map_y(minmax, y)?})))).flatten()
		.tuple_windows().for_each(|(s0, s1)| s0.zip(s1).for_each(|line| if let (Some(p0), Some(p1)) = line { self::line(&mut target, p0, p1, fg) }));
}
}
