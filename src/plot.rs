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
		if x < target.size.x && y < target.size.y { target[xy{x,y}].saturating_add( (opacity*color).into() ); }
	}
	let (i0, intery) = {
		let xend = f32::round(p0.x);
		let yend = p0.y + gradient * (xend - p0.x);
    let xgap = 1. - fract(p0.x + 1./2.);
		blend(target, xend as u32, yend as u32, color, (1.-fract(yend)) * xgap * opacity, transpose);
    blend(target, xend as u32, yend as u32+1, color, fract(yend) * xgap * opacity, transpose);
    (xend as i32, yend + gradient)
	};
	let i1 = {
		let xend = f32::round(p1.x);
    let yend = p1.y + gradient * (xend - p1.x);
    let xgap = fract(p1.x + 1./2.);
    blend(target, xend as u32, yend as u32, color, (1.-fract(yend)) * xgap * opacity, transpose);
    blend(target, xend as u32, yend as u32+1, color, fract(yend) * xgap * opacity, transpose);
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

#[allow(non_camel_case_types)] #[derive(PartialEq,Clone,Copy,PartialOrd,Debug,serde::Deserialize)] pub struct real(pub f32);
impl std::ops::Neg for real { type Output = Self; fn neg(self) -> Self { real(-self.0) } }
impl std::ops::Neg for &real { type Output = real; fn neg(self) -> real { real(-self.0) } }
impl std::ops::Add<Self> for real { type Output = Self; fn add(self, b: Self) -> Self { real(self.0.add(b.0)) } }
impl std::ops::AddAssign<Self> for real { fn add_assign(&mut self, b: Self) { self.0.add_assign(b.0) } }
impl std::iter::Sum<Self> for real { fn sum<I:Iterator<Item=Self>>(iter: I) -> Self { real(iter.map(|real(f)| f).sum()) } }
impl std::iter::Sum<&'t Self> for real { fn sum<I:Iterator<Item=&'t Self>>(iter: I) -> Self { real(iter.map(|real(f)| f).sum()) } }
impl std::ops::Sub<Self> for real { type Output = Self; fn sub(self, b: Self) -> Self { real(self.0.sub(b.0)) } }
impl std::ops::Sub<&real> for real { type Output = Self; fn sub(self, b: &real) -> Self { real(self.0.sub(b.0)) } }
impl std::ops::Sub<real> for &real { type Output = real; fn sub(self, b: real) -> Self::Output { real(self.0.sub(b.0)) } }
impl std::ops::Sub<&real> for &real { type Output = real; fn sub(self, b: &real) -> Self::Output { real(self.0.sub(b.0)) } }
impl std::ops::Mul<Self> for real { type Output = Self; fn mul(self, b: Self) -> Self { real(self.0.mul(b.0)) } }
impl std::ops::Mul<&real> for real { type Output = Self; fn mul(self, b: &real) -> Self { real(self.0.mul(b.0)) } }
impl std::ops::Mul<real> for &real { type Output = real; fn mul(self, b: real) -> Self::Output { real(self.0.mul(b.0)) } }
impl std::ops::Mul<&real> for &real { type Output = real; fn mul(self, b: &real) -> Self::Output { real(self.0.mul(b.0)) } }
impl num::Zero for real { const ZERO : Self = real(num::Zero::ZERO); }
impl std::ops::Div<Self> for real { type Output = Self; fn div(self, b: Self) -> Self { assert!(!num::IsZero::is_zero(&b)); real(self.0.div(b.0)) } }
impl std::ops::Div<real> for &real { type Output = real; fn div(self, b: real) -> Self::Output { assert!(!num::IsZero::is_zero(&b)); real(self.0.div(b.0)) } }
impl std::ops::Div<&real> for real { type Output = Self; fn div(self, b: &real) -> Self { assert!(!num::IsZero::is_zero(b)); real(self.0.div(b.0)) } }
impl std::cmp::Eq for real {}
impl std::cmp::Ord for real { fn cmp(&self, other: &Self) -> std::cmp::Ordering { self.partial_cmp(other).unwrap() } }
impl real {
	pub fn recip(self) -> real { real(f32::recip(self.0)) }
	pub fn exp(self) -> real { real(f32::exp(self.0)) }
	pub fn pow(self, n: real) -> real { real(f32::powf(self.0, n.0)) }
	pub fn powi(self, n: i32) -> real { real(f32::powi(self.0, n)) }
	pub fn ln(self) -> real { real(f32::ln(self.0)) }
	pub fn log10(self) -> real { real(f32::log10(self.0)) }
}

type Frame = (real, Box<[Box<[real]>]>);
pub struct Plot<'t, I: Iterator<Item=Frame>>{
	pub keys: &'t [&'t [&'t str]],
	pub values: Vec<Frame>,
	pub source: I
}
impl<'t, I: Iterator<Item=Frame>> ui::widget::Widget for Plot<'t, I> {
	#[fehler::throws(anyhow::Error)] fn paint(&mut self, mut target: &mut ui::widget::Target) {
		let Self{keys, values, source} = self;
		dbg!(keys);
		if let Some(v) = source.next() { values.push(v); }
		let x_minmax = vector::minmax(values.iter().map(|&(x,_)| x)).unwrap();
		let mut serie_of_sets = values.iter().map(|(_,sets)| sets);
		let mut sets_minmax = serie_of_sets.next().unwrap().iter().map(|set| vector::minmax(set.iter().copied()).unwrap()).collect::<Box<[_]>>();
		for sets in serie_of_sets { for (minmax, set) in sets_minmax.iter_mut().zip(sets.iter()) { *minmax = minmax.minmax(vector::minmax(set.iter().copied()).unwrap()) } }
		let size = target.size;
		use itertools::Itertools;
		values.iter().map(|(x, sets)| {
			let map = |vector::MinMax{min,max} : vector::MinMax<real>, v: real| -> f32 { (v-min / (max-min)).0 };
			sets.iter().zip(sets_minmax.iter()).map(move |(set, &minmax)| set.iter().map(move |&y| xy{x: map(x_minmax, *x), y: map(minmax, y)}*size.into())).flatten()
		}).flatten()
		.tuple_windows().for_each(|(p0, p1)| self::line(&mut target, p0, p1, num::zero(), 1.));
	}
}
