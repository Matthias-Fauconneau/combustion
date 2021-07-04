#![feature(trait_alias)]#![allow(non_snake_case)]
use {num::sq, iter::{IntoIterator, ConstRange, IntoConstSizeIterator, Dot}};
pub fn evaluate<const N: usize>(P: &[f64; N], x: f64) -> f64 { Dot::dot(P, ConstRange.map(|k| x.powi(k as i32))) }

pub trait Vector<const N: usize> = IntoConstSizeIterator<N>+IntoIterator<Item=f64>;
pub fn weighted_regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>, w: impl Vector<N>) -> [f64; D] {
	use nalgebra::{DMatrix, DVector, SVD};
	let ref w = DVector::from_iterator(N, w.into_iter());
	let ref x = x.collect();
	let A = DMatrix::from_iterator(N, D, (0..D).map(|k| x.iter().zip(w.iter()).map(move |(x, w)| w*x.powi(k as i32))).flatten());
	let b = DVector::from_iterator(N, y.into_iter().zip(w.iter()).map(|(y, w)| w*y));
	SVD::new(A, true, true).solve(&b, f64::EPSILON).unwrap().as_slice().try_into().unwrap()
}

// Regression with 1/y² weights (towards relative vertical error)
pub fn regression<const D: usize, const N: usize>(x: impl Vector<N>, y: impl Vector<N>+Clone) -> [f64; D] { weighted_regression(x, y.clone(), y.map(|y| 1./sq(y))) }
pub fn fit<T: Vector<N>+Copy, X: Fn(f64)->f64, Y: Fn(f64)->f64+Copy, const D: usize, const N: usize>(t: T, x: X, y: Y) -> [f64; D] { regression(t.map(x), t.map(y)) }
