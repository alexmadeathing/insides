use criterion::{black_box, criterion_group, criterion_main, Criterion};

use dilate::DilationMethod;
use insides::{Hilbert, Expand, Encoding, Coord, Index};

pub trait NumTraits {
    fn from_usize(value: usize) -> Self;
}

macro_rules! impl_num_traits {
    ($($t:ty),+) => {$(
        impl NumTraits for $t {
            fn from_usize(value: usize) -> Self {
                value as Self
            }
        }
    )+};
}

impl_num_traits!(u8, u16, u32, u64, u128);

#[inline]
fn bench_index_to_coords<E, const D: usize>()
where
    E: Encoding<D>,
    E::Index: NumTraits,
{
    for i in 0..1000 {
        black_box(E::from_index(black_box(E::Index::from_usize(i))).coords());
    }
}

type HilbertD2 = Hilbert::<Expand<u8, 2>, 2>;
type HilbertD3 = Hilbert::<Expand<u8, 3>, 3>;
type HilbertD4 = Hilbert::<Expand<u8, 4>, 4>;
type HilbertD5 = Hilbert::<Expand<u8, 5>, 5>;
type HilbertD6 = Hilbert::<Expand<u8, 6>, 6>;
type HilbertD7 = Hilbert::<Expand<u8, 7>, 7>;
type HilbertD8 = Hilbert::<Expand<u8, 8>, 8>;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Hilbert D2 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD2, 2>));
    c.bench_function("Hilbert D3 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD3, 3>));
    c.bench_function("Hilbert D4 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD4, 4>));
    c.bench_function("Hilbert D5 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD5, 5>));
    c.bench_function("Hilbert D6 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD6, 6>));
    c.bench_function("Hilbert D7 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD7, 7>));
    c.bench_function("Hilbert D8 index_to_coords", |b| b.iter(bench_index_to_coords::<HilbertD8, 8>));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);