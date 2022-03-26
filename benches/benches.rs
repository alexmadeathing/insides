use criterion::{black_box, criterion_group, criterion_main, Criterion};

use insides::{Hilbert, Expand, Encoding};

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

pub fn criterion_benchmark(c: &mut Criterion) {
    let coord_bits: usize = 8;
    let coord_length: usize = 1 << coord_bits;
    let index_length: usize = coord_length * coord_length;

    c.bench_function("insides: 2d hilbert index to coords", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(Hilbert::<Expand<u32, 2>, 2>::from_index(black_box(i as u64)).coords());
            }
        })
    });

    c.bench_function("fast_hilbert: 2d hilbert index to coords", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(fast_hilbert::h2xy::<u32>(black_box(i as u64)));
            }
        })
    });

    c.bench_function("hilbert_curve: 2d hilbert index to coords", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_curve::convert_1d_to_2d(black_box(i), black_box(coord_length)));
            }
        })
    });

    c.bench_function("hilbert_2d: 2d hilbert index to coords", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_2d::h2xy_discrete(black_box(i), black_box(coord_bits), black_box(hilbert_2d::Variant::Hilbert)));
            }
        })
    });

    c.bench_function("hilbert_index: 2d hilbert index to coords", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_index::FromHilbertIndex::<2>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
            }
        })
    });

    c.bench_function("hilbert: 2d hilbert index to coords", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(2)));
            }
        })
    });

    c.bench_function("insides: 2d hilbert coords to index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(Hilbert::<Expand<u32, 2>, 2>::from_coords(black_box([x as u32, y as u32])).index());
                }
            }
        })
    });

    c.bench_function("fast_hilbert: 2d hilbert coords to index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(fast_hilbert::xy2h::<u32>(black_box(x as u32), black_box(y as u32)));
                }
            }
        })
    });

    c.bench_function("hilbert_curve: 2d hilbert coords to index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(hilbert_curve::convert_2d_to_1d(black_box(x), black_box(y), black_box(coord_length)));
                }
            }
        })
    });

    c.bench_function("hilbert_2d: 2d hilbert coords to index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(hilbert_2d::xy2h_discrete(black_box(x), black_box(y), black_box(coord_bits), black_box(hilbert_2d::Variant::Hilbert)));
                }
            }
        })
    });

    c.bench_function("hilbert_index: 2d hilbert coords to index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y]), black_box(coord_bits)));
                }
            }
        })
    });

    c.bench_function("hilbert: 2d hilbert coords to index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32])).hilbert_transform(black_box(coord_bits)));
                }
            }
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(2000);
    targets = criterion_benchmark
);
criterion_main!(benches);
