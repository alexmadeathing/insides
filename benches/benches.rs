use criterion::{black_box, criterion_group, criterion_main, Criterion};

use insides::{Hilbert, Expand, SpaceFillingCurve};

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

fn benchmark_2d(c: &mut Criterion) {
    let coord_bits: usize = 8;
    let coord_length: usize = 1 << coord_bits;
    let index_length: usize = coord_length * coord_length;

    c.bench_function("2d hilbert index to coords: insides (u8)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(Hilbert::<Expand<u8, 2>, 2>::from_index(black_box(i as u16)).coords());
            }
        })
    });

    c.bench_function("2d hilbert coords to index: insides (u8)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(Hilbert::<Expand<u8, 2>, 2>::from_coords(black_box([x as u8, y as u8])).index());
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: insides (u16)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(Hilbert::<Expand<u16, 2>, 2>::from_index(black_box(i as u32)).coords());
            }
        })
    });

    c.bench_function("2d hilbert coords to index: insides (u16)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(Hilbert::<Expand<u16, 2>, 2>::from_coords(black_box([x as u16, y as u16])).index());
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: insides (u32)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(Hilbert::<Expand<u32, 2>, 2>::from_index(black_box(i as u64)).coords());
            }
        })
    });

    c.bench_function("2d hilbert coords to index: insides (u32)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(Hilbert::<Expand<u32, 2>, 2>::from_coords(black_box([x as u32, y as u32])).index());
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: insides (u64)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(Hilbert::<Expand<u64, 2>, 2>::from_index(black_box(i as u128)).coords());
            }
        })
    });

    c.bench_function("2d hilbert coords to index: insides (u64)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(Hilbert::<Expand<u64, 2>, 2>::from_coords(black_box([x as u64, y as u64])).index());
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: fast_hilbert (u8)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(fast_hilbert::h2xy::<u8>(black_box(i as u16)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: fast_hilbert (u8)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(fast_hilbert::xy2h::<u8>(black_box(x as u8), black_box(y as u8)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: fast_hilbert (u16)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(fast_hilbert::h2xy::<u16>(black_box(i as u32)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: fast_hilbert (u16)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(fast_hilbert::xy2h::<u16>(black_box(x as u16), black_box(y as u16)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: fast_hilbert (u32)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(fast_hilbert::h2xy::<u32>(black_box(i as u64)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: fast_hilbert (u32)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(fast_hilbert::xy2h::<u32>(black_box(x as u32), black_box(y as u32)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: fast_hilbert (u64)", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(fast_hilbert::h2xy::<u64>(black_box(i as u128)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: fast_hilbert (u64)", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(fast_hilbert::xy2h::<u64>(black_box(x as u64), black_box(y as u64)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: hilbert_curve", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_curve::convert_1d_to_2d(black_box(i), black_box(coord_length)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: hilbert_curve", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(hilbert_curve::convert_2d_to_1d(black_box(x), black_box(y), black_box(coord_length)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: hilbert_2d", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_2d::h2xy_discrete(black_box(i), black_box(coord_bits), black_box(hilbert_2d::Variant::Hilbert)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: hilbert_2d", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(hilbert_2d::xy2h_discrete(black_box(x), black_box(y), black_box(coord_bits), black_box(hilbert_2d::Variant::Hilbert)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: hilbert_index", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_index::FromHilbertIndex::<2>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: hilbert_index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y]), black_box(coord_bits)));
                }
            }
        })
    });

    c.bench_function("2d hilbert index to coords: hilbert", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(2)));
            }
        })
    });

    c.bench_function("2d hilbert coords to index: hilbert", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32])).hilbert_transform(black_box(coord_bits)));
                }
            }
        })
    });
}

fn benchmark_3d(c: &mut Criterion) {
    let coord_bits: usize = 5;
    let coord_length: usize = 1 << coord_bits;
    let index_length: usize = coord_length * coord_length * coord_length;

    c.bench_function("3d hilbert index to coords: insides", |b| {
        b.iter(|| {
            for i in 0..index_length {
                // Perform twice to balance with 2D and 4D
                black_box(Hilbert::<Expand<u16, 3>, 3>::from_index(black_box(i as u64)).coords());
                black_box(Hilbert::<Expand<u16, 3>, 3>::from_index(black_box(i as u64)).coords());
            }
        })
    });

    c.bench_function("3d hilbert coords to index: insides", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    for z in 0..coord_length {
                        // Perform twice to balance with 2D and 4D
                        black_box(Hilbert::<Expand<u16, 3>, 3>::from_coords(black_box([x as u16, y as u16, z as u16])).index());
                        black_box(Hilbert::<Expand<u16, 3>, 3>::from_coords(black_box([x as u16, y as u16, z as u16])).index());
                    }
                }
            }
        })
    });

    c.bench_function("3d hilbert index to coords: hilbert_index", |b| {
        b.iter(|| {
            for i in 0..index_length {
                // Perform twice to balance with 2D and 4D
                black_box(hilbert_index::FromHilbertIndex::<3>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
                black_box(hilbert_index::FromHilbertIndex::<3>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
            }
        })
    });

    c.bench_function("3d hilbert coords to index: hilbert_index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    for z in 0..coord_length {
                        // Perform twice to balance with 2D and 4D
                        black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z]), black_box(coord_bits)));
                        black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z]), black_box(coord_bits)));
                    }
                }
            }
        })
    });

    c.bench_function("3d hilbert index to coords: hilbert", |b| {
        b.iter(|| {
            for i in 0..index_length {
                // Perform twice to balance with 2D and 4D
                black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(3)));
                black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(3)));
            }
        })
    });

    c.bench_function("3d hilbert coords to index: hilbert", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    for z in 0..coord_length {
                        // Perform twice to balance with 2D and 4D
                        black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32])).hilbert_transform(black_box(coord_bits)));
                        black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32])).hilbert_transform(black_box(coord_bits)));
                    }
                }
            }
        })
    });
}

fn benchmark_4d(c: &mut Criterion) {
    let coord_bits: usize = 4;
    let coord_length: usize = 1 << coord_bits;
    let index_length: usize = coord_length * coord_length * coord_length * coord_length;

    c.bench_function("4d hilbert index to coords: insides", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(Hilbert::<Expand<u16, 4>, 4>::from_index(black_box(i as u64)).coords());
            }
        })
    });

    c.bench_function("4d hilbert coords to index: insides", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    for z in 0..coord_length {
                        for w in 0..coord_length {
                            black_box(Hilbert::<Expand<u16, 4>, 4>::from_coords(black_box([x as u16, y as u16, z as u16, w as u16])).index());
                        }
                    }
                }
            }
        })
    });

    c.bench_function("4d hilbert index to coords: hilbert_index", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert_index::FromHilbertIndex::<4>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
            }
        })
    });

    c.bench_function("4d hilbert coords to index: hilbert_index", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    for z in 0..coord_length {
                        for w in 0..coord_length {
                            black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z, w]), black_box(coord_bits)));
                        }
                    }
                }
            }
        })
    });

    c.bench_function("4d hilbert index to coords: hilbert", |b| {
        b.iter(|| {
            for i in 0..index_length {
                black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(4)));
            }
        })
    });

    c.bench_function("4d hilbert coords to index: hilbert", |b| {
        b.iter(|| {
            for x in 0..coord_length {
                for y in 0..coord_length {
                    for z in 0..coord_length {
                        for w in 0..coord_length {
                            black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32, w as u32])).hilbert_transform(black_box(coord_bits)));
                        }
                    }
                }
            }
        })
    });
}

pub fn criterion_benchmark(c: &mut Criterion) {
    benchmark_2d(c);
    benchmark_3d(c);
    benchmark_4d(c);
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(1000);
    targets = criterion_benchmark
);
criterion_main!(benches);
