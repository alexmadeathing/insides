use criterion::{black_box, criterion_group, criterion_main, Criterion};

use insides::{Fixed, Morton, Hilbert, SpaceFillingCurve};

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

fn benchmark_coords_and_indices_2d<FIToC: Fn(usize) -> RA, FCToI: Fn(usize, usize) -> RB, RA, RB>(
    c: &mut Criterion,
    crate_name: &str,
    coord_length: usize,
    index_to_coords: FIToC,
    coords_to_index: FCToI,
) {
    let index_length = coord_length * coord_length;
    c.bench_function(
        format!("2d index to coords: {}", crate_name).as_str(),
        |b| {
            b.iter(|| {
                for i in 0..index_length {
                    index_to_coords(i);
                }
            })
        },
    );

    c.bench_function(
        format!("2d coords to index: {}", crate_name).as_str(),
        |b| {
            b.iter(|| {
                for x in 0..coord_length {
                    for y in 0..coord_length {
                        coords_to_index(x, y);
                    }
                }
            })
        },
    );
}

fn benchmark_coords_and_indices_3d<FIToC: Fn(usize), FCToI: Fn(usize, usize, usize)>(
    c: &mut Criterion,
    crate_name: &str,
    coord_length: usize,
    index_to_coords: FIToC,
    coords_to_index: FCToI,
) {
    let index_length = coord_length * coord_length * coord_length;
    c.bench_function(
        format!("3d index to coords: {}", crate_name).as_str(),
        |b| {
            b.iter(|| {
                for i in 0..index_length {
                    index_to_coords(i);
                }
            })
        },
    );

    c.bench_function(
        format!("3d coords to index: {}", crate_name).as_str(),
        |b| {
            b.iter(|| {
                for x in 0..coord_length {
                    for y in 0..coord_length {
                        for z in 0..coord_length {
                            coords_to_index(x, y, z);
                        }
                    }
                }
            })
        },
    );
}

fn benchmark_coords_and_indices_4d<FIToC: Fn(usize), FCToI: Fn(usize, usize, usize, usize)>(
    c: &mut Criterion,
    crate_name: &str,
    coord_length: usize,
    index_to_coords: FIToC,
    coords_to_index: FCToI,
) {
    let index_length = coord_length * coord_length * coord_length * coord_length;
    c.bench_function(
        format!("4d index to coords: {}", crate_name).as_str(),
        |b| {
            b.iter(|| {
                for i in 0..index_length {
                    index_to_coords(i);
                }
            })
        },
    );

    c.bench_function(
        format!("4d coords to index: {}", crate_name).as_str(),
        |b| {
            b.iter(|| {
                for x in 0..coord_length {
                    for y in 0..coord_length {
                        for z in 0..coord_length {
                            for w in 0..coord_length {
                                coords_to_index(x, y, z, w);
                            }
                        }
                    }
                }
            })
        },
    );
}

fn benchmark_morton_2d(c: &mut Criterion) {
    let coord_bits: usize = 8;
    let coord_length: usize = 1 << coord_bits;

    benchmark_coords_and_indices_2d(
        c,
        "insides (morton u16)",
        coord_length,
        |i| black_box(Morton::<Fixed<u16, 2>, 2>::from_index(black_box(i as u16)).coords()),
        |x, y| {
            black_box(
                Morton::<Fixed<u16, 2>, 2>::from_coords(black_box([x as u16, y as u16])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "insides (morton u32)",
        coord_length,
        |i| black_box(Morton::<Fixed<u32, 2>, 2>::from_index(black_box(i as u32)).coords()),
        |x, y| {
            black_box(
                Morton::<Fixed<u32, 2>, 2>::from_coords(black_box([x as u32, y as u32])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "insides (morton u64)",
        coord_length,
        |i| black_box(Morton::<Fixed<u64, 2>, 2>::from_index(black_box(i as u64)).coords()),
        |x, y| {
            black_box(
                Morton::<Fixed<u64, 2>, 2>::from_coords(black_box([x as u64, y as u64])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "insides (morton u128)",
        coord_length,
        |i| black_box(Morton::<Fixed<u128, 2>, 2>::from_index(black_box(i as u128)).coords()),
        |x, y| {
            black_box(
                Morton::<Fixed<u128, 2>, 2>::from_coords(black_box([x as u128, y as u128])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "morton (u32)",
        coord_length,
        |i| black_box(morton::deinterleave_morton(black_box(i as u32))),
        |x, y| black_box(morton::interleave_morton(black_box(x as u16), black_box(y as u16))),
    );
}

fn benchmark_morton_3d(c: &mut Criterion) {
    let coord_bits: usize = 5;
    let coord_length: usize = 1 << coord_bits;

    benchmark_coords_and_indices_3d(
        c,
        "insides (morton u16)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u16, 3>, 3>::from_index(black_box(i as u16)).coords());
            black_box(Morton::<Fixed<u16, 3>, 3>::from_index(black_box(i as u16)).coords());
        },
        |x, y, z| {
            black_box(Morton::<Fixed<u16, 3>, 3>::from_coords(black_box([x as u16, y as u16, z as u16])).index());
            black_box(Morton::<Fixed<u16, 3>, 3>::from_coords(black_box([x as u16, y as u16, z as u16])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "insides (morton u32)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u32, 3>, 3>::from_index(black_box(i as u32)).coords());
            black_box(Morton::<Fixed<u32, 3>, 3>::from_index(black_box(i as u32)).coords());
        },
        |x, y, z| {
            black_box(Morton::<Fixed<u32, 3>, 3>::from_coords(black_box([x as u32, y as u32, z as u32])).index());
            black_box(Morton::<Fixed<u32, 3>, 3>::from_coords(black_box([x as u32, y as u32, z as u32])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "insides (morton u64)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u64, 3>, 3>::from_index(black_box(i as u64)).coords());
            black_box(Morton::<Fixed<u64, 3>, 3>::from_index(black_box(i as u64)).coords());
        },
        |x, y, z| {
            black_box(Morton::<Fixed<u64, 3>, 3>::from_coords(black_box([x as u64, y as u64, z as u64])).index());
            black_box(Morton::<Fixed<u64, 3>, 3>::from_coords(black_box([x as u64, y as u64, z as u64])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "insides (morton u128)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u128, 3>, 3>::from_index(black_box(i as u128)).coords());
            black_box(Morton::<Fixed<u128, 3>, 3>::from_index(black_box(i as u128)).coords());
        },
        |x, y, z| {
            black_box(Morton::<Fixed<u128, 3>, 3>::from_coords(black_box([x as u128, y as u128, z as u128])).index());
            black_box(Morton::<Fixed<u128, 3>, 3>::from_coords(black_box([x as u128, y as u128, z as u128])).index());
        },
    );
}

fn benchmark_morton_4d(c: &mut Criterion) {
    let coord_bits: usize = 4;
    let coord_length: usize = 1 << coord_bits;

    benchmark_coords_and_indices_4d(
        c,
        "insides (morton u16)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u16, 4>, 4>::from_index(black_box(i as u16)).coords());
            black_box(Morton::<Fixed<u16, 4>, 4>::from_index(black_box(i as u16)).coords());
        },
        |x, y, z, w| {
            black_box(Morton::<Fixed<u16, 4>, 4>::from_coords(black_box([x as u16, y as u16, z as u16, w as u16])).index());
            black_box(Morton::<Fixed<u16, 4>, 4>::from_coords(black_box([x as u16, y as u16, z as u16, w as u16])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "insides (morton u32)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u32, 4>, 4>::from_index(black_box(i as u32)).coords());
            black_box(Morton::<Fixed<u32, 4>, 4>::from_index(black_box(i as u32)).coords());
        },
        |x, y, z, w| {
            black_box(Morton::<Fixed<u32, 4>, 4>::from_coords(black_box([x as u32, y as u32, z as u32, w as u32])).index());
            black_box(Morton::<Fixed<u32, 4>, 4>::from_coords(black_box([x as u32, y as u32, z as u32, w as u32])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "insides (morton u64)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u64, 4>, 4>::from_index(black_box(i as u64)).coords());
            black_box(Morton::<Fixed<u64, 4>, 4>::from_index(black_box(i as u64)).coords());
        },
        |x, y, z, w| {
            black_box(Morton::<Fixed<u64, 4>, 4>::from_coords(black_box([x as u64, y as u64, z as u64, w as u64])).index());
            black_box(Morton::<Fixed<u64, 4>, 4>::from_coords(black_box([x as u64, y as u64, z as u64, w as u64])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "insides (morton u128)",
        coord_length,
        |i| {
            black_box(Morton::<Fixed<u128, 4>, 4>::from_index(black_box(i as u128)).coords());
            black_box(Morton::<Fixed<u128, 4>, 4>::from_index(black_box(i as u128)).coords());
        },
        |x, y, z, w| {
            black_box(Morton::<Fixed<u128, 4>, 4>::from_coords(black_box([x as u128, y as u128, z as u128, w as u128])).index());
            black_box(Morton::<Fixed<u128, 4>, 4>::from_coords(black_box([x as u128, y as u128, z as u128, w as u128])).index());
        },
    );
}

fn benchmark_hilbert_2d(c: &mut Criterion) {
    let coord_bits: usize = 8;
    let coord_length: usize = 1 << coord_bits;

    benchmark_coords_and_indices_2d(
        c,
        "insides (hilbert u16)",
        coord_length,
        |i| black_box(Hilbert::<Fixed<u16, 2>, 2>::from_index(black_box(i as u16)).coords()),
        |x, y| {
            black_box(
                Hilbert::<Fixed<u16, 2>, 2>::from_coords(black_box([x as u16, y as u16])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "insides (hilbert u32)",
        coord_length,
        |i| black_box(Hilbert::<Fixed<u32, 2>, 2>::from_index(black_box(i as u32)).coords()),
        |x, y| {
            black_box(
                Hilbert::<Fixed<u32, 2>, 2>::from_coords(black_box([x as u32, y as u32])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "insides (hilbert u64)",
        coord_length,
        |i| black_box(Hilbert::<Fixed<u64, 2>, 2>::from_index(black_box(i as u64)).coords()),
        |x, y| {
            black_box(
                Hilbert::<Fixed<u64, 2>, 2>::from_coords(black_box([x as u64, y as u64])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "insides (hilbert u128)",
        coord_length,
        |i| black_box(Hilbert::<Fixed<u128, 2>, 2>::from_index(black_box(i as u128)).coords()),
        |x, y| {
            black_box(
                Hilbert::<Fixed<u128, 2>, 2>::from_coords(black_box([x as u128, y as u128])).index(),
            )
        },
    );

    benchmark_coords_and_indices_2d(
        c,
        "fast_hilbert (u16)",
        coord_length,
        |i| black_box(fast_hilbert::h2xy::<u8>(black_box(i as u16), black_box(coord_bits as u8))),
        |x, y| black_box(fast_hilbert::xy2h::<u8>(black_box(x as u8), black_box(y as u8), black_box(coord_bits as u8))),
    );

    benchmark_coords_and_indices_2d(
        c,
        "fast_hilbert (u32)",
        coord_length,
        |i| black_box(fast_hilbert::h2xy::<u16>(black_box(i as u32), black_box(coord_bits as u8))),
        |x, y| black_box(fast_hilbert::xy2h::<u16>(black_box(x as u16), black_box(y as u16), black_box(coord_bits as u8))),
    );

    benchmark_coords_and_indices_2d(
        c,
        "fast_hilbert (u64)",
        coord_length,
        |i| black_box(fast_hilbert::h2xy::<u32>(black_box(i as u64), black_box(coord_bits as u8))),
        |x, y| black_box(fast_hilbert::xy2h::<u32>(black_box(x as u32), black_box(y as u32), black_box(coord_bits as u8))),
    );

    benchmark_coords_and_indices_2d(
        c,
        "fast_hilbert (u128)",
        coord_length,
        |i| black_box(fast_hilbert::h2xy::<u64>(black_box(i as u128), black_box(coord_bits as u8))),
        |x, y| black_box(fast_hilbert::xy2h::<u64>(black_box(x as u64), black_box(y as u64), black_box(coord_bits as u8))),
    );

    benchmark_coords_and_indices_2d(
        c,
        "hilbert_curve (usize)",
        coord_length,
        |i| black_box(hilbert_curve::convert_1d_to_2d(black_box(i), black_box(coord_length))),
        |x, y| black_box(hilbert_curve::convert_2d_to_1d(black_box(x), black_box(y), black_box(coord_length))),
    );

    benchmark_coords_and_indices_2d(
        c,
        "hilbert_2d (usize)",
        coord_length,
        |i| black_box(hilbert_2d::h2xy_discrete(black_box(i), black_box(coord_bits), black_box(hilbert_2d::Variant::Hilbert))),
        |x, y| black_box(hilbert_2d::xy2h_discrete(black_box(x), black_box(y), black_box(coord_bits), black_box(hilbert_2d::Variant::Hilbert))),
    );

    benchmark_coords_and_indices_2d(
        c,
        "hilbert_index (usize)",
        coord_length,
        |i| black_box(hilbert_index::FromHilbertIndex::<2>::from_hilbert_index(black_box(&i), black_box(coord_bits))),
        |x, y| black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y]), black_box(coord_bits))),
    );

    // This is causing compile errors
    /*
    benchmark_coords_and_indices_2d(
        c,
        "hilbert (BigUInt)",
        coord_length,
        |i| black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(2))),
        |x, y| black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32])).hilbert_transform(black_box(coord_bits))),
    );
    */
}

fn benchmark_hilbert_3d(c: &mut Criterion) {
    let coord_bits: usize = 5;
    let coord_length: usize = 1 << coord_bits;

    benchmark_coords_and_indices_3d(
        c,
        "insides (hilbert u16)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u16, 3>, 3>::from_index(black_box(i as u16)).coords());
            black_box(Hilbert::<Fixed<u16, 3>, 3>::from_index(black_box(i as u16)).coords());
        },
        |x, y, z| {
            black_box(Hilbert::<Fixed<u16, 3>, 3>::from_coords(black_box([x as u16, y as u16, z as u16])).index());
            black_box(Hilbert::<Fixed<u16, 3>, 3>::from_coords(black_box([x as u16, y as u16, z as u16])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "insides (hilbert u32)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u32, 3>, 3>::from_index(black_box(i as u32)).coords());
            black_box(Hilbert::<Fixed<u32, 3>, 3>::from_index(black_box(i as u32)).coords());
        },
        |x, y, z| {
            black_box(Hilbert::<Fixed<u32, 3>, 3>::from_coords(black_box([x as u32, y as u32, z as u32])).index());
            black_box(Hilbert::<Fixed<u32, 3>, 3>::from_coords(black_box([x as u32, y as u32, z as u32])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "insides (hilbert u64)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u64, 3>, 3>::from_index(black_box(i as u64)).coords());
            black_box(Hilbert::<Fixed<u64, 3>, 3>::from_index(black_box(i as u64)).coords());
        },
        |x, y, z| {
            black_box(Hilbert::<Fixed<u64, 3>, 3>::from_coords(black_box([x as u64, y as u64, z as u64])).index());
            black_box(Hilbert::<Fixed<u64, 3>, 3>::from_coords(black_box([x as u64, y as u64, z as u64])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "insides (hilbert u128)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u128, 3>, 3>::from_index(black_box(i as u128)).coords());
            black_box(Hilbert::<Fixed<u128, 3>, 3>::from_index(black_box(i as u128)).coords());
        },
        |x, y, z| {
            black_box(Hilbert::<Fixed<u128, 3>, 3>::from_coords(black_box([x as u128, y as u128, z as u128])).index());
            black_box(Hilbert::<Fixed<u128, 3>, 3>::from_coords(black_box([x as u128, y as u128, z as u128])).index());
        },
    );

    benchmark_coords_and_indices_3d(
        c,
        "hilbert_index (usize)",
        coord_length,
        |i| {
            black_box(hilbert_index::FromHilbertIndex::<3>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
            black_box(hilbert_index::FromHilbertIndex::<3>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
        },
        |x, y, z| {
            black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z]), black_box(coord_bits)));
            black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z]), black_box(coord_bits)));
        },
    );

    // This is causing compile errors
    /*
    benchmark_coords_and_indices_3d(
        c,
        "hilbert (BigUInt)",
        coord_length,
        |i| {
            black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(3)));
            black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(3)));
        },
        |x, y, z| {
            black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32])).hilbert_transform(black_box(coord_bits)));
            black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32])).hilbert_transform(black_box(coord_bits)));
        },
    );
    */
}

fn benchmark_hilbert_4d(c: &mut Criterion) {
    let coord_bits: usize = 4;
    let coord_length: usize = 1 << coord_bits;

    benchmark_coords_and_indices_4d(
        c,
        "insides (hilbert u16)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u16, 4>, 4>::from_index(black_box(i as u16)).coords());
            black_box(Hilbert::<Fixed<u16, 4>, 4>::from_index(black_box(i as u16)).coords());
        },
        |x, y, z, w| {
            black_box(Hilbert::<Fixed<u16, 4>, 4>::from_coords(black_box([x as u16, y as u16, z as u16, w as u16])).index());
            black_box(Hilbert::<Fixed<u16, 4>, 4>::from_coords(black_box([x as u16, y as u16, z as u16, w as u16])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "insides (hilbert u32)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u32, 4>, 4>::from_index(black_box(i as u32)).coords());
            black_box(Hilbert::<Fixed<u32, 4>, 4>::from_index(black_box(i as u32)).coords());
        },
        |x, y, z, w| {
            black_box(Hilbert::<Fixed<u32, 4>, 4>::from_coords(black_box([x as u32, y as u32, z as u32, w as u32])).index());
            black_box(Hilbert::<Fixed<u32, 4>, 4>::from_coords(black_box([x as u32, y as u32, z as u32, w as u32])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "insides (hilbert u64)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u64, 4>, 4>::from_index(black_box(i as u64)).coords());
            black_box(Hilbert::<Fixed<u64, 4>, 4>::from_index(black_box(i as u64)).coords());
        },
        |x, y, z, w| {
            black_box(Hilbert::<Fixed<u64, 4>, 4>::from_coords(black_box([x as u64, y as u64, z as u64, w as u64])).index());
            black_box(Hilbert::<Fixed<u64, 4>, 4>::from_coords(black_box([x as u64, y as u64, z as u64, w as u64])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "insides (hilbert u128)",
        coord_length,
        |i| {
            black_box(Hilbert::<Fixed<u128, 4>, 4>::from_index(black_box(i as u128)).coords());
            black_box(Hilbert::<Fixed<u128, 4>, 4>::from_index(black_box(i as u128)).coords());
        },
        |x, y, z, w| {
            black_box(Hilbert::<Fixed<u128, 4>, 4>::from_coords(black_box([x as u128, y as u128, z as u128, w as u128])).index());
            black_box(Hilbert::<Fixed<u128, 4>, 4>::from_coords(black_box([x as u128, y as u128, z as u128, w as u128])).index());
        },
    );

    benchmark_coords_and_indices_4d(
        c,
        "hilbert_index (usize)",
        coord_length,
        |i| {
            black_box(hilbert_index::FromHilbertIndex::<4>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
            black_box(hilbert_index::FromHilbertIndex::<4>::from_hilbert_index(black_box(&i), black_box(coord_bits)));
        },
        |x, y, z, w| {
            black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z, w]), black_box(coord_bits)));
            black_box(hilbert_index::ToHilbertIndex::to_hilbert_index(&black_box([x, y, z, w]), black_box(coord_bits)));
        },
    );

    // This is causing compile errors
    /*
    benchmark_coords_and_indices_4d(
        c,
        "hilbert (BigUInt)",
        coord_length,
        |i| {
            black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(4)));
            black_box(hilbert::Point::new_from_hilbert_index(0, &black_box(num::BigUint::from(i)), black_box(coord_bits), black_box(4)));
        },
        |x, y, z, w| {
            black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32, w as u32])).hilbert_transform(black_box(coord_bits)));
            black_box(black_box(hilbert::Point::new(0, &[x as u32, y as u32, z as u32, w as u32])).hilbert_transform(black_box(coord_bits)));
        },
    );
    */
}

pub fn criterion_benchmark(c: &mut Criterion) {
    benchmark_morton_2d(c);
    benchmark_morton_3d(c);
    benchmark_morton_4d(c);

    benchmark_hilbert_2d(c);
    benchmark_hilbert_3d(c);
    benchmark_hilbert_4d(c);
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(1000);
    targets = criterion_benchmark
);
criterion_main!(benches);
