use dilate::*;

use super::{CurveCoord, CurveIndex, SpaceFillingCurve};

use crate::{internal::NumTraits, Morton};

mod lut_small {
    use super::NumTraits;

    #[cfg(feature = "lut_small_d2")]
    mod d2 {
        pub const FLIP_TX_INTO_SHL_LUT: [[u8; 4]; 2] = [
            [00, 00, 00, 03, ],
            [00, 00, 00, 03, ],
        ];

        pub const SHR_INTO_GRAY_INV_LUT: [[u8; 4]; 2] = [
            [00, 01, 03, 02, ],
            [00, 03, 01, 02, ],
        ];

        pub const GRAY_INTO_SHL_LUT: [[u8; 4]; 2] = [
            [00, 01, 03, 02, ],
            [00, 02, 03, 01, ],
        ];

        pub const ROTATION_LUT: [[u8; 4]; 2] = [
            [01, 00, 00, 01, ],
            [00, 01, 01, 00, ],
        ];
    }

    #[cfg(feature = "lut_small_d3")]
    mod d3 {
        pub const FLIP_TX_INTO_SHL_LUT: [[u8; 8]; 3] = [
            [00, 00, 00, 03, 03, 06, 06, 05, ],
            [00, 00, 00, 06, 06, 05, 05, 03, ],
            [00, 00, 00, 05, 05, 03, 03, 06, ],
        ];

        pub const SHR_INTO_GRAY_INV_LUT: [[u8; 8]; 3] = [
            [00, 01, 03, 02, 07, 06, 04, 05, ],
            [00, 07, 01, 06, 03, 04, 02, 05, ],
            [00, 03, 07, 04, 01, 02, 06, 05, ],
        ];

        pub const GRAY_INTO_SHL_LUT: [[u8; 8]; 3] = [
            [00, 01, 03, 02, 06, 07, 05, 04, ],
            [00, 02, 06, 04, 05, 07, 03, 01, ],
            [00, 04, 05, 01, 03, 07, 06, 02, ],
        ];

        pub const ROTATION_LUT: [[u8; 8]; 3] = [
            [01, 02, 02, 00, 00, 02, 02, 01, ],
            [02, 00, 00, 01, 01, 00, 00, 02, ],
            [00, 01, 01, 02, 02, 01, 01, 00, ],
        ];
    }

    #[cfg(feature = "lut_small_d4")]
    mod d4 {
        pub const FLIP_TX_INTO_SHL_LUT: [[u8; 16]; 4] = [
            [00, 00, 00, 03, 03, 06, 06, 05, 05, 12, 12, 15, 15, 10, 10, 09, ],
            [00, 00, 00, 06, 06, 12, 12, 10, 10, 09, 09, 15, 15, 05, 05, 03, ],
            [00, 00, 00, 12, 12, 09, 09, 05, 05, 03, 03, 15, 15, 10, 10, 06, ],
            [00, 00, 00, 09, 09, 03, 03, 10, 10, 06, 06, 15, 15, 05, 05, 12, ],
        ];

        pub const SHR_INTO_GRAY_INV_LUT: [[u8; 16]; 4] = [
            [00, 01, 03, 02, 07, 06, 04, 05, 15, 14, 12, 13, 08, 09, 11, 10, ],
            [00, 15, 01, 14, 03, 12, 02, 13, 07, 08, 06, 09, 04, 11, 05, 10, ],
            [00, 07, 15, 08, 01, 06, 14, 09, 03, 04, 12, 11, 02, 05, 13, 10, ],
            [00, 03, 07, 04, 15, 12, 08, 11, 01, 02, 06, 05, 14, 13, 09, 10, ],
        ];

        pub const GRAY_INTO_SHL_LUT: [[u8; 16]; 4] = [
            [00, 01, 03, 02, 06, 07, 05, 04, 12, 13, 15, 14, 10, 11, 09, 08, ],
            [00, 02, 06, 04, 12, 14, 10, 08, 09, 11, 15, 13, 05, 07, 03, 01, ],
            [00, 04, 12, 08, 09, 13, 05, 01, 03, 07, 15, 11, 10, 14, 06, 02, ],
            [00, 08, 09, 01, 03, 11, 10, 02, 06, 14, 15, 07, 05, 13, 12, 04, ],
        ];

        pub const ROTATION_LUT: [[u8; 16]; 4] = [
            [01, 02, 02, 03, 03, 02, 02, 00, 00, 02, 02, 03, 03, 02, 02, 01, ],
            [02, 03, 03, 00, 00, 03, 03, 01, 01, 03, 03, 00, 00, 03, 03, 02, ],
            [03, 00, 00, 01, 01, 00, 00, 02, 02, 00, 00, 01, 01, 00, 00, 03, ],
            [00, 01, 01, 02, 02, 01, 01, 03, 03, 01, 01, 02, 02, 01, 01, 00, ],
        ];
    }

    #[inline(always)]
    pub fn has_lut<const D: usize>() -> bool {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => true,
            #[cfg(feature = "lut_small_d3")]
            3 => true,
            #[cfg(feature = "lut_small_d4")]
            4 => true,
            _ => false,
        }
    }

    #[inline(always)]
    fn flip_tx_into_shl_lut<T: NumTraits, const D: usize>(_index: T, _rotation: usize) -> T {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => T::from_u8(d2::FLIP_TX_INTO_SHL_LUT[_rotation][_index.to_usize()]),
            #[cfg(feature = "lut_small_d3")]
            3 => T::from_u8(d3::FLIP_TX_INTO_SHL_LUT[_rotation][_index.to_usize()]),
            #[cfg(feature = "lut_small_d4")]
            4 => T::from_u8(d4::FLIP_TX_INTO_SHL_LUT[_rotation][_index.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn shr_into_gray_inv_lut<T: NumTraits, const D: usize>(_index: T, _rotation: usize) -> T {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => T::from_u8(d2::SHR_INTO_GRAY_INV_LUT[_rotation][_index.to_usize()]),
            #[cfg(feature = "lut_small_d3")]
            3 => T::from_u8(d3::SHR_INTO_GRAY_INV_LUT[_rotation][_index.to_usize()]),
            #[cfg(feature = "lut_small_d4")]
            4 => T::from_u8(d4::SHR_INTO_GRAY_INV_LUT[_rotation][_index.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn gray_into_shl_lut<T: NumTraits, const D: usize>(_index: T, _rotation: usize) -> T {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => T::from_u8(d2::GRAY_INTO_SHL_LUT[_rotation][_index.to_usize()]),
            #[cfg(feature = "lut_small_d3")]
            3 => T::from_u8(d3::GRAY_INTO_SHL_LUT[_rotation][_index.to_usize()]),
            #[cfg(feature = "lut_small_d4")]
            4 => T::from_u8(d4::GRAY_INTO_SHL_LUT[_rotation][_index.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn rotation_lut<T: NumTraits, const D: usize>(_index: T, _rotation: usize) -> usize {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => d2::ROTATION_LUT[_rotation][_index.to_usize()] as usize,
            #[cfg(feature = "lut_small_d3")]
            3 => d3::ROTATION_LUT[_rotation][_index.to_usize()] as usize,
            #[cfg(feature = "lut_small_d4")]
            4 => d4::ROTATION_LUT[_rotation][_index.to_usize()] as usize,
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn walk_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        order: usize,
        mut f: F,
    ) -> (usize, T)
    where
        F: FnMut(T, usize, usize, T) -> T,
    {
        let lower_mask = T::one().shl(D).sub(T::one());

        // We iterate in reverse because higher orders affect the rotation of lower orders
        let mut flip = T::zero();
        let mut rotation = (D - 1) - (order % D);
        for shift in (0..order).rev().map(|i| i * D) {
            // Extract a portion of the index
            let index_partial = index.shr(shift).bit_and(lower_mask);

            // Context specific code generates the hilbert partial
            let hilbert_partial = f(index_partial, shift, rotation, flip);

            // Calculate the next axis flip flags and rotation amounts resulting from this order
            flip = flip.bit_xor(flip_tx_into_shl_lut::<_, D>(hilbert_partial, rotation));
            rotation = rotation_lut::<_, D>(hilbert_partial, rotation);
        }
        (rotation, flip)
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        debug_assert!(has_lut::<D>());

        let order = (T::bits() - morton_index.lz() + D - 1) / D;

        let mut hilbert_index = T::zero();
        walk_transforms::<T, _, D>(
            morton_index,
            order,
            |morton_partial, shift, rotation, flip| {
                let hilbert_partial = shr_into_gray_inv_lut::<_, D>(flip.bit_xor(morton_partial), rotation);
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        debug_assert!(has_lut::<D>());

        let order = (T::bits() - hilbert_index.lz() + D - 1) / D;

        let mut morton_index = T::zero();
        walk_transforms::<T, _, D>(
            hilbert_index,
            order,
            |hilbert_partial, shift, rotation, flip| {
                let morton_partial = flip.bit_xor(gray_into_shl_lut::<_, D>(hilbert_partial, rotation));
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
    }
}

mod lut_large {
    use super::NumTraits;

    #[cfg(feature = "lut_large_d2")]
    mod d2 {
        pub const STRIDE: usize = 4;
        pub const TRANSFORM_LUT: [[u8; 256]; 4] = include!("lut/hilbert_d2_transform_lut.in");
        pub const MORTON_LUT: [[u8; 256]; 4] = include!("lut/hilbert_d2_morton_lut.in");
        pub const HILBERT_LUT: [[u8; 256]; 4] = include!("lut/hilbert_d2_hilbert_lut.in");
    }
    
    #[cfg(feature = "lut_large_d3")]
    mod d3 {
        pub const STRIDE: usize = 2;
        pub const TRANSFORM_LUT: [[u8; 64]; 12] = include!("lut/hilbert_d3_transform_lut.in");
        pub const MORTON_LUT: [[u8; 64]; 12] = include!("lut/hilbert_d3_morton_lut.in");
        pub const HILBERT_LUT: [[u8; 64]; 12] = include!("lut/hilbert_d3_hilbert_lut.in");
    }
    
    #[cfg(feature = "lut_large_d4")]
    mod d4 {
        pub const STRIDE: usize = 2;
        pub const TRANSFORM_LUT: [[u8; 256]; 32] = include!("lut/hilbert_d4_transform_lut.in");
        pub const MORTON_LUT: [[u8; 256]; 32] = include!("lut/hilbert_d4_morton_lut.in");
        pub const HILBERT_LUT: [[u8; 256]; 32] = include!("lut/hilbert_d4_hilbert_lut.in");
    }

    #[inline(always)]
    pub fn has_lut<const D: usize>() -> bool {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => true,
            #[cfg(feature = "lut_large_d3")]
            3 => true,
            #[cfg(feature = "lut_large_d4")]
            4 => true,
            _ => false,
        }
    }

    #[inline(always)]
    fn stride<const D: usize>() -> usize {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => d2::STRIDE,
            #[cfg(feature = "lut_large_d3")]
            3 => d3::STRIDE,
            #[cfg(feature = "lut_large_d4")]
            4 => d4::STRIDE,
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn lower_mask<T: NumTraits, const D: usize>() -> T {
        (0..(stride::<D>() * D)).fold(T::zero(), |m, _| m.shl(1).bit_or(T::one()))
    }

    #[inline(always)]
    fn combined_transform_lut<T: NumTraits, const D: usize>(
        _transform: usize,
        _hilbert_partial: T,
    ) -> usize {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => d2::TRANSFORM_LUT[_transform][_hilbert_partial.to_usize()] as usize,
            #[cfg(feature = "lut_large_d3")]
            3 => d3::TRANSFORM_LUT[_transform][_hilbert_partial.to_usize()] as usize,
            #[cfg(feature = "lut_large_d4")]
            4 => d4::TRANSFORM_LUT[_transform][_hilbert_partial.to_usize()] as usize,
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn combined_morton_lut<T: NumTraits, const D: usize>(
        _transform: usize,
        _hilbert_partial: T,
    ) -> T {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => T::from_u8(d2::MORTON_LUT[_transform][_hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d3")]
            3 => T::from_u8(d3::MORTON_LUT[_transform][_hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d4")]
            4 => T::from_u8(d4::MORTON_LUT[_transform][_hilbert_partial.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn combined_hilbert_lut<T: NumTraits, const D: usize>(
        _transform: usize,
        _hilbert_partial: T,
    ) -> T {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => T::from_u8(d2::HILBERT_LUT[_transform][_hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d3")]
            3 => T::from_u8(d3::HILBERT_LUT[_transform][_hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d4")]
            4 => T::from_u8(d4::HILBERT_LUT[_transform][_hilbert_partial.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn walk_combined_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        order: usize,
        mut f: F,
    ) -> T
    where
        F: FnMut(usize, T) -> (T, T),
    {
        let stride_bits = stride::<D>() * D;

        let output = T::zero();
        let transform = order * stride::<D>() % D;
        let state = (output, transform);

        (0..order).rev().map(|i| i * stride_bits).fold(state, |state, shift| {
            let (hilbert_partial, out_partial) = f(state.1, index.shr(shift).bit_and(lower_mask::<T, D>()));
            (state.0.bit_or(out_partial.shl(shift)), combined_transform_lut::<T, D>(state.1, hilbert_partial))
        }).0
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        debug_assert!(has_lut::<D>());

        let stride_bits = stride::<D>() * D;
        let order = ((T::bits() + stride_bits - 1) - morton_index.lz()) / stride_bits;

        walk_combined_transforms::<T, _, D>(
            morton_index,
            order,
            |transform, morton_partial| {
                let hilbert_partial = combined_hilbert_lut::<T, D>(transform, morton_partial);
                (hilbert_partial, hilbert_partial)
            },
        )
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        debug_assert!(has_lut::<D>());

        let stride_bits = stride::<D>() * D;
        let order = ((T::bits() + stride_bits - 1) - hilbert_index.lz()) / stride_bits;

        walk_combined_transforms::<T, _, D>(
            hilbert_index,
            order,
            |transform, hilbert_partial| {
                let morton_partial = combined_morton_lut::<T, D>(transform, hilbert_partial);
                (hilbert_partial, morton_partial)
            },
        )
    }
}

mod explicit {
    use super::NumTraits;
    
    #[inline(always)]
    pub fn gray<T: NumTraits>(index: T) -> T {
        index.bit_xor(index.shr(1))
    }
    
    #[inline(always)]
    pub fn gray_inverse<T: NumTraits, const D: usize>(index: T) -> T {
        // Since we know that only a small subset of bits of the integer are used,
        // we can perform a minimal gray code inverse rather than a complete one
        (1..D).fold(index, |acc, shift| acc.bit_xor(index.shr(shift)))
    }
    
    #[inline(always)]
    pub fn shl_cyclic<T: NumTraits, const D: usize>(index: T, amount: usize) -> T {
        let lower_mask: T = <T as NumTraits>::one().shl(D).sub(NumTraits::one());
        lower_mask.bit_and(index.shl(amount).bit_or(index.shr(D - amount)))
    }
    
    #[inline(always)]
    pub fn shr_cyclic<T: NumTraits, const D: usize>(index: T, amount: usize) -> T {
        let lower_mask: T = <T as NumTraits>::one().shl(D).sub(NumTraits::one());
        lower_mask.bit_and(index.shr(amount).bit_or(index.shl(D - amount)))
    }
    
    #[inline(always)]
    pub fn rotation_transform(index: usize) -> usize {
        if index > 0 {
            index.wrapping_add((index & 0x1).wrapping_sub(1)).trailing_ones() as usize
        } else {
            0
        }
    }
    
    #[inline(always)]
    pub fn flip_transform<T: NumTraits>(index: T) -> T {
        if index == T::zero() {
            T::zero()
        } else {
            gray(index.sub(T::one()).bit_and(T::bit_not(T::one())))
        }
    }

    #[inline(always)]
    pub fn next_flip<T: NumTraits, const D: usize>(index: T, flip: T, rotation: usize) -> T {
        flip.bit_xor(shl_cyclic::<_, D>(flip_transform(index), rotation))
    }

    #[inline(always)]
    pub fn next_rotation<T: NumTraits, const D: usize>(index: T, rotation: usize) -> usize {
        (rotation_transform(index.to_usize()) + rotation + 1) % D
    }

    #[inline(always)]
    pub fn walk_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        min_order: usize,
        mut f: F,
    ) -> (usize, T)
    where
        F: FnMut(T, usize, usize, T) -> T,
    {
        let lower_mask = T::one().shl(D).sub(T::one());
        let order = (T::bits() - index.lz() + D - 1) / D;

        let mut flip = T::zero();
        let mut rotation = (D - 1) - (order % D);
        for shift in (min_order..order).rev().map(|i| i * D) {
            let index_partial = index.shr(shift).bit_and(lower_mask);
            let hilbert_partial = f(index_partial, shift, rotation, flip);
            flip = next_flip::<T, D>(hilbert_partial, flip, rotation);
            rotation = next_rotation::<T, D>(hilbert_partial, rotation);
        }
        (rotation, flip)
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        let mut hilbert_index = T::zero();
        walk_transforms::<_, _, D>(
            morton_index,
            0,
            |morton_partial, shift, rotation, flip| {
                let hilbert_partial = gray_inverse::<_, D>(shr_cyclic::<_, D>(
                    flip.bit_xor(morton_partial),
                    rotation,
                ));
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        let mut morton_index = T::zero();
        walk_transforms::<_, _, D>(
            hilbert_index,
            0,
            |hilbert_partial, shift, rotation, flip| {
                let morton_partial = flip.bit_xor(shl_cyclic::<_, D>(
                    gray(hilbert_partial),
                    rotation,
                ));
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
    }
}

#[inline(always)]
fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
    if lut_large::has_lut::<D>() {
        return lut_large::morton_to_hilbert::<_, D>(morton_index);
    }

    if lut_small::has_lut::<D>() {
        return lut_small::morton_to_hilbert::<_, D>(morton_index);
    }

    explicit::morton_to_hilbert::<_, D>(morton_index)
}

#[inline(always)]
fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
    if lut_large::has_lut::<D>() {
        return lut_large::hilbert_to_morton::<_, D>(hilbert_index);
    }

    if lut_small::has_lut::<D>() {
        return lut_small::hilbert_to_morton::<_, D>(hilbert_index);
    }

    explicit::hilbert_to_morton::<_, D>(hilbert_index)
}

/// A Hilbert encoded space filling curve implementation
///
/// [Hilbert encoding](https://en.wikipedia.org/wiki/Hilbert_curve)
/// is a fractal space filling algorithm which maps a multidimensional
/// set of coordinates to one dimension and vice versa.
///
/// Where the Morton curve features better CPU performance, the Hilbert curve
/// offers excellent spatial coherence, making it an ideal choice for memory
/// access bound applications.
/// 
/// # Self Similarity
/// The Hilbert implementation featured in this crate is self-similar
/// in that lower magnitude indices will always reference the same
/// coordinates regardless how deep the curve is iterated. There is no
/// "depth" or "level" parameter.
/// 
/// # Lookup Tables
/// Two lookup table (LUT) options are provided, which greatly improve
/// performance (for dimensions 2, 3, and 4). These may be toggled on or off
/// using the cargo features:
/// 
/// `lut_small_d2`, `lut_small_d3`, `lut_small_d4`
/// 
/// ...and...
/// 
/// `lut_large_d2`, `lut_large_d3`, `lut_large_d4`.
/// 
/// The `lut_large_*`
/// variants favour cpu performance over memory performance, the `lut_small_*`
/// variants strike a balance between CPU and memory performance, and omitting
/// cargo features will use no LUT at all. See Performance section below.
/// 
/// # Examples
/// ```rust
/// use insides::*;
///
/// let location = Hilbert::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
///
/// assert_eq!(location.index(), 36);
/// assert_eq!(location.coords(), [1, 2, 3]);
/// ```
/// 
/// # Performance
/// When using the `lut_large_d2` feature, the Hilbert implementation in
/// insides is on a par with the current fastest crate,
/// [Fast Hilbert](https://crates.io/crates/fast_hilbert).
/// 
/// A full breakdown of performance for each feature and similar libraries is
/// shown below. Tests were performed using all types supported by each crate
/// and times averaged. Times are specific to the machine running the benchmarks;
/// your results may vary.
/// 
/// ### 2D Benchmark
/// Benchmark performed on a 256x256 grid of 2D coordinates.
/// 
/// | Crate (and Features)                  | Index to Coords Time | Coords to Index Time |
/// | ------------------------------------- | -------------------- | -------------------- |
/// | insides (with feature `lut_large_d2`) | **373μs**            | 416μs                |
/// | fast_hilbert                          | 462μs                | **405μs**            |
/// | hilbert_2d                            | 533μs                | 697μs                |
/// | hilbert_curve                         | 616μs                | 643μs                |
/// | insides (with feature `lut_small_d2`) | 822μs                | 1.1ms                |
/// | insides                               | 1.5ms                | 2.2ms                |
/// | hilbert_index                         | 1.8ms                | 2.6ms                |
/// | hilbert                               | 18.1ms               | 14.3ms               |
///
/// ### 3D Benchmark
/// Benchmark performed on a 32x32x32 grid of 3D coordinates (performed twice to balance with other tables).
/// 
/// | Crate (and Features)                  | Index to Coords Time | Coords to Index Time |
/// | ------------------------------------- | -------------------- | -------------------- |
/// | insides (with feature `lut_large_d3`) | **600μs**            | **705μs**            |
/// | insides (with feature `lut_small_d3`) | 844μs                | 1.1ms                |
/// | insides                               | 1.2ms                | 2.0ms                |
/// | hilbert_index                         | 1.4ms                | 1.9ms                |
/// | hilbert                               | 17.1ms               | 14.3ms               |
///
/// ### 4D Benchmark
/// Benchmark performed on a 16x16x16x16 grid of 4D coordinates.
/// 
/// | Crate (and Features)                  | Index to Coords Time | Coords to Index Time |
/// | ------------------------------------- | -------------------- | -------------------- |
/// | insides (with feature `lut_large_d4`) | **394μs**            | **470μs**            |
/// | insides (with feature `lut_small_d4`) | 656μs                | 798μs                |
/// | hilbert_index                         | 1.1ms                | 1.1ms                |
/// | insides                               | 955μs                | 1.5ms                |
/// | hilbert                               | 17.4ms               | 14.5ms               |
/// 
// Until we have complex generic constants, we have to pass D in here (needed by coords)
// Waiting on: https://github.com/rust-lang/rust/issues/76560
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hilbert<DM, const D: usize>(pub(crate) DM::Dilated)
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex;
    
impl<DM, const D: usize> SpaceFillingCurve<D> for Hilbert<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex,
{
    type Coord = DM::Undilated;
    type Index = DM::Dilated;
    const COORD_BITS: usize = DM::UNDILATED_BITS;
    const COORD_MAX: Self::Coord = DM::UNDILATED_MAX;
    const INDEX_BITS: usize = DM::DILATED_BITS;
    const INDEX_MAX: Self::Index = DM::DILATED_MASK;
    const D: usize = D;

    #[inline(always)]
    fn from_index(index: Self::Index) -> Self {
        debug_assert!(
            index <= Self::INDEX_MAX,
            "Parameter 'index' exceeds maximum"
        );
        Self(index)
    }

    #[inline(always)]
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        let morton_index = Morton::<DM, D>::from_coords(coords).index();
        if morton_index.fits_in_usize() { // Is this really needed?
            Self(Self::Index::from_usize(morton_to_hilbert::<usize, D>(morton_index.to_usize())))
        } else {
            Self(morton_to_hilbert::<Self::Index, D>(morton_index))
        }
    }

    #[inline(always)]
    fn coords(&self) -> [Self::Coord; D] {
        if self.0.fits_in_usize() {
            Morton::<DM, D>::from_index(Self::Index::from_usize(hilbert_to_morton::<usize, D>(self.0.to_usize()))).coords()
        } else {
            Morton::<DM, D>::from_index(morton_to_hilbert::<Self::Index, D>(self.0)).coords()
        }
    }

    #[inline(always)]
    fn index(&self) -> Self::Index {
        self.0
    }
}

/// A Hilbert encoded space filling curve implementation using the Expand dilation method
// Not available until: https://github.com/rust-lang/rust/issues/112792
//pub type HilbertExpand<T, const D: usize> = Hilbert<Expand<T, D>, D> where T: CurveIndex;

/// A Hilbert encoded space filling curve implementation using the Expand dilation method
// Not available until: https://github.com/rust-lang/rust/issues/112792
//pub type HilbertFixed<T, const D: usize> = Hilbert<Fixed<T, D>, D> where T: CurveIndex;

#[cfg(test)]
mod tests {
    extern crate std;

    use std::{fs::File, io::Write, collections::HashMap};

    use super::*;

    use crate::internal::tests::test_curve;

    #[test]
    #[ignore]
    fn generate_lut() {
        generate_lut_small::<2>();
        generate_lut_small::<3>();
        generate_lut_small::<4>();
        generate_lut_large::<2>();
        generate_lut_large::<3>();
        generate_lut_large::<4>();
    }

    fn generate_lut_small<const D: usize>() {
        let num_elems: usize = 1 << D;

        println!("#[cfg(feature = \"lut_small_d{D}\")]");
        println!("mod d{D} {{");

        println!("    pub const SHR_INTO_GRAY_INV_LUT: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::gray_inverse::<_, D>(explicit::shr_cyclic::<_, D>(i, r))));
            println!("],");
        }
        println!("    ];\n");

        println!("    pub const GRAY_INTO_SHL_LUT: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::shl_cyclic::<_, D>(explicit::gray(i), r)));
            println!("],");
        }
        println!("    ];\n");

        println!("    pub const FLIP_TX_INTO_SHL_LUT: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::shl_cyclic::<_, D>(explicit::flip_transform(i), r)));
            println!("],");
        }
        println!("    ];\n");

        println!("    pub const ROTATION_LUT: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", (explicit::rotation_transform(i) + r + 1) % D));
            println!("],");
        }
        println!("    ];");

        println!("}}\n");
    }

    fn generate_lut_large<const D: usize>() {
        let num_children: usize = 1 << D;
        let num_axis_flips: usize = 1 << (D - 1);
        let num_rotations: usize = D;
        let max_rotation: usize = num_rotations - 1;
        let num_transforms: usize = num_rotations * num_axis_flips;
        let stride = (1usize..).into_iter().take_while(|i| num_children.pow(*i as u32) <= 256).last().unwrap();
        let num_combined_children = num_children.pow(stride as u32);
        let lower_mask: usize = num_children - 1;

        // Enumerate and map all possible axis flips
        let mut flip_ids = HashMap::<usize, usize>::default();
        let mut flip_axes = Vec::default();
        for flip in (0..num_children).map(|i| explicit::flip_transform(i)) {
            if !flip_ids.contains_key(&flip) {
                flip_ids.insert(flip, flip_axes.len());
                flip_axes.push(flip);
            }
        }
        assert_eq!(flip_ids.len(), num_axis_flips);
        assert_eq!(flip_axes.len(), num_axis_flips);

        // Build raw transforms
        let transforms: Vec<Vec<(usize, usize)>> = (0..num_transforms).map(|tx_id| {
            // Get current flip and rotation from transform
            let (flip, rotation) = (flip_axes[tx_id / num_rotations], max_rotation - (tx_id % num_rotations));

            (0..num_children).map(|index| {
                // Save hilbert to morton transformer for this transform and index
                let morton = flip ^ explicit::shl_cyclic::<_, D>(explicit::gray(index), rotation);

                // To find what the next transform would be, we flip and rotate again
                let flip = flip ^ explicit::shl_cyclic::<_, D>(explicit::flip_transform(index), rotation);
                let rotation = (explicit::rotation_transform(index) + rotation + 1) % D;
                
                // Next transform id is a magic combination of flip and rotation
                let next_tx_id = flip_ids.get(&flip).unwrap() * num_rotations + (max_rotation - rotation);
                (next_tx_id, morton)
            }).collect()
        }).collect();

        // Build combined transforms
        let mut combined_transforms: Vec<Vec<(usize, usize, usize)>> = (0..num_transforms).map(|tx_id| {
            (0..num_combined_children).map(|combined_partial| {
                let mut combined_tx: usize = tx_id;
                let mut combined_morton: usize = 0;

                // Iterates backwards through each element stride would step over
                for shift in (0..stride).into_iter().rev().map(|s| s * D) {
                    // Extract a portion of the combined partial (which by magic happens to be the hilbert partial)
                    let partial = (combined_partial >> shift) & lower_mask;

                    // Build a combined hilbert to morton transformer
                    combined_morton |= transforms[combined_tx][partial].1 << shift;

                    // Move to the next transform
                    combined_tx = transforms[combined_tx][partial].0;
                }
                (combined_tx, combined_morton, 0)
            }).collect()
        }).collect();

        // Build morton to hilbert transformer
        // It's slightly easier to do this as a separate step, since we know for sure what the values should be
        for tx_id in 0..num_transforms {
            for combined_partial in 0..num_combined_children {
                let morton = combined_transforms[tx_id][combined_partial].1;
                combined_transforms[tx_id][morton].2 = combined_partial;
            }
        }

        // Output Rust code to std out
        // Also output LUTs to dedicated files (since they're rather large)
        println!("#[cfg(feature = \"lut_large_d{D}\")]");
        println!("mod d{D} {{");
        println!("    pub const STRIDE: usize = {stride};");
        let transform_lut_file = format!("hilbert_d{D}_transform_lut.in");
        match write_lut_file(transform_lut_file.as_str(), &combined_transforms, |(t, _, _)| *t) {
            Err(e) => panic!("Failed to write transform lut: {e}"),
            _ => {}
        }
        println!("    pub const TRANSFORM_LUT: [[u8; {num_combined_children}]; {num_transforms}] = include!(\"lut/{transform_lut_file}\");");

        let morton_lut_file = format!("hilbert_d{D}_morton_lut.in");
        match write_lut_file(morton_lut_file.as_str(), &combined_transforms, |(_, m, _)| *m) {
            Err(e) => panic!("Failed to write morton lut: {e}"),
            _ => {}
        }
        println!("    pub const MORTON_LUT: [[u8; {num_combined_children}]; {num_transforms}] = include!(\"lut/{morton_lut_file}\");");

        let hilbert_lut_file = format!("hilbert_d{D}_hilbert_lut.in");
        match write_lut_file(hilbert_lut_file.as_str(), &combined_transforms, |(_, _, h)| *h) {
            Err(e) => panic!("Failed to write hilbert lut: {e}"),
            _ => {}
        }
        println!("    pub const HILBERT_LUT: [[u8; {num_combined_children}]; {num_transforms}] = include!(\"lut/{hilbert_lut_file}\");");
        println!("}}\n");
    }

    fn write_lut_file<T, F>(file: &str, transforms: &Vec<Vec<T>>, element: F) -> std::io::Result<()>
    where
        F: Fn(&T) -> usize
    {
        let file = format!("./src/lut/{file}");
        let mut output = File::create(file)?;
        writeln!(output, "[")?;
        for i in transforms {
            write!(output, "    [\n        ")?;
            let mut num_in_line = 0;
            for j in i {
                if num_in_line >= 32 {
                    write!(output, "\n        ")?;
                    num_in_line = 0;
                }
                write!(output, "{:0>3}, ", element(j))?;
                num_in_line += 1;
            }
            writeln!(output, "\n    ],")?;
        }
        writeln!(output, "]")?;
        Ok(())
    }

    macro_rules! hilbert_expand {
        ($t:ty, $d:literal) => {
            Hilbert::<Expand<$t, $d>, $d>
        };
    }

    macro_rules! hilbert_fixed {
        ($t:ty, $d:literal) => {
            Hilbert::<Fixed<$t, $d>, $d>
        };
    }

    macro_rules! hilbert_tests {
        ($curve:path, $t:ty, $($d:literal),+) => {
            test_curve!($curve, false, $t, $($d),+);
//            test_curve_siblings!($curve, false, $t, $($d),+);
//            test_curve_neighbours!($curve, false, $t, $($d),+);
        };
    }

    hilbert_tests!(hilbert_expand, u8, 2, 3, 4, 5, 6, 7, 8);
    hilbert_tests!(hilbert_expand, u16, 2, 3, 4, 5, 6, 7, 8);
    hilbert_tests!(hilbert_expand, u32, 2, 3, 4);
    hilbert_tests!(hilbert_expand, u64, 2);

    hilbert_tests!(hilbert_fixed, u8, 2, 3, 4);
    hilbert_tests!(hilbert_fixed, u16, 2, 3, 4, 5, 6, 7, 8);
    hilbert_tests!(hilbert_fixed, u32, 2, 3, 4, 5, 6, 7, 8);
    hilbert_tests!(hilbert_fixed, u64, 2, 3, 4, 5, 6, 7, 8);
    hilbert_tests!(hilbert_fixed, u128, 2, 3, 4, 5, 6, 7, 8);
}
