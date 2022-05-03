// ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)
//
// Copyright © 2022 Alex Blunt (alexmadeathing)
//
// This is anti-capitalist software, released for free use by individuals and
// organizations that do not operate by capitalist principles.
//
// Permission is hereby granted, free of charge, to any person or organization
// (the "User") obtaining a copy of this software and associated documentation
// files (the "Software"), to use, copy, modify, merge, distribute, and/or sell
// copies of the Software, subject to the following conditions:
//
// 1. The above copyright notice and this permission notice shall be included in
// all copies or modified versions of the Software.
//
// 2. The User is one of the following:
//   a. An individual person, laboring for themselves
//   b. A non-profit organization
//   c. An educational institution
//   d. An organization that seeks shared profit for all of its members, and
//      allows non-members to set the cost of their labor
//
// 3. If the User is an organization with owners, then all owners are workers
// and all workers are owners with equal equity and/or equal vote.
//
// 4. If the User is an organization, then the User is not law enforcement or
// military, or working for or under either.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY
// KIND, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use dilate::*;

use super::{CurveCoord, CurveIndex, Neighbours, QueryDirection, Siblings, SpaceFillingCurve};

use crate::{internal::NumTraits, Morton};

mod lut_small {
    use super::NumTraits;

    #[cfg(feature = "lut_small_d2")]
    mod d2 {
        pub const AXIS_TX_INTO_SHL_LUT: [[u8; 4]; 2] = [
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

        pub const ROTATION_INCREMENT_LUT: [[u8; 4]; 2] = [
            [01, 00, 00, 01, ],
            [00, 01, 01, 00, ],
        ];
    }

    #[cfg(feature = "lut_small_d3")]
    mod d3 {
        pub const AXIS_TX_INTO_SHL_LUT: [[u8; 8]; 3] = [
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

        pub const ROTATION_INCREMENT_LUT: [[u8; 8]; 3] = [
            [01, 02, 02, 00, 00, 02, 02, 01, ],
            [02, 00, 00, 01, 01, 00, 00, 02, ],
            [00, 01, 01, 02, 02, 01, 01, 00, ],
        ];
    }

    #[cfg(feature = "lut_small_d4")]
    mod d4 {
        pub const AXIS_TX_INTO_SHL_LUT: [[u8; 16]; 4] = [
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

        pub const ROTATION_INCREMENT_LUT: [[u8; 16]; 4] = [
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
    fn axis_tx_into_shl_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> T {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => T::from_u8(d2::AXIS_TX_INTO_SHL_LUT[rotate_amount][index.to_usize()]),
            #[cfg(feature = "lut_small_d3")]
            3 => T::from_u8(d3::AXIS_TX_INTO_SHL_LUT[rotate_amount][index.to_usize()]),
            #[cfg(feature = "lut_small_d4")]
            4 => T::from_u8(d4::AXIS_TX_INTO_SHL_LUT[rotate_amount][index.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn shr_into_gray_inv_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> T {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => T::from_u8(d2::SHR_INTO_GRAY_INV_LUT[rotate_amount][index.to_usize()]),
            #[cfg(feature = "lut_small_d3")]
            3 => T::from_u8(d3::SHR_INTO_GRAY_INV_LUT[rotate_amount][index.to_usize()]),
            #[cfg(feature = "lut_small_d4")]
            4 => T::from_u8(d4::SHR_INTO_GRAY_INV_LUT[rotate_amount][index.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn gray_into_shl_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> T {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => T::from_u8(d2::GRAY_INTO_SHL_LUT[rotate_amount][index.to_usize()]),
            #[cfg(feature = "lut_small_d3")]
            3 => T::from_u8(d3::GRAY_INTO_SHL_LUT[rotate_amount][index.to_usize()]),
            #[cfg(feature = "lut_small_d4")]
            4 => T::from_u8(d4::GRAY_INTO_SHL_LUT[rotate_amount][index.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn rotation_increment_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> usize {
        match D {
            #[cfg(feature = "lut_small_d2")]
            2 => d2::ROTATION_INCREMENT_LUT[rotate_amount][index.to_usize()] as usize,
            #[cfg(feature = "lut_small_d3")]
            3 => d3::ROTATION_INCREMENT_LUT[rotate_amount][index.to_usize()] as usize,
            #[cfg(feature = "lut_small_d4")]
            4 => d4::ROTATION_INCREMENT_LUT[rotate_amount][index.to_usize()] as usize,
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn walk_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        min_order: usize,
        mut f: F,
    ) -> (usize, T)
    where
        F: FnMut(T, usize, usize, T) -> T,
    {
        let lower_mask = T::one().shl(D).sub(T::one());

        // A subset of highest bits of index will be 0, therefore, we can skip a certain
        // number of iterations in the following loop. Note that skipping iterations
        // still requires us to rotate the pattern as if we had iterated all orders, see
        // rotate_amount below.
        let order = (T::bits() - index.lz() + D - 1) / D;

        // We iterate in reverse because higher orders affect the rotation of lower orders
        let mut flip_axes = T::zero();
        let mut rotate_amount = (D - 1) - (order % D);
        for shift in (min_order..order).rev().map(|i| i * D) {
            // Extract a portion of the index
            let index_partial = index.shr(shift).bit_and(lower_mask);

            // Context specific code generates the hilbert partial
            let hilbert_partial = f(index_partial, shift, rotate_amount, flip_axes);

            // Calculate the next axis flip flags and rotation amounts resulting from this order
            flip_axes = flip_axes.bit_xor(axis_tx_into_shl_lut::<_, D>(hilbert_partial, rotate_amount));
            rotate_amount = rotation_increment_lut::<_, D>(hilbert_partial, rotate_amount);
        }
        (rotate_amount, flip_axes)
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        debug_assert!(has_lut::<D>());
        let mut hilbert_index = T::zero();
        walk_transforms::<_, _, D>(
            morton_index,
            0,
            |morton_partial, shift, rotate_amount, flip_axes| {
                let hilbert_partial = shr_into_gray_inv_lut::<_, D>(flip_axes.bit_xor(morton_partial), rotate_amount);
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        debug_assert!(has_lut::<D>());
        let mut morton_index = T::zero();
        walk_transforms::<_, _, D>(
            hilbert_index,
            0,
            |hilbert_partial, shift, rotate_amount, flip_axes| {
                let morton_partial = flip_axes.bit_xor(gray_into_shl_lut::<_, D>(hilbert_partial, rotate_amount));
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
    fn combined_transform_lut<T: NumTraits, const D: usize>(
        transform: usize,
        hilbert_partial: T,
    ) -> usize {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => d2::TRANSFORM_LUT[transform][hilbert_partial.to_usize()] as usize,
            #[cfg(feature = "lut_large_d3")]
            3 => d3::TRANSFORM_LUT[transform][hilbert_partial.to_usize()] as usize,
            #[cfg(feature = "lut_large_d4")]
            4 => d4::TRANSFORM_LUT[transform][hilbert_partial.to_usize()] as usize,
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn combined_morton_lut<T: NumTraits, const D: usize>(
        transform: usize,
        hilbert_partial: T,
    ) -> T {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => T::from_u8(d2::MORTON_LUT[transform][hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d3")]
            3 => T::from_u8(d3::MORTON_LUT[transform][hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d4")]
            4 => T::from_u8(d4::MORTON_LUT[transform][hilbert_partial.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn combined_hilbert_lut<T: NumTraits, const D: usize>(
        transform: usize,
        hilbert_partial: T,
    ) -> T {
        match D {
            #[cfg(feature = "lut_large_d2")]
            2 => T::from_u8(d2::HILBERT_LUT[transform][hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d3")]
            3 => T::from_u8(d3::HILBERT_LUT[transform][hilbert_partial.to_usize()]),
            #[cfg(feature = "lut_large_d4")]
            4 => T::from_u8(d4::HILBERT_LUT[transform][hilbert_partial.to_usize()]),
            _ => unimplemented!(),
        }
    }

    #[inline(always)]
    fn walk_combined_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        mut f: F,
    ) -> usize
    where
        F: FnMut(usize, T, usize) -> T,
    {
        let stride_bits = stride::<D>() * D;
        let lower_mask = (0..stride_bits).fold(T::zero(), |lm, _| lm.shl(1).bit_or(T::one()));

        let order = ((T::bits() + stride_bits - 1) - index.lz()) / stride_bits;

        let mut transform = order * stride::<D>() % D;
        for shift in (0..order).rev().map(|i| i * stride_bits) {
            let hilbert_partial = f(transform, index.shr(shift).bit_and(lower_mask), shift);
            transform = combined_transform_lut::<_, D>(transform, hilbert_partial);
        }
        transform
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        debug_assert!(has_lut::<D>());
        let mut hilbert_index = T::zero();
        walk_combined_transforms::<_, _, D>(
            morton_index,
            |transform, morton_partial, shift| {
                let hilbert_partial = combined_hilbert_lut::<_, D>(transform, morton_partial);
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        debug_assert!(has_lut::<D>());
        let mut morton_index = T::zero();
        walk_combined_transforms::<_, _, D>(
            hilbert_index,
            |transform, hilbert_partial, shift| {
                let morton_partial = combined_morton_lut::<_, D>(transform, hilbert_partial);
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
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
    pub fn axis_transform<T: NumTraits>(index: T) -> T {
        if index == T::zero() {
            T::zero()
        } else {
            gray(index.sub(T::one()).bit_and(T::bit_not(T::one())))
        }
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
        // A subset of highest bits of index will be 0, therefore, we can skip a certain
        // number of iterations in the following loop. Note that skipping iterations
        // still requires us to rotate the pattern as if we had iterated all orders, see
        // rotate_amount below.
        let order = (T::bits() - index.lz() + D - 1) / D;

        // We iterate in reverse because higher orders affect the rotation of lower orders
        let mut flip_axes = T::zero();
        let mut rotate_amount = (D - 1) - (order % D);
        for shift in (min_order..order).rev().map(|i| i * D) {
            // Extract a portion of the index
            let index_partial = index.shr(shift).bit_and(lower_mask);

            // Context specific code generates the hilbert partial
            let hilbert_partial = f(index_partial, shift, rotate_amount, flip_axes);

            // Calculate the next axis flip flags and rotation amounts resulting from this order
            flip_axes = flip_axes.bit_xor(shl_cyclic::<_, D>(
                axis_transform(hilbert_partial),
                rotate_amount,
            ));
            rotate_amount =
                (rotation_transform(hilbert_partial.to_usize()) + rotate_amount + 1) % D;
        }
        (rotate_amount, flip_axes)
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        let mut hilbert_index = T::zero();
        walk_transforms::<_, _, D>(
            morton_index,
            0,
            |morton_partial, shift, rotate_amount, flip_axes| {
                let hilbert_partial = gray_inverse::<_, D>(shr_cyclic::<_, D>(
                    flip_axes.bit_xor(morton_partial),
                    rotate_amount,
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
            |hilbert_partial, shift, rotate_amount, flip_axes| {
                let morton_partial = flip_axes.bit_xor(shl_cyclic::<_, D>(
                    gray(hilbert_partial),
                    rotate_amount,
                ));
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
    }
}

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

    #[inline(never)] // Temporary
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        let morton_index = Morton::<DM, D>::from_coords(coords).index();

        if lut_large::has_lut::<D>() {
            return Self(lut_large::morton_to_hilbert::<_, D>(morton_index));
        }

        if lut_small::has_lut::<D>() {
            return Self(lut_small::morton_to_hilbert::<_, D>(morton_index));
        }

        Self(explicit::morton_to_hilbert::<_, D>(morton_index))
    }

    #[inline(never)] // Temporary
    fn coords(&self) -> [Self::Coord; D] {
        if lut_large::has_lut::<D>() {
            return Morton::<DM, D>::from_index(lut_large::hilbert_to_morton::<_, D>(self.0)).coords();
        }

        if lut_small::has_lut::<D>() {
            return  Morton::<DM, D>::from_index(lut_small::hilbert_to_morton::<_, D>(self.0)).coords();
        }

        Morton::<DM, D>::from_index(explicit::hilbert_to_morton::<_, D>(self.0)).coords()
    }

    #[inline(always)]
    fn index(&self) -> Self::Index {
        self.0
    }
}

impl<DM, const D: usize> Siblings<D> for Hilbert<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex,
{
    fn sibling_on_axis_toggle(&self, axis: usize) -> Self {
        todo!();
    }

    fn sibling_on_axis(&self, axis: usize, direction: QueryDirection) -> Self {
        todo!()
    }

    fn sibling_from_bits(&self, axis_bits: Self::Index) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    extern crate std;

    use std::{fs::File, io::Write};

    use dilate::DilationMethod;
    use super::*;

    use crate::internal::tests::{test_curve, test_curve_siblings};

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

    #[test]
    fn test_gray() {
//        for i in 0..10000usize {
//            let g = super::gray(i);
//            assert_eq!(super::gray_inverse::<_, { usize::BITS as usize }>(g), i);
//        }
    }

    #[test]
    fn moo() {
        let expected_index = 32;
        let coords = Hilbert::<dilate::Fixed<u8, 4>, 4>::from_index(expected_index).coords();
        let wrong = Hilbert::<dilate::Fixed<u8, 4>, 4>::from_coords(coords);
        assert_eq!(wrong.index(), expected_index);
    }

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

        println!("    pub const AXIS_TX_INTO_SHL_LUT: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::shl_cyclic::<_, D>(explicit::axis_transform(i), r)));
            println!("],");
        }
        println!("    ];\n");

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

        println!("    pub const ROTATION_INCREMENT_LUT: [[u8; {num_elems}]; {D}] = [");
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

        // Enumerate all possible axis flips
        // This forces the ordering of flip IDs to be known values
        let mut flip_axes_map: Vec<Option<usize>> = vec![None; num_children];
        let mut axis_flip_count = 0;
        for i in 0..num_children {
            let flip_axes = explicit::axis_transform(i);
            let flip_id = flip_axes_map[flip_axes].unwrap_or_else(|| {
                let a = axis_flip_count;
                axis_flip_count += 1;
                a
            });
            flip_axes_map[flip_axes] = Some(flip_id);
        }
        assert_eq!(axis_flip_count, num_axis_flips);

        // Enumerate all possible transforms
        // This forces the order of transform indices to be known values
        let mut transform_map: Vec<Vec<Option<usize>>> = vec![vec![None; num_rotations]; num_children];
        let mut transform_count = 0;
        for i in 0..num_children {
            for rotate_amount in 0..D {
                // Find initial flip and rotate for index partial i (note xor of previous flip_axes is not required here, but previous rotation is)
                let flip_axes = explicit::shl_cyclic::<_, D>(
                    explicit::axis_transform(i),
                    rotate_amount,
                );
                let rotate_amount = (explicit::rotation_transform(i) + rotate_amount + 1) % D;

                let flip_id = flip_axes_map[flip_axes].unwrap();

                // Flip rotation here avoids flip at runtime
                let flipped_rotation = max_rotation - rotate_amount;
                let tx_id = flip_id * num_rotations + flipped_rotation;

                // Store ID of this transform
                if transform_map[flip_axes][rotate_amount].is_none() {
                    transform_count += 1;
                }
                transform_map[flip_axes][rotate_amount] = Some(tx_id);
            }
        }
        assert_eq!(transform_count, num_transforms);

        // Build transforms
        let mut transforms: Vec<Vec<(usize, usize, usize)>> = vec![vec![(0, 0, 0); num_children]; num_transforms];
        for rotate_amount in 0..D {
            for i in 0..num_children {
                // Find initial flip and rotate for index partial i (note xor of previous flip_axes is not required here, but previous rotation is)
                let flip_axes = explicit::shl_cyclic::<_, D>(
                    explicit::axis_transform(i),
                    rotate_amount,
                );
                let rotate_amount = (explicit::rotation_transform(i) + rotate_amount + 1) % D;

                // ID of this transform
                let tx_id = transform_map[flip_axes][rotate_amount].unwrap();

                // Within this transform, it's possible to decend to a number of other transforms depending on the index partial j
                for j in 0..num_children {
                    // Store resultant morton partial and hilbert partial (to coords and from coords)
                    transforms[tx_id][j].1 = flip_axes ^ explicit::shl_cyclic::<_, D>(explicit::gray(j), rotate_amount);
                    transforms[tx_id][j].2 = explicit::gray_inverse::<_, D>(explicit::shr_cyclic::<_, D>(flip_axes ^ j, rotate_amount));

                    // To find what the next transform would be, we flip and rotate again
                    let flip_axes = flip_axes ^ explicit::shl_cyclic::<_, D>(explicit::axis_transform(j), rotate_amount);
                    let rotate_amount = (explicit::rotation_transform(j) + rotate_amount + 1) % D;

                    // Cache next transform ID
                    transforms[tx_id][j].0 = transform_map[flip_axes][rotate_amount].unwrap();
                }
            }
        }

        // Build combined transforms
        let mut combined_transforms: Vec<Vec<usize>> = vec![vec![0; num_combined_children]; num_transforms];
        let mut combined_mortons: Vec<Vec<usize>> = vec![vec![0; num_combined_children]; num_transforms];
        for t in 0..num_transforms {
            for combined_partial in 0..num_combined_children {
                let combined_tx = &mut combined_transforms[t][combined_partial];
                let combined_morton = &mut combined_mortons[t][combined_partial];
//                let (combined_tx, combined_morton, combined_hilbert) = &mut combined_transforms[t][combined_partial];
                *combined_tx = t;

                for s in (0..stride).into_iter().rev().map(|s| s * D) {
                    let partial = (combined_partial >> s) & lower_mask;
                    *combined_morton = *combined_morton | (transforms[*combined_tx][partial].1 << s);
//                    *combined_hilbert = *combined_hilbert | (transforms[*combined_tx][partial].2 << s);
                    *combined_tx = transforms[*combined_tx][partial].0;
                }
            }
        }

        // Build inverse combined transform
        let mut combined_hilberts: Vec<Vec<usize>> = vec![vec![0; num_combined_children]; num_transforms];
        for t in 0..num_transforms {
            for combined_partial in 0..num_combined_children {
                combined_hilberts[t][combined_mortons[t][combined_partial]] = combined_partial;
            }
        }

        println!("#[cfg(feature = \"lut_large_d{D}\")]");
        println!("mod d{D} {{");

//        println!("    pub const TRANSFORM_LUT: [[u8; {num_children}]; {num_transforms}] = [");
//        for i in 0..num_transforms {
//            print!("        [");
//            for j in 0..(num_children - 1) {
//                let lut_data = transforms[i][j];
//                print!("{:0>2}, ", lut_data.0);
//            }
//            let lut_data = transforms[i][num_children - 1];
//            println!("{:0>2}],", lut_data.0);
//        }
//        println!("    ];\n");
//        println!("    pub const MORTON_LUT: [[u8; {num_children}]; {num_transforms}] = [");
//        for i in 0..num_transforms {
//            print!("        [");
//            for j in 0..(num_children - 1) {
//                let lut_data = transforms[i][j];
//                print!("{:0>2}, ", lut_data.1);
//            }
//            let lut_data = transforms[i][num_children - 1];
//            println!("{:0>2}],", lut_data.1);
//        }
//        println!("    ];\n");
//        println!("    pub const HILBERT_LUT: [[u8; {num_children}]; {num_transforms}] = [");
//        for i in 0..num_transforms {
//            print!("        [");
//            for j in 0..(num_children - 1) {
//                let lut_data = transforms[i][j];
//                print!("{:0>2}, ", lut_data.2);
//            }
//            let lut_data = transforms[i][num_children - 1];
//            println!("{:0>2}],", lut_data.2);
//        }
//        println!("    ];");

        println!("    pub const STRIDE: usize = {stride};");

        let transform_lut_file = format!("hilbert_d{D}_transform_lut.in");
        match write_lut_file(transform_lut_file.as_str(), &combined_transforms, |t| *t) {
            Err(e) => panic!("Failed to write transform lut: {e}"),
            _ => {}
        }
        println!("    pub const TRANSFORM_LUT: [[u8; {num_combined_children}]; {num_transforms}] = include!(\"lut/{transform_lut_file}\");");

        let morton_lut_file = format!("hilbert_d{D}_morton_lut.in");
        match write_lut_file(morton_lut_file.as_str(), &combined_mortons, |m| *m) {
            Err(e) => panic!("Failed to write morton lut: {e}"),
            _ => {}
        }
        println!("    pub const MORTON_LUT: [[u8; {num_combined_children}]; {num_transforms}] = include!(\"lut/{morton_lut_file}\");");

        let hilbert_lut_file = format!("hilbert_d{D}_hilbert_lut.in");
        match write_lut_file(hilbert_lut_file.as_str(), &combined_hilberts, |h| *h) {
            Err(e) => panic!("Failed to write hilbert lut: {e}"),
            _ => {}
        }
        println!("    pub const HILBERT_LUT: [[u8; {num_combined_children}]; {num_transforms}] = include!(\"lut/{hilbert_lut_file}\");");

//        println!("    pub const COMBINED_TRANSFORM_LUT: [[u8; {num_combined_children}]; {num_transforms}] = [");
//        for i in 0..num_transforms {
//            print!("        [");
//            for j in 0..(num_combined_children - 1) {
//                let lut_data = combined_transforms[i][j];
//                print!("{:0>2}, ", lut_data.0);
//            }
//            let lut_data = combined_transforms[i][num_combined_children - 1];
//            println!("{:0>2}],", lut_data.0);
//        }
//        println!("    ];\n");
//        println!("    pub const COMBINED_MORTON_LUT: [[u8; {num_combined_children}]; {num_transforms}] = [");
//        for i in 0..num_transforms {
//            print!("        [");
//            for j in 0..(num_combined_children - 1) {
//                let lut_data = combined_transforms[i][j];
//                print!("{:0>2}, ", lut_data.1);
//            }
//            let lut_data = combined_transforms[i][num_combined_children - 1];
//            println!("{:0>2}],", lut_data.1);
//        }
//        println!("    ];\n");
//        println!("    pub const COMBINED_HILBERT_LUT: [[u8; {num_combined_children}]; {num_transforms}] = [");
//        for i in 0..num_transforms {
//            print!("        [");
//            for j in 0..(num_combined_children - 1) {
//                let lut_data = combined_transforms[i][j];
//                print!("{:0>2}, ", lut_data.2);
//            }
//            let lut_data = combined_transforms[i][num_combined_children - 1];
//            println!("{:0>2}],", lut_data.2);
//        }
//        println!("    ];");

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

    test_curve!(hilbert_expand, true, u8, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_expand, true, u16, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_expand, true, u32, 2, 3, 4);
    test_curve!(hilbert_expand, true, u64, 2);

    test_curve!(hilbert_fixed, true, u8, 2, 3, 4);
    test_curve!(hilbert_fixed, true, u16, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_fixed, true, u32, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_fixed, true, u64, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_fixed, true, u128, 2, 3, 4, 5, 6, 7, 8);

    //    test_curve_siblings!(hilbert_expand, true, u8, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_expand, true, u16, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_expand, true, u32, 2, 3, 4);
    //    test_curve_siblings!(hilbert_expand, true, u64, 2);
    //
    //    test_curve_siblings!(hilbert_fixed, true, u8, 2, 3, 4);
    //    test_curve_siblings!(hilbert_fixed, true, u16, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_fixed, true, u32, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_fixed, true, u64, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_fixed, true, u128, 2, 3, 4, 5, 6, 7, 8);
}
