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

use super::{Coord, SpaceFillingCurve, Index, Neighbours, QueryDirection, Siblings};

use crate::{internal::NumTraits, Morton};

#[cfg(feature = "lut16")]
const GRAY_TAB: [usize; 16] = [
    00, 01, 03, 02, 06, 07, 05, 04, 12, 13, 15, 14, 10, 11, 09, 08,
];

#[cfg(feature = "lut32")]
const GRAY_TAB: [usize; 32] = [
    00, 01, 03, 02, 06, 07, 05, 04, 12, 13, 15, 14, 10, 11, 09, 08,
    24, 25, 27, 26, 30, 31, 29, 28, 20, 21, 23, 22, 18, 19, 17, 16,
];

#[cfg(feature = "lut64")]
const GRAY_TAB: [usize; 64] = [
    00, 01, 03, 02, 06, 07, 05, 04, 12, 13, 15, 14, 10, 11, 09, 08,
    24, 25, 27, 26, 30, 31, 29, 28, 20, 21, 23, 22, 18, 19, 17, 16,
    48, 49, 51, 50, 54, 55, 53, 52, 60, 61, 63, 62, 58, 59, 57, 56,
    40, 41, 43, 42, 46, 47, 45, 44, 36, 37, 39, 38, 34, 35, 33, 32,
];

#[cfg(feature = "lut16")]
const GRAY_INV_TAB: [usize; 16] = [
    00, 01, 03, 02, 07, 06, 04, 05, 15, 14, 12, 13, 08, 09, 11, 10,
];

#[cfg(feature = "lut32")]
const GRAY_INV_TAB: [usize; 32] = [
    00, 01, 03, 02, 07, 06, 04, 05, 15, 14, 12, 13, 08, 09, 11, 10,
    31, 30, 28, 29, 24, 25, 27, 26, 16, 17, 19, 18, 23, 22, 20, 21,
];

#[cfg(feature = "lut64")]
const GRAY_INV_TAB: [usize; 64] = [
    00, 01, 03, 02, 07, 06, 04, 05, 15, 14, 12, 13, 08, 09, 11, 10,
    31, 30, 28, 29, 24, 25, 27, 26, 16, 17, 19, 18, 23, 22, 20, 21,
    63, 62, 60, 61, 56, 57, 59, 58, 48, 49, 51, 50, 55, 54, 52, 53,
    32, 33, 35, 34, 39, 38, 36, 37, 47, 46, 44, 45, 40, 41, 43, 42,
];

#[cfg(feature = "lut16")]
const ROTATION_TRANSFORM_TAB: [usize; 16] = [
    0, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 4,
];

#[cfg(feature = "lut32")]
const ROTATION_TRANSFORM_TAB: [usize; 32] = [
    0, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 4,
    4, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 5,
];

#[cfg(feature = "lut64")]
const ROTATION_TRANSFORM_TAB: [usize; 64] = [
    0, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 4,
    4, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 5,
    5, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 4,
    4, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1, 6,
];

#[cfg(feature = "lut16")]
const AXIS_TRANSFORM_TAB: [usize; 16] = [
    00, 00, 00, 03, 03, 06, 06, 05, 05, 12, 12, 15, 15, 10, 10, 09,
];

#[cfg(feature = "lut32")]
const AXIS_TRANSFORM_TAB: [usize; 32] = [
    00, 00, 00, 03, 03, 06, 06, 05, 05, 12, 12, 15, 15, 10, 10, 09,
    09, 24, 24, 27, 27, 30, 30, 29, 29, 20, 20, 23, 23, 18, 18, 17,
];

#[cfg(feature = "lut64")]
const AXIS_TRANSFORM_TAB: [usize; 64] = [
    00, 00, 00, 03, 03, 06, 06, 05, 05, 12, 12, 15, 15, 10, 10, 09,
    09, 24, 24, 27, 27, 30, 30, 29, 29, 20, 20, 23, 23, 18, 18, 17,
    17, 48, 48, 51, 51, 54, 54, 53, 53, 60, 60, 63, 63, 58, 58, 57,
    57, 40, 40, 43, 43, 46, 46, 45, 45, 36, 36, 39, 39, 34, 34, 33,
];

#[inline]
fn gray<T>(i: T) -> T where T: NumTraits {
    #[cfg(any(feature = "lut16", feature = "lut32", feature = "lut64"))]
    if i < NumTraits::from_usize(GRAY_TAB.len()) {
        return NumTraits::from_usize(GRAY_TAB[i.to_usize()]);
    }
    i.bit_xor(i.shr(NumTraits::one()))
}

#[inline]
fn gray_inverse<T, const D: usize>(i: T) -> T where T: NumTraits {
    #[cfg(any(feature = "lut16", feature = "lut32", feature = "lut64"))]
    if i < NumTraits::from_usize(GRAY_INV_TAB.len()) {
        return NumTraits::from_usize(GRAY_INV_TAB[i.to_usize()]);
    }
    (1..D).fold(i, |acc, shift| acc.bit_xor(i.shr(shift)))
}

#[inline]
fn shl_cyclic<T, const D: usize>(i: T, a: usize, lower_mask: T) -> T where T: NumTraits {
    NumTraits::bit_and(lower_mask, i.shl(a).bit_or(i.shr(D - a)))
}

#[inline]
fn shr_cyclic<T, const D: usize>(i: T, a: usize, lower_mask: T) -> T where T: NumTraits {
    NumTraits::bit_and(lower_mask, i.shr(a).bit_or(i.shl(D - a)))
}

#[inline]
fn rotation_transform(i: usize) -> usize {
    #[cfg(any(feature = "lut16", feature = "lut32", feature = "lut64"))]
    if i < ROTATION_TRANSFORM_TAB.len() {
        return ROTATION_TRANSFORM_TAB[i];
    }
    if i > 0 {
        i.wrapping_add((i & 0x1).wrapping_sub(1)).trailing_ones() as usize
    } else {
        0
    }
}

#[inline]
fn axis_transform<T>(i: T) -> T where T: NumTraits {
    #[cfg(any(feature = "lut16", feature = "lut32", feature = "lut64"))]
    if i < NumTraits::from_usize(AXIS_TRANSFORM_TAB.len()) {
        return NumTraits::from_usize(AXIS_TRANSFORM_TAB[i.to_usize()]);
    }
    if i == NumTraits::zero() {
        NumTraits::zero()
    } else {
        gray(NumTraits::bit_and(i.sub(NumTraits::one()), NumTraits::bit_not(NumTraits::one())))
    }
}

// Until we have complex generic constants, we have to pass D in here (needed by coords)
// Waiting on: https://github.com/rust-lang/rust/issues/76560
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hilbert<DM, const D: usize>(pub(crate) DM::Dilated)
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: Coord,
    DM::Dilated: Index;

impl<DM, const D: usize> SpaceFillingCurve<D> for Hilbert<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: Coord,
    DM::Dilated: Index,
{
    type Coord = DM::Undilated;
    type Index = DM::Dilated;
    const COORD_MAX: Self::Coord = DM::UNDILATED_MAX;
    const INDEX_MAX: Self::Index = DM::DILATED_MASK;
    const D: usize = D;

    #[inline]
    fn from_index(index: Self::Index) -> Self {
        debug_assert!(
            index <= Self::INDEX_MAX,
            "Parameter 'index' exceeds maximum"
        );
        Self(index)
    }

    fn from_coords(coords: [Self::Coord; D]) -> Self {
        let morton_index = Morton::<DM, D>::from_coords(coords).index();

        let lower_mask = Self::Index::one().shl(D).sub(Self::Index::one());

        let order = (Self::Index::bits() - morton_index.lz() + D - 1) / D;
        let skip_orders = DM::UNDILATED_BITS - order;

        let mut hilbert_index = Self::Index::zero();
        let mut flip_axes = Self::Index::zero();
        let mut rotate_amount = (skip_orders + 1) % D;
        for i in (0..order).rev() {
            let morton_partial = morton_index.shr(i * D).bit_and(lower_mask);
            let raw_index = gray_inverse::<_, D>(shr_cyclic::<_, D>(flip_axes.bit_xor(morton_partial), rotate_amount, lower_mask));
            hilbert_index = hilbert_index.bit_or(raw_index.shl(i * D));
            flip_axes = flip_axes.bit_xor(shl_cyclic::<_, D>(axis_transform(raw_index), rotate_amount, lower_mask));
            rotate_amount = (rotation_transform(raw_index.to_usize()) + rotate_amount + 1) % D;
        }
        Self(hilbert_index)
    }

    fn coords(&self) -> [Self::Coord; D] {
        let lower_mask = Self::Index::one().shl(D).sub(Self::Index::one());

        let order = (Self::Index::bits() - self.0.lz() + D - 1) / D;
        let skip_orders = DM::UNDILATED_BITS - order;

        let mut morton_index = Self::Index::zero();
        let mut flip_axes = Self::Index::zero();
        let mut rotate_amount = (skip_orders + 1) % D;
        for i in (0..order).rev() {
            let raw_index = self.0.shr(i * D).bit_and(lower_mask);
            let morton_partial = flip_axes.bit_xor(shl_cyclic::<_, D>(gray(raw_index), rotate_amount, lower_mask));
            morton_index = morton_index.bit_or(morton_partial.shl(i * D));
            flip_axes = flip_axes.bit_xor(shl_cyclic::<_, D>(axis_transform(raw_index), rotate_amount, lower_mask));
            rotate_amount = (rotation_transform(raw_index.to_usize()) + rotate_amount + 1) % D;
        }
        Morton::<DM, D>::from_index(morton_index).coords()
    }

    #[inline]
    fn index(&self) -> Self::Index {
        self.0
    }
}

#[cfg(test)]
mod tests {
    extern crate std;

    use crate::internal::tests::test_curve;

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
        for i in 0..10000usize {
            let g = super::gray(i);
            assert_eq!(super::gray_inverse::<_, { usize::BITS as usize }>(g), i);
        }
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
}
