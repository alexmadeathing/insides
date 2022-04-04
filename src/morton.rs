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

use super::{CurveCoord, CurveIndex, SpaceFillingCurve, Neighbours, QueryDirection, Siblings};

use crate::internal::NumTraits;

/// A Morton encoded space filling curve implementation
/// 
/// Morton encoding, also known as a
/// [Z-order curve](https://en.wikipedia.org/wiki/Z-order_curve), is a space
/// filling algorithm which maps a multidimensional set of coordinates to one
/// dimension, achieved by interleaving the bit sequence of each coordinate
/// value.
/// 
/// Whilst other encoding methods may exhibit better spatial locality (such as
/// the Hilbert curve), the Morton curve offers excellent CPU performance,
/// since most behaviours can be reduced to a simple set of bitwise operations,
/// making it an ideal choice for applications such and quad trees and octrees.
/// 
/// # Examples
/// ```rust
/// use insides::*;
/// 
/// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
/// 
/// assert_eq!(location.index(), 0b110101);
/// assert_eq!(location.coords(), [1, 2, 3]);
/// ```
// Until we have complex generic constants, we have to pass D in here (needed by coords)
// Waiting on: https://github.com/rust-lang/rust/issues/76560
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Morton<DM, const D: usize>(pub(crate) DM::Dilated)
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex;

impl<DM, const D: usize> SpaceFillingCurve<D> for Morton<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex,
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

    #[inline]
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        debug_assert!(
            *coords.iter().max().unwrap() <= Self::COORD_MAX,
            "Parameter 'coords' contains a value which exceeds maximum"
        );
        let mut v = DM::Dilated::zero();
        for (axis, coord) in coords.into_iter().enumerate() {
            v = v.bit_or(DM::dilate(coord).value().shl(axis));
        }
        Self(v)
    }

    #[inline]
    fn coords(&self) -> [Self::Coord; D] {
        crate::internal::array_from_fn::<_, _, D>(|i| DilatedInt::<DM>::new(self.0.shr(i).bit_and(DM::DILATED_MAX)).undilate())
    }

    #[inline]
    fn index(&self) -> Self::Index {
        self.0
    }
}

impl<DM, const D: usize> Siblings<D> for Morton<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex,
{
    #[inline]
    fn sibling_on_axis(&self, axis: usize) -> Self {
        Self(self.0.bit_xor(DM::Dilated::one().shl(axis)))
    }

    #[inline]
    fn sibling_or_same_on_axis(&self, axis: usize, direction: QueryDirection) -> Self {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");
        let lower_axis_mask = DM::Dilated::one().shl(axis);
        let search_mask = match direction {
            QueryDirection::Positive => lower_axis_mask,
            QueryDirection::Negative => DM::Dilated::zero(),
        };
        Self(
            self.0
                .bit_and(lower_axis_mask.bit_not())
                .bit_or(search_mask),
        )
    }
}

impl<DM, const D: usize> Neighbours<D> for Morton<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex,
{
    #[inline]
    fn neighbour_on_axis(&self, axis: usize, direction: QueryDirection) -> Option<Self> {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");

        let coord = DilatedInt::<DM>::new(self.0.shr(axis).bit_and(DM::DILATED_MAX));
        let coord = match direction {
            QueryDirection::Positive => if coord.value() < DM::DILATED_MAX { Some(coord.add_one()) } else { None },
            QueryDirection::Negative => if coord.value() > NumTraits::zero() { Some(coord.sub_one()) } else { None },
        };
        coord.map(|coord| {
            let index = self.0.bit_and(DM::DILATED_MAX.shl(axis).bit_not());
            Self(index.bit_or(coord.value().shl(axis)))
        })
    }

    #[inline]
    fn neighbour_on_axis_wrapping(&self, axis: usize, direction: QueryDirection) -> Self {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");

        let coord = DilatedInt::<DM>::new(self.0.shr(axis).bit_and(DM::DILATED_MAX));
        let coord = match direction {
            QueryDirection::Positive => coord.add_one(),
            QueryDirection::Negative => coord.sub_one(),
        };
        let index = self.0.bit_and(DM::DILATED_MAX.shl(axis).bit_not());
        Self(index.bit_or(coord.value().shl(axis)))
    }
}

#[cfg(test)]
mod tests {
    use crate::internal::tests::{test_curve, test_curve_siblings, test_curve_neighbours};

    macro_rules! morton_expand {
        ($t:ty, $d:literal) => {
            Morton::<Expand<$t, $d>, $d>
        };
    }

    macro_rules! morton_fixed {
        ($t:ty, $d:literal) => {
            Morton::<Fixed<$t, $d>, $d>
        };
    }

    macro_rules! morton_tests {
        ($curve:path, $t:ty, $($d:literal),+) => {
            test_curve!($curve, false, $t, $($d),+);
            test_curve_siblings!($curve, false, $t, $($d),+);
            test_curve_neighbours!($curve, false, $t, $($d),+);
        };
    }

    morton_tests!(morton_expand, u8, 2, 3, 4, 5, 6, 7, 8);
    morton_tests!(morton_expand, u16, 2, 3, 4, 5, 6, 7, 8);
    morton_tests!(morton_expand, u32, 2, 3, 4);
    morton_tests!(morton_expand, u64, 2);

    morton_tests!(morton_fixed, u8, 2, 3, 4);
    morton_tests!(morton_fixed, u16, 2, 3, 4, 5, 6, 7, 8);
    morton_tests!(morton_fixed, u32, 2, 3, 4, 5, 6, 7, 8);
    morton_tests!(morton_fixed, u64, 2, 3, 4, 5, 6, 7, 8);
    morton_tests!(morton_fixed, u128, 2, 3, 4, 5, 6, 7, 8);
}
