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

use super::{Coord, Index, SpaceFillingCurve, Neighbours, QueryDirection, Siblings};

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
    DM::Undilated: Coord,
    DM::Dilated: Index;

impl<DM, const D: usize> SpaceFillingCurve<D> for Morton<DM, D>
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
    DM::Undilated: Coord,
    DM::Dilated: Index,
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
    DM::Undilated: Coord,
    DM::Dilated: Index,
{
    #[inline]
    fn neighbour_on_axis(&self, axis: usize, direction: QueryDirection) -> Self {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");

        // This is a faster eqivalent of converting to coords, adding or subtracting one, then converting back
        // It's faster because it bypasses the undilation and dilation stage and instead uses dilated arithmetic
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
    extern crate std;

    macro_rules! expand {
        ($t:ty, $d:literal) => {
            Expand::<$t, $d>
        };
    }

    macro_rules! fixed {
        ($t:ty, $d:literal) => {
            Fixed::<$t, $d>
        };
    }

    macro_rules! test_morton_curve {
        ($name:path, $dm_macro:ty, $t:ty, $($d:literal),+) => {$(
            paste::paste!{
                mod [< $name _ $dm_macro _ $t _d $d >] {
                    use super::*;
                    use super::super::*;

                    use std::vec::Vec;

                    type TestedCurve = Morton::<$dm_macro!($t, $d), $d>;

                    type Index = <$dm_macro!($t, $d) as DilationMethod>::Dilated;

                    const COORD_BITS: usize = <$dm_macro!($t, $d)>::UNDILATED_BITS;
                    const COORD_MAX: $t = <$dm_macro!($t, $d)>::UNDILATED_MAX;
                    const INDEX_MAX: Index = <$dm_macro!($t, $d)>::DILATED_MASK;

                    fn generate_coords_test_cases(num_test_cases: usize) -> Vec<([$t; $d], Index)> {
                        let mut test_cases = Vec::new();
                        let mut prev_coords = [0; $d];
                        for tc in 0..num_test_cases {
                            let mut coords = [0; $d];
                            for i in 0..$d {
                                // The exact coord doesn't matter here
                                // We just need to generate an assorted set of numbers
                                // I'm xor'ing with the previous test case just to spice it up a bit
                                coords[i] = ((((tc * $d) + i) ^ tc) as $t & COORD_MAX) ^ prev_coords[i];
                                prev_coords[i] = coords[i];
                            }
                            test_cases.push((coords, coords_to_index(coords)));
                        }

                        // Also test 0 and max
                        test_cases.push(([0; $d], 0));
                        test_cases.push(([COORD_MAX; $d], INDEX_MAX));

                        test_cases
                    }

                    // Convert coords to morton index
                    // Since we know how morton codes should work, we can reasonably generate them
                    fn coords_to_index(coords: [$t; $d]) -> Index {
                        let mut index = 0;
                        // This is a low performance, but easy to read way to make a morton index
                        for bit in 0..COORD_BITS {
                            for i in 0..$d {
                                let coord_bit = ((coords[i] >> bit) & 0x1) as Index;
                                let shift = bit * $d + i;
                                index |= coord_bit << shift;
                            }
                        }
                        index
                    }

                    #[test]
                    #[should_panic(expected = "Parameter 'index' exceeds maximum")]
                    #[allow(arithmetic_overflow)]
                    fn from_index_too_large_panics() {
                        if INDEX_MAX != Index::MAX {
                            TestedCurve::from_index(INDEX_MAX + 1);
                        } else {
                            panic!("Parameter 'index' exceeds maximum");
                        }
                    }

                    #[test]
                    fn from_index_stores_unmodified_index() {
                        assert_eq!(TestedCurve::from_index(0).0, 0);
                        assert_eq!(TestedCurve::from_index(0b10101).0, 0b10101);
                        assert_eq!(TestedCurve::from_index(INDEX_MAX).0, INDEX_MAX);
                    }

                    #[test]
                    fn index_retrieves_unmodified_index() {
                        assert_eq!(Morton::<$dm_macro!($t, $d), $d>(0).index(), 0);
                        assert_eq!(Morton::<$dm_macro!($t, $d), $d>(0b10101).index(), 0b10101);
                        assert_eq!(Morton::<$dm_macro!($t, $d), $d>(INDEX_MAX).index(), INDEX_MAX);
                    }

                    #[test]
                    #[allow(arithmetic_overflow)]
                    fn encoding_from_coords_too_large_panics() {
                        // Expand dilation methods don't have a max and won't panic
                        if COORD_MAX != <$t>::MAX {
                            for i in 0..$d {
                                let mut coords = [0; $d];
                                coords[i] = COORD_MAX + 1;

                                // We are testing each component of coords here, so can't easily use should_panic
                                // We'll emulate it at a basic level instead
                                let result = std::panic::catch_unwind(|| TestedCurve::from_coords(coords));
                                if !result.is_err() {
                                    panic!("Test did not panic as expected");
                                }
                            }
                        }
                    }

                    #[test]
                    fn encoding_from_coords_is_correct() {
                        for (coords, index) in generate_coords_test_cases(64) {
                            std::println!("Test case: coords = {:?}, index = 0b{:b}", coords, index);
                            std::println!("Result: 0b{:b}", TestedCurve::from_coords(coords).0);
                            assert_eq!(TestedCurve::from_coords(coords).0, index);
                        }
                    }

                    #[test]
                    fn encoding_coords_is_correct() {
                        for (coords, index) in generate_coords_test_cases(64) {
                            std::println!("Test case: coords = {:?}, index = 0b{:b}", coords, index);
                            std::println!("Result: {:?}", Morton::<$dm_macro!($t, $d), $d>(index).coords());

                            // Aww why can't I init the type alias as a tuple struct would?
                            assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).coords(), coords);
                        }
                    }

                    #[test]
                    fn sibling_on_axis_is_correct() {
                        for (coords, index) in generate_coords_test_cases(64) {
                            for axis in 0..$d {
                                let mut coords = coords;

                                // The sibling in this axis is either 1 above or 1 below the current coord
                                coords[axis] = if coords[axis] % 2 == 0 { coords[axis] + 1 } else { coords[axis] - 1 };
                                let expect_index = coords_to_index(coords);

                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).sibling_on_axis(axis).0, expect_index);
                            }
                        }
                    }

                    #[test]
                    fn sibling_or_same_on_axis_is_correct() {
                        for (coords, index) in generate_coords_test_cases(64) {
                            for axis in 0..$d {
                                // Negative search direction
                                let mut sib_coords = coords;
                                sib_coords[axis] = if sib_coords[axis] % 2 == 1 { sib_coords[axis] - 1 } else { sib_coords[axis] };
                                let expect_index = coords_to_index(sib_coords);
                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).sibling_or_same_on_axis(axis, QueryDirection::Negative).0, expect_index);

                                // Positive search direction
                                let mut sib_coords = coords;
                                sib_coords[axis] = if sib_coords[axis] % 2 == 0 { sib_coords[axis] + 1 } else { sib_coords[axis] };
                                let expect_index = coords_to_index(sib_coords);
                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).sibling_or_same_on_axis(axis, QueryDirection::Positive).0, expect_index);
                            }
                        }
                    }

                    #[test]
                    fn neighbour_on_axis_is_correct() {
                        for (coords, index) in generate_coords_test_cases(64) {
                            for axis in 0..$d {
                                // Negative search direction
                                let mut nei_coords = coords;
                                nei_coords[axis] = nei_coords[axis].wrapping_sub(1) & <TestedCurve as SpaceFillingCurve<$d>>::COORD_MAX;
                                let expect_index = coords_to_index(nei_coords);
                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).neighbour_on_axis(axis, QueryDirection::Negative).0, expect_index);

                                // Positive search direction
                                let mut nei_coords = coords;
                                nei_coords[axis] = nei_coords[axis].wrapping_add(1) & <TestedCurve as SpaceFillingCurve<$d>>::COORD_MAX;
                                let expect_index = coords_to_index(nei_coords);
                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).neighbour_on_axis(axis, QueryDirection::Positive).0, expect_index);
                            }
                        }
                    }
                }
            }
        )+}
    }

    use crate::internal::tests::test_curve;

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

    test_curve!(morton_expand, false, u8, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(morton_expand, false, u16, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(morton_expand, false, u32, 2, 3, 4);
    test_curve!(morton_expand, false, u64, 2);

    test_curve!(morton_fixed, false, u8, 2, 3, 4);
    test_curve!(morton_fixed, false, u16, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(morton_fixed, false, u32, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(morton_fixed, false, u64, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(morton_fixed, false, u128, 2, 3, 4, 5, 6, 7, 8);

//    test_morton_curve!(morton_expand, expand, u8, 2, 3, 4, 5, 6, 7, 8);
//    test_morton_curve!(morton_expand, expand, u16, 2, 3, 4, 5, 6, 7, 8);
//    test_morton_curve!(morton_expand, expand, u32, 2, 3, 4);
//    test_morton_curve!(morton_expand, expand, u64, 2);
//
//    test_morton_curve!(morton_fixed, fixed, u8, 2, 3, 4);
//    test_morton_curve!(morton_fixed, fixed, u16, 2, 3, 4, 5, 6, 7, 8);
//    test_morton_curve!(morton_fixed, fixed, u32, 2, 3, 4, 5, 6, 7, 8);
//    test_morton_curve!(morton_fixed, fixed, u64, 2, 3, 4, 5, 6, 7, 8);
//    test_morton_curve!(morton_fixed, fixed, u128, 2, 3, 4, 5, 6, 7, 8);
}
