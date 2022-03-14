use std::{
    ops::{Add, BitAnd, BitOr, BitXor, Not, Shl, Shr, Sub, BitOrAssign},
};

use dilate::*;

use super::{SearchDirection, Encoding, Siblings, Neighbours};

pub trait MortonIndex:
    DilatableType
    + Add<Output = Self>
    + Sub<Output = Self>
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + BitOr<Output = Self>
    + BitOrAssign
    + BitAnd<Output = Self>
    + BitXor<Output = Self>
    + Not<Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
}

macro_rules! impl_morton_index {
    ($($t:ty),+) => {$(
        impl MortonIndex for $t {
            fn zero() -> Self {
                0
            }
            fn one() -> Self {
                1
            }
        }
    )+};
}
impl_morton_index!(u8, u16, u32, u64, u128, usize);

// Until we have complex generic constants, we have to pass D in here (needed by coords)
// Waiting on: https://github.com/rust-lang/rust/issues/76560
pub struct Morton<DM, const D: usize>(DM::Dilated) where DM: DilationMethod;

impl<DM, const D: usize> Encoding<D> for Morton<DM, D>
where
    DM: DilationMethod,
    DM::Undilated: Ord,
    DM::Dilated: MortonIndex,
{
    type Coord = DM::Undilated;
    const COORD_MAX: Self::Coord = DM::UNDILATED_MAX;

    #[inline]
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        debug_assert!(*coords.iter().max().unwrap() <= Self::COORD_MAX, "Parameter 'coords' contains a value which exceeds maximum");
        let mut v = DM::Dilated::zero();
        for (i, c) in coords.into_iter().enumerate() {
            v |= DM::dilate(c).0 << i;
        }
        Self(v)
    }

    #[inline]
    fn coords(&self) -> [Self::Coord; D] {
        let mut i = 0;
        [(); D].map(|_| {
            let coord = DilatedInt::<DM>((self.0 >> i) & DM::DILATED_MAX).undilate();
            i += 1;
            coord
        })
    }
}

impl<DM, const D: usize> Siblings for Morton<DM, D>
where
    DM: DilationMethod,
    DM::Dilated: MortonIndex,
{
    #[inline]
    fn sibling_on_axis(&self, axis: usize) -> Self {
        Self(self.0 ^ (DM::Dilated::one() << axis))
    }

    #[inline]
    fn sibling_or_self_on_axis(&self, axis: usize, search_direction: SearchDirection) -> Self {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");
        let lower_axis_mask = DM::Dilated::one() << axis;
        let search_mask = match search_direction {
            SearchDirection::Positive => lower_axis_mask,
            SearchDirection::Negative => DM::Dilated::zero(),
        };
        Self(self.0 & !lower_axis_mask | search_mask)
    }
}

impl<DM, const D: usize> Neighbours for Morton<DM, D>
where
    DM: DilationMethod,
    DM::Dilated: MortonIndex,
    DilatedInt<DM>: AddOne + SubOne,
{
    #[inline]
    fn neighbour_on_axis(&self, axis: usize, search_direction: SearchDirection) -> Self {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");
        // This needs proving
        let coord = DilatedInt::<DM>((self.0 >> axis) & DM::DILATED_MAX);
        let coord = match search_direction {
            SearchDirection::Positive => coord.add_one(),
            SearchDirection::Negative => coord.sub_one(),
        };
        let index = self.0 & !(DM::DILATED_MAX << axis);
        Self(index | (coord.0 << axis))

        //let mut coords = Self::coords(index);
        //coords[axis] = match search_direction {
        //    SearchDirection::Positive => coords[axis] + DM::Undilated::one(),
        //    SearchDirection::Negative => coords[axis] - DM::Undilated::one(),
        //};
        //Self::from_coords(coords)
    }
}

#[cfg(test)]
mod tests {
    use std::panic::catch_unwind;

    use dilate::*;
    use serde_json::{Value, Map};
    use crate::test_data::*;

    lazy_static::lazy_static! {
//        pub static ref COORDS_TEST_DATA: serde_json::Value = import("morton_curve", "coords.json");
    }

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


//    fn test_cases(path: &str) -> &Vec<Value> {
//        println!("Using test cases: {path}");
//        COORDS_TEST_DATA.pointer(path).expect(format!("Could not find {path}").as_str()).as_array().expect(format!("Expected {path} to be an array").as_str())
//    }
//
//    fn test_case_coords<const N: usize>(test_case: &Map<String, Value>, coord_max: u128, index_max: u128) -> [u128; N] {
//        str_num_array_to_array::<N>(test_case.get("coords").expect("Expected coords array in test case"), coord_max, index_max)
//    }
//
//    fn test_case_index(test_case: &Map<String, Value>, coord_max: u128, index_max: u128) -> u128 {
//        str_num(test_case.get("index").expect("Expected index in test case"), coord_max, index_max)
//    }

    macro_rules! test_morton_curve {
        ($name:path, $dm_macro:ty, $t:ty, $($d:literal),+) => {$(
            paste::paste!{
                mod [< $name _ $t _d $d >] {
                    use super::*;
                    use super::super::*;

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

                            // Convert coords to morton index
                            // This is a low performance, but easy to read way to make a morton index
                            let mut index = 0;
                            for bit in 0..COORD_BITS {
                                for i in 0..$d {
                                    let coord_bit = ((coords[i] >> bit) & 0x1) as Index;
                                    let shift = bit * $d + i;
                                    index |= coord_bit << shift;
                                }
                            }

                            test_cases.push((coords, index));
                        }

                        // Also test 0 and max
                        test_cases.push(([0; $d], 0));
                        test_cases.push(([COORD_MAX; $d], INDEX_MAX));

                        test_cases
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
                        for t in generate_coords_test_cases(64) {
                            println!("Test case: coords = {:?}, index = 0b{:b}", t.0, t.1);
                            println!("Result: 0b{:b}", TestedCurve::from_coords(t.0).0);
                            assert_eq!(TestedCurve::from_coords(t.0).0, t.1);
                        }
                    }

                    #[test]
                    fn encoding_coords_is_correct() {
                        for t in generate_coords_test_cases(64) {
                            println!("Test case: coords = {:?}, index = 0b{:b}", t.0, t.1);
                            println!("Result: {:?}", Morton::<$dm_macro!($t, $d), $d>(t.1).coords());

                            // Aww why can't I init the type alias as a tuple struct would?
                            assert_eq!(Morton::<$dm_macro!($t, $d), $d>(t.1).coords(), t.0);
                        }
                    }
                }
            }
        )+}
    }

    test_morton_curve!(morton_expand, expand, u8, 1, 2, 3, 4, 5, 6, 7, 8);
    test_morton_curve!(morton_expand, expand, u16, 1, 2, 3, 4, 5, 6, 7, 8);
    test_morton_curve!(morton_expand, expand, u32, 1, 2, 3, 4);
    test_morton_curve!(morton_expand, expand, u64, 1, 2);

    test_morton_curve!(morton_fixed, fixed, u8, 1, 2, 3, 4);
    test_morton_curve!(morton_fixed, fixed, u16, 1, 2, 3, 4, 5, 6, 7, 8);
    test_morton_curve!(morton_fixed, fixed, u32, 1, 2, 3, 4, 5, 6, 7, 8);
    test_morton_curve!(morton_fixed, fixed, u64, 1, 2, 3, 4, 5, 6, 7, 8);
    test_morton_curve!(morton_fixed, fixed, u128, 1, 2, 3, 4, 5, 6, 7, 8);
}
