// Can be replaced with core::array::from_fn when stabilised
// https://github.com/rust-lang/rust/pull/94119
#[inline]
pub fn array_from_fn<F, T, const N: usize>(mut cb: F) -> [T; N]
where
    F: FnMut(usize) -> T,
{
    let mut idx = 0;
    [(); N].map(|_| {
        let res = cb(idx);
        idx += 1;
        res
    })
}

pub trait NumTraits: dilate::DilatableType + Ord {
    // Add methods as needed
    fn zero() -> Self;
    fn one() -> Self;
    fn bits() -> usize;
    fn lz(self) -> usize;
    fn add(self, rhs: Self) -> Self;
    fn sub(self, rhs: Self) -> Self;
    fn shl(self, amount: usize) -> Self;
    fn shr(self, amount: usize) -> Self;
    fn bit_not(self) -> Self;
    fn bit_and(self, rhs: Self) -> Self;
    fn bit_or(self, rhs: Self) -> Self;
    fn bit_xor(self, rhs: Self) -> Self;
    fn from_u8(value: u8) -> Self;
    fn from_usize(value: usize) -> Self;
    fn to_usize(self) -> usize;

    #[cfg(test)]
    fn max_value() -> Self;

    #[cfg(test)]
    fn distance_sqr(self, rhs: Self) -> usize;
}

macro_rules! impl_num_traits {
    ($($t:ty),+) => {$(
        impl NumTraits for $t {
            #[inline]
            fn zero() -> Self {
                0
            }

            #[inline]
            fn one() -> Self {
                1
            }

            #[inline]
            fn bits() -> usize {
                Self::BITS as usize
            }

            #[inline]
            fn lz(self) -> usize {
                self.leading_zeros() as usize
            }

            #[inline]
            fn add(self, rhs: Self) -> Self {
                self + rhs
            }

            #[inline]
            fn sub(self, rhs: Self) -> Self {
                self - rhs
            }

            #[inline]
            fn shl(self, amount: usize) -> Self {
                self << amount
            }

            #[inline]
            fn shr(self, amount: usize) -> Self {
                self >> amount
            }

            #[inline]
            fn bit_not(self) -> Self {
                !self
            }

            #[inline]
            fn bit_and(self, rhs: Self) -> Self {
                self & rhs
            }

            #[inline]
            fn bit_or(self, rhs: Self) -> Self {
                self | rhs
            }

            #[inline]
            fn bit_xor(self, rhs: Self) -> Self {
                self ^ rhs
            }

            #[inline]
            fn from_u8(value: u8) -> Self {
                value as Self
            }

            #[inline]
            fn from_usize(value: usize) -> Self {
                value as Self
            }

            #[inline]
            fn to_usize(self) -> usize {
                self as usize
            }

            #[inline]
            #[cfg(test)]
            fn max_value() -> Self {
                <$t>::MAX
            }

            #[inline]
            #[cfg(test)]
            fn distance_sqr(self, rhs: Self) -> usize {
                let a = self as i32;
                let b = rhs as i32;
                let dist = b - a;
                (dist * dist) as usize
            }
        }
    )+};
}

impl_num_traits!(u8, u16, u32, u64, u128, usize);

#[cfg(test)]
pub(crate) mod tests {
    extern crate std;

    use crate::internal::NumTraits;
    use crate::{Siblings, SpaceFillingCurve};
    use core::{hash::Hash, panic::RefUnwindSafe};
    use std::{collections::HashSet, fmt::Debug, marker::PhantomData, panic::catch_unwind};

    pub struct SpaceFillingCurveTester<SFC, const D: usize>(PhantomData<SFC>)
    where
        SFC: SpaceFillingCurve<D>;

    impl<SFC, const D: usize> SpaceFillingCurveTester<SFC, D>
    where
        SFC: SpaceFillingCurve<D>,
    {
        pub fn from_index_too_large_panics() {
            if SFC::INDEX_MAX != NumTraits::max_value() {
                SFC::from_index(SFC::INDEX_MAX.add(NumTraits::one()));
            } else {
                panic!("Parameter 'index' exceeds maximum");
            }
        }

        pub fn from_index_stores_unmodified_index()
        where
            <SFC as SpaceFillingCurve<D>>::Index: Debug,
        {
            assert_eq!(
                SFC::from_index(SFC::Index::zero()).index(),
                SFC::Index::zero()
            );
            assert_eq!(
                SFC::from_index(SFC::Index::from_usize(0b10101)).index(),
                SFC::Index::from_usize(0b10101)
            );
            assert_eq!(SFC::from_index(SFC::INDEX_MAX).index(), SFC::INDEX_MAX);
        }

        pub fn from_coords_too_large_panics()
        where
            <SFC as SpaceFillingCurve<D>>::Coord: RefUnwindSafe,
        {
            // Expand dilation methods don't have a max and won't panic
            if SFC::COORD_MAX != NumTraits::max_value() {
                for i in 0..D {
                    let mut coords = [NumTraits::zero(); D];
                    coords[i] = SFC::COORD_MAX.add(NumTraits::one());
                    // We are testing each component of coords here, so can't easily use should_panic
                    // We'll emulate it at a basic level instead
                    let result = catch_unwind(|| SFC::from_coords(coords));
                    if !result.is_err() {
                        panic!("Test did not panic as expected");
                    }
                }
            }
        }

        pub fn from_and_to_coords_are_correct(require_adjacent: bool)
        where
            <SFC as SpaceFillingCurve<D>>::Coord: Hash + Debug,
        {
            let max_num_indices = 100000;
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(max_num_indices);

            let mut coord_visit = HashSet::new();
            let mut prev_coords: Option<[SFC::Coord; D]> = None;

            for i in 0..num_indices {
                // Transform index to coords
                let coords = SFC::from_index(NumTraits::from_usize(i)).coords();

                // Record visits to this coord (should never be hit twice)
                assert_eq!(coord_visit.contains(&coords), false);
                coord_visit.insert(coords);

                // Optionally compare distance_sqr to previous coord (should always be 1)
                if require_adjacent {
                    if let Some(prev_coords) = prev_coords {
                        let distance_sqr = coords
                            .iter()
                            .zip(prev_coords.iter())
                            .fold(0, |sum, (&a, &b)| sum.add(b.distance_sqr(a)));
                        assert_eq!(distance_sqr, 1);
                    }
                    prev_coords = Some(coords);
                }

                // Transform coords to index (should be the original index)
                assert_eq!(SFC::from_coords(coords).index().to_usize(), i);
            }
        }

        pub fn sibling_on_axis_is_correct()
        where
            SFC: Siblings<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let max_num_indices = 100000;
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(max_num_indices);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for axis in 0..D {
                    let mut sib_coords = coords;
                    sib_coords[axis] =
                        if sib_coords[axis].bit_and(NumTraits::one()) == NumTraits::zero() {
                            sib_coords[axis].add(NumTraits::one())
                        } else {
                            sib_coords[axis].sub(NumTraits::one())
                        };
                    if sfc.sibling_on_axis(axis).coords() != sib_coords {
                        //std::println!("{:?}", );
                    }
                    assert_eq!(sfc.sibling_on_axis(axis).coords(), sib_coords);
                }
            }
        }
    }

    macro_rules! test_curve {
        ($curve:path, $require_adjacent:literal, $t:ty, $($d:literal),+) => {$(
            paste::paste!{
                mod [< $curve _ $t _d $d >] {
                    use super::super::*;
                    use crate::internal::tests::SpaceFillingCurveTester;

                    type TestedCurve = $curve!($t, $d);

                    #[test]
                    #[should_panic(expected = "Parameter 'index' exceeds maximum")]
                    fn from_index_too_large_panics() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::from_index_too_large_panics();
                    }

                    #[test]
                    fn from_index_stores_unmodified_index() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::from_index_stores_unmodified_index();
                    }

                    #[test]
                    #[allow(arithmetic_overflow)]
                    fn from_coords_too_large_panics() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::from_coords_too_large_panics();
                    }

                    #[test]
                    fn from_and_to_coords_are_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::from_and_to_coords_are_correct($require_adjacent);
                    }
                }
            }
        )+}
    }
    pub(crate) use test_curve;

    macro_rules! test_curve_siblings {
        ($curve:path, $require_adjacent:literal, $t:ty, $($d:literal),+) => {$(
            paste::paste!{
                mod [< $curve _siblings_ $t _d $d >] {
                    use super::super::*;
                    use crate::internal::tests::SpaceFillingCurveTester;

                    type TestedCurve = $curve!($t, $d);
                    #[test]
                    fn sibling_on_axis_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::sibling_on_axis_is_correct();
                    }

                    #[test]
                    fn sibling_or_same_on_axis_is_correct() {
//                        for (coords, index) in generate_coords_test_cases(64) {
//                            for axis in 0..$d {
//                                // Negative search direction
//                                let mut sib_coords = coords;
//                                sib_coords[axis] = if sib_coords[axis] % 2 == 1 { sib_coords[axis] - 1 } else { sib_coords[axis] };
//                                let expect_index = coords_to_index(sib_coords);
//                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).sibling_or_same_on_axis(axis, QueryDirection::Negative).0, expect_index);
//
//                                // Positive search direction
//                                let mut sib_coords = coords;
//                                sib_coords[axis] = if sib_coords[axis] % 2 == 0 { sib_coords[axis] + 1 } else { sib_coords[axis] };
//                                let expect_index = coords_to_index(sib_coords);
//                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).sibling_or_same_on_axis(axis, QueryDirection::Positive).0, expect_index);
//                            }
//                        }
                    }

                    #[test]
                    fn neighbour_on_axis_is_correct() {
//                        for (coords, index) in generate_coords_test_cases(64) {
//                            for axis in 0..$d {
//                                // Negative search direction
//                                let mut nei_coords = coords;
//                                nei_coords[axis] = nei_coords[axis].wrapping_sub(1) & <TestedCurve as SpaceFillingCurve<$d>>::COORD_MAX;
//                                let expect_index = coords_to_index(nei_coords);
//                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).neighbour_on_axis(axis, QueryDirection::Negative).0, expect_index);
//
//                                // Positive search direction
//                                let mut nei_coords = coords;
//                                nei_coords[axis] = nei_coords[axis].wrapping_add(1) & <TestedCurve as SpaceFillingCurve<$d>>::COORD_MAX;
//                                let expect_index = coords_to_index(nei_coords);
//                                assert_eq!(Morton::<$dm_macro!($t, $d), $d>(index).neighbour_on_axis(axis, QueryDirection::Positive).0, expect_index);
//                            }
//                        }
                    }
                }
            }
        )+}
    }
    pub(crate) use test_curve_siblings;
}
