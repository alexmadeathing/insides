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

// It's much faster to process numbers using the word size of the processor
// This trait helps convert between outer type and a suitable word size inner type
pub trait InnerRep {
    type Inner: NumTraits;

    fn from_inner(value: Self::Inner) -> Self;
    fn to_inner(self) -> Self::Inner;
}

macro_rules! impl_inner_register {
    ($t:ty, $reg:ty) => {
        impl InnerRep for $t {
            type Inner = $reg;

            #[inline(always)]
            fn from_inner(value: Self::Inner) -> Self { value as Self }
        
            #[inline(always)]
            fn to_inner(self) -> Self::Inner { self as Self::Inner }
        }
    };
}

impl_inner_register!(u8, usize);

impl_inner_register!(u16, usize);

#[cfg(target_pointer_width = "16")]
impl_inner_register!(u32, u32);

#[cfg(any(target_pointer_width = "32", target_pointer_width = "64"))]
impl_inner_register!(u32, usize);

#[cfg(any(target_pointer_width = "16", target_pointer_width = "32"))]
impl_inner_register!(u64, u64);

#[cfg(target_pointer_width = "64")]
impl_inner_register!(u64, usize);

impl_inner_register!(u128, u128);

impl_inner_register!(usize, usize);

pub trait NumTraits: InnerRep + Copy + Ord {
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
    use crate::{Siblings, Neighbours, SpaceFillingCurve, QueryDirection};
    use core::{hash::Hash, panic::RefUnwindSafe};
    use std::{collections::HashSet, fmt::Debug, marker::PhantomData, panic::catch_unwind};

    const MAX_TESTED_INDICES: usize = 100000;

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
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES);

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

        pub fn sibling_on_axis_toggle_is_correct()
        where
            SFC: Siblings<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / D);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for axis in 0..D {
                    let mut expected_coords = coords;
                    expected_coords[axis] =
                        if coords[axis].bit_and(NumTraits::one()) == NumTraits::zero() {
                            coords[axis].add(NumTraits::one())
                        } else {
                            coords[axis].sub(NumTraits::one())
                        };
                    assert_eq!(sfc.sibling_on_axis_toggle(axis).coords(), expected_coords);
                }
            }
        }

        pub fn sibling_on_axis_is_correct()
        where
            SFC: Siblings<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / D);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for axis in 0..D {
                    let mut expected_coords = coords;
                    expected_coords[axis] =
                        if coords[axis].bit_and(NumTraits::one()) == NumTraits::zero() {
                            coords[axis].add(NumTraits::one())
                        } else {
                            coords[axis]
                        };
                    assert_eq!(sfc.sibling_on_axis(axis, QueryDirection::Positive).coords(), expected_coords);

                    expected_coords[axis] =
                        if coords[axis].bit_and(NumTraits::one()) == NumTraits::one() {
                            coords[axis].sub(NumTraits::one())
                        } else {
                            coords[axis]
                        };
                    assert_eq!(sfc.sibling_on_axis(axis, QueryDirection::Negative).coords(), expected_coords);
                }
            }
        }

        pub fn sibling_from_bits_is_correct()
        where
            SFC: Siblings<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_siblings = 1 << D;
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / num_siblings);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for n in 0..num_siblings {
                    let expected_coords = array_from_fn::<_, _, D>(|i| coords[i].bit_and(SFC::Coord::one().bit_not()).bit_or(NumTraits::from_usize((n >> i) & 0x1)));
                    assert_eq!(sfc.sibling_from_bits(NumTraits::from_usize(n)).coords(), expected_coords)
                }
            }
        }

        pub fn sibling_on_axes_is_correct()
        where
            SFC: Siblings<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_siblings = 1 << D;
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / num_siblings);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for n in 0..num_siblings {
                    let axes = array_from_fn::<_, _, D>(|i| if n >> i & 0x1 != 0 { QueryDirection::Positive } else { QueryDirection::Negative });
                    let expected_coords = array_from_fn::<_, _, D>(|i| coords[i].bit_and(SFC::Coord::one().bit_not()).bit_or(NumTraits::from_usize((n >> i) & 0x1)));
                    assert_eq!(sfc.sibling_on_axes(axes).coords(), expected_coords)
                }
            }
        }

        pub fn neighbour_on_axis_is_correct()
        where
            SFC: Neighbours<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / D);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for axis in 0..D {
                    let expected = if coords[axis] != SFC::COORD_MAX {
                        let mut nbr_coords = coords;
                        nbr_coords[axis] = coords[axis].add(NumTraits::one());
                        Some(SFC::from_coords(nbr_coords))
                    } else {
                        None
                    };
                    assert_eq!(sfc.neighbour_on_axis(axis, QueryDirection::Positive), expected);

                    let expected = if coords[axis] != NumTraits::zero() {
                        let mut nbr_coords = coords;
                        nbr_coords[axis] = coords[axis].sub(NumTraits::one());
                        Some(SFC::from_coords(nbr_coords))
                    } else {
                        None
                    };
                    assert_eq!(sfc.neighbour_on_axis(axis, QueryDirection::Negative), expected);
                }
            }
        }

        pub fn neighbour_on_axis_wrapping_is_correct()
        where
            SFC: Neighbours<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / D);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();
                for axis in 0..D {
                    let mut expected_coords = coords;
                    expected_coords[axis] =
                        if coords[axis] != SFC::COORD_MAX {
                            coords[axis].add(NumTraits::one())
                        } else {
                            NumTraits::zero()
                        };
                    assert_eq!(sfc.neighbour_on_axis_wrapping(axis, QueryDirection::Positive).coords(), expected_coords);

                    expected_coords[axis] =
                        if coords[axis] != NumTraits::zero() {
                            coords[axis].sub(NumTraits::one())
                        } else {
                            SFC::COORD_MAX
                        };
                    assert_eq!(sfc.neighbour_on_axis_wrapping(axis, QueryDirection::Negative).coords(), expected_coords);
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
                    fn sibling_on_axis_toggle_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::sibling_on_axis_toggle_is_correct();
                    }

                    #[test]
                    fn sibling_on_axis_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::sibling_on_axis_is_correct();
                    }

                    #[test]
                    fn sibling_from_bits_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::sibling_from_bits_is_correct();
                    }

                    #[test]
                    fn sibling_on_axes_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::sibling_on_axes_is_correct();
                    }
                }
            }
        )+}
    }
    pub(crate) use test_curve_siblings;

    macro_rules! test_curve_neighbours {
        ($curve:path, $require_adjacent:literal, $t:ty, $($d:literal),+) => {$(
            paste::paste!{
                mod [< $curve _neighbours_ $t _d $d >] {
                    use super::super::*;
                    use crate::internal::tests::SpaceFillingCurveTester;

                    type TestedCurve = $curve!($t, $d);

                    #[test]
                    fn neighbour_on_axis_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::neighbour_on_axis_is_correct();
                    }

                    #[test]
                    fn neighbour_on_axis_wrapping_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::neighbour_on_axis_wrapping_is_correct();
                    }
                }
            }
        )+}
    }
    pub(crate) use test_curve_neighbours;

    use super::array_from_fn;
}
