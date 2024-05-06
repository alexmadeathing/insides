use core::array::from_fn;
use dilate::DilatableType;
use dilate::DilatedInt;
use dilate::Fixed;
use dilate::DilateFixed;

pub(crate) trait Sealed: DilatableType {
    // Main entrypoints
    fn morton_encode<const D: usize>(coords: [Self; D]) -> Self;
    fn morton_decode<const D: usize>(self, dilated_max: Self) -> [Self; D];

    // Utility methods
    fn _morton_encode_d2(coords: [Self; 2]) -> Self;
    fn _morton_encode_d3(coords: [Self; 3]) -> Self;
    fn _morton_encode_d4(coords: [Self; 4]) -> Self;
    fn _morton_decode_d2(self, dilated_max: Self) -> [Self; 2];
    fn _morton_decode_d3(self, dilated_max: Self) -> [Self; 3];
    fn _morton_decode_d4(self, dilated_max: Self) -> [Self; 4];
}

macro_rules! impl_sealed_common {
    () => {
        #[inline(always)]
        fn morton_encode<const D: usize>(coords: [Self; D]) -> Self {
            match D {
//                2 => self._morton_encode_d2([coords[0], coords[1]]),
//                3 => self._morton_encode_d3([coords[0], coords[1], coords[2]]),
//                4 => self._morton_encode_d4([coords[0], coords[1], coords[2], coords[3]]),
                _ => {
                    // TODO test performance against while loop
                    coords.into_iter().enumerate().fold(0, |v, (i, c)| {
                        v | c.dilate_fixed::<D>().value() << i
                    })
                },
            }
        }

        #[inline(always)]
        fn morton_decode<const D: usize>(self, dilated_max: Self) -> [Self; D] {
            match D {
//                2 => self._morton_decode_d2(),
//                3 => self._morton_decode_d3(),
//                4 => self._morton_decode_d4(),
                _ => {
                    from_fn::<_, D, _>(|i| {
                        DilatedInt::<Fixed<Self, D>>::new(self >> i & dilated_max).undilate()
                    })
                },
            }
        }
    };
}

impl Sealed for u8 {
    impl_sealed_common!();
}

pub trait NumTraits: Copy + Ord {
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
    fn fits_in_usize(self) -> bool;
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
            #[inline(always)]
            fn zero() -> Self {
                0
            }

            #[inline(always)]
            fn one() -> Self {
                1
            }

            #[inline(always)]
            fn bits() -> usize {
                Self::BITS as usize
            }

            #[inline(always)]
            fn lz(self) -> usize {
                self.leading_zeros() as usize
            }

            #[inline(always)]
            fn add(self, rhs: Self) -> Self {
                self + rhs
            }

            #[inline(always)]
            fn sub(self, rhs: Self) -> Self {
                self - rhs
            }

            #[inline(always)]
            fn shl(self, amount: usize) -> Self {
                self << amount
            }

            #[inline(always)]
            fn shr(self, amount: usize) -> Self {
                self >> amount
            }

            #[inline(always)]
            fn bit_not(self) -> Self {
                !self
            }

            #[inline(always)]
            fn bit_and(self, rhs: Self) -> Self {
                self & rhs
            }

            #[inline(always)]
            fn bit_or(self, rhs: Self) -> Self {
                self | rhs
            }

            #[inline(always)]
            fn bit_xor(self, rhs: Self) -> Self {
                self ^ rhs
            }

            #[inline(always)]
            fn fits_in_usize(self) -> bool {
                core::mem::size_of::<Self>() <= core::mem::size_of::<usize>() || self <= (usize::MAX as Self)
            }

            #[inline(always)]
            fn from_usize(value: usize) -> Self {
                value as Self
            }

            #[inline(always)]
            fn to_usize(self) -> usize {
                self as usize
            }

            #[inline(always)]
            #[cfg(test)]
            fn max_value() -> Self {
                <$t>::MAX
            }

            #[inline(always)]
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

pub trait MortonEncode<const D: usize>: Copy {
    fn morton_encode(coords: [Self; D]) -> Self;
}

impl MortonEncode<2> for u8 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 2]) -> Self {
        let dilated = ((coords[1] as u16) << 4 | coords[0] as u16).dilate_fixed::<2>().value();
        (dilated >> 7 | dilated) as u8
    }
}

impl MortonEncode<3> for u8 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 3]) -> Self {
        let dilated = ((coords[2] as u32) << 4 | (coords[1] as u32) << 2 | coords[0] as u32).dilate_fixed::<3>().value();
        (dilated >> 10 | dilated >> 5 | dilated) as u8 & 0x3F
    }
}

impl MortonEncode<4> for u8 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 4]) -> Self {
        let dilated = ((coords[3] as u32) << 6 | (coords[2] as u32) << 4 | (coords[1] as u32) << 2 | coords[0] as u32).dilate_fixed::<4>().value();
        (dilated >> 21 | dilated >> 14 | dilated >> 7 | dilated) as u8
    }
}

impl MortonEncode<2> for u16 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 2]) -> Self {
        let dilated = ((coords[1] as u32) << 8 | coords[0] as u32).dilate_fixed::<2>().value();
        (dilated >> 15 | dilated) as u16
    }
}

impl MortonEncode<3> for u16 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 3]) -> Self {
        // These coords are intentionally not packed contiguously so that we can avoid a bitmask when recombining
        let dilated = ((coords[2] as u64) << 12 | (coords[1] as u64) << 6 | coords[0] as u64).dilate_fixed::<3>().value();
        (dilated >> 34 | dilated >> 17 | dilated) as u16
    }
}

impl MortonEncode<4> for u16 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 4]) -> Self {
        let dilated = ((coords[3] as u64) << 12 | (coords[2] as u64) << 8 | (coords[1] as u64) << 4 | coords[0] as u64).dilate_fixed::<4>().value();
        (dilated >> 45 | dilated >> 30 | dilated >> 15 | dilated) as u16
    }
}

impl MortonEncode<2> for u32 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 2]) -> Self {
        let dilated = ((coords[1] as u64) << 16 | coords[0] as u64).dilate_fixed::<2>().value();
        (dilated >> 32 << 1 | dilated) as u32
    }
}

impl MortonEncode<3> for u32 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 3]) -> Self {
        // These coords are intentionally not packed contiguously so that we can avoid a bitmask when recombining
        let dilated_xy = ((coords[1] as u64) << 11 | coords[0] as u64).dilate_fixed::<3>().value();
        let result_xy = (dilated_xy >> 32 | dilated_xy) as u32;
    
        let dilated_z = coords[2].dilate_fixed::<3>().value();
        let result_z = dilated_z << 2;
    
        result_xy | result_z
    }
}

impl MortonEncode<4> for u32 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 4]) -> Self {
        let dilated_xy = ((coords[1] as u64) << 8 | coords[0] as u64).dilate_fixed::<4>().value();
        let result_xy = (dilated_xy >> 31 | dilated_xy) as u32;
    
        let dilated_zw = ((coords[3] as u64) << 8 | coords[2] as u64).dilate_fixed::<4>().value();
        let result_zw = (dilated_zw >> 29 | dilated_zw << 2) as u32;
    
        result_xy | result_zw
    }
}

impl MortonEncode<2> for u64 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 2]) -> Self {
        0
    }
}

impl MortonEncode<3> for u64 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 3]) -> Self {
        0
    }
}

impl MortonEncode<4> for u64 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 4]) -> Self {
        0
    }
}

impl MortonEncode<2> for u128 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 2]) -> Self {
        0
    }
}

impl MortonEncode<3> for u128 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 3]) -> Self {
        0
    }
}

impl MortonEncode<4> for u128 {
    #[inline(always)]
    fn morton_encode(coords: [Self; 4]) -> Self {
        0
    }
}

impl MortonEncode<2> for usize {
    #[inline(always)]
    fn morton_encode(coords: [Self; 2]) -> Self {
        #[cfg(target_pointer_width = "8")]
        return u8::morton_encode(coords.map(|c| c as u8)) as Self;
        #[cfg(target_pointer_width = "16")]
        return u16::morton_encode(coords.map(|c| c as u16)) as Self;
        #[cfg(target_pointer_width = "32")]
        return u32::morton_encode(coords.map(|c| c as u32)) as Self;
        #[cfg(target_pointer_width = "64")]
        return u64::morton_encode(coords.map(|c| c as u64)) as Self;
    }
}

impl MortonEncode<3> for usize {
    #[inline(always)]
    fn morton_encode(coords: [Self; 3]) -> Self {
        #[cfg(target_pointer_width = "8")]
        return u8::morton_encode(coords.map(|c| c as u8)) as Self;
        #[cfg(target_pointer_width = "16")]
        return u16::morton_encode(coords.map(|c| c as u16)) as Self;
        #[cfg(target_pointer_width = "32")]
        return u32::morton_encode(coords.map(|c| c as u32)) as Self;
        #[cfg(target_pointer_width = "64")]
        return u64::morton_encode(coords.map(|c| c as u64)) as Self;
    }
}

impl MortonEncode<4> for usize {
    #[inline(always)]
    fn morton_encode(coords: [Self; 4]) -> Self {
        #[cfg(target_pointer_width = "8")]
        return u8::morton_encode(coords.map(|c| c as u8)) as Self;
        #[cfg(target_pointer_width = "16")]
        return u16::morton_encode(coords.map(|c| c as u16)) as Self;
        #[cfg(target_pointer_width = "32")]
        return u32::morton_encode(coords.map(|c| c as u32)) as Self;
        #[cfg(target_pointer_width = "64")]
        return u64::morton_encode(coords.map(|c| c as u64)) as Self;
    }
}

#[cfg(test)]
pub(crate) mod tests {
    extern crate std;

    use crate::internal::NumTraits;
    use crate::{Neighbours, QueryDirection, Siblings, SpaceFillingCurve};
    use core::{hash::Hash, panic::RefUnwindSafe};
    use std::array::from_fn;
    use std::println;
    use std::{collections::HashSet, fmt::Debug, marker::PhantomData, panic::catch_unwind};

    const MAX_TESTED_INDICES: usize = 100000;
    const MAX_PERMUTATIONS: usize = 10000;

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
                if SFC::from_coords(coords).index().to_usize() != i {
                    println!("OH :(");
                }
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
                    let expected_coords = from_fn::<_, D, _>(|i| coords[i].bit_and(SFC::Coord::one().bit_not()).bit_or(NumTraits::from_usize((n >> i) & 0x1)));
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
                    let axes = from_fn::<_, D, _>(|i| if n >> i & 0x1 != 0 { QueryDirection::Positive } else { QueryDirection::Negative });
                    let expected_coords = from_fn::<_, D, _>(|i| coords[i].bit_and(SFC::Coord::one().bit_not()).bit_or(NumTraits::from_usize((n >> i) & 0x1)));
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
        
        pub fn neighbour_on_corner_is_correct()
        where
            SFC: Neighbours<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / D);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();

                let num_corners = (1usize << D).min(MAX_PERMUTATIONS);
                for corner in 0..num_corners {
                    let test_axes = core::array::from_fn(|i| if (corner >> i) & 0x1 != 0 { QueryDirection::Positive } else { QueryDirection::Negative });

                    let mut expected_coords = coords;
                    let mut valid = true;
                    for axis in 0..D {
                        if matches!(test_axes[axis], QueryDirection::Positive) {
                            if coords[axis] < SFC::COORD_MAX {
                                expected_coords[axis] = coords[axis].add(NumTraits::one());
                            } else {
                                valid = false;
                            }
                        } else {
                            if coords[axis] > NumTraits::zero() {
                                expected_coords[axis] = coords[axis].sub(NumTraits::one());
                            } else {
                                valid = false;
                            }
                        }
                    }
                    let expected = valid.then(|| SFC::from_coords(expected_coords));
                    assert_eq!(sfc.neighbour_on_corner(test_axes), expected);
                }
            }
        }
        
        pub fn neighbour_on_corner_wrapping_is_correct()
        where
            SFC: Neighbours<D>,
            <SFC as SpaceFillingCurve<D>>::Coord: Debug,
        {
            let num_indices: usize = SFC::INDEX_MAX.to_usize().min(MAX_TESTED_INDICES / D);
            for i in 0..num_indices {
                // It's a shame that we have to rely on other SFC methods to test this trait... not sure of a better solution yet
                let sfc = SFC::from_index(NumTraits::from_usize(i));
                let coords = sfc.coords();

                let num_corners = (1usize << D).min(MAX_PERMUTATIONS);
                for corner in 0..num_corners {
                    let test_axes = core::array::from_fn(|i| if (corner >> i) & 0x1 != 0 { QueryDirection::Positive } else { QueryDirection::Negative });

                    let mut expected_coords = coords;
                    for axis in 0..D {
                        if matches!(test_axes[axis], QueryDirection::Positive) {
                            if coords[axis] < SFC::COORD_MAX {
                                expected_coords[axis] = coords[axis].add(NumTraits::one());
                            } else {
                                expected_coords[axis] = NumTraits::zero();
                            }
                        } else {
                            if coords[axis] > NumTraits::zero() {
                                expected_coords[axis] = coords[axis].sub(NumTraits::one());
                            } else {
                                expected_coords[axis] = SFC::COORD_MAX;
                            }
                        }
                    }
                    let expected = SFC::from_coords(expected_coords);
                    assert_eq!(sfc.neighbour_on_corner_wrapping(test_axes), expected);
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

                    #[test]
                    fn neighbour_on_corner_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::neighbour_on_corner_is_correct();
                    }

                    #[test]
                    fn neighbour_on_corner_wrapping_is_correct() {
                        SpaceFillingCurveTester::<TestedCurve, $d>::neighbour_on_corner_wrapping_is_correct();
                    }
                }
            }
        )+}
    }
    pub(crate) use test_curve_neighbours;
}
