use std::{
    marker::PhantomData,
    ops::{Add, BitAnd, BitOr, BitXor, Not, Shl, Shr, Sub},
};

use dilate::*;

use super::{SearchDirection, Coords, Siblings, Neighbours};

pub trait MortonIndex:
    DilatableType
    + Add<Output = Self>
    + Sub<Output = Self>
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + BitOr<Output = Self>
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
pub struct MortonCurve<DM, const D: usize> where DM: DilationMethod {
    marker: PhantomData<DM>,
}

impl<DM, const D: usize> Coords<D> for MortonCurve<DM, D>
    where DM: DilationMethod,
    DM::Dilated: MortonIndex,
{
    type Coord = DM::Undilated;
    type Index = DM::Dilated;

    #[inline]
    fn from_coords(coords: [Self::Coord; D]) -> Self::Index {
        let mut v = Self::Index::zero();
        for (i, c) in coords.into_iter().enumerate() {
            v = v | (DM::dilate(c).0 << i);
        }
        v
    }

    #[inline]
    fn coords(index: Self::Index) -> [Self::Coord; D] {
        let mut i = 0;
        [(); D].map(|_| {
            let coord = DilatedInt::<DM>((index >> i) & DM::DILATED_MAX).undilate();
            i += 1;
            coord
        })
    }
}

impl<DM, const D: usize> Siblings for MortonCurve<DM, D>
where
    DM: DilationMethod,
    DM::Dilated: MortonIndex,
{
    type Index = DM::Dilated;

    #[inline]
    fn sibling_on_axis(index: DM::Dilated, axis: usize) -> DM::Dilated {
        index ^ (DM::Dilated::one() << axis)
    }

    #[inline]
    fn sibling_or_self_on_axis(index: DM::Dilated, axis: usize, search_direction: SearchDirection) -> DM::Dilated {
        let lower_axis_mask = DM::Dilated::one() << axis;
        let search_mask = match search_direction {
            SearchDirection::Positive => lower_axis_mask,
            SearchDirection::Negative => DM::Dilated::zero(),
        };
        index & !lower_axis_mask | search_mask
    }
}

impl<DM, const D: usize> Neighbours for MortonCurve<DM, D>
where
    DM: DilationMethod,
    DM::Dilated: MortonIndex,
    DilatedInt<DM>: AddOne + SubOne,
{
    type Index = DM::Dilated;

    #[inline]
    fn neighbour_on_axis(index: DM::Dilated, axis: usize, search_direction: SearchDirection) -> DM::Dilated {
        // This needs proving
        let coord = DilatedInt::<DM>((index >> axis) & DM::DILATED_MAX);
        let coord = match search_direction {
            SearchDirection::Positive => coord.add_one(),
            SearchDirection::Negative => coord.sub_one(),
        };
        let index = index & !(DM::DILATED_MAX << axis);
        index | (coord.0 << axis)

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
}
