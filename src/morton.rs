use dilate::*;

use super::{CurveCoord, CurveIndex, Neighbours, QueryDirection, Siblings, SpaceFillingCurve};

use crate::internal::{array_from_fn, NumTraits, MortonEncode};

/// A Morton encoded space filling curve implementation
///
/// Morton encoding, also known as a
/// [Z-order curve](https://en.wikipedia.org/wiki/Z-order_curve), is a fractal
/// space filling algorithm which maps a multidimensional set of coordinates to
/// one dimension and vice versa, achieved by interleaving the bit sequence of
/// each coordinate value.
///
/// Whilst other encoding methods may exhibit better spatial locality (such as
/// the Hilbert curve), the Morton curve offers excellent CPU performance,
/// since the behaviour can be reduced to a simple set of bitwise operations,
/// making it an ideal choice for CPU bound applications.
///
/// # Self Similarity
/// The Morton implementation featured in this crate is self-similar
/// in that lower magnitude indices will always reference the same
/// coordinates regardless how deep the curve is iterated. There is no
/// "depth" or "level" parameter.
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
    const COORD_BITS: usize = DM::UNDILATED_BITS;
    const COORD_MAX: Self::Coord = DM::UNDILATED_MAX;
    const INDEX_BITS: usize = DM::DILATED_BITS;
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

    #[inline(never)]
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        debug_assert!(
            *coords.iter().max().unwrap() <= Self::COORD_MAX,
            "Parameter 'coords' contains a value which exceeds maximum"
        );

        if D <= 4 {
            let coords = coords.map(|c| DM::to_dilated(c));
            Self(match D {
                2 => <Self::Index as MortonEncode<2>>::morton_encode([coords[0], coords[1]]),
                3 => <Self::Index as MortonEncode<3>>::morton_encode([coords[0], coords[1], coords[2]]),
                4 => <Self::Index as MortonEncode<4>>::morton_encode([coords[0], coords[1], coords[2], coords[3]]),
                _ => unreachable!(),
            })
        } else {
            Self(
                coords
                    .into_iter()
                    .enumerate()
                    .fold(Self::Index::zero(), |v, (i, c)| {
                        v.bit_or(DM::dilate(c).value().shl(i))
                    }),
            )
        }
    }

    #[inline]
    fn coords(&self) -> [Self::Coord; D] {
        array_from_fn::<_, _, D>(|i| {
            DilatedInt::<DM>::new(self.0.shr(i).bit_and(DM::DILATED_MAX)).undilate()
        })
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
    fn sibling_on_axis_toggle(&self, axis: usize) -> Self {
        Self(self.0.bit_xor(Self::Index::one().shl(axis)))
    }

    #[inline]
    fn sibling_on_axis(&self, axis: usize, direction: QueryDirection) -> Self {
        debug_assert!(axis < D, "Parameter 'axis' exceeds maximum");
        let lower_axis_mask = Self::Index::one().shl(axis);
        let search_mask = match direction {
            QueryDirection::Positive => lower_axis_mask,
            QueryDirection::Negative => NumTraits::zero(),
        };
        Self(
            self.0
                .bit_and(lower_axis_mask.bit_not())
                .bit_or(search_mask),
        )
    }

    #[inline]
    fn sibling_from_bits(&self, axis_bits: Self::Index) -> Self {
        let lower_mask = Self::Index::one().shl(D).sub(NumTraits::one());
        debug_assert!(
            axis_bits <= lower_mask,
            "Paremeter 'axis_bits' contains set bits in positions beyond D"
        );
        Self(self.0.bit_and(lower_mask.bit_not()).bit_or(axis_bits))
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
        match direction {
            QueryDirection::Positive => {
                if coord.value() < DM::DILATED_MAX {
                    let index = self.0.bit_and(DM::DILATED_MAX.shl(axis).bit_not());
                    Some(Self(index.bit_or(coord.add_one().value().shl(axis))))
                } else {
                    None
                }
            }
            QueryDirection::Negative => {
                if coord.value() > NumTraits::zero() {
                    let index = self.0.bit_and(DM::DILATED_MAX.shl(axis).bit_not());
                    Some(Self(index.bit_or(coord.sub_one().value().shl(axis))))
                } else {
                    None
                }
            }
        }
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

    fn neighbour_on_corner(&self, axes: [QueryDirection; D]) -> Option<Self> {
        // TODO can we improve on this? Seems like a lot of processing...
        let mut index = Self::Index::zero();
        for (axis, direction) in axes.into_iter().enumerate() {
            let coord = DilatedInt::<DM>::new(self.0.shr(axis).bit_and(DM::DILATED_MAX));
            let coord = match direction {
                QueryDirection::Positive => {
                    if coord.value() >= DM::DILATED_MAX {
                        return None;
                    }
                    coord.add_one()
                }
                QueryDirection::Negative => {
                    if coord.value() == NumTraits::zero() {
                        return None;
                    }
                    coord.sub_one()
                }
            };
            index = index.bit_or(coord.value().shl(axis));
        }
        Some(Self(index))
    }

    fn neighbour_on_corner_wrapping(&self, axes: [QueryDirection; D]) -> Self {
        // TODO Can we improve on this? Seems like a lot of processing...
        let mut index = Self::Index::zero();
        for (axis, direction) in axes.into_iter().enumerate() {
            let coord = DilatedInt::<DM>::new(self.0.shr(axis).bit_and(DM::DILATED_MAX));
            let coord = match direction {
                QueryDirection::Positive => coord.add_one(),
                QueryDirection::Negative => coord.sub_one(),
            };
            index = index.bit_or(coord.value().shl(axis));
        }
        Self(index)
    }
}

/// A Morton encoded space filling curve implementation using the Expand dilation method
// Not available until: https://github.com/rust-lang/rust/issues/112792
//pub type MortonExpand<T, const D: usize> = Morton<Expand<T, D>, D> where T: CurveIndex;

/// A Morton encoded space filling curve implementation using the Fixed dilation method
// Not available until: https://github.com/rust-lang/rust/issues/112792
//pub type MortonFixed<T, const D: usize> = Morton<Fixed<T, D>, D> where T: CurveIndex;

#[cfg(test)]
mod tests {
    use crate::internal::tests::{test_curve, test_curve_neighbours, test_curve_siblings};

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
