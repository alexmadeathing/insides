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

use core::ops::Range;

use dilate::*;

use super::{Coord, Index, Encoding, Neighbours, QueryDirection, Siblings, Morton};

use crate::internal::NumTraits;

// Until we have complex generic constants, we have to pass D in here (needed by coords)
// Waiting on: https://github.com/rust-lang/rust/issues/76560
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hilbert<DM, const D: usize>(DM::Dilated)
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: Coord,
    DM::Dilated: Index;

impl<DM, const D: usize> Encoding<D> for Hilbert<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: Coord,
    DM::Dilated: Index,
    Morton<DM, D>: Encoding<D, Coord = DM::Undilated, Index = DM::Dilated>,
{
    type Coord = DM::Undilated;
    type Index = DM::Dilated;
    const COORD_MAX: Self::Coord = DM::UNDILATED_MAX;
    const INDEX_MAX: Self::Index = DM::DILATED_MASK;

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

        todo!();
    }

    #[inline]
    fn coords(&self) -> [Self::Coord; D] {
//        let rotate = |i, a| {  i.shl(a) }

        let mut morton_index = Self::Index::zero();
        let lower_mask = Self::Index::one().shl(D).sub(Self::Index::one());
//        let num_bits = (Self::Index::bits() - self.0.ls()) / D;
        for i in 0..DM::UNDILATED_BITS {
            let sub_index = self.0.shr(i * D).bit_and(lower_mask);
            let gray_index = sub_index.bit_xor(sub_index.shr(1));

            let rotate_amount = i.rem_euclid(D);
            let rotated_index = gray_index.shl(rotate_amount).bit_or(gray_index.shr(D - rotate_amount)).bit_and(lower_mask);

//            (i * D) + gray_index

            morton_index = morton_index.bit_or(rotated_index.shl(i * D));

        }

        Morton::<DM, D>::from_index(morton_index).coords()
    }

    #[inline]
    fn index(&self) -> Self::Index {
        self.0
    }
}
