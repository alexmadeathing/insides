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

use super::{Coord, Encoding, Index, Neighbours, QueryDirection, Siblings};

use crate::{internal::NumTraits, Morton};

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
        let lower_mask = Self::Index::one().shl(D).sub(Self::Index::one());

        let gray = |i: Self::Index| i.bit_xor(i.shr(1));

        let d = |i: usize| (i + ((i & 0x1) - 1)).trailing_ones() as usize;

        // TODO Can we remove this branch?
        let e = |i: Self::Index| {
            if i == NumTraits::zero() {
                NumTraits::zero()
            } else {
                gray(i.sub(NumTraits::one()).bit_and(Self::Index::one().bit_not()))
            }
        };

        let rotate =
            |i: Self::Index, a| lower_mask.bit_and(i.shl(a).bit_or(i.shr(D - a)));

        let order = (Self::Index::bits() - self.0.lz() + D - 1) / D;
        let skip_orders = DM::UNDILATED_BITS - order;

        let mut morton_index = Self::Index::zero();
        let mut flip_axes = Self::Index::zero();
        let mut rotate_amount = (skip_orders + 1) % D;
        for i in (0..order).rev() {
            let raw_index = self.0.shr(i * D).bit_and(lower_mask);
            let index = flip_axes.bit_xor(rotate(gray(raw_index), rotate_amount));
            morton_index = morton_index.bit_or(index.shl(i * D));
            flip_axes = flip_axes.bit_xor(rotate(e(raw_index), rotate_amount));
            rotate_amount = (d(raw_index.to_usize()) + rotate_amount + 1) % D;
        }
        Morton::<DM, D>::from_index(morton_index).coords()
    }

    #[inline]
    fn index(&self) -> Self::Index {
        self.0
    }
}
