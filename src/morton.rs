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

use super::{CurveCoord, CurveIndex, Neighbours, QueryDirection, Siblings, SpaceFillingCurve};

use crate::internal::{array_from_fn, NumTraits};
/*
struct PageMapping {
    page_bits: usize,
    data_bits: usize,
    data_count: usize,
    data_cursor: usize,
}

impl PageMapping {
    #[inline(always)]
    fn new<T, const N: usize>(page_bits: usize, data_bits: usize, data_count: usize) -> Self {
        Self { page_bits, data_bits, data_count, data_cursor: 0 }
    }

    #[inline(always)]
    fn load_page<F: Fn(usize) -> usize>(&mut self, get_bits: F) -> Option<usize> {
        let total_data_bits = self.data_bits * self.data_count;
        if self.data_cursor < total_data_bits {
            let num_to_copy = (total_data_bits - self.data_cursor).min(self.page_bits);
            let data_cursor_max = self.data_cursor + num_to_copy;

            let mut page = 0;
            while self.data_cursor < data_cursor_max {
                page |= get_bits()
            }
            self.data_cursor
            let page = get_page(self.load_cursor, self.load_cursor + num_to_copy);
            self.load_cursor += num_to_copy;
            Some(page)
        } else {
            None
        }
    }
}

struct PageFill {
    page_bits: usize,
    load_bits: usize,
    store_bits: usize,
    load_cursor: usize,
    store_cursor: usize,
}

impl PageFill {
    #[inline(always)]
    fn new(page_bits: usize, load_bits: usize, store_bits: usize) -> Self {
        Self { page_bits, load_bits, store_bits, load_cursor: 0, store_cursor: 0 }
    }

    #[inline(always)]
    fn load_page<F: Fn(usize, usize) -> usize>(&mut self, get_page: F) -> Option<usize> {
        if self.load_cursor < self.load_bits {
            let num_to_copy = (self.load_bits - self.load_cursor).min(self.page_bits);
            let page = get_page(self.load_cursor, self.load_cursor + num_to_copy);
            self.load_cursor += num_to_copy;
            Some(page)
        } else {
            None
        }
    }
}

fn test<DM: DilationMethod, const D: usize>(coords: [DM::Undilated; D])
where
    DM::Undilated: CurveCoord,
{
    let page_bits = core::mem::size_of::<usize>() / D;
    let coord_bits = DM::UNDILATED_BITS;
    let page_fill = PageFill::new(page_bits, coord_bits * D, 0);
    while let Some(page) = page_fill.load_page(|cursor_min, cursor_max| {
        let mut page = 0;
        let mut cursor = cursor_min;
        while cursor < cursor_max {
            let inner_cursor = cursor % coord_bits;
            let bits_remaining = coord_bits - inner_cursor;
            let mask = (0..bits_remaining).fold(0, |m, _| (m << 1) | 1);
            page |= coords[cursor / coord_bits].shr(inner_cursor).to_usize() & mask;
    
        }
    }) {
        let dilated = Fixed::<usize, D>::dilate(page);
    }
}
*/
/*
struct PageFill {
    num_data_bits: usize,
    num_page_bits: usize,
}

impl PageFill {
    pub const fn new(num_data_bits: usize, num_page_bits: usize) -> Self {
        Self {
            num_data_bits,
            num_page_bits,
        }
    }

    pub fn each_page(&self) {
        let mut page_index = 0;
        while page_index < self.num_pages() {
            let (data_index, data_cursor, mut total_num_to_copy) = self.page(page_index);
            while total_num_to_copy > 0 {
                let num_to_copy = total_num_to_copy.min()
            }
    
        }
    }



    const fn num_pages(&self) -> usize {
        self.num_data_bits / self.num_page_bits + 1
    }

    const fn page(&self, index: usize) -> (usize, usize, usize) {
        let shared_cursor = self.num_page_bits * index;
        let data_index = shared_cursor / self.num_data_bits;
        let data_cursor = shared_cursor % self.num_data_bits;
        (data_index, data_cursor, self.num_page_bits)
    }

    pub const fn page<F: (usize) -> usize>(&mut self, overlap: F) -> Option<usize> {
        if self.data_cursor < self.num_data_bits {
            let mut page = 0;
            let mut data_cursor = self.data_cursor;
            let mut page_cursor = 0;
            while data_cursor < self.num_data_bits && page_cursor < self.num_page_bits {
                let num_bits_to_copy = (self.num_data_bits - data_cursor).min(self.num_page_bits - page_cursor);
                let copy_mask = (1 << num_bits_to_copy) - 1;
                page |= (overlap(data_cursor) & copy_mask) << page_cursor;
                page_cursor += num_bits_to_copy;
                data_cursor += num_bits_to_copy;
            }
            Some(page)
        } else {
            None
        }
    }

    pub const fn has_page(&self) -> bool {
        self.data_cursor < self.num_data_bits
    }

    pub const fn has_overlap(&self) -> bool {
        self.page_cursor < self.num_page_bits && self.data_cursor < self.num_data_bits
    }*/

/*    #[inline(always)]
    const fn find_overlaps<F: FnMut(usize, usize, usize, usize)>(mut data_index: usize, mut data_cursor: usize, data_bits: usize, overlap: F) -> (usize, usize) {
        let mut page_cursor = 0;

        while page_cursor < PageBits && data_index < D {
            let num_bits_to_copy = (data_bits - data_cursor).min(PageBits - page_cursor);
            overlap(data_index, data_cursor, page_cursor, num_bits_to_copy);
            page_cursor += num_bits_to_copy;
            data_cursor += num_bits_to_copy;
            if data_cursor == data_bits {
                data_cursor = 0;
                data_index += 1;
            }
        }

        (data_index, data_cursor)
    }

    #[inline(always)]
    pub fn page(&self) -> usize {
        0
//        let mut page = 0;
//        if Self::find_overlaps(self.data_index, self.data_cursor, InBits, |data_index, data_cursor, page_cursor, num_bits_to_copy| {
//            let copy_mask = (0..num_bits_to_copy).fold(0, |m| (m << 1) | 1);
//            page |= ((self.data[data_index].shr(data_cursor).to_usize() & copy_mask) << page_cursor);
//        }) != (self.data_index, self.data_cursor) {
//            Some(page)
//        } else {
//            None
//        }
    }

    #[inline(always)]
    pub const fn next(mut self) -> Option<Self> {
        if self.first_iteration {
            self.first_iteration = false;
            Some(self)
        } else {
            let (data_index, data_cursor) = Self::find_overlaps(self.data_index, self.data_cursor, InBits, |_, _, _, _| { });
            if (data_index, data_cursor) != (self.data_index, self.data_cursor) {
                Some(Self { data: self.data, data_index, data_cursor, first_iteration: false })
            } else {
                None
            }
        }
    }*/
//}
/*
trait MortonDilate: Sized {
    fn morton_dilate_d2(coords: [Self; 2]) -> Self;
    fn morton_dilate_d3(coords: [Self; 3]) -> Self;
    fn morton_dilate_d4(coords: [Self; 4]) -> Self;
    fn morton_dilate_dn<const D: usize>(coords: [Self; D]) -> Self;
}

impl MortonDilate for u8 {
    fn morton_dilate_d2(coords: [Self; 2]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 2>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[1] as u16) << 4) | (coords[0] as u16);
        v = v.dilate_fixed::<2>().value();
        ((v >> 7) | v) as Self
    }

    fn morton_dilate_d3(coords: [Self; 3]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 3>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[2] as u32) << 4) | ((coords[1] as u32) << 2) | (coords[0] as u32);
        v = v.dilate_fixed::<3>().value();
        (((v >> 10) | (v >> 5) | v) as Self) & 0x3F
    }

    fn morton_dilate_d4(coords: [Self; 4]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 4>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[3] as u32) << 6) | ((coords[2] as u32) << 4) | ((coords[1] as u32) << 2) | (coords[0] as u32);
        v = v.dilate_fixed::<4>().value();
        ((v >> 21) | (v >> 14) | (v >> 7) | v) as Self
    }

    fn morton_dilate_dn<const D: usize>(coords: [Self; D]) -> Self {
        unimplemented!()
    }
}

impl MortonDilate for u16 {
    fn morton_dilate_d2(coords: [Self; 2]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 2>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[1] as u32) << 8) | (coords[0] as u32);
        v = v.dilate_fixed::<2>().value();
        ((v >> 15) | v) as Self
    }

    fn morton_dilate_d3(coords: [Self; 3]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 3>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[2] as u64) << 10) | ((coords[1] as u64) << 5) | (coords[0] as u64);
        v = v.dilate_fixed::<3>().value();
        (((v >> 28) | (v >> 14) | v) as Self) & 0x7FFF
    }

    fn morton_dilate_d4(coords: [Self; 4]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 4>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[3] as u64) << 12) | ((coords[2] as u64) << 8) | ((coords[1] as u64) << 4) | (coords[0] as u64);
        v = v.dilate_fixed::<4>().value();
        ((v >> 45) | (v >> 30) | (v >> 15) | v) as Self
    }

    fn morton_dilate_dn<const D: usize>(coords: [Self; D]) -> Self {
        unimplemented!()
//        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, D>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
//        let mut v = (0..D).fold(0, |v, i| v | ((coords[i] as u16) << (i * Fixed::<Self, D>::UNDILATED_BITS)));
//
//        
//        usize::BITS / D
//        v = v.dilate_fixed::<D>().value();
//
//        todo!();
    }
}

impl MortonDilate for u32 {
    fn morton_dilate_d2(coords: [Self; 2]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 2>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[1] as u64) << 16) | (coords[0] as u64);
        v = v.dilate_fixed::<2>().value();
        ((v >> 31) | v) as Self
    }

    fn morton_dilate_d3(coords: [Self; 3]) -> Self {
        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 3>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
        let mut v = ((coords[2] as u64) << 20) | ((coords[1] as u64) << 10) | (coords[0] as u64);
        v = v.dilate_fixed::<3>().value();
        (((v >> 28) | (v >> 14) | v) as Self) & 0x7FFF
    }

    fn morton_dilate_d4(coords: [Self; 4]) -> Self {
        unimplemented!()
        // Compress bits into separate word size dilations
        // Dilate
        // Uncompress bits
        // Merge

//        let fill_state = PageFill::<instride, outstride>(coords);
//        while let Some(page) = fill_state.next_page() {
//            match page.bits_occupied {
//                1..8 => { /* 8bit dilate */},
//                9..16 => { /* 16bit dilate */},
//                #[cfg(any(target_pointer_width = "32", target_pointer_width = "64"))]
//                17..32 => { /* 32bit dilate */},
//                #[cfg(any(target_pointer_width = "64"))]
//                33..64 => { /* 64bit dilate */},
//                _ => unreachable!(),
//            }
//        }

//        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, 4>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
//        let mut v = ((coords[3] as u64) << 12) | ((coords[2] as u64) << 8) | ((coords[1] as u64) << 4) | (coords[0] as u64);
//        v = v.dilate_fixed::<4>().value();
//        ((v >> 45) | (v >> 30) | (v >> 15) | v) as Self
    }

    fn morton_dilate_dn<const D: usize>(coords: [Self; D]) -> Self {
        unimplemented!()
//        debug_assert!(coords.iter().fold(0, |m, c| m | *c) <= Fixed::<Self, D>::UNDILATED_MAX, "Attempting to dilate a value which exceeds maximum (See DilationMethod::UNDILATED_MAX)");
//        let mut v = (0..D).fold(0, |v, i| v | ((coords[i] as u16) << (i * Fixed::<Self, D>::UNDILATED_BITS)));
//
//        
//        usize::BITS / D
//        v = v.dilate_fixed::<D>().value();
//
//        todo!();
    }
}*/

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

//trait MortonMethod<const D: usize> {
//    type Coord;
//    type Index;
//    const D: usize;
//
//    fn coords_to_morton(coords: [Self::Coord; D]) -> Self::Index;
//    fn morton_to_coords(morton: Self::Index) -> [Self::Coord; D];
//}

//impl MortonMethod<2> for Fixed<u8, 2> {
//    type Coord = Self::Coord;
//    type Index = Self::Index;
//    const D: usize = 2;
//
//    fn coords_to_morton(coords: [Self::Coord; Self::D]) -> Self::Index {
//
//    }
//
//    fn morton_to_coords(morton: Self::Index) -> [Self::Coord; Self::D] {
//    }
//}
//
//macro_rules! impl_regular_morton_method {
//    ($dm:ty) => {
//        impl MortonMethod<<$dm>::D> for $dm {
//            type Coord = Self::Coord;
//            type Index = Self::Index;
//            const D: usize = Self::D;
//        
//            fn coords_to_morton(coords: [Self::Coord; Self::D]) -> Self::Index {
//        
//            }
//        
//            fn morton_to_coords(morton: Self::Index) -> [Self::Coord; Self::D] {
//            }
//        }
//    };
//}

const fn test<const N: usize>() -> usize {
    N
}

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

    #[inline]
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        debug_assert!(
            *coords.iter().max().unwrap() <= Self::COORD_MAX,
            "Parameter 'coords' contains a value which exceeds maximum"
        );

        let _ = test::<{DM::UNDILATED_BITS * D}>();

//        let page_fill = PageFill::new(DM::UNDILATED_BITS * D, Fixed::<usize, D>::UNDILATED_BITS);
//        while page_fill.has_page() {
//            
//        }

//        PageFill::<Self::Coord, DM::UNDILATED_BITS, DM::DILATED_BITS, Fixed::<usize, D>::UNDILATED_BITS, D>::from_data(data)

        

        Self(
            coords
                .into_iter()
                .enumerate()
                .fold(Self::Index::zero(), |v, (i, c)| {
                    v.bit_or(DM::dilate(c).value().shl(i))
                }),
        )
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
        let coord = match direction {
            QueryDirection::Positive => {
                if coord.value() < DM::DILATED_MAX {
                    Some(coord.add_one())
                } else {
                    None
                }
            }
            QueryDirection::Negative => {
                if coord.value() > NumTraits::zero() {
                    Some(coord.sub_one())
                } else {
                    None
                }
            }
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
