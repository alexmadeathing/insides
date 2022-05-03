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

use crate::{internal::NumTraits, Morton};

#[cfg(all(feature = "lut_small", not(feature = "lut_large")))]
mod lut_small {
    use super::NumTraits;

    pub const MAX_LUT_D: usize = 4;

    #[inline(always)]
    pub fn axis_tx_into_shl_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> T {
//        #[cfg(feature = "lut_d2")]
        if D == 2 {
            const AXIS_TX_INTO_SHL_LUT_D2: [[u8; 4]; 2] = [
                [00, 00, 00, 03, ],
                [00, 00, 00, 03, ],
            ];
            return T::from_u8(AXIS_TX_INTO_SHL_LUT_D2[rotate_amount][index.to_usize()]);
        }

//        #[cfg(feature = "lut_d3")]
        if D == 3 {
            const AXIS_TX_INTO_SHL_LUT_D3: [[u8; 8]; 3] = [
                [00, 00, 00, 03, 03, 06, 06, 05, ],
                [00, 00, 00, 06, 06, 05, 05, 03, ],
                [00, 00, 00, 05, 05, 03, 03, 06, ],
            ];
            return T::from_u8(AXIS_TX_INTO_SHL_LUT_D3[rotate_amount][index.to_usize()]);
        }

//        #[cfg(feature = "lut_d4")]
        if D == 4 {
            const AXIS_TX_INTO_SHL_LUT_D4: [[u8; 16]; 4] = [
                [00, 00, 00, 03, 03, 06, 06, 05, 05, 12, 12, 15, 15, 10, 10, 09, ],
                [00, 00, 00, 06, 06, 12, 12, 10, 10, 09, 09, 15, 15, 05, 05, 03, ],
                [00, 00, 00, 12, 12, 09, 09, 05, 05, 03, 03, 15, 15, 10, 10, 06, ],
                [00, 00, 00, 09, 09, 03, 03, 10, 10, 06, 06, 15, 15, 05, 05, 12, ],
            ];
            return T::from_u8(AXIS_TX_INTO_SHL_LUT_D4[rotate_amount][index.to_usize()]);
        }
        unimplemented!()
    }

    #[inline(always)]
    pub fn shr_into_gray_inv_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> T {
//        #[cfg(feature = "lut_d2")]
        if D == 2 {
            const SHR_INTO_GRAY_INV_LUT_D2: [[u8; 4]; 2] = [
                [00, 01, 03, 02, ],
                [00, 03, 01, 02, ],
            ];
            return T::from_u8(SHR_INTO_GRAY_INV_LUT_D2[rotate_amount][index.to_usize()]);
        }

//        #[cfg(feature = "lut_d3")]
        if D == 3 {
            const SHR_INTO_GRAY_INV_LUT_D3: [[u8; 8]; 3] = [
                [00, 01, 03, 02, 07, 06, 04, 05, ],
                [00, 07, 01, 06, 03, 04, 02, 05, ],
                [00, 03, 07, 04, 01, 02, 06, 05, ],
            ];
            return T::from_u8(SHR_INTO_GRAY_INV_LUT_D3[rotate_amount][index.to_usize()]);
        }

//        #[cfg(feature = "lut_d4")]
        if D == 4 {
            const SHR_INTO_GRAY_INV_LUT_D4: [[u8; 16]; 4] = [
                [00, 01, 03, 02, 07, 06, 04, 05, 15, 14, 12, 13, 08, 09, 11, 10, ],
                [00, 15, 01, 14, 03, 12, 02, 13, 07, 08, 06, 09, 04, 11, 05, 10, ],
                [00, 07, 15, 08, 01, 06, 14, 09, 03, 04, 12, 11, 02, 05, 13, 10, ],
                [00, 03, 07, 04, 15, 12, 08, 11, 01, 02, 06, 05, 14, 13, 09, 10, ],
            ];
            return T::from_u8(SHR_INTO_GRAY_INV_LUT_D4[rotate_amount][index.to_usize()]);
        }
        unimplemented!()
    }

    #[inline(always)]
    pub fn gray_into_shl_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> T {
//        #[cfg(feature = "lut_d2")]
        if D == 2 {
            const GRAY_INTO_SHL_LUT_D2: [[u8; 4]; 2] = [
                [00, 01, 03, 02, ],
                [00, 02, 03, 01, ],
            ];
            return T::from_u8(GRAY_INTO_SHL_LUT_D2[rotate_amount][index.to_usize()]);
        }

//        #[cfg(feature = "lut_d3")]
        if D == 3 {
            const GRAY_INTO_SHL_LUT_D3: [[u8; 8]; 3] = [
                [00, 01, 03, 02, 06, 07, 05, 04, ],
                [00, 02, 06, 04, 05, 07, 03, 01, ],
                [00, 04, 05, 01, 03, 07, 06, 02, ],
            ];
            return T::from_u8(GRAY_INTO_SHL_LUT_D3[rotate_amount][index.to_usize()]);
        }

//        #[cfg(feature = "lut_d4")]
        if D == 4 {
            const GRAY_INTO_SHL_LUT_D4: [[u8; 16]; 4] = [
                [00, 01, 03, 02, 06, 07, 05, 04, 12, 13, 15, 14, 10, 11, 09, 08, ],
                [00, 02, 06, 04, 12, 14, 10, 08, 09, 11, 15, 13, 05, 07, 03, 01, ],
                [00, 04, 12, 08, 09, 13, 05, 01, 03, 07, 15, 11, 10, 14, 06, 02, ],
                [00, 08, 09, 01, 03, 11, 10, 02, 06, 14, 15, 07, 05, 13, 12, 04, ],
            ];
            return T::from_u8(GRAY_INTO_SHL_LUT_D4[rotate_amount][index.to_usize()]);
        }
        unimplemented!()
    }

    #[inline(always)]
    pub fn rotation_increment_lut<T: NumTraits, const D: usize>(index: T, rotate_amount: usize) -> usize {
//        #[cfg(feature = "lut_d2")]
        if D == 2 {
            const ROTATION_INCREMENT_LUT_D2: [[u8; 4]; 2] = [
                [01, 00, 00, 01, ],
                [00, 01, 01, 00, ],
            ];
            return ROTATION_INCREMENT_LUT_D2[rotate_amount][index.to_usize()] as usize;
        }

//        #[cfg(feature = "lut_d3")]
        if D == 3 {
            const ROTATION_INCREMENT_LUT_D3: [[u8; 8]; 3] = [
                [01, 02, 02, 00, 00, 02, 02, 01, ],
                [02, 00, 00, 01, 01, 00, 00, 02, ],
                [00, 01, 01, 02, 02, 01, 01, 00, ],
            ];
            return ROTATION_INCREMENT_LUT_D3[rotate_amount][index.to_usize()] as usize;
        }

//        #[cfg(feature = "lut_d4")]
        if D == 4 {
            const ROTATION_INCREMENT_LUT_D4: [[u8; 16]; 4] = [
                [01, 02, 02, 03, 03, 02, 02, 00, 00, 02, 02, 03, 03, 02, 02, 01, ],
                [02, 03, 03, 00, 00, 03, 03, 01, 01, 03, 03, 00, 00, 03, 03, 02, ],
                [03, 00, 00, 01, 01, 00, 00, 02, 02, 00, 00, 01, 01, 00, 00, 03, ],
                [00, 01, 01, 02, 02, 01, 01, 03, 03, 01, 01, 02, 02, 01, 01, 00, ],
            ];
            return ROTATION_INCREMENT_LUT_D4[rotate_amount][index.to_usize()] as usize;
        }
        unimplemented!()
    }

    #[inline(always)]
    pub fn rotation_transform_lut(index: usize) -> usize {
        const ROTATION_TRANSFORM_LUT: [u8; 16] = [
            00, 01, 01, 02, 02, 01, 01, 03, 03, 01, 01, 02, 02, 01, 01, 04,
        ];

        ROTATION_TRANSFORM_LUT[index] as usize
    }

    #[inline(always)]
    pub fn walk_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        min_order: usize,
        mut f: F,
    ) -> (usize, T)
    where
        F: FnMut(T, usize, usize, T) -> T,
    {
        let lower_mask = T::one().shl(D).sub(T::one());

        // A subset of highest bits of index will be 0, therefore, we can skip a certain
        // number of iterations in the following loop. Note that skipping iterations
        // still requires us to rotate the pattern as if we had iterated all orders, see
        // rotate_amount below.
        let order = (T::bits() - index.lz() + D - 1) / D;

        // We iterate in reverse because higher orders affect the rotation of lower orders
        let mut flip_axes = T::zero();
        let mut rotate_amount = (D - 1) - (order % D);
        for shift in (min_order..order).rev().map(|i| i * D) {
            // Extract a portion of the index
            let index_partial = index.shr(shift).bit_and(lower_mask);

            // Context specific code generates the hilbert partial
            let hilbert_partial = f(index_partial, shift, rotate_amount, flip_axes);

            // Calculate the next axis flip flags and rotation amounts resulting from this order
            flip_axes = flip_axes.bit_xor(axis_tx_into_shl_lut::<_, D>(hilbert_partial, rotate_amount));
            rotate_amount = rotation_increment_lut::<_, D>(hilbert_partial, rotate_amount);
        }
        (rotate_amount, flip_axes)
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        debug_assert!(D <= MAX_LUT_D);
        let mut hilbert_index = T::zero();
        walk_transforms::<_, _, D>(
            morton_index,
            0,
            |morton_partial, shift, rotate_amount, flip_axes| {
                let hilbert_partial = shr_into_gray_inv_lut::<_, D>(flip_axes.bit_xor(morton_partial), rotate_amount);
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        debug_assert!(D <= MAX_LUT_D);
        let mut morton_index = T::zero();
        walk_transforms::<_, _, D>(
            hilbert_index,
            0,
            |hilbert_partial, shift, rotate_amount, flip_axes| {
                let morton_partial = flip_axes.bit_xor(gray_into_shl_lut::<_, D>(hilbert_partial, rotate_amount));
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
    }
}

#[cfg(feature = "lut_large")]
mod lut_large {
    use super::NumTraits;
    
    pub const MAX_LUT_D: usize = 4;

    const TRANSFORM_LUT_D2: [[u8; 4]; 4] = [
        [01, 00, 00, 03],
        [00, 01, 01, 02],
        [03, 02, 02, 01],
        [02, 03, 03, 00],
    ];
        
    const MORTON_LUT_D2: [[u8; 4]; 4] = [
        [00, 02, 03, 01],
        [00, 01, 03, 02],
        [03, 01, 00, 02],
        [03, 02, 00, 01],
    ];
    
    const HILBERT_LUT_D2: [[u8; 4]; 4] = [
        [00, 03, 01, 02],
        [00, 01, 03, 02],
        [02, 01, 03, 00],
        [02, 03, 01, 00],
    ];
    
    const TRANSFORM_LUT_D3: [[u8; 8]; 12] = [
        [02, 01, 01, 09, 09, 04, 04, 08],
        [00, 02, 02, 07, 07, 11, 11, 03],
        [01, 00, 00, 05, 05, 06, 06, 10],
        [05, 04, 04, 06, 06, 01, 01, 11],
        [03, 05, 05, 10, 10, 08, 08, 00],
        [04, 03, 03, 02, 02, 09, 09, 07],
        [08, 07, 07, 03, 03, 10, 10, 02],
        [06, 08, 08, 01, 01, 05, 05, 09],
        [07, 06, 06, 11, 11, 00, 00, 04],
        [11, 10, 10, 00, 00, 07, 07, 05],
        [09, 11, 11, 04, 04, 02, 02, 06],
        [10, 09, 09, 08, 08, 03, 03, 01],
    ];
    
    const MORTON_LUT_D3: [[u8; 8]; 12] = [
        [00, 04, 05, 01, 03, 07, 06, 02],
        [00, 02, 06, 04, 05, 07, 03, 01],
        [00, 01, 03, 02, 06, 07, 05, 04],
        [03, 07, 06, 02, 00, 04, 05, 01],
        [03, 01, 05, 07, 06, 04, 00, 02],
        [03, 02, 00, 01, 05, 04, 06, 07],
        [06, 02, 03, 07, 05, 01, 00, 04],
        [06, 04, 00, 02, 03, 01, 05, 07],
        [06, 07, 05, 04, 00, 01, 03, 02],
        [05, 01, 00, 04, 06, 02, 03, 07],
        [05, 07, 03, 01, 00, 02, 06, 04],
        [05, 04, 06, 07, 03, 02, 00, 01],
    ];
    
    const HILBERT_LUT_D3: [[u8; 8]; 12] = [
        [00, 03, 07, 04, 01, 02, 06, 05],
        [00, 07, 01, 06, 03, 04, 02, 05],
        [00, 01, 03, 02, 07, 06, 04, 05],
        [04, 07, 03, 00, 05, 06, 02, 01],
        [06, 01, 07, 00, 05, 02, 04, 03],
        [02, 03, 01, 00, 05, 04, 06, 07],
        [06, 05, 01, 02, 07, 04, 00, 03],
        [02, 05, 03, 04, 01, 06, 00, 07],
        [04, 05, 07, 06, 03, 02, 00, 01],
        [02, 01, 05, 06, 03, 00, 04, 07],
        [04, 03, 05, 02, 07, 00, 06, 01],
        [06, 07, 05, 04, 01, 00, 02, 03],
    ];
    
    const TRANSFORM_LUT_D4: [[u8; 16]; 32] = [
        [03, 02, 02, 29, 29, 06, 06, 24, 24, 10, 10, 21, 21, 14, 14, 19],
        [00, 03, 03, 18, 18, 31, 31, 13, 13, 07, 07, 22, 22, 27, 27, 08],
        [01, 00, 00, 11, 11, 16, 16, 26, 26, 28, 28, 23, 23, 12, 12, 05],
        [02, 01, 01, 04, 04, 09, 09, 15, 15, 17, 17, 20, 20, 25, 25, 30],
        [07, 06, 06, 25, 25, 02, 02, 28, 28, 14, 14, 17, 17, 10, 10, 23],
        [04, 07, 07, 22, 22, 27, 27, 09, 09, 03, 03, 18, 18, 31, 31, 12],
        [05, 04, 04, 15, 15, 20, 20, 30, 30, 24, 24, 19, 19, 08, 08, 01],
        [06, 05, 05, 00, 00, 13, 13, 11, 11, 21, 21, 16, 16, 29, 29, 26],
        [11, 10, 10, 21, 21, 14, 14, 16, 16, 02, 02, 29, 29, 06, 06, 27],
        [08, 11, 11, 26, 26, 23, 23, 05, 05, 15, 15, 30, 30, 19, 19, 00],
        [09, 08, 08, 03, 03, 24, 24, 18, 18, 20, 20, 31, 31, 04, 04, 13],
        [10, 09, 09, 12, 12, 01, 01, 07, 07, 25, 25, 28, 28, 17, 17, 22],
        [15, 14, 14, 17, 17, 10, 10, 20, 20, 06, 06, 25, 25, 02, 02, 31],
        [12, 15, 15, 30, 30, 19, 19, 01, 01, 11, 11, 26, 26, 23, 23, 04],
        [13, 12, 12, 07, 07, 28, 28, 22, 22, 16, 16, 27, 27, 00, 00, 09],
        [14, 13, 13, 08, 08, 05, 05, 03, 03, 29, 29, 24, 24, 21, 21, 18],
        [19, 18, 18, 13, 13, 22, 22, 08, 08, 26, 26, 05, 05, 30, 30, 03],
        [16, 19, 19, 02, 02, 15, 15, 29, 29, 23, 23, 06, 06, 11, 11, 24],
        [17, 16, 16, 27, 27, 00, 00, 10, 10, 12, 12, 07, 07, 28, 28, 21],
        [18, 17, 17, 20, 20, 25, 25, 31, 31, 01, 01, 04, 04, 09, 09, 14],
        [23, 22, 22, 09, 09, 18, 18, 12, 12, 30, 30, 01, 01, 26, 26, 07],
        [20, 23, 23, 06, 06, 11, 11, 25, 25, 19, 19, 02, 02, 15, 15, 28],
        [21, 20, 20, 31, 31, 04, 04, 14, 14, 08, 08, 03, 03, 24, 24, 17],
        [22, 21, 21, 16, 16, 29, 29, 27, 27, 05, 05, 00, 00, 13, 13, 10],
        [27, 26, 26, 05, 05, 30, 30, 00, 00, 18, 18, 13, 13, 22, 22, 11],
        [24, 27, 27, 10, 10, 07, 07, 21, 21, 31, 31, 14, 14, 03, 03, 16],
        [25, 24, 24, 19, 19, 08, 08, 02, 02, 04, 04, 15, 15, 20, 20, 29],
        [26, 25, 25, 28, 28, 17, 17, 23, 23, 09, 09, 12, 12, 01, 01, 06],
        [31, 30, 30, 01, 01, 26, 26, 04, 04, 22, 22, 09, 09, 18, 18, 15],
        [28, 31, 31, 14, 14, 03, 03, 17, 17, 27, 27, 10, 10, 07, 07, 20],
        [29, 28, 28, 23, 23, 12, 12, 06, 06, 00, 00, 11, 11, 16, 16, 25],
        [30, 29, 29, 24, 24, 21, 21, 19, 19, 13, 13, 08, 08, 05, 05, 02],
    ];
    
    const MORTON_LUT_D4: [[u8; 16]; 32] = [
        [00, 08, 09, 01, 03, 11, 10, 02, 06, 14, 15, 07, 05, 13, 12, 04],
        [00, 04, 12, 08, 09, 13, 05, 01, 03, 07, 15, 11, 10, 14, 06, 02],
        [00, 02, 06, 04, 12, 14, 10, 08, 09, 11, 15, 13, 05, 07, 03, 01],
        [00, 01, 03, 02, 06, 07, 05, 04, 12, 13, 15, 14, 10, 11, 09, 08],
        [03, 11, 10, 02, 00, 08, 09, 01, 05, 13, 12, 04, 06, 14, 15, 07],
        [03, 07, 15, 11, 10, 14, 06, 02, 00, 04, 12, 08, 09, 13, 05, 01],
        [03, 01, 05, 07, 15, 13, 09, 11, 10, 08, 12, 14, 06, 04, 00, 02],
        [03, 02, 00, 01, 05, 04, 06, 07, 15, 14, 12, 13, 09, 08, 10, 11],
        [06, 14, 15, 07, 05, 13, 12, 04, 00, 08, 09, 01, 03, 11, 10, 02],
        [06, 02, 10, 14, 15, 11, 03, 07, 05, 01, 09, 13, 12, 08, 00, 04],
        [06, 04, 00, 02, 10, 08, 12, 14, 15, 13, 09, 11, 03, 01, 05, 07],
        [06, 07, 05, 04, 00, 01, 03, 02, 10, 11, 09, 08, 12, 13, 15, 14],
        [05, 13, 12, 04, 06, 14, 15, 07, 03, 11, 10, 02, 00, 08, 09, 01],
        [05, 01, 09, 13, 12, 08, 00, 04, 06, 02, 10, 14, 15, 11, 03, 07],
        [05, 07, 03, 01, 09, 11, 15, 13, 12, 14, 10, 08, 00, 02, 06, 04],
        [05, 04, 06, 07, 03, 02, 00, 01, 09, 08, 10, 11, 15, 14, 12, 13],
        [12, 04, 05, 13, 15, 07, 06, 14, 10, 02, 03, 11, 09, 01, 00, 08],
        [12, 08, 00, 04, 05, 01, 09, 13, 15, 11, 03, 07, 06, 02, 10, 14],
        [12, 14, 10, 08, 00, 02, 06, 04, 05, 07, 03, 01, 09, 11, 15, 13],
        [12, 13, 15, 14, 10, 11, 09, 08, 00, 01, 03, 02, 06, 07, 05, 04],
        [15, 07, 06, 14, 12, 04, 05, 13, 09, 01, 00, 08, 10, 02, 03, 11],
        [15, 11, 03, 07, 06, 02, 10, 14, 12, 08, 00, 04, 05, 01, 09, 13],
        [15, 13, 09, 11, 03, 01, 05, 07, 06, 04, 00, 02, 10, 08, 12, 14],
        [15, 14, 12, 13, 09, 08, 10, 11, 03, 02, 00, 01, 05, 04, 06, 07],
        [10, 02, 03, 11, 09, 01, 00, 08, 12, 04, 05, 13, 15, 07, 06, 14],
        [10, 14, 06, 02, 03, 07, 15, 11, 09, 13, 05, 01, 00, 04, 12, 08],
        [10, 08, 12, 14, 06, 04, 00, 02, 03, 01, 05, 07, 15, 13, 09, 11],
        [10, 11, 09, 08, 12, 13, 15, 14, 06, 07, 05, 04, 00, 01, 03, 02],
        [09, 01, 00, 08, 10, 02, 03, 11, 15, 07, 06, 14, 12, 04, 05, 13],
        [09, 13, 05, 01, 00, 04, 12, 08, 10, 14, 06, 02, 03, 07, 15, 11],
        [09, 11, 15, 13, 05, 07, 03, 01, 00, 02, 06, 04, 12, 14, 10, 08],
        [09, 08, 10, 11, 15, 14, 12, 13, 05, 04, 06, 07, 03, 02, 00, 01],
    ];
    
    const HILBERT_LUT_D4: [[u8; 16]; 32] = [
        [00, 03, 07, 04, 15, 12, 08, 11, 01, 02, 06, 05, 14, 13, 09, 10],
        [00, 07, 15, 08, 01, 06, 14, 09, 03, 04, 12, 11, 02, 05, 13, 10],
        [00, 15, 01, 14, 03, 12, 02, 13, 07, 08, 06, 09, 04, 11, 05, 10],
        [00, 01, 03, 02, 07, 06, 04, 05, 15, 14, 12, 13, 08, 09, 11, 10],
        [04, 07, 03, 00, 11, 08, 12, 15, 05, 06, 02, 01, 10, 09, 13, 14],
        [08, 15, 07, 00, 09, 14, 06, 01, 11, 12, 04, 03, 10, 13, 05, 02],
        [14, 01, 15, 00, 13, 02, 12, 03, 09, 06, 08, 07, 10, 05, 11, 04],
        [02, 03, 01, 00, 05, 04, 06, 07, 13, 12, 14, 15, 10, 11, 09, 08],
        [08, 11, 15, 12, 07, 04, 00, 03, 09, 10, 14, 13, 06, 05, 01, 02],
        [14, 09, 01, 06, 15, 08, 00, 07, 13, 10, 02, 05, 12, 11, 03, 04],
        [02, 13, 03, 12, 01, 14, 00, 15, 05, 10, 04, 11, 06, 09, 07, 08],
        [04, 05, 07, 06, 03, 02, 00, 01, 11, 10, 08, 09, 12, 13, 15, 14],
        [12, 15, 11, 08, 03, 00, 04, 07, 13, 14, 10, 09, 02, 01, 05, 06],
        [06, 01, 09, 14, 07, 00, 08, 15, 05, 02, 10, 13, 04, 03, 11, 12],
        [12, 03, 13, 02, 15, 00, 14, 01, 11, 04, 10, 05, 08, 07, 09, 06],
        [06, 07, 05, 04, 01, 00, 02, 03, 09, 08, 10, 11, 14, 15, 13, 12],
        [14, 13, 09, 10, 01, 02, 06, 05, 15, 12, 08, 11, 00, 03, 07, 04],
        [02, 05, 13, 10, 03, 04, 12, 11, 01, 06, 14, 09, 00, 07, 15, 08],
        [04, 11, 05, 10, 07, 08, 06, 09, 03, 12, 02, 13, 00, 15, 01, 14],
        [08, 09, 11, 10, 15, 14, 12, 13, 07, 06, 04, 05, 00, 01, 03, 02],
        [10, 09, 13, 14, 05, 06, 02, 01, 11, 08, 12, 15, 04, 07, 03, 00],
        [10, 13, 05, 02, 11, 12, 04, 03, 09, 14, 06, 01, 08, 15, 07, 00],
        [10, 05, 11, 04, 09, 06, 08, 07, 13, 02, 12, 03, 14, 01, 15, 00],
        [10, 11, 09, 08, 13, 12, 14, 15, 05, 04, 06, 07, 02, 03, 01, 00],
        [06, 05, 01, 02, 09, 10, 14, 13, 07, 04, 00, 03, 08, 11, 15, 12],
        [12, 11, 03, 04, 13, 10, 02, 05, 15, 08, 00, 07, 14, 09, 01, 06],
        [06, 09, 07, 08, 05, 10, 04, 11, 01, 14, 00, 15, 02, 13, 03, 12],
        [12, 13, 15, 14, 11, 10, 08, 09, 03, 02, 00, 01, 04, 05, 07, 06],
        [02, 01, 05, 06, 13, 14, 10, 09, 03, 00, 04, 07, 12, 15, 11, 08],
        [04, 03, 11, 12, 05, 02, 10, 13, 07, 00, 08, 15, 06, 01, 09, 14],
        [08, 07, 09, 06, 11, 04, 10, 05, 15, 00, 14, 01, 12, 03, 13, 02],
        [14, 15, 13, 12, 09, 08, 10, 11, 01, 00, 02, 03, 06, 07, 05, 04],
    ];

    #[inline(always)]
    pub fn transform_lut<T: NumTraits, const D: usize>(
        transform: usize,
        hilbert_partial: T,
    ) -> usize {
        (match D {
            2 => TRANSFORM_LUT_D2[transform][hilbert_partial.to_usize()],
            3 => TRANSFORM_LUT_D3[transform][hilbert_partial.to_usize()],
            4 => TRANSFORM_LUT_D4[transform][hilbert_partial.to_usize()],
            _ => unimplemented!(),
        }) as usize
    }

    #[inline(always)]
    pub fn morton_lut<T: NumTraits, const D: usize>(
        transform: usize,
        hilbert_partial: T,
    ) -> T {
        T::from_u8(match D {
            2 => MORTON_LUT_D2[transform][hilbert_partial.to_usize()],
            3 => MORTON_LUT_D3[transform][hilbert_partial.to_usize()],
            4 => MORTON_LUT_D4[transform][hilbert_partial.to_usize()],
            _ => unimplemented!(),
        })
    }

    #[inline(always)]
    pub fn hilbert_lut<T: NumTraits, const D: usize>(
        transform: usize,
        morton_partial: T,
    ) -> T {
        T::from_u8(match D {
            2 => HILBERT_LUT_D2[transform][morton_partial.to_usize()],
            3 => HILBERT_LUT_D3[transform][morton_partial.to_usize()],
            4 => HILBERT_LUT_D4[transform][morton_partial.to_usize()],
            _ => unimplemented!(),
        })
    }

    #[inline(always)]
    pub fn walk_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        min_order: usize,
        mut f: F,
    ) -> usize
    where
        F: FnMut(usize, T, usize) -> T,
    {
        let lower_mask = T::one().shl(D).sub(T::one());
        let order = ((T::bits() + D - 1) - index.lz()) / D;
        let mut transform = order % D;
        for shift in (min_order..order).rev().map(|i| i * D) {
            let hilbert_partial = f(transform, index.shr(shift).bit_and(lower_mask), shift);
            transform = transform_lut::<_, D>(transform, hilbert_partial);
        }
        transform
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        debug_assert!(D <= MAX_LUT_D);
        let mut hilbert_index = T::zero();
        walk_transforms::<_, _, D>(
            morton_index,
            0,
            |transform, morton_partial, shift| {
                let hilbert_partial = hilbert_lut::<_, D>(transform, morton_partial);
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        debug_assert!(D <= MAX_LUT_D);
        let mut morton_index = T::zero();
        walk_transforms::<_, _, D>(
            hilbert_index,
            0,
            |transform, hilbert_partial, shift| {
                let morton_partial = morton_lut::<_, D>(transform, hilbert_partial);
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
    }
}

mod explicit {
    use super::NumTraits;
    
    #[inline(always)]
    pub fn gray<T: NumTraits>(index: T) -> T {
        index.bit_xor(index.shr(1))
    }
    
    #[inline(always)]
    pub fn gray_inverse<T: NumTraits, const D: usize>(index: T) -> T {
        // Since we know that only a small subset of bits of the integer are used,
        // we can perform a minimal gray code inverse rather than a complete one
        (1..D).fold(index, |acc, shift| acc.bit_xor(index.shr(shift)))
    }
    
    #[inline(always)]
    pub fn shl_cyclic<T: NumTraits, const D: usize>(index: T, amount: usize) -> T {
        let lower_mask: T = <T as NumTraits>::one().shl(D).sub(NumTraits::one());
        lower_mask.bit_and(index.shl(amount).bit_or(index.shr(D - amount)))
    }
    
    #[inline(always)]
    pub fn shr_cyclic<T: NumTraits, const D: usize>(index: T, amount: usize) -> T {
        let lower_mask: T = <T as NumTraits>::one().shl(D).sub(NumTraits::one());
        lower_mask.bit_and(index.shr(amount).bit_or(index.shl(D - amount)))
    }
    
    #[inline(always)]
    pub fn rotation_transform(index: usize) -> usize {
        if index > 0 {
            index.wrapping_add((index & 0x1).wrapping_sub(1)).trailing_ones() as usize
        } else {
            0
        }
    }
    
    #[inline(always)]
    pub fn axis_transform<T: NumTraits>(index: T) -> T {
        if index == T::zero() {
            T::zero()
        } else {
            gray(index.sub(T::one()).bit_and(T::bit_not(T::one())))
        }
    }

    #[inline(always)]
    pub fn walk_transforms<T: NumTraits, F, const D: usize>(
        index: T,
        min_order: usize,
        mut f: F,
    ) -> (usize, T)
    where
        F: FnMut(T, usize, usize, T) -> T,
    {
        let lower_mask = T::one().shl(D).sub(T::one());
        // A subset of highest bits of index will be 0, therefore, we can skip a certain
        // number of iterations in the following loop. Note that skipping iterations
        // still requires us to rotate the pattern as if we had iterated all orders, see
        // rotate_amount below.
        let order = (T::bits() - index.lz() + D - 1) / D;

        // We iterate in reverse because higher orders affect the rotation of lower orders
        let mut flip_axes = T::zero();
        let mut rotate_amount = (D - 1) - (order % D);
        for shift in (min_order..order).rev().map(|i| i * D) {
            // Extract a portion of the index
            let index_partial = index.shr(shift).bit_and(lower_mask);

            // Context specific code generates the hilbert partial
            let hilbert_partial = f(index_partial, shift, rotate_amount, flip_axes);

            // Calculate the next axis flip flags and rotation amounts resulting from this order
            flip_axes = flip_axes.bit_xor(shl_cyclic::<_, D>(
                axis_transform(hilbert_partial),
                rotate_amount,
            ));
            rotate_amount =
                (rotation_transform(hilbert_partial.to_usize()) + rotate_amount + 1) % D;
        }
        (rotate_amount, flip_axes)
    }

    #[inline(always)]
    pub fn morton_to_hilbert<T: NumTraits, const D: usize>(morton_index: T) -> T {
        let mut hilbert_index = T::zero();
        walk_transforms::<_, _, D>(
            morton_index,
            0,
            |morton_partial, shift, rotate_amount, flip_axes| {
                let hilbert_partial = gray_inverse::<_, D>(shr_cyclic::<_, D>(
                    flip_axes.bit_xor(morton_partial),
                    rotate_amount,
                ));
                hilbert_index = hilbert_index.bit_or(hilbert_partial.shl(shift));
                hilbert_partial
            },
        );
        hilbert_index
    }

    #[inline(always)]
    pub fn hilbert_to_morton<T: NumTraits, const D: usize>(hilbert_index: T) -> T {
        let mut morton_index = T::zero();
        walk_transforms::<_, _, D>(
            hilbert_index,
            0,
            |hilbert_partial, shift, rotate_amount, flip_axes| {
                let morton_partial = flip_axes.bit_xor(shl_cyclic::<_, D>(
                    gray(hilbert_partial),
                    rotate_amount,
                ));
                morton_index = morton_index.bit_or(morton_partial.shl(shift));
                hilbert_partial
            },
        );
        morton_index
    }
}

// Until we have complex generic constants, we have to pass D in here (needed by coords)
// Waiting on: https://github.com/rust-lang/rust/issues/76560
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Hilbert<DM, const D: usize>(pub(crate) DM::Dilated)
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex;

impl<DM, const D: usize> SpaceFillingCurve<D> for Hilbert<DM, D>
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

    #[inline(always)]
    fn from_index(index: Self::Index) -> Self {
        debug_assert!(
            index <= Self::INDEX_MAX,
            "Parameter 'index' exceeds maximum"
        );
        Self(index)
    }

    #[inline(never)] // Temporary
    fn from_coords(coords: [Self::Coord; D]) -> Self {
        let morton_index = Morton::<DM, D>::from_coords(coords).index();

        #[cfg(feature = "lut_large")]
        if D <= lut_large::MAX_LUT_D {
            return Self(lut_large::morton_to_hilbert::<_, D>(morton_index));
        }

        #[cfg(feature = "lut_small")]
        if D <= lut_small::MAX_LUT_D {
            return  Self(lut_small::morton_to_hilbert::<_, D>(morton_index));
        }

        Self(explicit::morton_to_hilbert::<_, D>(morton_index))
    }

    #[inline(never)] // Temporary
    fn coords(&self) -> [Self::Coord; D] {
        #[cfg(feature = "lut_large")]
        if D <= lut_large::MAX_LUT_D {
            return Morton::<DM, D>::from_index(lut_large::hilbert_to_morton::<_, D>(self.0)).coords();
        }

        #[cfg(feature = "lut_small")]
        if D <= lut_small::MAX_LUT_D {
            return  Morton::<DM, D>::from_index(lut_small::hilbert_to_morton::<_, D>(self.0)).coords();
        }

        Morton::<DM, D>::from_index(explicit::hilbert_to_morton::<_, D>(self.0)).coords()
    }

    #[inline(always)]
    fn index(&self) -> Self::Index {
        self.0
    }
}

impl<DM, const D: usize> Siblings<D> for Hilbert<DM, D>
where
    // When https://github.com/rust-lang/rust/issues/52662 is available, we can clean this up
    DM: DilationMethod,
    DM::Undilated: CurveCoord,
    DM::Dilated: CurveIndex,
{
    fn sibling_on_axis_toggle(&self, axis: usize) -> Self {
        todo!();
    }

    fn sibling_on_axis(&self, axis: usize, direction: QueryDirection) -> Self {
        todo!()
    }

    fn sibling_from_bits(&self, axis_bits: Self::Index) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    extern crate std;

    use dilate::DilationMethod;
    use super::*;

    use crate::internal::tests::{test_curve, test_curve_siblings};

    macro_rules! hilbert_expand {
        ($t:ty, $d:literal) => {
            Hilbert::<Expand<$t, $d>, $d>
        };
    }

    macro_rules! hilbert_fixed {
        ($t:ty, $d:literal) => {
            Hilbert::<Fixed<$t, $d>, $d>
        };
    }

    #[test]
    fn test_gray() {
//        for i in 0..10000usize {
//            let g = super::gray(i);
//            assert_eq!(super::gray_inverse::<_, { usize::BITS as usize }>(g), i);
//        }
    }

    #[test]
    #[ignore]
    fn generate_lut() {
        generate_lut_small();
        //generate_lut_for_dimension::<dilate::Expand<u8, 2>, 2>();
        //generate_lut_for_dimension::<dilate::Expand<u8, 3>, 3>();
        //generate_lut_for_dimension::<dilate::Expand<u8, 4>, 4>();
    }

    fn generate_lut_small() {
        let max_d: usize = 4;
        let num_elems: usize = 1 << max_d;

        println!("    pub const MAX_LUT_D: usize = {max_d};\n");
        
        println!("    const GRAY_LUT: [u8; {num_elems}] = [");
        print!("        ");
        (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::gray(i)));
        println!("\n    ];\n");

        println!("    const GRAY_INV_LUT: [u8; {num_elems}] = [");
        print!("        ");
        (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::gray_inverse::<_, 4>(i)));
        println!("\n    ];\n");

        generate_shlr_luts::<2>();

        generate_shlr_luts::<3>();

        generate_shlr_luts::<4>();

        println!("    const ROTATION_TRANSFORM_LUT: [u8; {num_elems}] = [");
        print!("        ");
        (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::rotation_transform(i)));
        println!("\n    ];\n");

        println!("    const AXIS_TRANSFORM_LUT: [u8; {num_elems}] = [");
        print!("        ");
        (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::axis_transform(i)));
        println!("\n    ];");
        println!("}}\n");
    }

    fn generate_shlr_luts<const D: usize>() {
        let num_elems: usize = 1 << D;

        println!("    const AXIS_TX_INTO_SHL_LUT_D{D}: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::shl_cyclic::<_, D>(explicit::axis_transform(i), r)));
            print!("],\n");
        }
        println!("    ];\n");

        println!("    const SHR_INTO_GRAY_INV_LUT_D{D}: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::gray_inverse::<_, D>(explicit::shr_cyclic::<_, D>(i, r))));
            print!("],\n");
        }
        println!("    ];\n");

        println!("    const GRAY_INTO_SHL_LUT_D{D}: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", explicit::shl_cyclic::<_, D>(explicit::gray(i), r)));
            print!("],\n");
        }
        println!("    ];\n");

        println!("    const ROTATION_INCREMENT_LUT_D{D}: [[u8; {num_elems}]; {D}] = [");
        for r in 0..D {
            print!("        [");
            (0..num_elems).into_iter().for_each(|i| print!("{:0>2}, ", (explicit::rotation_transform(i) + r + 1) % D));
            print!("],\n");
        }
        println!("    ];\n");
    }

    fn generate_lut_for_dimension<DM, const D: usize>()
    where
        DM: DilationMethod,
        DM::Undilated: CurveCoord,
        DM::Dilated: CurveIndex + std::fmt::Display,
    {
        match D {
            2 => _generate_lut_for_dimension::<DM, D, { 1 << 2 }, { 2 * (1 << (2 - 1)) }>(),
            3 => _generate_lut_for_dimension::<DM, D, { 1 << 3 }, { 3 * (1 << (3 - 1)) }>(),
            4 => _generate_lut_for_dimension::<DM, D, { 1 << 4 }, { 4 * (1 << (4 - 1)) }>(),
            _ => unimplemented!(),
        }
    }

    fn _generate_lut_for_dimension<DM, const D: usize, const NUM_CHILDREN: usize, const NUM_TRANSFORMS: usize>()
    where
        DM: DilationMethod,
        DM::Undilated: CurveCoord,
        DM::Dilated: CurveIndex + std::fmt::Display,
    {
        let mut num_transforms = 0;
        let mut transform_map = [[None; D]; NUM_CHILDREN];
        let mut transforms = [[(0, DM::Dilated::zero(), DM::Dilated::zero()); NUM_CHILDREN]; NUM_TRANSFORMS];

        let mut num_flip_axes = 0;
        let mut flip_axes_map = [None; NUM_CHILDREN];

        let extra_debug = false;

        let num_rotations = D;
        let max_rotation = num_rotations - 1;
        let num_axis_flips = 1 << (D - 1);

        if extra_debug {
            println!("");
            println!("D: {}", D);
            println!("num_rotations: {}", num_rotations);
            println!("num_axis_flips: {}", num_axis_flips);
        }

        // Enumerate all possible axis flips
        // This forces the ordering of flip IDs to be known values
        for i in 0..NUM_CHILDREN {
            let flip_axes = explicit::axis_transform(DM::Dilated::from_usize(i));
            let flip_id = flip_axes_map[flip_axes.to_usize()].unwrap_or_else(|| {
                if extra_debug {
                    println!("flip_axes {num_flip_axes:?}: hilbert_partial = {i:?}, flip_axes = {flip_axes:?}");
                }
                num_flip_axes += 1;
                num_flip_axes - 1
            });
            flip_axes_map[flip_axes.to_usize()] = Some(flip_id);
        }

        // Enumerate all possible transforms
        // This forces the order of transform indices to be known values
        for i in 0..NUM_CHILDREN {
            for rotate_amount in 0..D {
                // Find initial flip and rotate for index partial i (note xor of previous flip_axes is not required here, but previous rotation is)
                let flip_axes = explicit::shl_cyclic::<_, D>(
                    explicit::axis_transform(DM::Dilated::from_usize(i)),
                    rotate_amount,
                );
                let rotate_amount = (explicit::rotation_transform(i) + rotate_amount + 1) % D;

                let flip_id = flip_axes_map[flip_axes.to_usize()].unwrap();

                // Flip rotation here avoids flip at runtime
                let flipped_rotation = max_rotation - rotate_amount;
                let tx_id = flip_id * num_rotations + flipped_rotation;

                if extra_debug {
                    println!("Transform {tx_id:?}: rotate_amount = {rotate_amount:?}, flip_axes = {flip_axes:?}");
                }

                // Store ID of this transform
                if transform_map[flip_axes.to_usize()][rotate_amount].is_none() {
                    num_transforms += 1;
                }
                transform_map[flip_axes.to_usize()][rotate_amount] = Some(tx_id);
            }
        }
        assert_eq!(num_transforms, NUM_TRANSFORMS);

        // Iterate all possible transforms
        for rotate_amount in 0..D {
            for i in 0..NUM_CHILDREN {
                // Find initial flip and rotate for index partial i (note xor of previous flip_axes is not required here, but previous rotation is)
                let flip_axes = explicit::shl_cyclic::<_, D>(
                    explicit::axis_transform(DM::Dilated::from_usize(i)),
                    rotate_amount,
                );
                let rotate_amount = (explicit::rotation_transform(i) + rotate_amount + 1) % D;

                // ID of this transform
                let tx_id = transform_map[flip_axes.to_usize()][rotate_amount].unwrap();

                // Within this transform, it's possible to decend to a number of other transforms depending on the index partial j
                for j in 0..NUM_CHILDREN {
                    // Store resultant morton partial and hilbert partial (to coords and from coords)
                    transforms[tx_id][j].1 = flip_axes.bit_xor(explicit::shl_cyclic::<_, D>(
                        explicit::gray(DM::Dilated::from_usize(j)),
                        rotate_amount,
                    ));
                    transforms[tx_id][j].2 = explicit::gray_inverse::<_, D>(explicit::shr_cyclic::<_, D>(
                        flip_axes.bit_xor(DM::Dilated::from_usize(j)),
                        rotate_amount,
                    ));

                    // To find what the next transform would be, we flip and rotate again
                    let flip_axes = flip_axes.bit_xor(explicit::shl_cyclic::<_, D>(
                        explicit::axis_transform(DM::Dilated::from_usize(j)),
                        rotate_amount,
                    ));
                    let rotate_amount = (explicit::rotation_transform(j) + rotate_amount + 1) % D;

                    // Cache next transform ID
                    transforms[tx_id][j].0 =
                        transform_map[flip_axes.to_usize()][rotate_amount].unwrap();
                }
            }
        }

        if extra_debug {
            println!("num_transforms: {}\n", num_transforms);
        }

        let md = false;
        if md {
            print!("|    hibert_index -> | ");
            for i in 0..NUM_CHILDREN {
                print!("{:2} | ", i);
            }
            println!("");
            print!("| ------------------ | ");
            for _ in 0..NUM_CHILDREN {
                print!("-- | ");
            }
            println!("");
            for i in 0..num_transforms {
                print!("| transform_index {:2} | ", i);
                for j in 0..NUM_CHILDREN {
                    let lut_data = transforms[i][j];
                    print!("{:2} | ", lut_data.0 / num_axis_flips);
                }
                println!("");
            }
        } else {
            println!(
                "pub const TRANSFORM_LUT_D{}: [[u8; {}]; {}] = [",
                D, NUM_CHILDREN, NUM_TRANSFORMS
            );
            for i in 0..num_transforms {
                print!("    [");
                for j in 0..(NUM_CHILDREN - 1) {
                    let lut_data = transforms[i][j];
                    print!("{:0>2}, ", lut_data.0);
                }
                let lut_data = transforms[i][NUM_CHILDREN - 1];
                println!("{:0>2}],", lut_data.0);
            }
            println!("];\n");
            println!(
                "pub const MORTON_LUT_D{}: [[u8; {}]; {}] = [",
                D, NUM_CHILDREN, NUM_TRANSFORMS
            );
            for i in 0..num_transforms {
                print!("    [");
                for j in 0..(NUM_CHILDREN - 1) {
                    let lut_data = transforms[i][j];
                    print!("{:0>2}, ", lut_data.1);
                }
                let lut_data = transforms[i][NUM_CHILDREN - 1];
                println!("{:0>2}],", lut_data.1);
            }
            println!("];\n");
            println!(
                "pub const HILBERT_LUT_D{}: [[u8; {}]; {}] = [",
                D, NUM_CHILDREN, NUM_TRANSFORMS
            );
            for i in 0..num_transforms {
                print!("    [");
                for j in 0..(NUM_CHILDREN - 1) {
                    let lut_data = transforms[i][j];
                    print!("{:0>2}, ", lut_data.2);
                }
                let lut_data = transforms[i][NUM_CHILDREN - 1];
                println!("{:0>2}],", lut_data.2);
            }
            println!("];\n");
        }
    }

    test_curve!(hilbert_expand, true, u8, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_expand, true, u16, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_expand, true, u32, 2, 3, 4);
    test_curve!(hilbert_expand, true, u64, 2);

    test_curve!(hilbert_fixed, true, u8, 2, 3, 4);
    test_curve!(hilbert_fixed, true, u16, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_fixed, true, u32, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_fixed, true, u64, 2, 3, 4, 5, 6, 7, 8);
    test_curve!(hilbert_fixed, true, u128, 2, 3, 4, 5, 6, 7, 8);

    //    test_curve_siblings!(hilbert_expand, true, u8, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_expand, true, u16, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_expand, true, u32, 2, 3, 4);
    //    test_curve_siblings!(hilbert_expand, true, u64, 2);
    //
    //    test_curve_siblings!(hilbert_fixed, true, u8, 2, 3, 4);
    //    test_curve_siblings!(hilbert_fixed, true, u16, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_fixed, true, u32, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_fixed, true, u64, 2, 3, 4, 5, 6, 7, 8);
    //    test_curve_siblings!(hilbert_fixed, true, u128, 2, 3, 4, 5, 6, 7, 8);
}
