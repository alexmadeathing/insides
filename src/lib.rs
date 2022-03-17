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

#![no_std]
#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]
#![deny(rustdoc::invalid_rust_codeblocks)]

//! A compact, high performance space filling curve library for Rust.
//! 
//! TODO Write intro
//! 
//! # Examples
//! ```
//! use insides::*;
//! 
//! let xyz = [1, 2, 3];
//! 
//! // Create a 3D morton location using u16 coordinate indices
//! // At the moment, we have to specify the number of dimensions twice, sorry!
//! // (this will change with improvements to Rust const generics)
//! let location = Morton::<Expand<u16, 3>, 3>::from_coords(xyz);
//! 
//! // Access raw morton location index
//! // Useful as a map key or array index!
//! // In this case, we're using a morton curve, so we know the underlying bit
//! // layout. If we used a different type of curve, the index would be
//! // different.
//! assert_eq!(location.index(), 0b110101);
//! 
//! // Get neighbour in the negative X axis
//! let neighbour_a = location.neighbour_on_axis(0, QueryDirection::Negative);
//! assert_eq!(neighbour_a.coords(), [0, 2, 3]);
//! 
//! // Get neighbour in the positive Y axis
//! let neighbour_b = location.neighbour_on_axis(1, QueryDirection::Positive);
//! assert_eq!(neighbour_b.coords(), [1, 3, 3]);
//! 
//! // Get sibling on Z axis
//! // Its use is subtle and may not be apparent at first glance. Imagine the
//! // world is divided into a grid where each cell contains 8 siblings (in the
//! // 3D case). This method operates upon that grid. This is mostly useful
//! // when applied to tree-like structures.
//! let sibling_a = location.sibling_on_axis(2);
//! assert_eq!(sibling_a.coords(), [1, 2, 2]);
//! 
//! // Calling again from the sibling gets our original location
//! let sibling_b = sibling_a.sibling_on_axis(2);
//! assert_eq!(sibling_b, location);
//! 
//! // This gets, either the sibling, or the same location, on the Z axis,
//! // depending on which is further in the query direction.
//! // This is useful when you want to find the locations on one edge of a cell.
//! let sibling_or_same = location.sibling_or_same_on_axis(2, QueryDirection::Positive);
//! assert_eq!(sibling_or_same.coords(), [1, 2, 3]);
//! ```

pub mod morton_curve;

pub use dilate::{Expand, Fixed};
pub use morton_curve::{Morton, MortonIndex};

#[derive(Clone, Copy, Debug)]
pub enum QueryDirection {
    Positive,
    Negative,
}

// I'd really love to get rid of the generic parameter here but I think it's waiting on:
// https://github.com/rust-lang/rust/issues/76560
pub trait Encoding<const D: usize> {
    type Coord;
    const COORD_MAX: Self::Coord;

    fn from_coords(coords: [Self::Coord; D]) -> Self;
    fn coords(&self) -> [Self::Coord; D];
}

pub trait Siblings {
    fn sibling_on_axis(&self, axis: usize) -> Self;
    fn sibling_or_same_on_axis(&self, axis: usize, direction: QueryDirection) -> Self;
}

pub trait Neighbours {
    fn neighbour_on_axis(&self, axis: usize, direction: QueryDirection) -> Self;
}
