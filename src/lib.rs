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

//#![no_std]
#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]
#![deny(rustdoc::invalid_rust_codeblocks)]

//! A compact, high performance space filling curve library for Rust.
//! 
//! This library provides an abstract interface to generalise interactions with
//! space filling curves as well as a morton curve implementation and
//! supporting manipulation methods.
//! 
//! We currently only support Morton Encoding, but the interface will support
//! a Hilbert implementation - and we have plans to include it soon.
//!
//! # Morton Encoding 
//! Morton encoding, also known as a
//! [Z-order curve](https://en.wikipedia.org/wiki/Z-order_curve), is a space
//! filling algorithm which maps a multidimensional set of coordinates to one
//! dimension, achieved by interleaving the bit sequence of each coordinate
//! value.
//! 
//! Whilst other encoding methods may exhibit better spatial locality (such as
//! the Hilbert curve), the Morton curve offers excellent CPU performance,
//! since most behaviours can be reduced to a simple set of bitwise operations -
//! making it an ideal choice for applications such and quad trees and octrees.
//! 
//! # Examples
//! ```rust
//! use insides::*;
//! 
//! // Create a 3D morton location using u16 coordinate indices
//! // At the moment, we have to specify the number of dimensions twice, sorry!
//! // (this will change with improvements to Rust const generics)
//! let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
//! 
//! assert_eq!(location.index(), 0b110101);
//! assert_eq!(location.coords(), [1, 2, 3]);
//! ```

#[doc(hidden)]
pub use dilate::*;

/// Morton curve implementation
pub mod morton;
pub use crate::morton::Morton;

/// Hilbert curve implementation
pub mod hilbert;
pub use crate::hilbert::Hilbert;

pub(crate) mod internal;

/// Trait wrapper for coordinates
pub trait Coord: internal::NumTraits {}

impl Coord for u8 {}
impl Coord for u16 {}
impl Coord for u32 {}
impl Coord for u64 {}
impl Coord for u128 {}
impl Coord for usize {}

/// Trait wrapper for index used by [Hilbert]
pub trait Index: internal::NumTraits {}

impl Index for u8 {}
impl Index for u16 {}
impl Index for u32 {}
impl Index for u64 {}
impl Index for u128 {}
impl Index for usize {}

/// Direction to search within an axis when making queries
#[derive(Clone, Copy, Debug)]
pub enum QueryDirection {
    /// Search in a positive direction along an axis relative to the source location
    Positive,

    /// Search in a negative direction along an axis relative to the source location
    Negative,
}

/// Provides conversion to and from coordinates
// I'd really love to get rid of the generic parameter here but I think it's waiting on:
// https://github.com/rust-lang/rust/issues/76560
pub trait SpaceFillingCurve<const D: usize> {
    /// Coordinate type
    type Coord: Coord;

    /// Index type
    type Index: Index;

    /// Maximum value for a single coordinate
    const COORD_MAX: Self::Coord;

    /// Maximum value for an encoded index
    const INDEX_MAX: Self::Index;

    /// Number of dimensions
    const D: usize;

    /// Create a new location from a pre-encoded curve index
    /// 
    /// This function is provided for situations where an existing
    /// encoded index is available and you wish to use the various
    /// manipulation methods provided by this library.
    /// 
    /// # Panics
    /// Panics if parameter 'index' is greater than Self::INDEX_MAX
    /// 
    /// # Examples
    /// ```rust
    /// use insides::*;
    /// 
    /// let location = Morton::<Expand<u16, 3>, 3>::from_index(0b110101);
    /// 
    /// assert_eq!(location.index(), 0b110101);
    /// assert_eq!(location.coords(), [1, 2, 3]);
    /// ```
    fn from_index(index: Self::Index) -> Self;

    /// Encode curve location from coordinates
    /// 
    /// This function converts raw coordinates to an encoded curve
    /// representation. The resultant object contains the raw index and may
    /// provide various additional manipulation methods via [Siblings] and
    /// [Neighbours].
    /// 
    /// # Panics
    /// Panics if coordinate is greater than [SpaceFillingCurve::COORD_MAX].
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
    fn from_coords(coords: [Self::Coord; D]) -> Self;

    /// Decode curve index into coordinates
    /// 
    /// This method converts an encoded curve location back into raw
    /// coordinates.
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
    fn coords(&self) -> [Self::Coord; D];

    /// Access curve location index
    /// 
    /// This method retrieves the encoded curve location index.
    /// 
    /// # Examples
    /// ```rust
    /// use insides::*;
    /// 
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    /// 
    /// assert_eq!(location.index(), 0b110101);
    /// ```
    fn index(&self) -> Self::Index;
}

/// Provides methods for retrieving adacent locations within a cluster
pub trait Siblings<const D: usize>: SpaceFillingCurve<D> {
    /// Get sibling location on axis within local cluster
    /// 
    /// The sibling is considered to be within the local cluster this location
    /// resides within, where a cluster is the collection of 2<sup>D</sup>
    /// adjacent locations and the cluster origin is a vector of even numbers.
    /// 
    /// # Panics
    /// Panics if parameter 'axis' is greater than or equal to D
    /// 
    /// # Examples
    /// ```rust
    /// use insides::*;
    /// 
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    /// 
    /// assert_eq!(location.sibling_on_axis(0).coords(), [0, 2, 3]);
    /// assert_eq!(location.sibling_on_axis(1).coords(), [1, 3, 3]);
    /// assert_eq!(location.sibling_on_axis(2).coords(), [1, 2, 2]);
    /// ```
    fn sibling_on_axis(&self, axis: usize) -> Self;

    /// Get sibling or same location on axis within local cluster
    /// 
    /// This method gets either the sibling, or the same location, depending on
    /// which is further in the query direction along a given axis.
    /// 
    /// The sibling is considered to be within the local cluster this location
    /// resides within, where a cluster is the collection of 2<sup>D</sup>
    /// adjacent locations and the cluster origin is a vector of even numbers.
    /// 
    /// # Panics
    /// Panics if parameter 'axis' is greater than or equal to D
    /// 
    /// # Examples
    /// ```rust
    /// use insides::*;
    /// 
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    /// 
    /// assert_eq!(location.sibling_or_same_on_axis(0, QueryDirection::Negative).coords(), [0, 2, 3]);
    /// assert_eq!(location.sibling_or_same_on_axis(1, QueryDirection::Negative).coords(), [1, 2, 3]);
    /// assert_eq!(location.sibling_or_same_on_axis(2, QueryDirection::Negative).coords(), [1, 2, 2]);
    /// 
    /// assert_eq!(location.sibling_or_same_on_axis(0, QueryDirection::Positive).coords(), [1, 2, 3]);
    /// assert_eq!(location.sibling_or_same_on_axis(1, QueryDirection::Positive).coords(), [1, 3, 3]);
    /// assert_eq!(location.sibling_or_same_on_axis(2, QueryDirection::Positive).coords(), [1, 2, 3]);
    /// ```
    fn sibling_or_same_on_axis(&self, axis: usize, direction: QueryDirection) -> Self;
}

/// Provides methods for retrieving adacent locations
pub trait Neighbours<const D: usize>: SpaceFillingCurve<D> {
    /// Get neighbour location on axis
    /// 
    /// This method gets the neighbour location in a direction along on an axis.
    /// It is equvalent to adding 1 to or subtracting 1 from the encoded
    /// coordinate.
    /// 
    /// Unlike the [Siblings] implementation, methods in the Neighbours trait
    /// may cross cluster boundaries.
    /// 
    /// # Panics
    /// Panics if parameter 'axis' is greater than or equal to D
    /// 
    /// # Examples
    /// ```rust
    /// use insides::*;
    /// 
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    /// 
    /// assert_eq!(location.neighbour_on_axis(0, QueryDirection::Negative).coords(), [0, 2, 3]);
    /// assert_eq!(location.neighbour_on_axis(1, QueryDirection::Negative).coords(), [1, 1, 3]);
    /// assert_eq!(location.neighbour_on_axis(2, QueryDirection::Negative).coords(), [1, 2, 2]);
    /// 
    /// assert_eq!(location.neighbour_on_axis(0, QueryDirection::Positive).coords(), [2, 2, 3]);
    /// assert_eq!(location.neighbour_on_axis(1, QueryDirection::Positive).coords(), [1, 3, 3]);
    /// assert_eq!(location.neighbour_on_axis(2, QueryDirection::Positive).coords(), [1, 2, 4]);
    /// ```
    fn neighbour_on_axis(&self, axis: usize, direction: QueryDirection) -> Self;
}
