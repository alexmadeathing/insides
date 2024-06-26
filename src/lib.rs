#![no_std]
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

use core::{fmt::Debug, hash::Hash};

#[doc(hidden)]
pub use dilate::*;

/// Morton curve implementation
pub mod morton;
pub use crate::morton::Morton;

/// Hilbert curve implementation
pub mod hilbert;
pub use crate::hilbert::Hilbert;

pub(crate) mod internal;
use internal::NumTraits;

/// Trait wrapper for coordinates
pub trait CurveCoord: dilate::DilatableType + NumTraits {}

impl CurveCoord for u8 {}
impl CurveCoord for u16 {}
impl CurveCoord for u32 {}
impl CurveCoord for u64 {}
impl CurveCoord for u128 {}
impl CurveCoord for usize {}

/// Trait wrapper for indices
pub trait CurveIndex: dilate::DilatableType + NumTraits {}

impl CurveIndex for u8 {}
impl CurveIndex for u16 {}
impl CurveIndex for u32 {}
impl CurveIndex for u64 {}
impl CurveIndex for u128 {}
impl CurveIndex for usize {}

/// Direction to search within an axis when making queries
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QueryDirection {
    /// Search in a negative direction along an axis relative to the source location
    Negative,

    /// Search in a positive direction along an axis relative to the source location
    Positive,
}

impl QueryDirection {
    #[inline]
    fn mix<I>(self, index: I, axis: usize) -> I
    where
        I: NumTraits,
    {
        if self == Self::Positive {
            index.bit_or(I::one().shl(axis))
        } else {
            index
        }
    }
}

/// Provides conversion to and from coordinates
// I'd really love to get rid of the generic parameter here but I think it's waiting on:
// https://github.com/rust-lang/rust/issues/76560
pub trait SpaceFillingCurve<const D: usize>: Sized + Ord + Copy + Default + Debug + Hash {
    /// Coordinate type
    type Coord: CurveCoord;

    /// Index type
    type Index: CurveIndex;

    /// Number of bits available in a single coordinate
    const COORD_BITS: usize;

    /// Maximum value for a single coordinate
    const COORD_MAX: Self::Coord;

    /// Number of bits available in an encoded index
    const INDEX_BITS: usize;

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
    /// assert_eq!(location.sibling_on_axis_toggle(0).coords(), [0, 2, 3]);
    /// assert_eq!(location.sibling_on_axis_toggle(1).coords(), [1, 3, 3]);
    /// assert_eq!(location.sibling_on_axis_toggle(2).coords(), [1, 2, 2]);
    /// ```
    fn sibling_on_axis_toggle(&self, axis: usize) -> Self;

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
    /// assert_eq!(location.sibling_on_axis(0, QueryDirection::Negative).coords(), [0, 2, 3]);
    /// assert_eq!(location.sibling_on_axis(1, QueryDirection::Negative).coords(), [1, 2, 3]);
    /// assert_eq!(location.sibling_on_axis(2, QueryDirection::Negative).coords(), [1, 2, 2]);
    ///
    /// assert_eq!(location.sibling_on_axis(0, QueryDirection::Positive).coords(), [1, 2, 3]);
    /// assert_eq!(location.sibling_on_axis(1, QueryDirection::Positive).coords(), [1, 3, 3]);
    /// assert_eq!(location.sibling_on_axis(2, QueryDirection::Positive).coords(), [1, 2, 3]);
    /// ```
    fn sibling_on_axis(&self, axis: usize, direction: QueryDirection) -> Self;

    /// Get sibling from axis bits within local cluster
    ///
    /// This method gets, within the local cluster, the sibling identified by a set
    /// of D axis bits. Where the least significant bit represents the X (or 0th) axis and
    /// each successive more significant bit represents the next higher axis. A value of 1 in
    /// an axis bit identifies the sibling in the positive direction along that axis
    /// and a value of 0 identifies the sibling in the negative direction.
    ///
    /// This combination of bits is sufficient to identify any sibling within the
    /// local cluster.
    ///
    /// The local cluster is the grouping of 2<sup>D</sup> adjacent locations within
    /// the coordinate space, where the cluster origin lies at even coordinates.
    ///
    /// # Panics
    /// Panics if parameter 'axis_bits' contains set bits in positions beyond D
    ///
    /// # Examples
    /// ```rust
    /// use insides::*;
    ///
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    ///
    /// assert_eq!(location.sibling_from_bits(0b000).coords(), [0, 2, 2]);
    /// assert_eq!(location.sibling_from_bits(0b001).coords(), [1, 2, 2]);
    /// assert_eq!(location.sibling_from_bits(0b010).coords(), [0, 3, 2]);
    /// assert_eq!(location.sibling_from_bits(0b011).coords(), [1, 3, 2]);
    /// assert_eq!(location.sibling_from_bits(0b100).coords(), [0, 2, 3]);
    /// assert_eq!(location.sibling_from_bits(0b101).coords(), [1, 2, 3]);
    /// assert_eq!(location.sibling_from_bits(0b110).coords(), [0, 3, 3]);
    /// assert_eq!(location.sibling_from_bits(0b111).coords(), [1, 3, 3]);
    /// ```
    fn sibling_from_bits(&self, axis_bits: Self::Index) -> Self;

    /// Get sibling from axes within local cluster
    ///
    /// This method gets, within the local cluster, the sibling identified by a set
    /// of axis query directions. Where the first query direction in the `axes`
    /// array identifies the sibling on the X (or 0th) axis and each successive
    /// query direction identifies the sibling on the next higher axis.
    ///
    /// This combination of query directions is sufficient to identify any sibling
    /// within the local cluster.
    ///
    /// The local cluster is the grouping of 2<sup>D</sup> adjacent locations within
    /// the coordinate space, where the cluster origin lies at even coordinates.
    ///
    /// # Examples
    /// ```rust
    /// use insides::*;
    /// use insides::QueryDirection::{Negative, Positive};
    ///
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    ///
    /// assert_eq!(location.sibling_on_axes([Negative, Negative, Negative]).coords(), [0, 2, 2]);
    /// assert_eq!(location.sibling_on_axes([Negative, Negative, Positive]).coords(), [0, 2, 3]);
    /// assert_eq!(location.sibling_on_axes([Negative, Positive, Negative]).coords(), [0, 3, 2]);
    /// assert_eq!(location.sibling_on_axes([Negative, Positive, Positive]).coords(), [0, 3, 3]);
    /// assert_eq!(location.sibling_on_axes([Positive, Negative, Negative]).coords(), [1, 2, 2]);
    /// assert_eq!(location.sibling_on_axes([Positive, Negative, Positive]).coords(), [1, 2, 3]);
    /// assert_eq!(location.sibling_on_axes([Positive, Positive, Negative]).coords(), [1, 3, 2]);
    /// assert_eq!(location.sibling_on_axes([Positive, Positive, Positive]).coords(), [1, 3, 3]);
    /// ```
    #[inline]
    fn sibling_on_axes(&self, axes: [QueryDirection; D]) -> Self {
        self.sibling_from_bits(
            (0..D)
                .into_iter()
                .fold(Self::Index::zero(), |b, i| axes[i].mix(b, i)),
        )
    }
}

/// Provides methods for retrieving adacent locations
pub trait Neighbours<const D: usize>: SpaceFillingCurve<D> {
    /// Get neighbour location on axis
    ///
    /// This method gets the neighbour location in a direction along an axis.
    /// It is equvalent to adding 1 to or subtracting 1 from the encoded
    /// coordinate.
    ///
    /// If the neighbour would cause a coordinate to wrap, this method returns
    /// None instead. For a wrapping version, please see
    /// [neighbour_on_axis_wrapping()](Neighbours::neighbour_on_axis_wrapping()).
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
    /// assert_eq!(location.neighbour_on_axis(0, QueryDirection::Negative).unwrap().coords(), [0, 2, 3]);
    /// assert_eq!(location.neighbour_on_axis(1, QueryDirection::Negative).unwrap().coords(), [1, 1, 3]);
    /// assert_eq!(location.neighbour_on_axis(2, QueryDirection::Negative).unwrap().coords(), [1, 2, 2]);
    ///
    /// assert_eq!(location.neighbour_on_axis(0, QueryDirection::Positive).unwrap().coords(), [2, 2, 3]);
    /// assert_eq!(location.neighbour_on_axis(1, QueryDirection::Positive).unwrap().coords(), [1, 3, 3]);
    /// assert_eq!(location.neighbour_on_axis(2, QueryDirection::Positive).unwrap().coords(), [1, 2, 4]);
    ///
    /// let edge_location = Morton::<Expand<u16, 2>, 2>::from_coords([0, u16::MAX]);
    /// assert_eq!(edge_location.neighbour_on_axis(0, QueryDirection::Negative), None);
    /// assert_eq!(edge_location.neighbour_on_axis(1, QueryDirection::Positive), None);
    /// ```
    fn neighbour_on_axis(&self, axis: usize, direction: QueryDirection) -> Option<Self>;

    /// Get neighbour location on axis with wrapping
    ///
    /// This method gets the neighbour location in a direction along an axis.
    /// It is equvalent to adding 1 to or subtracting 1 from the encoded
    /// coordinate.
    ///
    /// This method allows coordinates to wrap between zero and
    /// [COORD_MAX](SpaceFillingCurve::COORD_MAX). For a non-wrapping version,
    /// please see
    /// [neighbour_on_axis()](Neighbours::neighbour_on_axis()).
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
    /// assert_eq!(location.neighbour_on_axis_wrapping(0, QueryDirection::Negative).coords(), [0, 2, 3]);
    /// assert_eq!(location.neighbour_on_axis_wrapping(1, QueryDirection::Negative).coords(), [1, 1, 3]);
    /// assert_eq!(location.neighbour_on_axis_wrapping(2, QueryDirection::Negative).coords(), [1, 2, 2]);
    ///
    /// assert_eq!(location.neighbour_on_axis_wrapping(0, QueryDirection::Positive).coords(), [2, 2, 3]);
    /// assert_eq!(location.neighbour_on_axis_wrapping(1, QueryDirection::Positive).coords(), [1, 3, 3]);
    /// assert_eq!(location.neighbour_on_axis_wrapping(2, QueryDirection::Positive).coords(), [1, 2, 4]);
    ///
    /// let edge_location = Morton::<Expand<u16, 2>, 2>::from_coords([0, u16::MAX]);
    /// assert_eq!(edge_location.neighbour_on_axis_wrapping(0, QueryDirection::Negative).coords(), [u16::MAX, u16::MAX]);
    /// assert_eq!(edge_location.neighbour_on_axis_wrapping(1, QueryDirection::Positive).coords(), [0, 0]);
    /// ```
    fn neighbour_on_axis_wrapping(&self, axis: usize, direction: QueryDirection) -> Self;

    /// Get neighbour location on corner
    ///
    /// This method gets the neighbour location in a direction along each axis. This results finding the
    /// adjecant diagonal locations (the corners). It is equvalent to adding 1 to or subtracting 1 from
    /// each axis of the coordinates.
    ///
    /// If the neighbour would cause a coordinate to wrap, this method returns
    /// None instead. For a wrapping version, please see
    /// [neighbour_on_corner_wrapping()](Neighbours::neighbour_on_corner_wrapping()).
    ///
    /// Unlike the [Siblings] implementation, methods in the Neighbours trait
    /// may cross cluster boundaries.
    ///
    /// # Panics
    /// Panics if parameter 'axis' is greater than or equal to D
    ///
    /// # Examples
    /// ```rust
    /// use insides::{*, QueryDirection::*};
    ///
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    ///
    /// assert_eq!(location.neighbour_on_corner([Negative, Negative, Negative]).unwrap().coords(), [0, 1, 2]);
    /// assert_eq!(location.neighbour_on_corner([Positive, Positive, Positive]).unwrap().coords(), [2, 3, 4]);
    /// 
    /// assert_eq!(location.neighbour_on_corner([Positive, Negative, Negative]).unwrap().coords(), [2, 1, 2]);
    /// assert_eq!(location.neighbour_on_corner([Negative, Positive, Negative]).unwrap().coords(), [0, 3, 2]);
    /// assert_eq!(location.neighbour_on_corner([Negative, Negative, Positive]).unwrap().coords(), [0, 1, 4]);
    ///
    /// let edge_location = Morton::<Expand<u16, 2>, 2>::from_coords([0, u16::MAX]);
    /// assert_eq!(edge_location.neighbour_on_corner([Negative, Negative]), None);
    /// assert_eq!(edge_location.neighbour_on_corner([Positive, Positive]), None);
    /// ```
    fn neighbour_on_corner(&self, axes: [QueryDirection; D]) -> Option<Self>;

    /// Get neighbour location on corner with wrapping
    ///
    /// This method gets the neighbour location in a direction along each axis. This results finding the
    /// adjecant diagonal locations (the corners). It is equvalent to adding 1 to or subtracting 1 from
    /// each axis of the coordinates.
    ///
    /// This method allows coordinates to wrap between zero and
    /// [COORD_MAX](SpaceFillingCurve::COORD_MAX). For a non-wrapping version,
    /// please see
    /// [neighbour_on_corner()](Neighbours::neighbour_on_corner()).
    ///
    /// Unlike the [Siblings] implementation, methods in the Neighbours trait
    /// may cross cluster boundaries.
    ///
    /// # Panics
    /// Panics if parameter 'axis' is greater than or equal to D
    ///
    /// # Examples
    /// ```rust
    /// use insides::{*, QueryDirection::*};
    ///
    /// let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);
    ///
    /// assert_eq!(location.neighbour_on_corner_wrapping([Negative, Negative, Negative]).coords(), [0, 1, 2]);
    /// assert_eq!(location.neighbour_on_corner_wrapping([Positive, Positive, Positive]).coords(), [2, 3, 4]);
    /// 
    /// assert_eq!(location.neighbour_on_corner_wrapping([Positive, Negative, Negative]).coords(), [2, 1, 2]);
    /// assert_eq!(location.neighbour_on_corner_wrapping([Negative, Positive, Negative]).coords(), [0, 3, 2]);
    /// assert_eq!(location.neighbour_on_corner_wrapping([Negative, Negative, Positive]).coords(), [0, 1, 4]);
    ///
    /// let edge_location = Morton::<Expand<u16, 2>, 2>::from_coords([0, u16::MAX]);
    /// assert_eq!(edge_location.neighbour_on_corner_wrapping([Negative, Negative]).coords(), [u16::MAX, u16::MAX - 1]);
    /// assert_eq!(edge_location.neighbour_on_corner_wrapping([Positive, Positive]).coords(), [1, 0]);
    /// ```
    fn neighbour_on_corner_wrapping(&self, axes: [QueryDirection; D]) -> Self;
}
