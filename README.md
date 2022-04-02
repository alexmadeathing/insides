[![Crates.io](https://img.shields.io/crates/d/insides.svg)](https://crates.io/crates/insides)
[![Anti-Capitalist Software License (v 1.4)](https://img.shields.io/badge/license-Anti--Capitalist%20(v%201.4)-brightgreen)](LICENSE.md)
[![alexmadeathing](https://circleci.com/gh/alexmadeathing/insides.svg?style=shield)](https://app.circleci.com/pipelines/github/alexmadeathing/insides?filter=all)
[![rustc 1.51+]][Rust 1.51]

# WARNING
This library is in an alpha stage of development. Its interface may be subject to change.

For migration notes, please see: https://github.com/alexmadeathing/insides/releases

# insides
A compact, high performance space filling curve library for Rust.

This library provides an abstract interface to generalise interactions with
space filling curves as well as a morton curve implementation and
supporting manipulation methods.

# Supported Curves
* Morton - A morton, or Z-order, curve implementation

We currently only support Morton Encoding, but the interface will support
other curve implementations, such as Hilbert - and we have plans to include them soon.

# Features
* High performance - Ready to use in performance sensitive contexts
* N-dimensional - Suitable for multi-dimensional applications (up to 16 dimensions under certain conditions)
* Type safe - Multiple input types with known output types (supports `u8`, `u16`, `u32`, `u64`, `u128`, `usize`)
* `no_std` - Suitable for embedded devices (additional standard library features can be enabled via the `std` feature)
* Extensible - Flexible trait based implementation
* Minimal dependencies - Release build depends only on [dilate](https://github.com/alexmadeathing/dilate) for integer dilation

# Getting Started
First, link insides into your project's cargo.toml.

Check for the latest version at [crates.io](https://crates.io/crates/insides):
```toml
[dependencies]
insides = "0.1.0"
```

Next, import insides into your project and try out some of the features:

```rust
use insides::*;

// Create a 3D morton location using u16 coordinate indices
// At the moment, we have to specify the number of dimensions twice, sorry!
// (this will change with improvements to Rust const generics)
let location = Morton::<Expand<u16, 3>, 3>::from_coords([1, 2, 3]);

// Access raw morton location index
// Useful as a map key or array index!
// In this case, we're using a morton curve, so we know the underlying bit
// layout. If we used a different type of curve, the index would be
// different.
assert_eq!(location.index(), 0b110101);

// Get neighbour in the negative X axis
assert_eq!(location.neighbour_on_axis(0, QueryDirection::Negative).coords(), [0, 2, 3]);

// Get neighbour in the positive Y axis
assert_eq!(location.neighbour_on_axis(1, QueryDirection::Positive).coords(), [1, 3, 3]);

// Get sibling on Z axis
// Its use is subtle and may not be apparent at first glance. Imagine the
// world is divided into a grid where each cell contains a cluster of 8
// siblings (in the 3D case). This method operates upon that grid. This is
// mostly useful when applied to tree-like structures.
let sibling = location.sibling_on_axis(2);
assert_eq!(sibling.coords(), [1, 2, 2]);

// Calling again from the sibling gets our original location
assert_eq!(sibling.sibling_on_axis(2), location);

// This gets, either the sibling, or the same location, on the Z axis,
// depending on which is further in the query direction.
assert_eq!(location.sibling_or_same_on_axis(2, QueryDirection::Positive).coords(), [1, 2, 3]);
```

For more detailed info, please see the [code reference](https://docs.rs/insides/latest/insides/).

# Roadmap
Please refer to the [Roadmap to V1.0](https://github.com/alexmadeathing/insides/discussions/2) discussion.

# Contributing
Contributions are most welcome.

For bugs reports, please [submit a bug report](https://github.com/alexmadeathing/insides/issues/new?assignees=&labels=bug&template=bug_report.md&title=).

For feature requests, please [submit a feature request](https://github.com/alexmadeathing/insides/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).

If you have ideas and want to contribute directly, please start by creating an [Idea discussion](https://github.com/alexmadeathing/insides/discussions/new) in the discussions area. Allow others to comment prior to committing to doing the work. When all parties agree on the design, the work may begin. When your code is ready to be published, please submit a pull request referring back to your Idea discussion. We are unlikely to accept a pull request that has not gone through this process, unless it is for a very small change.

# License

insides is licensed under the [Anti-Capitalist Software License (v 1.4)](https://github.com/alexmadeathing/insides/blob/main/LICENSE.md). This means it is free and open source for use by individuals and organizations that do not operate by capitalist principles.

Unless explicitly stated, your contributions will be incorporated under this license.
