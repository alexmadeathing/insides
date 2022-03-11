[![Crates.io](https://img.shields.io/crates/d/insides.svg)](https://crates.io/crates/insides)
[![Anti-Capitalist Software License (v 1.4)](https://img.shields.io/badge/license-Anti--Capitalist%20(v%201.4)-brightgreen)](LICENSE.md)
[![alexmadeathing](https://circleci.com/gh/alexmadeathing/insides.svg?style=shield)](https://app.circleci.com/pipelines/github/alexmadeathing/insides?filter=all)

# WARNING
This library is in an alpha stage of development. Its interface may be subject to change.

For migration notes, please see: https://github.com/alexmadeathing/insides/releases

# insides
A compact, high performance space filling curve library for Rust.

# Supported Curves

# Goals
* High performance - Ready to use in performance sensitive contexts
* Multiple types - Supports `u8`, `u16`, `u32`, `u64`, `u128`, `usize` (signed versions not yet planned)
* N-dimensional - Suitable for multi-dimensional applications (up to 16 dimensions under certain conditions)
* Trait based implementation - Conforms to standard Rust implementation patterns
* Minimal dependencies - Depends only on [dilate](https://github.com/alexmadeathing/dilate) for integer dilation

# Getting Started
First, link insides into your project's cargo.toml.

Check for the latest version at [crates.io](https://crates.io/crates/insides):
```
[dependencies]
insides = "0.1.0"
```

Next, import insides into your project and try out some of the features:

```
use insides::*;
```

For more detailed info, please see the [code reference](https://docs.rs/insides/latest/insides/).

# Roadmap
Please refer to the [Roadmap to V1.0](https://github.com/alexmadeathing/insides/discussions/2) discussion.

# Contributing
Contributions are most welcome.

For bugs reports, please [submit a bug report](https://github.com/alexmadeathing/insides/issues/new?assignees=&labels=bug&template=bug_report.md&title=).

For feature requests, please [submit a feature request](https://github.com/alexmadeathing/insides/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).

If you have ideas and want to contribute directly, please start by creating an [Idea discussion](https://github.com/alexmadeathing/insides/discussions/new) in the discussions area. Allow others to comment prior to committing to doing the work. When all parties agree on the design, the work may begin. When your code is ready to be published, please submit a pull request referring back to your Idea discussion. We are unlikely to accept a pull request that has not gone through this process, unless it is for a very small change.

# References and Acknowledgments

# License

insides is licensed under the [Anti-Capitalist Software License (v 1.4)](https://github.com/alexmadeathing/insides/blob/main/LICENSE.md). This means it is free and open source for use by individuals and organizations that do not operate by capitalist principles.

Unless explicitly stated, your contributions will be incorporated under this license.
