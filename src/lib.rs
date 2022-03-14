pub mod morton_curve;

#[cfg(test)]
mod test_data;

pub use morton_curve::Morton;

#[derive(Clone, Copy, Debug)]
pub enum SearchDirection {
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
    fn sibling_or_self_on_axis(&self, axis: usize, search_direction: SearchDirection) -> Self;
}

pub trait Neighbours {
    fn neighbour_on_axis(&self, axis: usize, search_direction: SearchDirection) -> Self;
}
