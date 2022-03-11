pub mod morton_curve;

pub use morton_curve::MortonCurve;

#[derive(Clone, Copy, Debug)]
pub enum SearchDirection {
    Positive,
    Negative,
}

// I'd really love to get rid of the generic parameter here but I think it's waiting on:
// https://github.com/rust-lang/rust/issues/76560
pub trait Coords<const D: usize> {
    type Coord;
    type Index;

    fn from_coords(coords: [Self::Coord; D]) -> Self::Index;
    fn coords(index: Self::Index) -> [Self::Coord; D];
}

pub trait Siblings {
    type Index;

    fn sibling_on_axis(index: Self::Index, axis: usize) -> Self::Index;
    fn sibling_or_self_on_axis(index: Self::Index, axis: usize, search_direction: SearchDirection) -> Self::Index;
}

pub trait Neighbours {
    type Index;

    fn neighbour_on_axis(index: Self::Index, axis: usize, search_direction: SearchDirection) -> Self::Index;
}

#[cfg(test)]
mod tests {
}
