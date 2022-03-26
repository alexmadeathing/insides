
// Can be replaced with core::array::from_fn when stabilised
// https://github.com/rust-lang/rust/pull/94119
#[inline]
pub fn array_from_fn<F, T, const N: usize>(mut cb: F) -> [T; N]
where
    F: FnMut(usize) -> T,
{
    let mut idx = 0;
    [(); N].map(|_| {
        let res = cb(idx);
        idx += 1;
        res
    })
}

pub trait NumTraits: dilate::DilatableType + Ord {
    // Add methods as needed
    fn zero() -> Self;
    fn one() -> Self;
    fn bits() -> usize;
    fn lz(self) -> usize;
    fn add(self, rhs: Self) -> Self;
    fn sub(self, rhs: Self) -> Self;
    fn shl(self, amount: usize) -> Self;
    fn shr(self, amount: usize) -> Self;
    fn bit_not(self) -> Self;
    fn bit_and(self, rhs: Self) -> Self;
    fn bit_or(self, rhs: Self) -> Self;
    fn bit_xor(self, rhs: Self) -> Self;
    fn from_usize(value: usize) -> Self;
    fn to_usize(self) -> usize;
}

macro_rules! impl_num_traits {
    ($($t:ty),+) => {$(
        impl NumTraits for $t {
            #[inline]
            fn zero() -> Self {
                0
            }

            #[inline]
            fn one() -> Self {
                1
            }

            #[inline]
            fn bits() -> usize {
                Self::BITS as usize
            }

            #[inline]
            fn lz(self) -> usize {
                self.leading_zeros() as usize
            }

            #[inline]
            fn add(self, rhs: Self) -> Self {
                self + rhs
            }

            #[inline]
            fn sub(self, rhs: Self) -> Self {
                self - rhs
            }

            #[inline]
            fn shl(self, amount: usize) -> Self {
                self << amount
            }

            #[inline]
            fn shr(self, amount: usize) -> Self {
                self >> amount
            }

            #[inline]
            fn bit_not(self) -> Self {
                !self
            }

            #[inline]
            fn bit_and(self, rhs: Self) -> Self {
                self & rhs
            }

            #[inline]
            fn bit_or(self, rhs: Self) -> Self {
                self | rhs
            }

            #[inline]
            fn bit_xor(self, rhs: Self) -> Self {
                self ^ rhs
            }

            #[inline]
            fn from_usize(value: usize) -> Self {
                value as Self
            }

            #[inline]
            fn to_usize(self) -> usize {
                self as usize
            }
        }
    )+};
}

impl_num_traits!(u8, u16, u32, u64, u128, usize);
