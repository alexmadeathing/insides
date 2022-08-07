use insides::{Hilbert, Fixed, SpaceFillingCurve};

pub fn black_box<T>(dummy: T) -> T {
    unsafe {
        let ret = std::ptr::read_volatile(&dummy);
        std::mem::forget(dummy);
        ret
    }
}

fn main() {
    const D: usize = 4;
    let num: usize = 256usize.pow(D as u32);
    let mut total_coords = [0; D];
    let mut count = 0;
    for i in 0..num {
        let coords = Hilbert::<Fixed<u64, D>, D>::from_index(i as u64).coords();
        let mut any_non_zero = false;
        for j in 0..D {
            total_coords[j] = coords[j].max(total_coords[j]);
            any_non_zero = any_non_zero || coords[j] > 0;
        }
        if any_non_zero {
            count += 1;
        }
    }
    println!("Result: {total_coords:?} ({count})");
}
