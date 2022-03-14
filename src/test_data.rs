use std::path::PathBuf;

use serde_json::{Value, Map};

fn test_data_path(test_group: &str, file: &str) -> PathBuf {
    let mut path = std::env::current_exe().expect("Failed to get exe path");
    path.pop(); // exe
    path.pop(); // deps
    path.pop(); // debug/release
    path.pop(); // target
    path.push("test_data");
    path.push(String::from(test_group));
    path.push(String::from(file));
    path
}

pub fn import(test_group: &str, file: &str) -> Value {
    let file = std::fs::read_to_string(test_data_path(test_group, file)).expect("Failed to read test data file");
    let json: Value = serde_json::from_str(file.as_str()).unwrap_or_default();
    json
}

pub fn iter_object_array(object_array: &Vec<Value>) -> impl Iterator<Item = &Map<String, Value>> {
    object_array.iter().map(|tc| tc.as_object().expect("Expected object in object array"))
}

pub fn parse_radix(str: &str) -> (u32, &str) {
    assert!(!str.starts_with("0X"), "Number radices should be lower case - 0x{} instead of {}", str.strip_prefix("0X").unwrap(), str);
    assert!(!str.starts_with("0O"), "Number radices should be lower case - 0o{} instead of {}", str.strip_prefix("0O").unwrap(), str);
    assert!(!str.starts_with("0B"), "Number radices should be lower case - 0b{} instead of {}", str.strip_prefix("0B").unwrap(), str);
    if str.starts_with("0x") {
        (16, str.strip_prefix("0x").unwrap())
    } else if str.starts_with("0o") {
        (8, str.strip_prefix("0o").unwrap())
    } else if str.starts_with("0b") {
        (2, str.strip_prefix("0b").unwrap())
    } else {
        (10, str)
    }
}

pub fn parse_str_num(str: &str, coord_max: u128, index_max: u128) -> u128 {
    if str == "coord_max" {
        coord_max
    }
    else if str == "index_max" {
        index_max
    } else {
        let (radix, remaining) = parse_radix(str);
        let value = u128::from_str_radix(remaining, radix).expect("Failed to parse number expressed as string: {str}");
        value
    }
}

pub fn str_num_array_to_array<const N: usize>(array: &Value, coord_max: u128, index_max: u128) -> [u128; N] {
    assert!(array.is_array(), "Expected array containing {N} numbers expressed as strings");
    let array = array.as_array().unwrap();
    assert!(array.len() == N, "Expected array containing {N} numbers expressed as strings");

    let mut i = 0;
    [(); N].map(|_| {
        let r = str_num(&array[i], coord_max, index_max);
        i += 1;
        r
    })
}

pub fn str_num(str_num: &Value, coord_max: u128, index_max: u128) -> u128 {
    assert!(str_num.is_string(), "Numbers should be expressed as strings to allow for large integers and different radices");
    let str_num = str_num.as_str().unwrap();
    parse_str_num(str_num, coord_max, index_max)
}
