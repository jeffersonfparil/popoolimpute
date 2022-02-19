// // fn closure_test<F>(f:F) where F:Fn(&str) {
// pub fn closure_test<F:Fn(&str)>(f:F) {
//     f("Impution of genotype data from pool");
//     f("sequencing allele frequency information.");
// }

use std::fs::File;
use std::io::{Result, Lines, BufRead, BufReader};
use std::path::Path;

pub fn stream_file<P>(filename:P) -> Result<Lines<BufReader<File>>>
where P:AsRef<Path>, {
    let file = File::open(filename).expect("File not found!");
    Ok(BufReader::new(file).lines())
}

// pub fn print_line_by_line(stream:Result<Lines<BufReader<File>>>){
//     for line in stream.unwrap() {
//         println!("{}", line.unwrap());
//     }
// }

// pub fn print_column_by_column_per_line(stream:Result<Lines<BufReader<File>>>){
//     for line in stream.unwrap() {
//         for column in line.unwrap().split_whitespace() {
//             println!("{}", column);
//         }
//     }
// }

fn parse_read(read:&str, reference_allele:char, quality:&str, minimum_quality:u8) -> (u64, u64, u64, u64, u64, u64, u64) {
    // qualities
    let mut qual:Vec<u8> = vec![];
    let vec_qual:Vec<char> = quality.chars().collect();
    for q in vec_qual.iter() {
        qual.push(*q as u8 - 33);
    }
    // reads
    let vec:Vec<char> = read.chars().collect();
    let mut A:u64 = 0;
    let mut T:u64 = 0;
    let mut C:u64 = 0;
    let mut G:u64 = 0;
    let mut I:u64 = 0;
    let mut D:u64 = 0;
    let mut N:u64 = 0;
    let mut i:usize = 0;
    let mut q:u8;
    let mut read_start_marker:bool = false;
    let mut read_regex_marker:bool = false;
    let mut insertion_deletion_char:char = '+';
    let mut insertion_deletion_count:u32 = 0;
    let mut insertion_deletion_count_place_value:u32 = 1;

    let mut n_read_start = 0;
    let mut n_read_end = 0;

    for x in vec.iter(){
        q = qual[i];
        if *x == '^' {
            // println!("x = {}", *x);
            read_start_marker = true;
            n_read_start = n_read_start + 1;
            continue;
        } else if read_start_marker {
            read_start_marker = false;
            continue;
        } else if *x == '$' {
            n_read_end = n_read_end + 1;
            continue;
        } else if q >= minimum_quality {
            if (*x=='.') | (*x==',') {
                match reference_allele {
                    'A' => A = A + 1,
                    'T' => T = T + 1,
                    'C' => C = C + 1,
                    'G' => G = G + 1,
                    _ => N = N + 1
                }
            } else if (*x=='+') | (*x=='-') {
                insertion_deletion_char = *x;
                read_regex_marker = true;
                continue;
            } else if read_regex_marker {
                if let Some(d) = x.to_digit(10) {
                    insertion_deletion_count = (insertion_deletion_count*insertion_deletion_count_place_value)
                                                + d;
                    insertion_deletion_count_place_value = insertion_deletion_count_place_value * 10;
                } else {
                    match insertion_deletion_char {
                        '+' => I = I + (insertion_deletion_count as u64),
                        '-' => D = D + (insertion_deletion_count as u64),
                        _ => N = N + (insertion_deletion_count as u64),
                    }
                    read_regex_marker = false;
                    insertion_deletion_count_place_value = 1;
                }
            } else {
                match *x {
                    'A' | 'a' => A = A + 1,
                    'T' | 't' => T = T + 1,
                    'C' | 'c'  => C = C + 1,
                    'G' | 'g' => G = G + 1,
                    _ => N = N + 1
                }
            }
            i = i + 1;
        } else {
            i = i + 1;
        }
    }
    println!("Read starts: {}", n_read_start);
    println!("Read ends: {}", n_read_end);
    println!("i went up to: {}", i);
    println!("qual length: {}", qual.iter().count());
    println!("read length: {}", vec.iter().count());
    (A, T, C, G, I, D, N)
}


pub fn parse_line(stream:Result<Lines<BufReader<File>>>){
    let mut vec: Vec<String>;
    for line in stream.unwrap() {
        vec = line
                .unwrap()
                .split_whitespace()
                .map(str::to_string)
                .collect();
        // println!("{}", vec.iter().count());
        // println!("{}", vec[4]);
        // let chars: Vec<char> = vec[4].chars().collect();
        // println!("{}", chars[4]);
        let (a, t, c, g, i, d, n) = parse_read(&vec[4], 'A', &vec[5], 0);
        println!("{}", n);
    }
}