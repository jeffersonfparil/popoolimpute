// // fn closure_test<F>(f:F) where F:Fn(&str) {
// pub fn closure_test<F:Fn(&str)>(f:F) {
//     f("Impution of genotype data from pool");
//     f("sequencing allele frequency information.");
// }

use std::fs::File;
use std::io::{Result, Lines, BufRead, BufReader};
use std::path::Path;

pub fn stream_file<P>(filename:P) -> Result<Lines<BufReader<File>>> where P:AsRef<Path> {
    let file = File::open(filename).expect("File not found!");
    Ok(BufReader::new(file).lines())
}

struct Locus {
    chromosome: String,
    position: u64,
    reference_allele: char,
    ATCGIDN_counts: Vec<u64>
}

impl Locus {

    fn parse_1_read_column(read: &str, reference_allele:char, quality:&str, minimum_quality:u8) -> (u64, u64, u64, u64, u64, u64, u64) {
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
        let mut i:u64 = 0;
        let mut q:u8;
        let mut read_start_marker:bool = false;
        let mut read_regex_marker:bool = false;
        let mut insertion_deletion_char:char = '+';
        let mut insertion_deletion_count:u32 = 0;
        let mut insertion_deletion_count_place_value:u32 = 1;
        // for debugging
        let mut n_read_start = 0;
        let mut n_read_end = 0;
        // iterate across reads
        // println!("qual length: {}", qual.iter().count());
        // println!("read length: {}", vec.iter().count());
        for x in vec.iter(){
            q = qual[i as usize];
            if *x == '^' {
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
                    i = i + 1; // conitnue to the next quality score
                } else if (*x=='+') | (*x=='-') {
                    insertion_deletion_char = *x;
                    read_regex_marker = true;
                    continue;
                } else if read_regex_marker {
                    if let Some(d) = x.to_digit(10) {
                        insertion_deletion_count = (insertion_deletion_count*insertion_deletion_count_place_value)
                                                    + d;
                        insertion_deletion_count_place_value = insertion_deletion_count_place_value * 10;
                    } else if insertion_deletion_count > 1 {
                        match insertion_deletion_char {
                            '+' => I = I + 1,
                            '-' => D = D + 1,
                            _ => N = N + 1
                        }
                        insertion_deletion_count = insertion_deletion_count - 1;
                    } else {
                        match insertion_deletion_char {
                            '+' => I = I + 1,
                            '-' => D = D + 1,
                            _ => N = N + 1
                        }
                        read_regex_marker = false;
                        insertion_deletion_count = 0;
                        insertion_deletion_count_place_value = 0;
                    }
                    continue;
                } else {
                    match *x {
                        'A' | 'a' => A = A + 1,
                        'T' | 't' => T = T + 1,
                        'C' | 'c' => C = C + 1,
                        'G' | 'g' => G = G + 1,
                        '*'       => D = D + 1,
                        _ => N = N + 1
                    }
                }
            }
        }
        // debugging
        // println!("Read starts: {}", n_read_start);
        // println!("Read ends: {}", n_read_end);
        // println!("i went up to: {}", i);
        // println!("qual length: {}", qual.iter().count());
        // println!("read length: {}", vec.iter().count());
        // output tupe of counts
        (A, T, C, G, I, D, N)
    }

    fn new(vec: Vec<String>, min_base_quality: u8) -> Self {
        let mut chromosome: String = vec[0].to_string();
        let mut position: u64 = vec[1].parse().unwrap();
        let mut reference_allele: char = vec[2].chars().next().expect("We're expecting a single character as reference allele!");
        let p: usize = vec.iter().count();
        let mut out: Vec<u64> = vec![];
        for i in (4..p).step_by(3) {
            // println!("{}", i);
            let (A, T, C, G, I, D, N) = Locus::parse_1_read_column(&vec[i], reference_allele, &vec[i+1], min_base_quality);
            out.push(A); out.push(T); out.push(C); out.push(G); out.push(I); out.push(D); out.push(N);
        }
        Self{chromosome: chromosome, position: position, reference_allele: reference_allele, ATCGIDN_counts: out}
    }
}

pub fn iterate_across_lines(stream:Result<Lines<BufReader<File>>>){
    let min_base_quality: u8 = 20;
    let window_size: u64 = 1;
    let mut vec: Vec<String>;
    let mut window: Vec<Locus> = vec![];
    let mut i: u64 = 0;
    for line in stream.unwrap() {
        vec = line
                .unwrap()
                .split_whitespace()
                .map(str::to_string)
                .collect();
        let locus: Locus = Locus::new(vec, min_base_quality);
        // debug
        println!("{}", window[0].chromosome);
        println!("{}", window[0].position);
        println!("{}-{}-{}-{}-{}-{}-{}<--->{}-{}-{}-{}-{}-{}-{}",
                 window[0].ATCGIDN_counts[0], window[0].ATCGIDN_counts[1], window[0].ATCGIDN_counts[2], window[0].ATCGIDN_counts[3], window[0].ATCGIDN_counts[4], window[0].ATCGIDN_counts[5], window[0].ATCGIDN_counts[6],
                 window[0].ATCGIDN_counts[7], window[0].ATCGIDN_counts[8], window[0].ATCGIDN_counts[9], window[0].ATCGIDN_counts[10], window[0].ATCGIDN_counts[11], window[0].ATCGIDN_counts[12], window[0].ATCGIDN_counts[13]
                );
        // fill-up window
        if window.len() < window_size as usize {
            window.push(locus);
        } else {
            i = i + 1;

        }
    }
}
