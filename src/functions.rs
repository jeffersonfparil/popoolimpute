// // fn closure_test<F>(f:F) where F:Fn(&str) {
// pub fn closure_test<F:Fn(&str)>(f:F) {
//     f("Impution of genotype data from pool");
//     f("sequencing allele frequency information.");
// }


use std::fs::File;
use std::io::{Result, Lines, BufRead, BufReader};
use std::path::Path;
use ndarray::{Array1, Array2, Axis, arr1, arr2};

pub fn stream_file<P>(filename:P) -> Result<Lines<BufReader<File>>> where P:AsRef<Path> {
    let file = File::open(filename).expect("File not found!");
    Ok(BufReader::new(file).lines())
}

pub fn count_columns(mut lines: Lines<BufReader<File>>) -> usize {
    let line: String = lines.next().unwrap().unwrap();
    let vec_string: Vec<String> = line.split_whitespace()
                    .map(str::to_string).collect();
    let n: usize = vec_string.len() as usize;
    n
}

trait PileupLocus {
    fn parse_vec_string(&self, min_base_quality: u8) -> (String, u64, char, Vec<u64>);
}

impl PileupLocus for Vec<String> {

    fn parse_vec_string(&self, min_base_quality: u8) -> (String, u64, char, Vec<u64>) {
        let chromosome: String = self[0].to_string();
        let position: u64 = self[1].parse().unwrap();
        let reference_allele: char = self[2].chars().next().expect("We're expecting a single character as reference allele!");
        let p: usize = self.iter().count();
        let mut counts: Vec<u64> = vec![];

        for i in (4..p).step_by(3) {
            // println!("{}", i);
            let reads = &self[i];
            let qualities = &self[i+1];
            // reads
            let vec:Vec<char> = reads.chars().collect();
            let mut A:u64 = 0;
            let mut T:u64 = 0;
            let mut C:u64 = 0;
            let mut G:u64 = 0;
            let mut I:u64 = 0;
            let mut D:u64 = 0;
            let mut N:u64 = 0;
            let mut idx_qual:u64 = 0;
            let mut q:u8;
            let mut read_start_marker:bool = false;
            let mut read_regex_marker:bool = false;
            let mut insertion_deletion_char:char = '+';
            let mut insertion_deletion_count:u32 = 0;
            let mut insertion_deletion_count_place_value:u32 = 1;
            // for debugging
            let mut n_read_start = 0;
            let mut n_read_end = 0;
            // qualities
            let mut qual:Vec<u8> = vec![];
            let vec_qual:Vec<char> = qualities.chars().collect();
            for q in vec_qual.iter() {
                qual.push(*q as u8 - 33);
            }
            // iterate across reads
            // println!("qual length: {}", qual.iter().count());
            // println!("read length: {}", vec.iter().count());
            for x in vec.iter(){
                q = qual[idx_qual as usize];
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
                } else if q >= min_base_quality {
                    if (*x=='.') | (*x==',') {
                        match reference_allele {
                            'A' => A = A + 1,
                            'T' => T = T + 1,
                            'C' => C = C + 1,
                            'G' => G = G + 1,
                            _ => N = N + 1
                        }
                        idx_qual = idx_qual + 1; // conitnue to the next quality score
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


            counts.push(A); counts.push(T); counts.push(C); counts.push(G); counts.push(I); counts.push(D); counts.push(N);
        }
        (chromosome, position, reference_allele, counts)
    }

}

struct Window {
    chromosomes: Vec<String>,
    positions: Vec<u64>,
    reference_alleles: Vec<char>,
    count_matrix: Array2<f64> // dimensions: p loci x n pools
}

impl Window {
    fn dim(&self) -> (usize, usize){
        self.count_matrix.dim()
    }

    // fn slide(&mut self, chromosome: String, position: u64, reference_allele: char, counts: Vec<u64>) -> i64 {
    //     let (p, n) = self.dim();
    //     let mut chromosomes: Vec<String> = vec![chromosome]; chromosomes.extend(self.chromosomes[1..p-1].to_vec());
    //     let mut positions: Vec<u64> = vec![position]; positions.extend(self.positions[1..p-1].to_vec());
    //     let mut reference_alleles: Vec<char> = vec![reference_allele]; reference_alleles.extend(self.reference_alleles[1..p-1].to_vec());
    //     let out: Window = Window {
    //         chromosomes: chromosomes,
    //         positions: positions,
    //         reference_alleles: reference_alleles,
    //         count_matrix: self.count_matrix[1..p-1].append(count_matrix),
    //     }
    //     0
    // }
}

pub fn iterate_across_lines(stream:Result<Lines<BufReader<File>>>, n: usize){
    const NUMBER_OF_ALLELES: usize = 7;
    let min_base_quality: u8 = 20;
    let window_size: usize = 3;
    let mut locus_counter: u64 = 0;
    let mut vec: Vec<String>;
    let mut chromosomes: Vec<String> = Vec::new(); let chromosomes_ref = &mut chromosomes;
    let mut positions: Vec<u64> = Vec::new(); let positions_ref = &mut positions;
    let mut reference_alleles: Vec<char> = Vec::new(); let reference_alleles_ref = &mut reference_alleles;
    let mut data: Vec<u64> = Vec::new(); let data_ref = &mut data; // create a reference to the vector so we can mutate it within our for-loop
    let mut count_matrix = Array2::<f64>::zeros((window_size, n)); let count_matrix_ref = &mut count_matrix;
    for line in stream.unwrap() {
        locus_counter += 1;
        println!("{}: {} > {}", n, locus_counter, window_size);
        vec = line
                .unwrap()
                .split_whitespace()
                .map(str::to_string)
                .collect();
        let (chromosome, position, reference_allele, counts) = vec.parse_vec_string(min_base_quality);
        chromosomes_ref.extend([chromosome]);
        positions_ref.extend([position]);
        reference_alleles_ref.extend([reference_allele]);
        // data_ref.extend(counts);
        // println!("{:?}", data_ref);
        if locus_counter < window_size as u64 {
            for i in 0..n {
                count_matrix_ref.slice_mut(s![locus_counter as usize, i]).fill(counts[i] as f64);
                let x = count_matrix_ref[[locus_counter as usize, i]];
                println!("{}", x);
            }
        }
        if locus_counter == window_size as u64 {
            let t = data_ref.len();
            let n = window_size;
            let p = t / n;
            // println!("n={}; p={}", n, p);
            // count_matrix_ref = Array2::from_shape_vec((p, n), data_ref.clone()).unwrap();
            // println!("{:?}", chromosomes_ref);
            // x

        }

        if locus_counter > window_size as u64 {
            // encapsulate, impute and slide
            // encapsulate
            // let mut window = Window {
            //     chromosomes: chromosome,
            //     positions: Vec<u64>,
            //     reference_alleles: Vec<char>,
            //     count_matrix: Array2<f64> // dimensions: p loci x n pools
            // };

            // impute
            // slide
            println!("HERE");
            println!("{:?}", count_matrix_ref);
            println!("{:?}", chromosomes_ref);
            println!("{:?}", positions_ref);
            println!("{:?}", reference_alleles_ref);
        };
    }

}
