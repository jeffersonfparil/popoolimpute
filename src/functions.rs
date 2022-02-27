// // fn closure_test<F>(f:F) where F:Fn(&str) {
// pub fn closure_test<F:Fn(&str)>(f:F) {
//     f("Impution of genotype data from pool");
//     f("sequencing allele frequency information.");
// }


use std::fs::File;
use std::io::{Result, Lines, BufRead, BufReader};
use std::path::Path;
use ndarray::{Array1, Array2, Axis, arr1, arr2};
use std::mem::replace;

const NUMBER_OF_ALLELES: usize = 7;


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

struct Window<'a> {
    chromosomes: &'a mut Vec<String>,
    positions: &'a mut Vec<u64>,
    reference_alleles: &'a mut Vec<char>,
    count_matrix: &'a mut Array2<f64> // dimensions: n (loci x 7 alleles) x p pools
}

impl Window<'_> {
    fn dim(&self) -> (usize, usize, usize){
        let (n, p) = self.count_matrix.dim();
        let n_loci: usize = n / NUMBER_OF_ALLELES;
        (n, p, n_loci)
    }

    fn slide(&mut self, chromosome: String, position: u64, reference_allele: char, counts: Vec<u64>, debug: bool) -> &mut Self {
        let (n, p, n_loci) = self.dim();
        
        for i in 1..n_loci {
            self.chromosomes[i-1] = self.chromosomes[i].clone(); // we need to clone the chromosome name
            self.positions[i-1] = self.positions[i];
            self.reference_alleles[i-1] = self.reference_alleles[i];
        }
        self.chromosomes[n_loci-1] = chromosome;
        self.positions[n_loci-1] = position;
        self.reference_alleles[n_loci-1] = reference_allele;

        // let mut count_matrix: Array2<f64> = Array2::<f64>::zeros((n, p));
        // let count_matrix_ref = &mut count_matrix;
        let X = self.count_matrix.clone();
        for locus_counter in 1..n_loci {
            for i in 0..NUMBER_OF_ALLELES {
                for j in 0..p {
                    self.count_matrix
                    .slice_mut(s![(NUMBER_OF_ALLELES*(locus_counter-1))+i, j])
                    .fill(X[[(NUMBER_OF_ALLELES*(locus_counter-0))+i, j]] as f64);
                }
            }
        }

        for i in 0..NUMBER_OF_ALLELES {
            for j in 0..p {
                self.count_matrix
                .slice_mut(s![(NUMBER_OF_ALLELES*(n_loci-1))+i, j])
                .fill(counts[(j*NUMBER_OF_ALLELES)+i] as f64);
            }
        }

        if debug {
            println!("############################");
            println!("n={}; p={}", n, p);
            println!("chr={:?}", self.chromosomes);
            println!("pos={:?}", self.positions);
            println!("ref={:?}", self.reference_alleles);
            println!("ORIG={:?}", X.slice(s![n-15..n,0..5]));
            println!("NEW={:?}", self.count_matrix.slice(s![n-15..n,0..5]));
        }

        self
    }

    fn find_coordinates(&self) -> (Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>) {
        let (n, p, n_loci) = self.dim();

        let mut missing_loci = vec![];
        let mut missing_pools = vec![];
        let mut present_loci = vec![];
        let mut present_pools = vec![];
        // find indices of missing loci and pools
        for j in 0..p {
            // iterating across loci because we have an inner loop to represent the 7 alleles per locus structure of the count matrix
            for i in 0..n_loci {
                let mut sum: f64 = 0.0;
                for k in 0..NUMBER_OF_ALLELES {
                    sum += self.count_matrix[[(i*NUMBER_OF_ALLELES)+k,j]];
                }
                if sum == 0.0 {
                    missing_loci.push(i);
                    missing_pools.push(j);
                }
            }
        }
        // find indices of missing loci
        for i in 0..n_loci {
            let mut bool_not_missing: bool = false;
            for j in missing_loci.iter() {
                if i == *j {
                    bool_not_missing = false;
                    break;
                } else {
                    bool_not_missing = true;
                }
            }
            if bool_not_missing {
                present_loci.push(i);
            }
        }
        // find indices of missing pools
        for i in 0..p {
            let mut bool_not_missing: bool = false;
            for j in missing_pools.iter() {
                if i == *j {
                    bool_not_missing = false;
                    break;
                } else {
                    bool_not_missing = true;
                }
            }
            if bool_not_missing {
                present_pools.push(i);
            }
        }

        (missing_loci, missing_pools, present_loci, present_pools)
    }

    // fn find_response_and_explanatory_matrices(&self) -> i32 {
    //     let (n, p, n_loci) = self.dim();
    //     let (missing_loci, missing_pools, present_loci, present_pools) = self.find_coordinates();
    //     let sorted_unique_pool_iterator = missing_pools.clone();
    //     sorted_unique_pool_iterator.sort();
    //     sorted_unique_pool_iterator.dedup();
    //     for i in sorted_unique_pool_iterator.iter() {
    //         let y: Array2<f64>;
    //         for j in missing_pools.iter() {
    //             if i == j {
    //                 y = self.count_matrix.slice(s![.., i]);
    //             }
    //         }
    //     }


    //     0
    // }

}

pub fn iterate_across_lines(stream:Result<Lines<BufReader<File>>>, n: usize){
    let min_base_quality: u8 = 20;
    let window_size: usize = 3;
    let mut locus_counter: u64 = 0;
    let mut vec: Vec<String>;
    let mut chromosomes: Vec<String> = Vec::new(); let chromosomes_ref = &mut chromosomes;
    let mut positions: Vec<u64> = Vec::new(); let positions_ref = &mut positions;
    let mut reference_alleles: Vec<char> = Vec::new(); let reference_alleles_ref = &mut reference_alleles;
    let mut count_matrix = Array2::<f64>::zeros((window_size*NUMBER_OF_ALLELES as usize, (n-3)/3 as usize)); let count_matrix_ref = &mut count_matrix;
    for line in stream.unwrap() {
        vec = line
                .unwrap()
                .split_whitespace()
                .map(str::to_string)
                .collect();
        let (chromosome, position, reference_allele, counts) = vec.parse_vec_string(min_base_quality);
        
        // if locus_counter == 1 {
        //     println!("{:?}", counts);
        // }

        // data_ref.extend(counts);
        // println!("{:?}", data_ref);
        if locus_counter < window_size as u64 {
            chromosomes_ref.extend([chromosome]);
            positions_ref.extend([position]);
            reference_alleles_ref.extend([reference_allele]);
            for i in 0..NUMBER_OF_ALLELES {
                for j in 0..(n-3)/3 {
                    count_matrix_ref.slice_mut(s![(NUMBER_OF_ALLELES*(locus_counter as usize))+i as usize, j]).fill(counts[(j*NUMBER_OF_ALLELES)+i] as f64);
                    // let x = count_matrix_ref[[locus_counter as usize, i]];
                    // println!("{}", x);
                }
            }
        } else {
            // encapsulate, impute and slide
            // encapsulate
            let mut window = Window {
                chromosomes: chromosomes_ref,
                positions: positions_ref,
                reference_alleles: reference_alleles_ref,
                count_matrix: count_matrix_ref // dimensions: p loci x n pools
            };

            // impute
            let (missing_loci, missing_pools, present_loci, present_pools) = window.find_coordinates();
            println!("missing_loci: {:?}", missing_loci);
            println!("missing_pools: {:?}", missing_pools);
            println!("present_loci: {:?}", present_loci);
            println!("present_pools: {:?}", present_pools);

            // let (response_idx, explanatory_idx) = window.find_response_and_explanatory_matrices();
            // println!("Response columns: {:?}", response_idx);
            // println!("Explanatory columns: {:?}", explanatory_idx);

            // slide
            window.slide(chromosome, position, reference_allele, counts, false);
 
        }
        // update counter
        locus_counter += 1;

    }

}
