#[macro_use(s)]
extern crate ndarray;

mod functions;

// fn main() {
//     println!("Building PoPoolImpute on rust");
// 	functions::closure_test(|x| println!("{}", x));   
// }

fn main() {
    let stream = functions::stream_file("res/test.pileup");
    let n: usize = functions::count_columns(stream.unwrap());

    let stream = functions::stream_file("res/test.pileup");
    functions::iterate_across_lines(stream, n);

}
