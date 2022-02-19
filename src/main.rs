mod functions;

// fn main() {
//     println!("Building PoPoolImpute on rust");
// 	functions::closure_test(|x| println!("{}", x));   
// }

fn main() {
    let stream = functions::stream_file("res/test.pileup");
    // functions::print_line_by_line(stream);
    // functions::print_column_by_column_per_line(stream);
    functions::parse_line(stream);

}
