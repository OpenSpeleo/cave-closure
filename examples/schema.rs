fn main() {
    let schema = if std::env::args().nth(1).as_deref() == Some("result") {
        cave_closure::json::result_schema()
    } else {
        cave_closure::json::input_schema()
    };
    println!("{}", serde_json::to_string_pretty(&schema).unwrap());
}
