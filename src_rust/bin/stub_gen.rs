fn main() -> pyo3_stub_gen::Result<()> {
    let mut info = cave_closure::stub_info()?;
    let mut args = std::env::args_os().skip(1);
    if let Some(arg) = args.next() {
        if arg != "--output-dir" {
            return Err(std::io::Error::other("usage: stub_gen [--output-dir PATH]").into());
        }
        info.python_root = args
            .next()
            .ok_or_else(|| std::io::Error::other("--output-dir requires a path"))?
            .into();
        if args.next().is_some() {
            return Err(std::io::Error::other("usage: stub_gen [--output-dir PATH]").into());
        }
    }
    info.generate()?;
    Ok(())
}
