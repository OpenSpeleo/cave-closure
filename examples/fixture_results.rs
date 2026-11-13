//! Compute the persisted fixture corpus with every implemented linear solver.
//! Independent Python references, rather than these outputs, define the optimum.
use cave_closure::{
    LinearBackend, SolveConfig, SolveControl, SolveReport, Workspace, json, prepare,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, error::Error, fs, path::Path};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Corpus {
    schema_version: u32,
    cases: Vec<Case>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    name: String,
    input: String,
    reference: String,
    results: BTreeMap<String, String>,
}

#[derive(Serialize)]
struct SolverResult {
    schema_version: u32,
    solver: &'static str,
    reduce: bool,
    report: SolveReport,
}

fn member(directory: &Path, filename: &str) -> Result<std::path::PathBuf, Box<dyn Error>> {
    if Path::new(filename)
        .file_name()
        .and_then(|name| name.to_str())
        != Some(filename)
    {
        return Err(format!("expected a corpus filename, got {filename:?}").into());
    }
    Ok(directory.join(filename))
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err("usage: fixture_results FIXTURE_DIRECTORY OUTPUT_DIRECTORY".into());
    }
    let fixtures = Path::new(&args[0]);
    let output = Path::new(&args[1]);
    let corpus: Corpus = serde_json::from_str(&fs::read_to_string(fixtures.join("corpus.json"))?)?;
    if corpus.schema_version != 1 || corpus.cases.is_empty() {
        return Err("unsupported or empty fixture corpus".into());
    }
    fs::create_dir_all(output)?;
    let solvers = [
        ("reduced_pcg", LinearBackend::ReducedPcg, true),
        ("unreduced_pcg", LinearBackend::ReducedPcg, false),
        ("sparse_qr", LinearBackend::SparseQr, false),
        ("sparse_cholesky", LinearBackend::SparseCholesky, false),
    ];
    for case in corpus.cases {
        // Validate every manifest path, including the independently generated reference.
        member(fixtures, &case.reference)?;
        let original = json::from_json(&fs::read_to_string(member(fixtures, &case.input)?)?)?;
        if case.results.len() != solvers.len() {
            return Err(format!("{}: incomplete solver manifest", case.name).into());
        }
        for (solver, backend, reduce) in solvers {
            let mut problem = original.clone();
            problem.algorithm.backend = backend;
            let started = std::time::Instant::now();
            let report = prepare(&problem, &Default::default())?.solve(
                &SolveConfig {
                    reduce,
                    ..Default::default()
                },
                &mut Workspace::default(),
                &SolveControl::default(),
            )?;
            report
                .require_converged()
                .map_err(|error| format!("{} / {solver}: {error}", case.name))?;
            // The core's serializer explicitly rejects non-finite diagnostics/results.
            report.to_json()?;
            let result = SolverResult {
                schema_version: 1,
                solver,
                reduce,
                report,
            };
            fs::write(
                member(output, &case.results[solver])?,
                serde_json::to_string(&result)?,
            )?;
            println!(
                "{} / {solver}: converged ({:.3}s)",
                case.name,
                started.elapsed().as_secs_f64()
            );
        }
    }
    Ok(())
}
