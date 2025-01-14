mod cli;
use cli::Cli;

fn main() {
    let cli = Cli::from_args();
    dbg!(&cli);
    if let Err(e) = cli.validate_parameters() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
