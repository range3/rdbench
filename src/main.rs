use rdbench::{args, driver::Driver, model::GrayScottFactory, Result};

fn main() -> Result<()> {
    let universe = mpi::initialize().unwrap();
    let args = args::parse();
    let model = GrayScottFactory::create(&universe, &args)?;
    let mut driver = Driver::new(model, &args)?;
    driver.run()?;
    Ok(())
}
