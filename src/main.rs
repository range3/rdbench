use mpi::traits::*;
use rdbench::args;
use rdbench::model;

fn main() {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let args = args::parse();
    let mut model = model::GrayScottFactory::create(&universe, &args).unwrap();

    for i in 0..world.size() {
        world.barrier();
        if i == world.rank() {
            dbg!(&model.u());
        }
    }
    model.exchange_halos();
    for i in 0..world.size() {
        world.barrier();
        if i == world.rank() {
            dbg!(&model.u());
        }
    }
}
