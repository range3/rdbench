use mpi::traits::*;
use rdbench::args;
use rdbench::model;

fn main() {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let args = args::parse();
    let model = model::GrayScottFactory::create(&universe, &args).unwrap();

    // if world.rank() == 0 {
    //     dbg!(&model.params);
    //     dbg!(&args);
    // }

    for i in 0..world.size() {
        world.barrier();
        if i == world.rank() {
            dbg!(&model.domain);
        }
    }
}
