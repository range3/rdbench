use super::{Domain, Parameters};
use mpi::topology::CartesianCommunicator;

pub struct GrayScott {
    pub cart_comm: CartesianCommunicator,
    pub domain: Domain,
    pub params: Parameters,
}

impl GrayScott {
    // pub fn new(params: Parameters, sz_tile_x: usize, sz_tile_y: usize) -> Self {}
}
