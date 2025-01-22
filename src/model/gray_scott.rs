use super::{Domain, Parameters};
use mpi::topology::CartesianCommunicator;

pub struct GrayScott {
    pub cart_comm: CartesianCommunicator,
    pub domain: Domain,
    pub params: Parameters,
}

impl GrayScott {
    pub fn new(cart_comm: CartesianCommunicator, domain: Domain, params: Parameters) -> Self {
        Self {
            cart_comm,
            domain,
            params,
        }
    }
}
