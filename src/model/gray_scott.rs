use super::{Domain, Parameters};
use mpi::{datatype::UserDatatype, topology::CartesianCommunicator, traits::Equivalence, Count};
use ndarray::Array2;

pub struct GrayScott {
    pub cart_comm: CartesianCommunicator,
    pub domain: Domain,
    pub params: Parameters,
    u: Array2<f64>,
    v: Array2<f64>,
    u_next: Array2<f64>,
    v_next: Array2<f64>,
    halo_type: UserDatatype,
}

impl GrayScott {
    pub fn new(cart_comm: CartesianCommunicator, domain: Domain, params: Parameters) -> Self {
        let shape = domain.local_shape_with_halo();
        let halo_type = Self::create_halo_type(&domain);
        Self {
            cart_comm,
            domain,
            params,
            u: Array2::ones(shape),
            v: Array2::zeros(shape),
            u_next: Array2::zeros(shape),
            v_next: Array2::zeros(shape),
            halo_type,
        }
    }

    fn create_halo_type(domain: &Domain) -> UserDatatype {
        UserDatatype::vector(
            domain.ny as Count,
            1,
            domain.nx_with_halo() as Count,
            &f64::equivalent_datatype(),
        )
    }
}
