use super::{Domain, GrayScott, Parameters};
use crate::args::Args;
use crate::error::{Error, Result};
use mpi::environment::Universe;
use mpi::ffi;
use mpi::topology::CartesianCommunicator;
use mpi::traits::Communicator;

pub struct GrayScottFactory;

impl GrayScottFactory {
    pub fn create(universe: &Universe, args: &Args) -> Result<GrayScott> {
        let params = Parameters::new(
            args.param_du,
            args.param_dv,
            args.param_f,
            args.param_k,
            args.param_dt,
        );
        let cart_comm =
            Self::create_cart_comm(&universe.world(), [args.nr_tiles_y, args.nr_tiles_x])?;
        let domain = Domain::from_cart_comm(&cart_comm, [args.sz_tile_y, args.sz_tile_x]);
        Ok(GrayScott {
            cart_comm,
            domain,
            params,
        })
    }

    fn create_cart_comm(
        comm: &impl Communicator,
        dims: [usize; 2],
    ) -> Result<CartesianCommunicator> {
        let dims = Self::create_dims(comm.size(), &dims.map(|d| d as i32))?;
        comm.create_cartesian_communicator(&dims[..], &[true; 2], true)
            .ok_or(Error::invalid_domain(
                comm.size(),
                dims,
                "Failed to create Cartesian communicator",
            ))
    }

    fn create_dims(np: i32, dims: &[i32]) -> Result<Vec<i32>> {
        if np <= 0 {
            return Err(Error::InvalidParameters(
                "Number of processes must be positive".to_string(),
            ));
        }
        if dims.is_empty() {
            return Err(Error::InvalidParameters(
                "Dimensions must not be empty".to_string(),
            ));
        }
        if dims.iter().any(|&d| d < 0) {
            return Err(Error::InvalidParameters(
                "Dimensions must be positive".to_string(),
            ));
        }

        let mut dims = dims.to_vec();
        let ndims = dims.len() as i32;

        let result = unsafe { ffi::MPI_Dims_create(np, ndims, dims.as_mut_ptr()) };

        if result != ffi::MPI_SUCCESS as i32 {
            return Err(Error::mpi_error(result, "Failed to create dimensions"));
        }

        Ok(dims)
    }
}
