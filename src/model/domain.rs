use mpi::{topology::CartesianCommunicator, traits::Communicator};

#[derive(Debug)]
pub struct Domain {
    pub total_nx: usize,
    pub total_ny: usize,
    pub nx: usize,
    pub ny: usize,
    pub start_x: usize,
    pub start_y: usize,
}

impl Domain {
    pub fn new(dims: [usize; 2], coords: [usize; 2], tile_size: [usize; 2]) -> Self {
        Self {
            total_nx: dims[1] * tile_size[1],
            total_ny: dims[0] * tile_size[0],
            nx: tile_size[1],
            ny: tile_size[0],
            start_x: coords[1] * tile_size[1],
            start_y: coords[0] * tile_size[0],
        }
    }

    pub fn from_cart_comm(comm: &CartesianCommunicator, tile_size: [usize; 2]) -> Self {
        Self::new(
            comm.get_layout()
                .dims
                .try_into()
                .map(|arr: [i32; 2]| [arr[0] as usize, arr[1] as usize])
                .unwrap(),
            comm.rank_to_coordinates(comm.rank())
                .try_into()
                .map(|arr: [i32; 2]| [arr[0] as usize, arr[1] as usize])
                .unwrap(),
            tile_size,
        )
    }

    pub fn nx_with_halo(&self) -> usize {
        self.nx + 2
    }
    pub fn ny_with_halo(&self) -> usize {
        self.ny + 2
    }
    pub fn size_with_halo(&self) -> usize {
        self.nx_with_halo() * self.ny_with_halo()
    }
    pub fn size(&self) -> usize {
        self.nx * self.ny
    }
    pub fn total_size(&self) -> usize {
        self.total_nx * self.total_ny
    }
}
