use super::{Domain, Parameters};
use mpi::{
    datatype::{MutView, UserDatatype, View},
    // request::{LocalScope, RequestCollection},
    topology::CartesianCommunicator,
    traits::{Communicator, Destination, Equivalence, Source},
    Count,
};
use ndarray::{/*prelude::*,*/ Array2};
use std::{cell::UnsafeCell, ops::Deref};

struct Compass<T> {
    north: T,
    south: T,
    east: T,
    west: T,
}

struct Neighbors(Compass<Option<i32>>);

impl Neighbors {
    fn from_cart_comm(cart_comm: &CartesianCommunicator) -> Self {
        let (north, south) = cart_comm.shift(0, 1);
        let (west, east) = cart_comm.shift(1, 1);
        Self(Compass {
            north,
            south,
            east,
            west,
        })
    }
}

impl Deref for Neighbors {
    type Target = Compass<Option<i32>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Clone, Copy)]
enum Tag {
    S2N = 0,
    N2S = 1,
    E2W = 2,
    W2E = 3,
}

impl Tag {
    fn with_offset(self, offset: i32) -> i32 {
        offset + self as i32
    }
}

struct HaloDatatype {
    vertical: UserDatatype,
    horizontal: UserDatatype,
}

impl HaloDatatype {
    fn new(domain: &Domain) -> Self {
        let vertical = UserDatatype::vector(
            domain.ny as Count,
            1,
            domain.nx_with_halo() as Count,
            &f64::equivalent_datatype(),
        );
        let horizontal = UserDatatype::contiguous(domain.nx as Count, &f64::equivalent_datatype());
        Self {
            vertical,
            horizontal,
        }
    }
}

pub struct GrayScott {
    pub cart_comm: CartesianCommunicator,
    pub domain: Domain,
    pub params: Parameters,
    u: UnsafeCell<Array2<f64>>,
    v: UnsafeCell<Array2<f64>>,
    u_next: UnsafeCell<Array2<f64>>,
    v_next: UnsafeCell<Array2<f64>>,
    halo_type: HaloDatatype,
    neighbors: Neighbors,
}

impl GrayScott {
    pub fn new(cart_comm: CartesianCommunicator, domain: Domain, params: Parameters) -> Self {
        let shape = domain.local_shape_with_halo();
        Self {
            params,
            u: UnsafeCell::new(Array2::ones(shape)),
            v: UnsafeCell::new(Array2::zeros(shape)),
            u_next: UnsafeCell::new(Array2::zeros(shape)),
            v_next: UnsafeCell::new(Array2::zeros(shape)),
            halo_type: HaloDatatype::new(&domain),
            neighbors: Neighbors::from_cart_comm(&cart_comm),
            cart_comm,
            domain,
        }
    }

    pub fn u(&self) -> &Array2<f64> {
        unsafe { &*self.u.get() }
    }
    pub fn v(&self) -> &Array2<f64> {
        unsafe { &*self.v.get() }
    }

    fn create_recv_views<'a, 'b>(
        &'a self,
        tile: &'b UnsafeCell<Array2<f64>>,
    ) -> Compass<MutView<'a, 'b, UserDatatype, [f64]>> {
        let tile = unsafe { &mut *tile.get() };
        let ptr = tile.as_mut_ptr();
        let HaloDatatype {
            vertical: vtype,
            horizontal: htype,
        } = &self.halo_type;
        let domain = &self.domain;
        Compass {
            north: unsafe {
                MutView::with_count_and_datatype(
                    std::slice::from_raw_parts_mut(ptr.add(domain.disp(0, 1)), domain.nx),
                    1,
                    htype,
                )
            },
            south: unsafe {
                MutView::with_count_and_datatype(
                    std::slice::from_raw_parts_mut(
                        ptr.add(domain.disp(domain.ny_with_halo() - 1, 1)),
                        domain.nx,
                    ),
                    1,
                    htype,
                )
            },
            east: unsafe {
                let disp = domain.disp(1, domain.nx_with_halo() - 1);
                MutView::with_count_and_datatype(
                    std::slice::from_raw_parts_mut(ptr.add(disp), tile.len() - disp),
                    1,
                    vtype,
                )
            },
            west: unsafe {
                let disp = domain.disp(1, 0);
                MutView::with_count_and_datatype(
                    std::slice::from_raw_parts_mut(ptr.add(disp), tile.len() - disp),
                    1,
                    vtype,
                )
            },
        }
    }

    fn create_send_views<'a, 'b>(
        &'a self,
        tile: &'b UnsafeCell<Array2<f64>>,
    ) -> Compass<View<'a, 'b, UserDatatype, [f64]>> {
        let domain = &self.domain;
        let HaloDatatype {
            vertical: vtype,
            horizontal: htype,
        } = &self.halo_type;
        let tile = unsafe { &*tile.get() };
        Compass {
            north: unsafe {
                View::with_count_and_datatype(
                    &tile.as_slice_memory_order().unwrap()[domain.disp(1, 1)..],
                    1,
                    htype,
                )
            },
            south: unsafe {
                View::with_count_and_datatype(
                    &tile.as_slice_memory_order().unwrap()
                        [domain.disp(domain.ny_with_halo() - 2, 1)..],
                    1,
                    htype,
                )
            },
            east: unsafe {
                View::with_count_and_datatype(
                    &tile.as_slice_memory_order().unwrap()
                        [domain.disp(1, domain.nx_with_halo() - 2)..],
                    1,
                    vtype,
                )
            },
            west: unsafe {
                View::with_count_and_datatype(
                    &tile.as_slice_memory_order().unwrap()[domain.disp(1, 1)..],
                    1,
                    vtype,
                )
            },
        }
    }

    pub fn exchange_halos(&mut self) {
        let mut u_recv_views = self.create_recv_views(&self.u);
        let u_send_views = self.create_send_views(&self.u);
        let mut v_recv_views = self.create_recv_views(&self.v);
        let v_send_views = self.create_send_views(&self.v);

        mpi::request::multiple_scope(8, |scope, coll| {
            for (recv_views, ofs) in [(&mut u_recv_views, 0), (&mut v_recv_views, 100)] {
                for (neighbor, view, tag) in [
                    (self.neighbors.north, &mut recv_views.north, Tag::N2S),
                    (self.neighbors.south, &mut recv_views.south, Tag::S2N),
                    (self.neighbors.east, &mut recv_views.east, Tag::E2W),
                    (self.neighbors.west, &mut recv_views.west, Tag::W2E),
                ] {
                    if let Some(rank) = neighbor {
                        coll.add(
                            self.cart_comm
                                .process_at_rank(rank)
                                .immediate_receive_into_with_tag(scope, view, tag.with_offset(ofs)),
                        );
                    }
                }
            }

            for (send_views, ofs) in [(&u_send_views, 0), (&v_send_views, 100)] {
                for (neighbor, view, tag) in [
                    (self.neighbors.north, &send_views.north, Tag::S2N),
                    (self.neighbors.south, &send_views.south, Tag::N2S),
                    (self.neighbors.east, &send_views.east, Tag::W2E),
                    (self.neighbors.west, &send_views.west, Tag::E2W),
                ] {
                    if let Some(rank) = neighbor {
                        self.cart_comm
                            .process_at_rank(rank)
                            .send_with_tag(view, tag.with_offset(ofs));
                    }
                }
            }

            coll.wait_all(&mut vec![]);
        });
    }

    pub fn compute_next_state(&mut self) {
        let (f, k, dt, du, dv) = (
            self.params.f,
            self.params.k,
            self.params.dt,
            self.params.du,
            self.params.dv,
        );
        let u = unsafe { &*self.u.get() };
        let v = unsafe { &*self.v.get() };
        let u_next = unsafe { &mut *self.u_next.get() };
        let v_next = unsafe { &mut *self.v_next.get() };

        let (ny, nx) = (self.domain.ny, self.domain.nx);

        for y in 1..ny + 1 {
            for x in 1..nx + 1 {
                let cur_u = u[[y, x]];
                let cur_v = v[[y, x]];
                let laplacian_u =
                    u[[y - 1, x]] + u[[y, x - 1]] + u[[y, x + 1]] + u[[y + 1, x]] - 4.0 * cur_u;
                let laplacian_v =
                    v[[y - 1, x]] + v[[y, x - 1]] + v[[y, x + 1]] + v[[y + 1, x]] - 4.0 * cur_v;
                let diffusion_u = du * laplacian_u;
                let diffusion_v = dv * laplacian_v;
                let uv2 = cur_u * cur_v * cur_v;
                let react_u = -uv2 + f * (1.0 - cur_u);
                let react_v = uv2 - (f + k) * cur_v;
                u_next[[y, x]] = cur_u + dt * (diffusion_u + react_u);
                v_next[[y, x]] = cur_v + dt * (diffusion_v + react_v);
            }
        }
    }

    pub fn swap_buffers(&mut self) {
        std::mem::swap(&mut self.u, &mut self.u_next);
        std::mem::swap(&mut self.v, &mut self.v_next);
    }

    pub fn step(&mut self) {
        self.exchange_halos();
        self.compute_next_state();
        self.swap_buffers();
    }
}
