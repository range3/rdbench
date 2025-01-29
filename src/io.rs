use crate::{
    args::{Args, FileLayout},
    model::{Domain, FieldType},
    mpi::{
        datatype,
        io::{File, Mode},
    },
    Error, Result,
};
use mpi::{
    datatype::{UserDatatype, View},
    traits::{Communicator, Equivalence},
};
use ndarray::Array2;

fn create_tile_type(domain: &Domain) -> Result<UserDatatype> {
    let sizes = [domain.ny_with_halo() as i32, domain.nx_with_halo() as i32];
    let subsizes = [domain.ny as i32, domain.nx as i32];
    let starts = [1, 1];
    datatype::create_subarray(&sizes, &subsizes, &starts, &f64::equivalent_datatype())
}

pub enum IOStrategyEnum {
    Canonical(CanonicalIO),
    Log(LogIO),
}

pub trait IOStrategy {
    fn write<C: Communicator>(
        &self,
        comm: &C,
        data: &Array2<f64>,
        domain: &Domain,
        field_type: FieldType,
        idx: usize,
    ) -> Result<()>;
}

impl IOStrategy for IOStrategyEnum {
    fn write<C: Communicator>(
        &self,
        comm: &C,
        data: &Array2<f64>,
        domain: &Domain,
        field_type: FieldType,
        idx: usize,
    ) -> Result<()> {
        match self {
            Self::Canonical(strategy) => strategy.write(comm, data, domain, field_type, idx),
            Self::Log(strategy) => strategy.write(comm, data, domain, field_type, idx),
        }
    }
}

pub trait FileIOStrategy: IOStrategy {
    fn file_name(&self, domain: &Domain, field_type: FieldType, idx: usize) -> String {
        let type_str = match field_type {
            FieldType::U => "u",
            FieldType::V => "v",
        };
        format!(
            "{}-{}x{}-{}x{}-{:06}.bin",
            type_str, domain.total_nx, domain.total_ny, domain.nx, domain.ny, idx
        )
    }
}

pub struct CanonicalIO {
    output_prefix: String,
    collective: bool,
    sync: bool,
    file_type: UserDatatype,
    tile_type: UserDatatype,
}

impl CanonicalIO {
    pub fn new(
        output_prefix: String,
        collective: bool,
        sync: bool,
        domain: &Domain,
    ) -> Result<Self> {
        let sizes = [domain.total_ny as i32, domain.total_nx as i32];
        let subsizes = [domain.ny as i32, domain.nx as i32];
        let starts = [domain.start_y as i32, domain.start_x as i32];
        let file_type =
            datatype::create_subarray(&sizes, &subsizes, &starts, &f64::equivalent_datatype())?;

        Ok(Self {
            output_prefix,
            collective,
            sync,
            file_type,
            tile_type: create_tile_type(domain)?,
        })
    }
}

impl FileIOStrategy for CanonicalIO {}

impl IOStrategy for CanonicalIO {
    fn write<C: Communicator>(
        &self,
        comm: &C,
        data: &Array2<f64>,
        domain: &Domain,
        field_type: FieldType,
        idx: usize,
    ) -> Result<()> {
        let path = format!(
            "{}{}",
            self.output_prefix,
            self.file_name(domain, field_type, idx)
        );
        let mode = Mode::CREATE | Mode::WRONLY | Mode::UNIQUE_OPEN;
        let file = File::open(comm, &path, mode)?;
        file.set_atomicity(false)?;
        file.set_view(0, &f64::equivalent_datatype(), &self.file_type, "native")?;

        let tile_view = unsafe {
            View::with_count_and_datatype(
                data.as_slice_memory_order()
                    .ok_or_else(|| Error::invalid_data("Failed to get slice from Array2"))?,
                1,
                &self.tile_type,
            )
        };

        if self.collective {
            file.write_at_all(0, &tile_view)?;
        } else {
            file.write_at(0, &tile_view)?;
        }

        if self.sync {
            file.sync()?;
        }

        Ok(())
    }
}

pub struct LogIO {
    output_prefix: String,
    collective: bool,
    sync: bool,
    tile_type: UserDatatype,
}

impl LogIO {
    pub fn new(
        output_prefix: String,
        collective: bool,
        sync: bool,
        domain: &Domain,
    ) -> Result<Self> {
        Ok(Self {
            output_prefix,
            collective,
            sync,
            tile_type: create_tile_type(domain)?,
        })
    }
}

impl FileIOStrategy for LogIO {}

impl IOStrategy for LogIO {
    fn write<C: Communicator>(
        &self,
        comm: &C,
        data: &Array2<f64>,
        domain: &Domain,
        field_type: FieldType,
        idx: usize,
    ) -> Result<()> {
        let path = format!(
            "{}{}",
            self.output_prefix,
            self.file_name(domain, field_type, idx)
        );
        let mode = Mode::CREATE | Mode::WRONLY | Mode::UNIQUE_OPEN;
        let file = File::open(comm, &path, mode)?;
        file.set_atomicity(false)?;

        let offset = (domain.size() * comm.rank() as usize * std::mem::size_of::<f64>()) as i64;
        let tile_view = unsafe {
            View::with_count_and_datatype(
                data.as_slice_memory_order()
                    .ok_or_else(|| Error::invalid_data("Failed to get slice from Array2"))?,
                1,
                &self.tile_type,
            )
        };

        if self.collective {
            file.write_at_all(offset, &tile_view)?;
        } else {
            file.write_at(offset, &tile_view)?;
        }

        if self.sync {
            file.sync()?;
        }

        Ok(())
    }
}

pub fn create_io_strategy(args: &Args, domain: &Domain) -> Result<IOStrategyEnum> {
    match args.file_layout {
        FileLayout::Canonical => {
            let strategy =
                CanonicalIO::new(args.output.clone(), args.collective, !args.nosync, domain)?;
            Ok(IOStrategyEnum::Canonical(strategy))
        }
        FileLayout::Log => {
            let strategy = LogIO::new(args.output.clone(), args.collective, !args.nosync, domain)?;
            Ok(IOStrategyEnum::Log(strategy))
        }
    }
}
