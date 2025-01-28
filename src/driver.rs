use mpi::traits::Communicator;

use crate::{
    args::{Args, IoField},
    io::{self, IOStrategy},
    model::{FieldType, GrayScott},
    Result,
};

pub struct Driver<'a> {
    model: GrayScott,
    io_strategy: Box<dyn IOStrategy>,
    args: &'a Args,
}

impl<'a> Driver<'a> {
    pub fn new(model: GrayScott, args: &'a Args) -> Result<Self> {
        let io_strategy = io::create_io_strategy(args, &model.domain)?;
        Ok(Self {
            model,
            io_strategy,
            args,
        })
    }

    pub fn run(&mut self) -> Result<()> {
        let mut ckpt_idx = 0;
        let steps = self.args.steps;
        let interval = self.args.interval;

        if self.args.init_output {
            self.ckpt(ckpt_idx, 0)?;
            ckpt_idx += 1;
        }

        for step in 1..=steps {
            self.model.step();

            if interval != 0 && step % interval == 0 {
                self.ckpt(ckpt_idx, step)?;
                ckpt_idx += 1;
            }
        }

        Ok(())
    }

    fn ckpt(&self, idx: usize, step: usize) -> Result<()> {
        let comm = &self.model.cart_comm;

        if self.args.io_field.should_io_u() {
            self.model
                .checkpoint(&*self.io_strategy, idx, FieldType::U)?;
        }
        if self.args.io_field.should_io_v() {
            self.model
                .checkpoint(&*self.io_strategy, idx, FieldType::V)?;
        }

        if self.args.verbose && comm.rank() == 0 {
            if self.args.io_field != IoField::None {
                println!("Step: {}, Checkpoint: {}", step, idx);
            }
        }

        Ok(())
    }
}
