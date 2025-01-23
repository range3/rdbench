use super::Domain;
use ndarray::Array2;

pub trait Filler {
    fn fill(&self, tile: &mut Array2<f64>, domain: Domain);
}

pub struct ConstantFiller {
    value: f64,
}

impl ConstantFiller {
    pub fn new(value: f64) -> Self {
        Self { value }
    }
}

impl Filler for ConstantFiller {
    fn fill(&self, tile: &mut Array2<f64>, _domain: Domain) {
        tile.fill(self.value);
    }
}

pub struct CenterBlockFiller {
    value: f64,
    block_size_x: usize,
    block_size_y: usize,
}

impl CenterBlockFiller {
    pub fn new(value: f64, block_size_x: usize, block_size_y: usize) -> Self {
        Self {
            value,
            block_size_x,
            block_size_y,
        }
    }
}

impl Filler for CenterBlockFiller {
    fn fill(&self, tile: &mut Array2<f64>, domain: Domain) {
        let block_size_x = std::cmp::min(self.block_size_x, domain.total_nx);
        let block_size_y = std::cmp::min(self.block_size_y, domain.total_ny);
        let start_x = domain.start_x;
        let start_y = domain.start_y;
        let block_start_x = domain.total_nx / 2 - block_size_x / 2;
        let block_start_y = domain.total_ny / 2 - block_size_y / 2;
        let block_end_x = block_start_x + block_size_x;
        let block_end_y = block_start_y + block_size_y;

        for y in 0..domain.ny {
            for x in 0..domain.nx {
                if block_start_x <= start_x + x
                    && start_x + x < block_end_x
                    && block_start_y <= start_y + y
                    && start_y + y < block_end_y
                {
                    tile[[y + 1, x + 1]] = self.value;
                }
            }
        }
    }
}
