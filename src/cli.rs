use clap::Parser;
use clap::ValueEnum;

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Benchmark for MPI-IO using Gray-Scott reaction-diffusion system"
)]
pub struct Cli {
    // Tile partitioning settings
    #[arg(long, default_value = "0")]
    pub nr_tiles_x: usize,
    #[arg(long, default_value = "0")]
    pub nr_tiles_y: usize,
    #[arg(long, default_value = "1024")]
    pub sz_tile_x: usize,
    #[arg(long, default_value = "1024")]
    pub sz_tile_y: usize,

    // Simulation settings
    #[arg(long, default_value = "20000")]
    pub steps: usize,
    #[arg(long, default_value = "200")]
    pub interval: usize,
    #[arg(long)]
    pub init_output: bool,
    #[arg(long)]
    pub verbose: bool,

    // I/O settings
    #[arg(long, default_value = "out/o_")]
    pub output: String,
    #[arg(long, value_enum, default_value = "canonical")]
    pub file_layout: FileLayout,
    #[arg(long, value_enum, default_value = "v")]
    pub io_field: IoField,
    #[arg(long)]
    pub collective: bool,
    #[arg(long)]
    pub nosync: bool,
    #[arg(long)]
    pub validate: bool,
    #[arg(long)]
    pub prettify: bool,

    // Gray-Scott model parameters
    #[arg(long, default_value = "0.04")]
    pub param_f: f64,
    #[arg(long, default_value = "0.06075")]
    pub param_k: f64,
    #[arg(long, default_value = "0.2")]
    pub param_dt: f64,
    #[arg(long, default_value = "0.1")]
    pub param_du: f64,
    #[arg(long, default_value = "0.05")]
    pub param_dv: f64,
}

impl Cli {
    pub fn from_args() -> Self {
        Self::parse()
    }

    pub fn validate_parameters(&self) -> Result<(), String> {
        if self.param_f < 0.0
            || self.param_k < 0.0
            || self.param_dt <= 0.0
            || self.param_du <= 0.0
            || self.param_dv <= 0.0
        {
            return Err("Model parameters must be positive".to_string());
        }
        Ok(())
    }

    #[allow(dead_code)]
    pub fn nr_files(&self) -> usize {
        if self.interval == 0 {
            return 0;
        }
        let output_count = match self.io_field {
            IoField::U | IoField::V => 1,
            IoField::Both => 2,
            IoField::None => 0,
        };

        (self.steps / self.interval + if self.init_output { 1 } else { 0 }) * output_count
    }
}

#[derive(Parser, Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum FileLayout {
    Canonical,
    Log,
}

#[derive(Parser, Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum IoField {
    U,
    V,
    Both,
    None,
}

impl IoField {
    #[allow(dead_code)]
    pub fn should_io_u(&self) -> bool {
        matches!(self, IoField::U | IoField::Both)
    }

    #[allow(dead_code)]
    pub fn should_io_v(&self) -> bool {
        matches!(self, IoField::V | IoField::Both)
    }
}
