use crate::error::{Error, Result};

#[derive(Debug, Clone, Copy)]
pub struct Parameters {
    pub du: f64, // diffusion rate of u (activator)
    pub dv: f64, // diffusion rate of v (inhibitor)
    pub f: f64,  // feed rate
    pub k: f64,  // kill rate
    pub dt: f64, // time step
                 // Note: dx = dy = 1.0 (grid spacing is fixed)
}

impl Parameters {
    pub fn new(du: f64, dv: f64, f: f64, k: f64, dt: f64) -> Self {
        Self { du, dv, f, k, dt }
    }

    pub fn default() -> Self {
        Self {
            du: 0.1,
            dv: 0.05,
            f: 0.04,
            k: 0.06075,
            dt: 0.2,
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.du <= 0.0 || self.dv <= 0.0 {
            return Err(Error::InvalidParameters(
                "Diffusion rates must be > 0".to_string(),
            ));
        }
        if self.f < 0.0 {
            return Err(Error::InvalidParameters(
                "Feed rate must be >= 0".to_string(),
            ));
        }
        if self.k < 0.0 {
            return Err(Error::InvalidParameters(
                "Kill rate must be >= 0".to_string(),
            ));
        }
        if self.dt <= 0.0 {
            return Err(Error::InvalidParameters(
                "Time step must be > 0".to_string(),
            ));
        }

        let max_stable_dt = 0.5 / self.du.max(self.dv);
        if self.dt > max_stable_dt {
            return Err(Error::InvalidParameters(format!(
                "Time step exceeds stability limit (max: {})",
                max_stable_dt
            )));
        }

        Ok(())
    }

    pub fn warn(&self) {
        if self.f > 0.1 || self.k > 0.1 {
            eprintln!(
                "Warning: F and k values above 0.1 might not produce typical Turing patterns"
            );
        }

        if self.du <= self.dv {
            eprintln!("Warning: Du should typically be less than Dv for pattern formation");
        }
    }
}
