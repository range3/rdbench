#[derive(Debug)]
pub enum Error {
    InvalidParameters(String),
    InvalidDomain {
        np: i32,
        dims: Vec<i32>,
        msg: String,
    },
    MpiError {
        code: i32,
        msg: String,
    },
}

impl Error {
    pub fn invalid_parameters(msg: &str) -> Self {
        Error::InvalidParameters(msg.to_string())
    }

    pub fn invalid_domain(np: i32, dims: Vec<i32>, msg: &str) -> Self {
        Error::InvalidDomain {
            np,
            dims,
            msg: msg.to_string(),
        }
    }

    pub fn mpi_error(code: i32, msg: &str) -> Self {
        Error::MpiError {
            code,
            msg: msg.to_string(),
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidParameters(msg) => write!(f, "Invalid parameters: {}", msg),
            Error::InvalidDomain { np, dims, msg } => {
                write!(f, "Invalid domain (np={}, dims={:?}): {}", np, dims, msg)
            }
            Error::MpiError { code, msg } => write!(f, "MPI error: ({}) {}", code, msg),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::InvalidParameters(_) => None,
            Error::InvalidDomain { .. } => None,
            Error::MpiError { .. } => None,
        }
    }
}
