use crate::{Error, Result};
use mpi::{
    datatype::{UncommittedUserDatatype, UserDatatype},
    ffi,
    raw::traits::FromRaw,
    traits::UncommittedDatatype,
    Count,
};
use std::mem::MaybeUninit;
use std::os::raw::c_int;

pub fn create_subarray<D>(
    sizes: &[Count],
    subsizes: &[Count],
    starts: &[Count],
    oldtype: &D,
) -> Result<UserDatatype>
where
    D: UncommittedDatatype,
{
    if sizes.len() != subsizes.len() || sizes.len() != starts.len() {
        return Err(Error::invalid_parameters(
            "Sizes, subsizes, and starts must have the same length",
        ));
    }
    unsafe {
        let mut newtype = MaybeUninit::uninit();
        ffi::MPI_Type_create_subarray(
            sizes.len() as Count,
            sizes.as_ptr(),
            subsizes.as_ptr(),
            starts.as_ptr(),
            ffi::MPI_ORDER_C as c_int,
            oldtype.as_raw(),
            newtype.as_mut_ptr(),
        );
        Ok(UncommittedUserDatatype::from_raw(newtype.assume_init()).commit())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mpi::datatype::UserDatatype;
    use mpi::traits::*;

    #[test]
    fn test_create_subarray() {
        let _universe = mpi::initialize().unwrap();
        let sizes = [4, 4];
        let subsizes = [2, 2];
        let starts = [1, 1];
        let oldtype = UserDatatype::contiguous(4, &f64::equivalent_datatype());
        let _newtype = create_subarray(&sizes, &subsizes, &starts, &oldtype).unwrap();
    }
}
