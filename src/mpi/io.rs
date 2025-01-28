use mpi::{
    ffi,
    raw::traits::*,
    traits::{Buffer, BufferMut, Communicator, Datatype},
};
use std::{mem::MaybeUninit, os::raw::c_int};

type Offset = ffi::MPI_Offset;

#[derive(Debug, Clone, Copy)]
pub struct Mode(c_int);

impl Mode {
    pub const CREATE: Mode = Mode(ffi::MPI_MODE_CREATE as c_int);
    pub const RDONLY: Mode = Mode(ffi::MPI_MODE_RDONLY as c_int);
    pub const WRONLY: Mode = Mode(ffi::MPI_MODE_WRONLY as c_int);
    pub const RDWR: Mode = Mode(ffi::MPI_MODE_RDWR as c_int);
    pub const DELETE_ON_CLOSE: Mode = Mode(ffi::MPI_MODE_DELETE_ON_CLOSE as c_int);
    pub const UNIQUE_OPEN: Mode = Mode(ffi::MPI_MODE_UNIQUE_OPEN as c_int);
    pub const EXCL: Mode = Mode(ffi::MPI_MODE_EXCL as c_int);
    pub const APPEND: Mode = Mode(ffi::MPI_MODE_APPEND as c_int);
    pub const SEQUENTIAL: Mode = Mode(ffi::MPI_MODE_SEQUENTIAL as c_int);
}

unsafe impl AsRaw for Mode {
    type Raw = c_int;
    fn as_raw(&self) -> Self::Raw {
        self.0
    }
}

impl std::ops::BitOr for Mode {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Mode(self.0 | rhs.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mpi::ffi;

    #[test]
    fn test_mode() {
        let mode = Mode::RDONLY | Mode::WRONLY;
        assert_eq!(mode.0 as u32, ffi::MPI_MODE_RDONLY | ffi::MPI_MODE_WRONLY);
    }
}

pub struct File {
    handle: ffi::MPI_File,
}

unsafe impl AsRaw for File {
    type Raw = ffi::MPI_File;
    fn as_raw(&self) -> Self::Raw {
        self.handle
    }
}

impl FromRaw for File {
    unsafe fn from_raw(handle: Self::Raw) -> Self {
        Self { handle }
    }
}

unsafe impl MatchesRaw for File {}

impl Drop for File {
    fn drop(&mut self) {
        unsafe {
            ffi::MPI_File_close(&mut self.handle);
        }
    }
}

impl File {
    pub fn open<P: AsRef<std::path::Path>>(
        comm: &impl Communicator,
        filename: P,
        amode: Mode,
    ) -> crate::Result<Self> {
        let mut handle = MaybeUninit::uninit();
        let filename = std::ffi::CString::new(filename.as_ref().to_str().unwrap()).unwrap();
        unsafe {
            let ec = ffi::MPI_File_open(
                comm.as_raw(),
                filename.as_ptr(),
                amode.as_raw(),
                ffi::RSMPI_INFO_NULL,
                handle.as_mut_ptr(),
            );
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_open"));
            }
            Ok(Self {
                handle: handle.assume_init(),
            })
        }
    }

    pub fn write_at<B>(&self, offset: Offset, buf: &B) -> crate::Result<()>
    where
        B: Buffer,
    {
        unsafe {
            let ec = ffi::MPI_File_write_at(
                self.handle,
                offset,
                buf.pointer(),
                buf.count(),
                buf.as_datatype().as_raw(),
                ffi::RSMPI_STATUS_IGNORE,
            );
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_write_at"));
            }
            Ok(())
        }
    }

    pub fn read_at<B>(&self, offset: Offset, buf: &mut B) -> crate::Result<()>
    where
        B: BufferMut,
    {
        unsafe {
            let ec = ffi::MPI_File_read_at(
                self.handle,
                offset,
                buf.pointer_mut(),
                buf.count(),
                buf.as_datatype().as_raw(),
                ffi::RSMPI_STATUS_IGNORE,
            );
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_read_at"));
            }
            Ok(())
        }
    }

    pub fn write_at_all<B>(&self, offset: Offset, buf: &B) -> crate::Result<()>
    where
        B: Buffer,
    {
        unsafe {
            let ec = ffi::MPI_File_write_at_all(
                self.handle,
                offset,
                buf.pointer(),
                buf.count(),
                buf.as_datatype().as_raw(),
                ffi::RSMPI_STATUS_IGNORE,
            );
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_write_at_all"));
            }
            Ok(())
        }
    }

    pub fn read_at_all<B>(&self, offset: Offset, buf: &mut B) -> crate::Result<()>
    where
        B: BufferMut,
    {
        unsafe {
            let ec = ffi::MPI_File_read_at_all(
                self.handle,
                offset,
                buf.pointer_mut(),
                buf.count(),
                buf.as_datatype().as_raw(),
                ffi::RSMPI_STATUS_IGNORE,
            );
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_read_at_all"));
            }
            Ok(())
        }
    }

    pub fn set_atomicity(&self, flag: bool) -> crate::Result<()> {
        unsafe {
            let ec = ffi::MPI_File_set_atomicity(self.handle, flag as c_int);
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_set_atomicity"));
            }
            Ok(())
        }
    }

    pub fn sync(&self) -> crate::Result<()> {
        unsafe {
            let ec = ffi::MPI_File_sync(self.handle);
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_sync"));
            }
            Ok(())
        }
    }

    pub fn set_view<E, D>(
        &self,
        disp: Offset,
        etype: &E,
        filetype: &D,
        datarep: &str,
    ) -> crate::Result<()>
    where
        E: Datatype,
        D: Datatype,
    {
        let datarep = std::ffi::CString::new(datarep).unwrap();
        unsafe {
            let ec = ffi::MPI_File_set_view(
                self.handle,
                disp,
                etype.as_raw(),
                filetype.as_raw(),
                datarep.as_ptr(),
                ffi::RSMPI_INFO_NULL,
            );
            if ec != ffi::MPI_SUCCESS as c_int {
                return Err(crate::Error::mpi_error(ec, "MPI_File_set_view"));
            }
            Ok(())
        }
    }
}
