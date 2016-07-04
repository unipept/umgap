
use std::error;
use std::fmt;
use std::io;
use std::result;

extern crate csv;

extern crate fst;


pub type Result<T> = result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Csv(csv::Error),
    Fst(fst::Error),
    Io(io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Csv(ref err) => err.fmt(f),
            Error::Io(ref err) => err.fmt(f),
            Error::Fst(ref err) => err.fmt(f),
        }
    }
}

impl error::Error for Error {
    fn description(&self) -> &str {
        match *self {
            Error::Csv(ref err) => err.description(),
            Error::Io(ref err) => err.description(),
            Error::Fst(ref err) => err.description(),
        }
    }
}

impl From<csv::Error> for Error {
    fn from(err: csv::Error) -> Error {
        Error::Csv(err)
    }
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {
        Error::Io(err)
    }
}

impl From<fst::Error> for Error {
    fn from(err: fst::Error) -> Error {
        Error::Fst(err)
    }
}
