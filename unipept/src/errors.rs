//! Defines the possible errors in Unipept.

use std::error;
use std::fmt;
use std::io;
use std::result;
use std::num;

extern crate csv;

extern crate fst;

/// Represents a result with a possible Error.
pub type Result<T> = result::Result<T, Error>;

/// Represents an error in Unipept.
#[allow(missing_docs)]
#[derive(Debug)]
pub enum Error {
    Csv(csv::Error),
    Fst(fst::Error),
    Io(io::Error),
    Str(&'static str),
    Parse(num::ParseIntError)
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Csv(ref err) => err.fmt(f),
            Error::Io(ref err) => err.fmt(f),
            Error::Fst(ref err) => err.fmt(f),
            Error::Str(ref err) => write!(f, "{}", err),
            Error::Parse(ref err) => err.fmt(f),
        }
    }
}

impl error::Error for Error {
    fn description(&self) -> &str {
        match *self {
            Error::Csv(ref err) => err.description(),
            Error::Io(ref err) => err.description(),
            Error::Fst(ref err) => err.description(),
            Error::Str(ref err) => err,
            Error::Parse(ref err) => err.description(),
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

impl From<&'static str> for Error {
    fn from(err: &'static str) -> Error {
        Error::Str(err)
    }
}

impl From<num::ParseIntError> for Error {
    fn from(err: num::ParseIntError) -> Error {
        Error::Parse(err)
    }
}
