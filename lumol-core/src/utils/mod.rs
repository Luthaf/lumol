// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Various internal utilities, which do not have there own module

#[macro_use]
mod macros;

mod thread_vec;
pub use self::thread_vec::ThreadLocalVec;

mod thread_array;
pub use self::thread_array::ThreadLocalArray2;

#[cfg(test)]
mod xyz;
#[cfg(test)]
pub use self::xyz::system_from_xyz;
