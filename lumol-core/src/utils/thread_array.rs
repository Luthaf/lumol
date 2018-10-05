// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::cell::{RefCell, RefMut};
use std::ops::AddAssign;

use thread_local::CachedThreadLocal;
use ndarray::Zip;

use crate::Array2;

/// A collection of vectors, one by thread using this struct. All the vectors
/// are wrapped in `RefCell` to ensure mutability of the values. The thread
/// creating the `ThreadLocalVec` will get faster access to the underlying data.
pub struct ThreadLocalArray2<T> where T: Send {
    inner: CachedThreadLocal<RefCell<Array2<T>>>,
    size: (usize, usize)
}

impl<T: Send + Default + Clone> ThreadLocalArray2<T> {
    /// Create a new `ThreadLocalArray2` with the given size, initializing the
    /// values with `T::default`.
    pub fn with_size(size: (usize, usize)) -> Self {
        let inner = CachedThreadLocal::new();
        // Set the current thread as owner of the data
        let _ = inner.get_or(|| Box::new(RefCell::new(Array2::default(size))));
        ThreadLocalArray2 {
            inner: inner,
            size: size,
        }
    }

    /// Mutably borrow the thread local vector if it already exists, or create
    /// it and then borrow it.
    pub fn borrow_mut(&self) -> RefMut<Array2<T>> {
        self.inner
            .get_or(|| Box::new(RefCell::new(Array2::default(self.size))))
            .borrow_mut()
    }
}

impl<T: Send> ThreadLocalArray2<T> {
    /// Get an iterator over all the vectors created by the different threads
    pub fn into_iter(self) -> impl Iterator<Item = Array2<T>> {
        self.inner
            .into_iter()
            .map(|cell| cell.into_inner())
    }

    /// Sum the values from all the vectors created by the different threads in
    /// the `output` buffer
    pub fn sum_into(self, output: &mut Array2<T>) where T: AddAssign + Copy {
        for local in self.into_iter() {
            Zip::from(&mut **output).and(&*local).apply(|a, &b| {
                *a += b;
            });
        }
    }
}
