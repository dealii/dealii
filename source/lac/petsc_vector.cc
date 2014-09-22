// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/petsc_vector.h>

#ifdef DEAL_II_WITH_PETSC

#  include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{


  Vector::Vector ()
  {
    Vector::create_vector (0);
  }



  Vector::Vector (const size_type n)
  {
    Vector::create_vector (n);
  }



  Vector::Vector (const Vector &v)
    :
    VectorBase ()
  {
    // first create a dummy vector, then copy
    // over the other one
    Vector::create_vector (1);
    Vector::operator = (v);
  }



  Vector::Vector (const MPI::Vector &v)
  {
    // first create a dummy vector, then copy
    // over the other one
    Vector::create_vector (1);
    Vector::operator = (v);
  }



  void
  Vector::reinit (const size_type n,
                  const bool      fast)
  {
    // only do something if the sizes
    // mismatch
    if (size() != n)
      {
        // FIXME: I'd like to use this here,
        // but somehow it leads to odd errors
        // somewhere down the line in some of
        // the tests:
//         const int ierr = VecSetSizes (vector, n, n);
//         AssertThrow (ierr == 0, ExcPETScError(ierr));

        // so let's go the slow way:
        if (attained_ownership)
          {
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
            int ierr = VecDestroy (vector);
#else
            int ierr = VecDestroy (&vector);
#endif
            AssertThrow (ierr == 0, ExcPETScError(ierr));
          }

        create_vector (n);
      }

    // finally clear the new vector if so
    // desired
    if (fast == false)
      *this = 0;
  }



  void
  Vector::reinit (const Vector &v,
                  const bool    fast)
  {
    reinit (v.size(), fast);
  }



  void
  Vector::create_vector (const size_type n)
  {
    const int ierr
      = VecCreateSeq (PETSC_COMM_SELF, n, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    attained_ownership = true;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
