//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_vector.h>

#include <cmath>

#ifdef DEAL_II_USE_PETSC


namespace PETScWrappers
{


  Vector::Vector ()
  {
    Vector::create_vector (0);
  }



  Vector::Vector (const unsigned int n)
  {
    Vector::create_vector (n);
  }

  

  Vector::Vector (const VectorBase &v)
  {
    Vector::create_vector (v.size());
    VectorBase::operator = (v);
  }

  
  void
  Vector::create_vector (const unsigned int n,
                         const unsigned int local_size)
  {
    Assert (local_size < n, ExcIndexRange (local_size, 0, n));
    
    const int ierr
      = VecCreateSeq (PETSC_COMM_SELF, n, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
