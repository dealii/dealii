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
    const int n = 0;
    const int ierr
      = VecCreateSeq (PETSC_COMM_SELF, n, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  Vector::Vector (const unsigned int n)
  {
    const int ierr
      = VecCreateSeq (PETSC_COMM_SELF, n, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  

  Vector::Vector (const VectorBase &v)
  {
    int ierr
      = VecCreateSeq (PETSC_COMM_SELF, v.size(), &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    VectorBase::operator = (v);
  }

  
  void
  Vector::create_vector (const unsigned int n)
  {
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
