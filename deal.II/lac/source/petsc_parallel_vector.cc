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


#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_vector.h>

#include <cmath>
#include <algorithm>

#ifdef DEAL_II_USE_PETSC


namespace PETScWrappers
{
  namespace MPI
  {

    Vector::Vector ()
    {
                                       // this is an invalid empty vector, so we
                                       // can just as well create a sequential
                                       // one to avoid all the overhead incurred
                                       // by parallelism
      const int n = 0;
      const int ierr
        = VecCreateSeq (PETSC_COMM_SELF, n, &vector);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }



    Vector::Vector (const MPI_Comm    &communicator,
                    const unsigned int n,
                    const unsigned int local_size)
                    :
                    communicator (communicator)
    {
      Vector::create_vector (n, local_size);
    }

  

    Vector::Vector (const MPI_Comm    &communicator,
                    const VectorBase  &v,
                    const unsigned int local_size)
                    :
                    communicator (communicator)
    {
      Vector::create_vector (v.size(), local_size);

      VectorBase::operator = (v);
    }

  

    void
    Vector::reinit (const MPI_Comm    &comm,
                    const unsigned int n,
                    const unsigned int local_sz,
                    const bool         fast)
    {
      communicator = comm;
      
                                       // only do something if the sizes
                                       // mismatch
      if ((size() != n) || (local_size() != local_sz))
        {
                                           // FIXME: I'd like to use this here,
                                           // but somehow it leads to odd errors
                                           // somewhere down the line in some of
                                           // the tests:
//         const int ierr = VecSetSizes (vector, n, n);
//         AssertThrow (ierr == 0, ExcPETScError(ierr));

                                           // so let's go the slow way:
          int ierr;
          ierr = VecDestroy (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          create_vector (n, local_sz);
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
      communicator = v.communicator;
      
      reinit (communicator, v.size(), v.local_size(), fast);
    }



    Vector &
    Vector::operator = (const PETScWrappers::Vector &v)
    {
                                       // first flush buffers
      compress ();

      int ierr;

                                       // get a pointer to the local memory of
                                       // this vector
      PetscScalar *dest_array;
      ierr = VecGetArray (vector, &dest_array);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

                                       // then also a pointer to the source
                                       // vector
      PetscScalar *src_array;
      ierr = VecGetArray (static_cast<const Vec &>(v), &src_array);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

                                       // then copy:
      const std::pair<unsigned int, unsigned int>
        local_elements = local_range ();
      std::copy (src_array + local_elements.first,
                 src_array + local_elements.second,
		 dest_array);

                                       // finally restore the arrays
      ierr = VecRestoreArray (vector, &dest_array);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = VecRestoreArray (static_cast<const Vec &>(v), &src_array);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      return *this;
    }


    void
    Vector::create_vector (const unsigned int  n,
                           const unsigned int  local_size)
    {
      Assert (local_size <= n, ExcIndexRange (local_size, 0, n));

      const int ierr
        = VecCreateMPI (communicator, local_size, PETSC_DETERMINE,
                        &vector);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      Assert (size() == n, ExcInternalError());
    }

  }

}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
