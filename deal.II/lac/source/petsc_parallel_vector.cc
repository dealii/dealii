//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2006, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_parallel_vector.h>

#ifdef DEAL_II_USE_PETSC

#  include <lac/petsc_vector.h>
#  include <cmath>
#  include <algorithm>

DEAL_II_NAMESPACE_OPEN

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



    Vector::Vector (const MPI_Comm     &communicator,
		    const IndexSet &  local,
		    const IndexSet & ghost)
                    :
                    communicator (communicator)
    {
      Assert(local.is_contiguous(), ExcNotImplemented());
      
      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);
      
				       //possible optmization: figure out if
				       //there are ghost indices (collective
				       //operation!) and then create a
				       //non-ghosted vector.
//      Vector::create_vector (local.size(), local.n_elements());
      
      Vector::create_vector(local.size(), local.n_elements(), ghost_set);    
    }
    


    void
    Vector::reinit (const MPI_Comm    &comm,
                    const unsigned int n,
                    const unsigned int local_sz,
                    const bool         fast)
    {
      communicator = comm;
      
                                       // only do something if the sizes
                                       // mismatch (may not be true for every proc)
      
      int k_global, k = ((size() != n) || (local_size() != local_sz));
      MPI_Allreduce (&k, &k_global, 1,
	MPI_INT, MPI_LOR, communicator);

      if (k_global)
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


    
    void
    Vector::reinit (const MPI_Comm     &comm,
		    const IndexSet &  local,
		    const IndexSet & ghost)
    {
      communicator = comm;

      Assert(local.is_contiguous(), ExcNotImplemented());
      
      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);

      create_vector(local.size(), local.n_elements(), ghost_set); 
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

      Assert (size() == n,
	      ExcDimensionMismatch (size(), n));
    }


    
    void
    Vector::create_vector (const unsigned int  n,
                           const unsigned int  local_size,
			   const IndexSet & ghostnodes)
    {
      Assert (local_size <= n, ExcIndexRange (local_size, 0, n));
      ghosted = true;
      ghost_indices = ghostnodes;
      
				       //64bit indices won't work yet:
      Assert (sizeof(unsigned int)==sizeof(PetscInt), ExcInternalError());

      
      std::vector<unsigned int> ghostindices;
      ghostnodes.fill_index_vector(ghostindices);
      
      const PetscInt * ptr= (const PetscInt*)(&(ghostindices[0]));

      int ierr
	= VecCreateGhost(communicator,
			 local_size,
			 PETSC_DETERMINE,
			 ghostindices.size(),
			 ptr,
			 &vector);
      
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      Assert (size() == n,
	      ExcDimensionMismatch (size(), n));

#if DEBUG
				       // test ghost allocation in debug mode

#ifdef PETSC_USE_64BIT_INDICES
      PetscInt
#else
	int
#endif
	begin, end;

      ierr = VecGetOwnershipRange (vector, &begin, &end);

      Assert(local_size==(unsigned int)(end-begin), ExcInternalError());

      Vec l;
      ierr = VecGhostGetLocalForm(vector, &l);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      PetscInt lsize;
      ierr = VecGetSize(l, &lsize);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
	      
      ierr = VecGhostRestoreLocalForm(vector, &l);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
      
      Assert( lsize==end-begin+(PetscInt)ghost_indices.n_elements() ,ExcInternalError());

#endif

      
    }

    

  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
