// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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

#include <deal.II/lac/petsc_parallel_vector.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_vector.h>
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
      ghosted = false;
    }



    Vector::Vector (const MPI_Comm &communicator,
                    const size_type n,
                    const size_type local_size)
      :
      communicator (communicator)
    {
      Vector::create_vector (n, local_size);
    }



    Vector::Vector (const MPI_Comm   &communicator,
                    const VectorBase  &v,
                    const size_type   local_size)
      :
      communicator (communicator)
    {
      Vector::create_vector (v.size(), local_size);

      VectorBase::operator = (v);
    }



    Vector::Vector (const MPI_Comm     &communicator,
                    const IndexSet   &local,
                    const IndexSet &ghost)
      :
      communicator (communicator)
    {
      Assert(local.is_contiguous(), ExcNotImplemented());

      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);

      Vector::create_vector(local.size(), local.n_elements(), ghost_set);
    }

    Vector::Vector (const IndexSet   &local,
                    const IndexSet &ghost,
                    const MPI_Comm     &communicator)
      :
      communicator (communicator)
    {
      Assert(local.is_contiguous(), ExcNotImplemented());

      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);

      Vector::create_vector(local.size(), local.n_elements(), ghost_set);
    }


    Vector::Vector (const IndexSet   &local,
                    const MPI_Comm     &communicator)
      :
      communicator (communicator)
    {
      Assert(local.is_contiguous(), ExcNotImplemented());
      Vector::create_vector(local.size(), local.n_elements());
    }


    Vector::Vector (const MPI_Comm     &communicator,
                    const IndexSet   &local)
      :
      communicator (communicator)
    {
      Assert(local.is_contiguous(), ExcNotImplemented());
      Vector::create_vector(local.size(), local.n_elements());
    }

    void
    Vector::reinit (const MPI_Comm  &comm,
                    const size_type  n,
                    const size_type  local_sz,
                    const bool       fast)
    {
      communicator = comm;

      // only do something if the sizes
      // mismatch (may not be true for every proc)

      int k_global, k = ((size() != n) || (local_size() != local_sz));
      MPI_Allreduce (&k, &k_global, 1,
                     MPI_INT, MPI_LOR, communicator);

      if (k_global || has_ghost_elements())
        {
          // FIXME: I'd like to use this here,
          // but somehow it leads to odd errors
          // somewhere down the line in some of
          // the tests:
//         const int ierr = VecSetSizes (vector, n, n);
//         AssertThrow (ierr == 0, ExcPETScError(ierr));

          // so let's go the slow way:
          int ierr;

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
          ierr = VecDestroy (vector);
#else
          ierr = VecDestroy (&vector);
#endif

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
      if (v.has_ghost_elements())
        {
          reinit (v.locally_owned_elements(), v.ghost_indices, v.communicator);
          if (!fast)
            {
              int ierr = VecSet(vector, 0.0);
              AssertThrow (ierr == 0, ExcPETScError(ierr));
            }
        }
      else
        reinit (v.communicator, v.size(), v.local_size(), fast);
    }



    void
    Vector::reinit (const MPI_Comm     &comm,
                    const IndexSet   &local,
                    const IndexSet &ghost)
    {
      reinit(local, ghost, comm);
    }

    void
    Vector::reinit (const IndexSet   &local,
                    const IndexSet &ghost,
                    const MPI_Comm     &comm)
    {
      int ierr;
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      ierr = VecDestroy (vector);
#else
      ierr = VecDestroy (&vector);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      communicator = comm;

      Assert(local.is_contiguous(), ExcNotImplemented());

      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);

      create_vector(local.size(), local.n_elements(), ghost_set);
    }

    void
    Vector::reinit (const MPI_Comm     &comm,
                    const IndexSet   &local)
    {
      reinit(local, comm);
    }

    void
    Vector::reinit (const IndexSet &local,
                    const MPI_Comm &comm)
    {
      int ierr;
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      ierr = VecDestroy (vector);
#else
      ierr = VecDestroy (&vector);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      communicator = comm;

      Assert(local.is_contiguous(), ExcNotImplemented());
      Assert(local.size()>0, ExcMessage("can not create vector of size 0."));
      create_vector(local.size(), local.n_elements());
    }


    Vector &
    Vector::operator = (const PETScWrappers::Vector &v)
    {
      Assert(last_action==VectorOperation::unknown,
             ExcMessage("Call to compress() required before calling operator=."));
      //TODO [TH]: can not access v.last_action here. Implement is_compressed()?
      //Assert(v.last_action==VectorOperation::unknown,
      //    ExcMessage("Call to compress() required before calling operator=."));
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
      const std::pair<size_type, size_type>
      local_elements = local_range ();
      std::copy (src_array + local_elements.first,
                 src_array + local_elements.second,
                 dest_array);

      // finally restore the arrays
      ierr = VecRestoreArray (vector, &dest_array);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = VecRestoreArray (static_cast<const Vec &>(v), &src_array);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      if (has_ghost_elements())
        update_ghost_values();
      return *this;
    }


    void
    Vector::create_vector (const size_type n,
                           const size_type local_size)
    {
      Assert (local_size <= n, ExcIndexRange (local_size, 0, n));
      ghosted = false;

      const int ierr
        = VecCreateMPI (communicator, local_size, PETSC_DETERMINE,
                        &vector);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      Assert (size() == n,
              ExcDimensionMismatch (size(), n));
    }



    void
    Vector::create_vector (const size_type n,
                           const size_type local_size,
                           const IndexSet &ghostnodes)
    {
      Assert (local_size <= n, ExcIndexRange (local_size, 0, n));
      ghosted = true;
      ghost_indices = ghostnodes;

      std::vector<size_type> ghostindices;
      ghostnodes.fill_index_vector(ghostindices);

      const PetscInt *ptr
        = (ghostindices.size() > 0
           ?
           (const PetscInt *)(&(ghostindices[0]))
           :
           0);

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
      {
        // test ghost allocation in debug mode
        PetscInt begin, end;

        ierr = VecGetOwnershipRange (vector, &begin, &end);

        Assert(local_size==(size_type)(end-begin), ExcInternalError());

        Vec l;
        ierr = VecGhostGetLocalForm(vector, &l);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        PetscInt lsize;
        ierr = VecGetSize(l, &lsize);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        ierr = VecGhostRestoreLocalForm(vector, &l);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        Assert (lsize==end-begin+(PetscInt)ghost_indices.n_elements(),
                ExcInternalError());
      }
#endif


      // in PETSc versions up to 3.5, VecCreateGhost zeroed out the locally
      // owned vector elements but forgot about the ghost elements. we need to
      // do this ourselves
      //
      // see https://code.google.com/p/dealii/issues/detail?id=233
#if DEAL_II_PETSC_VERSION_LT(3,6,0)
      PETScWrappers::MPI::Vector zero;
      zero.reinit (communicator, this->size(), local_size);
      *this = zero;
#endif

    }



    bool
    Vector::all_zero() const
    {
      unsigned int has_nonzero = VectorBase::all_zero()?0:1;
#ifdef DEAL_II_WITH_MPI
      // in parallel, check that the vector
      // is zero on _all_ processors.
      unsigned int num_nonzero = Utilities::MPI::sum(has_nonzero, communicator);
      return num_nonzero == 0;
#else
      return has_nonzero == 0;
#endif
    }


    void
    Vector::print (std::ostream      &out,
                   const unsigned int precision,
                   const bool         scientific,
                   const bool         across) const
    {
      AssertThrow (out, ExcIO());

      // get a representation of the vector and
      // loop over all the elements
      PetscScalar *val;
      PetscInt    nlocal, istart, iend;

      int ierr = VecGetArray (vector, &val);

      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = VecGetLocalSize (vector, &nlocal);

      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = VecGetOwnershipRange (vector, &istart, &iend);

      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // save the state of out stream
      std::ios::fmtflags old_flags = out.flags();
      unsigned int old_precision = out.precision (precision);

      out.precision (precision);
      if (scientific)
        out.setf (std::ios::scientific, std::ios::floatfield);
      else
        out.setf (std::ios::fixed, std::ios::floatfield);

      for ( unsigned int i = 0;
            i < Utilities::MPI::n_mpi_processes(communicator);
            i++)
        {
          // This is slow, but most likely only used to debug.
          MPI_Barrier(communicator);
          if (i == Utilities::MPI::this_mpi_process(communicator))
            {
              if (across)
                {
                  out << "[Proc" << i << " " << istart << "-" << iend-1 << "]" << ' ';
                  for (PetscInt i=0; i<nlocal; ++i)
                    out << val[i] << ' ';
                }
              else
                {
                  out << "[Proc " << i << " " << istart << "-" << iend-1 << "]" << std::endl;
                  for (PetscInt i=0; i<nlocal; ++i)
                    out << val[i] << std::endl;
                }
              out << std::endl;
            }
        }
      // reset output format
      out.flags (old_flags);
      out.precision(old_precision);

      // restore the representation of the
      // vector
      ierr = VecRestoreArray (vector, &val);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      AssertThrow (out, ExcIO());
    }

  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
