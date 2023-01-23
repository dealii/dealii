// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/mpi.h>

#include <deal.II/lac/petsc_vector.h>

#ifdef DEAL_II_WITH_PETSC

#  include <algorithm>
#  include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
  {
    Vector::Vector()
    {
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      Vector::create_vector(MPI_COMM_SELF, 0, 0);
    }



    Vector::Vector(const MPI_Comm &communicator,
                   const size_type n,
                   const size_type locally_owned_size)
    {
      Vector::create_vector(communicator, n, locally_owned_size);
    }



    Vector::Vector(const IndexSet &local,
                   const IndexSet &ghost,
                   const MPI_Comm &communicator)
    {
      Assert(local.is_ascending_and_one_to_one(communicator),
             ExcNotImplemented());

      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);

      Vector::create_vector(communicator,
                            local.size(),
                            local.n_elements(),
                            ghost_set);
    }



    Vector::Vector(const Vector &v)
      : VectorBase()
    {
      if (v.has_ghost_elements())
        Vector::create_vector(v.get_mpi_communicator(),
                              v.size(),
                              v.locally_owned_size(),
                              v.ghost_indices);
      else
        Vector::create_vector(v.get_mpi_communicator(),
                              v.size(),
                              v.locally_owned_size());

      this->operator=(v);
    }



    Vector::Vector(const IndexSet &local, const MPI_Comm &communicator)
    {
      Assert(local.is_ascending_and_one_to_one(communicator),
             ExcNotImplemented());
      Vector::create_vector(communicator, local.size(), local.n_elements());
    }



    Vector &
    Vector::operator=(const Vector &v)
    {
      // make sure left- and right-hand side of the assignment are
      // compress()'ed:
      Assert(v.last_action == VectorOperation::unknown,
             internal::VectorReference::ExcWrongMode(VectorOperation::unknown,
                                                     v.last_action));
      Assert(last_action == VectorOperation::unknown,
             internal::VectorReference::ExcWrongMode(VectorOperation::unknown,
                                                     last_action));

      // if the vectors have different sizes,
      // then first resize the present one
      if (size() != v.size())
        {
          if (v.has_ghost_elements())
            reinit(v.locally_owned_elements(),
                   v.ghost_indices,
                   v.get_mpi_communicator());
          else
            reinit(v.get_mpi_communicator(),
                   v.size(),
                   v.locally_owned_size(),
                   true);
        }

      PetscErrorCode ierr = VecCopy(v.vector, vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      if (has_ghost_elements())
        {
          ierr = VecGhostUpdateBegin(vector, INSERT_VALUES, SCATTER_FORWARD);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          ierr = VecGhostUpdateEnd(vector, INSERT_VALUES, SCATTER_FORWARD);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
        }
      return *this;
    }



    void
    Vector::clear()
    {
      VectorBase::clear();

      create_vector(MPI_COMM_SELF, 0, 0);
    }



    void
    Vector::reinit(const MPI_Comm &communicator,
                   const size_type n,
                   const size_type local_sz,
                   const bool      omit_zeroing_entries)
    {
      // only do something if the sizes
      // mismatch (may not be true for every proc)

      int k_global, k = ((size() != n) || (locally_owned_size() != local_sz));
      {
        const int ierr =
          MPI_Allreduce(&k, &k_global, 1, MPI_INT, MPI_LOR, communicator);
        AssertThrowMPI(ierr);
      }

      if (k_global || has_ghost_elements())
        {
          // FIXME: I'd like to use this here,
          // but somehow it leads to odd errors
          // somewhere down the line in some of
          // the tests:
          //         const PetscErrorCode ierr = VecSetSizes (vector, n, n);
          //         AssertThrow (ierr == 0, ExcPETScError(ierr));

          // so let's go the slow way:

          const PetscErrorCode ierr = VecDestroy(&vector);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          create_vector(communicator, n, local_sz);
        }

      // finally clear the new vector if so
      // desired
      if (omit_zeroing_entries == false)
        *this = 0;
    }



    void
    Vector::reinit(const Vector &v, const bool omit_zeroing_entries)
    {
      if (v.has_ghost_elements())
        {
          reinit(v.locally_owned_elements(),
                 v.ghost_indices,
                 v.get_mpi_communicator());
          if (!omit_zeroing_entries)
            {
              const PetscErrorCode ierr = VecSet(vector, 0.0);
              AssertThrow(ierr == 0, ExcPETScError(ierr));
            }
        }
      else
        reinit(v.get_mpi_communicator(),
               v.size(),
               v.locally_owned_size(),
               omit_zeroing_entries);
    }



    void
    Vector::reinit(const IndexSet &local,
                   const IndexSet &ghost,
                   const MPI_Comm &comm)
    {
      const PetscErrorCode ierr = VecDestroy(&vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      Assert(local.is_ascending_and_one_to_one(comm), ExcNotImplemented());

      IndexSet ghost_set = ghost;
      ghost_set.subtract_set(local);

      create_vector(comm, local.size(), local.n_elements(), ghost_set);
    }

    void
    Vector::reinit(const IndexSet &local, const MPI_Comm &comm)
    {
      const PetscErrorCode ierr = VecDestroy(&vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      Assert(local.is_ascending_and_one_to_one(comm), ExcNotImplemented());
      Assert(local.size() > 0, ExcMessage("can not create vector of size 0."));
      create_vector(comm, local.size(), local.n_elements());
    }

    void
    Vector::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
    {
      this->reinit(partitioner->locally_owned_range(),
                   partitioner->ghost_indices(),
                   partitioner->get_mpi_communicator());
    }


    void
    Vector::create_vector(const MPI_Comm &communicator,
                          const size_type n,
                          const size_type locally_owned_size)
    {
      (void)n;
      AssertIndexRange(locally_owned_size, n + 1);
      ghosted = false;

      const PetscErrorCode ierr = VecCreateMPI(communicator,
                                               locally_owned_size,
                                               PETSC_DETERMINE,
                                               &vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      Assert(size() == n, ExcDimensionMismatch(size(), n));
    }



    void
    Vector::create_vector(const MPI_Comm &communicator,
                          const size_type n,
                          const size_type locally_owned_size,
                          const IndexSet &ghostnodes)
    {
      (void)n;
      AssertIndexRange(locally_owned_size, n + 1);
      ghosted       = true;
      ghost_indices = ghostnodes;

      std::vector<size_type> ghostindices;
      ghostnodes.fill_index_vector(ghostindices);

      const PetscInt *ptr =
        (ghostindices.size() > 0 ?
           reinterpret_cast<const PetscInt *>(ghostindices.data()) :
           nullptr);

      PetscErrorCode ierr = VecCreateGhost(communicator,
                                           locally_owned_size,
                                           PETSC_DETERMINE,
                                           ghostindices.size(),
                                           ptr,
                                           &vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      Assert(size() == n, ExcDimensionMismatch(size(), n));

#  if DEBUG
      {
        // test ghost allocation in debug mode
        PetscInt begin, end;

        ierr = VecGetOwnershipRange(vector, &begin, &end);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        AssertDimension(locally_owned_size,
                        static_cast<size_type>(end - begin));

        Vec l;
        ierr = VecGhostGetLocalForm(vector, &l);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscInt lsize;
        ierr = VecGetSize(l, &lsize);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        ierr = VecGhostRestoreLocalForm(vector, &l);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        AssertDimension(lsize,
                        end - begin +
                          static_cast<PetscInt>(ghost_indices.n_elements()));
      }
#  endif
    }



    bool
    Vector::all_zero() const
    {
      unsigned int has_nonzero = VectorBase::all_zero() ? 0 : 1;
#  ifdef DEAL_II_WITH_MPI
      // in parallel, check that the vector
      // is zero on _all_ processors.
      unsigned int num_nonzero =
        Utilities::MPI::sum(has_nonzero, this->get_mpi_communicator());
      return num_nonzero == 0;
#  else
      return has_nonzero == 0;
#  endif
    }


    void
    Vector::print(std::ostream &     out,
                  const unsigned int precision,
                  const bool         scientific,
                  const bool         across) const
    {
      AssertThrow(out.fail() == false, ExcIO());

      // get a representation of the vector and
      // loop over all the elements
      const PetscScalar *val;
      PetscInt           nlocal, istart, iend;

      PetscErrorCode ierr = VecGetArrayRead(vector, &val);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      ierr = VecGetLocalSize(vector, &nlocal);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      ierr = VecGetOwnershipRange(vector, &istart, &iend);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      // save the state of out stream
      std::ios::fmtflags old_flags     = out.flags();
      unsigned int       old_precision = out.precision(precision);

      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      // let each processor produce its output in turn. this requires
      // synchronizing output between processors using a barrier --
      // which is clearly slow, but nobody is going to print a whole
      // matrix this way on a regular basis for production runs, so
      // the slowness of the barrier doesn't matter
      MPI_Comm communicator = this->get_mpi_communicator();
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(communicator);
           i++)
        {
          const int mpi_ierr = MPI_Barrier(communicator);
          AssertThrowMPI(mpi_ierr);

          if (i == Utilities::MPI::this_mpi_process(communicator))
            {
              if (across)
                {
                  out << "[Proc" << i << " " << istart << "-" << iend - 1 << "]"
                      << ' ';
                  for (PetscInt i = 0; i < nlocal; ++i)
                    out << val[i] << ' ';
                }
              else
                {
                  out << "[Proc " << i << " " << istart << "-" << iend - 1
                      << "]" << std::endl;
                  for (PetscInt i = 0; i < nlocal; ++i)
                    out << val[i] << std::endl;
                }
              out << std::endl;
            }
        }
      // reset output format
      out.flags(old_flags);
      out.precision(old_precision);

      // restore the representation of the
      // vector
      ierr = VecRestoreArrayRead(vector, &val);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      AssertThrow(out.fail() == false, ExcIO());
    }

  } // namespace MPI

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
