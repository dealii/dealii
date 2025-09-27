// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/trilinos_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/mpi.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/read_write_vector.h>
#  include <deal.II/lac/trilinos_index_access.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

#  include <boost/io/ios_state.hpp>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_Export.h>
#  include <Epetra_Import.h>
#  include <Epetra_Vector.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <cmath>
#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
#  ifndef DOXYGEN
  namespace internal
  {
    VectorReference::operator TrilinosScalar() const
    {
      AssertIndexRange(index, vector.size());

      // Trilinos allows for vectors to be referenced by the [] or ()
      // operators but only () checks index bounds. We check these bounds by
      // ourselves, so we can use []. Note that we can only get local values.

      const TrilinosWrappers::types::int_type local_index =
        vector.vector->Map().LID(
          static_cast<TrilinosWrappers::types::int_type>(index));
#    ifndef DEAL_II_WITH_64BIT_INDICES
      Assert(local_index >= 0,
             MPI::Vector::ExcAccessToNonLocalElement(
               index,
               vector.vector->Map().NumMyElements(),
               vector.vector->Map().MinMyGID(),
               vector.vector->Map().MaxMyGID()));
#    else
      Assert(local_index >= 0,
             MPI::Vector::ExcAccessToNonLocalElement(
               index,
               vector.vector->Map().NumMyElements(),
               vector.vector->Map().MinMyGID64(),
               vector.vector->Map().MaxMyGID64()));
#    endif


      return (*(vector.vector))[0][local_index];
    }
  } // namespace internal
#  endif

  namespace MPI
  {
    Vector::Vector()
      : last_action(Zero)
      , compressed(true)
      , has_ghosts(false)
      , vector(new Epetra_FEVector(
          Epetra_Map(0, 0, 0, Utilities::Trilinos::comm_self())))
    {}



    Vector::Vector(const IndexSet &parallel_partitioning,
                   const MPI_Comm  communicator)
      : Vector()
    {
      reinit(parallel_partitioning, communicator);
    }



    Vector::Vector(const Vector &v)
      : Vector()
    {
      has_ghosts     = v.has_ghosts;
      vector         = std::make_unique<Epetra_FEVector>(*v.vector);
      owned_elements = v.owned_elements;
    }



    Vector::Vector(Vector &&v) // NOLINT
      : Vector()
    {
      // initialize a minimal, valid object and swap
      static_cast<EnableObserverPointer &>(*this) =
        static_cast<EnableObserverPointer &&>(v);
      swap(v);
    }



    Vector::Vector(const IndexSet &parallel_partitioner,
                   const Vector   &v,
                   const MPI_Comm  communicator)
      : Vector()
    {
      AssertThrow(parallel_partitioner.size() ==
                    static_cast<size_type>(
                      TrilinosWrappers::n_global_elements(v.vector->Map())),
                  ExcDimensionMismatch(parallel_partitioner.size(),
                                       TrilinosWrappers::n_global_elements(
                                         v.vector->Map())));

      vector = std::make_unique<Epetra_FEVector>(
        parallel_partitioner.make_trilinos_map(communicator, true));
      reinit(v, false, true);
    }



    Vector::Vector(const IndexSet &local,
                   const IndexSet &ghost,
                   const MPI_Comm  communicator)
      : Vector()
    {
      reinit(local, ghost, communicator, false);
    }



    void
    Vector::clear()
    {
      // When we clear the vector, reset the pointer and generate an empty
      // vector.
      Epetra_Map map(0, 0, Epetra_MpiComm(MPI_COMM_SELF));

      has_ghosts  = false;
      vector      = std::make_unique<Epetra_FEVector>(map);
      last_action = Zero;
    }



    void
    Vector::reinit(const IndexSet &parallel_partitioner,
                   const MPI_Comm  communicator,
                   const bool /*omit_zeroing_entries*/)
    {
      nonlocal_vector.reset();

      const bool overlapping =
        !parallel_partitioner.is_ascending_and_one_to_one(communicator);

      Epetra_Map map =
        parallel_partitioner.make_trilinos_map(communicator, overlapping);

      vector = std::make_unique<Epetra_FEVector>(map);

      has_ghosts = vector->Map().UniqueGIDs() == false;

      // If the IndexSets are overlapping, we don't really know
      // which process owns what. So we decide that no process
      // owns anything in that case. In particular asking for
      // the locally owned elements is not allowed.
      if (has_ghosts)
        {
          owned_elements.clear();
          owned_elements.set_size(0);
        }
      else
        owned_elements = parallel_partitioner;

      if constexpr (running_in_debug_mode())
        {
          const size_type n_elements_global =
            Utilities::MPI::sum(owned_elements.n_elements(), communicator);

          Assert(has_ghosts || n_elements_global == size(), ExcInternalError());
        }

      last_action = Zero;
    }



    void
    Vector::reinit(const Vector &v,
                   const bool    omit_zeroing_entries,
                   const bool    allow_different_maps)
    {
      nonlocal_vector.reset();

      // In case we do not allow to have different maps, this call means that
      // we have to reset the vector. So clear the vector, initialize our map
      // with the map in v, and generate the vector.
      if (allow_different_maps == false)
        {
          // check equality for MPI communicators: We can only choose the fast
          // version in case the underlying Epetra_MpiComm object is the same,
          // otherwise we might access an MPI_Comm object that has been
          // deleted
          const Epetra_MpiComm *my_comm =
            dynamic_cast<const Epetra_MpiComm *>(&vector->Comm());
          const Epetra_MpiComm *v_comm =
            dynamic_cast<const Epetra_MpiComm *>(&v.vector->Comm());
          const bool same_communicators =
            my_comm != nullptr && v_comm != nullptr &&
            my_comm->DataPtr() == v_comm->DataPtr();
          if (!same_communicators ||
              vector->Map().SameAs(v.vector->Map()) == false)
            {
              vector      = std::make_unique<Epetra_FEVector>(v.vector->Map());
              has_ghosts  = v.has_ghosts;
              last_action = Zero;
              owned_elements = v.owned_elements;
            }
          else if (omit_zeroing_entries == false)
            {
              // old and new vectors have exactly the same map, i.e. size and
              // parallel distribution
              int ierr = vector->GlobalAssemble(last_action);
              Assert(ierr == 0, ExcTrilinosError(ierr));

              ierr = vector->PutScalar(0.0);
              Assert(ierr == 0, ExcTrilinosError(ierr));

              last_action = Zero;
            }
        }

      // Otherwise, we have to check that the two vectors are already of the
      // same size, create an object for the data exchange and then insert all
      // the data. The first assertion is only a check whether the user knows
      // what they are doing.
      else
        {
          Assert(omit_zeroing_entries == false,
                 ExcMessage(
                   "It is not possible to exchange data with the "
                   "option 'omit_zeroing_entries' set, which would not write "
                   "elements."));

          AssertThrow(size() == v.size(),
                      ExcDimensionMismatch(size(), v.size()));

          Epetra_Import data_exchange(vector->Map(), v.vector->Map());

          const int ierr = vector->Import(*v.vector, data_exchange, Insert);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));

          last_action = Insert;
        }
      if constexpr (running_in_debug_mode())
        {
          const Epetra_MpiComm *comm_ptr =
            dynamic_cast<const Epetra_MpiComm *>(&(v.vector->Comm()));
          Assert(comm_ptr != nullptr, ExcInternalError());
          const size_type n_elements_global =
            Utilities::MPI::sum(owned_elements.n_elements(), comm_ptr->Comm());
          Assert(has_ghosts || n_elements_global == size(), ExcInternalError());
        }
    }



    void
    Vector::reinit(const MPI::BlockVector &v, const bool import_data)
    {
      nonlocal_vector.reset();
      owned_elements.clear();
      owned_elements.set_size(v.size());

      // In case we do not allow to have different maps, this call means that
      // we have to reset the vector. So clear the vector, initialize our map
      // with the map in v, and generate the vector.
      if (v.n_blocks() == 0)
        return;

      // create a vector that holds all the elements contained in the block
      // vector. need to manually create an Epetra_Map.
      size_type n_elements = 0, added_elements = 0, block_offset = 0;
      for (size_type block = 0; block < v.n_blocks(); ++block)
        n_elements += v.block(block).vector->Map().NumMyElements();
      std::vector<TrilinosWrappers::types::int_type> global_ids(n_elements, -1);
      for (size_type block = 0; block < v.n_blocks(); ++block)
        {
          TrilinosWrappers::types::int_type *glob_elements =
            TrilinosWrappers::my_global_elements(
              v.block(block).trilinos_partitioner());
          size_type vector_size = v.block(block).vector->Map().NumMyElements();
          for (size_type i = 0; i < vector_size; ++i)
            global_ids[added_elements++] = glob_elements[i] + block_offset;
          owned_elements.add_indices(v.block(block).owned_elements,
                                     block_offset);
          block_offset += v.block(block).size();
        }

      Assert(n_elements == added_elements, ExcInternalError());
      Epetra_Map new_map(v.size(),
                         n_elements,
                         global_ids.data(),
                         0,
                         v.block(0).trilinos_partitioner().Comm());

      auto actual_vec = std::make_unique<Epetra_FEVector>(new_map);

      TrilinosScalar *entries = (*actual_vec)[0];
      for (size_type block = 0; block < v.n_blocks(); ++block)
        {
          v.block(block).trilinos_vector().ExtractCopy(entries, 0);
          entries += v.block(block).vector->Map().NumMyElements();
        }

      if (import_data == true)
        {
          AssertThrow(static_cast<size_type>(TrilinosWrappers::global_length(
                        *actual_vec)) == v.size(),
                      ExcDimensionMismatch(TrilinosWrappers::global_length(
                                             *actual_vec),
                                           v.size()));

          Epetra_Import data_exchange(vector->Map(), actual_vec->Map());

          const int ierr = vector->Import(*actual_vec, data_exchange, Insert);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));

          last_action = Insert;
        }
      else
        vector = std::move(actual_vec);
      if constexpr (running_in_debug_mode())
        {
          const Epetra_MpiComm *comm_ptr =
            dynamic_cast<const Epetra_MpiComm *>(&(vector->Comm()));
          Assert(comm_ptr != nullptr, ExcInternalError());
          const size_type n_elements_global =
            Utilities::MPI::sum(owned_elements.n_elements(), comm_ptr->Comm());

          Assert(has_ghosts || n_elements_global == size(), ExcInternalError());
        }
    }



    void
    Vector::reinit(const IndexSet &locally_owned_entries,
                   const IndexSet &ghost_entries,
                   const MPI_Comm  communicator,
                   const bool      vector_writable)
    {
      nonlocal_vector.reset();
      owned_elements = locally_owned_entries;
      if (vector_writable == false)
        {
          IndexSet parallel_partitioner = locally_owned_entries;
          parallel_partitioner.add_indices(ghost_entries);
          Epetra_Map map =
            parallel_partitioner.make_trilinos_map(communicator, true);
          vector = std::make_unique<Epetra_FEVector>(map);
        }
      else
        {
          Epetra_Map map =
            locally_owned_entries.make_trilinos_map(communicator, true);
          Assert(map.IsOneToOne(),
                 ExcMessage("A writable vector must not have ghost entries in "
                            "its parallel partitioning"));

          if (vector->Map().SameAs(map) == false)
            vector = std::make_unique<Epetra_FEVector>(map);
          else
            {
              const int ierr = vector->PutScalar(0.);
              Assert(ierr == 0, ExcTrilinosError(ierr));
            }

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);
          if (Utilities::MPI::n_mpi_processes(communicator) > 1)
            {
              Epetra_Map nonlocal_map =
                nonlocal_entries.make_trilinos_map(communicator, true);
              nonlocal_vector =
                std::make_unique<Epetra_MultiVector>(nonlocal_map, 1);
            }
        }

      has_ghosts = vector->Map().UniqueGIDs() == false;

      last_action = Zero;

      if constexpr (running_in_debug_mode())
        {
          const size_type n_elements_global =
            Utilities::MPI::sum(owned_elements.n_elements(), communicator);

          Assert(has_ghosts || n_elements_global == size(), ExcInternalError());
        }
    }



    void
    Vector::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const bool                                                make_ghosted,
      const bool                                                vector_writable)
    {
      if (make_ghosted)
        {
          Assert(partitioner->ghost_indices_initialized(),
                 ExcMessage("You asked to create a ghosted vector, but the "
                            "partitioner does not provide ghost indices."));

          this->reinit(partitioner->locally_owned_range(),
                       partitioner->ghost_indices(),
                       partitioner->get_mpi_communicator(),
                       vector_writable);
        }
      else
        {
          this->reinit(partitioner->locally_owned_range(),
                       partitioner->get_mpi_communicator());
        }
    }



    Vector &
    Vector::operator=(const Vector &v)
    {
      Assert(vector.get() != nullptr,
             ExcMessage("Vector is not constructed properly."));

      // check equality for MPI communicators to avoid accessing a possibly
      // invalid MPI_Comm object
      const Epetra_MpiComm *my_comm =
        dynamic_cast<const Epetra_MpiComm *>(&vector->Comm());
      const Epetra_MpiComm *v_comm =
        dynamic_cast<const Epetra_MpiComm *>(&v.vector->Comm());
      const bool same_communicators = my_comm != nullptr && v_comm != nullptr &&
                                      my_comm->DataPtr() == v_comm->DataPtr();
      // Need to ask MPI whether the communicators are the same. We would like
      // to use the following checks but currently we cannot make sure the
      // memory of my_comm is not stale from some MPI_Comm_free
      // somewhere. This can happen when a vector lives in GrowingVectorMemory
      // data structures. Thus, the following code is commented out.
      //
      // if (my_comm != nullptr &&
      //     v_comm != nullptr &&
      //     my_comm->DataPtr() != v_comm->DataPtr())
      //  {
      //    int communicators_same = 0;
      //    const int ierr = MPI_Comm_compare (my_comm->GetMpiComm(),
      //                                       v_comm->GetMpiComm(),
      //                                       &communicators_same);
      //    AssertThrowMPI(ierr);
      //    if (!(communicators_same == MPI_IDENT ||
      //          communicators_same == MPI_CONGRUENT))
      //      same_communicators = false;
      //    else
      //      same_communicators = true;
      //  }

      // distinguish three cases. First case: both vectors have the same
      // layout (just need to copy the local data, not reset the memory and
      // the underlying Epetra_Map). The third case means that we have to
      // rebuild the calling vector.
      if (same_communicators && v.vector->Map().SameAs(vector->Map()))
        {
          *vector = *v.vector;
          if (v.nonlocal_vector.get() != nullptr)
            nonlocal_vector =
              std::make_unique<Epetra_MultiVector>(v.nonlocal_vector->Map(), 1);
          last_action = Zero;
        }
      // Second case: vectors have the same global
      // size, but different parallel layouts (and
      // one of them a one-to-one mapping). Then we
      // can call the import/export functionality.
      else if (size() == v.size() &&
               (v.vector->Map().UniqueGIDs() || vector->Map().UniqueGIDs()))
        {
          reinit(v, false, true);
        }
      // Third case: Vectors do not have the same
      // size.
      else
        {
          vector         = std::make_unique<Epetra_FEVector>(*v.vector);
          last_action    = Zero;
          has_ghosts     = v.has_ghosts;
          owned_elements = v.owned_elements;
        }

      if (v.nonlocal_vector.get() != nullptr)
        nonlocal_vector =
          std::make_unique<Epetra_MultiVector>(v.nonlocal_vector->Map(), 1);

      return *this;
    }



    Vector &
    Vector::operator=(Vector &&v) noexcept
    {
      static_cast<EnableObserverPointer &>(*this) =
        static_cast<EnableObserverPointer &&>(v);
      swap(v);
      return *this;
    }



    template <typename number>
    Vector &
    Vector::operator=(const ::dealii::Vector<number> &v)
    {
      Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

      // this is probably not very efficient but works. in particular, we could
      // do better if we know that number==TrilinosScalar because then we could
      // elide the copying of elements
      //
      // let's hope this isn't a particularly frequent operation
      std::pair<size_type, size_type> local_range = this->local_range();
      for (size_type i = local_range.first; i < local_range.second; ++i)
        (*vector)[0][i - local_range.first] = v(i);

      return *this;
    }



    void
    Vector::import_nonlocal_data_for_fe(const TrilinosWrappers::SparseMatrix &m,
                                        const Vector                         &v)
    {
      Assert(m.trilinos_matrix().Filled() == true,
             ExcMessage("Matrix is not compressed. "
                        "Cannot find exchange information!"));
      Assert(v.vector->Map().UniqueGIDs() == true,
             ExcMessage("The input vector has overlapping data, "
                        "which is not allowed."));

      if (vector->Map().SameAs(m.trilinos_matrix().ColMap()) == false)
        vector =
          std::make_unique<Epetra_FEVector>(m.trilinos_matrix().ColMap());

      Epetra_Import data_exchange(vector->Map(), v.vector->Map());
      const int     ierr = vector->Import(*v.vector, data_exchange, Insert);

      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      last_action = Insert;
    }


    void
    Vector::import_elements(const LinearAlgebra::ReadWriteVector<double> &rwv,
                            const VectorOperation::values operation)
    {
      Assert(
        this->size() == rwv.size(),
        ExcMessage(
          "Both vectors need to have the same size for import_elements() to work!"));
      // TODO: a generic import_elements() function should handle any kind of
      // data layout in ReadWriteVector, but this function is of limited use as
      // this class will (hopefully) be retired eventually.
      Assert(this->locally_owned_elements() == rwv.get_stored_elements(),
             ExcNotImplemented());

      if (operation == VectorOperation::insert)
        {
          for (const auto idx : this->locally_owned_elements())
            (*this)[idx] = rwv[idx];
        }
      else if (operation == VectorOperation::add)
        {
          for (const auto idx : this->locally_owned_elements())
            (*this)[idx] += rwv[idx];
        }
      else
        AssertThrow(false, ExcNotImplemented());

      this->compress(operation);
    }


    void
    Vector::compress(VectorOperation::values given_last_action)
    {
      Assert(has_ghost_elements() == false,
             ExcMessage(
               "Calling compress() is only useful if a vector "
               "has been written into, but this is a vector with ghost "
               "elements and consequently is read-only. It does "
               "not make sense to call compress() for such "
               "vectors."));

      // Select which mode to send to Trilinos. Note that we use last_action if
      // available and ignore what the user tells us to detect wrongly mixed
      // operations. Typically given_last_action is only used on machines that
      // do not execute an operation (because they have no own cells for
      // example).
      Epetra_CombineMode mode = last_action;
      if (last_action == Zero)
        {
          if (given_last_action == VectorOperation::add)
            mode = Add;
          else if (given_last_action == VectorOperation::insert)
            mode = Insert;
          else
            Assert(
              false,
              ExcMessage(
                "compress() can only be called with VectorOperation add, insert, or unknown"));
        }
      else
        {
          Assert(
            ((last_action == Add) &&
             (given_last_action == VectorOperation::add)) ||
              ((last_action == Insert) &&
               (given_last_action == VectorOperation::insert)),
            ExcMessage(
              "The last operation on the Vector and the given last action in the compress() call do not agree!"));
        }


      if constexpr (running_in_debug_mode())
        {
          // check that every process has decided to use the same mode. This
          // will otherwise result in undefined behavior in the call to
          // GlobalAssemble().
          const double          double_mode = mode;
          const Epetra_MpiComm *comm_ptr = dynamic_cast<const Epetra_MpiComm *>(
            &(trilinos_partitioner().Comm()));
          Assert(comm_ptr != nullptr, ExcInternalError());

          const Utilities::MPI::MinMaxAvg result =
            Utilities::MPI::min_max_avg(double_mode, comm_ptr->GetMpiComm());
          Assert(result.max == result.min,
                 ExcMessage(
                   "Not all processors agree whether the last operation on "
                   "this vector was an addition or a set operation. This will "
                   "prevent the compress() operation from succeeding."));
        }

      // Now pass over the information about what we did last to the vector.
      if (nonlocal_vector.get() == nullptr || mode != Add)
        {
          const auto ierr = vector->GlobalAssemble(mode);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        }
      else
        {
          Epetra_Export exporter(nonlocal_vector->Map(), vector->Map());

          int ierr = vector->Export(*nonlocal_vector, exporter, mode);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));

          ierr = nonlocal_vector->PutScalar(0.);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        }
      last_action = Zero;

      compressed = true;
    }



    TrilinosScalar
    Vector::operator()(const size_type index) const
    {
      // Extract local indices in the vector.
      TrilinosWrappers::types::int_type trilinos_i = vector->Map().LID(
        static_cast<TrilinosWrappers::types::int_type>(index));
      TrilinosScalar value = 0.;

      // If the element is not present on the current processor, we can't
      // continue. This is the main difference to the el() function.
      if (trilinos_i == -1)
        {
#  ifndef DEAL_II_WITH_64BIT_INDICES
          Assert(false,
                 ExcAccessToNonLocalElement(index,
                                            vector->Map().NumMyElements(),
                                            vector->Map().MinMyGID(),
                                            vector->Map().MaxMyGID()));
#  else
          Assert(false,
                 ExcAccessToNonLocalElement(index,
                                            vector->Map().NumMyElements(),
                                            vector->Map().MinMyGID64(),
                                            vector->Map().MaxMyGID64()));
#  endif
        }
      else
        value = (*vector)[0][trilinos_i];

      return value;
    }



    void
    Vector::add(const Vector &v, const bool allow_different_maps)
    {
      if (allow_different_maps == false)
        *this += v;
      else
        {
          Assert(!has_ghost_elements(), ExcGhostsPresent());
          AssertThrow(size() == v.size(),
                      ExcDimensionMismatch(size(), v.size()));

          Epetra_Import data_exchange(vector->Map(), v.vector->Map());
          int           ierr =
            vector->Import(*v.vector, data_exchange, Epetra_AddLocalAlso);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          last_action = Add;
        }
    }



    bool
    Vector::operator==(const Vector &v) const
    {
      Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));
      if (vector->Map().NumMyElements() != v.vector->Map().NumMyElements())
        return false;

      size_type vector_size = vector->Map().NumMyElements();
      for (size_type i = 0; i < vector_size; ++i)
        if ((*(v.vector))[0][i] != (*vector)[0][i])
          return false;

      return true;
    }



    bool
    Vector::operator!=(const Vector &v) const
    {
      Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

      return (!(*this == v));
    }



    bool
    Vector::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      TrilinosScalar       *start_ptr = (*vector)[0];
      const TrilinosScalar *ptr       = start_ptr,
                           *eptr = start_ptr + vector->Map().NumMyElements();
      unsigned int flag          = 0;
      while (ptr != eptr)
        {
          if (*ptr != 0)
            {
              flag = 1;
              break;
            }
          ++ptr;
        }

      // in parallel, check that the vector
      // is zero on _all_ processors.
      const Epetra_MpiComm *mpi_comm =
        dynamic_cast<const Epetra_MpiComm *>(&vector->Map().Comm());
      Assert(mpi_comm != nullptr, ExcInternalError());
      unsigned int num_nonzero = Utilities::MPI::sum(flag, mpi_comm->Comm());
      return num_nonzero == 0;
    }



    bool
    Vector::is_non_negative() const
    {
      // get a representation of the vector and
      // loop over all the elements
      TrilinosScalar       *start_ptr = (*vector)[0];
      const TrilinosScalar *ptr       = start_ptr,
                           *eptr = start_ptr + vector->Map().NumMyElements();
      unsigned int flag          = 0;
      while (ptr != eptr)
        {
          if (*ptr < 0.0)
            {
              flag = 1;
              break;
            }
          ++ptr;
        }

      // in parallel, check that the vector
      // is zero on _all_ processors.
      const auto max_n_negative =
        Utilities::MPI::max(flag, get_mpi_communicator());
      return max_n_negative == 0;
    }



    void
    Vector::print(std::ostream      &out,
                  const unsigned int precision,
                  const bool         scientific,
                  const bool         across) const
    {
      AssertThrow(out.fail() == false, ExcIO());
      boost::io::ios_flags_saver restore_flags(out);


      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      size_type vector_size = vector->Map().NumMyElements();
      if (size() != vector_size)
        {
          auto global_id = [&](const size_type index) {
            return gid(vector->Map(), index);
          };
          out << "size:" << size()
              << " locally_owned_size:" << vector->Map().NumMyElements() << " :"
              << std::endl;
          for (size_type i = 0; i < vector_size; ++i)
            out << "[" << global_id(i) << "]: " << (*(vector))[0][i]
                << std::endl;
        }
      else
        {
          TrilinosScalar *val;
          int             leading_dimension;
          int             ierr = vector->ExtractView(&val, &leading_dimension);

          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          if (across)
            for (size_type i = 0; i < size(); ++i)
              out << static_cast<double>(val[i]) << ' ';
          else
            for (size_type i = 0; i < size(); ++i)
              out << static_cast<double>(val[i]) << std::endl;
          out << std::endl;
        }

      AssertThrow(out.fail() == false, ExcIO());
    }



    void
    Vector::swap(Vector &v) noexcept
    {
      std::swap(last_action, v.last_action);
      std::swap(compressed, v.compressed);
      std::swap(has_ghosts, v.has_ghosts);
      std::swap(vector, v.vector);
      std::swap(nonlocal_vector, v.nonlocal_vector);
      std::swap(owned_elements, v.owned_elements);
    }



    std::size_t
    Vector::memory_consumption() const
    {
      // TODO[TH]: No accurate memory
      // consumption for Trilinos vectors
      // yet. This is a rough approximation with
      // one index and the value per local
      // entry.
      return sizeof(*this) +
             this->vector->Map().NumMyElements() *
               (sizeof(double) + sizeof(TrilinosWrappers::types::int_type));
    }

    // explicit instantiations
#  ifndef DOXYGEN
#    include "lac/trilinos_vector.inst"
#  endif
  } // namespace MPI
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
