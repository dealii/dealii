// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2016 by the deal.II authors
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

#ifndef dealii__la_parallel_vector_templates_h
#define dealii__la_parallel_vector_templates_h


#include <deal.II/base/config.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_operations_internal.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>


DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace distributed
  {

    template <typename Number>
    void
    Vector<Number>::clear_mpi_requests ()
    {
#ifdef DEAL_II_WITH_MPI
      for (size_type j=0; j<compress_requests.size(); j++)
        MPI_Request_free(&compress_requests[j]);
      compress_requests.clear();
      for (size_type j=0; j<update_ghost_values_requests.size(); j++)
        MPI_Request_free(&update_ghost_values_requests[j]);
      update_ghost_values_requests.clear();
#endif
    }



    template <typename Number>
    void
    Vector<Number>::resize_val (const size_type new_alloc_size)
    {
      if (new_alloc_size > allocated_size)
        {
          Assert (((allocated_size > 0 && val != 0) ||
                   val == 0), ExcInternalError());
          if (val != 0)
            free(val);

          Utilities::System::posix_memalign ((void **)&val, 64, sizeof(Number)*new_alloc_size);

          allocated_size = new_alloc_size;
        }
      else if (new_alloc_size == 0)
        {
          if (val != 0)
            free(val);
          val = 0;
          allocated_size = 0;
        }
      thread_loop_partitioner.reset(new ::dealii::parallel::internal::TBBPartitioner());
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const size_type size,
                            const bool      omit_zeroing_entries)
    {
      clear_mpi_requests();

      // check whether we need to reallocate
      resize_val (size);

      // delete previous content in import data
      if (import_data != 0)
        delete[] import_data;
      import_data = 0;

      // set partitioner to serial version
      partitioner.reset (new Utilities::MPI::Partitioner (size));

      // set entries to zero if so requested
      if (omit_zeroing_entries == false)
        this->operator = (Number());

      vector_is_ghosted = false;
    }



    template <typename Number>
    template <typename Number2>
    void
    Vector<Number>::reinit (const Vector<Number2> &v,
                            const bool             omit_zeroing_entries)
    {
      clear_mpi_requests();
      Assert (v.partitioner.get() != 0, ExcNotInitialized());

      // check whether the partitioners are
      // different (check only if the are allocated
      // differently, not if the actual data is
      // different)
      if (partitioner.get() != v.partitioner.get())
        {
          partitioner = v.partitioner;
          const size_type new_allocated_size = partitioner->local_size() +
                                               partitioner->n_ghost_indices();
          resize_val (new_allocated_size);
        }

      if (omit_zeroing_entries == false)
        this->operator= (Number());

      if (import_data != 0)
        {
          delete [] import_data;

          // do not reallocate import_data directly, but only upon request. It
          // is only used as temporary storage for compress() and
          // update_ghost_values, and we might have vectors where we never
          // call these methods and hence do not need to have the storage.
          import_data = 0;
        }

      vector_is_ghosted = false;

      thread_loop_partitioner = v.thread_loop_partitioner;
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const IndexSet &locally_owned_indices,
                            const IndexSet &ghost_indices,
                            const MPI_Comm  communicator)
    {
      // set up parallel partitioner with index sets and communicator
      std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> new_partitioner
      (new Utilities::MPI::Partitioner (locally_owned_indices,
                                        ghost_indices, communicator));
      reinit (new_partitioner);
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const IndexSet &locally_owned_indices,
                            const MPI_Comm  communicator)
    {
      // set up parallel partitioner with index sets and communicator
      std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> new_partitioner
      (new Utilities::MPI::Partitioner (locally_owned_indices,
                                        communicator));
      reinit (new_partitioner);
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_in)
    {
      clear_mpi_requests();
      partitioner = partitioner_in;

      // set vector size and allocate memory
      const size_type new_allocated_size = partitioner->local_size() +
                                           partitioner->n_ghost_indices();
      resize_val (new_allocated_size);

      // initialize to zero
      this->operator= (Number());

      if (import_data != 0)
        {
          delete [] import_data;

          // do not reallocate import_data directly, but only upon request. It
          // is only used as temporary storage for compress() and
          // update_ghost_values, and we might have vectors where we never
          // call these methods and hence do not need to have the storage.
          import_data = 0;
        }

      vector_is_ghosted = false;
    }



    template <typename Number>
    Vector<Number>::Vector ()
      :
      partitioner (new Utilities::MPI::Partitioner()),
      allocated_size (0),
      val (0),
      import_data (0)
    {
      reinit(0);
    }



    template <typename Number>
    Vector<Number>::Vector (const Vector<Number> &v)
      :
      Subscriptor(),
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false)
    {
      reinit (v, true);

      thread_loop_partitioner = v.thread_loop_partitioner;
      dealii::internal::Vector_copy<Number,Number> copier(v.val, val);
      internal::parallel_for(copier, partitioner->local_size(),
                             thread_loop_partitioner);

      zero_out_ghosts();
    }



    template <typename Number>
    Vector<Number>::Vector (const IndexSet &local_range,
                            const IndexSet &ghost_indices,
                            const MPI_Comm  communicator)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false)
    {
      reinit (local_range, ghost_indices, communicator);
    }



    template <typename Number>
    Vector<Number>::Vector (const IndexSet &local_range,
                            const MPI_Comm  communicator)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false)
    {
      reinit (local_range, communicator);
    }



    template <typename Number>
    Vector<Number>::Vector (const size_type size)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false)
    {
      reinit (size, false);
    }



    template <typename Number>
    Vector<Number>::
    Vector (const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false)
    {
      reinit (partitioner);
    }



    template <typename Number>
    inline
    Vector<Number>::~Vector ()
    {
      resize_val(0);

      if (import_data != 0)
        delete[] import_data;
      import_data = 0;

      clear_mpi_requests();
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator = (const Vector<Number> &c)
    {
#ifdef _MSC_VER
      return this->operator=<Number>(c);
#else
      return this->template operator=<Number>(c);
#endif
    }



    template <typename Number>
    template <typename Number2>
    inline
    Vector<Number> &
    Vector<Number>::operator = (const Vector<Number2> &c)
    {
      Assert (c.partitioner.get() != 0, ExcNotInitialized());

      // we update ghost values whenever one of the input or output vector
      // already held ghost values or when we import data from a vector with
      // the same local range but different ghost layout
      bool must_update_ghost_values = c.vector_is_ghosted;

      // check whether the two vectors use the same parallel partitioner. if
      // not, check if all local ranges are the same (that way, we can
      // exchange data between different parallel layouts). One variant which
      // is included here and necessary for compatibility with the other
      // distributed vector classes (Trilinos, PETSc) is the case when vector
      // c does not have any ghosts (constructed without ghost elements given)
      // but the current vector does: In that case, we need to exchange data
      // also when none of the two vector had updated its ghost values before.
      if (partitioner.get() == 0)
        reinit (c, true);
      else if (partitioner.get() != c.partitioner.get())
        {
          // local ranges are also the same if both partitioners are empty
          // (even if they happen to define the empty range as [0,0) or [c,c)
          // for some c!=0 in a different way).
          int local_ranges_are_identical =
            (partitioner->local_range() == c.partitioner->local_range() ||
             (partitioner->local_range().second == partitioner->local_range().first &&
              c.partitioner->local_range().second == c.partitioner->local_range().first));
          if ((c.partitioner->n_mpi_processes() > 1 &&
               Utilities::MPI::min(local_ranges_are_identical,
                                   c.partitioner->get_communicator()) == 0)
              ||
              !local_ranges_are_identical)
            reinit (c, true);
          else
            must_update_ghost_values |= vector_is_ghosted;

          must_update_ghost_values |=
            (c.partitioner->ghost_indices_initialized() == false &&
             partitioner->ghost_indices_initialized() == true);
        }
      else
        must_update_ghost_values |= vector_is_ghosted;

      thread_loop_partitioner = c.thread_loop_partitioner;
      dealii::internal::Vector_copy<Number,Number2> copier(c.val, val);
      internal::parallel_for(copier, partitioner->local_size(),
                             thread_loop_partitioner);

      if (must_update_ghost_values)
        update_ghost_values();
      else
        zero_out_ghosts();
      return *this;
    }



#ifdef DEAL_II_WITH_PETSC

    namespace petsc_helpers
    {
      template <typename PETSC_Number, typename Number>
      void copy_petsc_vector (const PETSC_Number *petsc_start_ptr,
                              const PETSC_Number *petsc_end_ptr,
                              Number *ptr)
      {
        std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
      }

      template <typename PETSC_Number, typename Number>
      void copy_petsc_vector (const std::complex<PETSC_Number> *petsc_start_ptr,
                              const std::complex<PETSC_Number> *petsc_end_ptr,
                              std::complex<Number> *ptr)
      {
        std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
      }

      template <typename PETSC_Number, typename Number>
      void copy_petsc_vector (const std::complex<PETSC_Number> * /*petsc_start_ptr*/,
                              const std::complex<PETSC_Number> * /*petsc_end_ptr*/,
                              Number * /*ptr*/)
      {
        AssertThrow(false, ExcMessage("Tried to copy complex -> real"));
      }
    }

    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator = (const PETScWrappers::MPI::Vector &petsc_vec)
    {
      // TODO: We would like to use the same compact infrastructure as for the
      // Trilinos vector below, but the interface through ReadWriteVector does
      // not support overlapping (ghosted) PETSc vectors, which we need for
      // backward compatibility.

      Assert(petsc_vec.locally_owned_elements() == locally_owned_elements(),
             StandardExceptions::ExcInvalidState());

      // get a representation of the vector and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec &>(petsc_vec), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      const size_type vec_size = local_size();
      petsc_helpers::copy_petsc_vector (start_ptr, start_ptr + vec_size, begin());

      // restore the representation of the vector
      ierr = VecRestoreArray (static_cast<const Vec &>(petsc_vec), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // spread ghost values between processes?
      if (vector_is_ghosted || petsc_vec.has_ghost_elements())
        update_ghost_values();

      // return a reference to this object per normal c++ operator overloading
      // semantics
      return *this;
    }

#endif



#ifdef DEAL_II_WITH_TRILINOS

    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator = (const TrilinosWrappers::MPI::Vector &trilinos_vec)
    {
#ifdef DEAL_II_WITH_MPI
      IndexSet combined_set = partitioner->locally_owned_range();
      combined_set.add_indices(partitioner->ghost_indices());
      ReadWriteVector<Number> rw_vector(combined_set);
      rw_vector.import(trilinos_vec, VectorOperation::insert);
      import(rw_vector, VectorOperation::insert);

      if (vector_is_ghosted || trilinos_vec.has_ghost_elements())
        update_ghost_values();
#else
      AssertThrow(false, ExcNotImplemented());
#endif

      return *this;
    }

#endif



    template <typename Number>
    void
    Vector<Number>::compress (::dealii::VectorOperation::values operation)
    {
      compress_start (0, operation);
      compress_finish(operation);
    }



    template <typename Number>
    void
    Vector<Number>::update_ghost_values () const
    {
      update_ghost_values_start ();
      update_ghost_values_finish ();
    }



    template <typename Number>
    void
    Vector<Number>::zero_out_ghosts ()
    {
      std::fill_n (&val[partitioner->local_size()],
                   partitioner->n_ghost_indices(),
                   Number());
      vector_is_ghosted = false;
    }



    template <typename Number>
    void
    Vector<Number>::compress_start (const unsigned int counter,
                                    ::dealii::VectorOperation::values operation)
    {
      (void)counter;
      (void)operation;
      Assert (vector_is_ghosted == false,
              ExcMessage ("Cannot call compress() on a ghosted vector"));

#ifdef DEAL_II_WITH_MPI
      // nothing to do for insert (only need to zero ghost entries in
      // compress_finish()). in debug mode we want to check consistency
      // of the inserted data, therefore the communication is still
      // initialized. Having different code in debug and optimized mode is
      // somewhat dangerous, but it really saves communication so it seems
      // still worthwhile.
#ifndef DEBUG
      if (operation == VectorOperation::insert)
        return;
#endif

      const Utilities::MPI::Partitioner &part = *partitioner;

      // nothing to do when we neither have import
      // nor ghost indices.
      if (part.n_ghost_indices()==0 && part.n_import_indices()==0)
        return;

      // make this function thread safe
      Threads::Mutex::ScopedLock lock (mutex);

      const unsigned int n_import_targets = part.import_targets().size();
      const unsigned int n_ghost_targets  = part.ghost_targets().size();

      // Need to send and receive the data. Use non-blocking communication,
      // where it is generally less overhead to first initiate the receive and
      // then actually send the data
      if (compress_requests.size() == 0)
        {
          // set channels in different range from update_ghost_values channels
          const unsigned int channel = counter + 400;
          unsigned int current_index_start = 0;
          compress_requests.resize (n_import_targets + n_ghost_targets);

          // allocate import_data in case it is not set up yet
          if (import_data == 0)
            import_data = new Number[part.n_import_indices()];
          for (unsigned int i=0; i<n_import_targets; i++)
            {
              AssertThrow (static_cast<size_type>(part.import_targets()[i].second)*
                           sizeof(Number) <
                           static_cast<size_type>(std::numeric_limits<int>::max()),
                           ExcMessage("Index overflow: Maximum message size in MPI is 2GB. "
                                      "The number of ghost entries times the size of 'Number' "
                                      "exceeds this value. This is not supported."));
              MPI_Recv_init (&import_data[current_index_start],
                             part.import_targets()[i].second*sizeof(Number),
                             MPI_BYTE,
                             part.import_targets()[i].first,
                             part.import_targets()[i].first +
                             part.n_mpi_processes()*channel,
                             part.get_communicator(),
                             &compress_requests[i]);
              current_index_start += part.import_targets()[i].second;
            }
          AssertDimension(current_index_start, part.n_import_indices());

          current_index_start = part.local_size();
          for (unsigned int i=0; i<n_ghost_targets; i++)
            {
              AssertThrow (static_cast<size_type>(part.ghost_targets()[i].second)*
                           sizeof(Number) <
                           static_cast<size_type>(std::numeric_limits<int>::max()),
                           ExcMessage("Index overflow: Maximum message size in MPI is 2GB. "
                                      "The number of ghost entries times the size of 'Number' "
                                      "exceeds this value. This is not supported."));
              MPI_Send_init (&this->val[current_index_start],
                             part.ghost_targets()[i].second*sizeof(Number),
                             MPI_BYTE,
                             part.ghost_targets()[i].first,
                             part.this_mpi_process() +
                             part.n_mpi_processes()*channel,
                             part.get_communicator(),
                             &compress_requests[n_import_targets+i]);
              current_index_start += part.ghost_targets()[i].second;
            }
          AssertDimension (current_index_start,
                           part.local_size()+part.n_ghost_indices());
        }

      AssertDimension(n_import_targets + n_ghost_targets,
                      compress_requests.size());
      if (compress_requests.size() > 0)
        {
          int ierr = MPI_Startall(compress_requests.size(),&compress_requests[0]);
          (void)ierr;
          Assert (ierr == MPI_SUCCESS, ExcInternalError());
        }
#endif
    }



    template <typename Number>
    void
    Vector<Number>::compress_finish (::dealii::VectorOperation::values operation)
    {
#ifdef DEAL_II_WITH_MPI

      // in optimized mode, no communication was started, so leave the
      // function directly (and only clear ghosts)
#ifndef DEBUG
      if (operation == VectorOperation::insert)
        {
          zero_out_ghosts();
          return;
        }
#endif

      const Utilities::MPI::Partitioner &part = *partitioner;

      // nothing to do when we neither have import nor ghost indices.
      if (part.n_ghost_indices()==0 && part.n_import_indices()==0)
        return;

      // make this function thread safe
      Threads::Mutex::ScopedLock lock (mutex);

      const unsigned int n_import_targets = part.import_targets().size();
      const unsigned int n_ghost_targets  = part.ghost_targets().size();

      if (operation != dealii::VectorOperation::insert)
        AssertDimension (n_ghost_targets+n_import_targets,
                         compress_requests.size());

      // first wait for the receive to complete
      if (compress_requests.size() > 0 && n_import_targets > 0)
        {
          int ierr = MPI_Waitall (n_import_targets, &compress_requests[0],
                                  MPI_STATUSES_IGNORE);
          (void)ierr;
          Assert (ierr == MPI_SUCCESS, ExcInternalError());

          Number *read_position = import_data;
          std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
          my_imports = part.import_indices().begin();

          // If the operation is no insertion, add the imported data to the
          // local values. For insert, nothing is done here (but in debug mode
          // we assert that the specified value is either zero or matches with
          // the ones already present
          if (operation != dealii::VectorOperation::insert)
            for ( ; my_imports!=part.import_indices().end(); ++my_imports)
              for (unsigned int j=my_imports->first; j<my_imports->second; j++)
                local_element(j) += *read_position++;
          else
            for ( ; my_imports!=part.import_indices().end(); ++my_imports)
              for (unsigned int j=my_imports->first; j<my_imports->second;
                   j++, read_position++)
                Assert(*read_position == Number() ||
                       std::abs(local_element(j) - *read_position) <=
                       std::abs(local_element(j)) * 1000. *
                       std::numeric_limits<real_type>::epsilon(),
                       ExcNonMatchingElements(*read_position, local_element(j),
                                              part.this_mpi_process()));
          AssertDimension(read_position-import_data,part.n_import_indices());
        }

      if (compress_requests.size() > 0 && n_ghost_targets > 0)
        {
          int ierr = MPI_Waitall (n_ghost_targets,
                                  &compress_requests[n_import_targets],
                                  MPI_STATUSES_IGNORE);
          (void)ierr;
          Assert (ierr == MPI_SUCCESS, ExcInternalError());
        }
      else
        AssertDimension (part.n_ghost_indices(), 0);

      zero_out_ghosts ();
#else
      (void)operation;
#endif
    }



    template <typename Number>
    void
    Vector<Number>::update_ghost_values_start (const unsigned int counter) const
    {
#ifdef DEAL_II_WITH_MPI
      const Utilities::MPI::Partitioner &part = *partitioner;

      // nothing to do when we neither have import nor ghost indices.
      if (part.n_ghost_indices()==0 && part.n_import_indices()==0)
        return;

      // make this function thread safe
      Threads::Mutex::ScopedLock lock (mutex);

      const unsigned int n_import_targets = part.import_targets().size();
      const unsigned int n_ghost_targets = part.ghost_targets().size();

      // Need to send and receive the data. Use non-blocking communication,
      // where it is generally less overhead to first initiate the receive and
      // then actually send the data
      if (update_ghost_values_requests.size() == 0)
        {
          size_type current_index_start = part.local_size();
          update_ghost_values_requests.resize (n_import_targets+n_ghost_targets);
          for (unsigned int i=0; i<n_ghost_targets; i++)
            {
              // allow writing into ghost indices even though we are in a
              // const function
              MPI_Recv_init (const_cast<Number *>(&val[current_index_start]),
                             part.ghost_targets()[i].second*sizeof(Number),
                             MPI_BYTE,
                             part.ghost_targets()[i].first,
                             part.ghost_targets()[i].first +
                             counter*part.n_mpi_processes(),
                             part.get_communicator(),
                             &update_ghost_values_requests[i]);
              current_index_start += part.ghost_targets()[i].second;
            }
          AssertDimension (current_index_start,
                           part.local_size()+part.n_ghost_indices());

          // allocate import_data in case it is not set up yet
          if (import_data == 0 && part.n_import_indices() > 0)
            import_data = new Number[part.n_import_indices()];
          current_index_start = 0;
          for (unsigned int i=0; i<n_import_targets; i++)
            {
              MPI_Send_init (&import_data[current_index_start],
                             part.import_targets()[i].second*sizeof(Number),
                             MPI_BYTE, part.import_targets()[i].first,
                             part.this_mpi_process() +
                             part.n_mpi_processes()*counter,
                             part.get_communicator(),
                             &update_ghost_values_requests[n_ghost_targets+i]);
              current_index_start += part.import_targets()[i].second;
            }
          AssertDimension (current_index_start, part.n_import_indices());
        }

      // copy the data that is actually to be send to the import_data field
      if (part.n_import_indices() > 0)
        {
          Assert (import_data != 0, ExcInternalError());
          Number *write_position = import_data;
          std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
          my_imports = part.import_indices().begin();
          for ( ; my_imports!=part.import_indices().end(); ++my_imports)
            for (unsigned int j=my_imports->first; j<my_imports->second; j++)
              *write_position++ = local_element(j);
        }

      AssertDimension (n_import_targets+n_ghost_targets,
                       update_ghost_values_requests.size());
      if (update_ghost_values_requests.size() > 0)
        {
          int ierr = MPI_Startall(update_ghost_values_requests.size(),
                                  &update_ghost_values_requests[0]);
          (void)ierr;
          Assert (ierr == MPI_SUCCESS, ExcInternalError());
        }
#else
      (void)counter;
#endif
    }



    template <typename Number>
    void
    Vector<Number>::update_ghost_values_finish () const
    {
#ifdef DEAL_II_WITH_MPI
      // wait for both sends and receives to complete, even though only
      // receives are really necessary. this gives (much) better performance
      AssertDimension (partitioner->ghost_targets().size() +
                       partitioner->import_targets().size(),
                       update_ghost_values_requests.size());
      if (update_ghost_values_requests.size() > 0)
        {
          // make this function thread safe
          Threads::Mutex::ScopedLock lock (mutex);

          int ierr = MPI_Waitall (update_ghost_values_requests.size(),
                                  &update_ghost_values_requests[0],
                                  MPI_STATUSES_IGNORE);
          (void)ierr;
          Assert (ierr == MPI_SUCCESS, ExcInternalError());
        }
#endif
      vector_is_ghosted = true;
    }



    template <typename Number>
    void
    Vector<Number>::import(const ReadWriteVector<Number>                  &V,
                           VectorOperation::values                         operation,
                           std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern)
    {
      // If no communication pattern is given, create one. Otherwise, use the
      // given one.
      std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> comm_pattern;
      if (communication_pattern.get() == NULL)
        {
          comm_pattern.reset(new Utilities::MPI::Partitioner(locally_owned_elements(),
                                                             V.get_stored_elements(),
                                                             get_mpi_communicator()));
        }
      else
        {
          comm_pattern =
            std_cxx11::dynamic_pointer_cast<const Utilities::MPI::Partitioner> (communication_pattern);
          AssertThrow(comm_pattern != NULL,
                      ExcMessage("The communication pattern is not of type "
                                 "Utilities::MPI::Partitioner."));
        }
      Vector<Number> tmp_vector(comm_pattern);

      // fill entries from ReadWriteVector into the distributed vector,
      // including ghost entries. this is not really efficient right now
      // because indices are translated twice, once by nth_index_in_set(i) and
      // once for operator() of tmp_vector
      const IndexSet &v_stored = V.get_stored_elements();
      for (size_type i=0; i<v_stored.n_elements(); ++i)
        tmp_vector(v_stored.nth_index_in_set(i)) = V.local_element(i);

      tmp_vector.compress(operation);

      dealii::internal::Vector_copy<Number,Number> copier(tmp_vector.val, val);
      internal::parallel_for(copier, partitioner->local_size(),
                             thread_loop_partitioner);
    }



    template <typename Number>
    void
    Vector<Number>::swap (Vector<Number> &v)
    {
#ifdef DEAL_II_WITH_MPI

#ifdef DEBUG
      if (Utilities::MPI::job_supports_mpi())
        {
          // make sure that there are not outstanding requests from updating
          // ghost values or compress
          int flag = 1;
          if (update_ghost_values_requests.size()>0)
            {
              int ierr = MPI_Testall (update_ghost_values_requests.size(),
                                      &update_ghost_values_requests[0],
                                      &flag, MPI_STATUSES_IGNORE);
              Assert (ierr == MPI_SUCCESS, ExcInternalError());
              Assert (flag == 1,
                      ExcMessage("MPI found unfinished update_ghost_values() requests"
                                 "when calling swap, which is not allowed"));
            }
          if (compress_requests.size()>0)
            {
              int ierr = MPI_Testall (compress_requests.size(), &compress_requests[0],
                                      &flag, MPI_STATUSES_IGNORE);
              Assert (ierr == MPI_SUCCESS, ExcInternalError());
              Assert (flag == 1,
                      ExcMessage("MPI found unfinished compress() requests "
                                 "when calling swap, which is not allowed"));
            }
        }
#endif

      std::swap (compress_requests, v.compress_requests);
      std::swap (update_ghost_values_requests, v.update_ghost_values_requests);
#endif

      std::swap (partitioner,       v.partitioner);
      std::swap (thread_loop_partitioner, v.thread_loop_partitioner);
      std::swap (allocated_size,    v.allocated_size);
      std::swap (val,               v.val);
      std::swap (import_data,       v.import_data);
      std::swap (vector_is_ghosted, v.vector_is_ghosted);
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator = (const Number s)
    {
      internal::Vector_set<Number> setter(s, val);

      internal::parallel_for(setter, partitioner->local_size(),
                             thread_loop_partitioner);

      // if we call Vector::operator=0, we want to zero out all the entries
      // plus ghosts.
      if (s==Number())
        zero_out_ghosts();

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator += (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_add_v<Number> vector_add(val, v.val);
      internal::parallel_for(vector_add, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator -= (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_subtract_v<Number> vector_subtract(val, v.val);
      internal::parallel_for(vector_subtract, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number>
    void
    Vector<Number>::add (const Number a)
    {
      AssertIsFinite(a);

      internal::Vectorization_add_factor<Number> vector_add(val, a);
      internal::parallel_for(vector_add, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::add (const Number a,
                         const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertIsFinite(a);
      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_add_av<Number> vector_add(val, v.val, a);
      internal::parallel_for(vector_add, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::add (const Number a,
                         const VectorSpaceVector<Number> &vv,
                         const Number b,
                         const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);
      Assert(dynamic_cast<const Vector<Number> *>(&ww)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &w = dynamic_cast<const Vector<Number> &>(ww);

      AssertIsFinite(a);
      AssertIsFinite(b);

      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());

      internal::Vectorization_add_avpbw<Number> vector_add(val, v.val, w.val, a, b);
      internal::parallel_for(vector_add, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::sadd (const Number x,
                          const Vector<Number> &v)
    {
      AssertIsFinite(x);
      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_sadd_xv<Number> vector_sadd(val, v.val, x);
      internal::parallel_for(vector_sadd, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::sadd (const Number x,
                          const Number a,
                          const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertIsFinite(x);
      AssertIsFinite(a);
      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_sadd_xav<Number> vector_sadd(val, v.val, a, x);
      internal::parallel_for(vector_sadd, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::sadd (const Number x,
                          const Number a,
                          const Vector<Number> &v,
                          const Number b,
                          const Vector<Number> &w)
    {
      AssertIsFinite(x);
      AssertIsFinite(a);
      AssertIsFinite(b);

      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());

      internal::Vectorization_sadd_xavbw<Number> vector_sadd(val, v.val, w.val,
                                                             x, a, b);
      internal::parallel_for(vector_sadd, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator *= (const Number factor)
    {
      AssertIsFinite(factor);
      internal::Vectorization_multiply_factor<Number> vector_multiply(val, factor);

      internal::parallel_for(vector_multiply, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator /= (const Number factor)
    {
      operator *= (static_cast<Number>(1.)/factor);
      return *this;
    }



    template <typename Number>
    void
    Vector<Number>::scale (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_scale<Number> vector_scale(val, v.val);
      internal::parallel_for(vector_scale, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::equ (const Number a,
                         const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertIsFinite(a);
      AssertDimension (local_size(), v.local_size());

      internal::Vectorization_equ_au<Number> vector_equ(val, v.val, a);
      internal::parallel_for(vector_equ, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number> &v,
                         const Number b,
                         const Vector<Number> &w)
    {
      AssertIsFinite(a);
      AssertIsFinite(b);

      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());

      internal::Vectorization_equ_aubv<Number> vector_equ(val, v.val, w.val, a, b);
      internal::parallel_for(vector_equ, partitioner->local_size(),
                             thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    bool
    Vector<Number>::all_zero_local () const
    {
      const size_type local_size = partitioner->local_size();
      for (size_type i=0; i<local_size; ++i)
        if (val[i] != Number(0))
          return false;
      return true;
    }



    template <typename Number>
    bool
    Vector<Number>::all_zero () const
    {
      // use int instead of bool. in order to make global reduction operations
      // work also when MPI_Init was not called, only call MPI_Allreduce
      // commands when there is more than one processor (note that reinit()
      // functions handle this case correctly through the job_supports_mpi()
      // query). this is the same in all the reduce functions below
      int local_result = -static_cast<int>(all_zero_local());
      if (partitioner->n_mpi_processes() > 1)
        return -Utilities::MPI::max(local_result,
                                    partitioner->get_communicator());
      else
        return -local_result;
    }



    template <typename Number>
    template <typename Number2>
    Number
    Vector<Number>::inner_product_local(const Vector<Number2> &v) const
    {
      if (PointerComparison::equal (this, &v))
        return norm_sqr_local();

      AssertDimension (partitioner->local_size(), v.partitioner->local_size());

      Number sum;
      internal::Dot<Number,Number2> dot(val, v.val);
      internal::parallel_reduce (dot, partitioner->local_size(), sum,
                                 thread_loop_partitioner);
      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number>
    Number
    Vector<Number>::operator * (const VectorSpaceVector<Number> &vv) const
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      Number local_result = inner_product_local(v);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::norm_sqr_local () const
    {
      real_type sum;
      internal::Norm2<Number,real_type> norm2(val);
      internal::parallel_reduce (norm2, partitioner->local_size(), sum,
                                 thread_loop_partitioner);
      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number>
    Number
    Vector<Number>::mean_value_local () const
    {
      Assert(size() != 0, ExcEmptyObject());

      if (partitioner->local_size() == 0)
        return Number();

      Number sum;
      internal::MeanValue<Number> mean(val);
      internal::parallel_reduce (mean, partitioner->local_size(), sum,
                                 thread_loop_partitioner);

      return sum / real_type(partitioner->local_size());
    }



    template <typename Number>
    Number
    Vector<Number>::mean_value () const
    {
      Number local_result = mean_value_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result *
                                    (real_type)partitioner->local_size(),
                                    partitioner->get_communicator())
               /(real_type)partitioner->size();
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm_local () const
    {
      real_type sum;
      internal::Norm1<Number, real_type> norm1(val);
      internal::parallel_reduce (norm1, partitioner->local_size(), sum,
                                 thread_loop_partitioner);

      return sum;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm () const
    {
      real_type local_result = l1_norm_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l2_norm () const
    {
      real_type local_result = norm_sqr_local();
      if (partitioner->n_mpi_processes() > 1)
        return std::sqrt(Utilities::MPI::sum(local_result,
                                             partitioner->get_communicator()));
      else
        return std::sqrt(local_result);
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::lp_norm_local (const real_type p) const
    {
      real_type sum;
      internal::NormP<Number, real_type> normp(val, p);
      internal::parallel_reduce (normp, partitioner->local_size(), sum,
                                 thread_loop_partitioner);
      return std::pow(sum, 1./p);
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::lp_norm (const real_type p) const
    {
      const real_type local_result = lp_norm_local(p);
      if (partitioner->n_mpi_processes() > 1)
        return std::pow (Utilities::MPI::sum(std::pow(local_result,p),
                                             partitioner->get_communicator()),
                         static_cast<real_type>(1.0/p));
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm_local () const
    {
      real_type max = 0.;

      const size_type local_size = partitioner->local_size();
      for (size_type i=0; i<local_size; ++i)
        max = std::max (numbers::NumberTraits<Number>::abs(val[i]), max);

      return max;
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm () const
    {
      const real_type local_result = linfty_norm_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::max (local_result,
                                    partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    Number
    Vector<Number>::add_and_dot_local(const Number          a,
                                      const Vector<Number> &v,
                                      const Vector<Number> &w)
    {
      const size_type vec_size = partitioner->local_size();
      AssertDimension (vec_size, v.local_size());
      AssertDimension (vec_size, w.local_size());

      Number sum;
      internal::AddAndDot<Number> adder(this->val, v.val, w.val, a);
      internal::parallel_reduce (adder, vec_size, sum, thread_loop_partitioner);
      AssertIsFinite(sum);
      return sum;
    }



    template <typename Number>
    Number
    Vector<Number>::add_and_dot (const Number                     a,
                                 const VectorSpaceVector<Number> &vv,
                                 const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);
      Assert(dynamic_cast<const Vector<Number> *>(&ww)!=NULL,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &w = dynamic_cast<const Vector<Number> &>(ww);

      Number local_result = add_and_dot_local(a, v, w);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::partitioners_are_compatible
    (const Utilities::MPI::Partitioner &part) const
    {
      return partitioner->is_compatible (part);
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::partitioners_are_globally_compatible
    (const Utilities::MPI::Partitioner &part) const
    {
      return partitioner->is_globally_compatible (part);
    }



    template <typename Number>
    std::size_t
    Vector<Number>::memory_consumption () const
    {
      std::size_t memory = sizeof(*this);
      memory += sizeof (Number) * static_cast<std::size_t>(allocated_size);

      // if the partitioner is shared between more processors, just count a
      // fraction of that memory, since we're not actually using more memory
      // for it.
      if (partitioner.use_count() > 0)
        memory += partitioner->memory_consumption()/partitioner.use_count()+1;
      if (import_data != 0)
        memory += (static_cast<std::size_t>(partitioner->n_import_indices())*
                   sizeof(Number));
      return memory;
    }



    template <typename Number>
    void
    Vector<Number>::print (std::ostream      &out,
                           const unsigned int precision,
                           const bool         scientific,
                           const bool         across) const
    {
      Assert (partitioner.get() !=0, ExcInternalError());
      AssertThrow (out, ExcIO());
      std::ios::fmtflags old_flags = out.flags();
      unsigned int old_precision = out.precision (precision);

      out.precision (precision);
      if (scientific)
        out.setf (std::ios::scientific, std::ios::floatfield);
      else
        out.setf (std::ios::fixed, std::ios::floatfield);

      // to make the vector write out all the information in order, use as
      // many barriers as there are processors and start writing when it's our
      // turn
#ifdef DEAL_II_WITH_MPI
      if (partitioner->n_mpi_processes() > 1)
        for (unsigned int i=0; i<partitioner->this_mpi_process(); i++)
          MPI_Barrier (partitioner->get_communicator());
#endif

      out << "Process #" << partitioner->this_mpi_process() << std::endl
          << "Local range: [" << partitioner->local_range().first << ", "
          << partitioner->local_range().second << "), global size: "
          << partitioner->size() << std::endl
          << "Vector data:" << std::endl;
      if (across)
        for (size_type i=0; i<partitioner->local_size(); ++i)
          out << local_element(i) << ' ';
      else
        for (size_type i=0; i<partitioner->local_size(); ++i)
          out << local_element(i) << std::endl;
      out << std::endl;

      if (vector_is_ghosted)
        {
          out << "Ghost entries (global index / value):" << std::endl;
          if (across)
            for (size_type i=0; i<partitioner->n_ghost_indices(); ++i)
              out << '(' << partitioner->ghost_indices().nth_index_in_set(i)
                  << '/' << local_element(partitioner->local_size()+i) << ") ";
          else
            for (size_type i=0; i<partitioner->n_ghost_indices(); ++i)
              out << '(' << partitioner->ghost_indices().nth_index_in_set(i)
                  << '/' << local_element(partitioner->local_size()+i) << ")"
                  << std::endl;
          out << std::endl;
        }
      out << std::flush;

#ifdef DEAL_II_WITH_MPI
      if (partitioner->n_mpi_processes() > 1)
        {
          MPI_Barrier (partitioner->get_communicator());

          for (unsigned int i=partitioner->this_mpi_process()+1;
               i<partitioner->n_mpi_processes(); i++)
            MPI_Barrier (partitioner->get_communicator());
        }
#endif

      AssertThrow (out, ExcIO());
      // reset output format
      out.flags (old_flags);
      out.precision(old_precision);
    }

  } // end of namespace distributed

} // end of namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
