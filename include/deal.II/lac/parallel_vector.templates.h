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

#ifndef dealii__parallel_vector_templates_h
#define dealii__parallel_vector_templates_h


#include <deal.II/base/config.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/vector_view.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>

DEAL_II_NAMESPACE_OPEN


namespace parallel
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
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const size_type size,
                            const bool      omit_zeroing_entries)
    {
      clear_mpi_requests();
      // check whether we need to reallocate
      resize_val (size);

      // reset vector view
      vector_view.reinit (size, val);

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
          vector_view.reinit (partitioner->local_size(), val);
        }
      else
        Assert (vector_view.size() == partitioner->local_size(),
                ExcInternalError());

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
      vector_view.reinit (partitioner->local_size(), val);

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



#ifdef DEAL_II_WITH_PETSC

    namespace internal
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
      Assert(petsc_vec.locally_owned_elements() == locally_owned_elements(),
             StandardExceptions::ExcInvalidState());

      // get a representation of the vector and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec &>(petsc_vec), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      const size_type vec_size = local_size();
      internal::copy_petsc_vector (start_ptr, start_ptr + vec_size, begin());

      // restore the representation of the vector
      ierr = VecRestoreArray (static_cast<const Vec &>(petsc_vec), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // spread ghost values between processes?
      if (vector_is_ghosted || petsc_vec.has_ghost_elements())
        update_ghost_values();

      // return a pointer to this object per normal c++ operator overloading
      // semantics
      return *this;
    }

#endif



#ifdef DEAL_II_WITH_TRILINOS

    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator = (const TrilinosWrappers::MPI::Vector &trilinos_vec)
    {
      if (trilinos_vec.has_ghost_elements() == false)
        {
          Assert(trilinos_vec.locally_owned_elements() == locally_owned_elements(),
                 StandardExceptions::ExcInvalidState());
        }
      else
        // ghosted trilinos vector must contain the local range of this vector
        // which is contiguous
        {
          Assert((trilinos_vec.locally_owned_elements() & locally_owned_elements())
                 == locally_owned_elements(),
                 StandardExceptions::ExcInvalidState());
        }

      // create on trilinos data
      const std::size_t start_index =
        trilinos_vec.vector_partitioner().NumMyElements() > 0 ?
        trilinos_vec.vector_partitioner().
        LID(static_cast<TrilinosWrappers::types::int_type>(this->local_range().first)) : 0;
      const VectorView<double> in_view (local_size(), trilinos_vec.begin()+start_index);
      static_cast<dealii::Vector<Number>&>(vector_view) =
        static_cast<const dealii::Vector<double>&>(in_view);

      // spread ghost values between processes?
      if (vector_is_ghosted || trilinos_vec.has_ghost_elements())
        update_ghost_values();

      // return a pointer to this object per normal c++ operator overloading
      // semantics
      return *this;
    }

#endif



    template <typename Number>
    void
    Vector<Number>::copy_from (const Vector<Number> &c,
                               const bool            call_update_ghost_values)
    {
      AssertDimension (local_range().first, c.local_range().first);
      AssertDimension (local_range().second, c.local_range().second);
      AssertDimension (vector_view.size(), c.vector_view.size());
      vector_view = c.vector_view;
      if (call_update_ghost_values == true)
        update_ghost_values();
      else
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

          Assert (part.local_size() == vector_view.size(), ExcInternalError());
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
                Assert(*read_position == 0. ||
                       std::abs(local_element(j) - *read_position) <=
                       std::abs(local_element(j)) * 1000. *
                       std::numeric_limits<Number>::epsilon(),
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
          Assert (part.local_size() == vector_view.size(),
                  ExcInternalError());
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
      std::swap (allocated_size,    v.allocated_size);
      std::swap (val,               v.val);
      std::swap (import_data,       v.import_data);
      std::swap (vector_is_ghosted, v.vector_is_ghosted);

      // vector view cannot be swapped so reset it manually (without touching
      // the vector elements)
      vector_view.reinit (partitioner->local_size(), val);
      v.vector_view.reinit (v.partitioner->local_size(), v.val);
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
