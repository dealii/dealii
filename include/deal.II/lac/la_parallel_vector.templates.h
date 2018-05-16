// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

#ifndef dealii_la_parallel_vector_templates_h
#define dealii_la_parallel_vector_templates_h


#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx14/memory.h>
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
        {
          const int ierr = MPI_Request_free(&compress_requests[j]);
          AssertThrowMPI(ierr);
        }
      compress_requests.clear();
      for (size_type j=0; j<update_ghost_values_requests.size(); j++)
        {
          const int ierr = MPI_Request_free(&update_ghost_values_requests[j]);
          AssertThrowMPI(ierr);
        }
      update_ghost_values_requests.clear();
#endif
    }



    template <typename Number>
    void
    Vector<Number>::resize_val (const size_type new_alloc_size)
    {
      if (new_alloc_size > allocated_size)
        {
          Assert (((allocated_size > 0 && values != nullptr) ||
                   values == nullptr), ExcInternalError());

          Number *new_val;
          Utilities::System::posix_memalign ((void **)&new_val, 64, sizeof(Number)*new_alloc_size);
          values.reset (new_val);

          allocated_size = new_alloc_size;
        }
      else if (new_alloc_size == 0)
        {
          values.reset();
          allocated_size = 0;
        }
      thread_loop_partitioner = std::make_shared<::dealii::parallel::internal::TBBPartitioner>();
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
      import_data.reset ();

      // set partitioner to serial version
      partitioner = std::make_shared<Utilities::MPI::Partitioner> (size);

      // set entries to zero if so requested
      if (omit_zeroing_entries == false)
        this->operator = (Number());
      else
        zero_out_ghosts();
    }



    template <typename Number>
    template <typename Number2>
    void
    Vector<Number>::reinit (const Vector<Number2> &v,
                            const bool             omit_zeroing_entries)
    {
      clear_mpi_requests();
      Assert (v.partitioner.get() != nullptr, ExcNotInitialized());

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
      else
        zero_out_ghosts();

      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.reset ();

      thread_loop_partitioner = v.thread_loop_partitioner;
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const IndexSet &locally_owned_indices,
                            const IndexSet &ghost_indices,
                            const MPI_Comm  communicator)
    {
      // set up parallel partitioner with index sets and communicator
      std::shared_ptr<const Utilities::MPI::Partitioner> new_partitioner
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
      std::shared_ptr<const Utilities::MPI::Partitioner> new_partitioner
      (new Utilities::MPI::Partitioner (locally_owned_indices,
                                        communicator));
      reinit (new_partitioner);
    }



    template <typename Number>
    void
    Vector<Number>::reinit (const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_in)
    {
      clear_mpi_requests();
      partitioner = partitioner_in;

      // set vector size and allocate memory
      const size_type new_allocated_size = partitioner->local_size() +
                                           partitioner->n_ghost_indices();
      resize_val (new_allocated_size);

      // initialize to zero
      this->operator= (Number());


      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.reset ();

      vector_is_ghosted = false;
    }



    template <typename Number>
    Vector<Number>::Vector ()
      :
      partitioner (new Utilities::MPI::Partitioner()),
      allocated_size (0),
      values (nullptr, &free)
    {
      reinit(0);
    }



    template <typename Number>
    Vector<Number>::Vector (const Vector<Number> &v)
      :
      Subscriptor(),
      allocated_size (0),
      values (nullptr, &free),
      vector_is_ghosted (false)
    {
      reinit (v, true);

      thread_loop_partitioner = v.thread_loop_partitioner;

      const size_type this_size = local_size();
      if (this_size>0)
        {
          dealii::internal::VectorOperations::Vector_copy<Number,Number> copier(v.values.get(), values.get());
          internal::VectorOperations::parallel_for(copier, 0, partitioner->local_size(),
                                                   thread_loop_partitioner);
        }
    }



    template <typename Number>
    Vector<Number>::Vector (const IndexSet &local_range,
                            const IndexSet &ghost_indices,
                            const MPI_Comm  communicator)
      :
      allocated_size (0),
      values (nullptr, &free),
      vector_is_ghosted (false)
    {
      reinit (local_range, ghost_indices, communicator);
    }



    template <typename Number>
    Vector<Number>::Vector (const IndexSet &local_range,
                            const MPI_Comm  communicator)
      :
      allocated_size (0),
      values (nullptr, &free),
      vector_is_ghosted (false)
    {
      reinit (local_range, communicator);
    }



    template <typename Number>
    Vector<Number>::Vector (const size_type size)
      :
      allocated_size (0),
      values (nullptr, &free),
      vector_is_ghosted (false)
    {
      reinit (size, false);
    }



    template <typename Number>
    Vector<Number>::
    Vector (const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      :
      allocated_size (0),
      values (nullptr, &free),
      vector_is_ghosted (false)
    {
      reinit (partitioner);
    }



    template <typename Number>
    inline
    Vector<Number>::~Vector ()
    {
      try
        {
          clear_mpi_requests();
        }
      catch (...)
        {}
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
      Assert (c.partitioner.get() != nullptr, ExcNotInitialized());

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
      if (partitioner.get() == nullptr)
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
                                   c.partitioner->get_mpi_communicator()) == 0)
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

      const size_type this_size = partitioner->local_size();
      if (this_size>0)
        {
          dealii::internal::VectorOperations::Vector_copy<Number,Number2> copier(c.values.get(), values.get());
          internal::VectorOperations::parallel_for(copier, 0, this_size,
                                                   thread_loop_partitioner);
        }

      if (must_update_ghost_values)
        update_ghost_values();
      else
        zero_out_ghosts();
      return *this;
    }



    template <typename Number>
    template <typename Number2>
    void
    Vector<Number>::copy_locally_owned_data_from(const Vector<Number2> &src)
    {
      AssertDimension(partitioner->local_size(), src.partitioner->local_size());
      if (partitioner->local_size() > 0)
        {
          dealii::internal::VectorOperations::Vector_copy<Number,Number2> copier(src.values.get(),
              values.get());
          internal::VectorOperations::parallel_for(copier, 0, partitioner->local_size(),
                                                   thread_loop_partitioner);
        }
    }




#ifdef DEAL_II_WITH_PETSC

    namespace petsc_helpers
    {
      template <typename PETSC_Number, typename Number>
      void
      copy_petsc_vector (const PETSC_Number *petsc_start_ptr,
                         const PETSC_Number *petsc_end_ptr,
                         Number *ptr)
      {
        std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
      }

      template <typename PETSC_Number, typename Number>
      void
      copy_petsc_vector (const std::complex<PETSC_Number> *petsc_start_ptr,
                         const std::complex<PETSC_Number> *petsc_end_ptr,
                         std::complex<Number> *ptr)
      {
        std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
      }

      template <typename PETSC_Number, typename Number>
      void
      copy_petsc_vector (const std::complex<PETSC_Number> * /*petsc_start_ptr*/,
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
      PetscErrorCode ierr = VecGetArray (static_cast<const Vec &>(petsc_vec),
                                         &start_ptr);
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
    Vector<Number>::zero_out_ghosts () const
    {
      if (values != nullptr)
        std::fill_n (values.get()+partitioner->local_size(),
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
      // make this function thread safe
      Threads::Mutex::ScopedLock lock (mutex);

      // allocate import_data in case it is not set up yet
      if (import_data == nullptr && partitioner->n_import_indices() > 0)
        import_data = std_cxx14::make_unique<Number[]>(partitioner->n_import_indices());

      partitioner->import_from_ghosted_array_start
      (operation, counter,
       ArrayView<Number>(values.get() + partitioner->local_size(),partitioner->n_ghost_indices()),
       ArrayView<Number>(import_data.get(), partitioner->n_import_indices()),
       compress_requests);
#endif
    }



    template <typename Number>
    void
    Vector<Number>::compress_finish (::dealii::VectorOperation::values operation)
    {
#ifdef DEAL_II_WITH_MPI
      vector_is_ghosted = false;

      // in order to zero ghost part of the vector, we need to call
      // import_from_ghosted_array_finish() regardless of
      // compress_requests.size() == 0

      // make this function thread safe
      Threads::Mutex::ScopedLock lock (mutex);

      Assert(partitioner->n_import_indices() == 0 ||
             import_data != nullptr,
             ExcNotInitialized());
      partitioner->import_from_ghosted_array_finish
      (operation,
       ArrayView<const Number>(import_data.get(), partitioner->n_import_indices()),
       ArrayView<Number>(values.get(), partitioner->local_size()),
       ArrayView<Number>(values.get() + partitioner->local_size(),partitioner->n_ghost_indices()),
       compress_requests);
#else
      (void)operation;
#endif
    }



    template <typename Number>
    void
    Vector<Number>::update_ghost_values_start (const unsigned int counter) const
    {
#ifdef DEAL_II_WITH_MPI
      // nothing to do when we neither have import nor ghost indices.
      if (partitioner->n_ghost_indices()==0 && partitioner->n_import_indices()==0)
        return;

      // make this function thread safe
      Threads::Mutex::ScopedLock lock (mutex);

      // allocate import_data in case it is not set up yet
      if (import_data == nullptr && partitioner->n_import_indices() > 0)
        import_data = std_cxx14::make_unique<Number[]>(partitioner->n_import_indices());

      partitioner->export_to_ghosted_array_start
      (counter,
       ArrayView<const Number>(values.get(), partitioner->local_size()),
       ArrayView<Number>(import_data.get(), partitioner->n_import_indices()),
       ArrayView<Number>(values.get() + partitioner->local_size(),partitioner->n_ghost_indices()),
       update_ghost_values_requests);

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

          partitioner->export_to_ghosted_array_finish
          (ArrayView<Number>(values.get() + partitioner->local_size(),partitioner->n_ghost_indices()),
           update_ghost_values_requests);
        }
#endif
      vector_is_ghosted = true;
    }



    template <typename Number>
    void
    Vector<Number>::import(const ReadWriteVector<Number>                  &V,
                           VectorOperation::values                         operation,
                           std::shared_ptr<const CommunicationPatternBase> communication_pattern)
    {
      IndexSet locally_owned_elem = locally_owned_elements();
      // If no communication pattern is given, create one. Otherwise, use the
      // given one.
      std::shared_ptr<const Utilities::MPI::Partitioner> comm_pattern;
      if (communication_pattern.get() == nullptr)
        {
          // Split the IndexSet of V in locally owned elements and ghost indices
          // then create the communication pattern
          IndexSet ghost_indices(V.get_stored_elements());
          ghost_indices.subtract_set(locally_owned_elem);
          IndexSet local_indices(V.get_stored_elements());
          local_indices.subtract_set(ghost_indices);
          comm_pattern = std::make_shared<Utilities::MPI::Partitioner>
                         (local_indices, ghost_indices, get_mpi_communicator());
        }
      else
        {
          comm_pattern =
            std::dynamic_pointer_cast<const Utilities::MPI::Partitioner> (communication_pattern);
          AssertThrow(comm_pattern != nullptr,
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

      // Copy the local elements of tmp_vector to the right place in val
      IndexSet tmp_index_set = tmp_vector.locally_owned_elements();
      for (size_type i=0; i<tmp_index_set.n_elements(); ++i)
        {
          values[locally_owned_elem.index_within_set(tmp_index_set.nth_index_in_set(i))] =
            tmp_vector.local_element(i);
        }
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
              const int ierr = MPI_Testall (update_ghost_values_requests.size(),
                                            update_ghost_values_requests.data(),
                                            &flag, MPI_STATUSES_IGNORE);
              AssertThrowMPI (ierr);
              Assert (flag == 1,
                      ExcMessage("MPI found unfinished update_ghost_values() requests"
                                 "when calling swap, which is not allowed"));
            }
          if (compress_requests.size()>0)
            {
              const int ierr = MPI_Testall (compress_requests.size(), compress_requests.data(),
                                            &flag, MPI_STATUSES_IGNORE);
              AssertThrowMPI (ierr);
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
      std::swap (values,               v.values);
      std::swap (import_data,       v.import_data);
      std::swap (vector_is_ghosted, v.vector_is_ghosted);
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator = (const Number s)
    {
      const size_type this_size = local_size();
      if (this_size>0)
        {
          internal::VectorOperations::Vector_set<Number> setter(s, values.get());

          internal::VectorOperations::parallel_for(setter, 0, this_size,
                                                   thread_loop_partitioner);
        }

      // if we call Vector::operator=0, we want to zero out all the entries
      // plus ghosts.
      if (s==Number())
        zero_out_ghosts();

      return *this;
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const VectorSpaceVector<Number> &V,
                           const bool omit_zeroing_entries)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&V)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);

      reinit(down_V, omit_zeroing_entries);
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator += (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertDimension (local_size(), v.local_size());

      internal::VectorOperations::Vectorization_add_v<Number> vector_add(values.get(), v.values.get());
      internal::VectorOperations::parallel_for(vector_add, 0, partitioner->local_size(),
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
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertDimension (local_size(), v.local_size());

      internal::VectorOperations::Vectorization_subtract_v<Number> vector_subtract(values.get(), v.values.get());
      internal::VectorOperations::parallel_for(vector_subtract, 0, partitioner->local_size(),
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

      internal::VectorOperations::Vectorization_add_factor<Number> vector_add(values.get(), a);
      internal::VectorOperations::parallel_for(vector_add, 0, partitioner->local_size(),
                                               thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::add_local (const Number a,
                               const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertIsFinite(a);
      AssertDimension (local_size(), v.local_size());

      // nothing to do if a is zero
      if (a == Number(0.))
        return;

      internal::VectorOperations::Vectorization_add_av<Number> vector_add(values.get(), v.values.get(), a);
      internal::VectorOperations::parallel_for(vector_add, 0, partitioner->local_size(),
                                               thread_loop_partitioner);
    }



    template <typename Number>
    void
    Vector<Number>::add (const Number a,
                         const VectorSpaceVector<Number> &vv)
    {
      add_local(a,vv);

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
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);
      Assert(dynamic_cast<const Vector<Number> *>(&ww)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &w = dynamic_cast<const Vector<Number> &>(ww);

      AssertIsFinite(a);
      AssertIsFinite(b);

      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());

      internal::VectorOperations::Vectorization_add_avpbw<Number> vector_add(values.get(), v.values.get(),
          w.values.get(), a, b);
      internal::VectorOperations::parallel_for(vector_add, 0, partitioner->local_size(),
                                               thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::add (const std::vector<size_type> &indices,
                         const std::vector<Number>    &values)
    {
      for (std::size_t i=0; i<indices.size(); ++i)
        {
          this->operator()(indices[i]) += values[i];
        }
    }




    template <typename Number>
    void
    Vector<Number>::sadd (const Number x,
                          const Vector<Number> &v)
    {
      AssertIsFinite(x);
      AssertDimension (local_size(), v.local_size());

      internal::VectorOperations::Vectorization_sadd_xv<Number> vector_sadd(values.get(), v.values.get(), x);
      internal::VectorOperations::parallel_for(vector_sadd, 0, partitioner->local_size(),
                                               thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    void
    Vector<Number>::sadd_local (const Number x,
                                const Number a,
                                const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertIsFinite(x);
      AssertIsFinite(a);
      AssertDimension (local_size(), v.local_size());

      internal::VectorOperations::Vectorization_sadd_xav<Number> vector_sadd(values.get(), v.values.get(), a, x);
      internal::VectorOperations::parallel_for(vector_sadd, 0, partitioner->local_size(),
                                               thread_loop_partitioner);
    }



    template <typename Number>
    void
    Vector<Number>::sadd (const Number x,
                          const Number a,
                          const VectorSpaceVector<Number> &vv)
    {
      sadd_local(x,a,vv);

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

      internal::VectorOperations::Vectorization_sadd_xavbw<Number> vector_sadd(values.get(), v.values.get(), w.values.get(),
          x, a, b);
      internal::VectorOperations::parallel_for(vector_sadd, 0, partitioner->local_size(),
                                               thread_loop_partitioner);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator *= (const Number factor)
    {
      AssertIsFinite(factor);
      internal::VectorOperations::Vectorization_multiply_factor<Number> vector_multiply(values.get(),
          factor);

      internal::VectorOperations::parallel_for(vector_multiply, 0, partitioner->local_size(),
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
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertDimension (local_size(), v.local_size());

      internal::VectorOperations::Vectorization_scale<Number> vector_scale(values.get(), v.values.get());
      internal::VectorOperations::parallel_for(vector_scale, 0, partitioner->local_size(),
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
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      AssertIsFinite(a);
      AssertDimension (local_size(), v.local_size());

      internal::VectorOperations::Vectorization_equ_au<Number> vector_equ(values.get(), v.values.get(), a);
      internal::VectorOperations::parallel_for(vector_equ, 0, partitioner->local_size(),
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

      internal::VectorOperations::Vectorization_equ_aubv<Number> vector_equ(values.get(), v.values.get(),
          w.values.get(), a, b);
      internal::VectorOperations::parallel_for(vector_equ, 0, partitioner->local_size(),
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
        if (values[i] != Number(0))
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
                                    partitioner->get_mpi_communicator());
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
      internal::VectorOperations::Dot<Number,Number2> dot(values.get(), v.values.get());
      internal::VectorOperations::parallel_reduce (dot, 0, partitioner->local_size(), sum,
                                                   thread_loop_partitioner);
      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number>
    Number
    Vector<Number>::operator * (const VectorSpaceVector<Number> &vv) const
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);

      Number local_result = inner_product_local(v);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::norm_sqr_local () const
    {
      real_type sum;
      internal::VectorOperations::Norm2<Number,real_type> norm2(values.get());
      internal::VectorOperations::parallel_reduce (norm2, 0, partitioner->local_size(), sum,
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
      internal::VectorOperations::MeanValue<Number> mean(values.get());
      internal::VectorOperations::parallel_reduce (mean, 0, partitioner->local_size(), sum,
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
                                    partitioner->get_mpi_communicator())
               /(real_type)partitioner->size();
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm_local () const
    {
      real_type sum;
      internal::VectorOperations::Norm1<Number, real_type> norm1(values.get());
      internal::VectorOperations::parallel_reduce (norm1, 0, partitioner->local_size(), sum,
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
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::norm_sqr () const
    {
      real_type local_result = norm_sqr_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l2_norm () const
    {
      return std::sqrt(norm_sqr());
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::lp_norm_local (const real_type p) const
    {
      real_type sum;
      internal::VectorOperations::NormP<Number, real_type> normp(values.get(), p);
      internal::VectorOperations::parallel_reduce (normp, 0, partitioner->local_size(), sum,
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
                                             partitioner->get_mpi_communicator()),
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
        max = std::max (numbers::NumberTraits<Number>::abs(values[i]), max);

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
                                    partitioner->get_mpi_communicator());
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
      internal::VectorOperations::AddAndDot<Number> adder(this->values.get(), v.values.get(), w.values.get(), a);
      internal::VectorOperations::parallel_reduce (adder, 0, vec_size, sum, thread_loop_partitioner);
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
      Assert(dynamic_cast<const Vector<Number> *>(&vv)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &v = dynamic_cast<const Vector<Number> &>(vv);
      Assert(dynamic_cast<const Vector<Number> *>(&ww)!=nullptr,
             ExcVectorTypeNotCompatible());
      const Vector<Number> &w = dynamic_cast<const Vector<Number> &>(ww);

      Number local_result = add_and_dot_local(a, v, w);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    partitioner->get_mpi_communicator());
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
      if (import_data != nullptr)
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
      Assert (partitioner.get() != nullptr, ExcInternalError());
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
          {
            const int ierr = MPI_Barrier (partitioner->get_mpi_communicator());
            AssertThrowMPI (ierr);
          }
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
          int ierr = MPI_Barrier (partitioner->get_mpi_communicator());
          AssertThrowMPI (ierr);

          for (unsigned int i=partitioner->this_mpi_process()+1;
               i<partitioner->n_mpi_processes(); i++)
            {
              ierr = MPI_Barrier (partitioner->get_mpi_communicator());
              AssertThrowMPI (ierr);
            }
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
