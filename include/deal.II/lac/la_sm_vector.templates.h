// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_la_sm_vector_templates_h
#define dealii_la_sm_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/cuda_size.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/la_sm_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_operations_internal.h>


DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace SharedMPI
  {
    namespace internal
    {
      template <typename Number, typename MemorySpaceType>
      struct la_parallel_vector_templates_functions
      {
        using size_type = types::global_dof_index;

        static void
        resize_val(const types::global_dof_index new_alloc_size,
                   types::global_dof_index &     allocated_size,
                   MemorySpaceData<Number> &     data,
                   const MPI_Comm &              comm_shared,
                   const bool                    contiguous_allocation_enabled)
        {
#ifndef DEAL_II_WITH_MPI
          Assert(false, ExcNeedsMPI());

          (void)new_alloc_size;
          (void)allocated_size;
          (void)data;
          (void)comm_shared;
#else
          // TODO: is assert fine?
          Assert(((allocated_size > 0 && data.values != nullptr) ||
                  data.values == nullptr),
                 ExcInternalError());

          allocated_size = new_alloc_size;

          const unsigned int size_sm =
            Utilities::MPI::n_mpi_processes(comm_shared);
          const unsigned int rank_sm =
            Utilities::MPI::this_mpi_process(comm_shared);

          MPI_Win *win       = new MPI_Win;
          Number * data_this = (Number *)malloc(0);
          data.others.resize(size_sm);

          MPI_Info info;
          MPI_Info_create(&info);

          if (contiguous_allocation_enabled)
            MPI_Info_set(info, "alloc_shared_noncontig", "false");
          else
            MPI_Info_set(info, "alloc_shared_noncontig", "true");

          const std::size_t align_by = 64;

          std::size_t s = ((new_alloc_size * sizeof(Number) + align_by - 1) /
                           sizeof(Number)) *
                          sizeof(Number);

          data.memory_constumption_values = s;

          MPI_Win_allocate_shared(
            s, sizeof(Number), info, comm_shared, data_this, win);

          for (unsigned int i = 0; i < size_sm; i++)
            {
              int      disp_unit;
              MPI_Aint ssize;
              MPI_Win_shared_query(
                *win, i, &ssize, &disp_unit, &data.others[i]);
            }

          Number *ptr_unaligned = data.others[rank_sm];
          Number *ptr_aligned   = ptr_unaligned;

          AssertThrow(std::align(align_by,
                                 new_alloc_size * sizeof(Number),
                                 reinterpret_cast<void *&>(ptr_aligned),
                                 s) != nullptr,
                      ExcNotImplemented());

          unsigned int              n_align_local = ptr_aligned - ptr_unaligned;
          std::vector<unsigned int> n_align_sm(size_sm);

          MPI_Allgather(&n_align_local,
                        1,
                        MPI_UNSIGNED,
                        n_align_sm.data(),
                        1,
                        MPI_UNSIGNED,
                        comm_shared);

          for (unsigned int i = 0; i < size_sm; i++)
            data.others[i] += n_align_sm[i];

          data.values     = {ptr_aligned,
                         [&data](Number *) { MPI_Win_free(data.values_win); }};
          data.values_win = win;
#endif
        }

        template <typename RealType>
        static void
        linfty_norm_local(const ::dealii::MemorySpace::MemorySpaceData<
                            Number,
                            ::dealii::MemorySpace::Host> &data,
                          const unsigned int              size,
                          RealType &                      max)
        {
          for (size_type i = 0; i < size; ++i)
            max =
              std::max(numbers::NumberTraits<Number>::abs(data.values[i]), max);
        }
      };
    } // namespace internal


    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::clear_mpi_requests()
    {
      // TODO: what is the use of this?
      // it should probably delegated to partitioner_old
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::resize_val(const size_type new_alloc_size,
                                                const MPI_Comm &comm_sm)
    {
      internal::la_parallel_vector_templates_functions<Number,
                                                       MemorySpaceType>::
        resize_val(new_alloc_size,
                   allocated_size,
                   data,
                   comm_sm,
                   this->partitioner->contiguous_allocation_enabled());

      thread_loop_partitioner =
        std::make_shared<::dealii::parallel::internal::TBBPartitioner>();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(const size_type size,
                                            const bool omit_zeroing_entries)
    {
      Assert(false, ExcNotImplemented());
      (void)size;
      (void)omit_zeroing_entries;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const Vector<Number2, MemorySpaceType> &v,
      const bool                              omit_zeroing_entries)
    {
      clear_mpi_requests();
      Assert(v.partitioner_old.get() != nullptr, ExcNotInitialized());

      // check whether the partitioners are
      // different (check only if the are allocated
      // differently, not if the actual data is
      // different)
      if (partitioner_old.get() != v.partitioner_old.get() ||
          setup_ghosts != v.setup_ghosts)
        {
          partitioner_old = v.partitioner_old;
          partitioner     = v.partitioner;
          setup_ghosts    = v.setup_ghosts;

          // AssertThrow(setup_ghosts,
          // dealii::StandardExceptions::ExcNotImplemented());

          const size_type new_allocated_size =
            partitioner->local_size() +
            (setup_ghosts ? partitioner->n_ghost_indices() : 0);
          resize_val(new_allocated_size,
                     partitioner->get_sm_mpi_communicator());
        }

      if (omit_zeroing_entries == false)
        this->operator=(Number());
      else if (setup_ghosts)
        zero_out_ghosts();

      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.clear(); // TODO

      thread_loop_partitioner = v.thread_loop_partitioner;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const IndexSet &ghost_indices,
      const MPI_Comm  communicator)
    {
      Assert(false, ExcNotImplemented());
      (void)locally_owned_indices;
      (void)ghost_indices;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const MPI_Comm  communicator)
    {
      Assert(false, ExcNotImplemented());
      (void)locally_owned_indices;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_old,
      const std::shared_ptr<const PartitionerBase> &            partitioner,
      const bool                                                setup_ghosts)
    {
      // AssertThrow(setup_ghosts,
      // dealii::StandardExceptions::ExcNotImplemented());

      clear_mpi_requests();
      this->partitioner_old = partitioner_old;
      this->partitioner     = partitioner;
      this->setup_ghosts    = setup_ghosts;

      // set vector size and allocate memory
      const size_type new_allocated_size =
        this->partitioner->local_size() +
        (setup_ghosts ? this->partitioner->n_ghost_indices() : 0);
      resize_val(new_allocated_size, partitioner->get_sm_mpi_communicator());

      // initialize to zero
      this->operator=(Number());


      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.clear();

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector()
      : partitioner_old(new Utilities::MPI::Partitioner())
      , allocated_size(0)
    {}



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const bool do_ghost_value_update,
                                            const bool do_compress)
      : partitioner_old(new Utilities::MPI::Partitioner())
      , allocated_size(0)
      , do_ghost_value_update(do_ghost_value_update)
      , do_compress(do_compress)
    {}



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const Vector<Number, MemorySpaceType> &v)
      : Subscriptor()
      , allocated_size(0)
      , vector_is_ghosted(false)
    {
      reinit(v, true);

      thread_loop_partitioner = v.thread_loop_partitioner;

      const size_type this_size = local_size();
      if (this_size > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number, MemorySpaceType>::copy(
              thread_loop_partitioner, partitioner->local_size(), v.data, data);
        }
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const IndexSet &ghost_indices,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)local_range;
      (void)ghost_indices;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)local_range;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const size_type size)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)size;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_old)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)partitioner_old;
    }



    template <typename Number, typename MemorySpaceType>
    inline Vector<Number, MemorySpaceType>::~Vector()
    {
      try
        {
          clear_mpi_requests();
        }
      catch (...)
        {}
    }



    template <typename Number, typename MemorySpaceType>
    inline Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator=(const Vector<Number, MemorySpaceType> &c)
    {
#ifdef _MSC_VER
      return this->operator=<Number>(c);
#else
      return this->template operator=<Number>(c);
#endif
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    inline Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator=(const Vector<Number2, MemorySpaceType> &c)
    {
      Assert(c.partitioner_old.get() != nullptr, ExcNotInitialized());

      // we update ghost values whenever one of the input or output vector
      // already held ghost values or when we import data from a vector with
      // the same local range but different ghost layout
      bool must_update_ghost_values = c.vector_is_ghosted;

      // check whether the two vectors use the same parallel partitioner_old. if
      // not, check if all local ranges are the same (that way, we can
      // exchange data between different parallel layouts). One variant which
      // is included here and necessary for compatibility with the other
      // distributed vector classes (Trilinos, PETSc) is the case when vector
      // c does not have any ghosts (constructed without ghost elements given)
      // but the current vector does: In that case, we need to exchange data
      // also when none of the two vector had updated its ghost values before.
      if (partitioner_old.get() == nullptr)
        reinit(c, true);
      else if (partitioner_old.get() != c.partitioner_old.get())
        {
          // local ranges are also the same if both partitioners are empty
          // (even if they happen to define the empty range as [0,0) or [c,c)
          // for some c!=0 in a different way).
          int local_ranges_are_identical =
            (partitioner_old->local_range() ==
               c.partitioner_old->local_range() ||
             (partitioner_old->local_range().second ==
                partitioner_old->local_range().first &&
              c.partitioner_old->local_range().second ==
                c.partitioner_old->local_range().first));
          if ((c.partitioner_old->n_mpi_processes() > 1 &&
               Utilities::MPI::min(local_ranges_are_identical,
                                   c.partitioner_old->get_mpi_communicator()) ==
                 0) ||
              !local_ranges_are_identical)
            reinit(c, true);
          else
            must_update_ghost_values |= vector_is_ghosted;

          must_update_ghost_values |=
            (c.partitioner_old->ghost_indices_initialized() == false &&
             partitioner_old->ghost_indices_initialized() == true);
        }
      else
        must_update_ghost_values |= vector_is_ghosted;

      thread_loop_partitioner = c.thread_loop_partitioner;

      const size_type this_size = partitioner->local_size();
      if (this_size > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number2, MemorySpaceType>::copy(
              thread_loop_partitioner, this_size, c.data, data);
        }

      if (setup_ghosts)
        {
          if (must_update_ghost_values)
            update_ghost_values();
          else
            zero_out_ghosts();
        }

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::copy_locally_owned_data_from(
      const Vector<Number2, MemorySpaceType> &src)
    {
      Assert(false, ExcNotImplemented());
      (void)src;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename MemorySpaceType2>
    void
    Vector<Number, MemorySpaceType>::import(
      const Vector<Number, MemorySpaceType2> &src,
      VectorOperation::values                 operation)
    {
      Assert(false, ExcNotImplemented());
      (void)src;
      (void)operation;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress(
      ::dealii::VectorOperation::values operation)
    {
      compress_start(0, operation);
      compress_finish(operation);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values() const
    {
      Assert(setup_ghosts, ExcNotImplemented());

      update_ghost_values_start();
      update_ghost_values_finish();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::zero_out_ghosts() const
    {
      Assert(setup_ghosts, ExcNotImplemented());

      if (data.values != nullptr)
        std::fill_n(data.values.get() + partitioner->local_size(),
                    partitioner->n_ghost_indices(),
                    Number());

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_start(
      const unsigned int                communication_channel,
      ::dealii::VectorOperation::values operation)
    {
      (void)operation;
      Assert(::dealii::VectorOperation::values::add == operation,
             ExcNotImplemented());
      Assert(vector_is_ghosted == false,
             ExcMessage("Cannot call compress() on a ghosted vector"));

      if (do_compress)
        partitioner->import_from_ghosted_array_start(operation,
                                                     communication_channel,
                                                     data.values.get(),
                                                     data.others,
                                                     import_data,
                                                     requests);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_finish(
      ::dealii::VectorOperation::values operation)
    {
      vector_is_ghosted = false;

      if (do_compress)
        partitioner->import_from_ghosted_array_finish(
          operation, data.values.get(), data.others, import_data, requests);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values_start(
      const unsigned int communication_channel) const
    {
      if (do_ghost_value_update)
        partitioner->export_to_ghosted_array_start(communication_channel,
                                                   data.values.get(),
                                                   data.others,
                                                   import_data,
                                                   requests);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values_finish() const
    {
      if (do_ghost_value_update)
        partitioner->export_to_ghosted_array_finish(data.values.get(),
                                                    data.others,
                                                    requests);
      vector_is_ghosted = true;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::import(
      const ReadWriteVector<Number> &                 V,
      VectorOperation::values                         operation,
      std::shared_ptr<const CommunicationPatternBase> communication_pattern)
    {
      Assert(false, ExcNotImplemented());
      (void)V;
      (void)operation;
      (void)communication_pattern;
    }

    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::swap(Vector<Number, MemorySpaceType> &v)
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator=(const Number s)
    {
      const size_type this_size = local_size();
      if (this_size > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number, MemorySpaceType>::set(
              thread_loop_partitioner, this_size, s, data);
        }

      // if we call Vector::operator=0, we want to zero out all the entries
      // plus ghosts.
      if (s == Number() && setup_ghosts)
        zero_out_ghosts();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(const VectorSpaceVector<Number> &V,
                                            const bool omit_zeroing_entries)
    {
      Assert(false, ExcNotImplemented());
      (void)V;
      (void)omit_zeroing_entries;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator+=(const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)vv;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator-=(const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      AssertDimension(local_size(), v.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::subtract_vector(
          thread_loop_partitioner, partitioner->local_size(), v.data, data);

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number a)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add_local(
      const Number                     a,
      const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      AssertIsFinite(a);
      AssertDimension(local_size(), v.local_size());

      // nothing to do if a is zero
      if (a == Number(0.))
        return;

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_av(
          thread_loop_partitioner, partitioner->local_size(), a, v.data, data);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number                     a,
                                         const VectorSpaceVector<Number> &vv)
    {
      add_local(a, vv);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number                     a,
                                         const VectorSpaceVector<Number> &vv,
                                         const Number                     b,
                                         const VectorSpaceVector<Number> &ww)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)vv;
      (void)b;
      (void)ww;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const std::vector<size_type> &indices,
                                         const std::vector<Number> &   values)
    {
      Assert(false, ExcNotImplemented());
      (void)indices;
      (void)values;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(
      const Number                           x,
      const Vector<Number, MemorySpaceType> &v)
    {
      Assert(false, ExcNotImplemented());
      (void)x;
      (void)v;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd_local(
      const Number                     x,
      const Number                     a,
      const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert((dynamic_cast<const VectorType *>(&vv) != nullptr),
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      AssertIsFinite(x);
      AssertIsFinite(a);
      AssertDimension(local_size(), v.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::sadd_xav(
          thread_loop_partitioner,
          partitioner->local_size(),
          x,
          a,
          v.data,
          data);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(const Number                     x,
                                          const Number                     a,
                                          const VectorSpaceVector<Number> &vv)
    {
      sadd_local(x, a, vv);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(
      const Number                           x,
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Number                           b,
      const Vector<Number, MemorySpaceType> &w)
    {
      Assert(false, ExcNotImplemented());
      (void)x;
      (void)a;
      (void)v;
      (void)b;
      (void)w;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator*=(const Number factor)
    {
      Assert(false, ExcNotImplemented());
      (void)factor;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator/=(const Number factor)
    {
      Assert(false, ExcNotImplemented());
      (void)factor;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::scale(const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      AssertDimension(local_size(), v.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::scale(
          thread_loop_partitioner, local_size(), v.data, data);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::equ(const Number                     a,
                                         const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      AssertIsFinite(a);
      AssertDimension(local_size(), v.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::equ_au(
          thread_loop_partitioner, partitioner->local_size(), a, v.data, data);


      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::equ(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Number                           b,
      const Vector<Number, MemorySpaceType> &w)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)v;
      (void)b;
      (void)w;
    }



    template <typename Number, typename MemorySpaceType>
    bool
    Vector<Number, MemorySpaceType>::all_zero() const
    {
      return (linfty_norm() == 0) ? true : false;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    Number
    Vector<Number, MemorySpaceType>::inner_product_local(
      const Vector<Number2, MemorySpaceType> &v) const
    {
      if (PointerComparison::equal(this, &v))
        return norm_sqr_local();

      AssertDimension(partitioner->local_size(), v.partitioner->local_size());

      return dealii::internal::VectorOperations::
        functions<Number, Number2, MemorySpaceType>::dot(
          thread_loop_partitioner, partitioner->local_size(), v.data, data);
    }



    template <typename Number, typename MemorySpaceType>
    Number Vector<Number, MemorySpaceType>::
           operator*(const VectorSpaceVector<Number> &vv) const
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert((dynamic_cast<const VectorType *>(&vv) != nullptr),
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      Number local_result = inner_product_local(v);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::norm_sqr_local() const
    {
      real_type sum;


      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::norm_2(
          thread_loop_partitioner, partitioner->local_size(), sum, data);

      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value_local() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l1_norm_local() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l1_norm() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::norm_sqr() const
    {
      real_type local_result = norm_sqr_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l2_norm() const
    {
      return std::sqrt(norm_sqr());
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::lp_norm_local(const real_type p) const
    {
      Assert(false, ExcNotImplemented());
      (void)p;
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::lp_norm(const real_type p) const
    {
      Assert(false, ExcNotImplemented());
      (void)p;
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::linfty_norm_local() const
    {
      real_type max = 0.;

      const size_type local_size = partitioner->local_size();
      internal::la_parallel_vector_templates_functions<
        Number,
        MemorySpaceType>::linfty_norm_local(data, local_size, max);

      return max;
    }



    template <typename Number, typename MemorySpaceType>
    inline typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::linfty_norm() const
    {
      const real_type local_result = linfty_norm_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::max(local_result,
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::add_and_dot_local(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Vector<Number, MemorySpaceType> &w)
    {
      const size_type vec_size = partitioner->local_size();
      AssertDimension(vec_size, v.local_size());
      AssertDimension(vec_size, w.local_size());

      Number sum = dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_and_dot(
          thread_loop_partitioner, vec_size, a, v.data, w.data, data);

      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::add_and_dot(
      const Number                     a,
      const VectorSpaceVector<Number> &vv,
      const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert((dynamic_cast<const VectorType *>(&vv) != nullptr),
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);
      Assert((dynamic_cast<const VectorType *>(&ww) != nullptr),
             ExcVectorTypeNotCompatible());
      const VectorType &w = dynamic_cast<const VectorType &>(ww);

      Number local_result = add_and_dot_local(a, v, w);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number, typename MemorySpaceType>
    inline bool
    Vector<Number, MemorySpaceType>::partitioners_are_compatible(
      const Utilities::MPI::Partitioner &part) const
    {
      return partitioner_old->is_compatible(part);
    }



    template <typename Number, typename MemorySpaceType>
    inline bool
    Vector<Number, MemorySpaceType>::partitioners_are_globally_compatible(
      const Utilities::MPI::Partitioner &part) const
    {
      Assert(false, ExcNotImplemented());
      return partitioner_old->is_globally_compatible(part);
    }



    template <typename Number, typename MemorySpaceType>
    std::size_t
    Vector<Number, MemorySpaceType>::memory_consumption() const
    {
      return data.memory_consumption() + import_data.memory_consumption();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::print(std::ostream &     out,
                                           const unsigned int precision,
                                           const bool         scientific,
                                           const bool         across) const
    {
      Assert(false, ExcNotImplemented());
      (void)out;
      (void)precision;
      (void)scientific;
      (void)across;
    }

  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
