// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_la_parallel_vector_templates_h
#define dealii_la_parallel_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_operations_internal.h>

#include <Kokkos_Core.hpp>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace distributed
  {
    namespace internal
    {
      // In the import_from_ghosted_array_finish we might need to calculate the
      // maximal and minimal value for the given number type, which is not
      // straightforward for complex numbers. Therefore, comparison of complex
      // numbers is prohibited and throws an exception.
      template <typename Number>
      Number
      get_min(const Number a, const Number b)
      {
        return std::min(a, b);
      }

      template <typename Number>
      std::complex<Number>
      get_min(const std::complex<Number> a, const std::complex<Number>)
      {
        AssertThrow(false,
                    ExcMessage("VectorOperation::min not "
                               "implemented for complex numbers"));
        return a;
      }

      template <typename Number>
      Number
      get_max(const Number a, const Number b)
      {
        return std::max(a, b);
      }

      template <typename Number>
      std::complex<Number>
      get_max(const std::complex<Number> a, const std::complex<Number>)
      {
        AssertThrow(false,
                    ExcMessage("VectorOperation::max not "
                               "implemented for complex numbers"));
        return a;
      }



      // Resize the underlying array on the host or on the device
      template <typename Number, typename MemorySpaceType>
      struct la_parallel_vector_templates_functions
      {
        static_assert(std::is_same_v<MemorySpaceType, MemorySpace::Host> ||
                        std::is_same_v<MemorySpaceType, MemorySpace::Default>,
                      "MemorySpace should be Host or Default");

        static void
        resize_val(
          const types::global_dof_index /*new_alloc_size*/,
          types::global_dof_index & /*allocated_size*/,
          ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpaceType>
            & /*data*/,
          const MPI_Comm /*comm_sm*/)
        {}

        static void
        import_elements(
          const ::dealii::LinearAlgebra::ReadWriteVector<Number> & /*V*/,
          VectorOperation::values /*operation*/,
          const std::shared_ptr<const ::dealii::Utilities::MPI::Partitioner> &
          /*communication_pattern*/,
          const IndexSet & /*locally_owned_elem*/,
          ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpaceType>
            & /*data*/)
        {}

        template <typename RealType>
        static void
        linfty_norm_local(
          const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpaceType>
            & /*data*/,
          const unsigned int /*size*/,
          RealType & /*max*/)
        {}
      };

      template <typename Number>
      struct la_parallel_vector_templates_functions<Number,
                                                    ::dealii::MemorySpace::Host>
      {
        using size_type = types::global_dof_index;

        static void
        resize_val(const types::global_dof_index new_alloc_size,
                   types::global_dof_index      &allocated_size,
                   ::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data,
                   const MPI_Comm comm_shared)
        {
          if (comm_shared == MPI_COMM_SELF)
            {
#if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
              Kokkos::resize(Kokkos::WithoutInitializing,
                             data.values,
                             new_alloc_size);
#else
              Kokkos::resize(data.values, new_alloc_size);
#endif

              allocated_size = new_alloc_size;

              data.values_sm = {
                ArrayView<const Number>(data.values.data(), new_alloc_size)};
            }
          else
            {
#ifdef DEAL_II_WITH_MPI
              allocated_size = new_alloc_size;

              const unsigned int size_sm =
                Utilities::MPI::n_mpi_processes(comm_shared);
              const unsigned int rank_sm =
                Utilities::MPI::this_mpi_process(comm_shared);

              MPI_Win mpi_window;
              Number *data_this;


              std::vector<Number *> others(size_sm);

              MPI_Info info;
              int      ierr = MPI_Info_create(&info);
              AssertThrowMPI(ierr);

              ierr = MPI_Info_set(info, "alloc_shared_noncontig", "true");
              AssertThrowMPI(ierr);

              const std::size_t align_by = 64;

              std::size_t s = ((new_alloc_size * sizeof(Number) + align_by) /
                               sizeof(Number)) *
                              sizeof(Number);

              ierr = MPI_Win_allocate_shared(
                s, sizeof(Number), info, comm_shared, &data_this, &mpi_window);
              AssertThrowMPI(ierr);

              for (unsigned int i = 0; i < size_sm; ++i)
                {
                  int        disp_unit;
                  MPI_Aint   ssize;
                  const auto ierr = MPI_Win_shared_query(
                    mpi_window, i, &ssize, &disp_unit, &others[i]);
                  AssertThrowMPI(ierr);
                }

              Number *ptr_unaligned = others[rank_sm];
              Number *ptr_aligned   = ptr_unaligned;

              AssertThrow(std::align(align_by,
                                     new_alloc_size * sizeof(Number),
                                     reinterpret_cast<void *&>(ptr_aligned),
                                     s) != nullptr,
                          ExcNotImplemented());

              unsigned int n_align_local = ptr_aligned - ptr_unaligned;
              std::vector<unsigned int> n_align_sm(size_sm);

              ierr = MPI_Allgather(&n_align_local,
                                   1,
                                   MPI_UNSIGNED,
                                   n_align_sm.data(),
                                   1,
                                   MPI_UNSIGNED,
                                   comm_shared);
              AssertThrowMPI(ierr);

              for (unsigned int i = 0; i < size_sm; ++i)
                others[i] += n_align_sm[i];

              std::vector<unsigned int> new_alloc_sizes(size_sm);

              ierr = MPI_Allgather(&new_alloc_size,
                                   1,
                                   MPI_UNSIGNED,
                                   new_alloc_sizes.data(),
                                   1,
                                   MPI_UNSIGNED,
                                   comm_shared);
              AssertThrowMPI(ierr);

              data.values_sm.resize(size_sm);
              for (unsigned int i = 0; i < size_sm; ++i)
                data.values_sm[i] =
                  ArrayView<const Number>(others[i], new_alloc_sizes[i]);

              data.values =
                Kokkos::View<Number *,
                             Kokkos::HostSpace,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                  ptr_aligned, new_alloc_size);

              // Kokkos will not free the memory because the memory is
              // unmanaged. Instead we use a shared pointer to take care of
              // that.
              data.values_sm_ptr = {ptr_aligned,
                                    [mpi_window](Number *) mutable {
                                      // note: we are creating here a copy of
                                      // the window other approaches led to
                                      // segmentation faults
                                      const auto ierr =
                                        MPI_Win_free(&mpi_window);
                                      AssertThrowMPI(ierr);
                                    }};

#else
              DEAL_II_ASSERT_UNREACHABLE();
#endif
            }
        }

        static void
        import_elements(
          const ::dealii::LinearAlgebra::ReadWriteVector<Number> &V,
          VectorOperation::values                                 operation,
          const std::shared_ptr<const ::dealii::Utilities::MPI::Partitioner>
                         &communication_pattern,
          const IndexSet &locally_owned_elem,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Host>
            &data)
        {
          Assert(
            (operation == VectorOperation::add) ||
              (operation == VectorOperation::insert),
            ExcMessage(
              "Only VectorOperation::add and VectorOperation::insert are allowed"));

          ::dealii::LinearAlgebra::distributed::
            Vector<Number, ::dealii::MemorySpace::Host>
              tmp_vector(communication_pattern);

          // fill entries from ReadWriteVector into the distributed vector,
          // including ghost entries. this is not really efficient right now
          // because indices are translated twice, once by nth_index_in_set(i)
          // and once for operator() of tmp_vector
          const IndexSet &v_stored = V.get_stored_elements();
          for (size_type i = 0; i < v_stored.n_elements(); ++i)
            tmp_vector(v_stored.nth_index_in_set(i)) = V.local_element(i);

          tmp_vector.compress(operation);

          // Copy the local elements of tmp_vector to the right place in val
          IndexSet tmp_index_set = tmp_vector.locally_owned_elements();
          if (operation == VectorOperation::add)
            {
              for (size_type i = 0; i < tmp_index_set.n_elements(); ++i)
                {
                  data.values[locally_owned_elem.index_within_set(
                    tmp_index_set.nth_index_in_set(i))] +=
                    tmp_vector.local_element(i);
                }
            }
          else
            {
              for (size_type i = 0; i < tmp_index_set.n_elements(); ++i)
                {
                  data.values[locally_owned_elem.index_within_set(
                    tmp_index_set.nth_index_in_set(i))] =
                    tmp_vector.local_element(i);
                }
            }
        }

        template <typename RealType>
        static void
        linfty_norm_local(const ::dealii::MemorySpace::MemorySpaceData<
                            Number,
                            ::dealii::MemorySpace::Host> &data,
                          const unsigned int              size,
                          RealType                       &max)
        {
          for (size_type i = 0; i < size; ++i)
            max =
              std::max(numbers::NumberTraits<Number>::abs(data.values[i]), max);
        }
      };

      template <typename Number>
      struct la_parallel_vector_templates_functions<
        Number,
        ::dealii::MemorySpace::Default>
      {
        using size_type = types::global_dof_index;

        static void
        resize_val(
          const types::global_dof_index new_alloc_size,
          types::global_dof_index      &allocated_size,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Default>
            &data,
          const MPI_Comm /*comm_sm*/)
        {
          static_assert(
            std::is_same_v<Number, float> || std::is_same_v<Number, double>,
            "Number should be float or double for Default memory space");

          if (new_alloc_size > allocated_size)
            {
              Assert(((allocated_size > 0 && data.values.size() != 0) ||
                      data.values.size() == 0),
                     ExcInternalError());

#if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
              Kokkos::resize(Kokkos::WithoutInitializing,
                             data.values,
                             new_alloc_size);
#else
              Kokkos::resize(data.values, new_alloc_size);
#endif

              allocated_size = new_alloc_size;
            }
          else if (new_alloc_size == 0)
            {
              Kokkos::resize(data.values, 0);
              allocated_size = 0;
            }
        }

        static void
        import_elements(
          const ReadWriteVector<Number> &V,
          VectorOperation::values        operation,
          std::shared_ptr<const Utilities::MPI::Partitioner>
                         &communication_pattern,
          const IndexSet &locally_owned_elem,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Default>
            &data)
        {
          Assert(
            (operation == VectorOperation::add) ||
              (operation == VectorOperation::insert),
            ExcMessage(
              "Only VectorOperation::add and VectorOperation::insert are allowed"));

          ::dealii::LinearAlgebra::distributed::
            Vector<Number, ::dealii::MemorySpace::Default>
              tmp_vector(communication_pattern);

          // fill entries from ReadWriteVector into the distributed vector,
          // including ghost entries. this is not really efficient right now
          // because indices are translated twice, once by nth_index_in_set(i)
          // and once for operator() of tmp_vector
          const IndexSet                   &v_stored = V.get_stored_elements();
          const size_type                   n_elements = v_stored.n_elements();
          Kokkos::DefaultHostExecutionSpace host_exec;
          Kokkos::View<size_type *, Kokkos::HostSpace> indices(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "indices"),
            n_elements);
          Kokkos::parallel_for(
            "dealii::import_elements: fill indices",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(host_exec,
                                                                   0,
                                                                   n_elements),
            KOKKOS_LAMBDA(size_type i) {
              indices[i] = communication_pattern->global_to_local(
                v_stored.nth_index_in_set(i));
            });
          host_exec.fence();

          // Move the indices to the device
          ::dealii::MemorySpace::Default::kokkos_space::execution_space exec;
          auto indices_dev = Kokkos::create_mirror_view_and_copy(
            ::dealii::MemorySpace::Default::kokkos_space{}, indices);

          // Move the data to the device
          Kokkos::View<Number *, Kokkos::HostSpace> V_view(V.begin(),
                                                           n_elements);
          auto V_dev = Kokkos::create_mirror_view_and_copy(
            ::dealii::MemorySpace::Default::kokkos_space{}, V_view);

          // Set the values in tmp_vector
          Kokkos::parallel_for(
            "dealii::import_elements: set values tmp_vector",
            Kokkos::RangePolicy<
              ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
              exec, 0, n_elements),
            KOKKOS_LAMBDA(size_type i) {
              tmp_vector(indices_dev(i)) = V_dev(i);
            });
          exec.fence();

          tmp_vector.compress(operation);

          // Copy the local elements of tmp_vector to the right place in val
          IndexSet        tmp_index_set  = tmp_vector.locally_owned_elements();
          const size_type tmp_n_elements = tmp_index_set.n_elements();
          Kokkos::realloc(indices, tmp_n_elements);
          Kokkos::parallel_for(
            "dealii::import_elements: copy local elements to val",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(host_exec,
                                                                   0,
                                                                   n_elements),
            KOKKOS_LAMBDA(size_type i) {
              indices[i] = locally_owned_elem.index_within_set(
                tmp_index_set.nth_index_in_set(i));
            });
          host_exec.fence();
          Kokkos::realloc(indices_dev, tmp_n_elements);
          Kokkos::deep_copy(indices_dev,
                            Kokkos::subview(indices,
                                            Kokkos::make_pair(size_type(0),
                                                              tmp_n_elements)));

          if (operation == VectorOperation::add)
            Kokkos::parallel_for(
              "dealii::import_elements: add values",
              Kokkos::RangePolicy<
                ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
                exec, 0, n_elements),
              KOKKOS_LAMBDA(size_type i) {
                data.values(indices_dev(i)) += tmp_vector(i);
              });
          else
            Kokkos::parallel_for(
              "dealii::import_elements: set values",
              Kokkos::RangePolicy<
                ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
                exec, 0, n_elements),
              KOKKOS_LAMBDA(size_type i) {
                data.values(indices_dev(i)) = tmp_vector(i);
              });
          exec.fence();
        }

        template <typename RealType>
        static void
        linfty_norm_local(const ::dealii::MemorySpace::MemorySpaceData<
                            Number,
                            ::dealii::MemorySpace::Default> &data,
                          const unsigned int                 size,
                          RealType                          &result)
        {
          static_assert(std::is_same_v<Number, RealType>,
                        "RealType should be the same type as Number");

          typename ::dealii::MemorySpace::Default::kokkos_space::execution_space
            exec;
          Kokkos::parallel_reduce(
            "dealii::linfty_norm_local",
            Kokkos::RangePolicy<
              ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
              exec, 0, size),
            KOKKOS_LAMBDA(size_type i, RealType & update) {
#if DEAL_II_KOKKOS_VERSION_GTE(3, 7, 0)
              update = Kokkos::fmax(update, Kokkos::abs(data.values(i)));
#else
              update = Kokkos::Experimental::fmax(
                update, Kokkos::Experimental::fabs(data.values(i)));
#endif
            },
            Kokkos::Max<RealType, Kokkos::HostSpace>(result));
        }
      };
    } // namespace internal


    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::clear_mpi_requests()
    {
#ifdef DEAL_II_WITH_MPI
      for (auto &compress_request : compress_requests)
        {
          const int ierr = MPI_Request_free(&compress_request);
          AssertThrowMPI(ierr);
        }
      compress_requests.clear();
      for (auto &update_ghost_values_request : update_ghost_values_requests)
        {
          const int ierr = MPI_Request_free(&update_ghost_values_request);
          AssertThrowMPI(ierr);
        }
      update_ghost_values_requests.clear();
#endif
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::resize_val(const size_type new_alloc_size,
                                                const MPI_Comm  comm_sm)
    {
      internal::la_parallel_vector_templates_functions<
        Number,
        MemorySpaceType>::resize_val(new_alloc_size,
                                     allocated_size,
                                     data,
                                     comm_sm);

      thread_loop_partitioner =
        std::make_shared<::dealii::parallel::internal::TBBPartitioner>();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(const size_type size,
                                            const bool omit_zeroing_entries)
    {
      clear_mpi_requests();

      // check whether we need to reallocate
      resize_val(size, comm_sm);

      // delete previous content in import data
      Kokkos::resize(import_data.values_host_buffer, 0);
      Kokkos::resize(import_data.values, 0);

      // set partitioner to serial version
      partitioner = std::make_shared<Utilities::MPI::Partitioner>(size);

      // set entries to zero if so requested
      if (omit_zeroing_entries == false)
        this->operator=(Number());
      else
        zero_out_ghost_values();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const types::global_dof_index local_size,
      const types::global_dof_index ghost_size,
      const MPI_Comm                comm,
      const MPI_Comm                comm_sm)
    {
      clear_mpi_requests();

      this->comm_sm = comm_sm;

      // check whether we need to reallocate
      resize_val(local_size + ghost_size, comm_sm);

      // delete previous content in import data
      Kokkos::resize(import_data.values_host_buffer, 0);
      Kokkos::resize(import_data.values, 0);

      // create partitioner
      partitioner = std::make_shared<Utilities::MPI::Partitioner>(local_size,
                                                                  ghost_size,
                                                                  comm);

      this->operator=(Number());
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const Vector<Number2, MemorySpaceType> &v,
      const bool                              omit_zeroing_entries)
    {
      clear_mpi_requests();
      Assert(v.partitioner.get() != nullptr, ExcNotInitialized());

      this->comm_sm = v.comm_sm;

      // check whether the partitioners are
      // different (check only if the are allocated
      // differently, not if the actual data is
      // different)
      if (partitioner.get() != v.partitioner.get())
        {
          partitioner = v.partitioner;
          const size_type new_allocated_size =
            partitioner->locally_owned_size() + partitioner->n_ghost_indices();
          resize_val(new_allocated_size, this->comm_sm);
        }

      if (omit_zeroing_entries == false)
        this->operator=(Number());
      else
        zero_out_ghost_values();

      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      Kokkos::resize(import_data.values_host_buffer, 0);
      Kokkos::resize(import_data.values, 0);

      thread_loop_partitioner = v.thread_loop_partitioner;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const IndexSet &ghost_indices,
      const MPI_Comm  communicator)
    {
      // set up parallel partitioner with index sets and communicator
      reinit(std::make_shared<Utilities::MPI::Partitioner>(
        locally_owned_indices, ghost_indices, communicator));
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const MPI_Comm  communicator)
    {
      // set up parallel partitioner with index sets and communicator
      reinit(
        std::make_shared<Utilities::MPI::Partitioner>(locally_owned_indices,
                                                      communicator));
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_in,
      const MPI_Comm                                            comm_sm)
    {
      clear_mpi_requests();

      this->comm_sm = comm_sm;

      // set vector size and allocate memory
      if (partitioner.get() != partitioner_in.get())
        {
          partitioner = partitioner_in;
          const size_type new_allocated_size =
            partitioner->locally_owned_size() + partitioner->n_ghost_indices();
          resize_val(new_allocated_size, comm_sm);
        }

      // initialize to zero
      *this = Number();


      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      Kokkos::resize(import_data.values_host_buffer, 0);
      Kokkos::resize(import_data.values, 0);

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_in,
      const bool /*make_ghosted*/,
      const MPI_Comm &comm_sm)
    {
      this->reinit(partitioner_in, comm_sm);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector()
      : partitioner(std::make_shared<Utilities::MPI::Partitioner>())
      , allocated_size(0)
      , comm_sm(MPI_COMM_SELF)
    {
      reinit(0);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const Vector<Number, MemorySpaceType> &v)
      : allocated_size(0)
      , vector_is_ghosted(false)
      , comm_sm(MPI_COMM_SELF)
    {
      reinit(v, true);

      thread_loop_partitioner = v.thread_loop_partitioner;

      copy_locally_owned_data_from(v);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector( // NOLINT
      Vector<Number, MemorySpaceType> &&v)
      : Vector()
    {
      static_cast<EnableObserverPointer &>(*this) =
        static_cast<EnableObserverPointer &&>(v);
      this->swap(v);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const IndexSet &ghost_indices,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
      , comm_sm(MPI_COMM_SELF)
    {
      reinit(local_range, ghost_indices, communicator);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
      , comm_sm(MPI_COMM_SELF)
    {
      reinit(local_range, communicator);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const size_type size)
      : allocated_size(0)
      , vector_is_ghosted(false)
      , comm_sm(MPI_COMM_SELF)
    {
      reinit(size, false);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      : allocated_size(0)
      , vector_is_ghosted(false)
      , comm_sm(MPI_COMM_SELF)
    {
      reinit(partitioner);
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
    Vector<Number, MemorySpaceType>::operator=(
      const Vector<Number, MemorySpaceType> &c)
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
    Vector<Number, MemorySpaceType>::operator=(
      const Vector<Number2, MemorySpaceType> &c)
    {
      Assert(c.partitioner.get() != nullptr, ExcNotInitialized());

      // we update ghost values whenever one of the input or output vector
      // already held ghost values or when we import data from a vector with
      // the same local range but different ghost layout
      bool must_update_ghost_values = c.vector_is_ghosted;

      this->comm_sm = c.comm_sm;

      // check whether the two vectors use the same parallel partitioner. if
      // not, check if all local ranges are the same (that way, we can
      // exchange data between different parallel layouts). One variant which
      // is included here and necessary for compatibility with the other
      // distributed vector classes (Trilinos, PETSc) is the case when vector
      // c does not have any ghosts (constructed without ghost elements given)
      // but the current vector does: In that case, we need to exchange data
      // also when none of the two vector had updated its ghost values before.
      if (partitioner.get() == nullptr)
        reinit(c, true);
      else if (partitioner.get() != c.partitioner.get())
        {
          // local ranges are also the same if both partitioners are empty
          // (even if they happen to define the empty range as [0,0) or [c,c)
          // for some c!=0 in a different way).
          int local_ranges_are_identical =
            (partitioner->local_range() == c.partitioner->local_range() ||
             (partitioner->local_range().second ==
                partitioner->local_range().first &&
              c.partitioner->local_range().second ==
                c.partitioner->local_range().first));
          if ((c.partitioner->n_mpi_processes() > 1 &&
               Utilities::MPI::min(local_ranges_are_identical,
                                   c.partitioner->get_mpi_communicator()) ==
                 0) ||
              !local_ranges_are_identical)
            reinit(c, true);
          else
            must_update_ghost_values |= vector_is_ghosted;

          must_update_ghost_values |=
            (c.partitioner->ghost_indices_initialized() == false &&
             partitioner->ghost_indices_initialized() == true);
        }
      else
        must_update_ghost_values |= vector_is_ghosted;

      thread_loop_partitioner = c.thread_loop_partitioner;

      copy_locally_owned_data_from(c);

      if (must_update_ghost_values)
        update_ghost_values();
      else
        zero_out_ghost_values();
      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::copy_locally_owned_data_from(
      const Vector<Number2, MemorySpaceType> &src)
    {
      AssertDimension(partitioner->locally_owned_size(),
                      src.partitioner->locally_owned_size());
      if (partitioner->locally_owned_size() > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number2, MemorySpaceType>::copy(
              thread_loop_partitioner,
              partitioner->locally_owned_size(),
              src.data,
              data);
        }
    }



    template <typename Number, typename MemorySpaceType>
    template <typename MemorySpaceType2>
    void
    Vector<Number, MemorySpaceType>::import_elements(
      const Vector<Number, MemorySpaceType2> &src,
      const VectorOperation::values           operation)
    {
      Assert(src.partitioner.get() != nullptr, ExcNotInitialized());
      Assert(partitioner->locally_owned_range() ==
               src.partitioner->locally_owned_range(),
             ExcMessage("Locally owned indices should be identical."));
      Assert(partitioner->ghost_indices() == src.partitioner->ghost_indices(),
             ExcMessage("Ghost indices should be identical."));
      ::dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::import_elements(
          thread_loop_partitioner, allocated_size, operation, src.data, data);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress(VectorOperation::values operation)
    {
      compress_start(0, operation);
      compress_finish(operation);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values() const
    {
      update_ghost_values_start();
      update_ghost_values_finish();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::zero_out_ghost_values() const
    {
      Kokkos::pair<size_type, size_type> range(
        partitioner->locally_owned_size(),
        partitioner->locally_owned_size() + partitioner->n_ghost_indices());
      if (data.values.size() > 0)
        Kokkos::deep_copy(Kokkos::subview(data.values, range), 0);

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_start(
      const unsigned int      communication_channel,
      VectorOperation::values operation)
    {
      AssertIndexRange(communication_channel, 200);
      Assert(vector_is_ghosted == false,
             ExcMessage("Cannot call compress() on a ghosted vector"));

#ifdef DEAL_II_WITH_MPI
      // make this function thread safe
      std::lock_guard<std::mutex> lock(mutex);

      // allocate import_data in case it is not set up yet
      if (partitioner->n_import_indices() > 0)
        {
#  if !defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
          if (std::is_same_v<MemorySpaceType, dealii::MemorySpace::Default>)
            {
              if (import_data.values_host_buffer.size() == 0)
#    if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
                Kokkos::resize(Kokkos::WithoutInitializing,
                               import_data.values_host_buffer,
                               partitioner->n_import_indices());
#    else
                Kokkos::resize(import_data.values_host_buffer,
                               partitioner->n_import_indices());
#    endif
            }
          else
#  endif
            {
              if (import_data.values.size() == 0)
#  if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
                Kokkos::resize(Kokkos::WithoutInitializing,
                               import_data.values,
                               partitioner->n_import_indices());
#  else
              Kokkos::resize(import_data.values,
                             partitioner->n_import_indices());
#  endif
            }
        }

#  if !defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
      if (std::is_same_v<MemorySpaceType, dealii::MemorySpace::Default>)
        {
          // Move the data to the host and then move it back to the
          // device. We use values to store the elements because the function
          // uses a view of the array and thus we need the data on the host to
          // outlive the scope of the function.
          data.values_host_buffer =
#    if DEAL_II_KOKKOS_VERSION_GTE(4, 0, 0)
            Kokkos::create_mirror_view_and_copy(Kokkos::SharedHostPinnedSpace{},
                                                data.values);
#    else
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                                data.values);
#    endif
          partitioner->import_from_ghosted_array_start(
            operation,
            communication_channel,
            ArrayView<Number, MemorySpace::Host>(
              data.values_host_buffer.data() +
                partitioner->locally_owned_size(),
              partitioner->n_ghost_indices()),
            ArrayView<Number, MemorySpace::Host>(
              import_data.values_host_buffer.data(),
              partitioner->n_import_indices()),
            compress_requests);
        }
      else
#  endif
        {
          partitioner->import_from_ghosted_array_start(
            operation,
            communication_channel,
            ArrayView<Number, MemorySpaceType>(
              data.values.data() + partitioner->locally_owned_size(),
              partitioner->n_ghost_indices()),
            ArrayView<Number, MemorySpaceType>(import_data.values.data(),
                                               partitioner->n_import_indices()),
            compress_requests);
        }
#else
      (void)communication_channel;
      (void)operation;
#endif
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_finish(
      VectorOperation::values operation)
    {
#ifdef DEAL_II_WITH_MPI
      vector_is_ghosted = false;

      // in order to zero ghost part of the vector, we need to call
      // import_from_ghosted_array_finish() regardless of
      // compress_requests.empty()

      // make this function thread safe
      std::lock_guard<std::mutex> lock(mutex);
#  if !defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
      if (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
        {
          Assert(partitioner->n_import_indices() == 0 ||
                   import_data.values_host_buffer.size() != 0,
                 ExcNotInitialized());
          partitioner
            ->import_from_ghosted_array_finish<Number, MemorySpace::Host>(
              operation,
              ArrayView<const Number, MemorySpace::Host>(
                import_data.values_host_buffer.data(),
                partitioner->n_import_indices()),
              ArrayView<Number, MemorySpace::Host>(
                data.values_host_buffer.data(),
                partitioner->locally_owned_size()),
              ArrayView<Number, MemorySpace::Host>(
                data.values_host_buffer.data() +
                  partitioner->locally_owned_size(),
                partitioner->n_ghost_indices()),
              compress_requests);

          // The communication is done on the host, so we need to
          // move the data back to the device.
          Kokkos::deep_copy(data.values, data.values_host_buffer);

          Kokkos::resize(data.values_host_buffer, 0);
        }
      else
#  endif
        {
          Assert(partitioner->n_import_indices() == 0 ||
                   import_data.values.size() != 0,
                 ExcNotInitialized());
          partitioner
            ->import_from_ghosted_array_finish<Number, MemorySpaceType>(
              operation,
              ArrayView<const Number, MemorySpaceType>(
                import_data.values.data(), partitioner->n_import_indices()),
              ArrayView<Number, MemorySpaceType>(
                data.values.data(), partitioner->locally_owned_size()),
              ArrayView<Number, MemorySpaceType>(
                data.values.data() + partitioner->locally_owned_size(),
                partitioner->n_ghost_indices()),
              compress_requests);
        }
#else
      (void)operation;
#endif
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values_start(
      const unsigned int communication_channel) const
    {
      AssertIndexRange(communication_channel, 200);
#ifdef DEAL_II_WITH_MPI
      // nothing to do when we neither have import nor ghost indices.
      if (partitioner->n_ghost_indices() == 0 &&
          partitioner->n_import_indices() == 0)
        return;

      // make this function thread safe
      std::lock_guard<std::mutex> lock(mutex);

      // allocate import_data in case it is not set up yet
      if (partitioner->n_import_indices() > 0)
        {
#  if !defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
          if (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
            {
              if (import_data.values_host_buffer.size() == 0)
#    if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
                Kokkos::resize(Kokkos::WithoutInitializing,
                               import_data.values_host_buffer,
                               partitioner->n_import_indices());
#    else
                Kokkos::resize(import_data.values_host_buffer,
                               partitioner->n_import_indices());
#    endif
            }
          else
#  endif
            {
              if (import_data.values.size() == 0)
#  if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
                Kokkos::resize(Kokkos::WithoutInitializing,
                               import_data.values,
                               partitioner->n_import_indices());
#  else
              Kokkos::resize(import_data.values,
                             partitioner->n_import_indices());
#  endif
            }
        }

#  if !defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
      if (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
        {
          // Move the data to the host and then move it back to the
          // device. We use values to store the elements because the function
          // uses a view of the array and thus we need the data on the host to
          // outlive the scope of the function.
          data.values_host_buffer =
#    if DEAL_II_KOKKOS_VERSION_GTE(4, 0, 0)
            Kokkos::create_mirror_view_and_copy(Kokkos::SharedHostPinnedSpace{},
                                                data.values);
#    else
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                                data.values);
#    endif

          partitioner->export_to_ghosted_array_start<Number, MemorySpace::Host>(
            communication_channel,
            ArrayView<const Number, MemorySpace::Host>(
              data.values_host_buffer.data(),
              partitioner->locally_owned_size()),
            ArrayView<Number, MemorySpace::Host>(
              import_data.values_host_buffer.data(),
              partitioner->n_import_indices()),
            ArrayView<Number, MemorySpace::Host>(
              data.values_host_buffer.data() +
                partitioner->locally_owned_size(),
              partitioner->n_ghost_indices()),
            update_ghost_values_requests);
        }
      else
#  endif
        {
          partitioner->export_to_ghosted_array_start<Number, MemorySpaceType>(
            communication_channel,
            ArrayView<const Number, MemorySpaceType>(
              data.values.data(), partitioner->locally_owned_size()),
            ArrayView<Number, MemorySpaceType>(import_data.values.data(),
                                               partitioner->n_import_indices()),
            ArrayView<Number, MemorySpaceType>(
              data.values.data() + partitioner->locally_owned_size(),
              partitioner->n_ghost_indices()),
            update_ghost_values_requests);
        }

#else
      (void)communication_channel;
#endif
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values_finish() const
    {
#ifdef DEAL_II_WITH_MPI
      // wait for both sends and receives to complete, even though only
      // receives are really necessary. this gives (much) better performance
      AssertDimension(partitioner->ghost_targets().size() +
                        partitioner->import_targets().size(),
                      update_ghost_values_requests.size());
      if (update_ghost_values_requests.size() > 0)
        {
          // make this function thread safe
          std::lock_guard<std::mutex> lock(mutex);

#  if !defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
          if (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
            {
              partitioner->export_to_ghosted_array_finish(
                ArrayView<Number, MemorySpace::Host>(
                  data.values_host_buffer.data() +
                    partitioner->locally_owned_size(),
                  partitioner->n_ghost_indices()),
                update_ghost_values_requests);

              // The communication is done on the host, so we need to
              // move the data back to the device.
              auto range = Kokkos::make_pair(partitioner->locally_owned_size(),
                                             partitioner->locally_owned_size() +
                                               partitioner->n_ghost_indices());
              Kokkos::deep_copy(Kokkos::subview(data.values, range),
                                Kokkos::subview(data.values_host_buffer,
                                                range));

              Kokkos::resize(data.values_host_buffer, 0);
            }
          else
#  endif
            {
              partitioner->export_to_ghosted_array_finish(
                ArrayView<Number, MemorySpaceType>(
                  data.values.data() + partitioner->locally_owned_size(),
                  partitioner->n_ghost_indices()),
                update_ghost_values_requests);
            }
        }

#endif
      vector_is_ghosted = true;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::import_elements(
      const ReadWriteVector<Number> &V,
      const VectorOperation::values  operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern)
    {
      // If no communication pattern is given, create one. Otherwise, use the
      // given one.
      std::shared_ptr<const Utilities::MPI::Partitioner> comm_pattern;
      if (communication_pattern.get() == nullptr)
        {
          // Split the IndexSet of V in locally owned elements and ghost indices
          // then create the communication pattern
          IndexSet locally_owned_elem = locally_owned_elements();
          IndexSet ghost_indices      = V.get_stored_elements();
          ghost_indices.subtract_set(locally_owned_elem);
          comm_pattern = std::make_shared<Utilities::MPI::Partitioner>(
            locally_owned_elem, ghost_indices, get_mpi_communicator());
        }
      else
        {
          comm_pattern =
            std::dynamic_pointer_cast<const Utilities::MPI::Partitioner>(
              communication_pattern);
          AssertThrow(comm_pattern != nullptr,
                      ExcMessage("The communication pattern is not of type "
                                 "Utilities::MPI::Partitioner."));
        }
      Vector<Number, ::dealii::MemorySpace::Host> tmp_vector(comm_pattern);

      data.copy_to(tmp_vector.begin(), locally_owned_size());

      // fill entries from ReadWriteVector into the distributed vector,
      // including ghost entries. this is not really efficient right now
      // because indices are translated twice, once by nth_index_in_set(i) and
      // once for operator() of tmp_vector
      const IndexSet &v_stored     = V.get_stored_elements();
      const size_type v_n_elements = v_stored.n_elements();
      switch (operation)
        {
          case VectorOperation::insert:
            {
              for (size_type i = 0; i < v_n_elements; ++i)
                tmp_vector(v_stored.nth_index_in_set(i)) = V.local_element(i);

              break;
            }
          case VectorOperation::add:
            {
              for (size_type i = 0; i < v_n_elements; ++i)
                tmp_vector(v_stored.nth_index_in_set(i)) += V.local_element(i);

              break;
            }
          case VectorOperation::min:
            {
              for (size_type i = 0; i < v_n_elements; ++i)
                tmp_vector(v_stored.nth_index_in_set(i)) =
                  internal::get_min(tmp_vector(v_stored.nth_index_in_set(i)),
                                    V.local_element(i));

              break;
            }
          case VectorOperation::max:
            {
              for (size_type i = 0; i < v_n_elements; ++i)
                tmp_vector(v_stored.nth_index_in_set(i)) =
                  internal::get_max(tmp_vector(v_stored.nth_index_in_set(i)),
                                    V.local_element(i));

              break;
            }
          default:
            {
              Assert(false, ExcMessage("This operation is not supported."));
            }
        }
      tmp_vector.compress(operation);

      data.copy_from(tmp_vector.begin(), locally_owned_size());
    }

    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::swap(
      Vector<Number, MemorySpaceType> &v) noexcept
    {
#ifdef DEAL_II_WITH_MPI

      if constexpr (running_in_debug_mode())
        {
          Assert(Utilities::MPI::job_supports_mpi() ||
                   (update_ghost_values_requests.empty() &&
                    compress_requests.empty()),
                 ExcInternalError());

          // make sure that there are not outstanding requests from updating
          // ghost values or compress
          if (update_ghost_values_requests.size() > 0)
            {
              int       flag = 1;
              const int ierr = MPI_Testall(update_ghost_values_requests.size(),
                                           update_ghost_values_requests.data(),
                                           &flag,
                                           MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
              Assert(flag == 1,
                     ExcMessage(
                       "MPI found unfinished update_ghost_values() requests "
                       "when calling swap, which is not allowed."));
            }
          if (compress_requests.size() > 0)
            {
              int       flag = 1;
              const int ierr = MPI_Testall(compress_requests.size(),
                                           compress_requests.data(),
                                           &flag,
                                           MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
              Assert(flag == 1,
                     ExcMessage("MPI found unfinished compress() requests "
                                "when calling swap, which is not allowed."));
            }
        }

      std::swap(compress_requests, v.compress_requests);
      std::swap(update_ghost_values_requests, v.update_ghost_values_requests);
      std::swap(comm_sm, v.comm_sm);
#endif

      std::swap(partitioner, v.partitioner);
      std::swap(thread_loop_partitioner, v.thread_loop_partitioner);
      std::swap(allocated_size, v.allocated_size);
      std::swap(data, v.data);
      std::swap(import_data, v.import_data);
      std::swap(vector_is_ghosted, v.vector_is_ghosted);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator=( // NOLINT
      Vector<Number, MemorySpaceType> &&v)
    {
      static_cast<EnableObserverPointer &>(*this) =
        static_cast<EnableObserverPointer &&>(v);
      this->swap(v);
      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator=(const Number s)
    {
      const size_type this_size = locally_owned_size();
      if (this_size > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number, MemorySpaceType>::set(
              thread_loop_partitioner, this_size, s, data);
        }

      // if we call Vector::operator=0, we want to zero out all the entries
      // plus ghosts.
      if (s == Number())
        zero_out_ghost_values();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator+=(
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_vector(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          v.data,
          data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator-=(
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::subtract_vector(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          v.data,
          data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number a)
    {
      AssertIsFinite(a);

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_factor(
          thread_loop_partitioner, partitioner->locally_owned_size(), a, data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add_local(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertIsFinite(a);
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      // nothing to do if a is zero
      if (a == Number(0.))
        return;

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_av(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          a,
          v.data,
          data);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &vv)
    {
      add_local(a, vv);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Number                           b,
      const Vector<Number, MemorySpaceType> &w)
    {
      AssertIsFinite(a);
      AssertIsFinite(b);

      AssertDimension(locally_owned_size(), v.locally_owned_size());
      AssertDimension(locally_owned_size(), w.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_avpbw(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          a,
          b,
          v.data,
          w.data,
          data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const std::vector<size_type> &indices,
                                         const std::vector<Number>    &values)
    {
      for (std::size_t i = 0; i < indices.size(); ++i)
        {
          this->operator()(indices[i]) += values[i];
        }
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(
      const Number                           x,
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertIsFinite(x);
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::sadd_xv(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          x,
          v.data,
          data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd_local(
      const Number                           x,
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertIsFinite(x);
      AssertIsFinite(a);
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::sadd_xav(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          x,
          a,
          v.data,
          data);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(
      const Number                           x,
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v)
    {
      sadd_local(x, a, v);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator*=(const Number factor)
    {
      AssertIsFinite(factor);

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::multiply_factor(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          factor,
          data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator/=(const Number factor)
    {
      operator*=(static_cast<Number>(1.) / factor);
      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::scale(
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::scale(
          thread_loop_partitioner, locally_owned_size(), v.data, data);

      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::equ(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v)
    {
      AssertIsFinite(a);
      AssertDimension(locally_owned_size(), v.locally_owned_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::equ_au(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          a,
          v.data,
          data);


      if (vector_is_ghosted)
        update_ghost_values();
      else
        assert_no_residual_content_in_ghost_region();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number,
           MemorySpaceType>::assert_no_residual_content_in_ghost_region() const
    {
      if constexpr (running_in_debug_mode())
        {
          // This should only be called for non-ghosted vectors
          Assert(!vector_is_ghosted, ExcInternalError());

          // Run a reduction over the ghost range only to find out whether some
          // entries are non-zero
          real_type sum = real_type();
          dealii::internal::VectorOperations::
            functions<Number, Number, MemorySpaceType>::norm_1(
              thread_loop_partitioner,
              partitioner->n_ghost_indices(),
              sum,
              data,
              partitioner->locally_owned_size());

          Assert(sum == real_type(),
                 ExcMessage(
                   "You called a vector space operation like add(), "
                   "scale(), operator* for a non-ghosted vector, which "
                   "will not update the content in the memory locations "
                   "reserved for ghost values. However, a non-zero "
                   "content was detected for some of those entries, which "
                   "can lead to an invalid state of the vector. Please "
                   "call Vector::compress(VectorOperation::add) or "
                   "Vector::zero_out_ghost_values() before calling a "
                   "vector space operation to avoid this problem."));
        }
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::extract_subvector_to(
      const ArrayView<const types::global_dof_index> &indices,
      const ArrayView<Number>                        &elements) const
    {
      AssertDimension(indices.size(), elements.size());
      for (unsigned int i = 0; i < indices.size(); ++i)
        {
          AssertIndexRange(indices[i], size());
          elements[i] = (*this)[indices[i]];
        }
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

      AssertDimension(partitioner->locally_owned_size(),
                      v.partitioner->locally_owned_size());

      return dealii::internal::VectorOperations::
        functions<Number, Number2, MemorySpaceType>::dot(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          v.data,
          data);
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::operator*(
      const Vector<Number, MemorySpaceType> &v) const
    {
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
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          sum,
          data);

      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value_local() const
    {
      Assert(size() != 0, ExcEmptyObject());

      if (partitioner->locally_owned_size() == 0)
        return Number();

      Number sum = ::dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::mean_value(
          thread_loop_partitioner, partitioner->locally_owned_size(), data);

      return sum / real_type(partitioner->locally_owned_size());
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value() const
    {
      Number local_result = mean_value_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result *
                                     static_cast<real_type>(
                                       partitioner->locally_owned_size()),
                                   partitioner->get_mpi_communicator()) /
               static_cast<real_type>(partitioner->size());
      else
        return local_result;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l1_norm_local() const
    {
      real_type sum;

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::norm_1(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          sum,
          data);

      return sum;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l1_norm() const
    {
      real_type local_result = l1_norm_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_mpi_communicator());
      else
        return local_result;
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
      real_type sum = 0.;

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::norm_p(
          thread_loop_partitioner,
          partitioner->locally_owned_size(),
          sum,
          p,
          data);

      return std::pow(sum, 1. / p);
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::lp_norm(const real_type p) const
    {
      const real_type local_result = lp_norm_local(p);
      if (partitioner->n_mpi_processes() > 1)
        return std::pow(
          Utilities::MPI::sum(std::pow(local_result, p),
                              partitioner->get_mpi_communicator()),
          static_cast<real_type>(1.0 / p));
      else
        return local_result;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::linfty_norm_local() const
    {
      real_type max = 0.;

      const size_type locally_owned_size = partitioner->locally_owned_size();
      internal::la_parallel_vector_templates_functions<
        Number,
        MemorySpaceType>::linfty_norm_local(data, locally_owned_size, max);

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
      const size_type vec_size = partitioner->locally_owned_size();
      AssertDimension(vec_size, v.locally_owned_size());
      AssertDimension(vec_size, w.locally_owned_size());

      Number sum = dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_and_dot(
          thread_loop_partitioner, vec_size, a, v.data, w.data, data);

      AssertIsFinite(sum);

      return sum;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::add_and_dot(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Vector<Number, MemorySpaceType> &w)
    {
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
      return partitioner->is_compatible(part);
    }



    template <typename Number, typename MemorySpaceType>
    inline bool
    Vector<Number, MemorySpaceType>::partitioners_are_globally_compatible(
      const Utilities::MPI::Partitioner &part) const
    {
      return partitioner->is_globally_compatible(part);
    }



    template <typename Number, typename MemorySpaceType>
    std::size_t
    Vector<Number, MemorySpaceType>::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      memory += sizeof(Number) * static_cast<std::size_t>(allocated_size);

      // if the partitioner is shared between more processors, just count a
      // fraction of that memory, since we're not actually using more memory
      // for it.
      if (partitioner.use_count() > 0)
        memory +=
          partitioner->memory_consumption() / partitioner.use_count() + 1;
      if (import_data.values_host_buffer.size() != 0 ||
          import_data.values.size() != 0)
        memory += (static_cast<std::size_t>(partitioner->n_import_indices()) *
                   sizeof(Number));
      return memory;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::print(std::ostream      &out,
                                           const unsigned int precision,
                                           const bool         scientific,
                                           const bool         across) const
    {
      Assert(partitioner.get() != nullptr, ExcInternalError());
      AssertThrow(out.fail() == false, ExcIO());
      std::ios::fmtflags old_flags     = out.flags();
      unsigned int       old_precision = out.precision(precision);

      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

        // to make the vector write out all the information in order, use as
        // many barriers as there are processors and start writing when it's our
        // turn
#ifdef DEAL_II_WITH_MPI
      if (partitioner->n_mpi_processes() > 1)
        for (unsigned int i = 0; i < partitioner->this_mpi_process(); ++i)
          {
            const int ierr = MPI_Barrier(partitioner->get_mpi_communicator());
            AssertThrowMPI(ierr);
          }
#endif

      std::vector<Number> stored_elements(allocated_size);
      data.copy_to(stored_elements.data(), allocated_size);

      out << "Process #" << partitioner->this_mpi_process() << std::endl
          << "Local range: [" << partitioner->local_range().first << ", "
          << partitioner->local_range().second
          << "), global size: " << partitioner->size() << std::endl
          << "Vector data:" << std::endl;
      if (across)
        for (size_type i = 0; i < partitioner->locally_owned_size(); ++i)
          out << stored_elements[i] << ' ';
      else
        for (size_type i = 0; i < partitioner->locally_owned_size(); ++i)
          out << stored_elements[i] << std::endl;
      out << std::endl;

      if (vector_is_ghosted)
        {
          out << "Ghost entries (global index / value):" << std::endl;
          if (across)
            for (size_type i = 0; i < partitioner->n_ghost_indices(); ++i)
              out << '(' << partitioner->ghost_indices().nth_index_in_set(i)
                  << '/'
                  << stored_elements[partitioner->locally_owned_size() + i]
                  << ") ";
          else
            for (size_type i = 0; i < partitioner->n_ghost_indices(); ++i)
              out << '(' << partitioner->ghost_indices().nth_index_in_set(i)
                  << '/'
                  << stored_elements[partitioner->locally_owned_size() + i]
                  << ')' << std::endl;
          out << std::endl;
        }
      out << std::flush;

#ifdef DEAL_II_WITH_MPI
      if (partitioner->n_mpi_processes() > 1)
        {
          int ierr = MPI_Barrier(partitioner->get_mpi_communicator());
          AssertThrowMPI(ierr);

          for (unsigned int i = partitioner->this_mpi_process() + 1;
               i < partitioner->n_mpi_processes();
               i++)
            {
              ierr = MPI_Barrier(partitioner->get_mpi_communicator());
              AssertThrowMPI(ierr);
            }
        }
#endif

      AssertThrow(out.fail() == false, ExcIO());
      // reset output format
      out.flags(old_flags);
      out.precision(old_precision);
    }

  } // end of namespace distributed
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
