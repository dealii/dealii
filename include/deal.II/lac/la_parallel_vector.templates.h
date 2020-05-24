// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

#ifndef dealii_la_parallel_vector_templates_h
#define dealii_la_parallel_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/cuda_size.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_operations_internal.h>

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
        static_assert(std::is_same<MemorySpaceType, MemorySpace::Host>::value ||
                        std::is_same<MemorySpaceType, MemorySpace::CUDA>::value,
                      "MemorySpace should be Host or CUDA");

        static void
        resize_val(
          const types::global_dof_index /*new_alloc_size*/,
          types::global_dof_index & /*allocated_size*/,
          ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpaceType>
            & /*data*/)
        {}

        static void
        import(
          const ::dealii::LinearAlgebra::ReadWriteVector<Number> & /*V*/,
          ::dealii::VectorOperation::values /*operation*/,
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
                   types::global_dof_index &     allocated_size,
                   ::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          if (new_alloc_size > allocated_size)
            {
              Assert(((allocated_size > 0 && data.values != nullptr) ||
                      data.values == nullptr),
                     ExcInternalError());

              Number *new_val;
              Utilities::System::posix_memalign(
                reinterpret_cast<void **>(&new_val),
                64,
                sizeof(Number) * new_alloc_size);
              data.values.reset(new_val);

              allocated_size = new_alloc_size;
            }
          else if (new_alloc_size == 0)
            {
              data.values.reset();
              allocated_size = 0;
            }
        }

        static void
        import(
          const ::dealii::LinearAlgebra::ReadWriteVector<Number> &V,
          ::dealii::VectorOperation::values                       operation,
          const std::shared_ptr<const ::dealii::Utilities::MPI::Partitioner>
            &             communication_pattern,
          const IndexSet &locally_owned_elem,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Host>
            &data)
        {
          Assert(
            (operation == ::dealii::VectorOperation::add) ||
              (operation == ::dealii::VectorOperation::insert),
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
                          RealType &                      max)
        {
          for (size_type i = 0; i < size; ++i)
            max =
              std::max(numbers::NumberTraits<Number>::abs(data.values[i]), max);
        }
      };

#ifdef DEAL_II_COMPILER_CUDA_AWARE
      template <typename Number>
      struct la_parallel_vector_templates_functions<Number,
                                                    ::dealii::MemorySpace::CUDA>
      {
        using size_type = types::global_dof_index;

        static void
        resize_val(const types::global_dof_index new_alloc_size,
                   types::global_dof_index &     allocated_size,
                   ::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          static_assert(
            std::is_same<Number, float>::value ||
              std::is_same<Number, double>::value,
            "Number should be float or double for CUDA memory space");

          if (new_alloc_size > allocated_size)
            {
              Assert(((allocated_size > 0 && data.values_dev != nullptr) ||
                      data.values_dev == nullptr),
                     ExcInternalError());

              Number *new_val_dev;
              Utilities::CUDA::malloc(new_val_dev, new_alloc_size);
              data.values_dev.reset(new_val_dev);

              allocated_size = new_alloc_size;
            }
          else if (new_alloc_size == 0)
            {
              data.values_dev.reset();
              allocated_size = 0;
            }
        }

        static void
        import(const ReadWriteVector<Number> &V,
               VectorOperation::values        operation,
               std::shared_ptr<const Utilities::MPI::Partitioner>
                               communication_pattern,
               const IndexSet &locally_owned_elem,
               ::dealii::MemorySpace::
                 MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          Assert(
            (operation == ::dealii::VectorOperation::add) ||
              (operation == ::dealii::VectorOperation::insert),
            ExcMessage(
              "Only VectorOperation::add and VectorOperation::insert are allowed"));

          ::dealii::LinearAlgebra::distributed::
            Vector<Number, ::dealii::MemorySpace::CUDA>
              tmp_vector(communication_pattern);

          // fill entries from ReadWriteVector into the distributed vector,
          // including ghost entries. this is not really efficient right now
          // because indices are translated twice, once by nth_index_in_set(i)
          // and once for operator() of tmp_vector
          const IndexSet &       v_stored   = V.get_stored_elements();
          const size_type        n_elements = v_stored.n_elements();
          std::vector<size_type> indices(n_elements);
          for (size_type i = 0; i < n_elements; ++i)
            indices[i] = communication_pattern->global_to_local(
              v_stored.nth_index_in_set(i));
          // Move the indices to the device
          size_type *indices_dev;
          ::dealii::Utilities::CUDA::malloc(indices_dev, n_elements);
          ::dealii::Utilities::CUDA::copy_to_dev(indices, indices_dev);
          // Move the data to the device
          Number *V_dev;
          ::dealii::Utilities::CUDA::malloc(V_dev, n_elements);
          cudaError_t cuda_error_code = cudaMemcpy(V_dev,
                                                   V.begin(),
                                                   n_elements * sizeof(Number),
                                                   cudaMemcpyHostToDevice);
          AssertCuda(cuda_error_code);

          // Set the values in tmp_vector
          const int n_blocks =
            1 + n_elements / (::dealii::CUDAWrappers::chunk_size *
                              ::dealii::CUDAWrappers::block_size);
          ::dealii::LinearAlgebra::CUDAWrappers::kernel::set_permutated<Number>
            <<<n_blocks, ::dealii::CUDAWrappers::block_size>>>(
              indices_dev, tmp_vector.begin(), V_dev, n_elements);

          tmp_vector.compress(operation);

          // Copy the local elements of tmp_vector to the right place in val
          IndexSet        tmp_index_set  = tmp_vector.locally_owned_elements();
          const size_type tmp_n_elements = tmp_index_set.n_elements();
          indices.resize(tmp_n_elements);
          for (size_type i = 0; i < tmp_n_elements; ++i)
            indices[i] = locally_owned_elem.index_within_set(
              tmp_index_set.nth_index_in_set(i));
          ::dealii::Utilities::CUDA::free(indices_dev);
          ::dealii::Utilities::CUDA::malloc(indices_dev, tmp_n_elements);
          ::dealii::Utilities::CUDA::copy_to_dev(indices, indices_dev);

          if (operation == VectorOperation::add)
            ::dealii::LinearAlgebra::CUDAWrappers::kernel::add_permutated<
              Number><<<n_blocks, ::dealii::CUDAWrappers::block_size>>>(
              indices_dev,
              data.values_dev.get(),
              tmp_vector.begin(),
              tmp_n_elements);
          else
            ::dealii::LinearAlgebra::CUDAWrappers::kernel::set_permutated<
              Number><<<n_blocks, ::dealii::CUDAWrappers::block_size>>>(
              indices_dev,
              data.values_dev.get(),
              tmp_vector.begin(),
              tmp_n_elements);

          ::dealii::Utilities::CUDA::free(indices_dev);
          ::dealii::Utilities::CUDA::free(V_dev);
        }

        template <typename RealType>
        static void
        linfty_norm_local(const ::dealii::MemorySpace::MemorySpaceData<
                            Number,
                            ::dealii::MemorySpace::CUDA> &data,
                          const unsigned int              size,
                          RealType &                      result)
        {
          static_assert(std::is_same<Number, RealType>::value,
                        "RealType should be the same type as Number");

          Number *    result_device;
          cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
          AssertCuda(error_code);
          error_code = cudaMemset(result_device, 0, sizeof(Number));

          const int n_blocks = 1 + size / (::dealii::CUDAWrappers::chunk_size *
                                           ::dealii::CUDAWrappers::block_size);
          ::dealii::LinearAlgebra::CUDAWrappers::kernel::reduction<
            Number,
            ::dealii::LinearAlgebra::CUDAWrappers::kernel::LInfty<Number>>
            <<<dim3(n_blocks, 1), dim3(::dealii::CUDAWrappers::block_size)>>>(
              result_device, data.values_dev.get(), size);

          // Copy the result back to the host
          error_code = cudaMemcpy(&result,
                                  result_device,
                                  sizeof(Number),
                                  cudaMemcpyDeviceToHost);
          AssertCuda(error_code);
          // Free the memory on the device
          error_code = cudaFree(result_device);
          AssertCuda(error_code);
        }
      };
#endif
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
    Vector<Number, MemorySpaceType>::resize_val(const size_type new_alloc_size)
    {
      internal::la_parallel_vector_templates_functions<
        Number,
        MemorySpaceType>::resize_val(new_alloc_size, allocated_size, data);

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
      resize_val(size);

      // delete previous content in import data
      import_data.values.reset();
      import_data.values_dev.reset();

      // set partitioner to serial version
      partitioner = std::make_shared<Utilities::MPI::Partitioner>(size);

      // set entries to zero if so requested
      if (omit_zeroing_entries == false)
        this->operator=(Number());
      else
        zero_out_ghosts();
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

      // check whether the partitioners are
      // different (check only if the are allocated
      // differently, not if the actual data is
      // different)
      if (partitioner.get() != v.partitioner.get())
        {
          partitioner = v.partitioner;
          const size_type new_allocated_size =
            partitioner->local_size() + partitioner->n_ghost_indices();
          resize_val(new_allocated_size);
        }

      if (omit_zeroing_entries == false)
        this->operator=(Number());
      else
        zero_out_ghosts();

      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.values.reset();
      import_data.values_dev.reset();

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
      std::shared_ptr<const Utilities::MPI::Partitioner> new_partitioner(
        new Utilities::MPI::Partitioner(locally_owned_indices,
                                        ghost_indices,
                                        communicator));
      reinit(new_partitioner);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const MPI_Comm  communicator)
    {
      // set up parallel partitioner with index sets and communicator
      std::shared_ptr<const Utilities::MPI::Partitioner> new_partitioner(
        new Utilities::MPI::Partitioner(locally_owned_indices, communicator));
      reinit(new_partitioner);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_in)
    {
      clear_mpi_requests();
      partitioner = partitioner_in;

      // set vector size and allocate memory
      const size_type new_allocated_size =
        partitioner->local_size() + partitioner->n_ghost_indices();
      resize_val(new_allocated_size);

      // initialize to zero
      this->operator=(Number());


      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.values.reset();
      import_data.values_dev.reset();

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector()
      : partitioner(new Utilities::MPI::Partitioner())
      , allocated_size(0)
    {
      reinit(0);
    }



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
      reinit(local_range, ghost_indices, communicator);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      reinit(local_range, communicator);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const size_type size)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      reinit(size, false);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      : allocated_size(0)
      , vector_is_ghosted(false)
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
      Assert(c.partitioner.get() != nullptr, ExcNotInitialized());

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

      const size_type this_size = partitioner->local_size();
      if (this_size > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number2, MemorySpaceType>::copy(
              thread_loop_partitioner, this_size, c.data, data);
        }

      if (must_update_ghost_values)
        update_ghost_values();
      else
        zero_out_ghosts();
      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::copy_locally_owned_data_from(
      const Vector<Number2, MemorySpaceType> &src)
    {
      AssertDimension(partitioner->local_size(), src.partitioner->local_size());
      if (partitioner->local_size() > 0)
        {
          dealii::internal::VectorOperations::
            functions<Number, Number2, MemorySpaceType>::copy(
              thread_loop_partitioner,
              partitioner->local_size(),
              src.data,
              data);
        }
    }



    template <typename Number, typename MemorySpaceType>
    template <typename MemorySpaceType2>
    void
    Vector<Number, MemorySpaceType>::import(
      const Vector<Number, MemorySpaceType2> &src,
      VectorOperation::values                 operation)
    {
      Assert(src.partitioner.get() != nullptr, ExcNotInitialized());
      Assert(partitioner->locally_owned_range() ==
               src.partitioner->locally_owned_range(),
             ExcMessage("Locally owned indices should be identical."));
      Assert(partitioner->ghost_indices() == src.partitioner->ghost_indices(),
             ExcMessage("Ghost indices should be identical."));
      ::dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::import(
          thread_loop_partitioner, allocated_size, operation, src.data, data);
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
      update_ghost_values_start();
      update_ghost_values_finish();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::zero_out_ghosts() const
    {
      if (data.values != nullptr)
        std::fill_n(data.values.get() + partitioner->local_size(),
                    partitioner->n_ghost_indices(),
                    Number());
#ifdef DEAL_II_COMPILER_CUDA_AWARE
      if (data.values_dev != nullptr)
        {
          const cudaError_t cuda_error_code =
            cudaMemset(data.values_dev.get() + partitioner->local_size(),
                       0,
                       partitioner->n_ghost_indices() * sizeof(Number));
          AssertCuda(cuda_error_code);
        }
#endif

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_start(
      const unsigned int                communication_channel,
      ::dealii::VectorOperation::values operation)
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
#  if defined(DEAL_II_COMPILER_CUDA_AWARE) && \
    defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
          if (std::is_same<MemorySpaceType, dealii::MemorySpace::CUDA>::value)
            {
              if (import_data.values_dev == nullptr)
                import_data.values_dev.reset(
                  Utilities::CUDA::allocate_device_data<Number>(
                    partitioner->n_import_indices()));
            }
          else
#  endif
            {
#  if !defined(DEAL_II_COMPILER_CUDA_AWARE) && \
    defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
              static_assert(
                std::is_same<MemorySpaceType, dealii::MemorySpace::Host>::value,
                "This code path should only be compiled for CUDA-aware-MPI for MemorySpace::Host!");
#  endif
              if (import_data.values == nullptr)
                {
                  Number *new_val;
                  Utilities::System::posix_memalign(
                    reinterpret_cast<void **>(&new_val),
                    64,
                    sizeof(Number) * partitioner->n_import_indices());
                  import_data.values.reset(new_val);
                }
            }
        }

#  if defined DEAL_II_COMPILER_CUDA_AWARE && \
    !defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
      if (std::is_same<MemorySpaceType, dealii::MemorySpace::CUDA>::value)
        {
          // Move the data to the host and then move it back to the
          // device. We use values to store the elements because the function
          // uses a view of the array and thus we need the data on the host to
          // outlive the scope of the function.
          Number *new_val;
          Utilities::System::posix_memalign(reinterpret_cast<void **>(&new_val),
                                            64,
                                            sizeof(Number) * allocated_size);

          data.values.reset(new_val);

          cudaError_t cuda_error_code =
            cudaMemcpy(data.values.get(),
                       data.values_dev.get(),
                       allocated_size * sizeof(Number),
                       cudaMemcpyDeviceToHost);
          AssertCuda(cuda_error_code);
        }
#  endif

#  if defined(DEAL_II_COMPILER_CUDA_AWARE) && \
    defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
      if (std::is_same<MemorySpaceType, dealii::MemorySpace::CUDA>::value)
        {
          partitioner->import_from_ghosted_array_start(
            operation,
            communication_channel,
            ArrayView<Number, MemorySpace::CUDA>(
              data.values_dev.get() + partitioner->local_size(),
              partitioner->n_ghost_indices()),
            ArrayView<Number, MemorySpace::CUDA>(
              import_data.values_dev.get(), partitioner->n_import_indices()),
            compress_requests);
        }
      else
#  endif
        {
          partitioner->import_from_ghosted_array_start(
            operation,
            communication_channel,
            ArrayView<Number, MemorySpace::Host>(
              data.values.get() + partitioner->local_size(),
              partitioner->n_ghost_indices()),
            ArrayView<Number, MemorySpace::Host>(
              import_data.values.get(), partitioner->n_import_indices()),
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
      ::dealii::VectorOperation::values operation)
    {
#ifdef DEAL_II_WITH_MPI
      vector_is_ghosted = false;

      // in order to zero ghost part of the vector, we need to call
      // import_from_ghosted_array_finish() regardless of
      // compress_requests.size() == 0

      // make this function thread safe
      std::lock_guard<std::mutex> lock(mutex);
#  if defined(DEAL_II_COMPILER_CUDA_AWARE) && \
    defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
      if (std::is_same<MemorySpaceType, MemorySpace::CUDA>::value)
        {
          Assert(partitioner->n_import_indices() == 0 ||
                   import_data.values_dev != nullptr,
                 ExcNotInitialized());
          partitioner
            ->import_from_ghosted_array_finish<Number, MemorySpace::CUDA>(
              operation,
              ArrayView<const Number, MemorySpace::CUDA>(
                import_data.values_dev.get(), partitioner->n_import_indices()),
              ArrayView<Number, MemorySpace::CUDA>(data.values_dev.get(),
                                                   partitioner->local_size()),
              ArrayView<Number, MemorySpace::CUDA>(
                data.values_dev.get() + partitioner->local_size(),
                partitioner->n_ghost_indices()),
              compress_requests);
        }
      else
#  endif
        {
          Assert(partitioner->n_import_indices() == 0 ||
                   import_data.values != nullptr,
                 ExcNotInitialized());
          partitioner
            ->import_from_ghosted_array_finish<Number, MemorySpace::Host>(
              operation,
              ArrayView<const Number, MemorySpace::Host>(
                import_data.values.get(), partitioner->n_import_indices()),
              ArrayView<Number, MemorySpace::Host>(data.values.get(),
                                                   partitioner->local_size()),
              ArrayView<Number, MemorySpace::Host>(
                data.values.get() + partitioner->local_size(),
                partitioner->n_ghost_indices()),
              compress_requests);
        }

#  if defined DEAL_II_COMPILER_CUDA_AWARE && \
    !defined  DEAL_II_MPI_WITH_CUDA_SUPPORT
      // The communication is done on the host, so we need to
      // move the data back to the device.
      if (std::is_same<MemorySpaceType, MemorySpace::CUDA>::value)
        {
          cudaError_t cuda_error_code =
            cudaMemcpy(data.values_dev.get(),
                       data.values.get(),
                       allocated_size * sizeof(Number),
                       cudaMemcpyHostToDevice);
          AssertCuda(cuda_error_code);

          data.values.reset();
        }
#  endif
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
#  if defined(DEAL_II_COMPILER_CUDA_AWARE) && \
    defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
          Assert(
            (std::is_same<MemorySpaceType, dealii::MemorySpace::CUDA>::value),
            ExcMessage(
              "Using MemorySpace::CUDA only allowed if the code is compiled with a CUDA compiler!"));
          if (import_data.values_dev == nullptr)
            import_data.values_dev.reset(
              Utilities::CUDA::allocate_device_data<Number>(
                partitioner->n_import_indices()));
#  else
#    ifdef DEAL_II_MPI_WITH_CUDA_SUPPORT
          static_assert(
            std::is_same<MemorySpaceType, dealii::MemorySpace::Host>::value,
            "This code path should only be compiled for CUDA-aware-MPI for MemorySpace::Host!");
#    endif
          if (import_data.values == nullptr)
            {
              Number *new_val;
              Utilities::System::posix_memalign(
                reinterpret_cast<void **>(&new_val),
                64,
                sizeof(Number) * partitioner->n_import_indices());
              import_data.values.reset(new_val);
            }
#  endif
        }

#  if defined DEAL_II_COMPILER_CUDA_AWARE && \
    !defined(DEAL_II_MPI_WITH_CUDA_SUPPORT)
      // Move the data to the host and then move it back to the
      // device. We use values to store the elements because the function
      // uses a view of the array and thus we need the data on the host to
      // outlive the scope of the function.
      Number *new_val;
      Utilities::System::posix_memalign(reinterpret_cast<void **>(&new_val),
                                        64,
                                        sizeof(Number) * allocated_size);

      data.values.reset(new_val);

      cudaError_t cuda_error_code = cudaMemcpy(data.values.get(),
                                               data.values_dev.get(),
                                               allocated_size * sizeof(Number),
                                               cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
#  endif

#  if !(defined(DEAL_II_COMPILER_CUDA_AWARE) && \
        defined(DEAL_II_MPI_WITH_CUDA_SUPPORT))
      partitioner->export_to_ghosted_array_start<Number, MemorySpace::Host>(
        communication_channel,
        ArrayView<const Number, MemorySpace::Host>(data.values.get(),
                                                   partitioner->local_size()),
        ArrayView<Number, MemorySpace::Host>(import_data.values.get(),
                                             partitioner->n_import_indices()),
        ArrayView<Number, MemorySpace::Host>(data.values.get() +
                                               partitioner->local_size(),
                                             partitioner->n_ghost_indices()),
        update_ghost_values_requests);
#  else
      partitioner->export_to_ghosted_array_start<Number, MemorySpace::CUDA>(
        communication_channel,
        ArrayView<const Number, MemorySpace::CUDA>(data.values_dev.get(),
                                                   partitioner->local_size()),
        ArrayView<Number, MemorySpace::CUDA>(import_data.values_dev.get(),
                                             partitioner->n_import_indices()),
        ArrayView<Number, MemorySpace::CUDA>(data.values_dev.get() +
                                               partitioner->local_size(),
                                             partitioner->n_ghost_indices()),
        update_ghost_values_requests);
#  endif

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

#  if !(defined(DEAL_II_COMPILER_CUDA_AWARE) && \
        defined(DEAL_II_MPI_WITH_CUDA_SUPPORT))
          partitioner->export_to_ghosted_array_finish(
            ArrayView<Number, MemorySpace::Host>(
              data.values.get() + partitioner->local_size(),
              partitioner->n_ghost_indices()),
            update_ghost_values_requests);
#  else
          partitioner->export_to_ghosted_array_finish(
            ArrayView<Number, MemorySpace::CUDA>(
              data.values_dev.get() + partitioner->local_size(),
              partitioner->n_ghost_indices()),
            update_ghost_values_requests);
#  endif
        }

#  if defined DEAL_II_COMPILER_CUDA_AWARE && \
    !defined  DEAL_II_MPI_WITH_CUDA_SUPPORT
      // The communication is done on the host, so we need to
      // move the data back to the device.
      if (std::is_same<MemorySpaceType, MemorySpace::CUDA>::value)
        {
          cudaError_t cuda_error_code =
            cudaMemcpy(data.values_dev.get() + partitioner->local_size(),
                       data.values.get() + partitioner->local_size(),
                       partitioner->n_ghost_indices() * sizeof(Number),
                       cudaMemcpyHostToDevice);
          AssertCuda(cuda_error_code);

          data.values.reset();
        }
#  endif

#endif
      vector_is_ghosted = true;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::import(
      const ReadWriteVector<Number> &                 V,
      VectorOperation::values                         operation,
      std::shared_ptr<const CommunicationPatternBase> communication_pattern)
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

      data.copy_to(tmp_vector.begin(), local_size());

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

      data.copy_from(tmp_vector.begin(), local_size());
    }

    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::swap(Vector<Number, MemorySpaceType> &v)
    {
#ifdef DEAL_II_WITH_MPI

#  ifdef DEBUG
      if (Utilities::MPI::job_supports_mpi())
        {
          // make sure that there are not outstanding requests from updating
          // ghost values or compress
          int flag = 1;
          if (update_ghost_values_requests.size() > 0)
            {
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
#  endif

      std::swap(compress_requests, v.compress_requests);
      std::swap(update_ghost_values_requests, v.update_ghost_values_requests);
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
      if (s == Number())
        zero_out_ghosts();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(const VectorSpaceVector<Number> &V,
                                            const bool omit_zeroing_entries)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &down_V = dynamic_cast<const VectorType &>(V);

      reinit(down_V, omit_zeroing_entries);
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator+=(const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);

      AssertDimension(local_size(), v.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_vector(
          thread_loop_partitioner, partitioner->local_size(), v.data, data);

      if (vector_is_ghosted)
        update_ghost_values();

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
      AssertIsFinite(a);

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_factor(
          thread_loop_partitioner, partitioner->local_size(), a, data);

      if (vector_is_ghosted)
        update_ghost_values();
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
      // Downcast. Throws an exception if invalid.
      using VectorType = Vector<Number, MemorySpaceType>;
      Assert(dynamic_cast<const VectorType *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &v = dynamic_cast<const VectorType &>(vv);
      Assert(dynamic_cast<const VectorType *>(&ww) != nullptr,
             ExcVectorTypeNotCompatible());
      const VectorType &w = dynamic_cast<const VectorType &>(ww);

      AssertIsFinite(a);
      AssertIsFinite(b);

      AssertDimension(local_size(), v.local_size());
      AssertDimension(local_size(), w.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::add_avpbw(
          thread_loop_partitioner,
          partitioner->local_size(),
          a,
          b,
          v.data,
          w.data,
          data);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const std::vector<size_type> &indices,
                                         const std::vector<Number> &   values)
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
      AssertDimension(local_size(), v.local_size());

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::sadd_xv(
          thread_loop_partitioner, partitioner->local_size(), x, v.data, data);

      if (vector_is_ghosted)
        update_ghost_values();
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
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator*=(const Number factor)
    {
      AssertIsFinite(factor);

      dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::multiply_factor(
          thread_loop_partitioner, partitioner->local_size(), factor, data);

      if (vector_is_ghosted)
        update_ghost_values();

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
      Assert(size() != 0, ExcEmptyObject());

      if (partitioner->local_size() == 0)
        return Number();

      Number sum = ::dealii::internal::VectorOperations::
        functions<Number, Number, MemorySpaceType>::mean_value(
          thread_loop_partitioner, partitioner->local_size(), data);

      return sum / real_type(partitioner->local_size());
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value() const
    {
      Number local_result = mean_value_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result * static_cast<real_type>(
                                                    partitioner->local_size()),
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
          thread_loop_partitioner, partitioner->local_size(), sum, data);

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
          thread_loop_partitioner, partitioner->local_size(), sum, p, data);

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
      if (import_data.values != nullptr || import_data.values_dev != nullptr)
        memory += (static_cast<std::size_t>(partitioner->n_import_indices()) *
                   sizeof(Number));
      return memory;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::print(std::ostream &     out,
                                           const unsigned int precision,
                                           const bool         scientific,
                                           const bool         across) const
    {
      Assert(partitioner.get() != nullptr, ExcInternalError());
      AssertThrow(out, ExcIO());
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
        for (unsigned int i = 0; i < partitioner->this_mpi_process(); i++)
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
        for (size_type i = 0; i < partitioner->local_size(); ++i)
          out << stored_elements[i] << ' ';
      else
        for (size_type i = 0; i < partitioner->local_size(); ++i)
          out << stored_elements[i] << std::endl;
      out << std::endl;

      if (vector_is_ghosted)
        {
          out << "Ghost entries (global index / value):" << std::endl;
          if (across)
            for (size_type i = 0; i < partitioner->n_ghost_indices(); ++i)
              out << '(' << partitioner->ghost_indices().nth_index_in_set(i)
                  << '/' << stored_elements[partitioner->local_size() + i]
                  << ") ";
          else
            for (size_type i = 0; i < partitioner->n_ghost_indices(); ++i)
              out << '(' << partitioner->ghost_indices().nth_index_in_set(i)
                  << '/' << stored_elements[partitioner->local_size() + i]
                  << ")" << std::endl;
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

      AssertThrow(out, ExcIO());
      // reset output format
      out.flags(old_flags);
      out.precision(old_precision);
    }

  } // end of namespace distributed
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
