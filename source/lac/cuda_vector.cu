// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#include <deal.II/base/cuda.h>
#include <deal.II/base/cuda_size.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_atomic.h>
#include <deal.II/lac/cuda_kernels.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <cmath>

#ifdef DEAL_II_WITH_CUDA

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    using ::dealii::CUDAWrappers::block_size;
    using ::dealii::CUDAWrappers::chunk_size;

    template <typename Number>
    Vector<Number>::Vector()
      : val(nullptr)
      , n_elements(0)
    {}



    template <typename Number>
    Vector<Number>::Vector(const Vector<Number> &V)
      : n_elements(V.n_elements)
    {
      // Allocate the memory
      cudaError_t error_code = cudaMalloc(&val, n_elements * sizeof(Number));
      AssertCuda(error_code);
      // Copy the values.
      error_code = cudaMemcpy(val,
                              V.val,
                              n_elements * sizeof(Number),
                              cudaMemcpyDeviceToDevice);
      AssertCuda(error_code);
    }



    template <typename Number>
    Vector<Number>::Vector(const size_type n)
      : val(nullptr)
      , n_elements(0)
    {
      reinit(n, false);
    }



    template <typename Number>
    Vector<Number>::~Vector()
    {
      if (val != nullptr)
        {
          cudaError_t error_code = cudaFree(val);
          AssertCuda(error_code);
          val        = nullptr;
          n_elements = 0;
        }
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const size_type n, const bool omit_zeroing_entries)
    {
      // Resize the underlying array if necessary
      if (n == 0)
        {
          if (val != nullptr)
            {
              cudaError_t error_code = cudaFree(val);
              AssertCuda(error_code);
              val = nullptr;
            }
        }
      else
        {
          if ((n_elements != n) && (val != nullptr))
            {
              cudaError_t error_code = cudaFree(val);
              AssertCuda(error_code);
            }

          cudaError_t error_code = cudaMalloc(&val, n * sizeof(Number));
          AssertCuda(error_code);

          // If necessary set the elements to zero
          if (omit_zeroing_entries == false)
            {
              cudaError_t error_code = cudaMemset(val, 0, n * sizeof(Number));
              AssertCuda(error_code);
            }
        }
      n_elements = n;
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const VectorSpaceVector<Number> &V,
                           const bool omit_zeroing_entries)
    {
      reinit(V.size(), omit_zeroing_entries);
    }



    template <typename Number>
    void
    Vector<Number>::import(const ReadWriteVector<Number> &V,
                           VectorOperation::values        operation,
                           std::shared_ptr<const CommunicationPatternBase>)
    {
      if (operation == VectorOperation::insert)
        {
          cudaError_t error_code = cudaMemcpy(val,
                                              V.begin(),
                                              n_elements * sizeof(Number),
                                              cudaMemcpyHostToDevice);
          AssertCuda(error_code);
        }
      else if (operation == VectorOperation::add)
        {
          // Create a temporary vector on the device
          Number *    tmp;
          cudaError_t error_code =
            cudaMalloc(&tmp, n_elements * sizeof(Number));
          AssertCuda(error_code);

          // Copy the vector from the host to the temporary vector on the device
          error_code = cudaMemcpy(tmp,
                                  V.begin(),
                                  n_elements * sizeof(Number),
                                  cudaMemcpyHostToDevice);
          AssertCuda(error_code);

          // Add the two vectors
          const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);

          kernel::vector_bin_op<Number, kernel::Binop_Addition>
            <<<n_blocks, block_size>>>(val, tmp, n_elements);
          // Check that the kernel was launched correctly
          AssertCuda(cudaGetLastError());
          // Check that there was no problem during the execution of the kernel
          AssertCuda(cudaDeviceSynchronize());

          // Delete the temporary vector
          error_code = cudaFree(tmp);
          AssertCuda(error_code);
        }
      else
        AssertThrow(false, ExcNotImplemented());
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator=(const Number s)
    {
      Assert(s == Number(), ExcMessage("Only 0 can be assigned to a vector."));
      (void)s;

      cudaError_t error_code = cudaMemset(val, 0, n_elements * sizeof(Number));
      AssertCuda(error_code);

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator*=(const Number factor)
    {
      AssertIsFinite(factor);
      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::vec_scale<Number>
        <<<n_blocks, block_size>>>(val, factor, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator/=(const Number factor)
    {
      AssertIsFinite(factor);
      Assert(factor != Number(0.), ExcZero());
      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::vec_scale<Number>
        <<<n_blocks, block_size>>>(val, 1. / factor, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator+=(const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If it fails, it throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements"));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);

      kernel::vector_bin_op<Number, kernel::Binop_Addition>
        <<<n_blocks, block_size>>>(val, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator-=(const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);

      kernel::vector_bin_op<Number, kernel::Binop_Subtraction>
        <<<n_blocks, block_size>>>(val, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Number Vector<Number>::operator*(const VectorSpaceVector<Number> &V) const
    {
      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements"));

      Number *    result_device;
      cudaError_t error_code =
        cudaMalloc(&result_device, n_elements * sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::double_vector_reduction<Number, kernel::DotProduct<Number>>
        <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                  val,
                                                  down_V.val,
                                                  static_cast<unsigned int>(
                                                    n_elements));

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result,
                              result_device,
                              sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result;
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number a)
    {
      AssertIsFinite(a);
      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::vec_add<Number><<<n_blocks, block_size>>>(val, a, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V)
    {
      AssertIsFinite(a);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::add_aV<Number><<<dim3(n_blocks, 1), dim3(block_size)>>>(
        val, a, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number                     a,
                        const VectorSpaceVector<Number> &V,
                        const Number                     b,
                        const VectorSpaceVector<Number> &W)
    {
      AssertIsFinite(a);
      AssertIsFinite(b);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements."));

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&W) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_W = dynamic_cast<const Vector<Number> &>(W);
      Assert(down_W.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::add_aVbW<Number><<<dim3(n_blocks, 1), dim3(block_size)>>>(
        val, a, down_V.val, b, down_W.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void
    Vector<Number>::sadd(const Number                     s,
                         const Number                     a,
                         const VectorSpaceVector<Number> &V)
    {
      AssertIsFinite(s);
      AssertIsFinite(a);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage(
               "Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::sadd<Number><<<dim3(n_blocks, 1), dim3(block_size)>>>(
        s, val, a, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void
    Vector<Number>::scale(const VectorSpaceVector<Number> &scaling_factors)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&scaling_factors) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_scaling_factors =
        dynamic_cast<const Vector<Number> &>(scaling_factors);
      Assert(down_scaling_factors.size() == this->size(),
             ExcMessage(
               "Cannot scale two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::scale<Number>
        <<<dim3(n_blocks, 1), dim3(block_size)>>>(val,
                                                  down_scaling_factors.val,
                                                  n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void
    Vector<Number>::equ(const Number a, const VectorSpaceVector<Number> &V)
    {
      AssertIsFinite(a);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(
        down_V.size() == this->size(),
        ExcMessage(
          "Cannot assign two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::equ<Number><<<dim3(n_blocks, 1), dim3(block_size)>>>(val,
                                                                   a,
                                                                   down_V.val,
                                                                   n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    bool
    Vector<Number>::all_zero() const
    {
      return (linfty_norm() == 0) ? true : false;
    }



    template <typename Number>
    typename Vector<Number>::value_type
    Vector<Number>::mean_value() const
    {
      Number *    result_device;
      cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::reduction<Number, kernel::ElemSum<Number>>
        <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                  val,
                                                  n_elements);

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result,
                              result_device,
                              sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result /
             static_cast<typename Vector<Number>::value_type>(n_elements);
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm() const
    {
      Number *    result_device;
      cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::reduction<Number, kernel::L1Norm<Number>>
        <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                  val,
                                                  n_elements);

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result,
                              result_device,
                              sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result;
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l2_norm() const
    {
      return std::sqrt((*this) * (*this));
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm() const
    {
      Number *    result_device;
      cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::reduction<Number, kernel::LInfty<Number>>
        <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                  val,
                                                  n_elements);

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result,
                              result_device,
                              sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result;
    }



    template <typename Number>
    Number
    Vector<Number>::add_and_dot(const Number                     a,
                                const VectorSpaceVector<Number> &V,
                                const VectorSpaceVector<Number> &W)
    {
      AssertIsFinite(a);

      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());
      Assert(dynamic_cast<const Vector<Number> *>(&W) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V and W. If it fails, throw an exceptiion.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage("Vector V has the wrong size."));
      const Vector<Number> &down_W = dynamic_cast<const Vector<Number> &>(W);
      Assert(down_W.size() == this->size(),
             ExcMessage("Vector W has the wrong size."));

      Number *    res_d;
      cudaError_t error_code = cudaMalloc(&res_d, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(res_d, 0., sizeof(Number));
      AssertCuda(error_code);

      const int n_blocks = 1 + (n_elements - 1) / (chunk_size * block_size);
      kernel::add_and_dot<Number><<<dim3(n_blocks, 1), dim3(block_size)>>>(
        res_d, val, down_V.val, down_W.val, a, n_elements);

      Number res;
      error_code =
        cudaMemcpy(&res, res_d, sizeof(Number), cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      error_code = cudaFree(res_d);

      return res;
    }



    template <typename Number>
    void
    Vector<Number>::print(std::ostream &     out,
                          const unsigned int precision,
                          const bool         scientific,
                          const bool) const
    {
      AssertThrow(out, ExcIO());
      std::ios::fmtflags old_flags     = out.flags();
      unsigned int       old_precision = out.precision(precision);

      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      out << "IndexSet: ";
      complete_index_set(n_elements).print(out);
      out << std::endl;

      // Copy the vector to the host
      std::vector<Number> cpu_val(n_elements);
      Utilities::CUDA::copy_to_host(val, cpu_val);
      for (unsigned int i = 0; i < n_elements; ++i)
        out << cpu_val[i] << std::endl;
      out << std::flush;

      AssertThrow(out, ExcIO());
      // reset output format
      out.flags(old_flags);
      out.precision(old_precision);
    }



    template <typename Number>
    std::size_t
    Vector<Number>::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      memory += sizeof(Number) * static_cast<std::size_t>(n_elements);

      return memory;
    }



    // Explicit Instanationation
    template class Vector<float>;
    template class Vector<double>;
  } // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
