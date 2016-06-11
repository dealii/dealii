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

#ifndef dealii__mpi_templates_h
#define dealii__mpi_templates_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/vector.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
#ifdef DEAL_II_WITH_MPI
      /**
       * Return the corresponding MPI data type id for the argument given.
       */
      MPI_Datatype mpi_type_id (const int *)
      {
        return MPI_INT;
      }


      MPI_Datatype mpi_type_id (const long int *)
      {
        return MPI_LONG;
      }


      MPI_Datatype mpi_type_id (const unsigned int *)
      {
        return MPI_UNSIGNED;
      }


      MPI_Datatype mpi_type_id (const unsigned long int *)
      {
        return MPI_UNSIGNED_LONG;
      }


      MPI_Datatype mpi_type_id (const unsigned long long int *)
      {
        return MPI_UNSIGNED_LONG_LONG;
      }


      MPI_Datatype mpi_type_id (const float *)
      {
        return MPI_FLOAT;
      }


      MPI_Datatype mpi_type_id (const double *)
      {
        return MPI_DOUBLE;
      }


      MPI_Datatype mpi_type_id (const long double *)
      {
        return MPI_LONG_DOUBLE;
      }
#endif

      template <typename T>
      void all_reduce (const MPI_Op      &mpi_op,
                       const T *const    values,
                       const MPI_Comm    &mpi_communicator,
                       T                 *output,
                       const std::size_t  size)
      {
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
            MPI_Allreduce (values != output
                           ?
                           // TODO This const_cast is only needed for older
                           // (e.g., openMPI 1.6, released in 2012)
                           // implementations of MPI-2. It is not needed as of
                           // MPI-3 and we should remove it at some point in
                           // the future.
                           const_cast<void *>(static_cast<const void *>(values))
                           :
                           MPI_IN_PLACE,
                           static_cast<void *>(output),
                           static_cast<int>(size),
                           internal::mpi_type_id(values),
                           mpi_op,
                           mpi_communicator);
          }
        else
#endif
          {
            (void)mpi_op;
            (void)mpi_communicator;
            for (std::size_t i=0; i<size; ++i)
              output[i] = values[i];
          }
      }

      template <typename T>
      void all_reduce (const MPI_Op                   &mpi_op,
                       const std::complex<T> *const    values,
                       const MPI_Comm                 &mpi_communicator,
                       std::complex<T>                *output,
                       const std::size_t               size)
      {
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
            T dummy_selector;
            MPI_Allreduce (values != output
                           ?
                           // TODO This const_cast is only needed for older
                           // (e.g., openMPI 1.6, released in 2012)
                           // implementations of MPI-2. It is not needed as of
                           // MPI-3 and we should remove it at some point in
                           // the future.
                           const_cast<void *>(static_cast<const void *>(values))
                           :
                           MPI_IN_PLACE,
                           static_cast<void *>(output),
                           static_cast<int>(size*2),
                           internal::mpi_type_id(&dummy_selector),
                           mpi_op,
                           mpi_communicator);
          }
        else
#endif
          {
            (void)mpi_op;
            (void)mpi_communicator;
            for (std::size_t i=0; i<size; ++i)
              output[i] = values[i];
          }
      }

      template <typename T>
      T all_reduce (const MPI_Op &mpi_op,
                    const T &t,
                    const MPI_Comm &mpi_communicator)
      {
        T output;
        all_reduce(mpi_op, &t, mpi_communicator, &output, 1);
        return output;
      }

      template <typename T>
      void all_reduce (const MPI_Op &mpi_op,
                       const std::vector<T> &values,
                       const MPI_Comm       &mpi_communicator,
                       std::vector<T>       &output)
      {
        Assert(values.size() == output.size(),
               ExcDimensionMismatch(values.size(), output.size()));
        all_reduce(mpi_op, &values[0], mpi_communicator, &output[0], values.size());
      }

      template <typename T>
      void all_reduce (const MPI_Op    &mpi_op,
                       const Vector<T> &values,
                       const MPI_Comm  &mpi_communicator,
                       Vector<T>  &output)
      {
        Assert(values.size() == output.size(),
               ExcDimensionMismatch(values.size(), output.size()));
        all_reduce(mpi_op, values.begin(), mpi_communicator, output.begin(), values.size());
      }
    }


    template <typename T>
    T sum (const T &t,
           const MPI_Comm &mpi_communicator)
    {
      return internal::all_reduce(MPI_SUM, t, mpi_communicator);
    }


    template <typename T>
    void sum (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &sums)
    {
      internal::all_reduce(MPI_SUM, values, mpi_communicator, sums);
    }

    template <typename T>
    void sum (const Vector<T> &values,
              const MPI_Comm &mpi_communicator,
              Vector<T> &sums)
    {
      internal::all_reduce(MPI_SUM, values, mpi_communicator, sums);
    }


    template <int rank, int dim, typename Number>
    Tensor<rank,dim,Number>
    sum (const Tensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator)
    {
      const unsigned int n_entries = Tensor<rank,dim,Number>::n_independent_components;
      Number entries[ Tensor<rank,dim,Number>::n_independent_components ];

      for (unsigned int i=0; i< n_entries; ++i)
        entries[i] = local[ local.unrolled_to_component_indices(i) ];

      Number global_entries[ Tensor<rank,dim,Number>::n_independent_components ];
      dealii::Utilities::MPI::sum( entries, mpi_communicator, global_entries );

      Tensor<rank,dim,Number> global;
      for (unsigned int i=0; i< n_entries; ++i)
        global[ global.unrolled_to_component_indices(i) ] = global_entries[i];

      return global;
    }

    template <int rank, int dim, typename Number>
    SymmetricTensor<rank,dim,Number>
    sum (const SymmetricTensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator)
    {
      const unsigned int n_entries = SymmetricTensor<rank,dim,Number>::n_independent_components;
      Number entries[ SymmetricTensor<rank,dim,Number>::n_independent_components ];

      for (unsigned int i=0; i< n_entries; ++i)
        entries[i] = local[ local.unrolled_to_component_indices(i) ];

      Number global_entries[ SymmetricTensor<rank,dim,Number>::n_independent_components ];
      dealii::Utilities::MPI::sum( entries, mpi_communicator, global_entries );

      SymmetricTensor<rank,dim,Number> global;
      for (unsigned int i=0; i< n_entries; ++i)
        global[ global.unrolled_to_component_indices(i) ] = global_entries[i];

      return global;
    }

    template <typename T>
    T max (const T &t,
           const MPI_Comm &mpi_communicator)
    {
      return internal::all_reduce(MPI_MAX, t, mpi_communicator);
    }


    template <typename T>
    void max (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &maxima)
    {
      internal::all_reduce(MPI_MAX, values, mpi_communicator, maxima);
    }


    template <typename T>
    T min (const T &t,
           const MPI_Comm &mpi_communicator)
    {
      return internal::all_reduce(MPI_MIN, t, mpi_communicator);
    }


    template <typename T>
    void min (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &minima)
    {
      internal::all_reduce(MPI_MIN, values, mpi_communicator, minima);
    }
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
