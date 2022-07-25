//-----------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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
//-----------------------------------------------------------

#include <deal.II/sundials/copy.h>

#ifdef DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_OPEN
namespace SUNDIALS
{
  namespace internal
  {
    /**
     * SUNDIALS provides different macros for getting the local length of a
     * vector for serial and parallel vectors (as well as various parallel
     * vectors that are not yet supported by deal.II). This function provides
     * a generic interface to both and does a (checked) conversion from long
     * int (the type SUNDIALS uses for lengths) to std::size_t.
     */
    inline std::size_t
    N_Vector_length(const N_Vector &vec)
    {
      const N_Vector_ID id     = N_VGetVectorID(vec);
      long int          length = -1;
      switch (id)
        {
          case SUNDIALS_NVEC_SERIAL:
            {
              length = NV_LENGTH_S(vec);
              break;
            }
#  ifdef DEAL_II_WITH_MPI
          case SUNDIALS_NVEC_PARALLEL:
            {
              length = NV_LOCLENGTH_P(vec);
              break;
            }
#  endif
          default:
            Assert(false, ExcNotImplemented());
        }

      Assert(length >= 0, ExcInternalError());
      return static_cast<std::size_t>(length);
    }


    void
    copy(LinearAlgebra::distributed::Vector<double> &dst, const N_Vector &src)
    {
      const IndexSet    is = dst.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(src));
      dst.zero_out_ghost_values();
      for (std::size_t i = 0; i < N; ++i)
        {
#  ifdef DEAL_II_WITH_MPI
          dst.local_element(i) = NV_Ith_P(src, i);
#  else
          dst.local_element(i)        = NV_Ith_S(src, i);
#  endif
        }
      dst.compress(VectorOperation::insert);
    }

    void
    copy(N_Vector &dst, const LinearAlgebra::distributed::Vector<double> &src)
    {
      const IndexSet    is = src.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(dst));
      for (std::size_t i = 0; i < N; ++i)
        {
#  ifdef DEAL_II_WITH_MPI
          NV_Ith_P(dst, i) = src.local_element(i);
#  else
          NV_Ith_S(dst, i)            = src.local_element(i);
#  endif
        }
    }

    void
    copy(LinearAlgebra::distributed::BlockVector<double> &dst,
         const N_Vector &                                 src)
    {
      const IndexSet    is = dst.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(src));
      dst.zero_out_ghost_values();
      for (std::size_t i = 0; i < N; ++i)
        {
#  ifdef DEAL_II_WITH_MPI
          dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
#  else
          dst[is.nth_index_in_set(i)] = NV_Ith_S(src, i);
#  endif
        }
      dst.compress(VectorOperation::insert);
    }

    void
    copy(N_Vector &                                             dst,
         const LinearAlgebra::distributed::BlockVector<double> &src)
    {
      IndexSet          is = src.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(dst));
      for (std::size_t i = 0; i < N; ++i)
        {
#  ifdef DEAL_II_WITH_MPI
          NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
#  else
          NV_Ith_S(dst, i)            = src[is.nth_index_in_set(i)];
#  endif
        }
    }

#  ifdef DEAL_II_WITH_MPI

#    ifdef DEAL_II_WITH_TRILINOS


    void
    copy(TrilinosWrappers::MPI::Vector &dst, const N_Vector &src)
    {
      const IndexSet    is = dst.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(src));
      for (std::size_t i = 0; i < N; ++i)
        {
          dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
        }
      dst.compress(VectorOperation::insert);
    }

    void
    copy(N_Vector &dst, const TrilinosWrappers::MPI::Vector &src)
    {
      const IndexSet    is = src.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(dst));
      for (std::size_t i = 0; i < N; ++i)
        {
          NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
        }
    }

    void
    copy(TrilinosWrappers::MPI::BlockVector &dst, const N_Vector &src)
    {
      const IndexSet    is = dst.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(src));
      for (std::size_t i = 0; i < N; ++i)
        {
          dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
        }
      dst.compress(VectorOperation::insert);
    }

    void
    copy(N_Vector &dst, const TrilinosWrappers::MPI::BlockVector &src)
    {
      IndexSet          is = src.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(dst));
      for (std::size_t i = 0; i < N; ++i)
        {
          NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
        }
    }

#    endif // DEAL_II_WITH_TRILINOS

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX

    void
    copy(PETScWrappers::MPI::Vector &dst, const N_Vector &src)
    {
      const IndexSet    is = dst.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(src));
      for (std::size_t i = 0; i < N; ++i)
        {
          dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
        }
      dst.compress(VectorOperation::insert);
    }

    void
    copy(N_Vector &dst, const PETScWrappers::MPI::Vector &src)
    {
      const IndexSet    is = src.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(dst));
      for (std::size_t i = 0; i < N; ++i)
        {
          NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
        }
    }

    void
    copy(PETScWrappers::MPI::BlockVector &dst, const N_Vector &src)
    {
      const IndexSet    is = dst.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(src));
      for (std::size_t i = 0; i < N; ++i)
        {
          dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
        }
      dst.compress(VectorOperation::insert);
    }

    void
    copy(N_Vector &dst, const PETScWrappers::MPI::BlockVector &src)
    {
      const IndexSet    is = src.locally_owned_elements();
      const std::size_t N  = is.n_elements();
      AssertDimension(N, N_Vector_length(dst));
      for (std::size_t i = 0; i < N; ++i)
        {
          NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
        }
    }

#      endif // PETSC_USE_COMPLEX
#    endif   // DEAL_II_WITH_PETSC

#  endif // mpi

    void
    copy(BlockVector<double> &dst, const N_Vector &src)
    {
      const std::size_t N = dst.size();
      AssertDimension(N_Vector_length(src), N);
      for (std::size_t i = 0; i < N; ++i)
        {
          dst[i] = NV_Ith_S(src, i);
        }
    }

    void
    copy(N_Vector &dst, const BlockVector<double> &src)
    {
      const std::size_t N = src.size();
      AssertDimension(N_Vector_length(dst), N);
      for (std::size_t i = 0; i < N; ++i)
        {
          NV_Ith_S(dst, i) = src[i];
        }
    }

    void
    copy(Vector<double> &dst, const N_Vector &src)
    {
      const std::size_t N = dst.size();
      AssertDimension(N_Vector_length(src), N);
      for (std::size_t i = 0; i < N; ++i)
        {
          dst[i] = NV_Ith_S(src, i);
        }
    }

    void
    copy(N_Vector &dst, const Vector<double> &src)
    {
      const std::size_t N = src.size();
      AssertDimension(N_Vector_length(dst), N);
      for (std::size_t i = 0; i < N; ++i)
        {
          NV_Ith_S(dst, i) = src[i];
        }
    }
  } // namespace internal
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#endif
