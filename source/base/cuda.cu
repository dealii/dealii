// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/vector_memory.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace CUDA
  {
    Handle::Handle()
    {
      cusolverStatus_t cusolver_error_code =
        cusolverDnCreate(&cusolver_dn_handle);
      AssertCusolver(cusolver_error_code);

      cusolver_error_code = cusolverSpCreate(&cusolver_sp_handle);
      AssertCusolver(cusolver_error_code);

      cusparseStatus_t cusparse_error_code = cusparseCreate(&cusparse_handle);
      AssertCusparse(cusparse_error_code);
    }



    Handle::~Handle()
    {
      cusolverStatus_t cusolver_error_code =
        cusolverDnDestroy(cusolver_dn_handle);
      AssertCusolver(cusolver_error_code);

      cusolver_error_code = cusolverSpDestroy(cusolver_sp_handle);
      AssertCusolver(cusolver_error_code);

      cusparseStatus_t cusparse_error_code = cusparseDestroy(cusparse_handle);
      AssertCusparse(cusparse_error_code);
    }
  } // namespace CUDA
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
