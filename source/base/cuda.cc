// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
