// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/lac/tensor_product_matrix.templates.h>

#include <array>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TensorProductMatrixSymmetricSum
  {
#include "lac/tensor_product_matrix.inst"

  }
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
