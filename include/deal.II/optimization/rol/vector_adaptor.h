// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_optimization_rol_vector_adaptor_h
#define dealii_optimization_rol_vector_adaptor_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_ROL
#  include <deal.II/trilinos/rol_adaptor.h>


DEAL_II_NAMESPACE_OPEN

namespace Rol
{
  /**
   * @copydoc TrilinosWrappers::ROLAdaptor
   *
   * @deprecated Use TrilinosWrappers::ROLAdaptor instead.
   */
  template <typename VectorType>
  using VectorAdaptor DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use dealii::TrilinosWrappers::ROLAdaptor instead.") =
    dealii::TrilinosWrappers::ROLAdaptor<VectorType>;
} // namespace Rol

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_ROL

#endif // dealii_optimization_rol_vector_adaptor_h
