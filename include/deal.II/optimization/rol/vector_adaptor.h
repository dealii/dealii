// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_optimization_rol_vector_adaptor_h
#define dealii_optimization_rol_vector_adaptor_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_ROL
#  include <deal.II/trilinos/rol_adaptor.h>

#endif // DEAL_II_TRILINOS_WITH_ROL

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_TRILINOS_WITH_ROL
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

#endif // DEAL_II_TRILINOS_WITH_ROL

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_optimization_rol_vector_adaptor_h
