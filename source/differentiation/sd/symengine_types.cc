// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_types.h>

#  include <symengine/basic.h>
#  include <symengine/dict.h>


DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;

    namespace types
    {
      namespace internal
      {
        bool
        ExpressionKeyLess::operator()(const SD::Expression &lhs,
                                      const SD::Expression &rhs) const
        {
          return SE::RCPBasicKeyLess()(lhs.get_RCP(), rhs.get_RCP());
        }
      } // namespace internal

    } // namespace types

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
