// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
