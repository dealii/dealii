// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_differentiation_sd_symengine_utilities_h
#define dealii_differentiation_sd_symengine_utilities_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_types.h>

#  include <symengine/basic.h>
#  include <symengine/dict.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;

    namespace Utilities
    {
      /**
       * Convert a map of Expressions to its SymEngine counterpart.
       */
      SE::map_basic_basic
      convert_expression_map_to_basic_map(
        const SD::types::substitution_map &substitution_map);

      /**
       * Convert a vector of Expressions to its SymEngine counterpart.
       */
      SE::vec_basic
      convert_expression_vector_to_basic_vector(
        const SD::types::symbol_vector &symbol_vector);

    } // namespace Utilities

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_utilities_h
