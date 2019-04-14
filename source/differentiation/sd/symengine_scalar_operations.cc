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


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    /* ------------------------- Symbol creation -----------------------*/


    Expression
    make_symbol(const std::string &symbol)
    {
      return Expression(symbol);
    }


    Expression
    make_symbolic_function(const std::string &             symbol,
                           const SD::types::symbol_vector &arguments)
    {
      return Expression(symbol, arguments);
    }


    Expression
    make_symbolic_function(const std::string &                symbol,
                           const SD::types::substitution_map &arguments)
    {
      return make_symbolic_function(symbol,
                                    SD::Utilities::extract_symbols(arguments));
    }


  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
