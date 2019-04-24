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

#  include <symengine/real_double.h>

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


    /* --------------------------- Differentiation ------------------------- */


    Expression
    differentiate(const Expression &func, const Expression &op)
    {
      return func.differentiate(op);
    }


    /* ------------------------ Symbolic map creation ----------------------*/


    namespace internal
    {
      bool
      is_valid_substitution_symbol(const SE::Basic &entry)
      {
        // Allow substitution of a symbol
        // It is pretty clear as to why this is wanted...
        if (SE::is_a<SE::Symbol>(entry))
          return true;

        // Allow substitution of a function symbol
        // If desired, we can transform general but undefined functional
        // relationships to an explicit form that is concrete. This is
        // required for a symbolic expression to be parsed by a Lambda or LLVM
        // optimizer.
        if (SE::is_a<SE::FunctionSymbol>(entry))
          return true;

        // Allow substitution of the explicit expression of the derivative of
        // an implicitly defined symbol (i.e. the result of the derivative of
        // a FunctionSymbol).
        if (SE::is_a<SE::Derivative>(entry))
          return true;

        // When performing tensor differentiation, one may end up with a
        // coefficient of one half due to symmetry operations, e.g.
        // 0.5*Derivative(Qi_00(C_11, C_00, C_01), C_01)
        // So we explicitly check for this exact case
        if (SE::is_a<SE::Mul>(entry))
          {
            const SE::Mul &entry_mul = SE::down_cast<const SE::Mul &>(entry);
            // Check that the factor is a half...
            if (SE::eq(*(entry_mul.get_coef()), *SE::real_double(0.5)))
              {
                // ...and that there is only one entry and that its a
                // Derivative type
                const SE::map_basic_basic &entry_mul_dict =
                  entry_mul.get_dict();
                if (entry_mul_dict.size() == 1 &&
                    SE::is_a<SE::Derivative>(*(entry_mul_dict.begin()->first)))
                  return true;
              }
          }

        return false;
      }

    } // namespace internal

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
