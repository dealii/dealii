// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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
    namespace SE = ::SymEngine;


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


    /* ---------------- Symbol map creation and manipulation --------------*/


#  ifndef DOXYGEN
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


      void
      set_value_in_symbol_map(
        types::substitution_map &                     substitution_map,
        const SymEngine::RCP<const SymEngine::Basic> &symbol,
        const SymEngine::RCP<const SymEngine::Basic> &value)
      {
        Assert(
          internal::is_valid_substitution_symbol(*symbol),
          ExcMessage(
            "Substitution with a number that does not represent a symbol or symbolic derivative"));

        auto it_sym = substitution_map.find(Expression(symbol));
        Assert(it_sym != substitution_map.end(),
               ExcMessage("Did not find this symbol in the map."));

        it_sym->second = Expression(value);
      }

    } // namespace internal


    void
    set_value_in_symbol_map(types::substitution_map &substitution_map,
                            const Expression &       symbol,
                            const Expression &       value)
    {
      internal::set_value_in_symbol_map(substitution_map,
                                        symbol.get_RCP(),
                                        value.get_RCP());
    }


    void
    set_value_in_symbol_map(types::substitution_map &      substitution_map,
                            const types::substitution_map &symbol_values)
    {
      for (const auto &entry : symbol_values)
        set_value_in_symbol_map(substitution_map, entry.first, entry.second);
    }
#  endif


    /* ------------------ Symbol substitution map creation ----------------*/


    types::substitution_map
    make_substitution_map(const Expression &symbol, const Expression &value)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbol, value);
      return substitution_map;
    }


#  ifndef DOXYGEN
    /* ---------------- Symbolic substitution map enlargement --------------*/


    void
    merge_substitution_maps(types::substitution_map &      symb_map_out,
                            const types::substitution_map &symb_map_in)
    {
      // Do this by hand so that we can perform some sanity checks
      for (const auto &entry : symb_map_in)
        {
          const typename types::substitution_map::const_iterator it_other =
            symb_map_out.find(entry.first);
          if (it_other == symb_map_out.end())
            symb_map_out.insert(std::make_pair(entry.first, entry.second));
          else
            {
              Assert(SE::eq(*(entry.second.get_RCP()),
                            *(it_other->second.get_RCP())),
                     ExcMessage("Key already in map, but values don't match"));
            }
        }
    }


    /* ---------------- Symbol substitution and evaluation --------------*/


    types::substitution_map
    resolve_explicit_dependencies(const types::substitution_map &symbol_values,
                                  const bool force_cyclic_dependency_resolution)
    {
      types::substitution_map symbol_values_resolved = symbol_values;
      const std::size_t       size                   = symbol_values.size();
      (void)size;
      for (auto &entry : symbol_values_resolved)
        {
          // Perform dictionary-based substitution to
          // resolve all explicit relations in a map.
          // Instead of checking by value (and thus having
          // to store a temporary value), we check to see
          // if the hash of the map entry changes.
          Expression & out = entry.second;
          SE::hash_t   hash_old;
          SE::hash_t   hash_new = out.get_RCP()->hash();
          unsigned int iter     = 0;
          do
            {
              // Write the substituted value straight back
              // into the map.
              if (force_cyclic_dependency_resolution)
                {
                  // Here we substitute in symbol_values_resolved instead of
                  // symbol_values, in order to break any cyclic dependencies.
                  // The earlier entries in the dictionary are in this way
                  // guaranteed to be resolved before any subsequent entries,
                  // thereby breaking the dependency loop.
                  out = out.substitute(symbol_values_resolved);
                }
              else
                {
                  out = out.substitute(symbol_values);
                }

              // Compute and store the hash of the new object
              hash_old = hash_new;
              hash_new = out.get_RCP()->hash();
              AssertThrow(
                iter < size,
                ExcMessage(
                  "Unresolvable cyclic dependency detected in substitution map."));
              ++iter;
            }
          while (hash_new != hash_old);
        }

      return symbol_values_resolved;
    }


    Expression
    substitute(const Expression &             expression,
               const types::substitution_map &substitution_map)
    {
      return expression.substitute(substitution_map);
    }
#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
