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
#  include <deal.II/differentiation/sd/symengine_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;

    namespace Utilities
    {
#  ifndef DOXYGEN

      SE::map_basic_basic
      convert_expression_map_to_basic_map(
        const SD::types::substitution_map &substitution_map)
      {
        SE::map_basic_basic sub_map;
        for (const auto &entry : substitution_map)
          sub_map[entry.first.get_RCP()] = entry.second.get_RCP();
        return sub_map;
      }



      SD::types::substitution_map
      convert_basic_map_to_expression_map(
        const SymEngine::map_basic_basic &substitution_map)
      {
        SD::types::substitution_map sub_map;
        for (const auto &entry : substitution_map)
          sub_map[Expression(entry.first)] = SD::Expression(entry.second);
        return sub_map;
      }



      SE::vec_basic
      convert_expression_vector_to_basic_vector(
        const SD::types::symbol_vector &symbol_vector)
      {
        SE::vec_basic symb_vec;
        symb_vec.reserve(symbol_vector.size());
        for (const auto &entry : symbol_vector)
          symb_vec.push_back(entry.get_RCP());
        return symb_vec;
      }



      SD::types::symbol_vector
      convert_basic_vector_to_expression_vector(
        const SE::vec_basic &symbol_vector)
      {
        SD::types::symbol_vector symb_vec;
        symb_vec.reserve(symbol_vector.size());
        for (const auto &entry : symbol_vector)
          symb_vec.push_back(SD::Expression(entry));
        return symb_vec;
      }



      SD::types::symbol_vector
      extract_symbols(const SD::types::substitution_map &substitution_values)
      {
        SD::types::symbol_vector symbols;
        symbols.reserve(substitution_values.size());

        for (const auto &substitution : substitution_values)
          symbols.push_back(substitution.first);

        return symbols;
      }



      std::vector<std::pair<SD::Expression, SD::Expression>>
      convert_basic_pair_vector_to_expression_pair_vector(
        const SymEngine::vec_pair &symbol_value_vector)
      {
        std::vector<std::pair<SD::Expression, SD::Expression>> symb_val_vec;
        symb_val_vec.reserve(symbol_value_vector.size());
        for (const auto &entry : symbol_value_vector)
          symb_val_vec.push_back(std::make_pair(SD::Expression(entry.first),
                                                SD::Expression(entry.second)));
        return symb_val_vec;
      }

#  endif // DOXYGEN

    } // namespace Utilities

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
