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

#ifndef dealii_differentiation_sd_symengine_utilities_h
#define dealii_differentiation_sd_symengine_utilities_h

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
    namespace Utilities
    {
      /**
       * Convert a map of Expressions to its SymEngine counterpart.
       */
      SymEngine::map_basic_basic
      convert_expression_map_to_basic_map(
        const SD::types::substitution_map &substitution_map);

      /**
       * Convert to a map of Expressions from its SymEngine counterpart.
       */
      SD::types::substitution_map
      convert_basic_map_to_expression_map(
        const SymEngine::map_basic_basic &substitution_map);

      /**
       * Convert a vector of Expressions to its SymEngine counterpart.
       */
      SymEngine::vec_basic
      convert_expression_vector_to_basic_vector(
        const SD::types::symbol_vector &symbol_vector);

      /**
       * Convert to a vector of Expressions from its SymEngine counterpart.
       */
      SD::types::symbol_vector
      convert_basic_vector_to_expression_vector(
        const SymEngine::vec_basic &symbol_vector);

      /**
       * Convert to a vector of pairs of Expressions from its SymEngine
       * counterpart.
       */
      std::vector<std::pair<Expression, Expression>>
      convert_basic_pair_vector_to_expression_pair_vector(
        const SymEngine::vec_pair &symbol_value_vector);

      /**
       * Extract the symbols (key entries) from a substitution map.
       *
       * @note It is guaranteed that the order of extraction of data into the
       * output vector is the same as that for extract_values().
       * That is to say that the unzipped key and value pairs as given by
       * extract_symbols() and extract_values() always have a 1:1
       * correspondence.
       */
      SD::types::symbol_vector
      extract_symbols(const SD::types::substitution_map &substitution_values);

      /**
       * Extract the values from a substitution map.
       * The value entries will be converted into the @p NumberType given
       * as a template parameter to this function via the @p ExpressionType.
       *
       * @note It is guaranteed that the order of extraction of data into the
       * output vector is the same as that for extract_symbols().
       * That is to say that the unzipped key and value pairs as given by
       * extract_symbols() and extract_values() always have a 1:1
       * correspondence.
       */
      template <typename NumberType, typename ExpressionType = SD::Expression>
      std::vector<NumberType>
      extract_values(const SD::types::substitution_map &substitution_values);

      /**
       * Print the key and value pairs stored in a substitution map.
       */
      template <typename StreamType>
      StreamType &
      print_substitution_map(
        StreamType &                       stream,
        const SD::types::substitution_map &symbol_value_map);

    } // namespace Utilities

  } // namespace SD
} // namespace Differentiation


/* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN


namespace Differentiation
{
  namespace SD
  {
    namespace Utilities
    {
      template <typename NumberType, typename ExpressionType>
      std::vector<NumberType>
      extract_values(const SD::types::substitution_map &substitution_values)
      {
        std::vector<NumberType> values;
        values.reserve(substitution_values.size());

        for (const auto &entry : substitution_values)
          values.push_back(
            static_cast<NumberType>(ExpressionType(entry.second)));

        return values;
      }


      template <typename StreamType>
      StreamType &
      print_substitution_map(
        StreamType &                       stream,
        const SD::types::substitution_map &symbol_value_map)
      {
        for (const auto &entry : symbol_value_map)
          stream << entry.first << " = " << entry.second << "\n";

        stream << std::flush;
        return stream;
      }

    } // namespace Utilities

  } // namespace SD
} // namespace Differentiation


#  endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_utilities_h
