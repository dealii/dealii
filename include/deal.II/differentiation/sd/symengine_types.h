// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_sd_symengine_types_h
#define dealii_differentiation_sd_symengine_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <boost/serialization/map.hpp>

#  include <map>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    // Forward declarations
    class Expression;


    namespace types
    {
      namespace internal
      {
        /**
         * A comparator for Expressions used as keys in maps.
         */
        struct ExpressionKeyLess
        {
          bool
          operator()(const SD::Expression &lhs,
                     const SD::Expression &rhs) const;
        };
      } // namespace internal

      /**
       * Type definition for a value substitution map.
       *
       * This serves the same purpose as a `SymEngine::map_basic_basic`, which
       * is equivalent to a `std::map<SymEngine::RCP<const SymEngine::Basic>,
       * SymEngine::RCP<const SymEngine::Basic>>`.
       */
      using substitution_map =
        std::map<SD::Expression, SD::Expression, internal::ExpressionKeyLess>;

      /**
       * Type definition for a vector of symbols.
       *
       * This serves the same purpose as a `SymEngine::vec_basic`, which is
       * equivalent to a `std::vector<SymEngine::RCP<const SymEngine::Basic>>`.
       */
      using symbol_vector = std::vector<SD::Expression>;

    } // namespace types

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE


#  ifndef DOXYGEN

// Add serialization capability for SD::types::internal::ExpressionKeyLess
// We need to define this so that we can use this comparator in maps that
// are to be serialized.
namespace boost
{
  namespace serialization
  {
    namespace SD = dealii::Differentiation::SD;

    template <typename Archive>
    void
    serialize(Archive & /*ar*/,
              SD::types::internal::ExpressionKeyLess & /*cmp*/,
              unsigned int /*version*/)
    {
      // Nothing to do.
    }
  } // namespace serialization
} // namespace boost

#  endif // DOXYGEN

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_types_h
