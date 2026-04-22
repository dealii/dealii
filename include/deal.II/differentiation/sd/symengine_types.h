// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_differentiation_sd_symengine_types_h
#define dealii_differentiation_sd_symengine_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <boost/serialization/map.hpp>

#  include <map>
#  include <vector>

#endif // DEAL_II_WITH_SYMENGINE

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SYMENGINE
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

#endif // DEAL_II_WITH_SYMENGINE

// Close the dealii namespace, but do not close the export{} block yet if we
// are building C++20 modules:
DEAL_II_NAMESPACE_CLOSE // Do not convert for module purposes

#ifdef DEAL_II_WITH_SYMENGINE

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

#endif // DEAL_II_WITH_SYMENGINE

// Re-open and then close the namespace, and this time also close the
// export{} block if we are building C++20 modules:
DEAL_II_NAMESPACE_OPEN // Do not convert for module purposes
  DEAL_II_NAMESPACE_CLOSE

#endif // dealii_differentiation_sd_symengine_types_h
