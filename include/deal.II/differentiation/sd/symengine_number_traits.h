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

#ifndef dealii_differentiation_sd_symengine_number_traits_h
#define dealii_differentiation_sd_symengine_number_traits_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <symengine/expression.h>

#  include <type_traits>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * symbolically differentiable number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     */
    template <typename NumberType>
    struct is_sd_number : std::false_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * SymEngine number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     */
    template <typename NumberType>
    struct is_symengine_number : std::false_type
    {};


    /*--- SymEngine Expression class ---*/


    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * symbolically differentiable number or not.
     * This is a specialization for the SymEngine Expression class.
     */
    template <>
    struct is_symengine_number<SymEngine::Expression> : std::true_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * symbolically differentiable number or not.
     * This is a specialization for the SymEngine Expression class.
     */
    template <>
    struct is_sd_number<SymEngine::Expression> : std::true_type
    {};

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_number_traits_h
