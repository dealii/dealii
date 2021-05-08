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

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_number_traits_h
