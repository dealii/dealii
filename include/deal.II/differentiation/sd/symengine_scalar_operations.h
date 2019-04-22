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

#ifndef dealii_differentiation_sd_symengine_scalar_operations_h
#define dealii_differentiation_sd_symengine_scalar_operations_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/numbers.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_types.h>

#  include <algorithm>
#  include <type_traits>
#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    /**
     * @name Symbolic variable creation
     */
    //@{

    /**
     * Return an Expression representing a scalar symbolic variable
     * with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"x"` then the scalar
     * symbolic variable that is returned represents the scalar $x$.
     *
     * @param[in] symbol An indentifier (or name) for the returned symbolic
     *            variable.
     * @return A scalar symbolic variable with the name @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    Expression
    make_symbol(const std::string &symbol);

    /**
     * Return an Expression representing a scalar symbolic function
     * with the identifier specified by @p symbol. The function's symbolic
     * dependencies are specified by the input @p arguments.
     *
     * For example, if the @p symbol is the string `"f"`, and the
     * arguments to the function that is generated  are the
     * symbolic variable `x` and the symbolic expression `y+z`, then the
     * generic symbolic function that is returned represents $f(x, y+z)$.
     *
     * @param[in] symbol An indentifier (or name) for the returned symbolic
     *            function.
     * @param[in] arguments A vector of input arguments to the returned
     *            symbolic function.
     * @return A generic symbolic function with the identifier @p symbolic and
     *         the number of input arguments equal to the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    Expression
    make_symbolic_function(const std::string &         symbol,
                           const types::symbol_vector &arguments);

    /**
     * Return an Expression representing a scalar symbolic function
     * with the identifier specified by @p symbol. The function's symbolic
     * dependencies are specified by the keys to the input @p arguments map;
     * the values stored in the map are ignored.
     *
     * For example, if the @p symbol is the string `"f"`, and the
     * arguments to the function that is generated are the
     * symbolic variable `x` and the symbolic expression `y+z`, then the
     * generic symbolic function that is returned represents $f(x, y+z)$.
     *
     * @param[in] symbol An indentifier (or name) for the returned symbolic
     *            function.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic function.
     * @return A generic symbolic function with the identifier @p symbolic and
     *         the number of input arguments equal to the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    Expression
    make_symbolic_function(const std::string &            symbol,
                           const types::substitution_map &arguments);

    //@}

    /**
     * @name Symbolic differentiation
     */
    //@{

    /**
     * Return the symbolic result of computing the partial derivative of the
     * scalar @p f with respect to the scalar @p x.
     * In most use cases the function or variable @f would be called the
     * dependent variable, and @p x the independent variable.
     *
     * @param[in] f A scalar symbolic function or (dependent) variable.
     * @param[in] x A scalar symbolic (independent) variable.
     * @return The symbolic function or variable representing the result
     *         $\frac{\partial f}{\partial x}$.
     */
    Expression
    differentiate(const Expression &f, const Expression &x);

    //@}

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
