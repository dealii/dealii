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

#ifndef dealii_differentiation_sd_symengine_tensor_operations_h
#define dealii_differentiation_sd_symengine_tensor_operations_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_types.h>

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
     * Return a vector of Expressions representing a vectorial symbolic
     * variable with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"v"` then the vectorial
     * symbolic variable that is returned represents the vector $v$. Each
     * component of $v$ is prefixed by the given @p symbol, and has a suffix that
     * indicates its component index.
     *
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the vector of returned
     *            symbolic variables.
     * @return A vector (a rank-1 tensor) of symbolic variables with the name of
     *         each individual component prefixed by @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int dim>
    Tensor<1, dim, Expression>
    make_vector_of_symbols(const std::string &symbol);

    /**
     * Return a tensor of Expressions representing a tensorial symbolic
     * variable with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"T"` then the tensorial
     * symbolic variable that is returned represents the vector $T$. Each
     * component of $T$ is prefixed by the given @p symbol, and has a suffix that
     * indicates its component indices.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic variables.
     * @return A tensor of symbolic variables with the name of each individual
     *         component prefixed by @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    make_tensor_of_symbols(const std::string &symbol);

    /**
     * Return a symmetric tensor of Expressions representing a tensorial
     * symbolic variable with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"S"` then the tensorial
     * symbolic variable that is returned represents the vector $S$. Each
     * component of $S$ is prefixed by the given @p symbol, and has a suffix that
     * indicates its component indices.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic variables.
     * @return A tensor of symbolic variables with the name of each individual
     *         component prefixed by @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    make_symmetric_tensor_of_symbols(const std::string &symbol);

    /**
     * Return a vector of Expression representing a vectorial symbolic
     * function with the identifier specified by @p symbol. The functions'
     * symbolic dependencies are specified by the keys to the input
     * @p arguments map; the values stored in the map are ignored.
     *
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the vector of returned
     *            symbolic functions.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic functions.
     * @return A vector (a rank-1 tensor) of generic symbolic functions with the
     *         name of each individual component prefixed by @p symbol, a suffix
     *         that indicates its component index, and the number of input
     *         arguments equal to the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int dim>
    Tensor<1, dim, Expression>
    make_vector_of_symbolic_functions(const std::string &            symbol,
                                      const types::substitution_map &arguments);

    /**
     * Return a tensor of Expression representing a tensorial symbolic
     * function with the identifier specified by @p symbol. The functions'
     * symbolic dependencies are specified by the keys to the input
     * @p arguments map; the values stored in the map are ignored.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic functions.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic functions.
     * @return A tensor of generic symbolic functions with the name of each
     *         individual component prefixed by @p symbol, a suffix that
     *         indicates its component indeices, and the number of input
     *         arguments equal to  the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    make_tensor_of_symbolic_functions(const std::string &            symbol,
                                      const types::substitution_map &arguments);

    /**
     * Return a symmetric tensor of Expression representing a tensorial
     * symbolic function with the identifier specified by @p symbol. The
     * functions' symbolic dependencies are specified by the keys to the input
     * @p arguments map; the values stored in the map are ignored.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic functions.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic functions.
     * @return A symmetric tensor of generic symbolic functions with the name of
     *         each individual component prefixed by @p symbol, a suffix that
     *         indicates its component indeices, and the number of input
     *         arguments equal to the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    make_symmetric_tensor_of_symbolic_functions(
      const std::string &            symbol,
      const types::substitution_map &arguments);

    //@}

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
