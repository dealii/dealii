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

#ifndef dealii_differentiation_sd_symengine_tensor_operations_h
#define dealii_differentiation_sd_symengine_tensor_operations_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
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
     * @param[in] symbol An identifier (or name) for the vector of returned
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
     * @param[in] symbol An identifier (or name) for the tensor of returned
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
     * @param[in] symbol An identifier (or name) for the tensor of returned
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
     * @param[in] symbol An identifier (or name) for the vector of returned
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
     * @param[in] symbol An identifier (or name) for the tensor of returned
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
     * @param[in] symbol An identifier (or name) for the tensor of returned
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

    /**
     * @name Symbolic differentiation
     */
    //@{

    /**
     * Return the symbolic result of computing the partial derivative of the
     * scalar @p f with respect to the tensor @p T.
     *
     * @param[in] f A scalar symbolic function or (dependent) expression.
     * @param[in] T A tensor of symbolic (independent) variables.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{T}}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Expression &f, const Tensor<rank, dim, Expression> &T);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * scalar @p f with respect to the symmetric tensor @p S.
     *
     * @param[in] f A scalar symbolic function or (dependent) expression.
     * @param[in] S A symmetric tensor of symbolic (independent) variables.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{S}}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Expression &                            f,
                  const SymmetricTensor<rank, dim, Expression> &S);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * rank-0 tensor (or scalar) @p f with respect to the tensor @p T.
     *
     * @param[in] f A rank-0 tensor symbolic function or (dependent) expression.
     * @param[in] T A tensor of symbolic (independent) variables.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{T}}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &   f,
                  const Tensor<rank, dim, Expression> &T);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * rank-0 tensor (or scalar) @p f with respect to the symmetric tensor @p S.
     *
     * @param[in] f A rank-0 tensor symbolic function or (dependent) expression.
     * @param[in] S A symmetric tensor of symbolic (independent) variables.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{S}}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &            f,
                  const SymmetricTensor<rank, dim, Expression> &S);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T with respect to the scalar @p x.
     *
     * @param[in] T A tensor of symbolic functions or (dependent) expressions.
     * @param[in] x A scalar symbolic (independent) variable.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{T}}{\partial x}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &T, const Expression &x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S with respect to the scalar @p x.
     *
     * @param[in] S A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] x A scalar symbolic (independent) variable.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{S}}{\partial x}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &S,
                  const Expression &                            x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T with respect to the rank-0 tensor @p x.
     *
     * @param[in] T A tensor of symbolic functions or (dependent) expressions.
     * @param[in] x A rank-0 tensor containing a symbolic (independent)
     *            variable.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{T}}{\partial x}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &T,
                  const Tensor<0, dim, Expression> &   x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S with respect to the rank-0 tensor @p x.
     *
     * @param[in] S A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] x A rank-0 tensor containing a symbolic (independent)
     *            variable.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{S}}{\partial x}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &S,
                  const Tensor<0, dim, Expression> &            x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T1 with respect to the tensor @p T2.
     *
     * @param[in] T1 A tensor of symbolic functions or (dependent) expressions.
     * @param[in] T2 A tensor of symbolic (independent) variables.
     * @return The tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{T}_{1}}{\partial
     * \mathbf{T}_{2}}$.
     */
    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &T1,
                  const Tensor<rank_2, dim, Expression> &T2);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S1 with respect to the symmetric tensor @p S2.
     *
     * @param[in] S1 A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] S2 A symmetric tensor of symbolic (independent)
     *            variables.
     * @return The symmetric tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{S}_{1}}{\partial
     * \mathbf{S}_{2}}$.
     */
    template <int rank_1, int rank_2, int dim>
    SymmetricTensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &S1,
                  const SymmetricTensor<rank_2, dim, Expression> &S2);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T with respect to the symmetric tensor @p S.
     *
     * @param[in] T A tensor of symbolic functions or (dependent) expressions.
     * @param[in] S A symmetric tensor of symbolic (independent)
     *            variables.
     * @return The tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{T}}{\partial \mathbf{S}}$.
     */
    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &         T,
                  const SymmetricTensor<rank_2, dim, Expression> &S);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S with respect to the tensor @p T.
     *
     * @param[in] S A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] T A tensor of symbolic (independent) variables.
     * @return The tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{S}}{\partial \mathbf{T}}$.
     */
    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &S,
                  const Tensor<rank_2, dim, Expression> &         T);

    //@}

    /**
     * @name Symbol map creation and manipulation
     */
    //@{

    /**
     * A convenience function for adding empty entries, with the key values
     * equal to the entries in the @p symbol_tensor, to the symbolic
     * map @p symbol_map.
     *
     * For more context which this function is used, see the other
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function.
     *
     * @tparam ignore_invalid_symbols See the other
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function for a detailed discussion on the role of this
     * template argument.
     *
     * @tparam SymbolicType A type that represents a symbolic variable.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p SymbolicType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              int rank,
              int dim,
              typename SymbolicType>
    void
    add_to_symbol_map(types::substitution_map &              symbol_map,
                      const Tensor<rank, dim, SymbolicType> &symbol_tensor);

    /**
     * A convenience function for adding empty entries, with the key values
     * equal to the entries in the @p symbol_tensor, to the symbolic
     * map @p symbol_map.
     *
     * For more context which this function is used, see the other
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function.
     *
     * @tparam ignore_invalid_symbols See the other
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function for a detailed discussion on the role of this
     * template argument.
     *
     * @tparam SymbolicType A type that represents a symbolic variable.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p SymbolicType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              int rank,
              int dim,
              typename SymbolicType>
    void
    add_to_symbol_map(
      types::substitution_map &                       symbol_map,
      const SymmetricTensor<rank, dim, SymbolicType> &symbol_tensor);

    /**
     * Find the input @p symbols in the @p substitution_map and set the entries
     * corresponding to the key values given by @p symbol_tensor to the values
     * given by @p value_tensor.
     *
     * This function may be used to safely transform an existing or null
     * symbolic map (one with uninitialized entries) into one that can be used
     * to conduct symbolic substitution operations (i.e., a substitution map).
     *
     * @tparam SymbolicType A type that represents a symbolic variable.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p SymbolicType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <int rank, int dim, typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &              substitution_map,
      const Tensor<rank, dim, SymbolicType> &symbol_tensor,
      const Tensor<rank, dim, ValueType> &   value_tensor);

    /**
     * Find the input @p symbols in the @p substitution_map and set the entries
     * corresponding to the key values given by @p symbol_tensor to the values
     * given by @p value_tensor.
     *
     * This function may be used to safely transform an existing or null
     * symbolic map (one with uninitialized entries) into one that can be used
     * to conduct symbolic substitution operations (i.e., a substitution map).
     *
     * @tparam SymbolicType A type that represents a symbolic variable.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p SymbolicType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <int rank, int dim, typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &                       substitution_map,
      const SymmetricTensor<rank, dim, SymbolicType> &symbol_tensor,
      const SymmetricTensor<rank, dim, ValueType> &   value_tensor);

    //@}

    /**
     * @name Symbol substitution map creation
     */
    //@{

    /**
     * Return a substitution map that has the entry keys given by the
     * @p symbol_tensor and the values given by the @p value_tensor. It is
     * expected that all key entries be valid symbols or symbolic expressions.
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &,const ValueType &)`
     * function.
     *
     * @tparam ExpressionType A type that represents a symbolic expression.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p ExpressionType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <int rank, int dim, typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const Tensor<rank, dim, ExpressionType> &symbol_tensor,
      const Tensor<rank, dim, ValueType> &     value_tensor);

    /**
     * Return a substitution map that has the entry keys given by the
     * @p symbol_tensor and the values given by the @p value_tensor. It is
     * expected that all key entries be valid symbols or symbolic expressions.
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &,const ValueType &)`
     * function.
     *
     * @tparam ExpressionType A type that represents a symbolic expression.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p ExpressionType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <int rank, int dim, typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const SymmetricTensor<rank, dim, ExpressionType> &symbol_tensor,
      const SymmetricTensor<rank, dim, ValueType> &     value_tensor);

    //@}

    /**
     * @name Symbol substitution map enlargement
     */
    //@{

    /**
     * A convenience function for adding an entry to the @p substitution_map.
     * The new entries will have the keys given in the @p symbol_tensor with
     * their paired values extracted from the corresponding elements of the
     * @p value_tensor.
     *
     * For more context which this function is used, see the other
     * `add_to_substitution_map(types::substitution_map &, const Expression &,
     * const Expression &)` function.
     *
     * @tparam ignore_invalid_symbols See the other
     * `add_to_substitution_map(types::substitution_map &, const Expression &,
     * const Expression &)` function for a detailed discussion on the role of
     * this template argument.
     *
     * @tparam ExpressionType A type that represents a symbolic expression.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p ExpressionType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <bool ignore_invalid_symbols = false,
              int  rank,
              int  dim,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                substitution_map,
      const Tensor<rank, dim, ExpressionType> &symbol_tensor,
      const Tensor<rank, dim, ValueType> &     value_tensor);

    /**
     * A convenience function for adding an entry to the @p substitution_map.
     * The new entries will have the keys given in the @p symbol_tensor with
     * their paired values extracted from the corresponding elements of the @p value_tensor.
     *
     * For more context which this function is used, see the other
     * `add_to_substitution_map(types::substitution_map &,const Expression &,
     * const Expression &)` function.
     *
     * @tparam ignore_invalid_symbols See the other
     * `add_to_substitution_map(types::substitution_map &, const Expression &,
     * const Expression &)` function for a detailed discussion on the role of
     * this template argument.
     *
     * @tparam ExpressionType A type that represents a symbolic expression.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, although if the @p ValueType is not supported
     *         by this class then a user-defined @p ExpressionType should be
     *         used.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <bool ignore_invalid_symbols = false,
              int  rank,
              int  dim,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                         substitution_map,
      const SymmetricTensor<rank, dim, ExpressionType> &symbol_tensor,
      const SymmetricTensor<rank, dim, ValueType> &     value_tensor);

    //@}

    /**
     * @name Symbol substitution and evaluation
     */
    //@{

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * tensor of symbolic expressions.
     * The symbols in the @p expression_tensor that correspond to the entry keys
     * of the @p substitution_map are substituted with the map entry's associated
     * value.
     * This substitution function may be used to give a set of symbolic
     * variables either a numeric interpretation or some symbolic definition.
     *
     * For more information regarding the performance of symbolic substitution,
     * and the outcome of evaluation using a substitution map with cyclic
     * dependencies, see the
     * `substitute(const Expression &, const types::substitution_map &)`
     * function.
     *
     * @note It is not required that all symbolic expressions be fully resolved
     * when using this function. In other words, partial substitutions are
     * valid.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    substitute(const Tensor<rank, dim, Expression> &expression_tensor,
               const types::substitution_map &      substitution_map);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * symmetric tensor of symbolic expressions.
     * The symbols in the @p expression_tensor that correspond to the entry keys
     * of the @p substitution_map are substituted with the map entry's associated
     * value.
     * This substitution function may be used to give a set of symbolic
     * variables either a numeric interpretation or some symbolic definition.
     *
     * For more information regarding the performance of symbolic substitution,
     * and the outcome of evaluation using a substitution map with cyclic
     * dependencies, see the
     * `substitute(const Expression &, const types::substitution_map &)`
     * function.
     *
     * @note It is not required that all symbolic expressions be fully resolved
     * when using this function. In other words, partial substitutions are
     * valid.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    substitute(const SymmetricTensor<rank, dim, Expression> &expression_tensor,
               const types::substitution_map &               substitution_map);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * tensor of symbolic expressions, and immediately evaluate the tensorial
     * result.
     * The symbols in the @p expression_tensor that correspond to the entry keys
     * of the @p substitution_map are substituted with the map entry's associated
     * value.
     * This substitution function is used to give a set of symbolic variables
     * a numeric interpretation with the returned result being of the type
     * specified by the @p ValueType template argument.
     *
     * For more information regarding the performance of symbolic substitution,
     * and the outcome of evaluation using a substitution map with cyclic
     * dependencies, see the
     * `substitute(const Expression &, const types::substitution_map &)`
     * function.
     *
     * @note It is required that all symbols in @p expression_tensor be
     * successfully resolved by the @p substitution_map.
     * If only partial substitution is performed, then an error is thrown.
     *
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. In the context of this particular
     *         function, this template parameter is typically arithmetic in
     *         nature.
     */
    template <typename ValueType, int rank, int dim>
    Tensor<rank, dim, ValueType>
    substitute_and_evaluate(
      const Tensor<rank, dim, Expression> &expression_tensor,
      const types::substitution_map &      substitution_map);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * symmetric tensor of symbolic expressions, and immediately evaluate the
     * tensorial result.
     * The symbols in the @p expression_tensor that correspond to the entry keys
     * of the @p substitution_map are substituted with the map entry's associated
     * value.
     * This substitution function is used to give a set of symbolic variables
     * a numeric interpretation, with the returned result being of the type
     * specified by the @p ValueType template argument.
     *
     * For more information regarding the performance of symbolic substitution,
     * and the outcome of evaluation using a substitution map with cyclic
     * dependencies, see the
     * `substitute(const Expression &, const types::substitution_map &)`
     * function.
     *
     * @note It is required that all symbols in @p expression_tensor be
     * successfully resolved by the @p substitution_map.
     * If only partial substitution is performed, then an error is thrown.
     *
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. In the context of this particular
     *         function, this template parameter is typically arithmetic in
     *         nature.
     */
    template <typename ValueType, int rank, int dim>
    SymmetricTensor<rank, dim, ValueType>
    substitute_and_evaluate(
      const SymmetricTensor<rank, dim, Expression> &expression_tensor,
      const types::substitution_map &               substitution_map);

    //@}

  } // namespace SD
} // namespace Differentiation


/* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN

namespace Differentiation
{
  namespace SD
  {
    /* ---------------- Symbolic differentiation --------------*/


    namespace internal
    {
      template <int dim>
      TableIndices<4>
      make_rank_4_tensor_indices(const unsigned int idx_i,
                                 const unsigned int idx_j)
      {
        const TableIndices<2> indices_i(
          SymmetricTensor<2, dim>::unrolled_to_component_indices(idx_i));
        const TableIndices<2> indices_j(
          SymmetricTensor<2, dim>::unrolled_to_component_indices(idx_j));
        return TableIndices<4>(indices_i[0],
                               indices_i[1],
                               indices_j[0],
                               indices_j[1]);
      }


      template <int rank_1, int rank_2>
      TableIndices<rank_1 + rank_2>
      concatenate_indices(const TableIndices<rank_1> &indices_1,
                          const TableIndices<rank_2> &indices_2)
      {
        TableIndices<rank_1 + rank_2> indices_out;
        for (unsigned int i = 0; i < rank_1; ++i)
          indices_out[i] = indices_1[i];
        for (unsigned int j = 0; j < rank_2; ++j)
          indices_out[rank_1 + j] = indices_2[j];
        return indices_out;
      }


      template <int rank>
      TableIndices<rank>
      transpose_indices(const TableIndices<rank> &indices)
      {
        return indices;
      }


      template <>
      inline TableIndices<2>
      transpose_indices(const TableIndices<2> &indices)
      {
        return TableIndices<2>(indices[1], indices[0]);
      }


      template <int rank, int dim, typename ValueType>
      bool
      is_symmetric_component(const TableIndices<rank> &,
                             const Tensor<rank, dim, ValueType> &)
      {
        return false;
      }


      template <int rank, int dim, typename ValueType>
      bool
      is_symmetric_component(const TableIndices<rank> &,
                             const SymmetricTensor<rank, dim, ValueType> &)
      {
        static_assert(
          rank == 0 || rank == 2,
          "Querying symmetric component for non rank-2 symmetric tensor index is not allowed.");
        return false;
      }


      template <int dim, typename ValueType>
      bool
      is_symmetric_component(const TableIndices<2> &table_indices,
                             const SymmetricTensor<2, dim, ValueType> &)
      {
        return table_indices[0] != table_indices[1];
      }


      template <int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<0, dim, ValueType>
      scalar_diff_tensor(const ValueType &                    func,
                         const TensorType<0, dim, ValueType> &op)
      {
        return differentiate(func, op);
      }


      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      scalar_diff_tensor(const ValueType &                       func,
                         const TensorType<rank, dim, ValueType> &op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(func, op[indices]);

            if (is_symmetric_component(indices, op))
              out[indices] *= 0.5;
          }
        return out;
      }


      // Specialization for rank-0 tensor
      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      tensor_diff_tensor(const TensorType<0, dim, ValueType> &   func,
                         const TensorType<rank, dim, ValueType> &op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(func, op[indices]);

            if (is_symmetric_component(indices, op))
              out[indices] *= 0.5;
          }
        return out;
      }


      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      tensor_diff_scalar(const TensorType<rank, dim, ValueType> &funcs,
                         const ValueType &                       op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(funcs[indices], op);
          }
        return out;
      }


      // Specialization for rank-0 tensor
      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      tensor_diff_tensor(const TensorType<rank, dim, ValueType> &funcs,
                         const TensorType<0, dim, ValueType> &   op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(funcs[indices], op);
          }
        return out;
      }


      // For either symmetric or normal tensors
      template <int rank_1,
                int rank_2,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank_1 + rank_2, dim, ValueType>
      tensor_diff_tensor(const TensorType<rank_1, dim, ValueType> &funcs,
                         const TensorType<rank_2, dim, ValueType> &op)
      {
        TensorType<rank_1 + rank_2, dim, ValueType> out;
        for (unsigned int i = 0; i < funcs.n_independent_components; ++i)
          {
            const TableIndices<rank_1> indices_i(
              funcs.unrolled_to_component_indices(i));
            for (unsigned int j = 0; j < op.n_independent_components; ++j)
              {
                const TableIndices<rank_2> indices_j(
                  op.unrolled_to_component_indices(j));
                const TableIndices<rank_1 + rank_2> indices_out =
                  concatenate_indices(indices_i, indices_j);

                out[indices_out] =
                  differentiate(funcs[indices_i], op[indices_j]);

                if (is_symmetric_component(indices_j, op))
                  out[indices_out] *= 0.5;
              }
          }
        return out;
      }


      // For mixed symmetric/standard tensors
      // The return type is always a standard tensor, since we cannot be sure
      // that any symmetries exist in either the function tensor or the
      // differential operator.
      template <int rank_1,
                int rank_2,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType_1,
                template <int, int, typename> class TensorType_2>
      Tensor<rank_1 + rank_2, dim, ValueType>
      tensor_diff_tensor(const TensorType_1<rank_1, dim, ValueType> &funcs,
                         const TensorType_2<rank_2, dim, ValueType> &op)
      {
        Tensor<rank_1 + rank_2, dim, ValueType> out;
        for (unsigned int i = 0; i < funcs.n_independent_components; ++i)
          {
            const TableIndices<rank_1> indices_i(
              funcs.unrolled_to_component_indices(i));
            for (unsigned int j = 0; j < op.n_independent_components; ++j)
              {
                const TableIndices<rank_2> indices_j(
                  op.unrolled_to_component_indices(j));
                const TableIndices<rank_1 + rank_2> indices_out =
                  concatenate_indices(indices_i, indices_j);

                out[indices_out] =
                  differentiate(funcs[indices_i], op[indices_j]);

                if (is_symmetric_component(indices_j, op))
                  out[indices_out] *= 0.5;

                // TODO: Implement for SymmetricTensor<4,dim,...>
                if (std::is_same<TensorType_1<rank_1, dim, ValueType>,
                                 SymmetricTensor<2, dim, ValueType>>::
                      value) // Symmetric function
                  {
                    const TableIndices<rank_1 + rank_2> indices_out_t =
                      concatenate_indices(transpose_indices(indices_i),
                                          indices_j);
                    out[indices_out_t] = out[indices_out];
                  }
                else if (std::is_same<TensorType_2<rank_2, dim, ValueType>,
                                      SymmetricTensor<2, dim, ValueType>>::
                           value) // Symmetric operator
                  {
                    const TableIndices<rank_1 + rank_2> indices_out_t =
                      concatenate_indices(indices_i,
                                          transpose_indices(indices_j));
                    out[indices_out_t] = out[indices_out];
                  }
                else
                  {
                    Assert(
                      false,
                      ExcMessage(
                        "Expect mixed tensor differentiation to have at least "
                        "one component stemming from a symmetric tensor."));
                  }
              }
          }
        return out;
      }

    } // namespace internal


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Expression &                   func,
                  const Tensor<rank, dim, Expression> &op)
    {
      return internal::scalar_diff_tensor(func, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Expression &                            func,
                  const SymmetricTensor<rank, dim, Expression> &op)
    {
      return internal::scalar_diff_tensor(func, op);
    }


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &   func,
                  const Tensor<rank, dim, Expression> &op)
    {
      return internal::scalar_diff_tensor(func, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &            func,
                  const SymmetricTensor<rank, dim, Expression> &op)
    {
      // Ensure that this returns a symmetric tensor by
      // invoking the scalar value function
      const Expression tmp = func;
      return internal::scalar_diff_tensor(tmp, op);
    }


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &symbol_tensor,
                  const Expression &                   op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &symbol_tensor,
                  const Tensor<0, dim, Expression> &   op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &symbol_tensor,
                  const Expression &                            op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &symbol_tensor,
                  const Tensor<0, dim, Expression> &            op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &symbol_tensor,
                  const Tensor<rank_2, dim, Expression> &op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    SymmetricTensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &symbol_tensor,
                  const SymmetricTensor<rank_2, dim, Expression> &op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &         symbol_tensor,
                  const SymmetricTensor<rank_2, dim, Expression> &op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &symbol_tensor,
                  const Tensor<rank_2, dim, Expression> &         op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    /* ---------------- Symbol map creation and manipulation --------------*/


    namespace internal
    {
      template <typename SymbolicType,
                typename ValueType,
                int rank,
                int dim,
                template <int, int, typename> class TensorType>
      void
      set_tensor_value_in_symbol_map(
        types::substitution_map &                  substitution_map,
        const TensorType<rank, dim, SymbolicType> &symbol_tensor,
        const TensorType<rank, dim, ValueType> &   value_tensor)
      {
        TensorType<rank, dim, Expression> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            set_value_in_symbol_map(substitution_map,
                                    symbol_tensor[indices],
                                    value_tensor[indices]);
          }
      }


      template <typename SymbolicType, typename ValueType, int dim>
      void
      set_tensor_value_in_symbol_map(
        types::substitution_map &                    substitution_map,
        const SymmetricTensor<4, dim, SymbolicType> &symbol_tensor,
        const SymmetricTensor<4, dim, ValueType> &   value_tensor)
      {
        SymmetricTensor<4, dim, Expression> out;
        for (unsigned int i = 0;
             i < SymmetricTensor<2, dim>::n_independent_components;
             ++i)
          for (unsigned int j = 0;
               j < SymmetricTensor<2, dim>::n_independent_components;
               ++j)
            {
              const TableIndices<4> indices =
                make_rank_4_tensor_indices<dim>(i, j);
              set_value_in_symbol_map(substitution_map,
                                      symbol_tensor[indices],
                                      value_tensor[indices]);
            }
      }
    } // namespace internal


    template <bool ignore_invalid_symbols,
              typename ValueType,
              int rank,
              int dim,
              typename SymbolicType>
    void
    add_to_symbol_map(types::substitution_map &              symbol_map,
                      const Tensor<rank, dim, SymbolicType> &symbol_tensor)
    {
      // Call the above function
      add_to_substitution_map<ignore_invalid_symbols>(
        symbol_map, symbol_tensor, Tensor<rank, dim, ValueType>());
    }


    template <bool ignore_invalid_symbols,
              typename ValueType,
              int rank,
              int dim,
              typename SymbolicType>
    void
    add_to_symbol_map(
      types::substitution_map &                       symbol_map,
      const SymmetricTensor<rank, dim, SymbolicType> &symbol_tensor)
    {
      // Call the above function
      add_to_substitution_map<ignore_invalid_symbols>(
        symbol_map, symbol_tensor, SymmetricTensor<rank, dim, ValueType>());
    }


    template <int rank, int dim, typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &              substitution_map,
      const Tensor<rank, dim, SymbolicType> &symbol_tensor,
      const Tensor<rank, dim, ValueType> &   value_tensor)
    {
      internal::set_tensor_value_in_symbol_map(substitution_map,
                                               symbol_tensor,
                                               value_tensor);
    }


    template <int rank, int dim, typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &                       substitution_map,
      const SymmetricTensor<rank, dim, SymbolicType> &symbol_tensor,
      const SymmetricTensor<rank, dim, ValueType> &   value_tensor)
    {
      internal::set_tensor_value_in_symbol_map(substitution_map,
                                               symbol_tensor,
                                               value_tensor);
    }


    /* ------------------ Symbol substitution map creation ----------------*/


    template <int rank, int dim, typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const Tensor<rank, dim, ExpressionType> &symbol_tensor,
      const Tensor<rank, dim, ValueType> &     value_tensor)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbol_tensor, value_tensor);
      return substitution_map;
    }


    template <int rank, int dim, typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const SymmetricTensor<rank, dim, ExpressionType> &symbol_tensor,
      const SymmetricTensor<rank, dim, ValueType> &     value_tensor)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbol_tensor, value_tensor);
      return substitution_map;
    }


    /* ---------------- Symbolic substitution map enlargement --------------*/


    namespace internal
    {
      template <int rank,
                int dim,
                typename ExpressionType,
                typename ValueType,
                template <int, int, typename> class TensorType>
      std::vector<std::pair<ExpressionType, ValueType>>
      make_tensor_entries_for_substitution_map(
        const TensorType<rank, dim, ExpressionType> &symbol_tensor,
        const TensorType<rank, dim, ValueType> &     value_tensor)
      {
        std::vector<std::pair<ExpressionType, ValueType>> symbol_values;
        for (unsigned int i = 0; i < symbol_tensor.n_independent_components;
             ++i)
          {
            const TableIndices<rank> indices(
              symbol_tensor.unrolled_to_component_indices(i));
            symbol_values.push_back(
              std::make_pair(symbol_tensor[indices], value_tensor[indices]));
          }
        return symbol_values;
      }


      template <int dim, typename ExpressionType, typename ValueType>
      std::vector<std::pair<ExpressionType, ValueType>>
      make_tensor_entries_for_substitution_map(
        const Tensor<0, dim, ExpressionType> &symbol_tensor,
        const Tensor<0, dim, ValueType> &     value_tensor)
      {
        const ExpressionType &expression = symbol_tensor;
        const ValueType &     value      = value_tensor;
        return {std::make_pair(expression, value)};
      }


      template <int dim, typename ExpressionType, typename ValueType>
      std::vector<std::pair<ExpressionType, ValueType>>
      make_tensor_entries_for_substitution_map(
        const SymmetricTensor<4, dim, ExpressionType> &symbol_tensor,
        const SymmetricTensor<4, dim, ValueType> &     value_tensor)
      {
        std::vector<std::pair<ExpressionType, ValueType>> symbol_values;
        for (unsigned int i = 0;
             i < SymmetricTensor<2, dim>::n_independent_components;
             ++i)
          for (unsigned int j = 0;
               j < SymmetricTensor<2, dim>::n_independent_components;
               ++j)
            {
              const TableIndices<4> indices =
                make_rank_4_tensor_indices<dim>(i, j);
              symbol_values.push_back(
                std::make_pair(symbol_tensor[indices], value_tensor[indices]));
            }
        return symbol_values;
      }
    } // namespace internal


    template <bool ignore_invalid_symbols,
              int  rank,
              int  dim,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                substitution_map,
      const Tensor<rank, dim, ExpressionType> &symbol_tensor,
      const Tensor<rank, dim, ValueType> &     value_tensor)
    {
      add_to_substitution_map<ignore_invalid_symbols>(
        substitution_map,
        internal::make_tensor_entries_for_substitution_map(symbol_tensor,
                                                           value_tensor));
    }


    template <bool ignore_invalid_symbols,
              int  rank,
              int  dim,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                         substitution_map,
      const SymmetricTensor<rank, dim, ExpressionType> &symbol_tensor,
      const SymmetricTensor<rank, dim, ValueType> &     value_tensor)
    {
      add_to_substitution_map<ignore_invalid_symbols>(
        substitution_map,
        internal::make_tensor_entries_for_substitution_map(symbol_tensor,
                                                           value_tensor));
    }


    /* ---------------- Symbol substitution and evaluation --------------*/


    namespace internal
    {
      template <int rank,
                int dim,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, Expression>
      substitute_tensor(
        const TensorType<rank, dim, Expression> &expression_tensor,
        const types::substitution_map &          substitution_map)
      {
        TensorType<rank, dim, Expression> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] =
              substitute(expression_tensor[indices], substitution_map);
          }
        return out;
      }


      template <int dim>
      Tensor<0, dim, Expression>
      substitute_tensor(const Tensor<0, dim, Expression> &expression_tensor,
                        const types::substitution_map &   substitution_map)
      {
        const Expression &expression = expression_tensor;
        return substitute(expression, substitution_map);
      }


      template <int dim>
      SymmetricTensor<4, dim, Expression>
      substitute_tensor(
        const SymmetricTensor<4, dim, Expression> &expression_tensor,
        const types::substitution_map &            substitution_map)
      {
        SymmetricTensor<4, dim, Expression> out;
        for (unsigned int i = 0;
             i < SymmetricTensor<2, dim>::n_independent_components;
             ++i)
          for (unsigned int j = 0;
               j < SymmetricTensor<2, dim>::n_independent_components;
               ++j)
            {
              const TableIndices<4> indices =
                make_rank_4_tensor_indices<dim>(i, j);
              out[indices] =
                substitute(expression_tensor[indices], substitution_map);
            }
        return out;
      }


      template <typename ValueType,
                int rank,
                int dim,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      substitute_and_evaluate_tensor(
        const TensorType<rank, dim, Expression> &expression_tensor,
        const types::substitution_map &          substitution_map)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] =
              substitute_and_evaluate<ValueType>(expression_tensor[indices],
                                                 substitution_map);
          }
        return out;
      }


      template <typename ValueType, int dim>
      Tensor<0, dim, ValueType>
      substitute_and_evaluate_tensor(
        const Tensor<0, dim, Expression> &expression_tensor,
        const types::substitution_map &   substitution_map)
      {
        const Expression &expression = expression_tensor;
        return substitute_and_evaluate<ValueType>(expression, substitution_map);
      }


      template <typename ValueType, int dim>
      SymmetricTensor<4, dim, ValueType>
      substitute_and_evaluate_tensor(
        const SymmetricTensor<4, dim, Expression> &expression_tensor,
        const types::substitution_map &            substitution_map)
      {
        SymmetricTensor<4, dim, ValueType> out;
        for (unsigned int i = 0;
             i < SymmetricTensor<2, dim>::n_independent_components;
             ++i)
          for (unsigned int j = 0;
               j < SymmetricTensor<2, dim>::n_independent_components;
               ++j)
            {
              const TableIndices<4> indices =
                make_rank_4_tensor_indices<dim>(i, j);
              out[indices] =
                substitute_and_evaluate<ValueType>(expression_tensor[indices],
                                                   substitution_map);
            }
        return out;
      }
    } // namespace internal


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    substitute(const Tensor<rank, dim, Expression> &expression_tensor,
               const types::substitution_map &      substitution_map)
    {
      return internal::substitute_tensor(expression_tensor, substitution_map);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    substitute(const SymmetricTensor<rank, dim, Expression> &expression_tensor,
               const types::substitution_map &               substitution_map)
    {
      return internal::substitute_tensor(expression_tensor, substitution_map);
    }


    template <typename ValueType, int rank, int dim>
    Tensor<rank, dim, ValueType>
    substitute_and_evaluate(
      const Tensor<rank, dim, Expression> &expression_tensor,
      const types::substitution_map &      substitution_map)
    {
      return internal::substitute_and_evaluate_tensor<ValueType>(
        expression_tensor, substitution_map);
    }


    template <typename ValueType, int rank, int dim>
    SymmetricTensor<rank, dim, ValueType>
    substitute_and_evaluate(
      const SymmetricTensor<rank, dim, Expression> &expression_tensor,
      const types::substitution_map &               substitution_map)
    {
      return internal::substitute_and_evaluate_tensor<ValueType>(
        expression_tensor, substitution_map);
    }



  } // namespace SD
} // namespace Differentiation

#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
