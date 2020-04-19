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

#  include <symengine/basic.h>
#  include <symengine/symengine_rcp.h>

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
     * @param[in] symbol An identifier (or name) for the returned symbolic
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
     * @param[in] symbol An identifier (or name) for the returned symbolic
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
     * @param[in] symbol An identifier (or name) for the returned symbolic
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
     * In most use cases the function or variable @p f would be called the
     * dependent variable, and @p x the independent variable.
     *
     * @param[in] f A scalar symbolic function or (dependent) expression.
     * @param[in] x A scalar symbolic (independent) variable.
     * @return The symbolic function or expression representing the result
     *         $\frac{\partial f}{\partial x}$.
     */
    Expression
    differentiate(const Expression &f, const Expression &x);

    //@}

    /**
     * @name Symbol map creation and manipulation
     */
    //@{

    namespace internal
    {
      /**
       * Return whether or not an @p entry is a valid symbol that
       * we can expect to perform substitution on.
       */
      bool
      is_valid_substitution_symbol(const SymEngine::Basic &entry);

      /**
       * A convenience function to set the @p value associated with
       * the @p symbol in the @p substitution_map.
       *
       * Using this function ensures that the @p symbol is one that is
       * valid specifically for the purpose of symbolic substitution.
       * It must therefore represent a symbol or symbolic derivative,
       * otherwise an error will be thrown.
       */
      void
      set_value_in_symbol_map(
        types::substitution_map &                     substitution_map,
        const SymEngine::RCP<const SymEngine::Basic> &symbol,
        const SymEngine::RCP<const SymEngine::Basic> &value);
    } // namespace internal

    /**
     * Return a symbolic map that has a single entry with the key given by
     * the @p symbol.
     * It is expected that all entries to be added to the symbolic map are
     * valid symbols or symbolic expressions.
     *
     * @tparam ignore_invalid_symbols See the
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function for a detailed discussion on the role of this
     * template argument.
     *
     * @tparam SymbolicType Any symbolic type that is understood by the
     *         add_to_symbol_map() functions. This includes individual
     *         Expression, std::vector<Expression>, as well as
     *         Tensors and SymmetricTensors of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              typename SymbolicType>
    types::substitution_map
    make_symbol_map(const SymbolicType &symbol);

    /**
     * Return a symbolic map that has the entry keys given by @p symbol and all
     * @p other_symbols.
     * It is expected that all entries to be added to the symbolic map are
     * valid symbols or symbolic expressions.
     *
     * With this function it is possible to construct a symbolic map from
     * different types. An example may be as follows:
     *
     * @code
     *   const types::substitution_map symbol_map
     *     = make_symbol_map(
     *         Expression(...),
     *         Tensor<1,dim,Expression>(...),
     *         SymmetricTensor<2,dim,Expression>(...));
     * @endcode
     *
     * @tparam ignore_invalid_symbols See the
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function for a detailed discussion on the role of this
     * template argument.
     *
     * @tparam SymbolicType Any symbolic type that is understood by the
     *         add_to_symbol_map() functions. This includes individual
     *         Expression, std::vector<Expression>, as well as
     *         Tensors and SymmetricTensors of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     * @tparam Args A type associated with the parameter pack that contains
     *         any number of other @p SymbolicTypes. All types held by the
     *         parameter pack share the same restriction as the @p SymbolicType
     *         documented above.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              typename SymbolicType,
              typename... Args>
    types::substitution_map
    make_symbol_map(const SymbolicType &symbol, const Args &... other_symbols);

    /**
     * A convenience function for adding an empty entry, with the key value
     * given by @p symbol, to the symbolic map @p symbol_map.
     *
     * This function is guaranteed to create an ordering that is identical
     * to the typical add_to_substitution_map() call that is used when
     * constructing a map to perform symbol substitution.
     * It exists primarily to create an initial map that can be used in the
     * optimize() call to a BatchOptimizer, specifically if the values that
     * are to be substituted into the map are not known at the time that the
     * symbols used to construct symbolic expressions are defined.
     * This helps one conform to the requirement that the arguments sent into
     * lambda and LLVM JIT compiled functions (created by optimizing symbolic
     * expressions) (i) be the same, and (ii) have a constant ordering.
     *
     * @tparam ignore_invalid_symbols A template parameter that enforces whether
     *         or not the @p symbol has to be a valid one or not. In the
     *         overwhelming majority of cases, the default value of
     *         <tt>false</tt> should be selected, with the result that an
     *         exception will be thrown if the input @p symbolic is, in fact,
     *         not a symbolic value or expression.
     *         An exceptional case is, for example, when performing symbolic
     *         assembly on a finite element level. When extracting the symbolic
     *         equivalent of the shape function gradients using FEExtractors,
     *         the returned tensor will have some <em>a priori</em> determined
     *         zero-valued components. These trivial components are not valid
     *         symbols (as they are not symbolic expressions), and we would
     *         typically wish to guard against their (erroneous) inclusion. In
     *         this scenario, for convenience, one could set
     *         @p ignore_invalid_symbols to <tt>true</tt> and these zero-valued
     *         entries would be skipped over and ignored.
     *
     * @note In this function, the @p ValueType is somewhat arbitrary as
     * it is only used to create default-constructed values as entries in
     * the map.
     */
    template <bool ignore_invalid_symbols = false, typename ValueType = double>
    void
    add_to_symbol_map(types::substitution_map &symbol_map,
                      const Expression &       symbol);

    /**
     * A convenience function for adding an empty entry, with the key value
     * given by @p symbol, to the symbolic map @p symbol_map.
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
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *         The required condition is fulfilled when the @p SymbolicType
     *         can be explicitly converted to a
     *         `const SymEngine::RCP<const SymEngine::Basic> &`.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              typename SymbolicType,
              typename T = typename std::enable_if<
                !std::is_base_of<Expression, SymbolicType>::value &&
                dealii::internal::is_explicitly_convertible<
                  SymbolicType,
                  const SymEngine::RCP<const SymEngine::Basic> &>::value>::type>
    void
    add_to_symbol_map(types::substitution_map &symbol_map,
                      const SymbolicType &     symbol);

    /**
     * A convenience function for adding empty entries, with the key values
     * equal to the entries in @p symbols, to the symbolic map @p symbol_map.
     * It is expected that all entries in the input vector @p symbols be
     * of @p SymbolicType, compatible with the other add_to_symbol_map()
     * functions.
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
     * @tparam SymbolicType Any symbolic type that is understood by the
     *         add_to_symbol_map() functions. This includes an individual
     *         Expression, as well as Tensors and SymmetricTensors of
     *         Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              typename SymbolicType>
    void
    add_to_symbol_map(types::substitution_map &        symbol_map,
                      const std::vector<SymbolicType> &symbols);

    /**
     * A convenience function for adding empty entries, with the key values
     * equal to the key entries in @p other_symbols, to the symbolic
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
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     */
    template <bool ignore_invalid_symbols = false, typename ValueType = double>
    void
    add_to_symbol_map(types::substitution_map &      symbol_map,
                      const types::substitution_map &other_symbols);

    /**
     * A convenience function for adding empty entries, with the key values
     * equal to the entries in @p symbol plus @p other_symbols, to the symbolic
     * map @p symbol_map.
     * It is expected that all entries in @p symbol and @p other_symbols be
     * of a @p SymbolicType, compatible with the other add_to_symbol_map()
     * functions.
     *
     * For more context which this function is used, see the other
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function.
     *
     * With this function it is possible to add entries from different types
     * to a symbolic map. An example may be as follows:
     *
     * @code
     *   types::substitution_map symbol_map = ...;
     *   add_to_symbol_map(
     *     symbol_map,
     *     Expression(...),
     *     Tensor<1,dim,Expression>(...),
     *     SymmetricTensor<2,dim,Expression>(...));
     * @endcode
     *
     * @tparam ignore_invalid_symbols See the other
     * `add_to_symbol_map(types::substitution_map &, const Expression &)`
     * function for a detailed discussion on the role of this
     * template argument.
     *
     * @tparam SymbolicType Any symbolic type that is understood by the
     *         add_to_symbol_map() functions. This includes individual
     *         Expression, std::vector<Expression>, as well as
     *         Tensors and SymmetricTensors of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. This @p ValueType is somewhat
     *         arbitrary as it is only used to create default-constructed
     *         values as entries in the map.
     * @tparam Args A type associated with the parameter pack that contains
     *         any number of other @p SymbolicTypes. All types held by the
     *         parameter pack share the same restriction as the @p SymbolicType
     *         documented above.
     */
    template <bool ignore_invalid_symbols = false,
              typename ValueType          = double,
              typename SymbolicType,
              typename... Args>
    void
    add_to_symbol_map(types::substitution_map &symbol_map,
                      const SymbolicType &     symbol,
                      const Args &... other_symbols);

    /**
     * Find the entry for @p symbol in the @p substitution_map and set its
     * corresponding @p value.
     *
     * This function may be used to safely transform an existing or null
     * symbolic map (one with uninitialized entries) into one that can be used
     * to conduct symbolic substitution operations (i.e., a substitution map).
     */
    void
    set_value_in_symbol_map(types::substitution_map &substitution_map,
                            const Expression &       symbol,
                            const Expression &       value);

    /**
     * Find the entry for @p symbol in the @p substitution_map and set its
     * corresponding @p value.
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
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *         The required condition is fulfilled when the @p SymbolicType
     *         can be explicitly converted to a
     *         `const SymEngine::RCP<const SymEngine::Basic> &`, and it is
     *         possible to construct an @p SymbolicType directly from the
     *         @p ValueType.
     */
    template <typename SymbolicType,
              typename ValueType,
              typename T = typename std::enable_if<
                dealii::internal::is_explicitly_convertible<
                  SymbolicType,
                  const SymEngine::RCP<const SymEngine::Basic> &>::value &&
                std::is_constructible<SymbolicType, ValueType>::value>::type>
    void
    set_value_in_symbol_map(types::substitution_map &substitution_map,
                            const SymbolicType &     symbol,
                            const ValueType &        value);

    /**
     * Find the entries for @p symbols in the @p substitution_map and set their
     * corresponding @p values.
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
     *         @p SymbolicType can be constructed from.
     */
    template <typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(types::substitution_map &        substitution_map,
                            const std::vector<SymbolicType> &symbols,
                            const std::vector<ValueType> &   values);

    /**
     * Find the entry for @p symbols in the @p substitution_map and set their
     * corresponding @p values. The modified symbol will have the key given by
     * the first element of @p symbol_value and the value given by its second
     * element.
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
     *         @p SymbolicType can be constructed from.
     */
    template <typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value);

    /**
     * Find the entries for @p symbols in the @p substitution_map and set their
     * corresponding @p values, followed by the same operation for the
     * @p other_symbol_values.
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
     *         @p SymbolicType can be constructed from.
     * @tparam Args A type associated with the parameter pack that contains
     *         any number of other pairs of @p SymbolicTypes and @p ValueTypes.
     *         All types held by the parameter pack share the same restriction
     *         as the @p SymbolicType and @p ValueType documented above.
     */
    template <typename SymbolicType, typename ValueType, typename... Args>
    void
    set_value_in_symbol_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value,
      const Args &... other_symbol_values);

    /**
     * Find the entries for @p symbols in the @p substitution_map and set their
     * corresponding @p values. The modified symbol will have the key given by
     * the first element of each paired entry in the @p symbol_values vector
     * and the value given by its respective second element.
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
     *         @p SymbolicType can be constructed from.
     */
    template <typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &                              substitution_map,
      const std::vector<std::pair<SymbolicType, ValueType>> &symbol_values);

    /**
     * Find the entries for @p symbols in the @p substitution_map and set their
     * corresponding @p values. The modified symbol will have the key given by
     * the each element the @p symbol_values map and the value given by its
     * respective mapped element.
     *
     * This function may be used to safely transform an existing or null
     * symbolic map (one with uninitialized entries) into one that can be used
     * to conduct symbolic substitution operations (i.e., a substitution map).
     */
    void
    set_value_in_symbol_map(types::substitution_map &      substitution_map,
                            const types::substitution_map &symbol_values);

    //@}

    /**
     * @name Symbol substitution map creation
     */
    //@{

    /**
     * Return a substitution map that has the entry key given by @p symbol
     * and the value given by @p value. It is expected that the key entry
     * be valid symbol or symbolic expression.
     *
     * The values that map to a @p symbol would typically be of an arithmetic
     * type. However, in some instances is may be useful to map a symbolic type
     * to another symbolic type (i.e. perform partial substitution). In such
     * a situation the resolve_explicit_dependencies() function may be useful
     * to simplify the final substitution map by resolving all explicit
     * interdependencies between entries in the substitution map.
     */
    types::substitution_map
    make_substitution_map(const Expression &symbol, const Expression &value);

    /**
     * Return a substitution map that has the entry key given by @p symbol
     * and the value given by @p value. It is expected that the key entry
     * be valid symbol or symbolic expression.
     *
     * The values that map to a @p symbol would typically be of a @p ValueType
     * (i.e., an arithmetic type).
     * However, in some instances is may be useful to map a symbolic type
     * to another symbolic type (i.e. perform partial substitution). In such
     * a situation the resolve_explicit_dependencies() function may be useful
     * to simplify the final substitution map by resolving all explicit
     * interdependencies between entries in the substitution map.
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
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *         The required condition is fulfilled when the @p ExpressionType
     *         can be explicitly converted to a
     *         `const SymEngine::RCP<const SymEngine::Basic> &`, and it is
     *         possible to construct an @p ExpressionType directly from the
     *         @p ValueType.
     */
    template <typename ExpressionType,
              typename ValueType,
              typename T = typename std::enable_if<
                dealii::internal::is_explicitly_convertible<
                  ExpressionType,
                  const SymEngine::RCP<const SymEngine::Basic> &>::value &&
                std::is_constructible<ExpressionType, ValueType>::value>::type>
    types::substitution_map
    make_substitution_map(const ExpressionType &symbol, const ValueType &value);

    /**
     * Return a substitution map that has the entry keys given by @p symbols
     * and the values given by @p values. It is expected that all key entries
     * be valid symbols or symbolic expressions.
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &, const ValueType &)`
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
    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(const std::vector<ExpressionType> &symbols,
                          const std::vector<ValueType> &     values);

    /**
     * Return a substitution map that has the key given by the first entry in
     * @p symbol_value, and the value of its second entry. It is expected that
     * the key entry be a valid symbol or symbolic expression.
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &, const ValueType &)`
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
    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const std::pair<ExpressionType, ValueType> &symbol_value);

    /**
     * Return a substitution map that has the keys given by the first entry of
     * each element of @p symbol_values, and the values given its second entry.
     * It is expected that all key entries be valid symbols or symbolic
     * expressions.
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &, const ValueType &)`
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
    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const std::vector<std::pair<ExpressionType, ValueType>> &symbol_values);

    /**
     * Return a substitution map that has the key given by the first entry in
     * @p symbol_value, and the value of its second entry, followed by the
     * addition of the @p other_symbol_values. It is expected that all key
     * entries be valid symbols or symbolic expressions.
     *
     * With this function it is possible to construct a symbolic substitution
     * map from different types, so long as there exists a
     * add_to_substitution_map() function with the signature corresponding to
     * the pair types. An example may be as follows:
     *
     * @code
     *   const types::substitution_map substitution_map
     *     = make_substitution_map(
     *         std::make_pair(Expression(...), 3),
     *         std::make_pair(Tensor<1,dim,Expression>(...),
     *                        Tensor<1,dim,float>(...)),
     *         std::make_pair(SymmetricTensor<2,dim,Expression>(...),
     *                        SymmetricTensor<2,dim,double>(...)));
     * @endcode
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &, const ValueType &)`
     * function.
     *
     * @tparam ExpressionType Any symbolic expression type that is understood
     *         by the make_substitution_map() functions. This includes
     *         individual Expression, as well as Tensors and SymmetricTensors
     *         of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     * @tparam Args A type associated with the parameter pack that contains
     *         any number of other @p ExpressionTypes. All types held by the
     *         parameter pack share the same restriction as the
     *         @p ExpressionType documented above.
     */
    template <typename ExpressionType, typename ValueType, typename... Args>
    types::substitution_map
    make_substitution_map(
      const std::pair<ExpressionType, ValueType> &symbol_value,
      const Args &... other_symbol_values);

    //@}

    /**
     * @name Symbol substitution map enlargement
     */
    //@{

    namespace internal
    {
      /**
       * A convenience function to add an entry to the @p substitution_map.
       * The new entry will have the key given by @p symbol with its paired
       * @p value. Such maps are required to perform substitution of symbolic
       * expressions, where all entries given as keys in the map are substituted
       * by their value counterparts.
       *
       * There are cases where it is convenient to simply ignore the fact that
       * we may be trying to add invalid symbols. For example, one may wish to
       * add a tensor where only some entries are symbols and the others are
       * zero'd. In this case we ensure that the user knows what they're doing
       * by forcing them to override this safety mechanism with a template
       * parameter. These types of functions are typically called in a manner
       * indicating that the user knows exactly what they're passing into it.
       *
       * @tparam ignore_invalid_symbols A template parameter that enforces whether
       *         or not the @p symbol has to be a valid one or not. In the
       *         overwhelming majority of cases, the default value of
       *         <tt>false</tt> should be selected, with the result that an
       *         exception will be thrown if the input @p symbolic is, in fact,
       *         not a symbolic value or expression.
       */
      template <bool ignore_invalid_symbols = false>
      void
      add_to_substitution_map(
        types::substitution_map &                     substitution_map,
        const SymEngine::RCP<const SymEngine::Basic> &symbol,
        const SymEngine::RCP<const SymEngine::Basic> &value);
    } // namespace internal

    /**
     * A convenience function to add an entry to the @p substitution_map.
     * The new entry will have the key given by @p symbol with its paired
     * @p value. Such maps are required to perform substitution of symbolic
     * expressions, with key entries being exchanges with their pair values.
     *
     * This function is guaranteed to create an ordering that is identical
     * to the typical add_to_symbol_map() call that may be used when initially
     * configuring a BatchOptimizer.
     * This helps one conform to the requirement that the arguments sent into
     * lambda and LLVM JIT compiled functions (created by optimizing symbolic
     * expressions) (i) be the same, and (ii) have a constant ordering.
     *
     * @tparam ignore_invalid_symbols A template parameter that enforces whether
     *         or not the @p symbol has to be a valid one or not. In the
     *         overwhelming majority of cases, the default value of
     *         <tt>false</tt> should be selected, with the result that an
     *         exception will be thrown if the input @p symbolic is, in fact,
     *         not a symbolic value or expression.
     *         An exceptional case is, for example, when performing symbolic
     *         assembly on a finite element level. When extracting the symbolic
     *         equivalent of the shape function gradients using FEExtractors,
     *         the returned tensor will have some <em>a priori</em> determined
     *         zero-valued components. These trivial components are not valid
     *         symbols (as they are not symbolic expressions), and we would
     *         typically wish to guard against their (erroneous) inclusion. In
     *         this scenario, for convenience, one could set
     *         @p ignore_invalid_symbols to <tt>true</tt> and these zero-valued
     *         entries would be skipped over and ignored.
     */
    template <bool ignore_invalid_symbols = false>
    void
    add_to_substitution_map(types::substitution_map &substitution_map,
                            const Expression &       symbol,
                            const Expression &       value);

    /**
     * A convenience function to add an entry to the @p substitution_map.
     * The new entry will have the key given by @p symbol with its paired
     * @p value.
     *
     * The @p ExpressionType will be used to convert the @p value to a
     * compatible SymEngine number type.
     * It is therefore required that the @p ExpressionType
     * 1. can be constructed from a @p ValueType, and that
     * 2. it is convertible to a `const SymEngine::RCP<const SymEngine::Basic>
     * &`.
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
              typename ExpressionType,
              typename ValueType,
              typename = typename std::enable_if<
                dealii::internal::is_explicitly_convertible<
                  ExpressionType,
                  const SymEngine::RCP<const SymEngine::Basic> &>::value &&
                std::is_constructible<ExpressionType, ValueType>::value>::type>
    void
    add_to_substitution_map(types::substitution_map &substitution_map,
                            const ExpressionType &   symbol,
                            const ValueType &        value);

    /**
     * A convenience function for adding multiple entries to the
     * @p substitution_map. The new entries will have the keys given in the
     * @p symbols vector, each of which will be paired index-wise with its
     * corresponding element in the @p values vector.
     *
     * The class represented by the @p ExpressionType template parameter will
     * be used to convert the p value to a compatible SymEngine number type.
     * It is therefore required that the @p ExpressionType
     * 1. can be constructed from a @p ValueType, and that
     * 2. it is convertible to a `const SymEngine::RCP<const SymEngine::Basic>
     * &`.
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
              typename ExpressionType,
              typename ValueType,
              typename = typename std::enable_if<
                dealii::internal::is_explicitly_convertible<
                  ExpressionType,
                  const SymEngine::RCP<const SymEngine::Basic> &>::value &&
                std::is_constructible<ExpressionType, ValueType>::value>::type>
    void
    add_to_substitution_map(types::substitution_map &          substitution_map,
                            const std::vector<ExpressionType> &symbols,
                            const std::vector<ValueType> &     values);

    /**
     * A convenience function for adding multiple entries to the
     * @p substitution_map. The new entries will have the keys given in the
     * @p symbols vector, each of which will be paired index-wise with its
     * corresponding element in the @p values vector. It is expected that there
     * are no duplicate entries between the two maps.
     *
     * For more context which this function is used, see the other
     * `add_to_substitution_map(types::substitution_map &, const Expression &,
     * const Expression &)` function.
     *
     * @tparam ignore_invalid_symbols See the other
     * `add_to_substitution_map(types::substitution_map &, const Expression &,
     * const Expression &)` function for a detailed discussion on the role of
     * this template argument.
     */
    template <bool ignore_invalid_symbols = false>
    void
    add_to_substitution_map(types::substitution_map &      substitution_map,
                            const types::substitution_map &symbol_values);

    /**
     * A convenience function to add an entry to the @p substitution_map.
     * The new entry will have the key given by the first element of
     * @p symbol_value and the value given by its second element.
     * It is expected that the key entry be a valid symbol or symbolic
     * expression, and that the paired @p symbol_value elements are compatible
     * with the other add_to_substitution_map() functions.
     *
     * The @p ExpressionType and its associated @p ValueType need not be scalar
     * types. So, for example, this function could be used to add tensor-valued
     * data to the map in the following way:
     *
     * @code
     *   types::substitution_map substitution_map = ...;
     *   add_to_substitution_map(
     *     substitution_map,
     *     std::make_pair(Tensor<1,dim,Expression>(...),
     *                    Tensor<1,dim,double>(...))
     *   );
     * @endcode
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
     * @tparam ExpressionType Any symbolic expression type that is understood
     *         by the add_to_substitution_map() functions. This includes
     *         individual Expression, as well as Tensors and SymmetricTensors
     *         of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <bool ignore_invalid_symbols = false,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                   substitution_map,
      const std::pair<ExpressionType, ValueType> &symbol_value);

    /**
     * A convenience function for adding multiple entries to the
     * @p substitution_map. The new entries will have the keys given by first
     * entry of each element of @p symbol_values, and the values given its
     * second entry. It is expected that the key entry be a valid symbols or
     * symbolic expressions, and that the paired @p symbol_value elements are
     * compatible with the other add_to_substitution_map() functions.
     *
     * The @p ExpressionType and its associated @p ValueType need not be scalar
     * types. So, for example, this function could be used to add tensor-valued
     * data to the map in the following way:
     *
     * @code
     *   types::substitution_map substitution_map = ...;
     *   using vector_entry_t = std::vector<std::pair<
     *     Tensor<1,dim,Expression>, Tensor<1,dim,double>
     *   >>;
     *   add_to_substitution_map(
     *     substitution_map,
     *     vector_entry_t{
     *       {Tensor<1,dim,Expression>(...), Tensor<1,dim,double>(...)},
     *       {Tensor<1,dim,Expression>(...), Tensor<1,dim,double>(...)}
     *     });
     * @endcode
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
     * @tparam ExpressionType Any symbolic expression type that is understood
     *         by the add_to_substitution_map() functions. This includes
     *         individual Expression, as well as Tensors and SymmetricTensors
     *         of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     */
    template <bool ignore_invalid_symbols = false,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                                substitution_map,
      const std::vector<std::pair<ExpressionType, ValueType>> &symbol_values);

    /**
     * A convenience function for adding multiple entries to the
     * @p substitution_map. The new entries will have the keys given by first
     * entry of each element of @p symbol_values, and the values given its
     * second entry, along with the addition of the @p other_symbol_values.
     * It is expected that the key entry be a valid symbols or symbolic
     * expressions, and that the paired @p symbol_value elements are compatible
     * with the other add_to_substitution_map() functions.
     *
     * With this function it is possible to construct a symbolic substitution
     * map from different types, so long as there exists a
     * add_to_substitution_map() function with the signature corresponding to
     * the pair types. An example may be as follows:
     *
     * @code
     *   types::substitution_map substitution_map = ...;
     *   add_to_substitution_map(
     *     substitution_map,
     *     std::make_pair(Expression(...), 3),
     *     std::make_pair(Tensor<1,dim,Expression>(...),
     *                    Tensor<1,dim,float>(...)),
     *     std::make_pair(SymmetricTensor<2,dim,Expression>(...),
     *                    SymmetricTensor<2,dim,double>(...)));
     * @endcode
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * `make_substitution_map(const Expression &, const ValueType &)`
     * function.
     *
     * @tparam ExpressionType Any symbolic expression type that is understood
     *         by the add_to_substitution_map() functions. This includes
     *         individual Expression, as well as Tensors and SymmetricTensors
     *         of Expressions.
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type or be a special type that a user-defined
     *         @p ExpressionType can be constructed from.
     * @tparam Args A type associated with the parameter pack that contains
     *         any number of other @p ExpressionTypes. All types held by the
     *         parameter pack share the same restriction as the
     *         @p ExpressionType documented above.
     */
    template <bool ignore_invalid_symbols = false,
              typename ExpressionType,
              typename ValueType,
              typename... Args>
    void
    add_to_substitution_map(
      types::substitution_map &                   substitution_map,
      const std::pair<ExpressionType, ValueType> &symbol_value,
      const Args &... other_symbol_values);

    /**
     * Concatenate two symbolic maps, merging a second map @p substitution_map_in
     * in-place into the initial and resultant map @p substitution_map_out. The map
     * @p substitution_map_out need not initially be empty.
     *
     * @note Duplicate symbols (keys) in the maps are permitted, so long as their
     * values are equal. If this is not the case then an error will be thrown.
     */
    void
    merge_substitution_maps(types::substitution_map &      substitution_map_out,
                            const types::substitution_map &substitution_map_in);

    /**
     * Concatenate multiple symbolic maps, merging the maps @p substitution_map_in
     * and @p other_substitution_maps_in, a collection of other maps,
     * in-place into the resultant map @p substitution_map_out. The map
     * @p substitution_map_out need not initially be empty.
     *
     * @note Duplicate symbols (keys) in the maps are permitted, so long as their
     * values are equal. If this is not the case then an error will be thrown.
     */
    template <typename... Args>
    void
    merge_substitution_maps(types::substitution_map &      substitution_map_out,
                            const types::substitution_map &substitution_map_in,
                            const Args &... other_substitution_maps_in);

    /**
     * Concatenate multiple symbolic maps, merging the maps @p substitution_map_in
     * and @p other_substitution_maps_in and returning the result.
     *
     * @note Duplicate symbols (keys) in the maps are permitted, so long as their
     * values are equal. If this is not the case then an error will be thrown.
     */
    template <typename... Args>
    types::substitution_map
    merge_substitution_maps(const types::substitution_map &substitution_map_in,
                            const Args &... other_substitution_maps_in);

    //@}

    /**
     * @name Symbol substitution and evaluation
     */
    //@{

    /**
     * Return a substitution map that has any explicit interdependencies between
     * the entries of the input @p substitution_map resolved.
     *
     * The @p force_cyclic_dependency_resolution flag exists to ensure, if
     * desired, that no cyclic dependencies can exist in the returned map.
     * If a cyclic dependency exists in the input substitution map,
     * @p substitution_map, then with this flag set to <tt>true</tt> the
     * dependency cycle is broken by a dictionary-ordered substitution.
     * For example, if the substitution map contains two entries
     * <tt>map["a"] -> "b"</tt> and <tt>map["b"] -> "a"</tt>, then the result
     * of calling this function would be a map with the elements
     * <tt>map["a"] -> "a"</tt> and <tt>map["b"] -> "a"</tt>.
     *
     * If one symbol is an explicit function of another, and it is desired that
     * all their values are completely resolved, then it may be necessary to
     * perform substitution a number of times before the result is finalized.
     * This function performs substitution sweeps for a set of symbolic
     * variables until all explicit relationships between the symbols in the map
     * have been resolved. Whether each entry returns a symbolic or real value
     * depends on the nature of the values stored in the substitution map. If
     * the values associated with a key are also symbolic then the returned
     * result may still be symbolic in nature. The terminal result of using the
     * input substitution map, @p symbol_values, is then guaranteed to be
     * rendered by a single substitition of the returned dependency-resolved
     * map.
     *
     * Example:
     * If <tt>map["a"] -> 1</tt> and
     * <tt>map["b"] -> "a"+ 2</tt>, then then the function $f(a,b(a)) = a+b$
     * will be evaluated and the result $f\vert_{a=1,b=a+2} = 3+a$ is determined
     * upon the completion of the first sweep. A second sweep is therefore
     * necessary to resolve the final symbol, and the returned value is
     * ultimately $f = [3+a]_{a=1} = 4$. By resolving the explicit relationships
     * between all symbols in the map, we determine that <tt>map["a"] -> 1</tt>
     * and <tt>map["b"] -> 1 + 2 = 3</tt> and thus, using only one substitution,
     * that $f = a+b = 1 + 3 = 4$.
     */
    types::substitution_map
    resolve_explicit_dependencies(
      const types::substitution_map &substitution_map,
      const bool force_cyclic_dependency_resolution = false);

    /**
     * Return a substitution map that has any explicit interdependencies between
     * the entries of the map generated by the paired elements in the
     * @p symbol_values vector resolved.
     * The @p force_cyclic_dependency_resolution exists to ensure, if desired,
     * that no cyclic dependencies can exist in the returned map.
     *
     * This function performs substitution sweeps for a set of symbolic
     * variables until all explicit relationships between the symbols in the map
     * have been resolved. Whether each entry returns a symbolic or real value
     * depends on the chosen ValueType and the values represented therein. If
     * the ValueType is another Expression, and they contain symbols then
     * the returned result may still be symbolic in nature.
     *
     * For an example of what this function does, see the documentation for the
     * other
     * `resolve_explicit_dependencies(const types::substitution_map &)`
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
    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    resolve_explicit_dependencies(
      const std::vector<std::pair<ExpressionType, ValueType>> &symbol_values,
      const bool force_cyclic_dependency_resolution = false);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * symbolic expression.
     * The symbols in the @p expression that correspond to the entry keys
     * of the @p substitution_map are substituted with the map entry's associated
     * value.
     * This substitution function may be used to give a set of symbolic
     * variables either a numeric interpretation or some symbolic definition.
     *
     * @note It is not required that all symbolic expressions be fully resolved
     * when using this function. In other words, partial substitutions are
     * valid.
     *
     * @note This function call is typically expensive, as by default it performs a
     * dictionary substitution for the symbols in the symbolic expression.
     * Should the numerical value of some symbolic expression be desired, then
     * this performance deficit may be mitigated through the use of the
     * BatchOptimizer class.
     * Situation dependent, the overhead of using a typical dictionary based
     * substitution may be on par with that of a substitution performed
     * using a BatchOptimizer. This is because there is an overhead to setting
     * up the optimizer, so this should be taken into consideration if
     * substitution is to occur for the given symbolic
     * expression only a few times.
     *
     * @note If the symbols stored in the map are explicitly dependent on one another,
     * then the returned result depends on the order in which the map is
     * traversed. It is recommended to first resolve all interdependencies in
     * the map using the resolve_explicit_dependencies() function.
     *
     * Examples:
     * <ol>
     *   <li>If <tt>map["a"] == 1</tt> and <tt>map["b"] == "a" + 2</tt>, then
     *   the function $f(a,b(a)) := a+b$ will be evaluated and the result
     *   $f\vert_{a=1,b=a+2} = 3+a$ is returned. This return is because the
     *   symbol "a" is substituted throughout the function first, and only
     *   then is the symbol "b(a)" substituted, by which time its explicit
     *   dependency on "a" cannot be resolved.</li>
     *
     *   <li>If <tt>map["a"] == "b"+2</tt> and <tt>map["b"] == 1</tt>, then
     *   the function $f(a(b),b): = a+b$ will be evaluated and the result
     *   $f\vert_{a=b+2, b} = [b+2+b]_{b=1} = 4$ is returned. This is because
     *   the explicitly dependent symbol "a(b)" is substituted first followed
     *   by the symbol "b".</li>
     * </ol>
     */
    Expression
    substitute(const Expression &             expression,
               const types::substitution_map &substitution_map);

    /**
     * Perform a substitution of the @p symbol into the given
     * @p expression. All matches are assigned
     * the corresponding @p value.
     * This substitution function may be used to give a set of symbolic
     * variables either a numeric interpretation or some symbolic definition.
     *
     * For more information regarding the performance of symbolic substitution,
     * see the other
     * `substitute(const Expression &, const types::substitution_map &)`
     * function.
     *
     * @note It is not required that all symbolic expressions be fully resolved
     * when using this function. In other words, partial substitutions are
     * valid.
     *
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. Although it is typically
     *         arithmetic in nature, it may also represent another symbolic
     *         expression type.
     */
    template <typename ValueType>
    Expression
    substitute(const Expression &expression,
               const Expression &symbol,
               const ValueType & value);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * symbolic expression.
     * The symbols in the @p expression that correspond to a matching entry key
     * of the @p symbol_values vector entry are substituted by the entry's associated
     * value.
     * This substitution function may be used to give a set of symbolic
     * variables either a numeric interpretation or some symbolic definition.
     *
     * For more information regarding the performance of symbolic substitution,
     * see the other
     * `substitute(const Expression &, const types::substitution_map &)`
     * function.
     *
     * @note It is not required that all symbolic expressions be fully resolved
     * when using this function. In other words, partial substitutions are
     * valid.
     *
     * @tparam ExpressionType A type that represents a symbolic expression.
     *         The Differentiation::SD::Expression class is often suitable for
     *         this purpose, but this may also represent Tensors and
     *         SymmetricTensors of Expressions.
     * @tparam Args Any symbolic type and value combination that is understood
     *         by the make_substitution_map() functions. This includes
     *         arguments involving individual Expressions,
     *         std::vector<Expression>, as well as Tensors and SymmetricTensors
     *         of Expressions.
     */
    template <typename ExpressionType, typename... Args>
    ExpressionType
    substitute(const ExpressionType &expression, const Args &... symbol_values);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * symbolic function, and immediately evaluate the result.
     * The symbols in the @p expression that correspond to the entry keys
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
     * @note It is required that all symbols in the @p expression be
     * successfully resolved by the @p substitution_map.
     * If only partial substitution is performed, then an error is thrown.
     *
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. In the context of this particular
     *         function, this template parameter is typically arithmetic in
     *         nature.
     */
    template <typename ValueType>
    ValueType
    substitute_and_evaluate(const Expression &             expression,
                            const types::substitution_map &substitution_map);

    /**
     * Perform a single substitution sweep of a set of symbols into the given
     * symbolic function, and immediately evaluate the result.
     * The symbols in the @p expression that correspond to the entry keys
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
     * @note It is required that all symbols in the @p expression be
     * successfully resolved by the substitution map that is generated with the
     * input collection of @p symbol_values.
     * If only partial substitution is performed, then an error is thrown.
     *
     * @tparam ValueType A type that corresponds to the @p value that the
     *         @p symbol is to represent. In the context of this particular
     *         function, this template parameter is typically arithmetic in
     *         nature.
     * @tparam Args Any symbolic type and value combination that is understood
     *         by the make_substitution_map() functions. This includes
     *         arguments involving individual Expressions,
     *         std::vector<Expression>, as well as Tensors and SymmetricTensors
     *         of Expressions.
     */
    template <typename ValueType, typename... Args>
    ValueType
    substitute_and_evaluate(const Expression &expression,
                            const Args &... symbol_values);

    //@}

  } // namespace SD
} // namespace Differentiation


/* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN


namespace Differentiation
{
  namespace SD
  {
    /* ---------------- Symbol map creation and manipulation --------------*/


    template <bool ignore_invalid_symbols,
              typename ValueType,
              typename SymbolicType>
    types::substitution_map
    make_symbol_map(const SymbolicType &symbol)
    {
      types::substitution_map symbol_map;
      add_to_symbol_map<ignore_invalid_symbols, ValueType>(symbol_map, symbol);
      return symbol_map;
    }


    template <bool ignore_invalid_symbols,
              typename ValueType,
              typename SymbolicType,
              typename... Args>
    types::substitution_map
    make_symbol_map(const SymbolicType &symbol, const Args &... other_symbols)
    {
      types::substitution_map symbol_map;
      add_to_symbol_map<ignore_invalid_symbols, ValueType>(symbol_map,
                                                           symbol,
                                                           other_symbols...);
      return symbol_map;
    }


    template <bool ignore_invalid_symbols, typename ValueType>
    void
    add_to_symbol_map(types::substitution_map &symbol_map,
                      const Expression &       symbol)
    {
      // Call the above function
      add_to_substitution_map<ignore_invalid_symbols>(
        symbol_map,
        symbol,
        dealii::internal::NumberType<ValueType>::value(0.0));
    }


    template <bool ignore_invalid_symbols,
              typename ValueType,
              typename SymbolicType,
              typename>
    void
    add_to_symbol_map(types::substitution_map &symbol_map,
                      const SymbolicType &     symbol)
    {
      // Call the above function
      using SE_RCP_Basic = const SymEngine::RCP<const SymEngine::Basic> &;
      add_to_substitution_map<ignore_invalid_symbols>(
        symbol_map,
        static_cast<SE_RCP_Basic>(symbol),
        dealii::internal::NumberType<ValueType>::value(0.0));
    }


    template <bool ignore_invalid_symbols,
              typename ValueType,
              typename SymbolicType>
    void
    add_to_symbol_map(types::substitution_map &        symbol_map,
                      const std::vector<SymbolicType> &symbols)
    {
      for (const auto &symbol : symbols)
        add_to_symbol_map<ignore_invalid_symbols, ValueType>(symbol_map,
                                                             symbol);
    }


    template <bool ignore_invalid_symbols, typename ValueType>
    void
    add_to_symbol_map(types::substitution_map &      symbol_map,
                      const types::substitution_map &other_symbols)
    {
      // We should be cautious as to whether or not the user has
      // hand-made the other_symbols map to be "merged" in.
      // So instead of blindly merging using the merge_substitution_maps()
      // function, we do it by hand and check for invalid symbols
      // if required.
      for (types::substitution_map::const_iterator it = other_symbols.begin();
           it != other_symbols.end();
           ++it)
        {
          Assert(symbol_map.find(it->first) == symbol_map.end(),
                 ExcMessage("Entry already exists in symbol map"));
          add_to_symbol_map<ignore_invalid_symbols, ValueType>(
            symbol_map, Expression(it->first));
        }
    }


    template <bool ignore_invalid_symbols,
              typename ValueType,
              typename SymbolicType,
              typename... Args>
    void
    add_to_symbol_map(types::substitution_map &symbol_map,
                      const SymbolicType &     symbol,
                      const Args &... other_symbols)
    {
      add_to_symbol_map<ignore_invalid_symbols, ValueType>(symbol_map, symbol);
      add_to_symbol_map<ignore_invalid_symbols, ValueType>(symbol_map,
                                                           other_symbols...);
    }


    template <typename SymbolicType, typename ValueType, typename>
    void
    set_value_in_symbol_map(types::substitution_map &substitution_map,
                            const SymbolicType &     symbol,
                            const ValueType &        value)
    {
      // Call the above function
      using SE_RCP_Basic = const SymEngine::RCP<const SymEngine::Basic> &;
      internal::set_value_in_symbol_map(substitution_map,
                                        static_cast<SE_RCP_Basic>(symbol),
                                        static_cast<SE_RCP_Basic>(
                                          SymbolicType(value)));
    }


    template <typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(types::substitution_map &        substitution_map,
                            const std::vector<SymbolicType> &symbols,
                            const std::vector<ValueType> &   values)
    {
      Assert(symbols.size() == values.size(),
             ExcDimensionMismatch(symbols.size(), values.size()));

      typename std::vector<SymbolicType>::const_iterator it_symb =
        symbols.begin();
      typename std::vector<ValueType>::const_iterator it_val = values.begin();
      for (; it_symb != symbols.end(); ++it_symb, ++it_val)
        {
          Assert(it_val != values.end(), ExcInternalError());
          set_value_in_symbol_map(substitution_map, *it_symb, *it_val);
        }
    }


    template <typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value)
    {
      set_value_in_symbol_map(substitution_map,
                              symbol_value.first,
                              symbol_value.second);
    }


    template <typename SymbolicType, typename ValueType, typename... Args>
    void
    set_value_in_symbol_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value,
      const Args &... other_symbol_values)
    {
      set_value_in_symbol_map(substitution_map, symbol_value);
      set_value_in_symbol_map(substitution_map, other_symbol_values...);
    }


    template <typename SymbolicType, typename ValueType>
    void
    set_value_in_symbol_map(
      types::substitution_map &                              substitution_map,
      const std::vector<std::pair<SymbolicType, ValueType>> &symbol_values)
    {
      // Call the above function
      for (const auto &entry : symbol_values)
        set_value_in_symbol_map(substitution_map, entry.first, entry.second);
    }


    /* ------------------ Symbol substitution map creation ----------------*/


    template <typename ExpressionType, typename ValueType, typename>
    types::substitution_map
    make_substitution_map(const ExpressionType &symbol, const ValueType &value)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbol, value);
      return substitution_map;
    }


    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(const std::vector<ExpressionType> &symbols,
                          const std::vector<ValueType> &     values)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbols, values);
      return substitution_map;
    }


    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const std::pair<ExpressionType, ValueType> &symbol_value)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbol_value);
      return substitution_map;
    }


    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    make_substitution_map(
      const std::vector<std::pair<ExpressionType, ValueType>> &symbol_values)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map, symbol_values);
      return substitution_map;
    }


    template <typename ExpressionType, typename ValueType, typename... Args>
    types::substitution_map
    make_substitution_map(
      const std::pair<ExpressionType, ValueType> &symbol_value,
      const Args &... other_symbol_values)
    {
      types::substitution_map substitution_map;
      add_to_substitution_map(substitution_map,
                              symbol_value,
                              other_symbol_values...);
      return substitution_map;
    }


    /* ---------------- Symbolic substitution map enlargement --------------*/

    namespace internal
    {
      /**
       * There are cases where it is convenient to simply ignore
       * the fact that we may be trying to add invalid symbols.
       * For example, one may wish to add a tensor where only
       * some entries are symbols and the others are zero'd.
       * In this case we ensure that the user knows what they're
       * doing by forcing them to override this safety mechanism
       * with a template parameter. These types of functions are
       * typically called in a manner indicating that the user
       * knows exactly what they're passing into it.
       */
      template <bool ignore_invalid_symbols>
      void
      add_to_substitution_map(
        types::substitution_map &                     substitution_map,
        const SymEngine::RCP<const SymEngine::Basic> &symbol,
        const SymEngine::RCP<const SymEngine::Basic> &value)
      {
        if (ignore_invalid_symbols == false)
          {
            Assert(
              internal::is_valid_substitution_symbol(*symbol),
              ExcMessage(
                "Substitution with a number that does not represent a symbol or symbolic derivative"));
            Assert(substitution_map.find(symbol) == substitution_map.end(),
                   ExcMessage("This symbol is already in the map."));
            substitution_map[Expression(symbol)] = Expression(value);
          }
        else
          {
            if (internal::is_valid_substitution_symbol(*symbol))
              {
                Assert(substitution_map.find(symbol) == substitution_map.end(),
                       ExcMessage("This symbol is already in the map."));
                substitution_map[Expression(symbol)] = Expression(value);
              }
          }
      }

    } // namespace internal


    template <bool ignore_invalid_symbols>
    void
    add_to_substitution_map(types::substitution_map &substitution_map,
                            const Expression &       symbol,
                            const Expression &       value)
    {
      internal::add_to_substitution_map<ignore_invalid_symbols>(
        substitution_map, symbol.get_RCP(), value.get_RCP());
    }


    template <bool ignore_invalid_symbols,
              typename ExpressionType,
              typename ValueType,
              typename>
    void
    add_to_substitution_map(types::substitution_map &substitution_map,
                            const ExpressionType &   symbol,
                            const ValueType &        value)
    {
      using SE_RCP_Basic = const SymEngine::RCP<const SymEngine::Basic> &;
      internal::add_to_substitution_map<ignore_invalid_symbols>(
        substitution_map,
        static_cast<SE_RCP_Basic>(symbol),
        static_cast<SE_RCP_Basic>(ExpressionType(value)));
    }


    template <bool ignore_invalid_symbols,
              typename ExpressionType,
              typename ValueType,
              typename>
    void
    add_to_substitution_map(types::substitution_map &          substitution_map,
                            const std::vector<ExpressionType> &symbols,
                            const std::vector<ValueType> &     values)
    {
      Assert(symbols.size() == values.size(),
             ExcMessage(
               "Vector of symbols and values must be of equal length."));

      typename types::symbol_vector::const_iterator   it_s = symbols.begin();
      typename std::vector<ValueType>::const_iterator it_v = values.begin();
      using SE_RCP_Basic = const SymEngine::RCP<const SymEngine::Basic> &;
      for (; it_s != symbols.end(); ++it_s, ++it_v)
        {
          Assert(it_v != values.end(), ExcInternalError());
          internal::add_to_substitution_map<ignore_invalid_symbols>(
            substitution_map,
            static_cast<SE_RCP_Basic>(*it_s),
            static_cast<SE_RCP_Basic>(ExpressionType(*it_v)));
        }
    }


    template <bool ignore_invalid_symbols,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                   substitution_map,
      const std::pair<ExpressionType, ValueType> &symbol_value)
    {
      add_to_substitution_map<ignore_invalid_symbols>(substitution_map,
                                                      symbol_value.first,
                                                      symbol_value.second);
    }


    template <bool ignore_invalid_symbols,
              typename ExpressionType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                                substitution_map,
      const std::vector<std::pair<ExpressionType, ValueType>> &symbol_values)
    {
      for (const auto &entry : symbol_values)
        {
          add_to_substitution_map<ignore_invalid_symbols>(substitution_map,
                                                          entry);
        }
    }


    template <bool ignore_invalid_symbols>
    void
    add_to_substitution_map(types::substitution_map &      substitution_map,
                            const types::substitution_map &symbol_values)
    {
      for (const auto &entry : symbol_values)
        {
          const SymEngine::RCP<const SymEngine::Basic> &symbol = entry.first;
          const SymEngine::RCP<const SymEngine::Basic> &value  = entry.second;
          internal::add_to_substitution_map<ignore_invalid_symbols>(
            substitution_map, symbol, value);
        }
    }


    template <bool ignore_invalid_symbols,
              typename ExpressionType,
              typename ValueType,
              typename... Args>
    void
    add_to_substitution_map(
      types::substitution_map &                   substitution_map,
      const std::pair<ExpressionType, ValueType> &symbol_value,
      const Args &... other_symbol_values)
    {
      add_to_substitution_map<ignore_invalid_symbols>(substitution_map,
                                                      symbol_value);
      add_to_substitution_map<ignore_invalid_symbols>(substitution_map,
                                                      other_symbol_values...);
    }


    template <typename... Args>
    types::substitution_map
    merge_substitution_maps(
      const types::substitution_map &substitution_map_in_1,
      const Args &... other_substitution_maps_in)
    {
      types::substitution_map substitution_map_out = substitution_map_in_1;
      merge_substitution_maps(substitution_map_out,
                              other_substitution_maps_in...);
      return substitution_map_out;
    }


    template <typename... Args>
    void
    merge_substitution_maps(
      types::substitution_map &      substitution_map_out,
      const types::substitution_map &substitution_map_in_1,
      const Args &... other_substitution_maps_in)
    {
      merge_substitution_maps(substitution_map_out, substitution_map_in_1);
      merge_substitution_maps(substitution_map_out,
                              other_substitution_maps_in...);
    }


    /* ---------------- Symbol substitution and evaluation --------------*/


    template <typename ExpressionType, typename ValueType>
    types::substitution_map
    resolve_explicit_dependencies(
      const std::vector<std::pair<ExpressionType, ValueType>> &symbol_values,
      const bool force_cyclic_dependency_resolution)
    {
      return resolve_explicit_dependencies(make_substitution_map(symbol_values),
                                           force_cyclic_dependency_resolution);
    }


    template <typename ValueType>
    Expression
    substitute(const Expression &expression,
               const Expression &symbol,
               const ValueType & value)
    {
      return expression.substitute(symbol, value);
    }


    template <typename ExpressionType, typename... Args>
    ExpressionType
    substitute(const ExpressionType &expression, const Args &... symbol_values)
    {
      // Call other function
      return substitute(expression, make_substitution_map(symbol_values...));
    }


    template <typename ValueType>
    ValueType
    substitute_and_evaluate(const Expression &             expression,
                            const types::substitution_map &substitution_map)
    {
      return expression.substitute_and_evaluate<ValueType>(substitution_map);
    }


    template <typename ValueType, typename... Args>
    ValueType
    substitute_and_evaluate(const Expression &expression,
                            const Args &... symbol_values)
    {
      // Call other function
      return substitute_and_evaluate<ValueType>(
        expression, make_substitution_map(symbol_values...));
    }

  } // namespace SD
} // namespace Differentiation


#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
