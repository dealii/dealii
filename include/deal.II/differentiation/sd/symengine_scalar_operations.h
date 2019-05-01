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
    namespace SE = ::SymEngine;

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
     * In most use cases the function or variable @f would be called the
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
     * @name Symbolic map creation
     */
    //@{

    namespace internal
    {
      /**
       * Return whether or not an @p entry is a valid symbol that
       * we can expect to perform substitution on.
       */
      bool
      is_valid_substitution_symbol(const SE::Basic &entry);
    } // namespace internal

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
     */
    template <typename ExpressionType = SD::Expression, typename ValueType>
    types::substitution_map
    make_substitution_map(const ExpressionType &symbol, const ValueType &value);

    /**
     * Return a substitution map that has the entry keys given by @p symbols
     * and the values given by @p values. It is expected that all key entries
     * be valid symbols or symbolic expressions.
     *
     * It is possible to map symbolic types to other symbolic types
     * using this function. For more details on this, see the other
     * \ref make_substitution_map(const Expression &,const ValueType &)
     * function.
     */
    template <typename ExpressionType = SD::Expression, typename ValueType>
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
     * \ref make_substitution_map(const Expression &,const ValueType &)
     * function.
     */
    template <typename ExpressionType = SD::Expression, typename ValueType>
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
     * \ref make_substitution_map(const Expression &,const ValueType &)
     * function.
     */
    template <typename ExpressionType = SD::Expression, typename ValueType>
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
     * \ref make_substitution_map(const Expression &,const ValueType &)
     * function.
     */
    template <typename SymbolicType, typename ValueType, typename... Args>
    types::substitution_map
    make_substitution_map(
      const std::pair<SymbolicType, ValueType> &symbol_value,
      const Args &... other_symbol_values);

    //@}

    /**
     * @name Symbolic substitution map enlargement
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
      add_to_substitution_map(types::substitution_map &       substitution_map,
                              const SE::RCP<const SE::Basic> &symbol,
                              const SE::RCP<const SE::Basic> &value);
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
     * 2. it is convertible to a `const SE::RCP<const SE::Basic> &`.
     *
     * For more context which this function is used, see the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &, const Expression &) function.
     *
     * @tparam ignore_invalid_symbols See the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &,const Expression &) function for a detailed
     * discussion on the role of this template argument.
     */
    template <bool ignore_invalid_symbols = false,
              typename ExpressionType     = SD::Expression,
              typename ValueType,
              typename = typename std::enable_if<
                dealii::internal::is_explicitly_convertible<
                  ExpressionType,
                  const SE::RCP<const SE::Basic> &>::value &&
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
     * 2. it is convertible to a `const SE::RCP<const SE::Basic> &`.
     *
     * For more context which this function is used, see the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &, const Expression &) function.
     *
     * @tparam ignore_invalid_symbols See the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &,const Expression &) function for a detailed
     * discussion on the role of this template argument.
     */
    template <bool ignore_invalid_symbols = false,
              typename ExpressionType     = SD::Expression,
              typename ValueType,
              typename = typename std::enable_if<
                dealii::internal::is_explicitly_convertible<
                  ExpressionType,
                  const SE::RCP<const SE::Basic> &>::value &&
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
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &, const Expression &) function.
     *
     * @tparam ignore_invalid_symbols See the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &,const Expression &) function for a detailed
     * discussion on the role of this template argument.
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
     * The @p SymbolicType and its associated @ValueType need not be scalar
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
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &, const Expression &) function.
     *
     * @tparam ignore_invalid_symbols See the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &,const Expression &) function for a detailed
     * discussion on the role of this template argument.
     */
    template <bool ignore_invalid_symbols = false,
              typename SymbolicType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value);

    /**
     * A convenience function for adding multiple entries to the
     * @p substitution_map. The new entries will have the keys given by first
     * entry of each element of @p symbol_values, and the values given its
     * second entry. It is expected that the key entry be a valid symbols or
     * symbolic expressions, and that the paired @p symbol_value elements are
     * compatible with the other add_to_substitution_map() functions.
     *
     * The @p SymbolicType and its associated @ValueType need not be scalar
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
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &, const Expression &) function.
     *
     * @tparam ignore_invalid_symbols See the other
     * \ref add_to_substitution_map(types::substitution_map &,const
     * Expression &,const Expression &) function for a detailed
     * discussion on the role of this template argument.
     */
    template <bool ignore_invalid_symbols = false,
              typename SymbolicType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                              substitution_map,
      const std::vector<std::pair<SymbolicType, ValueType>> &symbol_values);

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
     * \ref make_substitution_map(const Expression &,const ValueType &)
     * function.
     */
    template <bool ignore_invalid_symbols = false,
              typename SymbolicType,
              typename ValueType,
              typename... Args>
    void
    add_to_substitution_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value,
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

  } // namespace SD
} // namespace Differentiation


/* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN


namespace Differentiation
{
  namespace SD
  {
    /* ---------------- Symbolic substitution map creation --------------*/


    template <typename ValueType>
    types::substitution_map
    make_substitution_map(const Expression &symbol, const ValueType &value)
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


    template <typename SymbolicType, typename ValueType, typename... Args>
    types::substitution_map
    make_substitution_map(
      const std::pair<SymbolicType, ValueType> &symbol_value,
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
      add_to_substitution_map(types::substitution_map &       substitution_map,
                              const SE::RCP<const SE::Basic> &symbol,
                              const SE::RCP<const SE::Basic> &value)
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
      using SE_RCP_Basic = const SE::RCP<const SE::Basic> &;
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
      using SE_RCP_Basic = const SE::RCP<const SE::Basic> &;
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
              typename SymbolicType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value)
    {
      add_to_substitution_map<ignore_invalid_symbols>(substitution_map,
                                                      symbol_value.first,
                                                      symbol_value.second);
    }


    template <bool ignore_invalid_symbols,
              typename SymbolicType,
              typename ValueType>
    void
    add_to_substitution_map(
      types::substitution_map &                              substitution_map,
      const std::vector<std::pair<SymbolicType, ValueType>> &symbol_values)
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
          const SE::RCP<const SE::Basic> &symbol = entry.first;
          const SE::RCP<const SE::Basic> &value  = entry.second;
          internal::add_to_substitution_map<ignore_invalid_symbols>(
            substitution_map, symbol, value);
        }
    }


    template <bool ignore_invalid_symbols,
              typename SymbolicType,
              typename ValueType,
              typename... Args>
    void
    add_to_substitution_map(
      types::substitution_map &                 substitution_map,
      const std::pair<SymbolicType, ValueType> &symbol_value,
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

  } // namespace SD
} // namespace Differentiation


#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
