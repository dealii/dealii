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

#ifndef dealii_differentiation_sd_symengine_number_types_h
#define dealii_differentiation_sd_symengine_number_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

// Low level
#  include <symengine/basic.h>
#  include <symengine/dict.h>
#  include <symengine/symengine_exception.h>
#  include <symengine/symengine_rcp.h>

// Number types
#  include <symengine/expression.h>
#  include <symengine/integer.h>
#  include <symengine/logic.h>
#  include <symengine/number.h>
#  include <symengine/rational.h>

// Number operations
#  include <symengine/add.h>
#  include <symengine/functions.h>
#  include <symengine/mul.h>
#  include <symengine/pow.h>

// Evaluation
#  include <symengine/eval.h>
#  include <symengine/eval_arb.h>
#  include <symengine/eval_double.h>

// Differentiation
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/numbers.h>

#  include <deal.II/differentiation/sd/symengine_number_traits.h>
#  include <deal.II/differentiation/sd/symengine_types.h>

#  include <boost/serialization/split_member.hpp>

#  include <symengine/derivative.h>

#  include <algorithm>
#  include <memory>
#  include <sstream>
#  include <type_traits>
#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * An exception to indicate that the string sent to the SymEngine
     * parser is not valid.
     */
    DeclException1(
      ExcSymEngineParserError,
      std::string,
      << "The string '" << arg1
      << "' could not be parsed successfully. Are you sure that (1) it "
      << "consists of legitimate operations and syntax, and (2) you've "
      << "previously declared all symbolic variables that are present "
      << "in the expression?");

    //@}

    /**
     * A class to wrap SymEngine expressions.
     *
     * With this number type, SymEngine numbers can be used to perform scalar
     * and tensor mathematics in deal.II. It (or the SymEngine::Expression
     * class, of which it stores an instance) therefore forms the basis of
     * symbolic computation via SymEngine in deal.II. With it one can
     * perform symbolic differentiation and subsequent substitution with both
     * scalars and deal.II's native Tensor and SymmetricTensor types.
     *
     * The symbolic features that this class supports includes:
     * - expression parsing,
     * - comparison operations,
     * - logical operations,
     * - math operations,
     * - conditional expression construction,
     * - differentiation,
     * - substitution (partial and complete), and
     * - serialization.
     *
     * A simple example of how this class may be used is as follows:
     * @code
     *   // Constructing a symbolic expression:
     *   // This is a symbol, which we will treat as an argument to a symbolic
     *   // function.
     *   const Expression x("x");
     *   const Expression y("y");
     *   // This is a symbolic expression, which is an expression constructed
     *   // from individual symbols.
     *   const Expression f = (x + y)*(x + y);
     *
     *   // Value substitution
     *   types::substitution_map substitution_map;
     *   substitution_map[x] = Expression(1);
     *   substitution_map[y] = Expression(2.5);
     *   const double evaluated_f =
     *     f.substitute_and_evaluate<double>(substitution_map);
     *   // We could also have performed substitution of each individual
     *   // argument, if we wanted to. This means that one can partially
     *   // substitute an expression at any time.
     * @endcode
     *
     * A more intricate example of conditional evaluation is as follows:
     * @code
     *   // Construct symbolic expressions
     *   const Expression x("x");
     *   const Expression y("y");
     *   const Expression f_plus = (x + y)*(x + y);
     *   // Parsing expressions from a string is also possible. Its arguments
     *   // must have been previously declared through (and in scope).
     *   const Expression f_minus ("(x-y)*(x-y)", true);
     *
     *   // Constructing a conditional expression
     *   const SD_number_t f((x > Expression(0.0)), f_plus, f_minus);
     *
     *   // Value substitution
     *   types::substitution_map substitution_map;
     *   substitution_map[x] = Expression(1);
     *   substitution_map[y] = Expression(2.5);
     *   const double evaluated_f =
     *     f.substitute_and_evaluate<double>(substitution_map);
     *   // Since the substituted value for x was greater than zero, we expect
     *   // that the returned result now in evaluated_f was evaluated from
     *   // the function f_plus.
     * @endcode
     *
     * Lastly, here is an example using symbolic differentiation:
     * @code
     *   // Construct symbolic expressions
     *   const Expression x("x");
     *   const Expression f("x**2", true);
     *
     *   // Now perform differentiation. Specifically, we differentiate the
     *   // function "f" with respect to the symbolic variable "x".
     *   // The result should be the expression "2*x".
     *   const Expression df_dx = f.differentiate(x);
     *
     *   // Value substitution
     *   types::substitution_map substitution_map;
     *   substitution_map[x] = Expression(10.0);
     *   const double evaluated_df_dx =
     *     evaluated_df_dx.substitute_and_evaluate<double>(substitution_map);
     *   // We can expect the above to evaluate to "2*10" which is,
     *   // of course, the numeric value 20.
     * @endcode
     *
     * @author Jean-Paul Pelteret, 2019
     */
    class Expression
    {
    public:
      /**
       * @name Constructors
       */
      //@{

      /**
       * Default constructor.
       */
      Expression();

      /**
       * Constructor for boolean types.
       *
       * @note This constructor is marked as explicit so that there are no
       * potential ambiguities related to implicit conversions in either user
       * code or math functions that are loaded into the standard namespace.
       */
      explicit Expression(const bool value);

      /**
       * Constructor for arithmetic number types.
       *
       * @note This constructor is marked as explicit so that there are no
       * potential ambiguities related to implicit conversions in either user
       * code or math functions that are loaded into the standard namespace.
       */
      template <typename NumberType,
                typename = typename std::enable_if<
                  std::is_arithmetic<NumberType>::value>::type>
      explicit Expression(const NumberType &value);

      /**
       * Constructor for complex numbers templated on arithmetic number types.
       *
       * @note This constructor is marked as explicit so that there are no
       * potential ambiguities related to implicit conversions in either user
       * code or math functions that are loaded into the standard namespace.
       */
      template <typename NumberType,
                typename = typename std::enable_if<
                  std::is_arithmetic<NumberType>::value>::type>
      explicit Expression(const std::complex<NumberType> &value);

      /**
       * Constructor for integer types.
       */
      Expression(const SymEngine::integer_class &value);

      /**
       * Constructor for rational types.
       *
       * It is expected that both the @p numerator and @p denominator
       * be integral types.
       */
      template <typename NumberType,
                typename = typename std::enable_if<
                  std::is_integral<NumberType>::value>::type>
      Expression(const NumberType &numerator, const NumberType &denominator);

      /**
       * Constructor for rational types.
       */
      Expression(const SymEngine::rational_class &value);

      /**
       * Constructor for a piecewise defined function.
       *
       * The generated expression may be interpreted as the result of the
       * teniary operator, i.e.
       * <code>(condition ? expression_if_true : expression_if_false)</code>.
       *
       * The @p condition can be any expression that renders an expression that
       * is convertible to a SymEngine::Boolean operator. This includes:
       * - the logical operators of the deal.II Expression class
       * - SymEngine::boolean()
       * - SymEngine::contains()
       * - SymEngine::Eq()
       * - SymEngine::Ne()
       * - SymEngine::Ge()
       * - SymEngine::Gt()
       * - SymEngine::Le()
       * - SymEngine::Lt()
       * - SymEngine::logical_and()
       * - SymEngine::logical_nand()
       * - SymEngine::logical_or()
       * - SymEngine::logical_not()
       * - SymEngine::logical_nor()
       * - SymEngine::logical_xor()
       * - SymEngine::logical_xnor()
       * - ...
       *
       * An example of this constructor's use is as follows:
       * @code
       *   const Expression x("x");
       *   const Expression y("y");
       *
       *   // Construct a conditional expression using the symbolic variables.
       *   const Expression f ((x < Expression(0.0)), x+y, x-y);
       * @endcode
       */
      Expression(const Expression &condition,
                 const Expression &expression_if_true,
                 const Expression &expression_if_false);

      /**
       * Constructor for a piecewise defined function.
       *
       * The generated expression may be interpreted as the result of the
       * set of nested if-elseif-else statements, i.e. (in pseudo-code)
       * @code
       *   if (condition_expression[0].first == true)
       *     return condition_expression[0].second;
       *   else if (condition_expression[1].first == true)
       *     return condition_expression[1].second;
       *   else if (...)
       *     return ...;
       *   else
       *     return expression_otherwise;
       * @endcode
       * if the input vector has more than 2 elements.
       *
       * This variant takes the piecewise evaluated conditions and its results
       * as the first argument, and the default return value as the second
       * argument.
       */
      Expression(const std::vector<std::pair<Expression, Expression>>
                   &               condition_expression,
                 const Expression &expression_otherwise);

      /**
       * Constructor for a piecewise defined function.
       *
       * The generated expression may be interpreted as the result of the
       * set of nested if-elseif statements, i.e. (in pseudo-code)
       * @code
       *   if (condition_expression[0].first == true)
       *     return condition_expression[0].second;
       *   else if (condition_expression[1].first == true)
       *     return condition_expression[1].second;
       *   else if (...)
       *     return ...;
       * @endcode
       * if the input vector has more than 2 elements.
       *
       * This variant takes only the piecewise evaluated conditions and its
       * results. If none of the conditions are met upon evaluation then the
       * returned result will be NaN.
       */
      Expression(const std::vector<std::pair<Expression, Expression>>
                   &condition_expression);


      /**
       * Constructor for symbolic types.
       *
       * This constructor initializes a symbolic type with a character array
       * representing its symbolic value.
       */
      Expression(const char *symbol);

      /**
       * Constructor for symbolic types.
       *
       * This constructor initializes a symbolic type with a string
       * representing its symbolic value.
       * If the @p parse_as_expression flag is <tt>false</tt>, then the
       * @p symb_expr (potentially composed of multiple characters) will be
       * interpreted as a single symbol.
       * If the @p parse_as_expression flag is <tt>true</tt>, then the
       * @p symb_expr will be parsed as a symbolic expression (potentially
       * composed of multiple symbols, constants, etc.).
       */
      Expression(const std::string &symb_expr,
                 const bool         parse_as_expression = false);

      /**
       * Constructor for function symbol types.
       *
       * This constructor initializes a function symbol with a string
       * representing its symbolic name.
       */
      Expression(const std::string &         symbol_func,
                 const types::symbol_vector &arguments);

      /**
       * Copy constructor.
       */
      Expression(const Expression &rhs) = default;

      /**
       * Copy constructor.
       *
       * @note This constructor is marked as explicit to prevent any ambiguities
       * from implicit conversion when both the deal.II and SymEngine namespaces
       * are imported.
       */
      explicit Expression(const SymEngine::Expression &rhs);

      /**
       * Copy constructor.
       *
       * This allows us to create our class straight from the results of
       * SymEngine operations. This is especially important for operations like
       * "diff", because the returned result is not primitive, but rather a set
       * of compound operations.
       */
      Expression(const SymEngine::RCP<const SymEngine::Basic> &rhs);

      /**
       * Move constructor.
       */
      Expression(Expression &&rhs) = default;

      /**
       * Move constructor.
       *
       * This allows us to create our class straight from the results of
       * SymEngine operations. This is especially important for operations like
       * "diff", because the returned result is not primitive, but rather a set
       * of compound operations.
       */
      Expression(SymEngine::RCP<const SymEngine::Basic> &&rhs);

      /**
       * Destructor.
       */
      virtual ~Expression() = default;

      //@}

      /**
       * Utilities
       */
      //@{

      /**
       * Parse an expression from a string representing a symbolic @p expression.
       * This overwrites any existing value or expression that this object
       * represents.
       */
      Expression &
      parse(const std::string &expression);

      /**
       * Print the value stored by this object.
       *
       * Since the stored value could be one of a number of types, we leave
       * SymEngine to cast and output the correct representation of the data.
       */
      std::ostream &
      print(std::ostream &stream) const;

      /**
       * Save the value stored by this object to the @p stream.
       *
       * Each expression will be saved on a new line of the @p stream.
       */
      void
      save(std::ostream &stream) const;

      /**
       * Load the value stored in the @p stream into this object.
       *
       * It is expected that each expression appears on its own on single
       * line of @p stream.
       *
       * @note When loading a symbolic expression, it is imperative that
       * you first create or load all of the symbolic variables used in
       * the saved expression.
       */
      void
      load(std::istream &stream);

      /**
       * Write the data of this object to a stream for the purpose
       * of serialization.
       *
       * This effectively saves the value stored in this object to the
       * @p archive with the given @p version number.
       */
      template <class Archive>
      void
      save(Archive &archive, const unsigned int version) const;

      /**
       * Read the data of this object from a stream for the purpose
       * of serialization.
       *
       * This effectively loads into this object the value stored in of the
       * @p archive with the given @p version number.
       * In doing so, the previous contents of this object are thrown away.
       *
       * @note When deserializing a symbolic expression, it is imperative that
       * you first create or deserialize all of the symbolic variables used in
       * the serialized expression.
       */
      template <class Archive>
      void
      load(Archive &archive, const unsigned int version);

#  ifdef DOXYGEN
      /**
       * Write and read the data of this object from a stream for the purpose
       * of serialization.
       *
       * This effectively saves or loads the value stored into/out of the
       * @p archive with the given @p version number into this object.
       * If deserializing data, then the previous contents of this object
       * are thrown away.
       *
       * @note When deserializing a symbolic expression, it is imperative that
       * you first create or deserialize all of the symbolic variables used in
       * the serialized expression.
       */
      template <class Archive>
      void
      serialize(Archive &archive, const unsigned int version);
#  else
      // This macro defines the serialize() method that is compatible with
      // the templated save() and load() method that have been implemented.
      BOOST_SERIALIZATION_SPLIT_MEMBER()
#  endif

      //@}

      /**
       * @name Values
       */
      //@{

      /**
       * Return the value or expression that this class instance represents.
       */
      const SymEngine::Expression &
      get_expression() const;

      /**
       * Return the primitive SymEngine data type that stores the value or
       * expression represented by this object.
       */
      const SymEngine::Basic &
      get_value() const;

      /**
       * Return the pointer to the primitive SymEngine data type that stores
       * the value or expression represented by this object.
       */
      const SymEngine::RCP<const SymEngine::Basic> &
      get_RCP() const;

      //@}

      /**
       * @name Math and relational operators with (potentially) symbolic types
       */
      //@{

      /**
       * Assignment operator.
       *
       * Sets the data of this object's @p expression equal
       * to that of the @p rhs object.
       */
      Expression &
      operator=(const Expression &rhs);

      /**
       * Assignment operator.
       *
       * Sets the data of this object's @p expression equal
       * to that of the @p rhs object.
       */
      Expression &
      operator=(Expression &&rhs) noexcept;

      /**
       * Addition assignment.
       *
       * The @p rhs @p expression is added in-place to that of
       * this object's @p expression.
       */
      Expression &
      operator+=(const Expression &rhs);

      /**
       * Subtraction assignment.
       *
       * The @p rhs @p expression is subtracted in-place from that of
       * this object's @p expression.
       */
      Expression &
      operator-=(const Expression &rhs);

      /**
       * Multiplication assignment.
       *
       * This object's @p expression is multiplied in-place by that of
       * the @p rhs @p expression.
       */
      Expression &
      operator*=(const Expression &rhs);

      /**
       * Division assignment.
       *
       * This object's @p expression is divided in-place by that of
       * the @p rhs @p expression.
       */
      Expression &
      operator/=(const Expression &rhs);

      //@}

      /**
       * @name Math and relational operators with numeric types
       */
      //@{

      /**
       * Assignment operator.
       *
       * Set the data of this object's @p expression equal
       * to the numerical value of the @p rhs.
       */
      template <typename NumberType>
      Expression &
      operator=(const NumberType &rhs);

      /**
       * Negation operator.
       *
       * Return a the result of pre-multipying this object's @p expression
       * by <tt>-1</tt>.
       *
       * @note This operation is not performed in-place.
       */
      Expression
      operator-() const;

      /**
       * Addition assignment.
       *
       * The numerical value of the @p rhs is added in-place to that of
       * this object's @p expression.
       */
      template <typename NumberType>
      Expression &
      operator+=(const NumberType &rhs);

      /**
       * Subtraction assignment.
       *
       * The numerical value of the @p rhs is subtracted in-place from that of
       * this object's @p expression.
       */
      template <typename NumberType>
      Expression &
      operator-=(const NumberType &rhs);

      /**
       * Multiplication assignment.
       *
       * This object's @p expression is multiplied in-place by that of
       * the numerical value of the @p rhs.
       */
      template <typename NumberType>
      Expression &
      operator*=(const NumberType &rhs);

      /**
       * Division assignment.
       *
       * This object's @p expression is divided in-place by that of
       * the numerical value of the @p rhs.
       */
      template <typename NumberType>
      Expression &
      operator/=(const NumberType &rhs);

      //@}

      /**
       * @name Differentiation
       */
      //@{

      /**
       * Return the derivative of this object's @p expression
       * with respect to the given @p symbol.
       */
      Expression
      differentiate(const Expression &symbol) const;

      /**
       * Return the derivative of this object's @p expression
       * with respect to the given @p symbol.
       */
      Expression
      differentiate(
        const SymEngine::RCP<const SymEngine::Symbol> &symbol) const;

      /**
       * Return the derivative of this object's @p expression
       * with respect to the potential @p symbol.
       */
      Expression
      differentiate(const SymEngine::RCP<const SymEngine::Basic> &symbol) const;

      //@}

      /**
       * @name Dictionary-based substitution
       */
      //@{

      /**
       * Perform substitution of all symbols found in this object's @p expression
       * that match a key in the @p substitution_values map.
       *
       * @note The replacement value (the entry in the @p substitution_values
       * that is paired with a key) need not necessarily be numerical, but may
       * also be another symbolic type.
       *
       * @note With dictionary substitution, partial substitution is allowed
       * (i.e. an incomplete substitution map can be used and the return type
       * can be symbolic).
       */
      Expression
      substitute(const types::substitution_map &substitution_values) const;

      /**
       * Perform substitution of all symbols found in this object's @p expression
       * that match a key in the @p substitution_values map.
       *
       * This function is like the one above, but takes in a SymEngine map
       * (one that maps a `SymEngine::RCP<const SymEngine::Basic>` to another
       * `SymEngine::RCP<const SymEngine::Basic>`) as an argument.
       *
       * @note The replacement value (the entry in the @p substitution_values
       * that is paired with a key) need not necessarily be numerical, but may
       * also be another symbolic type.
       *
       * @note With dictionary substitution, partial substitution is allowed
       * (i.e. an incomplete substitution map can be used and the return type
       * can be symbolic).
       */
      Expression
      substitute(const SymEngine::map_basic_basic &substitution_values) const;

      /**
       * Perform substitution of all symbols found in this object's @p expression
       * that match the @p symbol. Each @p symbol will be substituted with
       * the given @p value.
       *
       * @note With dictionary substitution, partial substitution is allowed
       * (i.e. an incomplete substitution map can be used and the return type
       * can be symbolic).
       */
      Expression
      substitute(const Expression &symbol, const Expression &value) const;

      /**
       * Perform substitution of all symbols found in this object's @p expression
       * that match the @p symbol. Each @p symbol will be substituted with
       * the given @p value.
       *
       * @note With dictionary substitution, partial substitution is allowed
       * (i.e. an incomplete substitution map can be used and the return type
       * can be symbolic).
       */
      template <typename NumberType>
      Expression
      substitute(const Expression &symbol, const NumberType &value) const;

      /**
       * Full substitution and evaluation. This creates a Expression by
       * symbol substitution and then immediately computes its numerical value.
       *
       * @note All symbols must be resolved by the substitution map in order
       * for this function to return successfully.
       */
      template <typename ReturnType>
      ReturnType
      substitute_and_evaluate(
        const types::substitution_map &substitution_values) const;

      /**
       * Full substitution and evaluation. This creates a Expression by
       * symbol substitution and then immediately computes its numerical value.
       *
       * This function is like the one above, but takes in a SymEngine map
       * (one that maps a `SymEngine::RCP<const SymEngine::Basic>` to another
       * `SymEngine::RCP<const SymEngine::Basic>`) as an argument.
       *
       * @note All symbols must be resolved by the substitution map in order
       * for this function to return successfully.
       */
      template <typename ReturnType>
      ReturnType
      substitute_and_evaluate(
        const SymEngine::map_basic_basic &substitution_values) const;

      //@}

      /**
       * @name Conversion operators
       */
      //@{

      /**
       * Conversion operator for real integer or floating point values, and
       * complex integer or floating point values.
       *
       * @note This function is marked explicit so that the conversion must be
       * performed using a static_cast. In normal use, one would have expected
       * (Expression)*(double) --> (Expression)
       * If this function were not marked as explicit, then we could potentially
       * have
       * (Expression)*(double) --> (double)
       * So, to get out a value on needs to do the following:
       *
       * <code>
       *   const NumberType val = static_cast<NumberType>(Expression);
       * </code>
       *
       * or, probably less desirably,
       *
       * <code>
       *   const NumberType val = NumberType(Expression);
       * </code>
       *
       * @note If the underlying number is a custom type (i.e. encapsulated by a
       * NumberWrapper), then it is necessary to derive a new class from
       * Expression and define a specialized conversion operator
       * that calls an Evaluator that is specialized for this custom number
       * type. This could be achieved with an overriding conversion function
       * in the base class, for example:
       *
       * @code
       *   class MyNumber : public Expression
       *   {
       *     ...
       *     template <typename ResultType>
       *     explicit operator ResultType() const
       *     {
       *       if (this->get_value()->get_type_code() ==
       * SymEngine::NUMBER_WRAPPER)
       *       {
       *         // Implement custom evaluation function
       *         const ResultType result = ...;
       *         return result;
       *       }
       *       else
       *         // Call base class conversion operator
       *         return Expression::operator ResultType();
       *     }
       *   }
       * @endcode
       */
      template <typename ResultType>
      explicit operator ResultType() const;

      /**
       * Conversion operator that returns the value or expression that this
       * class instance represents.
       */
      explicit operator const SymEngine::Expression &() const;

      /**
       * Conversion operator that returns a SymEngine reference counted pointer
       * to the fundamental type.
       */
      operator const SymEngine::RCP<const SymEngine::Basic> &() const;

      //@}

    protected:
      /**
       * Return the value or expression that this class instance represents.
       */
      SymEngine::Expression &
      get_expression();

    private:
      /**
       * The value or expression that this instance of this class is to
       * represent.
       */
      SymEngine::Expression expression;
    };

    /**
     * @name Type traits
     */
    //@{

    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * symbolically differentiable number or not.
     * This is a specialization for the deal.II Expression class.
     *
     * @author Jean-Paul Pelteret, 2019
     */
    template <>
    struct is_symengine_number<Expression> : std::true_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * SymEngine number or not.
     * This is a specialization for the deal.II Expression class.
     *
     * @author Jean-Paul Pelteret, 2019
     */
    template <>
    struct is_sd_number<Expression> : std::true_type
    {};

    //@}

    /**
     * @name Bitwise operators
     */
    //@{

    /**
     * Bitwise left shift operator.
     *
     * This is used to output the @p expression to the input @p stream.
     */
    std::ostream &
    operator<<(std::ostream &stream, const Expression &expression);

    /**
     * Bitwise right shift operator.
     *
     * This is used to read in an @p expression from the input @p stream.
     */
    std::istream &
    operator>>(std::istream &stream, Expression &expression);

    //@}

    /**
     * @name Comparison operators
     */
    //@{

    /**
     * Equality operator.
     *
     * Return whether the @p lhs is equal to the @p rhs.
     */
    Expression
    operator==(const Expression &lhs, const Expression &rhs);

    /**
     * Non-equality operator.
     *
     * Return whether the @p lhs is not equal to the @p rhs.
     */
    Expression
    operator!=(const Expression &lhs, const Expression &rhs);

    /**
     * Less than operator.
     *
     * Return whether the @p lhs is less than the @p rhs.
     */
    Expression
    operator<(const Expression &lhs, const Expression &rhs);

    /**
     * Greater than operator.
     *
     * Return whether the @p lhs is greater than the @p rhs.
     */
    Expression
    operator>(const Expression &lhs, const Expression &rhs);

    /**
     * Less than or equals operator.
     *
     * Return whether the @p lhs is less than, or equal to, the @p rhs.
     */
    Expression
    operator<=(const Expression &lhs, const Expression &rhs);

    /**
     * Greater than or equals operator.
     *
     * Return whether the @p lhs is greater than, or equal to, the @p rhs.
     */
    Expression
    operator>=(const Expression &lhs, const Expression &rhs);

    //@}

    /**
     * @name Logical operators
     */
    //@{

    /**
     * Logical not operator.
     *
     * @note This operator can only be applied on boolean and conditional
     * expressions.
     */
    Expression operator!(const Expression &expression);

    /**
     * Logical and operator.
     *
     * @note This operator can only be applied on boolean and conditional
     * expressions.
     */
    Expression operator&(const Expression &lhs, const Expression &rhs);

    /**
     * Logical inclusive or operator.
     *
     * @note This operator can only be applied on boolean and conditional
     * expressions.
     */
    Expression
    operator|(const Expression &lhs, const Expression &rhs);

    /**
     * Logical exclusive or (xor) operator.
     *
     * @note This operator can only be applied on boolean and conditional
     * expressions.
     */
    Expression
    operator^(const Expression &lhs, const Expression &rhs);

    /**
     * And operator.
     *
     * @note This operator can only be applied on boolean and conditional
     * expressions.
     * This operator is a convenience wrapper for the logical and operator.
     */
    Expression
    operator&&(const Expression &lhs, const Expression &rhs);

    /**
     * Inclusive or operator.
     *
     * @note This operator can only be applied on boolean and conditional
     * expressions.
     * This operator is a convenience wrapper for the logical or operator.
     */
    Expression
    operator||(const Expression &lhs, const Expression &rhs);

    //@}

    /**
     * @name Mathematical operators
     */
    //@{

    /**
     * Addition operator.
     *
     * Return the result of adding the @p rhs to the @p lhs.
     */
    Expression
    operator+(Expression lhs, const Expression &rhs);

    /**
     * Subtraction operator.
     *
     * Return the result of subtracting the @p rhs from the @p lhs.
     */
    Expression
    operator-(Expression lhs, const Expression &rhs);

    /**
     * Multiplication operator.
     *
     * Return the result of multiplying the @p lhs by the @p rhs.
     */
    Expression operator*(Expression lhs, const Expression &rhs);

    /**
     * Division operator.
     *
     * Return the result of dividing the @p lhs by the @p rhs.
     */
    Expression
    operator/(Expression lhs, const Expression &rhs);

    /**
     * General addition operator.
     *
     * Return the result of adding the @p rhs to the @p lhs.
     * The @p lhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression
    operator+(const NumberType &lhs, const Expression &rhs)
    {
      return Expression(lhs) + rhs;
    }

    /**
     * General addition operator.
     *
     * Return the result of adding the @p rhs to the @p lhs.
     * The @p rhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression
    operator+(const Expression &lhs, const NumberType &rhs)
    {
      return lhs + Expression(rhs);
    }

    /**
     * General subtraction operator.
     *
     * Return the result of subtracting the @p rhs from the @p lhs.
     * The @p lhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression
    operator-(const NumberType &lhs, const Expression &rhs)
    {
      return Expression(lhs) - rhs;
    }

    /**
     * General subtraction operator.
     *
     * Return the result of subtracting the @p rhs from the @p lhs.
     * The @p rhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression
    operator-(const Expression &lhs, const NumberType &rhs)
    {
      return lhs - Expression(rhs);
    }

    /**
     * General multiplication operator.
     *
     * Return the result of multiplying the @p lhs by the @p rhs.
     * The @p lhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression operator*(const NumberType &lhs, const Expression &rhs)
    {
      return Expression(lhs) * rhs;
    }

    /**
     * General multiplication operator.
     *
     * Return the result of multiplying the @p lhs by the @p rhs.
     * The @p rhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression operator*(const Expression &lhs, const NumberType &rhs)
    {
      return lhs * Expression(rhs);
    }

    /**
     * General division operator.
     *
     * Return the result of dividing the @p lhs by the @p rhs.
     * The @p lhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression
    operator/(const NumberType &lhs, const Expression &rhs)
    {
      return Expression(lhs) / rhs;
    }

    /**
     * General division operator.
     *
     * Return the result of dividing the @p lhs by the @p rhs.
     * The @p lhs type may be any supported number type, and the result
     * is promoted to a Expression. The type conversion makes writing
     * scalar expressions using Expression more natural.
     */
    template <typename NumberType,
              typename = typename std::enable_if<
                std::is_constructible<Expression, NumberType>::value>::type>
    inline Expression
    operator/(const Expression &lhs, const NumberType &rhs)
    {
      return lhs / Expression(rhs);
    }

    //@}

  } // namespace SD
} // namespace Differentiation


/* ---------------------- inline and template functions -------------------- */


#  ifndef DOXYGEN


namespace Differentiation
{
  namespace SD
  {
    template <typename NumberType, typename>
    Expression::Expression(const NumberType &value)
      : expression(value)
    {}


    template <typename NumberType, typename>
    Expression::Expression(const std::complex<NumberType> &value)
      : expression(value)
    {}


    template <typename NumberType, typename>
    Expression::Expression(const NumberType &numerator,
                           const NumberType &denominator)
      : expression(
          SymEngine::Rational::from_two_ints(*SymEngine::integer(numerator),
                                             *SymEngine::integer(denominator)))
    {}


    template <class Archive>
    void
    Expression::save(Archive &ar, const unsigned int /*version*/) const
    {
      std::stringstream sstream;
      sstream << *this;
      const std::string expr = sstream.str();
      ar &              expr;
    }


    template <class Archive>
    void
    Expression::load(Archive &ar, const unsigned int /*version*/)
    {
      std::string expr;
      ar &        expr;
      parse(expr);
    }


    template <typename NumberType>
    Expression
    Expression::substitute(const Expression &symbol,
                           const NumberType &value) const
    {
      Assert(SymEngine::is_a<SymEngine::Symbol>(symbol.get_value()),
             ExcMessage(
               "Substitution with a number that does not represent a symbol."));

      types::substitution_map sub_vals;
      sub_vals[symbol] = Expression(value);
      return substitute(sub_vals);
    }


    template <typename ReturnType>
    ReturnType
    Expression::substitute_and_evaluate(
      const types::substitution_map &substitution_values) const
    {
      return static_cast<ReturnType>(substitute(substitution_values));
    }


    template <typename ReturnType>
    ReturnType
    Expression::substitute_and_evaluate(
      const SymEngine::map_basic_basic &substitution_values) const
    {
      return static_cast<ReturnType>(substitute(substitution_values));
    }


    template <typename NumberType>
    Expression &
    Expression::operator=(const NumberType &rhs)
    {
      *this = Expression(rhs);
      return *this;
    }


    template <typename NumberType>
    Expression &
    Expression::operator+=(const NumberType &rhs)
    {
      *this = Expression(SymEngine::add(get_RCP(), Expression(rhs).get_RCP()));
      return *this;
    }


    template <typename NumberType>
    Expression &
    Expression::operator-=(const NumberType &rhs)
    {
      *this = Expression(SymEngine::sub(get_RCP(), Expression(rhs).get_RCP()));
      return *this;
    }


    template <typename NumberType>
    Expression &
    Expression::operator*=(const NumberType &rhs)
    {
      *this = Expression(SymEngine::mul(get_RCP(), Expression(rhs).get_RCP()));
      return *this;
    }


    template <typename NumberType>
    Expression &
    Expression::operator/=(const NumberType &rhs)
    {
      *this = Expression(SymEngine::div(get_RCP(), Expression(rhs).get_RCP()));
      return *this;
    }


    template <typename ResultType>
    Expression::operator ResultType() const
    {
      return static_cast<ResultType>(get_expression());
    }

  } // namespace SD
} // namespace Differentiation


// Specializations from numbers.h
// These are required in order to make the SD::Expression class a compatible
// number for the Tensor and SymmetricTensor classes

// Forward declarations:
template <int rank_, int dim, typename Number>
class Tensor;

namespace internal
{
  template <>
  struct NumberType<Differentiation::SD::Expression>
  {
    static const Differentiation::SD::Expression &
    value(const Differentiation::SD::Expression &t)
    {
      return t;
    }

    template <
      typename T,
      typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    static Differentiation::SD::Expression
    value(const T &t)
    {
      return Differentiation::SD::Expression(t);
    }

    template <
      typename T,
      typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    static Differentiation::SD::Expression
    value(T &&t)
    {
      return Differentiation::SD::Expression(t);
    }

    template <
      typename T,
      typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    static Differentiation::SD::Expression
    value(const std::complex<T> &t)
    {
      return Differentiation::SD::Expression(t);
    }

    template <typename T, int dim>
    static Differentiation::SD::Expression
    value(const Tensor<0, dim, T> &t)
    {
      return Differentiation::SD::Expression(static_cast<T>(t));
    }

    template <typename T, int dim>
    static Differentiation::SD::Expression
    value(const Tensor<0, dim, std::complex<T>> &t)
    {
      return Differentiation::SD::Expression(static_cast<std::complex<T>>(t));
    }
  };
} // namespace internal


namespace numbers
{
  template <>
  inline bool
  value_is_zero(const Differentiation::SD::Expression &value)
  {
    if (SymEngine::is_a_Number(value.get_value()))
      {
        const SymEngine::RCP<const SymEngine::Number> number_rcp =
          SymEngine::rcp_static_cast<const SymEngine::Number>(value.get_RCP());
        return number_rcp->is_zero();
      }

    return false;
  }

  template <>
  inline bool
  values_are_equal(const Differentiation::SD::Expression &value_1,
                   const Differentiation::SD::Expression &value_2)
  {
    return (value_1.get_value().__cmp__(value_2.get_value()) == 0);
  }

  template <>
  inline bool
  value_is_less_than(const Differentiation::SD::Expression &value_1,
                     const Differentiation::SD::Expression &value_2)
  {
    return (value_1.get_value().__cmp__(value_2.get_value()) == -1);
  }
} // namespace numbers


#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_number_types_h
