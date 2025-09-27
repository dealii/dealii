// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_parser_h
#define dealii_function_parser_h


#include <deal.II/base/config.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <map>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename>
class Vector;
#endif

/**
 * This class implements a function object that gets its value by parsing a
 * string describing this function. It is a wrapper class for the muparser
 * library (see https://beltoforion.de/en/muparser/). This class lets you
 * evaluate strings such as "sqrt(1-x^2+y^2)" for given values of 'x' and 'y'.
 * Please refer to the muparser documentation for more information.  This class
 * is used in the step-33 and step-36 tutorial programs (the latter being much
 * simpler to understand).
 *
 * In addition to the built-in functions of muparser, namely
 * @code
 * sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh,
 * atan2, log2, log10, log, ln, exp, sqrt, sign, rint, abs, min, max, sum, avg
 * @endcode
 * this class also supports the following operations:
 * - <code>if(condition, then-value, else-value)</code>
 * - <code>|</code> and <code>&</code> (logic or and and)
 * - <code>int(x)</code>, <code>ceil(x)</code>, <code>floor(x)</code> (rounding)
 * - <code>cot(x)</code>, <code>csc(x)</code>, <code>sec(x)</code>
 * - <code>pow(x,n)</code>, <code>log(x)</code>
 * - <code>erfc(x)</code>
 * - <code>rand()</code>, <code>rand_seed(seed)</code>
 *
 * @note This class implements the list of functions just mentioned as
 *   user-defined functions by extending muparser. This means, in particular,
 *   that the `if(condition, then-value, else-value)` syntax evaluates all
 *   three arguments before determining whether the condition is true, and
 *   then discarding either the "then" or the "else" expressions. In almost
 *   all situations, this is not a problem except if the evaluation of
 *   one of the expressions throws a floating point exception in cases
 *   where it will later be discarded. (Assuming floating point exceptions
 *   are switched on, as is the default for deal.II in debug mode on most
 *   systems.) An example would be the expression `if(x>0, sqrt(x), 0)`
 *   which is mathematically well defined, but on systems where this is
 *   enabled will abort the program with a floating point exception when
 *   evaluated with a negative `x`. This is because the square root of
 *   `x` is computed before the `if` statement's condition is considered
 *   to determine whether the result should be the second or third
 *   argument. If this kind of behavior is a problem, you can resort to
 *   the muparser built-in syntax `(condition ? then-value : else-value)`,
 *   using the ternary syntax familiar to C++ programmers. If this
 *   syntax is used, muparser uses lazy evaluation in which only one of the
 *   branches is evaluated, depending on whether the `condition` is
 *   true or false.
 *
 * The following examples shows how to use this class:
 * @code
 * // set up problem:
 * std::string variables = "x,y";
 * std::string expression = "cos(x) + sqrt(y)";
 * std::map<std::string, double> constants;
 *
 * // FunctionParser with 2 variables and 1 component:
 * FunctionParser<2> fp(1);
 * fp.initialize(variables,
 *               expression,
 *               constants);
 *
 * // Point at which we want to evaluate the function
 * Point<2> point(0.0, 4.0);
 *
 * // evaluate the expression at 'point':
 * double result = fp.value(point);
 *
 * deallog << "Function '" << expression << "'"
 *         << " @ " << point
 *         << " is " << result << std::endl;
 * @endcode
 * The second example is a bit more complex:
 * @code
 * // Define some constants that will be used by the function parser
 * std::map<std::string, double> constants;
 * constants["pi"] = numbers::PI;
 *
 * // Define the variables that will be used inside the expressions
 * std::string variables = "x,y,z";
 *
 * // Define the expressions of the individual components of a
 * // vector valued function with two components:
 * std::vector<std::string> expressions(2);
 * expressions[0] = "sin(2*pi*x)+sinh(pi*z)";
 * expressions[1] = "sin(2*pi*y)*exp(x^2)";
 *
 * // function parser with 3 variables and 2 components
 * FunctionParser<3> vector_function(2);
 *
 * // And populate it with the newly created objects.
 * vector_function.initialize(variables,
 *                            expressions,
 *                            constants);
 *
 * // Point at which we want to evaluate the function
 * Point<3> point(0.0, 1.0, 1.0);
 *
 * // This Vector will store the result
 * Vector<double> result(2);
 *
 * // Fill 'result' by evaluating the function
 * vector_function.vector_value(point, result);
 *
 * // We can also only evaluate the 2nd component:
 * const double c = vector_function.value(point, 1);
 *
 * // Output the evaluated function
 * deallog << "Function '" << expressions[0] << ',' << expressions[1] << "'"
 *         << " at " << point
 *         << " is " << result << std::endl;
 * @endcode
 *
 * This class overloads the virtual methods value() and vector_value() of the
 * Function base class with the byte compiled versions of the expressions
 * given to the initialize() methods. Note that the class will not work unless
 * you first call the initialize() method that accepts the text description of
 * the function as an argument (among other things).
 *
 * The syntax to describe a function follows usual programming practice, and
 * is explained in detail at the homepage of the underlying muparser library
 * at https://beltoforion.de/en/muparser/.
 *
 * If you would like to check that muparser is parsing your functions correctly,
 * and to evaluate the functions at given parameter values, you may consider
 * running your expressions through pymuparser
 * (https://github.com/bobmyhill/pymuparser),
 * which can be installed using pip (python -m pip install pymuparser).
 * This module also allows users to define functions not included in MuParser,
 * such as the extended library provided by deal.II.
 *
 * For a wrapper of the FunctionParser class that supports ParameterHandler,
 * see Functions::ParsedFunction.
 *
 * Vector-valued functions can either be declared using strings where the
 * function components are separated by semicolons, or using a vector of
 * strings each defining one vector component.
 *
 * An example of time dependent scalar function is the following:
 * @code
 *    // Empty constants object
 *    std::map<std::string,double> constants;
 *
 *    // Variables that will be used inside the expressions
 *    std::string variables = "x,y,t";
 *
 *    // Define the expression of the scalar time dependent function.
 *    std::string expression = "exp(y*x)*exp(-t)";
 *
 *    // Generate an empty scalar function
 *    FunctionParser<2> function;
 *
 *    // And populate it with the newly created objects.
 *    function.initialize(variables,
 *                        expression,
 *                        constants,
 * // Treat the last variable ("t") as time.
 *                        true);
 * @endcode
 *
 * The following is another example of how to instantiate a vector valued
 * function by using a single string:
 * @code
 *    // Empty constants object
 *    std::map<std::string,double> constants;
 *
 *    // Variables that will be used inside the expressions
 *    std::string variables = "x,y";
 *
 *    // Define the expression of the vector valued  function.
 *    std::string expression = "cos(2*pi*x)*y^2; sin(2*pi*x)*exp(y)";
 *
 *    // Generate an empty vector valued function
 *    FunctionParser<2> function(2);
 *
 *    // And populate it with the newly created objects.
 *    function.initialize(variables,
 *                        expression,
 *                        constants);
 * @endcode
 *
 * @note The difference between this class and the SymbolicFunction class is
 * that the SymbolicFunction class allows to compute first and second order
 * derivatives (in a symbolic way), while this class computes first order
 * derivatives only, using finite differences. For complicated expressions,
 * this class is generally faster than SymbolicFunction.
 *
 * @ingroup functions
 */
template <int dim>
class FunctionParser
  : public AutoDerivativeFunction<dim>,
    protected internal::FunctionParser::ParserImplementation<dim, double>
{
public:
  /**
   * Constructor. Its arguments are the same of the base class Function, with
   * the additional parameter @p h, used for the computation of gradients
   * using finite differences. This object needs to be initialized with the
   * initialize() method before you can use it. If an attempt to use this
   * function is made before the initialize() method has been called, then an
   * exception is thrown.
   */
  FunctionParser(const unsigned int n_components = 1,
                 const double       initial_time = 0.0,
                 const double       h            = 1e-8);

  /**
   * Constructor for parsed functions. Takes directly a semi-colon separated
   * list of expressions (one for each component of the function), an optional
   * comma-separated list of constants, variable names and step size for the
   * computation of first order derivatives by finite differences.
   */
  FunctionParser(const std::string &expression,
                 const std::string &constants      = "",
                 const std::string &variable_names = default_variable_names() +
                                                     ",t",
                 const double h = 1e-8);

  /**
   * Copy constructor. Objects of this type can not be copied, and
   * consequently this constructor is deleted.
   */
  FunctionParser(const FunctionParser &) = delete;

  /**
   * Move constructor. Objects of this type can not be moved, and
   * consequently this constructor is deleted.
   */
  FunctionParser(FunctionParser &&) = delete;

  /**
   * Copy operator. Objects of this type can not be copied, and
   * consequently this operator is deleted.
   */
  FunctionParser &
  operator=(const FunctionParser &) = delete;

  /**
   * Move operator. Objects of this type can not be moved, and
   * consequently this operator is deleted.
   */
  FunctionParser &
  operator=(FunctionParser &&) = delete;

  /**
   * Type for the constant map. Used by the initialize() method.
   */
  using ConstMap = std::map<std::string, double>;

  /**
   * Initialize the object by setting the actual parsed functions.
   *
   * @param[in] vars a string with the variables, separated by commas, that will
   * be used by the expressions to be evaluated. Note that the variables can
   * have any name (of course different from the function names defined above!),
   * but the order IS important. The first variable will correspond to the first
   * component of the point in which the function is evaluated, the second
   * variable to the second component and so forth. If this function is also
   * time dependent, then it is necessary to specify it by setting the
   * <code>time_dependent</code> parameter to true. An exception is thrown if
   * the number of variables specified here is different from dim (if this
   * function is not time-dependent) or from dim+1 (if it is time-dependent).
   *
   * @param[in] expressions a list of strings containing the expressions that
   * will be byte compiled by the internal parser (muParser). Note that
   * the size of this vector must match exactly the number of components of
   * the FunctionParser, as declared in the constructor. If this is not the
   * case, an exception is thrown.
   *
   * @param[in] constants a map of constants used to pass any necessary constant
   * that we want to specify in our expressions (in the example above the
   * number pi). An expression is valid if and only if it contains only
   * defined variables and defined constants (other than the functions
   * specified above). If a constant is given whose name is not valid (eg:
   * <code>constants["sin"] = 1.5;</code>) an exception is thrown.
   *
   * @param[in] time_dependent If this is a time dependent function, then the
   * last variable declared in @p vars is assumed to be the time variable, and
   * FunctionTime::get_time() is used to initialize it when evaluating the
   * function. Naturally the number of variables parsed by initialize() in
   * this case is dim+1. The value of this parameter defaults to false, i.e.,
   * do not consider time.
   */
  virtual void
  initialize(const std::string              &vars,
             const std::vector<std::string> &expressions,
             const ConstMap                 &constants,
             const bool                      time_dependent = false) override;

  /**
   * Initialize the function. Same as above, but accepts a string rather than
   * a vector of strings. If this is a vector valued function, its components
   * are expected to be separated by a semicolon. An exception is thrown if
   * this method is called and the number of components successfully parsed
   * does not match the number of components of the base function.
   */
  void
  initialize(const std::string &vars,
             const std::string &expression,
             const ConstMap    &constants,
             const bool         time_dependent = false);

  /**
   * A function that returns default names for variables, to be used in the
   * first argument of the initialize() functions: it returns "x" in 1d, "x,y"
   * in 2d, and "x,y,z" in 3d.
   */
  static std::string
  default_variable_names();

  /**
   * Return the value of the function at the given point. Unless there is only
   * one component (i.e., the function is scalar), you should state the
   * component you want to have evaluated; it defaults to zero, i.e., the first
   * component.
   */
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * Return an array of function expressions (one per component), used to
   * initialize this function.
   */
  const std::vector<std::string> &
  get_expressions() const;

  /**
   * @addtogroup Exceptions
   * @{
   */
  DeclException2(ExcParseError,
                 int,
                 std::string,
                 << "Parsing Error at Column " << arg1
                 << ". The parser said: " << arg2);

  DeclException2(ExcInvalidExpressionSize,
                 int,
                 int,
                 << "The number of components (" << arg1
                 << ") is not equal to the number of expressions (" << arg2
                 << ").");

  /** @} */
};


template <int dim>
std::string
FunctionParser<dim>::default_variable_names()
{
  switch (dim)
    {
      case 1:
        return "x";
      case 2:
        return "x,y";
      case 3:
        return "x,y,z";
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return "";
}



DEAL_II_NAMESPACE_CLOSE

#endif
