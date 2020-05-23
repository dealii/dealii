// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_tensor_function_parser_h
#define dealii_tensor_function_parser_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/thread_local_storage.h>

#include <map>
#include <memory>
#include <vector>

namespace mu
{
  class Parser;
}

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename>
class Vector;
#endif

/**
 * This class implements a tensor function object that gets its value by parsing
 * a string describing this function. It is a wrapper class for the muparser
 * library (see http://muparser.beltoforion.de/). This class is essentially an
 * extension of the FunctionParser class to read in a TensorFunction. The class
 * reads in an expression of length dim<sup>rank</sup> (separated by a
 * semicolon) where the components of the tensor function are filled according
 * to the C++ convention (fastest index is the most right one).
 *
 * @note In contrast to the FunctionParser class the TensorFunctionParser class does not support
 * 	 automatic differentiation.
 *
 * A minimal example for the usage of the class would be:
 * @code
 * // set up time dependent tensor function:
 * const std::string variables = "x,y,t";
 * const std::string expression =
 *       "exp(-t)*cos(x+y);-sin(pi*x*y-t);sin(pi*x*y-t);exp(t)*cos(x+y)";
 * std::map<std::string,double> constants;
 * constants["pi"] = numbers::PI;
 *
 * // TensorFunctionParser with 2+1 variables (space + time) in 2D of rank 2.
 * // It is necessary to tell the parser that there is an additional variable
 * // to be taken into account (t).
 * TensorFunctionParser<2,2> tfp;
 * tfp.initialize(variables,
 *               expression,
 *               constants,
 *               true); // flag for time dependence
 *
 * // Point at which we want to evaluate the function
 * Point<2> point(0.0, 1.0);
 *
 * // evaluate the expression at 'point':
 * double result = tfp.value(point);
 *
 * deallog << "Function '" << expression << "'"
 *         << " @ " << point
 *         << " is: "
 *         << std::endl
 *         << result[0][0] << " " << result[0][1] << std::endl
 *         << result[1][0] << " " << result[1][1]
 *         << std::endl;
 * @endcode
 *
 * See also the documentation of the FunctionParser class.
 *
 * This class overloads the virtual method value() and value_list() of the
 * TensorFunction base class with the byte compiled versions of the expressions
 * given to the initialize() methods. Note that the class will not work unless
 * you first call the initialize() method that accepts the text description of
 * the function as an argument (among other things).
 *
 * The syntax to describe a function follows usual programming practice, and
 * is explained in detail at the homepage of the underlying muparser library
 * at http://muparser.beltoforion.de/ .
 *
 *
 * Vector-valued functions can either be declared using strings where the
 * function components are separated by semicolons, or using a vector of
 * strings each defining one vector component.
 *
 *
 * @ingroup functions
 * @author Konrad Simon, 2019
 */
template <int rank, int dim, typename Number = double>
class TensorFunctionParser : public TensorFunction<rank, dim, Number>
{
public:
  /**
   * Standard constructor. Only set initial time. This object needs to be
   * initialized with the initialize() method before you can use it. If an
   * attempt to use this function is made before the initialize() method has
   * been called, then an exception is thrown.
   */
  TensorFunctionParser(const double initial_time = 0.0);

  /**
   * Constructor for parsed functions. This object needs to be initialized
   * with the initialize() method before you can use it. If an attempt to
   * use this function is made before the initialize() method has been called,
   * then an exception is thrown.
   * Takes a semicolon separated list of expressions (one for each component
   * of the tensor function), an optional comma-separated list of constants.
   */
  TensorFunctionParser(
    const std::string &expression,
    const std::string &constants      = "",
    const std::string &variable_names = default_variable_names() + ",t");

  /**
   * Copy constructor. Objects of this type can not be copied, and
   * consequently this constructor is deleted.
   */
  TensorFunctionParser(const TensorFunctionParser &) = delete;

  /**
   * Move constructor. Objects of this type can not be moved, and
   * consequently this constructor is deleted.
   */
  TensorFunctionParser(TensorFunctionParser &&) = delete;

  /**
   * Destructor.
   */
  virtual ~TensorFunctionParser() override;

  /**
   * Copy operator. Objects of this type can not be copied, and
   * consequently this operator is deleted.
   */
  TensorFunctionParser &
  operator=(const TensorFunctionParser &) = delete;

  /**
   * Move operator. Objects of this type can not be moved, and
   * consequently this operator is deleted.
   */
  TensorFunctionParser &
  operator=(TensorFunctionParser &&) = delete;

  /**
   * Type for the constant map. Used by the initialize() method.
   */
  using ConstMap = std::map<std::string, double>;


  /**
   * Initialize the tensor function.  This method accepts the following
   * parameters:
   *
   * @param[in] vars A string with the variables that will be used by the
   * expressions to be evaluated. Note that the variables can have any name
   * (of course different from the function names defined above!), but the
   * order IS important. The first variable will correspond to the first
   * component of the point in which the function is evaluated, the second
   * variable to the second component and so forth. If this function is also
   * time dependent, then it is necessary to specify it by setting the
   * <code>time_dependent</code> parameter to true.  An exception is thrown if
   * the number of variables specified here is different from dim (if this
   * function is not time-dependent) or from dim+1 (if it is time-dependent).
   *
   * @param[in] expressions A vector of strings containing the expressions that
   * will be byte compiled by the internal parser (TensorFunctionParser). Note
   * that the size of this vector must match exactly the number of components of
   * the TensorFunctionParser, as declared in the constructor. If this is not
   * the case, an exception is thrown.
   *
   * @param[in] constants A map of constants used to pass any necessary constant
   * that we want to specify in our expressions (in the example above the
   * number pi). An expression is valid if and only if it contains only
   * defined variables and defined constants (other than the functions
   * specified above). If a constant is given whose name is not valid (eg:
   * <code>constants["sin"] = 1.5;</code>) an exception is thrown.
   *
   * @param[in] time_dependent If this is a time dependent function, then the
   * last variable declared in <b>vars</b> is assumed to be the time variable,
   * and this->get_time() is used to initialize it when evaluating the
   * function. Naturally the number of variables parsed by the initialize()
   * method in this case is dim+1. The value of this parameter defaults to
   * false, i.e. do not consider time.
   */
  void
  initialize(const std::string &             vars,
             const std::vector<std::string> &expressions,
             const ConstMap &                constants,
             const bool                      time_dependent = false);

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
             const ConstMap &   constants,
             const bool         time_dependent = false);

  /**
   * A function that returns default names for variables, to be used in the
   * first argument of the initialize() functions: it returns "x" in 1d, "x,y"
   * in 2d, and "x,y,z" in 3d.
   */
  static std::string
  default_variable_names();

  /**
   * Return the value of the tensor function at the given point.
   */
  virtual Tensor<rank, dim, Number>
  value(const Point<dim> &p) const override;

  /**
   * Return the value of the tensor function at the given point.
   */
  virtual void
  value_list(const std::vector<Point<dim>> &         p,
             std::vector<Tensor<rank, dim, Number>> &values) const override;

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

  //@}

private:
#ifdef DEAL_II_WITH_MUPARSER
  /**
   * Place for the variables for each thread
   */
  mutable Threads::ThreadLocalStorage<std::vector<double>> vars;

  /**
   * The muParser objects for each thread (and one for each component). We are
   * storing a unique_ptr so that we don't need to include the definition of
   * mu::Parser in this header.
   */
  mutable Threads::ThreadLocalStorage<std::vector<std::unique_ptr<mu::Parser>>>
    tfp;

  /**
   * An array to keep track of all the constants, required to initialize tfp in
   * each thread.
   */
  std::map<std::string, double> constants;

  /**
   * An array for the variable names, required to initialize tfp in each
   * thread.
   */
  std::vector<std::string> var_names;

  /**
   * Initialize tfp and vars on the current thread. This function may only be
   * called once per thread. A thread can test whether the function has
   * already been called by testing whether 'tfp.get().size()==0' (not
   * initialized) or >0 (already initialized).
   */
  void
  init_muparser() const;
#endif

  /**
   * An array of function expressions (one per component), required to
   * initialize tfp in each thread.
   */
  std::vector<std::string> expressions;

  /**
   * State of usability. This variable is checked every time the function is
   * called for evaluation. It's set to true in the initialize() methods.
   */
  bool initialized;

  /**
   * Number of variables. If this is also a function of time, then the number
   * of variables is dim+1, otherwise it is dim. In the case that this is a
   * time dependent function, the time is supposed to be the last variable. If
   * #n_vars is not identical to the number of the variables parsed by the
   * initialize() method, then an exception is thrown.
   */
  unsigned int n_vars;

  /**
   * Number of components is equal dim<sup>rank</sup>.
   */
  unsigned int n_components;
};


template <int rank, int dim, typename Number>
std::string
TensorFunctionParser<rank, dim, Number>::default_variable_names()
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
        Assert(false, ExcNotImplemented());
    }
  return "";
}



DEAL_II_NAMESPACE_CLOSE

#endif
