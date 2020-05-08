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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_symbolic_function_h
#define dealii_symbolic_function_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/sd.h>

#include <functional>
#include <iostream>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
namespace Functions
{
  template <int dim, typename RangeNumberType>
  class SymbolicFunction;
}
#endif

namespace Functions
{
#ifdef DEAL_II_WITH_SYMENGINE
  /**
   * A Function class that leverages symbolic differentiation to compute
   * gradients, Laplacians, Hessians, and time derivatives.
   *
   * This class can be used to define functions using methods provided by the
   * Differentiation::SD namespace. In particular, one can define a symbolic
   * evaluation point (the argument of the function), as well as a symbolic
   * expression.
   *
   * The symbolic gradients and the symbolic Hessians are computed at
   * construction time, and when a substitution in the symbolic functions is
   * requested by the user using the method update_user_substitution_map().
   *
   * Whenever one of the evaluation methods is called, a substitution is
   * attempted with the coordinate symbols argument replaced by the evaluation
   * point and the symbolic time replaced by the current time, as returned by
   * the get_time() method. The user has to make sure that at evaluation time
   * argument substitution provides a fully evaluated expression (i.e., no other
   * symbols are contained in the function expression, except numerical values),
   * or an exception will be thrown. Additional symbols can be partially
   * evaluated or substituted by storing them in a user supplied substitution
   * maps, that can be updated by calling update_user_substitution_map() or the
   * set_additional_function_arguments() methods.
   *
   * The simplest use case of this class is given in the following example:
   * @code
   * SymbolicFunction<2> fun("x^2+y; t*x*y");
   * fun.set_time(3.0);
   * Point<2> p(1.0, 2.0);
   *
   * auto a = fun.value(p, / * component * / 0); // a = 3.0
   * auto b = fun.value(p, / * component * / 1); // b = 6.0
   *
   * auto df_dt = fun.time_derivative();
   *
   * auto c = df_dt.value(p, / * component * / 0); // c = 0.0
   * auto d = df_dt.value(p, / * component * / 1); // d = 2.0
   * @endcode
   * where a Function with two components is defined using a string containing
   * their expressions separated by semicolons.
   *
   * A more involved example, that explicitly uses
   * Differentiation::SD::Expression objects, is given by
   * @code
   * using namespace Differentiation::SD;
   * // Create a position Tensor<1,2,Differentiation::SD::Expression>
   * // with symbols "x" and "y", and the symbol "t"
   * const auto x = SymbolicFunction<2>::get_default_coordinate_symbols();
   * const auto t = make_symbol("t");
   *
   * // Use directly x[0] (the symbol "x"), x[1] (the symbol "y"), and t
   * // (the symbol "t").
   * Expression f = std::sin(x[0])*std::cos(x[1])*std::sin(t);
   * // Alternatively, you can achieve the same result parsing a string:
   * // Expression f("sin(x)*cos(y)*sin(t)", true);
   * SymbolicFunction<2> function({f}, x);
   *
   * // Evaluate the function, its gradient, and its Laplacian
   * Point<2> p(1.0, 2.0);
   * auto fp = function.value(p);
   * auto gradfp = function.gradient(p);
   * auto lapfp = function.laplacian(p);
   *
   * // Evaluate the time derivative of the function, its gradient, and its
   * // Laplacian
   * auto time_derivative = function.time_derivative();
   * auto dt_fp = time_derivative.value(p);
   * auto dt_gradfp = time_derivative.gradient(p);
   * auto dt_lapfp = time_derivative.laplacian(p);
   * @endcode
   *
   * Partial substitution is possible (i.e., you can define the function using
   * additional symbols). However, as soon as you evaluate the function, you
   * have to make sure that all extraneous symbols (i.e., those not referring
   * to the spacial @p coordinate_symbols or to the @p time_symbol variable)
   * have been substituted with numerical values, or expressions of the spatial
   * or temporal argument, by calling the update_user_substitution_map() or the
   * set_additional_function_arguments() methods.
   *
   * If your function requires additional arguments to be evaluated, you can
   * specify them by calling the set_additional_function_arguments() method.
   *
   * If you call update_user_substitution_map() and
   * set_additional_function_arguments() with the same argument, the effect on
   * the function evaluation will be the same, however, the internal behaviour
   * and function derivatives will be different. The method
   * update_user_substitution_map() performs the substitution once (the first
   * time it is required), and then stores internally a copy of the resulting
   * expression, together with its derivatives (if required). These are then
   * used in all subsequent evaluations. Calling
   * set_additional_function_arguments() will evaluate the passed
   * substitution map on the fly during evaluation time, *after* all
   * derivatives have been computed.
   *
   * @note The difference between this class and the FunctionParser class is
   * that this class allows to compute first and second order derivatives (in a
   * symbolic way), while the FunctionParser class computes first order
   * derivatives only, using finite differences. For complicated expressions,
   * this class may be slower than the FunctionParser class.
   *
   * @ingroup functions
   * @author Luca Heltai 2019
   */
  template <int dim, typename RangeNumberType = double>
  class SymbolicFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * Constructor.
     *
     * The resulting Function object will have as many components as there
     * are entries in the vector of symbolic expressions @p function.
     *
     * The vector @p function should contain a list of symbolic expression
     * involving the coordinate symbols argument @p coordinate_symbols and
     * possibly the symbolic time argument @p time_symbol. It is possible to
     * define it in terms of other symbols, as long as the optional parameter
     * @p user_substitution_map replaces all symbols except
     * @p coordinate_symbols and @p time_symbol.
     * This is useful if, for example, you want to express formulas in terms of
     * material parameters that you want to name symbolically, rather than
     * through their numeric values when defining the formula, or when you want
     * to express your formula in terms of polar coordinates rather than
     * cartesian ones, and you want the symbolic engine to compute the
     * derivatives for you.
     * You may later update the symbol map contained in @p user_substitution_map
     * by calling update_user_substitution_map().
     *
     * @param function A vector of symbolic expressions of type
     * Differentiation::SD::Expression, representing the components of this
     * Function.
     *
     * @param coordinate_symbols A tensor of symbols representing coordinates,
     * used as input argument in the symbolic expressions contained in the
     * @p function vector. The default @p coordinate_symbols is a
     * Tensor<1,dim,Differentiation::SD::Expression>
     * containing the symbols "x" for `dim` equal to one, "x", "y" for `dim`
     * equal to two, and "x", "y", "z" for `dim` equal to three.
     *
     * @param time_symbol A symbolic variable representing time. It defaults
     * to a symbolic variable named "t".
     *
     * @param user_substitution_map Any other symbol that may be contained in
     * the symbolic function needs to be specified in this map. The map may be
     * empty, and the functions may still contain unevaluated symbols, provided
     * that you call update_user_substitution_map() and provide a replacement of
     * all symbols except @p coordinate_symbols and @p time_symbol before any
     * evaluation occurs.
     */
    SymbolicFunction(
      const std::vector<Differentiation::SD::Expression> &function,
      const Tensor<1, dim, Differentiation::SD::Expression>
        &coordinate_symbols = get_default_coordinate_symbols(),
      const Differentiation::SD::Expression &time_symbol =
        Differentiation::SD::make_symbol("t"),
      const Differentiation::SD::types::substitution_map
        &user_substitution_map = {});

    /**
     * Constructor that takes a single string that describes the function
     * expression as a semicolon separated list of expressions.
     *
     * The symbolic expression can use the default argument and the default
     * symbolic time variable, plus any additional symbols that you may
     * need, provided that you update the user substitution map that substitutes
     * all of them before you try to evaluate the function or its derivatives,
     * by calling update_user_substitution_map(), and that you provide all the
     * additional function arguments of your function using the method
     * set_additional_function_arguments().
     */
    SymbolicFunction(const std::string &expressions);

    /**
     * Store and apply the substitution map @p substitutions to each symbolic
     * component of this Function object.
     *
     * Notice that this method will trigger a recomputation of the
     * gradients, Hessians, and Laplacians of each component.
     */
    void
    update_user_substitution_map(
      const Differentiation::SD::types::substitution_map &substitutions);

    /**
     * Set the additional @p arguments to be substituted in next evaluation
     * step.
     *
     * Notice that the @p arguments are substituted *after* evaluating the
     * @p permanent_user_substitution_map, and after all derivatives are
     * computed. If the additional arguments you pass still depend on the
     * coordinate or time symbols, then evaluation of derivatives will result in
     * a partial derivative evaluation.
     *
     * This method provides a way to evaluate functions that depend on more
     * arguments than simply the coordinates and time. If you want to compute
     * the total derivative w.r.t. to complicated symbolic expressions, you
     * should call update_user_substitution_map() instead.
     */
    void
    set_additional_function_arguments(
      const Differentiation::SD::types::substitution_map &arguments);

    /**
     * Return a tensor of coordinate symbols that can be used to define the
     * expressions of this symbolic function object.
     *
     * The default argument is a Tensor<1,dim,Differentiation::SD::Expression>
     * containing the symbols "x" for `dim` equal to one, "x", "y" for `dim`
     * equal to two, and "x", "y", "z" for `dim` equal to three.
     */
    static Tensor<1, dim, Differentiation::SD::Expression>
    get_default_coordinate_symbols();

    /**
     * Get the actual arguments used for the coordinates in the symbolic
     * function. This object does not include any user-defined arguments.
     */
    const Tensor<1, dim, Differentiation::SD::Expression> &
    get_coordinate_symbols() const;

    /**
     * Get the actual symbolic time in use in this symbolic function.
     */
    const Differentiation::SD::Expression &
    get_time_symbol() const;

    /**
     * Get the actual symbolic expressions used in this symbolic function.
     */
    const std::vector<Differentiation::SD::Expression> &
    get_symbolic_function_expressions() const;

    /**
     * Get the currently stored @p user_substitution_map.
     */
    const Differentiation::SD::types::substitution_map &
    get_user_substitution_map() const;

    /**
     * Return a SymbolicFunction object that represents the time derivative of
     * this function. The spatial argument, the symbolic time, and the currently
     * stored user substitution map are forwarded to the new function.
     */
    SymbolicFunction<dim, RangeNumberType>
    time_derivative() const;

    // documentation inherited from the base class
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    // documentation inherited from the base class
    virtual Tensor<1, dim, RangeNumberType>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    // documentation inherited from the base class
    virtual RangeNumberType
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    // documentation inherited from the base class
    virtual SymmetricTensor<2, dim, RangeNumberType>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;

    /**
     * Print the stored arguments and function expression, as it would be
     * evaluated when calling the method value().
     */
    template <typename StreamType>
    StreamType &
    print(StreamType &out) const;

  private:
    /**
     * Return a substitution map that replaces the argument with the values of
     * @p point, the symbolic time with the value of this->get_time(), and any
     * additional arguments with the substitution map given by
     * @p additional_function_arguments.
     */
    Differentiation::SD::types::substitution_map
    create_evaluation_substitution_map(const Point<dim> &point) const;

    /**
     * Recompute the symbolic value of the function, applying the user
     * substitution map. This may be an expensive computation, and it is called
     * only if necessary.
     */
    void
    update_values() const;

    /**
     * Recompute the symbolic gradient of the function, applying the user
     * substitution map. This may be an expensive computation, and it is called
     * only if necessary.
     */
    void
    update_first_derivatives() const;

    /**
     * Recompute the symbolic Hessian and the symbolic Lapalacian of the
     * function. This may be an expensive computation, and it is called
     * only if necessary.
     */
    void
    update_second_derivatives() const;

    /**
     * The components of this symbolic function, before any subustitution took
     * place. This is immutable, and generated at construction time.
     *
     * Before any evaluation takes place, the @p user_substitution_map is
     * applied to this object, and the result is stored in the internal variable
     * function.
     *
     * During evaluation, the @p symbolic_coordinate, the @p symbolic_time, and
     * any remaining symbols are substituted with the input evaluation point,
     * the current time, and the content of @p additional_function_arguments.
     */
    const std::vector<Differentiation::SD::Expression> user_function;

    /**
     * Store the user substitution map used for expression substitutions. This
     * may be updated with a call to update_user_substitution_map(). Notice that
     * the function may still have unresolved symbols, provided that they are
     * resolved by a call to set_additional_function_arguments().
     */
    Differentiation::SD::types::substitution_map user_substitution_map;

    /**
     * Store a user substitution map used for additional argument
     * substitutions. This will be updated by a call to
     * set_additional_function_arguments().
     */
    Differentiation::SD::types::substitution_map additional_function_arguments;

    /**
     * The actual components of this symbolic function. This is obtained from
     * the @p user_function, after applying the @p user_substitution_map.
     */
    mutable std::vector<Differentiation::SD::Expression> function;

    /**
     * The gradients of each component of this symbolic function. This is
     * obtained by computing the symbolic gradient of the object @p function,
     * that is, after applying the @p user_substitution_map to @p user_function.
     */
    mutable std::vector<Tensor<1, dim, Differentiation::SD::Expression>>
      function_gradient;

    /**
     * The Hessians of each component of this symbolic function. This is
     * obtained by computing the symbolic Hessian of the object @p function,
     * that is, after applying the @p user_substitution_map to @p user_function.
     */
    mutable std::vector<Tensor<2, dim, Differentiation::SD::Expression>>
      function_hessian;

    /**
     * The Laplacians of each component of this symbolic function. This is
     * obtained by computing the symbolic Laplacian of the object @p function,
     * that is, after applying the @p user_substitution_map to @p user_function.
     */
    mutable std::vector<Differentiation::SD::Expression> function_laplacian;

    /**
     * The coordinate symbols argument of the function.
     */
    Tensor<1, dim, Differentiation::SD::Expression> coordinate_symbols;

    /**
     * The symbolic time argument of the function.
     */
    mutable Differentiation::SD::Expression time_symbol;
  };

  /**
   * Allow output using the bitwise left shift operator.
   */
  template <int dim, typename RangeNumberType>
  inline std::ostream &
  operator<<(std::ostream &out, const SymbolicFunction<dim, RangeNumberType> &f)
  {
    return f.print(out);
  }



  // Inline and template functions
  template <int dim, typename RangeNumberType>
  template <typename StreamType>
  StreamType &
  SymbolicFunction<dim, RangeNumberType>::print(StreamType &out) const
  {
    for (unsigned int i = 0; i < dim; ++i)
      out << coordinate_symbols[i] << ", ";
    for (const auto &argument_pair : additional_function_arguments)
      out << argument_pair.first << ", ";
    out << time_symbol << " -> " << user_function[0];
    for (unsigned int i = 1; i < user_function.size(); ++i)
      out << "; " << user_function[i];
    if (!user_substitution_map.empty())
      {
        out << " # ( ";
        std::string sep = "";
        for (const auto &substitution : user_substitution_map)
          {
            out << sep << substitution.first << " = " << substitution.second;
            sep = ", ";
          }
        out << " )";
      }
    return out;
  }
#else
  template <int dim, typename RangeNumberType = double>
  class SymbolicFunction
  {
  public:
    SymbolicFunction()
    {
      AssertThrow(
        false,
        ExcMessage(
          "This class is not available if you did not enable SymEngine "
          "when compiling deal.II."));
    }
  };
#endif
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
