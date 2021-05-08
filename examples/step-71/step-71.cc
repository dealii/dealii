/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Jean-Paul Pelteret, 2021
 */


// We start by including all the necessary deal.II header files and some C++
// related ones.
// This first header will give us access to a data structure that will allow
// us to store arbitrary data within it.
#include <deal.II/algorithms/general_data_storage.h>

// Next come some core classes, including one that provides an implementation
// for time-stepping.
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// Then some headers that define some useful coordinate transformations and
// kinematic relationships that are often found in nonlinear elasticity.
#include <deal.II/physics/transformations.h>
#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

// The following two headers provide all of the functionality that we need
// to perform automatic differentiation, and use the symbolic computer algebra
// system that deal.II can utilize. The headers of all automatic
// differentiation and symbolic differentiation wrapper classes, and any
// ancillary data structures that are required, are all collected inside these
// unifying headers.
#include <deal.II/differentiation/ad.h>
#include <deal.II/differentiation/sd.h>

// Including this header allows us the capability to write output to a
// file stream.
#include <fstream>


// As per usual, the entire tutorial program is defined within its own unique
// namespace.
namespace Step71
{
  using namespace dealii;

  // @sect3{An introductory example: The fundamentals of automatic and symbolic differentiation}

  // Automatic and symbolic differentiation have some magical and mystical
  // qualities. Although their use in a project can be beneficial for a
  // multitude of reasons, the barrier to understanding how to use these
  // frameworks or how they can be leveraged may exceed the patience of
  // the developer that is trying to (reliably) integrate them into their work.
  //
  // Although it is the wish of the author to successfully illustrate how these
  // tools can be integrated into workflows for finite element modelling, it
  // might be best to first take a step back and start right from the basics.
  // So to start off with, we'll first have a look at differentiating a "simple"
  // mathematical function using both frameworks, so that the fundamental
  // operations (both their sequence and function) can be firmly established and
  // understood with minimal complication. In the second part of this tutorial
  // we will put these fundamentals into practice and build on them further.
  //
  // Accompanying the description of the algorithmic steps to use the frameworks
  // will be a simplified view as to what they *might* be doing in the
  // background. This description will be very much one designed to aid
  // understanding, and the reader is encouraged to view the @ref auto_symb_diff
  // module documentation for a far more formal description into how these tools
  // actually work.
  //
  // @sect4{An analytical function}
  namespace SimpleExample
  {
    // In order to convince the reader that these tools are indeed useful in
    // practice, let us choose a function for which it is not too difficult to
    // compute the analytical derivatives by hand. It's just sufficiently
    // complicated to make you think about whether or not you truly want to go
    // through with this exercise, and might also make you question whether you
    // are completely sure that your calculations and implementation for its
    // derivatives are correct. The point, of course, is that differentiation of
    // functions is in a sense relatively formulaic and should be something
    // computers are good at -- if we could build on existing software that
    // understands the rules, we wouldn't have to bother with doing it
    // ourselves.
    //
    // We choose the two variable trigonometric function
    // $f(x,y) = \cos\left(\frac{y}{x}\right)$ for this purpose. Notice that
    // this function is templated on the number type. This is done because we
    // can often (but not always) use special auto-differentiable and symbolic
    // types as drop-in replacements for real or complex valued types, and these
    // will then perform some elementary calculations, such as evaluate a
    // function value along with its derivatives. We will exploit that property
    // and make sure that we need only define our function once, and then it can
    // be re-used in whichever context we wish to perform differential
    // operations on it.
    template <typename NumberType>
    NumberType f(const NumberType &x, const NumberType &y)
    {
      return std::cos(y / x);
    }

    // Rather than revealing this function's derivatives immediately, we'll
    // forward declare functions that return them and defer their definition to
    // later. As implied by the function names, they respectively return
    // the derivatives $\frac{df(x,y)}{dx}$:
    double df_dx(const double x, const double y);

    // $\frac{df(x,y)}{dy}$:
    double df_dy(const double x, const double y);

    // $\frac{d^{2}f(x,y)}{dx^{2}}$:
    double d2f_dx_dx(const double x, const double y);

    // $\frac{d^{2}f(x,y)}{dx dy}$:
    double d2f_dx_dy(const double x, const double y);

    // $\frac{d^{2}f(x,y)}{dy dx}$:
    double d2f_dy_dx(const double x, const double y);

    // and, lastly, $\frac{d^{2}f(x,y)}{dy^{2}}$:
    double d2f_dy_dy(const double x, const double y);


    // @sect4{Computing derivatives using automatic differentiation}

    // To begin, we'll use AD as the tool to automatically
    // compute derivatives for us. We will evaluate the function with the
    // arguments `x` and `y`, and expect the resulting value and all of the
    // derivatives to match to within the given tolerance.
    void
    run_and_verify_ad(const double x, const double y, const double tol = 1e-12)
    {
      // Our function $f(x,y)$ is a scalar-valued function, with arguments that
      // represent the typical input variables that one comes across in
      // algebraic calculations or tensor calculus. For this reason, the
      // Differentiation::AD::ScalarFunction class is the appropriate wrapper
      // class to use to do the computations that we require. (As a point of
      // comparison, if the function arguments represented finite element cell
      // degrees-of-freedom, we'd want to treat them differently.) The spatial
      // dimension of the problem is irrelevant since we have no vector- or
      // tensor-valued arguments to accommodate, so the `dim` template argument
      // is arbitrarily assigned a value of 1. The second template argument
      // stipulates which AD framework will be used (deal.II has support for
      // several external AD frameworks), and what the underlying number type
      // provided by this framework is to be used. This number type
      // influences the maximum order of the differential operation, and the
      // underlying algorithms that are used to compute them. Given its template
      // nature, this choice is a compile-time decision because many (but not
      // all) of the AD libraries exploit compile-time meta-programming to
      // implement these special number types in an efficient manner. The third
      // template parameter states what the result type is; in our case, we're
      // working with `double`s.
      constexpr unsigned int                     dim = 1;
      constexpr Differentiation::AD::NumberTypes ADTypeCode =
        Differentiation::AD::NumberTypes::sacado_dfad_dfad;
      using ADHelper =
        Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;

      // It is necessary that we pre-register with our @p ADHelper class how many
      // arguments (what we will call "independent variables") the function
      // $f(x,y)$ has. Those arguments are `x` and `y`, so obviously there
      // are two of them.
      constexpr unsigned int n_independent_variables = 2;

      // We now have sufficient information to create and initialize an
      // instance of the helper class. We can also get the concrete
      // number type that will be used in all subsequent calculations.
      // This is useful, because we can write everything from here on by
      // referencing this type, and if we ever want to change the framework
      // used, or number type (e.g., if we need more differential operations)
      // then we need only adjust the `ADTypeCode` template parameter.
      ADHelper ad_helper(n_independent_variables);
      using ADNumberType = typename ADHelper::ad_type;

      // The next step is to register the numerical values of the independent
      // variables with the helper class. This is done because the function
      // and its derivatives will be evaluated for exactly these arguments.
      // Since we register them in the order `{x,y}`, the variable `x` will
      // be assigned component number `0`, and `y` will be component `1`
      // -- a detail that will be used in the next few lines.
      ad_helper.register_independent_variables({x, y});

      // We now ask for the helper class to give to us the independent variables
      // with their auto-differentiable representation. These are termed
      // "sensitive variables", because from this point on any operations that
      // we do with the components `independent_variables_ad` are tracked and
      // recorded by the AD framework, and will be considered
      // when we ask for the derivatives of something that they're used to
      // compute. What the helper returns is a `vector` of auto-differentiable
      // numbers, but we can be sure that the zeroth element represents `x`
      // and the first element `y`. Just to make completely sure that there's
      // no ambiguity of what number type these variables are, we suffix all of
      // the auto-differentiable variables with `ad`.
      const std::vector<ADNumberType> independent_variables_ad =
        ad_helper.get_sensitive_variables();
      const ADNumberType &x_ad = independent_variables_ad[0];
      const ADNumberType &y_ad = independent_variables_ad[1];

      // We can immediately pass in our sensitive representation of the
      // independent variables to our templated function that computes
      // $f(x,y)$.
      // This also returns an auto-differentiable number.
      const ADNumberType f_ad = f(x_ad, y_ad);

      // So now the natural question to ask is what we have actually just
      // computed by passing these special `x_ad` and `y_ad` variables to the
      // function `f`, instead of the original `double` variables `x` and `y`?
      // In other words, how is all of this related to the computation of the
      // derivatives that we were wanting to determine? Or, more concisely: What
      // is so special about this returned `ADNumberType` object that gives it
      // the ability to magically return derivatives?
      //
      // In essence, how this *could* be done is the following:
      // This special number can be viewed as a data structure that stores the
      // function value, and the prescribed number of derivatives. For a
      // once-differentiable number expecting two arguments, it might look like
      // this:
      //
      // @code
      // struct ADNumberType
      // {
      //   double value;          // The value of the object
      //   double derivatives[2]; // Array of derivatives of the object with
      //                          // respect to x and y
      // };
      // @endcode
      //
      // For our independent variable `x_ad`, the starting value of `x_ad.value`
      // would simply be its assigned value (i.e., the real value of that this
      // variable represents). The derivative `x_ad.derivatives[0]` would be
      // initialized to `1`, since `x` is the zeroth independent variable and
      // $\frac{d(x)}{dx} = 1$. The derivative `x.derivatives[1]` would be
      // initialized to zero, since the first independent variable is `y` and
      // $\frac{d(x)}{dy} = 0$.
      //
      // For the function derivatives to be meaningful, we must assume that not
      // only is this function differentiable in an analytical sense, but that
      // it is also differentiable at the evaluation point `x,y`.
      // We can exploit both of these assumptions: when we use this number type
      // in mathematical operations, the AD framework *could*
      // overload the operations (e.g., `%operator+()`, `%operator*()` as well
      // as `%sin()`, `%exp()`, etc.) such that the returned result has the
      // expected value. At the same time, it would then compute the derivatives
      // through the knowledge of exactly what function is being overloaded and
      // rigorous application of the chain-rule. So, the `%sin()` function
      // (with its argument `a` itself being a function of the independent
      // variables `x` and `y`) *might* be defined as follows:
      //
      // @code
      // ADNumberType sin(const ADNumberType &a)
      // {
      //   ADNumberType output;
      //
      //   // For the input argument "a", "a.value" is simply its value.
      //   output.value = sin(a.value);
      //
      //   // We know that the derivative of sin(a) is cos(a), but we need
      //   // to also consider the chain rule and that the input argument
      //   // `a` is also differentiable with respect to the original
      //   // independent variables `x` and `y`. So `a.derivatives[0]`
      //   // and `a.derivatives[1]` respectively represent the partial
      //   // derivatives of `a` with respect to its inputs `x` and `y`.
      //   output.derivatives[0] = cos(a.value)*a.derivatives[0];
      //   output.derivatives[1] = cos(a.value)*a.derivatives[1];
      //
      //   return output;
      // };
      // @endcode
      //
      // All of that could of course also be done for second and even higher
      // order derivatives.
      //
      // So it is now clear that with the above representation the
      // `ADNumberType` is carrying around some extra data that represents the
      // various derivatives of differentiable functions with respect to the
      // original (sensitive) independent variables. It should therefore be
      // noted that there is computational overhead associated with using them
      // (as we compute extra functions when doing derivative computations) as
      // well as memory overhead in storing these results. So the prescribed
      // number of levels of differential operations should ideally be kept to a
      // minimum to limit computational cost. We could, for instance, have
      // computed the first derivatives ourself and then have used the
      // Differentiation::AD::VectorFunction helper class to determine the
      // gradient of the collection of dependent functions, which would be the
      // second derivatives of the original scalar function.
      //
      // It is also worth noting that because the chain rule is indiscriminately
      // applied and we only see the beginning and end-points of the calculation
      // `{x,y}` $\rightarrow$ `f(x,y)`, we will only ever be able to query
      // the total derivatives of `f`; the partial derivatives
      // (`a.derivatives[0]` and `a.derivatives[1]` in the above example) are
      // intermediate values and are hidden from us.

      // Okay, since we now at least have some idea as to exactly what `f_ad`
      // represents and what is encoded within it, let's put all of that to
      // some actual use. To gain access to those hidden derivative results,
      // we register the final result with the helper class. After this point,
      // we can no longer change the value of `f_ad` and have those changes
      // reflected in the results returned by the helper class.
      ad_helper.register_dependent_variable(f_ad);

      // The next step is to extract the derivatives (specifically, the function
      // gradient and Hessian). To do so we first create some temporary data
      // structures (with the result type `double`) to store the derivatives
      // (noting that all derivatives are returned at once, and not
      // individually)...
      Vector<double>     Df(ad_helper.n_dependent_variables());
      FullMatrix<double> D2f(ad_helper.n_dependent_variables(),
                             ad_helper.n_independent_variables());

      // ... and we then request that the helper class compute these
      // derivatives, and the function value itself. And that's it. We have
      // everything that we were aiming to get.
      const double computed_f = ad_helper.compute_value();
      ad_helper.compute_gradient(Df);
      ad_helper.compute_hessian(D2f);

      // We can convince ourselves that the AD framework is
      // correct by comparing it to the analytical solution. (Or, if you're
      // like the author, you'll be doing the opposite and will rather verify
      // that your implementation of the analytical solution is correct!)
      AssertThrow(std::abs(f(x, y) - computed_f) < tol,
                  ExcMessage(std::string("Incorrect value computed for f. ") +
                             std::string("Hand-calculated value: ") +
                             Utilities::to_string(f(x, y)) +
                             std::string(" ; ") +
                             std::string("Value computed by AD: ") +
                             Utilities::to_string(computed_f)));

      // Because we know the ordering of the independent variables, we know
      // which component of the gradient relates to which derivative...
      const double computed_df_dx = Df[0];
      const double computed_df_dy = Df[1];

      AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dx)));
      AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dy)));

      // ... and similar for the Hessian.
      const double computed_d2f_dx_dx = D2f[0][0];
      const double computed_d2f_dx_dy = D2f[0][1];
      const double computed_d2f_dy_dx = D2f[1][0];
      const double computed_d2f_dy_dy = D2f[1][1];

      AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dx_dx)));
      AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dx_dy)));
      AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dy_dx)));
      AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dy_dy)));
    }

    // That's pretty great. There wasn't too much work involved in computing
    // second-order derivatives of this trigonometric function.

    // @sect4{Hand-calculated derivatives of the analytical solution}

    // Since we now know how much "implementation effort" it takes to have the
    // AD framework compute those derivatives for us, let's
    // compare that to the same computed by hand and implemented in several
    // stand-alone functions.

    // Here are the two first derivatives of $f(x,y) =
    // \cos\left(\frac{y}{x}\right)$:
    //
    // $\frac{df(x,y)}{dx} = \frac{y}{x^2} \sin\left(\frac{y}{x}\right)$
    double df_dx(const double x, const double y)
    {
      Assert(x != 0.0, ExcDivideByZero());
      return y * std::sin(y / x) / (x * x);
    }

    // $\frac{df(x,y)}{dx} = -\frac{1}{x} \sin\left(\frac{y}{x}\right)$
    double df_dy(const double x, const double y)
    {
      return -std::sin(y / x) / x;
    }

    // And here are the four second derivatives of $f(x,y)$:
    //
    // $\frac{d^{2}f(x,y)}{dx^{2}} = -\frac{y}{x^4} (2x
    // \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
    double d2f_dx_dx(const double x, const double y)
    {
      return -y * (2 * x * std::sin(y / x) + y * std::cos(y / x)) /
             (x * x * x * x);
    }

    // $\frac{d^{2}f(x,y)}{dx dy} = \frac{1}{x^3} (x
    // \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
    double d2f_dx_dy(const double x, const double y)
    {
      return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x);
    }

    // $\frac{d^{2}f(x,y)}{dy dx} = \frac{1}{x^3} (x
    // \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$ (as
    // expected, on the basis of [Schwarz's
    // theorem](https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives))
    double d2f_dy_dx(const double x, const double y)
    {
      return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x);
    }

    // $\frac{d^{2}f(x,y)}{dy^{2}} = -\frac{1}{x^2}
    // \cos\left(\frac{y}{x}\right)$
    double d2f_dy_dy(const double x, const double y)
    {
      return -(std::cos(y / x)) / (x * x);
    }

    // Hmm... there's a lot of places in the above where we could have
    // introduced an error in the above, especially when it comes to applying
    // the chain rule. Although they're no silver bullet, at the very least
    // these AD frameworks can serve as a verification tool to make sure that we
    // haven't made any errors (either by calculation or by implementation) that
    // would negatively affect our results.

    // The point of this example of course is that we might have
    // chosen a relatively simple function $f(x,y)$ for which we can
    // hand-verify that the derivatives the AD framework computed is
    // correct. But the AD framework didn't care that the function was
    // simple: It could have been a much much more convoluted
    // expression, or could have depended on more than two variables,
    // and it would still have been able to compute the derivatives --
    // the only difference would have been that *we* wouldn't have
    // been able to come up with the derivatives any more to verify
    // correctness of the AD framework.



    // @sect4{Computing derivatives using symbolic differentiation}

    // We'll now repeat the same exercise using symbolic differentiation. The
    // term "symbolic differentiation" is a little bit misleading because
    // differentiation is just one tool that the Computer Algebra System (CAS)
    // (i.e., the symbolic framework) provides. Nevertheless, in the context
    // of finite element modeling and applications it is the most common use
    // of a CAS and will therefore be the one that we'll focus on.
    // Once more, we'll supply the argument values `x` and `y` with which to
    // evaluate our function $f(x,y) = \cos\left(\frac{y}{x}\right)$ and its
    // derivatives, and a tolerance with which to test the correctness of the
    // returned results.
    void
    run_and_verify_sd(const double x, const double y, const double tol = 1e-12)
    {
      // The first step that we need to take is to form the symbolic variables
      // that represent the function arguments that we wish to differentiate
      // with respect to. Again, these will be the independent variables for
      // our problem and as such are, in some sense, primitive variables that
      // have no dependencies on any other variable. We create these types of
      // (independent) variables by initializing a symbolic type
      // Differentiation::SD::Expression, which is a wrapper to a set of classes
      // used by the symbolic framework, with a unique identifier. On this
      // occasion it makes sense that this identifier, a `std::string`, be
      // simply `"x"` for the $x$ argument, and likewise `"y"` for the $y$
      // argument to the dependent function. Like before, we'll suffix symbolic
      // variable names with `sd` so that we can clearly see which variables are
      // symbolic (as opposed to numeric) in nature.
      const Differentiation::SD::Expression x_sd("x");
      const Differentiation::SD::Expression y_sd("y");

      // Using the templated function that computes $f(x,y)$, we can pass
      // these independent variables as arguments to the function. The returned
      // result will be another symbolic type that represents the sequence of
      // operations used to compute $\cos\left(\frac{y}{x}\right)$.
      const Differentiation::SD::Expression f_sd = f(x_sd, y_sd);

      // At this point it is legitimate to print out the expression `f_sd`, and
      // if we did so
      // @code
      // std::cout << "f(x,y) = " << f_sd << std::endl;
      // @endcode
      // we would see `f(x,y) = cos(y/x)` printed to the console.
      //
      // You might notice that we've constructed our symbolic function `f_sd`
      // with no context as to how we might want to use it: In contrast to the
      // AD approach shown above, what we were returned from calling
      // `f(x_sd, y_sd)` is not the evaluation of the function `f` at some
      // specific point, but is in fact a symbolic representation of the
      // evaluation at a generic, as yet undetermined, point. This is one of the
      // key points that makes symbolic frameworks (the CAS) different from
      // automatic differentiation frameworks. Each of the variables `x_sd` and
      // `y_sd`, and even the composite dependent function `f_sd`, are in some
      // sense respectively "placeholders" for numerical values and a
      // composition of operations. In fact, the individual components that are
      // used to compose the function are also placeholders. The sequence of
      // operations are encoded into in a tree-like data structure (conceptually
      // similar to an [abstract syntax
      // tree](https://en.wikipedia.org/wiki/Abstract_syntax_tree)).
      //
      // Once we form these data structures we can defer any operations that we
      // might want to do with them until some later time. Each of these
      // placeholders represents something, but we have the opportunity to
      // define or redefine what they represent at any convenient point in time.
      // So for this particular problem it makes sense that we somehow want to
      // associate "x" and "y" with *some* numerical value (with type yet to be
      // determined), but we could conceptually (and if it made sense) assign
      // the ratio "y/x" a value instead of the variables "x" and "y"
      // individually. We could also associate with "x" or "y" some other
      // symbolic function `g(a,b)`. Any of these operations involves
      // manipulating the recorded tree of operations, and substituting the
      // salient nodes on the tree (and that nodes' subtree) with something
      // else. The key word here is "substitution", and indeed there are many
      // functions in the Differentiation::SD namespace that have this word
      // in their names.
      //
      // This capability makes the framework entirely generic.
      // In the context of finite element simulations, the types of operations
      // that we would typically perform with our symbolic types are
      // function composition, differentiation, substitution (partial or
      // complete), and evaluation (i.e., conversion of the symbolic type to its
      // numerical counterpart). But should you need it, a CAS is often capable
      // of more than just this: It could be forming anti-derivatives
      // (integrals) of functions, perform simplifications on the expressions
      // that form a function (e.g., replace $(\sin a)^2 + (\cos a)^2$ by
      // $1$; or, more simply: if the function did an operation like `1+2`, a
      // CAS could replace it by `3`), and so forth: The *expression* that a
      // variable represents is obtained from how the function $f$ is
      // implemented, but a CAS can do with it whatever its functionality
      // happens to be.
      //
      // Specifically, to compute the symbolic representation of the first
      // derivatives of the dependent function with respect to its individual
      // independent variables, we use the
      // Differentiation::SD::Expression::differentiate() function with the
      // independent variable given as its argument. Each call will cause the
      // CAS to go through the tree of operations that compose `f_sd` and
      // differentiate each node of the expression tree with respect to the
      // given symbolic argument.
      const Differentiation::SD::Expression df_dx_sd = f_sd.differentiate(x_sd);
      const Differentiation::SD::Expression df_dy_sd = f_sd.differentiate(y_sd);

      // To compute the symbolic representation of the second derivatives, we
      // simply differentiate the first derivatives with respect to the
      // independent variables. So to compute a higher order derivative, we
      // first need to compute the lower order derivative.
      // (As the return type of the call to `differentiate()` is an expression,
      // we could in principal execute double differentiation directly from the
      // scalar by chaining two calls together. But this is unnecessary in this
      // particular case, since we have the intermediate results at hand.)
      const Differentiation::SD::Expression d2f_dx_dx_sd =
        df_dx_sd.differentiate(x_sd);
      const Differentiation::SD::Expression d2f_dx_dy_sd =
        df_dx_sd.differentiate(y_sd);
      const Differentiation::SD::Expression d2f_dy_dx_sd =
        df_dy_sd.differentiate(x_sd);
      const Differentiation::SD::Expression d2f_dy_dy_sd =
        df_dy_sd.differentiate(y_sd);
      // Printing the expressions for the first and second derivatives, as
      // computed by the CAS, using the statements
      // @code
      // std::cout << "df_dx_sd: " << df_dx_sd << std::endl;
      // std::cout << "df_dy_sd: " << df_dy_sd << std::endl;
      // std::cout << "d2f_dx_dx_sd: " << d2f_dx_dx_sd << std::endl;
      // std::cout << "d2f_dx_dy_sd: " << d2f_dx_dy_sd << std::endl;
      // std::cout << "d2f_dy_dx_sd: " << d2f_dy_dx_sd << std::endl;
      // std::cout << "d2f_dy_dy_sd: " << d2f_dy_dy_sd << std::endl;
      // @endcode
      // renders the following output:
      // @code{.sh}
      // df_dx_sd: y*sin(y/x)/x**2
      // df_dy_sd: -sin(y/x)/x
      // d2f_dx_dx_sd: -y**2*cos(y/x)/x**4 - 2*y*sin(y/x)/x**3
      // d2f_dx_dy_sd: sin(y/x)/x**2 + y*cos(y/x)/x**3
      // d2f_dy_dx_sd: sin(y/x)/x**2 + y*cos(y/x)/x**3
      // d2f_dy_dy_sd: -cos(y/x)/x**2
      // @endcode
      // This compares favorably to the analytical expressions for these
      // derivatives that were presented earlier.

      // Now that we have formed the symbolic expressions for the function and
      // its derivatives, we want to evaluate them for the numeric values for
      // the main function arguments `x` and `y`. To accomplish this, we
      // construct a *substitution map*, which maps the symbolic values to their
      // numerical counterparts.
      const Differentiation::SD::types::substitution_map substitution_map =
        Differentiation::SD::make_substitution_map(
          std::pair<Differentiation::SD::Expression, double>{x_sd, x},
          std::pair<Differentiation::SD::Expression, double>{y_sd, y});

      // The last step in the process is to convert all symbolic variables and
      // operations into numerical values, and produce the numerical result of
      // this operation. To do this we combine the substitution map with the
      // symbolic variable in the step we have already mentioned above:
      // "substitution".
      //
      // Once we pass this substitution map to the CAS, it will
      // substitute each instance of the symbolic variable (or, more generally,
      // sub-expression) with its numerical counterpart and then propagate these
      // results up the operation tree, simplifying each node on the tree if
      // possible. If the tree is reduced to a single value (i.e., we have
      // substituted all of the independent variables with their numerical
      // counterpart) then the evaluation is complete.
      //
      // Due to the strongly-typed nature of C++, we need to instruct the CAS to
      // convert its representation of the result into an intrinsic data type
      // (in this case a `double`). This is the "evaluation" step, and through
      // the template type we define the return type of this process.
      // Conveniently, these two steps can be done at once if we are certain
      // that we've performed a full substitution.
      const double computed_f =
        f_sd.substitute_and_evaluate<double>(substitution_map);

      AssertThrow(std::abs(f(x, y) - computed_f) < tol,
                  ExcMessage(std::string("Incorrect value computed for f. ") +
                             std::string("Hand-calculated value: ") +
                             Utilities::to_string(f(x, y)) +
                             std::string(" ; ") +
                             std::string("Value computed by AD: ") +
                             Utilities::to_string(computed_f)));

      // We can do the same for the first derivatives...
      const double computed_df_dx =
        df_dx_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_df_dy =
        df_dy_sd.substitute_and_evaluate<double>(substitution_map);

      AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dx)));
      AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dy)));

      // ... and the second derivatives.
      // Notice that we can reuse the same substitution map for each of these
      // operations because we wish to evaluate all of these functions for the
      // same values of `x` and `y`. Modifying the values in the substitution
      // map renders the result of same symbolic expression evaluated with
      // different values being assigned to the independent variables.
      // We could also happily have each variable represent a real value in
      // one pass, and a complex value in the next.
      const double computed_d2f_dx_dx =
        d2f_dx_dx_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_d2f_dx_dy =
        d2f_dx_dy_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_d2f_dy_dx =
        d2f_dy_dx_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_d2f_dy_dy =
        d2f_dy_dy_sd.substitute_and_evaluate<double>(substitution_map);

      AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dx_dx)));
      AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dx_dy)));
      AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dy_dx)));
      AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dy_dy)));
    }


    // @sect4{The SimpleExample::run() function}

    // The function used to drive these initial examples is straightforward.
    // We'll arbitrarily choose some values at which to evaluate the function
    // (although knowing that `x = 0` is not permissible), and then pass these
    // values to the functions that use the AD and SD frameworks.
    void run()
    {
      const double x = 1.23;
      const double y = 0.91;

      std::cout << "Simple example using automatic differentiation..."
                << std::endl;
      run_and_verify_ad(x, y);
      std::cout << "... all calculations are correct!" << std::endl;

      std::cout << "Simple example using symbolic differentiation."
                << std::endl;
      run_and_verify_sd(x, y);
      std::cout << "... all calculations are correct!" << std::endl;
    }

  } // namespace SimpleExample


  // @sect3{A more complex example: Using automatic and symbolic differentiation to compute derivatives at continuum points}

  // Now that we've introduced the principles behind automatic and symbolic
  // differentiation, we'll put them into action by formulating two coupled
  // magneto-mechanical constitutive laws: one that is rate-independent, and
  // another that exhibits rate-dependent behavior.
  //
  // As you will recall from the introduction, the material
  // constitutive laws we will consider are far more complicated than
  // the simple example above. This is not just because of the form of
  // the function $\psi_{0}$ that we will consider, but in particular
  // because $\psi_{0}$ doesn't just depend on two scalar variables, but
  // instead on a whole bunch of *tensors*, each with several
  // components. In some cases, these are *symmetric* tensors, for
  // which only a subset of components is in fact independent, and one has
  // to think about what it actually means to compute a derivative
  // such as $\frac{\partial\psi_{0}}{\partial \mathbf{C}}$ where $\mathbf
  // C$ is a symmetric tensor. How all of this will work will,
  // hopefully, become clear below. It will also become clear that
  // doing this by hand is going to be, at the very best, *exceedingly*
  // *tedious* and, at worst, riddled with hard-to-find bugs.
  namespace CoupledConstitutiveLaws
  {
    // @sect4{Constitutive parameters}

    // We start with a description of the various material parameters
    // that appear in the description of the energy function $\psi_{0}$.
    //
    // The ConstitutiveParameters class is used to hold these values.
    // Values for all parameters (both constitutive and rheological) are taken
    // from @cite Pelteret2018a, and are given values that produce a
    // constitutive response that is broadly representative of a real,
    // laboratory-made magneto-active polymer, though the specific values used
    // here are of no consequence to the purpose of this program of course.
    //
    // The first four constitutive parameters respectively represent
    // - the elastic shear modulus $\mu_{e}$,
    // - the elastic shear modulus at magnetic saturation $\mu_{e}^{\infty}$,
    // - the saturation magnetic field strength for the elastic shear
    //   modulus $h_{e}^{\text{sat}}$, and
    // - the Poisson ratio $\nu$.
    class ConstitutiveParameters : public ParameterAcceptor
    {
    public:
      ConstitutiveParameters();

      double mu_e       = 30.0e3;
      double mu_e_inf   = 250.0e3;
      double mu_e_h_sat = 212.2e3;
      double nu_e       = 0.49;

      // The next four, which only pertain to the rate-dependent material, are
      // parameters for
      // - the viscoelastic shear modulus $\mu_{v}$,
      // - the viscoelastic shear modulus at magnetic saturation
      // $\mu_{v}^{\infty}$,
      // - the saturation magnetic field strength for the viscoelastic
      //   shear modulus $h_{v}^{\text{sat}}$, and
      // - the characteristic relaxation time $\tau$.
      double mu_v       = 20.0e3;
      double mu_v_inf   = 35.0e3;
      double mu_v_h_sat = 92.84e3;
      double tau_v      = 0.6;

      // The last parameter is the relative magnetic permeability $\mu_{r}$.
      double mu_r = 6.0;

      bool initialized = false;
    };

    // The parameters are initialized through the ParameterAcceptor
    // framework, which is discussed in detail in step-60.
    ConstitutiveParameters::ConstitutiveParameters()
      : ParameterAcceptor("/Coupled Constitutive Laws/Constitutive Parameters/")
    {
      add_parameter("Elastic shear modulus", mu_e);
      add_parameter("Elastic shear modulus at magnetic saturation", mu_e_inf);
      add_parameter(
        "Saturation magnetic field strength for elastic shear modulus",
        mu_e_h_sat);
      add_parameter("Poisson ratio", nu_e);

      add_parameter("Viscoelastic shear modulus", mu_v);
      add_parameter("Viscoelastic shear modulus at magnetic saturation",
                    mu_v_inf);
      add_parameter(
        "Saturation magnetic field strength for viscoelastic shear modulus",
        mu_v_h_sat);
      add_parameter("Characteristic relaxation time", tau_v);

      add_parameter("Relative magnetic permeability", mu_r);

      parse_parameters_call_back.connect([&]() { initialized = true; });
    }


    // @sect4{Constitutive laws: Base class}

    // Since we'll be formulating two constitutive laws for the same class of
    // materials, it makes sense to define a base class that ensures a unified
    // interface to them.
    //
    // The class declaration starts with the constructor that will
    // accept the set of constitutive parameters that, in conjunction
    // with the material law itself, dictate the material response.
    template <int dim>
    class Coupled_Magnetomechanical_Constitutive_Law_Base
    {
    public:
      Coupled_Magnetomechanical_Constitutive_Law_Base(
        const ConstitutiveParameters &constitutive_parameters);

      // Instead of computing and returning the kinetic variables or their
      // linearization at will, we'll calculate and store these values within a
      // single method. These cached results will then be returned upon request.
      // We'll defer the precise explanation as to why we'd want to do this to
      // a later stage. What is important for now is to see that this function
      // accepts all of the field variables, namely the magnetic field vector
      // $\boldsymbol{\mathbb{H}}$ and right Cauchy-Green deformation tensor
      // $\mathbf{C}$, as well as the time discretizer. These, in addition to
      // the @p constitutive_parameters, are all the fundamental quantities that
      // are required to compute the material response.
      virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
                                        const Tensor<1, dim> &         H,
                                        const DiscreteTime &time) = 0;

      // The next few functions provide the interface to probe the material
      // response due subject to the applied deformation and magnetic loading.
      //
      // Since the class of materials can be expressed in terms of a free energy
      // $\psi_{0}$, we can compute that...
      virtual double get_psi() const = 0;

      // ... as well as the two kinetic quantities:
      // - the magnetic induction vector $\boldsymbol{\mathbb{B}}$, and
      // - the total Piola-Kirchhoff stress tensor $\mathbf{S}^{\text{tot}}$
      virtual Tensor<1, dim> get_B() const = 0;

      virtual SymmetricTensor<2, dim> get_S() const = 0;

      // ... and the linearization of the kinetic quantities, which are:
      // - the magnetostatic tangent tensor $\mathbb{D}$,
      // - the total referential magnetoelastic coupling tensor
      // $\mathfrak{P}^{\text{tot}}$, and
      // - the total referential elastic tangent tensor
      // $\mathcal{H}^{\text{tot}}$.
      virtual SymmetricTensor<2, dim> get_DD() const = 0;

      virtual Tensor<3, dim> get_PP() const = 0;

      virtual SymmetricTensor<4, dim> get_HH() const = 0;

      // We'll also define a method that provides a mechanism for this class
      // instance to do any additional tasks before moving on to the next
      // timestep. Again, the reason for doing this will become clear a little
      // later.
      virtual void update_end_of_timestep(){};

      // In the `protected` part of the class,
      // we store a reference to an instance of the constitutive parameters
      // that govern the material response.
      // For convenience, we also define some functions that return
      // various constitutive parameters (both explicitly defined, as well
      // as calculated).
      //
      // The parameters related to the elastic response of the material are,
      // in order:
      // - the elastic shear modulus,
      // - the elastic shear modulus at saturation magnetic field,
      // - the saturation magnetic field strength for the elastic shear
      //   modulus,
      // - the Poisson ratio,
      // - the Lam&eacute; parameter, and
      // - the bulk modulus.
    protected:
      const ConstitutiveParameters &constitutive_parameters;

      double get_mu_e() const;

      double get_mu_e_inf() const;

      double get_mu_e_h_sat() const;

      double get_nu_e() const;

      double get_lambda_e() const;

      double get_kappa_e() const;

      // The parameters related to the elastic response of the material are,
      // in order:
      // - the viscoelastic shear modulus,
      // - the viscoelastic shear modulus at magnetic saturation,
      // - the saturation magnetic field strength for the viscoelastic
      //   shear modulus, and
      // - the characteristic relaxation time.
      double get_mu_v() const;

      double get_mu_v_inf() const;

      double get_mu_v_h_sat() const;

      double get_tau_v() const;

      // The parameters related to the magnetic response of the material are,
      // in order:
      // - the relative magnetic permeability, and
      // - the magnetic permeability constant $\mu_{0}$ (not really a material
      // constant,
      //   but rather a universal constant that we'll group here for
      //   simplicity).
      //
      // We'll also implement a function that returns the
      // timestep size from the time discretizion.
      double get_mu_r() const;

      constexpr double get_mu_0() const;
      double           get_delta_t(const DiscreteTime &time) const;
    };



    // In the following, let us start by implementing the several
    // relatively trivial member functions of the class just defined:
    template <int dim>
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::
      Coupled_Magnetomechanical_Constitutive_Law_Base(
        const ConstitutiveParameters &constitutive_parameters)
      : constitutive_parameters(constitutive_parameters)
    {
      Assert(get_kappa_e() > 0, ExcInternalError());
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e() const
    {
      return constitutive_parameters.mu_e;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_inf() const
    {
      return constitutive_parameters.mu_e_inf;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_h_sat() const
    {
      return constitutive_parameters.mu_e_h_sat;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_nu_e() const
    {
      return constitutive_parameters.nu_e;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_lambda_e() const
    {
      return 2.0 * get_mu_e() * get_nu_e() / (1.0 - 2.0 * get_nu_e());
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_kappa_e() const
    {
      return (2.0 * get_mu_e() * (1.0 + get_nu_e())) /
             (3.0 * (1.0 - 2.0 * get_nu_e()));
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v() const
    {
      return constitutive_parameters.mu_v;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_inf() const
    {
      return constitutive_parameters.mu_v_inf;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_h_sat() const
    {
      return constitutive_parameters.mu_v_h_sat;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_tau_v() const
    {
      return constitutive_parameters.tau_v;
    }


    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_r() const
    {
      return constitutive_parameters.mu_r;
    }


    template <int dim>
    constexpr double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_0() const
    {
      return 4.0 * numbers::PI * 1e-7;
    }


    template <int dim>
    double Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_delta_t(
      const DiscreteTime &time) const
    {
      return time.get_previous_step_size();
    }


    // @sect4{Magnetoelastic constitutive law (using automatic differentiation)}

    // We'll begin by considering a non-dissipative material, namely one that
    // is governed by a magneto-hyperelastic constitutive law that exhibits
    // stiffening when immersed in a magnetic field. As described in
    // the introduction, the stored energy density function for such a material
    // might be given by
    // @f[
    //   \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // = \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //     \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //     \right]
    // + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    // - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //     \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //     \boldsymbol{\mathbb{H}} \right]
    // @f]
    // with
    // @f[
    //  f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    // = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //     \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //     \boldsymbol{\mathbb{H}}}
    //       {\left(h_{e}^{\text{sat}}\right)^{2}} \right) .
    // @f]
    //
    // Now on to the class that implements this behavior.
    // Since we expect that this class fully describes a single material, we'll
    // mark it as "final" so that the inheritance tree terminated here.
    // At the top of the class, we define the helper type that we will use in
    // the AD computations for our scalar energy density function. Note that we
    // expect it to return values of type `double`. We also have to specify the
    // number of spatial dimensions, `dim`, so that the link between vector,
    // tensor and symmetric tensor fields and the number of components that they
    // contain may be established. The concrete `ADTypeCode` used for the
    // ADHelper class will be provided as a template argument at the point where
    // this class is actually used.
    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    class Magnetoelastic_Constitutive_Law_AD final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
      using ADHelper =
        Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;
      using ADNumberType = typename ADHelper::ad_type;

    public:
      Magnetoelastic_Constitutive_Law_AD(
        const ConstitutiveParameters &constitutive_parameters);

      // Since the public interface to the base class is pure-`virtual`, here
      // we'll declare that this class will override all of these base class
      // methods.
      virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
                                        const Tensor<1, dim> &         H,
                                        const DiscreteTime &) override;

      virtual double get_psi() const override;

      virtual Tensor<1, dim> get_B() const override;

      virtual SymmetricTensor<2, dim> get_S() const override;

      virtual SymmetricTensor<2, dim> get_DD() const override;

      virtual Tensor<3, dim> get_PP() const override;

      virtual SymmetricTensor<4, dim> get_HH() const override;

      // In the `private` part of the class,
      // we need to define some extractors that will help us set independent
      // variables and later get the computed values related to the dependent
      // variables. If this class were to be used in the context of a finite
      // element problem, then each of these extractors is (most likely) related
      // to the gradient of a component of the solution field (in this case,
      // displacement and magnetic scalar potential). As you can probably infer
      // by now, here "C" denotes the right Cauchy-Green tensor and "H" denotes
      // the magnetic field vector.
    private:
      const FEValuesExtractors::Vector             H_components;
      const FEValuesExtractors::SymmetricTensor<2> C_components;

      // This is an instance of the automatic differentiation helper that
      // we'll set up to do all of the differential calculations related to
      // the constitutive law...
      ADHelper ad_helper;

      // ... and the following three member variables will store the output from
      // the
      // @p ad_helper. The @p ad_helper returns the derivatives with respect
      // to all field variables at once, so we'll retain the full gradient
      // vector and Hessian matrix. From that, we'll extract the individual
      // entries that we're actually interested in.
      double             psi;
      Vector<double>     Dpsi;
      FullMatrix<double> D2psi;
    };

    // When setting up the field component extractors, it is completely
    // arbitrary as to how they are ordered. But it is important that the
    // extractors do not have overlapping indices. The total number of
    // components of these extractors defines the number of independent
    // variables that the
    // @p ad_helper needs to track, and with respect to which we'll be taking
    // derivatives. The resulting data structures @p Dpsi and @p D2psi must also
    // be sized accordingly. Once the @p ad_helper is configured (its input
    // argument being the total number of components of $\mathbf{C}$ and
    // $\boldsymbol{\mathbb{H}}$), we can directly interrogate it as to how many
    // independent variables it uses.
    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::
      Magnetoelastic_Constitutive_Law_AD(
        const ConstitutiveParameters &constitutive_parameters)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , H_components(0)
      , C_components(Tensor<1, dim>::n_independent_components)
      , ad_helper(Tensor<1, dim>::n_independent_components +
                  SymmetricTensor<2, dim>::n_independent_components)
      , psi(0.0)
      , Dpsi(ad_helper.n_independent_variables())
      , D2psi(ad_helper.n_independent_variables(),
              ad_helper.n_independent_variables())
    {}

    // As stated before, due to the way that the automatic differentiation
    // libraries
    // work, the @p ad_helper will always returns the derivatives of the energy
    // density function with respect to all field variables simultaneously.
    // For this reason, it does not make sense to compute the derivatives in
    // the functions `get_B()`, `get_S()`, etc. because we'd be doing a lot of
    // extra computations that are then simply discarded. So, the best way to
    // deal with that is to have a single function call that does all of the
    // calculations up-front, and then we extract the stored data as its needed.
    // That's what we'll do in the `update_internal_data()` method. As the
    // material is rate-independent, we can ignore the DiscreteTime argument.
    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    void
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &)
    {
      Assert(determinant(C) > 0, ExcInternalError());

      // Since we reuse the @p ad_helper data structure at each time step,
      // we need to clear it of all stale information before use.
      ad_helper.reset();

      // The next step is to set the values for all field components.
      // These define the "point" around which we'll be computing the function
      // gradients and their linearization.
      // The extractors that we created before provide the association between
      // the fields and the registry within the @p ad_helper -- they'll be used
      // repeatedly to ensure that we have the correct interpretation of which
      // variable corresponds to which component of `H` or `C`.
      ad_helper.register_independent_variable(H, H_components);
      ad_helper.register_independent_variable(C, C_components);

      // Now that we've done the initial setup, we can retrieve the AD
      // counterparts of our fields. These are truly the independent variables
      // for the energy function, and are "sensitive" to the calculations that
      // are performed with them. Notice that the AD number are treated as a
      // special number type, and can be used in many templated classes (in this
      // example, as the scalar type for the Tensor and SymmetricTensor class).
      const Tensor<1, dim, ADNumberType> H_ad =
        ad_helper.get_sensitive_variables(H_components);
      const SymmetricTensor<2, dim, ADNumberType> C_ad =
        ad_helper.get_sensitive_variables(C_components);

      // We can also use them in many functions that are templated on the
      // scalar type. So, for these intermediate values that we require,
      // we can perform tensor operations and some mathematical functions.
      // The resulting type will also be an automatically differentiable
      // number, which encodes the operations performed in these functions.
      const ADNumberType det_F_ad = std::sqrt(determinant(C_ad));
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      AssertThrow(det_F_ad > ADNumberType(0.0),
                  ExcMessage("Volumetric Jacobian must be positive."));

      // Next we'll compute the scaling function that will cause the shear
      // modulus to change (increase) under the influence of a magnetic field...
      const ADNumberType f_mu_e_ad =
        1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
                std::tanh((2.0 * H_ad * H_ad) /
                          (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()));

      // ... and then we can define the material stored energy density function.
      // We'll see later that this example is sufficiently complex to warrant
      // the use of AD to, at the very least, verify an unassisted
      // implementation.
      const ADNumberType psi_ad =
        0.5 * this->get_mu_e() * f_mu_e_ad *
          (trace(C_ad) - dim - 2.0 * std::log(det_F_ad))                 //
        + this->get_lambda_e() * std::log(det_F_ad) * std::log(det_F_ad) //
        - 0.5 * this->get_mu_0() * this->get_mu_r() * det_F_ad *
            (H_ad * C_inv_ad * H_ad); //

      // The stored energy density function is, in fact, the dependent variable
      // for this problem, so as a final step in the  "configuration" phase,
      // we register its definition with the @p ad_helper.
      ad_helper.register_dependent_variable(psi_ad);

      // Finally, we can retrieve the resulting value of the stored energy
      // density function, as well as its gradient and Hessian with respect
      // to the input fields, and cache them.
      psi = ad_helper.compute_value();
      ad_helper.compute_gradient(Dpsi);
      ad_helper.compute_hessian(D2psi);
    }

    // The following few functions then allow for querying the so-stored value
    // of $\psi_{0}$, and to extract the desired components of the gradient
    // vector and Hessian matrix. We again make use of the extractors to express
    // which parts of the total gradient vector and Hessian matrix we wish to
    // retrieve. They only return the derivatives of the energy function, so
    // for our definitions of the kinetic variables and their linearization a
    // few more manipulations are required to form the desired result.
    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    double Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_psi() const
    {
      return psi;
    }


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    Tensor<1, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_B() const
    {
      const Tensor<1, dim> dpsi_dH =
        ad_helper.extract_gradient_component(Dpsi, H_components);
      return -dpsi_dH;
    }


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    SymmetricTensor<2, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_S() const
    {
      const SymmetricTensor<2, dim> dpsi_dC =
        ad_helper.extract_gradient_component(Dpsi, C_components);
      return 2.0 * dpsi_dC;
    }


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    SymmetricTensor<2, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_DD() const
    {
      const Tensor<2, dim> dpsi_dH_dH =
        ad_helper.extract_hessian_component(D2psi, H_components, H_components);
      return -symmetrize(dpsi_dH_dH);
    }

    // Note that for coupled terms the order of the extractor
    // arguments is especially important, as it dictates the order in which
    // the directional derivatives are taken. So, if we'd reversed the order
    // of the extractors in the call to `extract_hessian_component()` then we'd
    // actually have been retrieving part of $\left[ \mathfrak{P}^{\text{tot}}
    // \right]^{T}$.
    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    Tensor<3, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_PP() const
    {
      const Tensor<3, dim> dpsi_dC_dH =
        ad_helper.extract_hessian_component(D2psi, C_components, H_components);
      return -2.0 * dpsi_dC_dH;
    }


    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    SymmetricTensor<4, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_HH() const
    {
      const SymmetricTensor<4, dim> dpsi_dC_dC =
        ad_helper.extract_hessian_component(D2psi, C_components, C_components);
      return 4.0 * dpsi_dC_dC;
    }


    // @sect4{Magneto-viscoelastic constitutive law (using symbolic algebra and differentiation)}

    // The second material law that we'll consider will be one that represents
    // a magneto-viscoelastic material with a single dissipative mechanism.
    // We'll consider the free energy density function for such a material to
    // be defined as
    // @f{align*}{
    //   \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
    //   \right)
    // &= \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // + \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)
    // \\ \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // &= \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
    // \right)
    //     \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //     \right]
    // + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    // - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //     \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //     \boldsymbol{\mathbb{H}} \right]
    // \\ \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)
    // &= \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    // \right)
    //     \left[ \mathbf{C}_{v} : \left[
    //       \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //       \mathbf{C} \right] - d - \ln\left(
    //       \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    // @f}
    // with
    // @f[
    //   f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
    // = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //     \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //     \boldsymbol{\mathbb{H}}}
    //       {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    // @f]
    // @f[
    //   f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
    // = 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
    //     \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //     \boldsymbol{\mathbb{H}}}
    //       {\left(h_{v}^{\text{sat}}\right)^{2}} \right),
    // @f]
    // in conjunction with the evolution law for the internal viscous variable
    // @f[
    // \mathbf{C}_{v}^{(t)}
    // = \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
    //     \mathbf{C}_{v}^{(t-1)}
    //   + \frac{\Delta t}{\tau_{v}}
    //     \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //     \mathbf{C} \right]^{-1}
    //   \right]
    // @f]
    // that was discretized using a first-order backward difference
    // approximation.
    //
    // Again, let us see how this is implemented in a concrete class. Instead of
    // the AD framework used in the previous class, we will now utilize the SD
    // approach. To support this,
    // the class constructor accepts not only the @p constitutive_parameters,
    // but also two additional variables that will be used to initialize
    // a Differentiation::SD::BatchOptimizer. We'll give more context to this
    // later.
    template <int dim>
    class Magnetoviscoelastic_Constitutive_Law_SD final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
    public:
      Magnetoviscoelastic_Constitutive_Law_SD(
        const ConstitutiveParameters &               constitutive_parameters,
        const Differentiation::SD::OptimizerType     optimizer_type,
        const Differentiation::SD::OptimizationFlags optimization_flags);

      // Like for the automatic differentiation helper, the
      // Differentiation::SD::BatchOptimizer will return a collection of
      // results all at once. So, in order to do that just once, we'll utilize
      // a similar approach to before and do all of the expensive calculations
      // within the `update_internal_data()` function, and cache the results
      // for layer extraction.
      virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
                                        const Tensor<1, dim> &         H,
                                        const DiscreteTime &time) override;

      virtual double get_psi() const override;

      virtual Tensor<1, dim> get_B() const override;

      virtual SymmetricTensor<2, dim> get_S() const override;

      virtual SymmetricTensor<2, dim> get_DD() const override;

      virtual Tensor<3, dim> get_PP() const override;

      virtual SymmetricTensor<4, dim> get_HH() const override;

      // Since we're dealing with a rate dependent material, we'll have to
      // update the history variable at the appropriate time. That will be the
      // purpose of this function.
      virtual void update_end_of_timestep() override;

      // In the `private` part of the class, we will want to
      // keep track of the internal viscous deformation, so the following
      // two (real-valued, non-symbolic) member variables respectively hold
      // - the value of internal variable time step (and, if embedded within a
      //   nonlinear solver framework, Newton step), and
      // - the value of internal variable at the previous timestep.
      //
      // (We've labeled these variables "Q" so that they're easy to identify;
      // in a sea of calculations it is not necessarily easy to distinguish
      // `Cv` or `C_v` from `C`.)
    private:
      SymmetricTensor<2, dim> Q_t;
      SymmetricTensor<2, dim> Q_t1;

      // As we'll be using symbolic types, we'll need to define some symbolic
      // variables to use with the framework. (They are all suffixed with "SD"
      // to make it easy to distinguish the symbolic types or expressions from
      // real-valued types or scalars.) This can be done once up front
      // (potentially even as `static` variables) to minimize the overhead
      // associated with creating these variables. For the ultimate in generic
      // programming, we can even describe the constitutive parameters
      // symbolically, *potentially* allowing a single class instance to be
      // reused with different inputs for these values too.
      //
      // These are the symbolic scalars that represent the elastic, viscous, and
      // magnetic material parameters
      // (defined mostly in the same order as they appear in the @p ConstitutiveParameters
      // class). We also store a symbolic expression, @p delta_t_sd, that represents the
      // time step size):
      const Differentiation::SD::Expression mu_e_sd;
      const Differentiation::SD::Expression mu_e_inf_sd;
      const Differentiation::SD::Expression mu_e_h_sat_sd;
      const Differentiation::SD::Expression lambda_e_sd;
      const Differentiation::SD::Expression mu_v_sd;
      const Differentiation::SD::Expression mu_v_inf_sd;
      const Differentiation::SD::Expression mu_v_h_sat_sd;
      const Differentiation::SD::Expression tau_v_sd;
      const Differentiation::SD::Expression delta_t_sd;
      const Differentiation::SD::Expression mu_r_sd;

      // Next we define some tensorial symbolic variables that represent the
      // independent field variables, upon which the energy density function
      // is parameterized:
      const Tensor<1, dim, Differentiation::SD::Expression>          H_sd;
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_sd;

      // And similarly we have the symbolic representation of the internal
      // viscous variables (both its current value and its value at the
      // previous timestep):
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t_sd;
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t1_sd;

      // We should also store the definitions of the dependent expressions:
      // Although we'll only compute them once, we require them to retrieve
      // data from the @p optimizer that is declared below.
      // Furthermore, when serializing a material class like this one (not done
      // as a part of this tutorial) we'd either need to serialize these
      // expressions as well or we'd need to reconstruct them upon reloading.
      Differentiation::SD::Expression                          psi_sd;
      Tensor<1, dim, Differentiation::SD::Expression>          B_sd;
      SymmetricTensor<2, dim, Differentiation::SD::Expression> S_sd;
      SymmetricTensor<2, dim, Differentiation::SD::Expression> BB_sd;
      Tensor<3, dim, Differentiation::SD::Expression>          PP_sd;
      SymmetricTensor<4, dim, Differentiation::SD::Expression> HH_sd;

      // The next variable is then the optimizer that is used to evaluate the
      // dependent functions. More specifically, it provides the possibility to
      // accelerate the evaluation of the symbolic dependent expressions. This
      // is a vital tool, because the native evaluation of lengthy expressions
      // (using no method of acceleration, but rather direct evaluation directly
      // of the symbolic expressions) can be very slow. The
      // Differentiation::SD::BatchOptimizer class provides a mechanism by which
      // to transform the symbolic expression tree into another code path that,
      // for example, shares intermediate results between the various dependent
      // expressions (meaning that these intermediate values only get calculated
      // once per evaluation) and/or compiling the code using a just-in-time
      // compiler (thereby retrieving near-native performance for the evaluation
      // step).
      //
      // Performing this code transformation is very computationally expensive,
      // so we store the optimizer so that it is done just once per class
      // instance. This also further motivates the decision to make the
      // constitutive parameters themselves symbolic. We could then reuse a
      // single instance
      // of this @p optimizer across several materials (with the same energy
      // function, of course) and potentially multiple continuum points (if
      // embedded within a finite element simulation).
      //
      // As specified by the template parameter, the numerical result will be of
      // type <tt>double</tt>.
      Differentiation::SD::BatchOptimizer<double> optimizer;

      // During the evaluation phase, we must map the symbolic variables to
      // their real-valued counterparts. The next method will provide this
      // functionality.
      //
      // The final method of this class will configure the @p optimizer.
      Differentiation::SD::types::substitution_map
      make_substitution_map(const SymmetricTensor<2, dim> &C,
                            const Tensor<1, dim> &         H,
                            const double                   delta_t) const;

      void initialize_optimizer();
    };

    // As the resting deformation state is one at which the material is
    // considered to be completely relaxed, the internal viscous variables are
    // initialized with the identity tensor, i.e. $\mathbf{C}_{v} = \mathbf{I}$.
    // The various symbolic variables representing the constitutive parameters,
    // time step size, and field and internal variables all get a unique
    // identifier. The optimizer is passed the two parameters that declare which
    // optimization (acceleration) technique should be applied, as well as
    // which additional steps should be taken by the CAS to help improve
    // performance during evaluation.
    template <int dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::
      Magnetoviscoelastic_Constitutive_Law_SD(
        const ConstitutiveParameters &               constitutive_parameters,
        const Differentiation::SD::OptimizerType     optimizer_type,
        const Differentiation::SD::OptimizationFlags optimization_flags)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , Q_t(Physics::Elasticity::StandardTensors<dim>::I)
      , Q_t1(Physics::Elasticity::StandardTensors<dim>::I)
      , mu_e_sd("mu_e")
      , mu_e_inf_sd("mu_e_inf")
      , mu_e_h_sat_sd("mu_e_h_sat")
      , lambda_e_sd("lambda_e")
      , mu_v_sd("mu_v")
      , mu_v_inf_sd("mu_v_inf")
      , mu_v_h_sat_sd("mu_v_h_sat")
      , tau_v_sd("tau_v")
      , delta_t_sd("delta_t")
      , mu_r_sd("mu_r")
      , H_sd(Differentiation::SD::make_vector_of_symbols<dim>("H"))
      , C_sd(Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("C"))
      , Q_t_sd(
          Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t"))
      , Q_t1_sd(
          Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t1"))
      , optimizer(optimizer_type, optimization_flags)
    {
      initialize_optimizer();
    }

    // The substitution map simply pairs all of the following data together:
    // - the constitutive parameters (with values retrieved from the base
    // class),
    // - the time step size (with its value retrieved from the time
    // discretizer),
    // - the field values (with their values being prescribed by an external
    //   function that is calling into this @p Magnetoviscoelastic_Constitutive_Law_SD instance), and
    // - the current and previous internal viscous deformation (with their
    // values
    //   stored within this class instance).
    template <int dim>
    Differentiation::SD::types::substitution_map
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::make_substitution_map(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const double                   delta_t) const
    {
      return Differentiation::SD::make_substitution_map(
        std::make_pair(mu_e_sd, this->get_mu_e()),
        std::make_pair(mu_e_inf_sd, this->get_mu_e_inf()),
        std::make_pair(mu_e_h_sat_sd, this->get_mu_e_h_sat()),
        std::make_pair(lambda_e_sd, this->get_lambda_e()),
        std::make_pair(mu_v_sd, this->get_mu_v()),
        std::make_pair(mu_v_inf_sd, this->get_mu_v_inf()),
        std::make_pair(mu_v_h_sat_sd, this->get_mu_v_h_sat()),
        std::make_pair(tau_v_sd, this->get_tau_v()),
        std::make_pair(delta_t_sd, delta_t),
        std::make_pair(mu_r_sd, this->get_mu_r()),
        std::make_pair(H_sd, H),
        std::make_pair(C_sd, C),
        std::make_pair(Q_t_sd, Q_t),
        std::make_pair(Q_t1_sd, Q_t1));
    }

    // Due to the "natural" use of the symbolic expressions, much of the
    // procedure to configure the @p optimizer looks very similar to that which
    // is used to construct the automatic differentiation helper.
    // Nevertheless, we'll detail these steps again to highlight the differences
    // that underlie the two frameworks.
    //
    // The function starts with expressions that symbolically encode the
    // determinant of the deformation gradient (as expressed in terms of the
    // right Cauchy-Green deformation tensor, our primary field variable), as
    // well as the inverse of $\mathbf{C}$ itself:
    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law_SD<dim>::initialize_optimizer()
    {
      const Differentiation::SD::Expression det_F_sd =
        std::sqrt(determinant(C_sd));
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_inv_sd =
        invert(C_sd);

      // Next is the symbolic representation of the saturation function for
      // the elastic part of the free energy density function, followed by the
      // magnetoelastic contribution to the free energy density function.
      // This all has the same stucture as we'd seen previously.
      const Differentiation::SD::Expression f_mu_e_sd =
        1.0 +
        (mu_e_inf_sd / mu_e_sd - 1.0) *
          std::tanh((2.0 * H_sd * H_sd) / (mu_e_h_sat_sd * mu_e_h_sat_sd));

      const Differentiation::SD::Expression psi_ME_sd =
        0.5 * mu_e_sd * f_mu_e_sd *
          (trace(C_sd) - dim - 2.0 * std::log(det_F_sd)) +
        lambda_e_sd * std::log(det_F_sd) * std::log(det_F_sd) -
        0.5 * this->get_mu_0() * mu_r_sd * det_F_sd * (H_sd * C_inv_sd * H_sd);

      // In addition, we define the magneto-viscoelastic contribution to the
      // free energy density function. The first component required to implement
      // this is a scaling function that will cause the viscous shear modulus to
      // change (increase) under the influence of a magnetic field (see
      // @cite Pelteret2018a, equation 29). Thereafter we can compute the
      // dissipative component of the energy density function; its expression
      // is stated in @cite Pelteret2018a (equation 28), which is a
      // straight-forward extension of an energy density function formulated in
      // @cite Linder2011a (equation 46).
      const Differentiation::SD::Expression f_mu_v_sd =
        1.0 +
        (mu_v_inf_sd / mu_v_sd - 1.0) *
          std::tanh((2.0 * H_sd * H_sd) / (mu_v_h_sat_sd * mu_v_h_sat_sd));

      const Differentiation::SD::Expression psi_MVE_sd =
        0.5 * mu_v_sd * f_mu_v_sd *
        (Q_t_sd * (std::pow(det_F_sd, -2.0 / dim) * C_sd) - dim -
         std::log(determinant(Q_t_sd)));

      // From these building blocks, we can then define the material's total
      // free energy density function:
      psi_sd = psi_ME_sd + psi_MVE_sd;

      // As it stands, to the CAS the variable @p Q_t_sd appears
      // to be independent of @p C_sd. Our tensorial symbolic expression
      // @p Q_t_sd just has an identifier associated with it, and there is
      // nothing that links it to the other tensorial symbolic expression
      // @p C_sd. So any derivatives taken with respect to @p C_sd will ignore
      // this inherent dependence which, as we can see from the evolution law,
      // is in fact
      // $\mathbf{C}_{v} = \mathbf{C}_{v} \left( \mathbf{C}, t \right)$.
      // This means that deriving any function $f = f(\mathbf{C}, \mathbf{Q})$
      // with respect to  $\mathbf{C}$ will return partial derivatives
      // $\frac{\partial f(\mathbf{C}, \mathbf{Q})}{\partial \mathbf{C}}
      // \Big\vert_{\mathbf{Q}}$ as opposed to the total derivative
      // $\frac{d f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))}{d \mathbf{C}} =
      // \frac{\partial f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))}{\partial
      // \mathbf{C}} \Big\vert_{\mathbf{Q}} + \frac{\partial f(\mathbf{C},
      // \mathbf{Q}(\mathbf{C}))}{\partial \mathbf{Q}}
      // \Big\vert_{\mathbf{C}} : \frac{d \mathbf{Q}(\mathbf{C}))}{d
      // \mathbf{C}}$.
      //
      // By contrast, with the current AD libraries the total derivative would
      // always be returned. This implies that the computed kinetic variables
      // would be incorrect for this class of material model, making AD the
      // incorrect tool from which to derive (at the continuum point level) the
      // constitutive law for this dissipative material from an energy density
      // function.
      //
      // It is this specific level of control that characterizes a defining
      // difference difference between the SD and AD frameworks. In a few lines
      // we'll be manipulating the expression for the internal variable
      // @p Q_t_sd such that it produces the correct linearization.

      // But, first, we'll compute the symbolic expressions for the kinetic
      // variables, i.e., the magnetic induction vector and the Piola-Kirchhoff
      // stress tensor. The code that performs the differentiation quite closely
      // mimics the definition stated in the theory.
      B_sd = -Differentiation::SD::differentiate(psi_sd, H_sd);
      S_sd = 2.0 * Differentiation::SD::differentiate(psi_sd, C_sd);

      // Since the next step is to linearize the above, it is the appropriate
      // time to inform the CAS of the explicit dependency of @p Q_t_sd on @p C_sd,
      // i.e., state that $\mathbf{C}_{v} = \mathbf{C}_{v} \left( \mathbf{C}, t
      // \right)$. This means that all future differential operations made with
      // respect
      // to @p C_sd will take into account this dependence (i.e., compute total derivatives).
      // In other words, we will transform some expression such that their
      // intrinsic parameterization changes from $f(\mathbf{C}, \mathbf{Q})$
      // to $f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))$.
      //
      // To do this, we consider the time-discrete evolution law.
      // From that, we have the explicit expression for the internal
      // variable in terms of its history as well as the primary
      // field variable. That is what it described in this expression:
      const SymmetricTensor<2, dim, Differentiation::SD::Expression>
        Q_t_sd_explicit =
          (1.0 / (1.0 + delta_t_sd / tau_v_sd)) *
          (Q_t1_sd +
           (delta_t_sd / tau_v_sd * std::pow(det_F_sd, 2.0 / dim) * C_inv_sd));

      // Next we produce an intermediate substitution map, which will take
      // every instance of @p Q_t_sd (our identifier) found in an expression
      // and replace it with the full expression held in @p Q_t_sd_explicit.
      const Differentiation::SD::types::substitution_map
        substitution_map_explicit = Differentiation::SD::make_substitution_map(
          std::make_pair(Q_t_sd, Q_t_sd_explicit));

      // We can the perform this substitution on the two kinetic variables
      // and immediately differentiate the result that appears after that
      // substitution with the field variables. (If you'd like, this could
      // be split up into two steps with the intermediate results stored in
      // a temporary variable.) Again, if you overlook the "complexity"
      // generated by the substitution, these calls that linearize the kinetic
      // variables and produce the three tangent tensors quite closely resembles
      // what's stated in the theory.
      BB_sd = symmetrize(Differentiation::SD::differentiate(
        Differentiation::SD::substitute(B_sd, substitution_map_explicit),
        H_sd));
      PP_sd = -Differentiation::SD::differentiate(
        Differentiation::SD::substitute(S_sd, substitution_map_explicit), H_sd);
      HH_sd =
        2.0 *
        Differentiation::SD::differentiate(
          Differentiation::SD::substitute(S_sd, substitution_map_explicit),
          C_sd);

      // Now we need to tell the @p optimizer what entries we need to provide
      // numerical values for in order for it to successfully perform its
      // calculations. These essentially act as the input arguments to
      // all dependent functions that the @p optimizer must evaluate.
      // They are, collectively, the independent variables
      // for the problem, the history variables, the time step sie and the
      // constitutive parameters (since we've not hard encoded them in the
      // energy density function).
      //
      // So what we really want is to provide it a collection of
      // symbols, which one could accomplish in this way:
      // @code
      // optimizer.register_symbols(Differentiation::SD::make_symbol_map(
      //   mu_e_sd, mu_e_inf_sd, mu_e_h_sat_sd, lambda_e_sd,
      //   mu_v_sd, mu_v_inf_sd, mu_v_h_sat_sd, tau_v_sd,
      //   delta_t_sd, mu_r_sd,
      //   H_sd, C_sd,
      //   Q_t_sd, Q_t1_sd));
      // @endcode
      // But this is all actually already encoded as the keys of the
      // substitution map. Doing the above would also mean that we
      // need to manage the symbols in two places (here and when constructing
      // the substitution map), which is annoying and a potential source of
      // error if this material class is modified or extended.
      // Since we're not interested in the values at this point,
      // it is alright if the substitution map is filled with invalid data
      // for the values associated with each key entry.
      // So we'll simply create a fake substitution map, and extract the
      // symbols from that. Note that any substitution map passed to the
      // @p optimizer will have to, at the very least, contain entries for
      // these symbols.
      optimizer.register_symbols(
        Differentiation::SD::Utilities::extract_symbols(
          make_substitution_map({}, {}, 0)));

      // We then inform the optimizer of what values we want calculated, which
      // in our situation encompasses all of the dependent variables (namely
      // the energy density function and its various derivatives).
      optimizer.register_functions(psi_sd, B_sd, S_sd, BB_sd, PP_sd, HH_sd);

      // The last step is to finalize the optimizer. With this call it will
      // determine an equivalent code path that will evaluate all of the
      // dependent functions at once, but with less computational
      // cost than when evaluating the symbolic expression directly.
      // Note: This is an expensive call, so we want execute it as few times
      // as possible. We've done it in the constructor of our class, which
      // achieves the goal of being called only once per class instance.
      optimizer.optimize();
    }

    // Since the configuration of the @p optimizer was done up front, there's
    // very little to do each time we want to compute kinetic variables or
    // their linearization (derivatives).
    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &           time)
    {
      // To update the internal history variable, we first need to compute
      // a few fundamental quantities, which we've seen before.
      // We can also ask the time discretizer for the time step size that
      // was used to iterate from the previous time step to the current one.
      const double delta_t = this->get_delta_t(time);

      const double                  det_F = std::sqrt(determinant(C));
      const SymmetricTensor<2, dim> C_inv = invert(C);
      AssertThrow(det_F > 0.0,
                  ExcMessage("Volumetric Jacobian must be positive."));

      // Now we can update the (real valued) internal viscous deformation
      // tensor, as per the definition given by the evolution law in conjunction
      // with the chosen time discretization scheme.
      Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) *
            (Q_t1 + (delta_t / this->get_tau_v()) * std::pow(det_F, 2.0 / dim) *
                      C_inv);

      // Next we pass the optimizer the numeric values that we wish the
      // independent variables, time step size and (implicit to this call),
      // the constitutive parameters to represent.
      const auto substitution_map = make_substitution_map(C, H, delta_t);

      // When making this next call, the call path used to (numerically)
      // evaluate the dependent functions is quicker than dictionary
      // substitution.
      optimizer.substitute(substitution_map);
    }

    // Having called `update_internal_data()`, it is then valid to
    // extract data from the optimizer.
    // When doing the evaluation, we need the exact symbolic expressions of
    // the data to extracted from the optimizer. The implication of this
    // is that we needed to store the symbolic expressions of all dependent
    // variables for the lifetime of the optimizer (naturally, the same
    // is implied for the input variables).
    template <int dim>
    double Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_psi() const
    {
      return optimizer.evaluate(psi_sd);
    }


    template <int dim>
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_B() const
    {
      return optimizer.evaluate(B_sd);
    }


    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_S() const
    {
      return optimizer.evaluate(S_sd);
    }


    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_DD() const
    {
      return optimizer.evaluate(BB_sd);
    }


    template <int dim>
    Tensor<3, dim> Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_PP() const
    {
      return optimizer.evaluate(PP_sd);
    }


    template <int dim>
    SymmetricTensor<4, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_HH() const
    {
      return optimizer.evaluate(HH_sd);
    }

    // When moving forward in time, the "current" state of the internal variable
    // instantaneously defines the state at the "previous" timestep. As such, we
    // record value of history variable for use as the "past value" at the next
    // time step.
    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_end_of_timestep()
    {
      Q_t1 = Q_t;
    };


    // @sect3{A more complex example (continued): Parameters and hand-derived material classes}

    // Now that we've seen how the AD and SD frameworks can make light(er) work
    // of defining these constitutive laws, we'll implement the equivalent
    // classes by hand for the purpose of verification and to do some
    // preliminary benchmarking of the frameworks versus a native
    // implementation.
    //
    // At the expense of the author's sanity, what is documented below
    // (hopefully accurately) are the full definitions for the kinetic variables
    // and their tangents, as well as some intermediate computations. Since the
    // structure and design of the constitutive law classes has been outlined
    // earlier, we'll gloss over it and simply delineate between the various
    // stages of calculations in the `update_internal_data()` method definition.
    // It should be easy enough to link the derivative calculations (with their
    // moderately expressive variable names) to their documented definitions
    // that appear in the class descriptions.
    // We will, however, take the opportunity to present two different paradigms
    // for implementing constitutive law classes. The second will provide more
    // flexibility than the first (thereby making it more easily extensible,
    // in the author's opinion) at the expense of some performance.

    // @sect4{Magnetoelastic constitutive law (hand-derived)}

    // From the stored energy that, as mentioned earlier, is defined as
    // @f[
    //   \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // = \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //     \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //     \right]
    // + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    // - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //     \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //     \boldsymbol{\mathbb{H}} \right]
    // @f]
    // with
    // @f[
    //  f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    // = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //     \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //     \boldsymbol{\mathbb{H}}}
    //       {\left(h_{e}^{\text{sat}}\right)^{2}} \right) ,
    // \\ \text{det}(\mathbf{F}) = \sqrt{\text{det}(\mathbf{C})}
    // @f]
    // for this magnetoelastic material, the first derivatives that correspond
    // to the magnetic induction vector and total Piola-Kirchhoff stress
    // tensor are
    // @f[
    //  \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right)
    // \dealcoloneq - \frac{d \psi_{0}}{d \boldsymbol{\mathbb{H}}}
    // = - \frac{1}{2} \mu_{e} \left[ \text{tr}(\mathbf{C}) - d - 2 \ln
    // (\text{det}(\mathbf{F}))
    //       \right] \frac{d f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}}
    //       \right)}{d \boldsymbol{\mathbb{H}}}
    // + \mu_{0} \mu_{r} \text{det}(\mathbf{F}) \left[ \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}}
    //     \right]
    // @f]
    // @f{align}
    //  \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right)
    // \dealcoloneq 2 \frac{d \psi_{0} \left( \mathbf{C},
    // \boldsymbol{\mathbb{H}} \right)}{d \mathbf{C}}
    // &= \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //     \left[ \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
    //     - 2 \frac{1}{\text{det}(\mathbf{F})}
    //     \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \right]
    // + 4 \lambda_{e} \ln \left(\text{det}(\mathbf{F}) \right)
    //     \frac{1}{\text{det}(\mathbf{F})} \frac{d\,\text{det}(\mathbf{F})}{d
    //     \mathbf{C}}
    // - \mu_{0} \mu_{r} \left[
    //     \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //     \boldsymbol{\mathbb{H}} \right] \frac{d\,\text{det}(\mathbf{F})}{d
    //     \mathbf{C}} + \text{det}(\mathbf{F}) \frac{d \left[
    //     \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //     \boldsymbol{\mathbb{H}}
    //       \right]}{d \mathbf{C}} \right]
    // \\ &= \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //     \left[ \mathbf{I} - \mathbf{C}^{-1} \right]
    // + 2 \lambda_{e} \ln \left(\text{det}(\mathbf{F}) \right) \mathbf{C}^{-1}
    // - \mu_{0} \mu_{r} \left[
    //     \frac{1}{2}  \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //     \cdot \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F})
    //     \mathbf{C}^{-1}
    // - \text{det}(\mathbf{F})
    //     \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //       \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //       \right]
    // @f}
    // with
    // @f[
    //   \frac{d f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)}{d
    //   \boldsymbol{\mathbb{H}}}
    // = \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //   \text{sech}^{2} \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //   \boldsymbol{\mathbb{H}}}
    //     {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //   \left[ \frac{4} {\left(h_{e}^{\text{sat}}\right)^{2}}
    //   \boldsymbol{\mathbb{H}} \right]
    // @f]
    // @f[
    //   \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
    // = \mathbf{I}
    // \quad \text{(the second-order identity tensor)}
    // @f]
    // @f[
    //   \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}}
    // = \frac{1}{2} \text{det}(\mathbf{F}) \mathbf{C}^{-1}
    // @f]
    // @f[
    // \frac{d C^{-1}_{ab}}{d C_{cd}}
    // = - \text{sym} \left( C^{-1}_{ac} C^{-1}_{bd} \right)
    // = -\frac{1}{2} \left[ C^{-1}_{ac} C^{-1}_{bd} + C^{-1}_{ad} C^{-1}_{bc}
    // \right]
    // @f]
    // @f[
    //   \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //   \boldsymbol{\mathbb{H}} \right]}{d \mathbf{C}}
    // = - \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //   \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    // @f]
    // The use of the symmetry operator $\text{sym} \left( \bullet \right)$ in
    // the one derivation above helps to ensure that the resulting rank-4
    // tensor, which holds minor symmetries due to the symmetry of $\mathbf{C}$,
    // still maps rank-2 symmetric tensors to rank-2 symmetric tensors. See the
    // SymmetricTensor class documentation and the introduction to step-44 and
    // for further explanation as to what symmetry means in the context of
    // fourth-order tensors.
    //
    // The linearization of each of the kinematic variables with respect to
    // their arguments are
    // @f[
    // \mathbb{D} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // = \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}}
    // = - \frac{1}{2} \mu_{e} \left[ \text{tr}(\mathbf{C}) - d - 2 \ln
    // (\text{det}(\mathbf{F}))
    //     \right] \frac{d^{2} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}}
    //     \right)}{d \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}}
    // + \mu_{0} \mu_{r} \text{det}(\mathbf{F}) \mathbf{C}^{-1}
    // @f]
    // @f{align}
    // \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    // \right) = - \frac{d \mathbf{S}^{\text{tot}}}{d \boldsymbol{\mathbb{H}}}
    // &= - \mu_{e}
    //     \left[ \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
    //     - 2 \frac{1}{\text{det}(\mathbf{F})}
    //     \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \right]
    //       \otimes \frac{d f_{\mu_{e} \left( \boldsymbol{\mathbb{H}}
    //       \right)}}{d \boldsymbol{\mathbb{H}}}
    // + \mu_{0} \mu_{r} \left[
    //     \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \otimes
    //       \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //       \boldsymbol{\mathbb{H}}
    //         \right]}{d \boldsymbol{\mathbb{H}}} \right]
    // + \text{det}(\mathbf{F})
    //     \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //     \cdot \boldsymbol{\mathbb{H}}
    //       \right]}{d \mathbf{C} \otimes d \boldsymbol{\mathbb{H}}}
    // \\ &= - \mu_{e}
    //     \left[ \mathbf{I} - \mathbf{C}^{-1} \right] \otimes
    //       \frac{d f_{\mu_{e} \left( \boldsymbol{\mathbb{H}} \right)}}{d
    //       \boldsymbol{\mathbb{H}}}
    // + \mu_{0} \mu_{r} \left[
    //     \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes
    //       \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //       \right]
    // + \text{det}(\mathbf{F})
    //     \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //     \cdot \boldsymbol{\mathbb{H}}
    //       \right]}{d \mathbf{C} \otimes \mathbf{C} \boldsymbol{\mathbb{H}}}
    // @f}
    // @f{align}
    // \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    // \right) = 2 \frac{d \mathbf{S}^{\text{tot}}}{d \mathbf{C}}
    // &= 2 \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //     \left[ - \frac{d \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    //   + 4 \lambda_{e} \left[ \mathbf{C}^{-1} \otimes \left[
    //   \frac{1}{\text{det}(\mathbf{F})} \frac{d \, \text{det}(\mathbf{F})}{d
    //   \mathbf{C}} \right] + \ln \left(\text{det}(\mathbf{F}) \right) \frac{d
    //   \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    // \\ &- \mu_{0} \mu_{r}  \left[
    //  \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes \frac{d \left[
    //  \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}} \right]}{d \mathbf{C}}
    // + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}} \right] \mathbf{C}^{-1} \otimes \frac{d \,
    // \text{det}(\mathbf{F})}{d \mathbf{C}}
    // + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F}) \frac{d
    // \mathbf{C}^{-1}}{d \mathbf{C}}
    // \right]
    // \\ &+ 2 \mu_{0} \mu_{r} \left[ \left[
    //     \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //       \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //       \right] \otimes \frac{d \, \text{det}(\mathbf{F})}{d \mathbf{C}}
    //     - \text{det}(\mathbf{F})
    //     \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //     \cdot \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
    //     \mathbf{C}}
    // \right]
    // \\ &= 2 \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //     \left[ - \frac{d \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    //  + 4 \lambda_{e} \left[ \frac{1}{2} \mathbf{C}^{-1} \otimes
    //  \mathbf{C}^{-1} + \ln \left(\text{det}(\mathbf{F}) \right) \frac{d
    //  \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    // \\ &- \mu_{0} \mu_{r}  \left[
    //  - \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes \left[ \left[
    //  \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //   \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \right]
    // + \frac{1}{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F})  \mathbf{C}^{-1}
    // \otimes \mathbf{C}^{-1}
    // + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F}) \frac{d
    // \mathbf{C}^{-1}}{d \mathbf{C}}
    // \right]
    // \\ &+ 2 \mu_{0} \mu_{r} \left[ \frac{1}{2} \text{det}(\mathbf{F}) \left[
    //     \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //       \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //       \right] \otimes \mathbf{C}^{-1}
    //     - \text{det}(\mathbf{F})
    //     \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //     \cdot \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
    //     \mathbf{C}}
    // \right]
    // @f}
    // with
    // @f[
    //  \frac{d^{2} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)}{d
    //  \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}}
    // = -2 \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //   \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //   \boldsymbol{\mathbb{H}}}
    //     {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //   \text{sech}^{2} \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //   \boldsymbol{\mathbb{H}}}
    //     {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //   \left[ \frac{4} {\left(h_{e}^{\text{sat}}\right)^{2}} \mathbf{I}
    //   \right]
    // @f]
    // @f[
    // \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}}
    //         \right]}{d \boldsymbol{\mathbb{H}}}
    // = 2 \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
    // @f]
    // @f[
    // \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
    // \boldsymbol{\mathbb{H}}} \Rightarrow \frac{d^{2} \left[ \mathbb{H}_{e}
    // C^{-1}_{ef} \mathbb{H}_{f}
    //       \right]}{d C_{ab} d \mathbb{H}_{c}}
    // = - C^{-1}_{ac} C^{-1}_{be} \mathbb{H}_{e} - C^{-1}_{ae} \mathbb{H}_{e}
    // C^{-1}_{bc}
    // @f]
    // @f{align}
    // \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    // \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d \mathbf{C}}
    // &= -\frac{d \left[\left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
    // \right] \otimes
    //       \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
    //       \right]\right]}{d \mathbf{C}}
    // \\ \Rightarrow
    // \frac{d^{2} \left[ \mathbb{H}_{e} C^{-1}_{ef} \mathbb{H}_{f}
    //       \right]}{d C_{ab} d C_{cd}}
    // &= \text{sym} \left( C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{cf}
    // \mathbb{H}_{f} C^{-1}_{bd}
    //           + C^{-1}_{ce} \mathbb{H}_{e} C^{-1}_{bf} \mathbb{H}_{f}
    //           C^{-1}_{ad} \right)
    // \\ &= \frac{1}{2} \left[
    //      C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{cf} \mathbb{H}_{f} C^{-1}_{bd}
    //    + C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{df} \mathbb{H}_{f} C^{-1}_{bc}
    //    + C^{-1}_{ce} \mathbb{H}_{e} C^{-1}_{bf} \mathbb{H}_{f} C^{-1}_{ad}
    //    + C^{-1}_{be} \mathbb{H}_{e} C^{-1}_{df} \mathbb{H}_{f} C^{-1}_{ac}
    //   \right]
    // @f}
    //
    // Well, that escalated quickly -- although the the definition of $\psi_{0}$
    // and $f_{\mu_e}$ might have given some hints that the calculating
    // the kinetic fields and their linearization would take some effort, it is
    // likely that there's a little more complexity to the final definitions
    // that perhaps initially thought.
    // Knowing what we now do, it's probably fair to say that we really do not
    // want to compute first and second derivatives of these functions with
    // respect to their arguments -- regardless of well we did in calculus
    // classes, or how good a programmer we may be.
    //
    // In the class method definition where these are ultimately implemented,
    // we've composed these calculations slightly differently. Some intermediate
    // steps are also retained to give another perspective of how to
    // systematically compute the derivatives. Additionally, some calculations
    // are decomposed less or further to reuse some of the intermediate values
    // and, hopefully, aid the reader to follow the derivative operations.
    template <int dim>
    class Magnetoelastic_Constitutive_Law final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
    public:
      Magnetoelastic_Constitutive_Law(
        const ConstitutiveParameters &constitutive_parameters);

      virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
                                        const Tensor<1, dim> &         H,
                                        const DiscreteTime &) override;

      virtual double get_psi() const override;

      virtual Tensor<1, dim> get_B() const override;

      virtual SymmetricTensor<2, dim> get_S() const override;

      virtual SymmetricTensor<2, dim> get_DD() const override;

      virtual Tensor<3, dim> get_PP() const override;

      virtual SymmetricTensor<4, dim> get_HH() const override;

    private:
      double                  psi;
      Tensor<1, dim>          B;
      SymmetricTensor<2, dim> S;
      SymmetricTensor<2, dim> BB;
      Tensor<3, dim>          PP;
      SymmetricTensor<4, dim> HH;
    };


    template <int dim>
    Magnetoelastic_Constitutive_Law<dim>::Magnetoelastic_Constitutive_Law(
      const ConstitutiveParameters &constitutive_parameters)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , psi(0.0)
    {}

    // For this class's update method, we'll simply precompute a collection of
    // intermediate values (for function evaluations, derivative calculations,
    // and the like) and "manually" arrange them in the order that's required
    // to maximize their reuse. This means that we have to manage this
    // ourselves, and decide what values must be compute before others, all
    // while keeping some semblance of order or structure in the code itself.
    // It's effective, but perhaps a little tedious. It also doesn't do too much
    // to help future extension of the class, because all of these values remain
    // local to this single method.
    //
    // Interestingly, this basic technique of precomputing intermediate
    // expressions that are used in more than one place has a name:
    // [common subexpression elimination
    // (CSE)](https://en.wikipedia.org/wiki/Common_subexpression_elimination).
    // It is a strategy used by Computer Algebra Systems to reduce the
    // computational expense when they are tasked with evaluating similar
    // expressions.
    template <int dim>
    void Magnetoelastic_Constitutive_Law<dim>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &)
    {
      const double                  det_F = std::sqrt(determinant(C));
      const SymmetricTensor<2, dim> C_inv = invert(C);
      AssertThrow(det_F > 0.0,
                  ExcMessage("Volumetric Jacobian must be positive."));

      // The saturation function for the magneto-elastic energy.
      const double two_h_dot_h_div_h_sat_squ =
        (2.0 * H * H) / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat());
      const double tanh_two_h_dot_h_div_h_sat_squ =
        std::tanh(two_h_dot_h_div_h_sat_squ);

      const double f_mu_e =
        1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
                tanh_two_h_dot_h_div_h_sat_squ;

      // The first derivative of the saturation function, noting that
      // $\frac{d \tanh(x)}{dx} = \text{sech}^{2}(x)$.
      const double dtanh_two_h_dot_h_div_h_sat_squ =
        std::pow(1.0 / std::cosh(two_h_dot_h_div_h_sat_squ), 2.0);
      const Tensor<1, dim> dtwo_h_dot_h_div_h_sat_squ_dH =
        2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) * H;

      const Tensor<1, dim> df_mu_e_dH =
        (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
        (dtanh_two_h_dot_h_div_h_sat_squ * dtwo_h_dot_h_div_h_sat_squ_dH);

      // The second derivative of saturation function, noting that
      // $\frac{d \text{sech}^{2}(x)}{dx} = -2 \tanh(x) \text{sech}^{2}(x)$.
      const double d2tanh_two_h_dot_h_div_h_sat_squ =
        -2.0 * tanh_two_h_dot_h_div_h_sat_squ * dtanh_two_h_dot_h_div_h_sat_squ;
      const SymmetricTensor<2, dim> d2two_h_dot_h_div_h_sat_squ_dH_dH =
        2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) *
        Physics::Elasticity::StandardTensors<dim>::I;

      const SymmetricTensor<2, dim> d2f_mu_e_dH_dH =
        (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
        (d2tanh_two_h_dot_h_div_h_sat_squ *
           symmetrize(outer_product(dtwo_h_dot_h_div_h_sat_squ_dH,
                                    dtwo_h_dot_h_div_h_sat_squ_dH)) +
         dtanh_two_h_dot_h_div_h_sat_squ * d2two_h_dot_h_div_h_sat_squ_dH_dH);

      // Some intermediate quantities attained directly from the
      // field / kinematic variables.
      const double         log_det_F         = std::log(det_F);
      const double         tr_C              = trace(C);
      const Tensor<1, dim> C_inv_dot_H       = C_inv * H;
      const double         H_dot_C_inv_dot_H = H * C_inv_dot_H;

      // First derivatives of the intermediate quantities.
      const SymmetricTensor<2, dim> d_tr_C_dC =
        Physics::Elasticity::StandardTensors<dim>::I;
      const SymmetricTensor<2, dim> ddet_F_dC     = 0.5 * det_F * C_inv;
      const SymmetricTensor<2, dim> dlog_det_F_dC = 0.5 * C_inv;

      const Tensor<1, dim> dH_dot_C_inv_dot_H_dH = 2.0 * C_inv_dot_H;

      SymmetricTensor<4, dim> dC_inv_dC;
      for (unsigned int A = 0; A < dim; ++A)
        for (unsigned int B = A; B < dim; ++B)
          for (unsigned int C = 0; C < dim; ++C)
            for (unsigned int D = C; D < dim; ++D)
              dC_inv_dC[A][B][C][D] -=               //
                0.5 * (C_inv[A][C] * C_inv[B][D]     //
                       + C_inv[A][D] * C_inv[B][C]); //

      const SymmetricTensor<2, dim> dH_dot_C_inv_dot_H_dC =
        -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H));

      // Second derivatives of the intermediate quantities.
      const SymmetricTensor<4, dim> d2log_det_F_dC_dC = 0.5 * dC_inv_dC;

      const SymmetricTensor<4, dim> d2det_F_dC_dC =
        0.5 * (outer_product(C_inv, ddet_F_dC) + det_F * dC_inv_dC);

      const SymmetricTensor<2, dim> d2H_dot_C_inv_dot_H_dH_dH = 2.0 * C_inv;

      Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH;
      for (unsigned int A = 0; A < dim; ++A)
        for (unsigned int B = 0; B < dim; ++B)
          for (unsigned int C = 0; C < dim; ++C)
            d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -=
              C_inv[A][C] * C_inv_dot_H[B] + //
              C_inv_dot_H[A] * C_inv[B][C];  //

      SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC;
      for (unsigned int A = 0; A < dim; ++A)
        for (unsigned int B = A; B < dim; ++B)
          for (unsigned int C = 0; C < dim; ++C)
            for (unsigned int D = C; D < dim; ++D)
              d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] +=
                0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] +
                       C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] +
                       C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] +
                       C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]);

      // The stored energy density function.
      psi =
        (0.5 * this->get_mu_e() * f_mu_e) *
          (tr_C - dim - 2.0 * std::log(det_F)) +
        this->get_lambda_e() * (std::log(det_F) * std::log(det_F)) -
        (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H);

      // The kinetic quantities.
      B = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) *
            df_mu_e_dH //
          + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
              dH_dot_C_inv_dot_H_dH; //

      S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                        //
            (d_tr_C_dC - 2.0 * dlog_det_F_dC)                              //
          + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC) //
          - 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *            //
              (H_dot_C_inv_dot_H * ddet_F_dC                               //
               + det_F * dH_dot_C_inv_dot_H_dC);                           //

      // The linearization of the kinetic quantities.
      BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * //
             d2f_mu_e_dH_dH                                             //
           + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
               d2H_dot_C_inv_dot_H_dH_dH; //

      PP = -2.0 * (0.5 * this->get_mu_e()) *                                  //
             outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC),   //
                           df_mu_e_dH)                                        //
           +                                                                  //
           2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                //
             (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) //
              + det_F * d2H_dot_C_inv_dot_H_dC_dH);                           //

      HH =
        4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) //
        + 4.0 * this->get_lambda_e() *                                       //
            (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               //
             + 2.0 * log_det_F * d2log_det_F_dC_dC)                          //
        - 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                //
            (H_dot_C_inv_dot_H * d2det_F_dC_dC                               //
             + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)               //
             + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)               //
             + det_F * d2H_dot_C_inv_dot_H_dC_dC);                           //
    }

    template <int dim>
    double Magnetoelastic_Constitutive_Law<dim>::get_psi() const
    {
      return psi;
    }

    template <int dim>
    Tensor<1, dim> Magnetoelastic_Constitutive_Law<dim>::get_B() const
    {
      return B;
    }

    template <int dim>
    SymmetricTensor<2, dim> Magnetoelastic_Constitutive_Law<dim>::get_S() const
    {
      return S;
    }

    template <int dim>
    SymmetricTensor<2, dim> Magnetoelastic_Constitutive_Law<dim>::get_DD() const
    {
      return BB;
    }

    template <int dim>
    Tensor<3, dim> Magnetoelastic_Constitutive_Law<dim>::get_PP() const
    {
      return PP;
    }

    template <int dim>
    SymmetricTensor<4, dim> Magnetoelastic_Constitutive_Law<dim>::get_HH() const
    {
      return HH;
    }


    // @sect4{Magneto-viscoelastic constitutive law (hand-derived)}

    // As mentioned before, the free energy density function for the
    // magneto-viscoelastic material with one dissipative mechanism that we'll
    // be considering is defined as
    // @f[
    //   \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
    //   \right)
    // = \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // + \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)
    // @f]
    // @f[
    //   \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // = \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
    // \right)
    //     \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //     \right]
    // + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    // - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //     \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //     \boldsymbol{\mathbb{H}} \right]
    // @f]
    // @f[
    //   \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //   \boldsymbol{\mathbb{H}} \right)
    // = \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    // \right)
    //     \left[ \mathbf{C}_{v} : \left[
    //       \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //       \mathbf{C} \right] - d - \ln\left(
    //       \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    // @f]
    // with
    // @f[
    //   f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
    // = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //     \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //     \boldsymbol{\mathbb{H}}}
    //       {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    // @f]
    // @f[
    //   f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
    // = 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
    //     \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //     \boldsymbol{\mathbb{H}}}
    //       {\left(h_{v}^{\text{sat}}\right)^{2}} \right)
    // @f]
    // and the evolution law
    // @f[
    //  \dot{\mathbf{C}}_{v} \left( \mathbf{C} \right)
    // = \frac{1}{\tau} \left[
    //       \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //         \mathbf{C}\right]^{-1}
    //     - \mathbf{C}_{v} \right]
    // @f]
    // that itself is parameterized in terms of $\mathbf{C}$.
    // By design, the magnetoelastic part of the energy
    // $\psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)$
    // is identical to that of the magnetoelastic material presented earlier.
    // So, for the derivatives of the various contributions stemming from this
    // part of the energy, please refer to the previous section. We'll continue
    // to highlight the specific contributions from those terms by
    // superscripting the salient terms with $ME$, while contributions from the
    // magneto-viscoelastic component are superscripted with $MVE$.
    // Furthermore, the magnetic saturation function
    // $f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)$
    // for the damping term has the identical form as that of the elastic
    // term (i.e.,
    // $f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)$
    // ), and so the structure of its derivatives are identical to that
    // seen before; the only change is for the three constitutive parameters
    // that are now associated with the viscous shear modulus $\mu_{v}$ rather
    // than the elastic shear modulus $\mu_{e}$.
    //
    // For this magneto-viscoelastic material, the first derivatives that
    // correspond to the magnetic induction vector and total Piola-Kirchhoff
    // stress tensor are
    // @f[
    //  \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)
    // \dealcoloneq - \frac{\partial \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
    // \Big\vert_{\mathbf{C}, \mathbf{C}_{v}} \equiv
    // \boldsymbol{\mathbb{B}}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    // \right)
    // + \boldsymbol{\mathbb{B}}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) =  - \frac{d \psi_{0}^{ME} \left(
    // \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}}}
    //    - \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //    \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
    // @f]
    // @f[
    //  \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)
    // \dealcoloneq 2 \frac{\partial \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
    // \Big\vert_{\mathbf{C}_{v}, \boldsymbol{\mathbb{H}}} \equiv
    // \mathbf{S}^{\text{tot}, ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    // \right)
    // + \mathbf{S}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}}
    //     \right)
    // =  2 \frac{d \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    // \right)}{d \mathbf{C}}
    //  + 2 \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
    // @f]
    // with the viscous contributions being
    // @f[
    //   \boldsymbol{\mathbb{B}}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //   \boldsymbol{\mathbb{H}} \right)
    // = - \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
    // \Big\vert_{\mathbf{C}, \mathbf{C}_{v}} = - \frac{1}{2} \mu_{v}
    //     \left[ \mathbf{C}_{v} : \left[
    //       \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //       \mathbf{C} \right] - d - \ln\left(
    //       \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    //       \frac{\partial f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //       \right)}{\partial \boldsymbol{\mathbb{H}}}
    // @f]
    // @f[
    //   \mathbf{S}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //   \boldsymbol{\mathbb{H}}
    //     \right)
    // = 2 \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
    // \Big\vert_{\mathbf{C}_{v}, \boldsymbol{\mathbb{H}}} = \mu_{v}
    // f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //        \left[  \left[ \mathbf{C}_{v} : \mathbf{C} \right] \left[ -
    //        \frac{1}{d}
    //        \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C}^{-1} \right]
    //        + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C}_{v}
    //  \right]
    // @f]
    // and with
    // @f[
    // \frac{\partial f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    // \right)}{\partial \boldsymbol{\mathbb{H}}} \equiv \frac{d
    // f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)}{d
    // \boldsymbol{\mathbb{H}}} .
    // @f]
    // The time-discretized evolution law,
    // @f[
    // \mathbf{C}_{v}^{(t)} \left( \mathbf{C} \right)
    // = \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
    //     \mathbf{C}_{v}^{(t-1)}
    //   + \frac{\Delta t}{\tau_{v}}
    //     \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //     \mathbf{C} \right]^{-1}
    //   \right]
    // @f]
    // will also dictate how the linearization of the internal
    // variable with respect to the field variables is composed.
    //
    // Observe that in order to attain the *correct* expressions for the
    // magnetic induction vector and total Piola-Kirchhoff stress tensor for
    // this dissipative material, we must adhere strictly to the outcome of
    // applying the Coleman-Noll procedure: we must take *partial derivatives*
    // of the free energy density function with respect to the field variables.
    // (For our non-dissipative magnetoelastic material, taking either partial
    // or total derivatives would have had the same result, so there was no
    // need to draw your attention to this before.)
    // The crucial part of the operation is to freeze the internal variable
    // $\mathbf{C}_{v}^{(t)} \left( \mathbf{C} \right)$ while computing the
    // derivatives of $\psi_{0}^{MVE} \left( \mathbf{C},
    // \mathbf{C}_{v} \left( \mathbf{C} \right), \boldsymbol{\mathbb{H}}
    // \right)$ with respect to $\mathbf{C}$ -- the dependence of
    // $\mathbf{C}_{v}^{(t)}$ on $\mathbf{C}$ is not to be taken into account.
    // When deciding whether to use AD or SD to perform this task
    // the choice is clear -- only the symbolic framework provides a mechanism
    // to do this; as was mentioned before, AD can only return total derivatives
    // so it is unsuitable for the task.
    //
    // To wrap things up, we'll present the material tangents for this
    // rate-dependent coupled material. The linearization of both kinetic
    // variables with respect to their arguments are
    // @f[
    // \mathbb{D} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
    // \right) = \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}}
    // \equiv \mathbb{D}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // + \mathbb{D}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = \frac{d \boldsymbol{\mathbb{B}}^{ME}}{d
    // \boldsymbol{\mathbb{H}}}
    // + \frac{d \boldsymbol{\mathbb{B}}^{MVE}}{d \boldsymbol{\mathbb{H}}}
    // @f]
    // @f[
    // \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = - \frac{d \mathbf{S}^{\text{tot}}}{d
    // \boldsymbol{\mathbb{H}}} \equiv \mathfrak{P}^{\text{tot}, ME} \left(
    // \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    // + \mathfrak{P}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = - \frac{d \mathbf{S}^{\text{tot},
    // ME}}{d \boldsymbol{\mathbb{H}}}
    // - \frac{d \mathbf{S}^{\text{tot}, MVE}}{d \boldsymbol{\mathbb{H}}}
    // @f]
    // @f[
    // \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = 2 \frac{d \mathbf{S}^{\text{tot}}}{d
    // \mathbf{C}} \equiv \mathcal{H}^{\text{tot}, ME} \left( \mathbf{C},
    // \boldsymbol{\mathbb{H}} \right)
    // + \mathcal{H}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = 2 \frac{d \mathbf{S}^{\text{tot},
    // ME}}{d \mathbf{C}}
    // + 2 \frac{d \mathbf{S}^{\text{tot}, MVE}}{d \mathbf{C}}
    // @f]
    // where the tangents for the viscous contributions are
    // @f[
    // \mathbb{D}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = - \frac{1}{2} \mu_{v}
    //     \left[ \mathbf{C}_{v} : \left[
    //       \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //       \mathbf{C} \right] - d - \ln\left(
    //       \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    //       \frac{\partial^{2} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //       \right)}{\partial \boldsymbol{\mathbb{H}} \otimes
    //       d \boldsymbol{\mathbb{H}}}
    // @f]
    // @f[
    // \mathfrak{P}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right) = - \mu_{v}
    //        \left[  \left[ \mathbf{C}_{v} : \mathbf{C} \right] \left[ -
    //        \frac{1}{d}
    //        \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C}^{-1} \right]
    //        + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C}_{v}
    //  \right] \otimes \frac{d f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //  \right)}{d \boldsymbol{\mathbb{H}}}
    // @f]
    // @f{align}
    // \mathcal{H}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    // \boldsymbol{\mathbb{H}} \right)
    // &= 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //   \left[ - \frac{1}{d}
    //   \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //   \mathbf{C}^{-1} \right] \otimes
    //   \left[ \mathbf{C}_{v} + \mathbf{C} : \frac{d \mathbf{C}_{v}}{d
    //   \mathbf{C}} \right]
    // \\ &+ 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    // \left[ \mathbf{C}_{v} : \mathbf{C} \right]
    //   \left[
    //     \frac{1}{d^{2}}
    //     \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //     \mathbf{C}^{-1} \otimes \mathbf{C}^{-1}
    //     - \frac{1}{d}
    //     \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}} \frac{d
    //     \mathbf{C}^{-1}}{d \mathbf{C}}
    //   \right]
    // \\ &+ 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //   \left[
    //     -\frac{1}{d}
    //     \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //     \mathbf{C}_{v} \otimes \mathbf{C}^{-1}
    //     + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //     \frac{d \mathbf{C}_{v}}{d \mathbf{C}}
    //   \right]
    // @f}
    // with
    // @f[
    // \frac{\partial^{2} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    // \right)}{\partial \boldsymbol{\mathbb{H}} \otimes
    // d \boldsymbol{\mathbb{H}}} \equiv \frac{d^{2} f_{\mu_{v}^{MVE}} \left(
    // \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}} \otimes d
    // \boldsymbol{\mathbb{H}}}
    // @f]
    // and, from the evolution law,
    // @f[
    // \frac{d \mathbf{C}_{v}}{d \mathbf{C}}
    // \equiv \frac{d \mathbf{C}_{v}^{(t)}}{d \mathbf{C}}
    //  = \frac{\frac{\Delta t}{\tau_{v}} }{1 + \frac{\Delta t}{\tau_{v}}}
    //  \left[
    //     \frac{1}{d}
    //     \left[\text{det}\left(\mathbf{F}\right)\right]^{\frac{2}{d}}
    //     \mathbf{C}^{-1} \otimes \mathbf{C}^{-1}
    //    + \left[\text{det}\left(\mathbf{F}\right)\right]^{\frac{2}{d}} \frac{d
    //    \mathbf{C}^{-1}}{d \mathbf{C}}
    //   \right] .
    // @f]
    // Notice that just the last term of $\mathcal{H}^{\text{tot}, MVE}$
    // contains the tangent of the internal variable. The linearization of this
    // particular evolution law is linear. For an example of a nonlinear
    // evolution law, for which this linearization must be solved for in an
    // iterative manner, see @cite Koprowski-Theiss2011a.
    template <int dim>
    class Magnetoviscoelastic_Constitutive_Law final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
    public:
      Magnetoviscoelastic_Constitutive_Law(
        const ConstitutiveParameters &constitutive_parameters);

      virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
                                        const Tensor<1, dim> &         H,
                                        const DiscreteTime &time) override;

      virtual double get_psi() const override;

      virtual Tensor<1, dim> get_B() const override;

      virtual SymmetricTensor<2, dim> get_S() const override;

      virtual SymmetricTensor<2, dim> get_DD() const override;

      virtual Tensor<3, dim> get_PP() const override;

      virtual SymmetricTensor<4, dim> get_HH() const override;

      virtual void update_end_of_timestep() override;

    private:
      SymmetricTensor<2, dim> Q_t;
      SymmetricTensor<2, dim> Q_t1;

      double                  psi;
      Tensor<1, dim>          B;
      SymmetricTensor<2, dim> S;
      SymmetricTensor<2, dim> BB;
      Tensor<3, dim>          PP;
      SymmetricTensor<4, dim> HH;

      // A data structure that is used to store all intermediate calculations.
      // We'll see shortly precisely how this can be leveraged to make the part
      // of the code where we actually perform calculations clean and easy
      // (well, at least easier) to follow and maintain. But for now, we can say
      // that it will allow us to move the parts of the code where we compute
      // the derivatives of intermediate quantities away from where they are
      // used.
      mutable GeneralDataStorage cache;

      // The next two functions are used to update the state of the field and
      // internal variables, and will be called before we perform any
      // detailed calculations.
      void set_primary_variables(const SymmetricTensor<2, dim> &C,
                                 const Tensor<1, dim> &         H) const;

      void update_internal_variable(const DiscreteTime &time);

      // The remainder of the class interface is dedicated to methods that
      // are used to compute the components required to calculate the free
      // energy density function, and all of its derivatives:

      // The kinematic, or field, variables.
      const Tensor<1, dim> &get_H() const;

      const SymmetricTensor<2, dim> &get_C() const;

      // A generalized formulation for the saturation function, with the
      // required constitutive parameters passed as arguments to each function.
      double get_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      double get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      double get_f_mu(const double mu,
                      const double mu_inf,
                      const double mu_h_sat) const;

      // A generalized formulation for the first derivative of saturation
      // function, with the required constitutive parameters passed as arguments
      // to each function.
      double get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      Tensor<1, dim>
      get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const;

      Tensor<1, dim> get_df_mu_dH(const double mu,
                                  const double mu_inf,
                                  const double mu_h_sat) const;

      // A generalized formulation for the second derivative of saturation
      // function, with the required constitutive parameters passed as arguments
      // to each function.
      double get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      SymmetricTensor<2, dim>
      get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const;

      SymmetricTensor<2, dim> get_d2f_mu_dH_dH(const double mu,
                                               const double mu_inf,
                                               const double mu_h_sat) const;

      // Intermediate quantities attained directly from the
      // field / kinematic variables.
      const double &get_det_F() const;

      const SymmetricTensor<2, dim> &get_C_inv() const;

      const double &get_log_det_F() const;

      const double &get_trace_C() const;

      const Tensor<1, dim> &get_C_inv_dot_H() const;

      const double &get_H_dot_C_inv_dot_H() const;

      // First derivatives of the intermediate quantities.
      const SymmetricTensor<4, dim> &get_dC_inv_dC() const;

      const SymmetricTensor<2, dim> &get_d_tr_C_dC() const;

      const SymmetricTensor<2, dim> &get_ddet_F_dC() const;

      const SymmetricTensor<2, dim> &get_dlog_det_F_dC() const;

      const Tensor<1, dim> &get_dH_dot_C_inv_dot_H_dH() const;

      const SymmetricTensor<2, dim> &get_dH_dot_C_inv_dot_H_dC() const;

      // Derivative of internal variable with respect to field variables.
      // Notice that we only need this one derivative of the internal variable,
      // as this variable is only differentiated as part of the linearization
      // of the kinetic variables.
      const SymmetricTensor<4, dim> &
      get_dQ_t_dC(const DiscreteTime &time) const;

      // Second derivatives of the intermediate quantities.
      const SymmetricTensor<4, dim> &get_d2log_det_F_dC_dC() const;

      const SymmetricTensor<4, dim> &get_d2det_F_dC_dC() const;

      const SymmetricTensor<2, dim> &get_d2H_dot_C_inv_dot_H_dH_dH() const;

      const Tensor<3, dim> &get_d2H_dot_C_inv_dot_H_dC_dH() const;

      const SymmetricTensor<4, dim> &get_d2H_dot_C_inv_dot_H_dC_dC() const;
    };


    template <int dim>
    Magnetoviscoelastic_Constitutive_Law<
      dim>::Magnetoviscoelastic_Constitutive_Law(const ConstitutiveParameters
                                                   &constitutive_parameters)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , Q_t(Physics::Elasticity::StandardTensors<dim>::I)
      , Q_t1(Physics::Elasticity::StandardTensors<dim>::I)
      , psi(0.0)
    {}


    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &           time)
    {
      // Record the applied deformation state as well as the magnetic load.
      // Thereafter, update internal (viscous) variable based on new deformation
      // state.
      set_primary_variables(C, H);
      update_internal_variable(time);

      // Get the values for the elastic and viscous saturation function based
      // on the current magnetic field...
      const double f_mu_e = get_f_mu(this->get_mu_e(),
                                     this->get_mu_e_inf(),
                                     this->get_mu_e_h_sat());

      const double f_mu_v = get_f_mu(this->get_mu_v(),
                                     this->get_mu_v_inf(),
                                     this->get_mu_v_h_sat());

      // ... as well as their first derivatives...
      const Tensor<1, dim> df_mu_e_dH = get_df_mu_dH(this->get_mu_e(),
                                                     this->get_mu_e_inf(),
                                                     this->get_mu_e_h_sat());

      const Tensor<1, dim> df_mu_v_dH = get_df_mu_dH(this->get_mu_v(),
                                                     this->get_mu_v_inf(),
                                                     this->get_mu_v_h_sat());


      // ... and their second derivatives.
      const SymmetricTensor<2, dim> d2f_mu_e_dH_dH =
        get_d2f_mu_dH_dH(this->get_mu_e(),
                         this->get_mu_e_inf(),
                         this->get_mu_e_h_sat());

      const SymmetricTensor<2, dim> d2f_mu_v_dH_dH =
        get_d2f_mu_dH_dH(this->get_mu_v(),
                         this->get_mu_v_inf(),
                         this->get_mu_v_h_sat());

      // Intermediate quantities. Note that, since we're fetching these values
      // from a cache that has a lifetime that outlasts this function call, we
      // can alias the result rather than copying the value from the cache.
      const double &                 det_F = get_det_F();
      const SymmetricTensor<2, dim> &C_inv = get_C_inv();

      const double &log_det_F         = get_log_det_F();
      const double &tr_C              = get_trace_C();
      const double &H_dot_C_inv_dot_H = get_H_dot_C_inv_dot_H();

      // First derivatives of intermediate values, as well as the that of the
      // internal variable with respect to the right Cauchy-Green deformation
      // tensor.
      const SymmetricTensor<2, dim> &d_tr_C_dC     = get_d_tr_C_dC();
      const SymmetricTensor<2, dim> &ddet_F_dC     = get_ddet_F_dC();
      const SymmetricTensor<2, dim> &dlog_det_F_dC = get_dlog_det_F_dC();

      const SymmetricTensor<4, dim> &dQ_t_dC = get_dQ_t_dC(time);

      const Tensor<1, dim> &dH_dot_C_inv_dot_H_dH = get_dH_dot_C_inv_dot_H_dH();

      const SymmetricTensor<2, dim> &dH_dot_C_inv_dot_H_dC =
        get_dH_dot_C_inv_dot_H_dC();

      // Second derivatives of intermediate values.
      const SymmetricTensor<4, dim> &d2log_det_F_dC_dC =
        get_d2log_det_F_dC_dC();

      const SymmetricTensor<4, dim> &d2det_F_dC_dC = get_d2det_F_dC_dC();

      const SymmetricTensor<2, dim> &d2H_dot_C_inv_dot_H_dH_dH =
        get_d2H_dot_C_inv_dot_H_dH_dH();

      const Tensor<3, dim> &d2H_dot_C_inv_dot_H_dC_dH =
        get_d2H_dot_C_inv_dot_H_dC_dH();

      const SymmetricTensor<4, dim> &d2H_dot_C_inv_dot_H_dC_dC =
        get_d2H_dot_C_inv_dot_H_dC_dC();

      // Since the definitions of the linearizations become particularly
      // lengthy, we'll decompose the free energy density function into three
      // additive components:
      // - the "Neo-Hookean"-like term,
      // - the rate-dependent term, and
      // - the term that resembles that of the energy stored in the magnetic
      // field.
      //
      // To remain consistent, each of these contributions will be individually
      // added to the variables that we want to compute in that same order.
      //
      // So, first of all this is the energy density function itself:
      psi = (0.5 * this->get_mu_e() * f_mu_e) *
              (tr_C - dim - 2.0 * std::log(det_F)) +
            this->get_lambda_e() * (std::log(det_F) * std::log(det_F));
      psi += (0.5 * this->get_mu_v() * f_mu_v) *
             (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim -
              std::log(determinant(Q_t)));
      psi -=
        (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H);

      // ... followed by the magnetic induction vector and Piola-Kirchhoff
      // stress:
      B =
        -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * df_mu_e_dH;
      B -= (0.5 * this->get_mu_v()) *
           (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim -
            std::log(determinant(Q_t))) *
           df_mu_v_dH;
      B += 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
           dH_dot_C_inv_dot_H_dH;

      S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                         //
            (d_tr_C_dC - 2.0 * dlog_det_F_dC)                               //
          + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC); //
      S += 2.0 * (0.5 * this->get_mu_v() * f_mu_v) *
           ((Q_t * C) *
              ((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) * ddet_F_dC) +
            std::pow(det_F, -2.0 / dim) * Q_t);                // dC/dC = II
      S -= 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * //
           (H_dot_C_inv_dot_H * ddet_F_dC                      //
            + det_F * dH_dot_C_inv_dot_H_dC);                  //

      // ... and lastly the tangents due to the linearization of the kinetic
      // variables.
      BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) *
           d2f_mu_e_dH_dH;
      BB -= (0.5 * this->get_mu_v()) *
            (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim -
             std::log(determinant(Q_t))) *
            d2f_mu_v_dH_dH;
      BB += 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
            d2H_dot_C_inv_dot_H_dH_dH;

      PP = -2.0 * (0.5 * this->get_mu_e()) *
           outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC),
                         df_mu_e_dH);
      PP -= 2.0 * (0.5 * this->get_mu_v()) *
            outer_product(Tensor<2, dim>((Q_t * C) *
                                           ((-2.0 / dim) *
                                            std::pow(det_F, -2.0 / dim - 1.0) *
                                            ddet_F_dC) +
                                         std::pow(det_F, -2.0 / dim) * Q_t),
                          df_mu_v_dH);
      PP += 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *
            (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) +
             det_F * d2H_dot_C_inv_dot_H_dC_dH);

      HH =
        4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) //
        + 4.0 * this->get_lambda_e() *                                       //
            (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               //
             + 2.0 * log_det_F * d2log_det_F_dC_dC);                         //
      HH += 4.0 * (0.5 * this->get_mu_v() * f_mu_v) *
            (outer_product((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) *
                             ddet_F_dC,
                           C * dQ_t_dC + Q_t) +
             (Q_t * C) *
               (outer_product(ddet_F_dC,
                              (-2.0 / dim) * (-2.0 / dim - 1.0) *
                                std::pow(det_F, -2.0 / dim - 2.0) * ddet_F_dC) +
                ((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) *
                 d2det_F_dC_dC)) +
             outer_product(Q_t,
                           (-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) *
                             ddet_F_dC) +
             std::pow(det_F, -2.0 / dim) * dQ_t_dC);
      HH -= 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * //
            (H_dot_C_inv_dot_H * d2det_F_dC_dC                  //
             + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)  //
             + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)  //
             + det_F * d2H_dot_C_inv_dot_H_dC_dC);              //


      // Now that we're done using all of those temporary variables stored
      // in our cache, we can clear it out to free up some memory.
      cache.reset();
    }

    template <int dim>
    double Magnetoviscoelastic_Constitutive_Law<dim>::get_psi() const
    {
      return psi;
    }


    template <int dim>
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_B() const
    {
      return B;
    }


    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_S() const
    {
      return S;
    }


    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_DD() const
    {
      return BB;
    }


    template <int dim>
    Tensor<3, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_PP() const
    {
      return PP;
    }


    template <int dim>
    SymmetricTensor<4, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_HH() const
    {
      return HH;
    }


    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law<dim>::update_end_of_timestep()
    {
      Q_t1 = Q_t;
    };


    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_variable(
      const DiscreteTime &time)
    {
      const double delta_t = this->get_delta_t(time);

      Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) *
            (Q_t1 + (delta_t / this->get_tau_v()) *
                      std::pow(get_det_F(), 2.0 / dim) * get_C_inv());
    }

    // The next few functions implement the generalized formulation for the
    // saturation function, as well as its various derivatives.
    template <int dim>
    double
    Magnetoviscoelastic_Constitutive_Law<dim>::get_two_h_dot_h_div_h_sat_squ(
      const double mu_h_sat) const
    {
      const Tensor<1, dim> &H = get_H();
      return (2.0 * H * H) / (mu_h_sat * mu_h_sat);
    };


    template <int dim>
    double Magnetoviscoelastic_Constitutive_Law<
      dim>::get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
    {
      return std::tanh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat));
    };

    // A scaling function that will cause the shear modulus
    // to change (increase) under the influence of a magnetic
    // field.
    template <int dim>
    double Magnetoviscoelastic_Constitutive_Law<dim>::get_f_mu(
      const double mu,
      const double mu_inf,
      const double mu_h_sat) const
    {
      return 1.0 +
             (mu_inf / mu - 1.0) * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat);
    };

    // First derivative of scaling function
    template <int dim>
    double Magnetoviscoelastic_Constitutive_Law<
      dim>::get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
    {
      return std::pow(1.0 / std::cosh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat)),
                      2.0);
    };


    template <int dim>
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<
      dim>::get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const
    {
      return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) * get_H();
    };


    template <int dim>
    Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_df_mu_dH(
      const double mu,
      const double mu_inf,
      const double mu_h_sat) const
    {
      return (mu_inf / mu - 1.0) *
             (get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
              get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat));
    };


    template <int dim>
    double Magnetoviscoelastic_Constitutive_Law<
      dim>::get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
    {
      return -2.0 * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
             get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat);
    };


    template <int dim>
    SymmetricTensor<2, dim> Magnetoviscoelastic_Constitutive_Law<
      dim>::get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const
    {
      return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) *
             Physics::Elasticity::StandardTensors<dim>::I;
    };


    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2f_mu_dH_dH(
      const double mu,
      const double mu_inf,
      const double mu_h_sat) const
    {
      return (mu_inf / mu - 1.0) *
             (get_d2tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
                symmetrize(
                  outer_product(get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat),
                                get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat))) +
              get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
                get_d2two_h_dot_h_div_h_sat_squ_dH_dH(mu_h_sat));
    };

    // For the cached calculation approach that we've adopted for this material
    // class, the root of all calculations are the field variables, and the
    // immutable ancillary data such as the constitutive parameters and time
    // step size. As such, we need to enter them into the cache in a different
    // manner to the other variables, since they are inputs that are prescribed
    // from outside the class itself. This function simply adds them to the
    // cache directly from the input arguments, checking that there is no
    // equivalent data there in the first place (we expect to call the
    // `update_internal_data()` method only once per time step, or Newton
    // iteration).
    template <int dim>
    void Magnetoviscoelastic_Constitutive_Law<dim>::set_primary_variables(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H) const
    {
      // Set value for $\boldsymbol{\mathbb{H}}$.
      const std::string name_H("H");
      Assert(!cache.stores_object_with_name(name_H),
             ExcMessage(
               "The primary variable has already been added to the cache."));
      cache.add_unique_copy(name_H, H);

      // Set value for $\mathbf{C}$.
      const std::string name_C("C");
      Assert(!cache.stores_object_with_name(name_C),
             ExcMessage(
               "The primary variable has already been added to the cache."));
      cache.add_unique_copy(name_C, C);
    }

    // After that, we can fetch them from the cache at any point in time.
    template <int dim>
    const Tensor<1, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_H() const
    {
      const std::string name("H");
      Assert(cache.stores_object_with_name(name),
             ExcMessage("Primary variables must be added to the cache."));
      return cache.template get_object_with_name<Tensor<1, dim>>(name);
    }

    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_C() const
    {
      const std::string name("C");
      Assert(cache.stores_object_with_name(name),
             ExcMessage("Primary variables must be added to the cache."));
      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }

    // With the primary variables guaranteed to be in the cache when we need
    // them, we can not compute all intermediate values (either directly, or
    // indirectly) from them.
    //
    // If the cache does not already store the value that we're looking for,
    // then we quickly calculate it, store it in the cache and return the value
    // just stored in the cache. That way we can return it as a reference and
    // avoid copying the object. The same goes for any values that a compound
    // function might depend on. Said another way, if there is a dependency
    // chain of calculations that come before the one that we're currently
    // interested in doing, then we're guaranteed to resolve the dependencies
    // before we proceed with using any of those values. Although there is a
    // cost to fetching data from the cache, the "resolved dependency" concept
    // might be sufficiently convenient to make it worth looking past the extra
    // cost. If these material laws are embedded within a finite element
    // framework, then the added cost might not even be noticeable.
    template <int dim>
    const double &Magnetoviscoelastic_Constitutive_Law<dim>::get_det_F() const
    {
      const std::string name("det_F");
      if (cache.stores_object_with_name(name) == false)
        {
          const double det_F = std::sqrt(determinant(get_C()));
          AssertThrow(det_F > 0.0,
                      ExcMessage("Volumetric Jacobian must be positive."));
          cache.add_unique_copy(name, det_F);
        }

      return cache.template get_object_with_name<double>(name);
    }


    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv() const
    {
      const std::string name("C_inv");
      if (cache.stores_object_with_name(name) == false)
        {
          cache.add_unique_copy(name, invert(get_C()));
        }

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }


    template <int dim>
    const double &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_log_det_F() const
    {
      const std::string name("log(det_F)");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, std::log(get_det_F()));

      return cache.template get_object_with_name<double>(name);
    }


    template <int dim>
    const double &Magnetoviscoelastic_Constitutive_Law<dim>::get_trace_C() const
    {
      const std::string name("trace(C)");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, trace(get_C()));

      return cache.template get_object_with_name<double>(name);
    }


    template <int dim>
    const Tensor<1, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv_dot_H() const
    {
      const std::string name("C_inv_dot_H");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, get_C_inv() * get_H());

      return cache.template get_object_with_name<Tensor<1, dim>>(name);
    }


    template <int dim>
    const double &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_H_dot_C_inv_dot_H() const
    {
      const std::string name("H_dot_C_inv_dot_H");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, get_H() * get_C_inv_dot_H());

      return cache.template get_object_with_name<double>(name);
    }


    template <int dim>
    const SymmetricTensor<4, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dQ_t_dC(
      const DiscreteTime &time) const
    {
      const std::string name("dQ_t_dC");
      if (cache.stores_object_with_name(name) == false)
        {
          const double  delta_t = this->get_delta_t(time);
          const double &det_F   = get_det_F();

          const SymmetricTensor<4, dim> dQ_t_dC =
            (1.0 / (1.0 + delta_t / this->get_tau_v())) *
            (delta_t / this->get_tau_v()) *
            ((2.0 / dim) * std::pow(det_F, 2.0 / dim - 1.0) *
               outer_product(get_C_inv(), get_ddet_F_dC()) +
             std::pow(det_F, 2.0 / dim) * get_dC_inv_dC());

          cache.add_unique_copy(name, dQ_t_dC);
        }

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<4, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dC_inv_dC() const
    {
      const std::string name("dC_inv_dC");
      if (cache.stores_object_with_name(name) == false)
        {
          const SymmetricTensor<2, dim> &C_inv = get_C_inv();
          SymmetricTensor<4, dim>        dC_inv_dC;

          for (unsigned int A = 0; A < dim; ++A)
            for (unsigned int B = A; B < dim; ++B)
              for (unsigned int C = 0; C < dim; ++C)
                for (unsigned int D = C; D < dim; ++D)
                  dC_inv_dC[A][B][C][D] -=               //
                    0.5 * (C_inv[A][C] * C_inv[B][D]     //
                           + C_inv[A][D] * C_inv[B][C]); //

          cache.add_unique_copy(name, dC_inv_dC);
        }

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d_tr_C_dC() const
    {
      const std::string name("d_tr_C_dC");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name,
                              Physics::Elasticity::StandardTensors<dim>::I);

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_ddet_F_dC() const
    {
      const std::string name("ddet_F_dC");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, 0.5 * get_det_F() * get_C_inv());

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dlog_det_F_dC() const
    {
      const std::string name("dlog_det_F_dC");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, 0.5 * get_C_inv());

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }


    template <int dim>
    const Tensor<1, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dH() const
    {
      const std::string name("dH_dot_C_inv_dot_H_dH");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, 2.0 * get_C_inv_dot_H());

      return cache.template get_object_with_name<Tensor<1, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dC() const
    {
      const std::string name("dH_dot_C_inv_dot_H_dC");
      if (cache.stores_object_with_name(name) == false)
        {
          const Tensor<1, dim> C_inv_dot_H = get_C_inv_dot_H();
          cache.add_unique_copy(
            name, -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H)));
        }

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<4, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2log_det_F_dC_dC() const
    {
      const std::string name("d2log_det_F_dC_dC");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, 0.5 * get_dC_inv_dC());

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<4, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2det_F_dC_dC() const
    {
      const std::string name("d2det_F_dC_dC");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name,
                              0.5 *
                                (outer_product(get_C_inv(), get_ddet_F_dC()) +
                                 get_det_F() * get_dC_inv_dC()));

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<2, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dH_dH()
      const
    {
      const std::string name("d2H_dot_C_inv_dot_H_dH_dH");
      if (cache.stores_object_with_name(name) == false)
        cache.add_unique_copy(name, 2.0 * get_C_inv());

      return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
    }


    template <int dim>
    const Tensor<3, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dH()
      const
    {
      const std::string name("d2H_dot_C_inv_dot_H_dC_dH");
      if (cache.stores_object_with_name(name) == false)
        {
          const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H();
          const SymmetricTensor<2, dim> &C_inv       = get_C_inv();

          Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH;
          for (unsigned int A = 0; A < dim; ++A)
            for (unsigned int B = 0; B < dim; ++B)
              for (unsigned int C = 0; C < dim; ++C)
                d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -=
                  C_inv[A][C] * C_inv_dot_H[B] + //
                  C_inv_dot_H[A] * C_inv[B][C];  //

          cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dH);
        }

      return cache.template get_object_with_name<Tensor<3, dim>>(name);
    }


    template <int dim>
    const SymmetricTensor<4, dim> &
    Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dC()
      const
    {
      const std::string name("d2H_dot_C_inv_dot_H_dC_dC");
      if (cache.stores_object_with_name(name) == false)
        {
          const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H();
          const SymmetricTensor<2, dim> &C_inv       = get_C_inv();

          SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC;
          for (unsigned int A = 0; A < dim; ++A)
            for (unsigned int B = A; B < dim; ++B)
              for (unsigned int C = 0; C < dim; ++C)
                for (unsigned int D = C; D < dim; ++D)
                  d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] +=
                    0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] +
                           C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] +
                           C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] +
                           C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]);

          cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dC);
        }

      return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
    }

    // @sect4{Rheological experiment parameters}

    // The @p RheologicalExperimentParameters class is used to drive the
    // numerical experiments that are to be conducted on the coupled materials
    // that we've implemented constitutive laws for.
    class RheologicalExperimentParameters : public ParameterAcceptor
    {
    public:
      RheologicalExperimentParameters();

      // These are  dimensions of the rheological specimen that is to be
      // simulated. They, effectively, define the measurement point for our
      // virtual experiment.
      double sample_radius = 0.01;
      double sample_height = 0.001;

      // The three steady-state loading parameters are respectively
      // - the axial stretch,
      // - the shear strain amplitude, and
      // - the axial magnetic field strength.
      double lambda_2 = 0.95;
      double gamma_12 = 0.05;
      double H_2      = 60.0e3;

      // Moreover, the parameters for the time-dependent rheological
      // loading conditions are
      // - the loading cycle frequency,
      // - the number of load cycles, and
      // - the number of discrete timesteps per cycle.
      double       frequency         = 1.0 / (2.0 * numbers::PI);
      unsigned int n_cycles          = 5;
      unsigned int n_steps_per_cycle = 2500;

      // We also declare some self-explanatory parameters related to output
      // data generated for the experiments conducted with rate-dependent and
      // rate-independent materials.
      bool        output_data_to_file = true;
      std::string output_filename_rd =
        "experimental_results-rate_dependent.csv";
      std::string output_filename_ri =
        "experimental_results-rate_independent.csv";

      // The next few functions compute time-related parameters for the
      // experiment...
      double start_time() const;

      double end_time() const;

      double delta_t() const;

      // ... while the following two prescribe the mechanical and magnetic
      // loading at any given time...
      Tensor<1, 3> get_H(const double time) const;

      Tensor<2, 3> get_F(const double time) const;

      // ... and this last one outputs the status of the experiment to the
      // console.
      bool print_status(const int step_number) const;

      bool initialized = false;
    };



    RheologicalExperimentParameters::RheologicalExperimentParameters()
      : ParameterAcceptor("/Coupled Constitutive Laws/Rheological Experiment/")
    {
      add_parameter("Experimental sample radius", sample_radius);
      add_parameter("Experimental sample radius", sample_height);

      add_parameter("Axial stretch", lambda_2);
      add_parameter("Shear strain amplitude", gamma_12);
      add_parameter("Axial magnetic field strength", H_2);

      add_parameter("Frequency", frequency);
      add_parameter("Number of loading cycles", n_cycles);
      add_parameter("Discretisation for each cycle", n_steps_per_cycle);

      add_parameter("Output experimental results to file", output_data_to_file);
      add_parameter("Output file name (rate dependent constitutive law)",
                    output_filename_rd);
      add_parameter("Output file name (rate independent constitutive law)",
                    output_filename_ri);

      parse_parameters_call_back.connect([&]() -> void { initialized = true; });
    }


    double RheologicalExperimentParameters::start_time() const
    {
      return 0.0;
    }


    double RheologicalExperimentParameters::end_time() const
    {
      return n_cycles / frequency;
    }


    double RheologicalExperimentParameters::delta_t() const
    {
      return (end_time() - start_time()) / (n_steps_per_cycle * n_cycles);
    }


    bool
    RheologicalExperimentParameters::print_status(const int step_number) const
    {
      return (step_number % (n_cycles * n_steps_per_cycle / 100)) == 0;
    }

    // The applied magnetic field is always aligned with the axis of rotation
    // of the rheometer's rotor.
    Tensor<1, 3> RheologicalExperimentParameters::get_H(const double) const
    {
      return Tensor<1, 3>({0.0, 0.0, H_2});
    }

    // The applied deformation (gradient) is computed based on the geometry
    // of the rheometer and the sample, the sampling point, and the experimental
    // parameters. From the displacement profile documented in the introduction,
    // the deformation gradient may be expressed in Cartesian coordinates as
    // @f[
    // \mathbf{F} = \begin{bmatrix}
    //    \frac{\cos\left(\alpha\right)}{\sqrt{\lambda_{3}}}
    // & -\frac{\sin\left(\alpha\right)}{\sqrt{\lambda_{3}}}
    // & -\tau R \sqrt{\lambda_{3}} \sin\left(\Theta + \alpha\right)
    // \\  \frac{\sin\left(\alpha\right)}{\sqrt{\lambda_{3}}}
    // & \frac{\cos\left(\alpha\right)}{\sqrt{\lambda_{3}}}
    // & -\tau R \sqrt{\lambda_{3}} \cos\left(\Theta + \alpha\right)
    // \\  0 & 0 & \lambda_{3}
    // \end{bmatrix}
    // @f]
    Tensor<2, 3> RheologicalExperimentParameters::get_F(const double time) const
    {
      AssertThrow((sample_radius > 0.0 && sample_height > 0.0),
                  ExcMessage("Non-physical sample dimensions"));
      AssertThrow(lambda_2 > 0.0,
                  ExcMessage("Non-physical applied axial stretch"));

      const double sqrt_lambda_2     = std::sqrt(lambda_2);
      const double inv_sqrt_lambda_2 = 1.0 / sqrt_lambda_2;

      const double alpha_max =
        std::atan(std::tan(gamma_12) * sample_height /
                  sample_radius); // Small strain approximation
      const double A       = sample_radius * alpha_max;
      const double w       = 2.0 * numbers::PI * frequency; // in rad /s
      const double gamma_t = A * std::sin(w * time);
      const double tau_t =
        gamma_t /
        (sample_radius * sample_height); // Torsion angle per unit length
      const double alpha_t = tau_t * lambda_2 * sample_height;

      Tensor<2, 3> F;
      F[0][0] = inv_sqrt_lambda_2 * std::cos(alpha_t);
      F[0][1] = -inv_sqrt_lambda_2 * std::sin(alpha_t);
      F[0][2] = -tau_t * sample_radius * sqrt_lambda_2 * std::sin(alpha_t);
      F[1][0] = inv_sqrt_lambda_2 * std::sin(alpha_t);
      F[1][1] = inv_sqrt_lambda_2 * std::cos(alpha_t);
      F[1][2] = tau_t * sample_radius * sqrt_lambda_2 * std::cos(alpha_t);
      F[2][0] = 0.0;
      F[2][1] = 0.0;
      F[2][2] = lambda_2;

      AssertThrow((F[0][0] > 0) && (F[1][1] > 0) && (F[2][2] > 0),
                  ExcMessage("Non-physical deformation gradient component."));
      AssertThrow(std::abs(determinant(F) - 1.0) < 1e-6,
                  ExcMessage("Volumetric Jacobian is not equal to unity."));

      return F;
    }

    // @sect4{Rheological experiment: Parallel plate rotational rheometer}

    // This is the function that will drive the numerical experiments.
    template <int dim>
    void run_rheological_experiment(
      const RheologicalExperimentParameters &experimental_parameters,
      Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
        &material_hand_calculated,
      Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
        &               material_assisted_computation,
      TimerOutput &     timer,
      const std::string filename)
    {
      // We can take the hand-implemented constitutive law and compare the
      // results that we attain with it to those that we get using AD or SD.
      // In this way, we can verify that they produce identical results (which
      // indicates that either both implementations have a high probability of
      // being correct, or that they're incorrect with identical flaws being
      // present in both). Either way, it is a decent sanity check for the
      // fully self-implemented variants and can certainly be used as a
      // debugging strategy when differences between the results are
      // detected).
      const auto check_material_class_results =
        [](
          const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &to_verify,
          const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &blessed,
          const double tol = 1e-6) {
          (void)to_verify;
          (void)blessed;
          (void)tol;

          Assert(std::abs(blessed.get_psi() - to_verify.get_psi()) < tol,
                 ExcMessage("No match for psi. Error: " +
                            Utilities::to_string(std::abs(
                              blessed.get_psi() - to_verify.get_psi()))));

          Assert((blessed.get_B() - to_verify.get_B()).norm() < tol,
                 ExcMessage("No match for B. Error: " +
                            Utilities::to_string(
                              (blessed.get_B() - to_verify.get_B()).norm())));
          Assert((blessed.get_S() - to_verify.get_S()).norm() < tol,
                 ExcMessage("No match for S. Error: " +
                            Utilities::to_string(
                              (blessed.get_S() - to_verify.get_S()).norm())));

          Assert((blessed.get_DD() - to_verify.get_DD()).norm() < tol,
                 ExcMessage("No match for BB. Error: " +
                            Utilities::to_string(
                              (blessed.get_DD() - to_verify.get_DD()).norm())));
          Assert((blessed.get_PP() - to_verify.get_PP()).norm() < tol,
                 ExcMessage("No match for PP. Error: " +
                            Utilities::to_string(
                              (blessed.get_PP() - to_verify.get_PP()).norm())));
          Assert((blessed.get_HH() - to_verify.get_HH()).norm() < tol,
                 ExcMessage("No match for HH. Error: " +
                            Utilities::to_string(
                              (blessed.get_HH() - to_verify.get_HH()).norm())));
        };

      // We'll be outputting the constitutive response of the material to file
      // for post-processing, so here we declare a `stream` that will act as
      // a buffer for this output. We'll use a simple CSV format for the
      // outputted results.
      std::ostringstream stream;
      stream
        << "Time;Axial magnetic field strength [A/m];Axial magnetic induction [T];Shear strain [%];Shear stress [Pa]\n";

      // Using the DiscreteTime class, we iterate through each timestep using
      // a fixed time step size.
      for (DiscreteTime time(experimental_parameters.start_time(),
                             experimental_parameters.end_time() +
                               experimental_parameters.delta_t(),
                             experimental_parameters.delta_t());
           time.is_at_end() == false;
           time.advance_time())
        {
          if (experimental_parameters.print_status(time.get_step_number()))
            std::cout << "Timestep = " << time.get_step_number()
                      << " @ time = " << time.get_current_time() << "s."
                      << std::endl;

          // We fetch and compute the loading to be applied to the material
          // at this time step...
          const Tensor<1, dim> H =
            experimental_parameters.get_H(time.get_current_time());
          const Tensor<2, dim> F =
            experimental_parameters.get_F(time.get_current_time());
          const SymmetricTensor<2, dim> C =
            Physics::Elasticity::Kinematics::C(F);

          // ... then we update the state of the materials...
          {
            TimerOutput::Scope timer_section(timer, "Hand calculated");
            material_hand_calculated.update_internal_data(C, H, time);
            material_hand_calculated.update_end_of_timestep();
          }

          {
            TimerOutput::Scope timer_section(timer, "Assisted computation");
            material_assisted_computation.update_internal_data(C, H, time);
            material_assisted_computation.update_end_of_timestep();
          }

          // ... and test for discrepencies between the two.
          check_material_class_results(material_hand_calculated,
                                       material_assisted_computation);

          if (experimental_parameters.output_data_to_file)
            {
              // The next thing that we will do is collect some results to
              // post-process. All quantities are in the "current configuration"
              // (rather than the "reference configuration", in which all
              // quantities computed by the constitutive laws are framed).
              const Tensor<1, dim> h =
                Physics::Transformations::Covariant::push_forward(H, F);
              const Tensor<1, dim> b =
                Physics::Transformations::Piola::push_forward(
                  material_hand_calculated.get_B(), F);
              const SymmetricTensor<2, dim> sigma =
                Physics::Transformations::Piola::push_forward(
                  material_hand_calculated.get_S(), F);
              stream << time.get_current_time() << ";" << h[2] << ";" << b[2]
                     << ";" << F[1][2] * 100.0 << ";" << sigma[1][2] << "\n";
            }
        }

      // Finally, we output the strain-stress and magnetic loading history to
      // file.
      if (experimental_parameters.output_data_to_file)
        {
          std::ofstream output(filename);
          output << stream.str();
        }
    };

    // @sect4{The CoupledConstitutiveLaws::run() function}

    // The purpose of this driver function is to read in all of the parameters
    // from file and, based off of that, create a representative instance of
    // each constitutive law and invoke the function that conducts a rheological
    // experiment with it.
    void run(int argc, char *argv[])
    {
      using namespace dealii;

      constexpr unsigned int dim = 3;

      const ConstitutiveParameters          constitutive_parameters;
      const RheologicalExperimentParameters experimental_parameters;

      std::string parameter_file;
      if (argc > 1)
        parameter_file = argv[1];
      else
        parameter_file = "parameters.prm";
      ParameterAcceptor::initialize(parameter_file, "used_parameters.prm");

      // We start the actual work by configuring and running the experiment
      // using our rate-independent constitutive law. The automatically
      // differentiable number type is hard-coded here, but with some clever
      // templating it is possible to select which framework to use at run time
      // (e.g., as selected through the parameter file). We'll simultaneously
      // perform the experiments with the counterpary material law that was
      // fully implemented by hand, and check what it computes against our
      // assisted implementation.
      {
        TimerOutput timer(std::cout,
                          TimerOutput::summary,
                          TimerOutput::wall_times);
        std::cout
          << "Coupled magnetoelastic constitutive law using automatic differentiation."
          << std::endl;

        constexpr Differentiation::AD::NumberTypes ADTypeCode =
          Differentiation::AD::NumberTypes::sacado_dfad_dfad;

        Magnetoelastic_Constitutive_Law<dim> material(constitutive_parameters);
        Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode> material_ad(
          constitutive_parameters);

        run_rheological_experiment(experimental_parameters,
                                   material,
                                   material_ad,
                                   timer,
                                   experimental_parameters.output_filename_ri);

        std::cout << "... all calculations are correct!" << std::endl;
      }

      // Next we do the same for the rate-dependent constitutive law.
      // The highest performance option is selected as default if SymEngine
      // is set up to use the LLVM just-in-time compiler which (in conjunction
      // with some aggressive compilation flags) produces the fastest code
      // evaluation path of all of the available option. As a fall-back, the
      // so called "lambda" optimizer (which only requires a C++11 compliant
      // compiler) will be selected. At the same time, we'll ask the CAS to
      // perform common subexpression elimination to minimize the number of
      // intermediate calculations used during evaluation.
      // We'll record how long it takes to execute the "initialization" step
      // inside the constructor for the SD implementation, as this is where the
      // abovementioned transformations occur.
      {
        TimerOutput timer(std::cout,
                          TimerOutput::summary,
                          TimerOutput::wall_times);
        std::cout
          << "Coupled magneto-viscoelastic constitutive law using symbolic differentiation."
          << std::endl;

#ifdef DEAL_II_SYMENGINE_WITH_LLVM
        std::cout << "Using LLVM optimizer." << std::endl;
        constexpr Differentiation::SD::OptimizerType optimizer_type =
          Differentiation::SD::OptimizerType::llvm;
        constexpr Differentiation::SD::OptimizationFlags optimization_flags =
          Differentiation::SD::OptimizationFlags::optimize_all;
#else
        std::cout << "Using lambda optimizer." << std::endl;
        constexpr Differentiation::SD::OptimizerType optimizer_type =
          Differentiation::SD::OptimizerType::lambda;
        constexpr Differentiation::SD::OptimizationFlags optimization_flags =
          Differentiation::SD::OptimizationFlags::optimize_cse;
#endif

        Magnetoviscoelastic_Constitutive_Law<dim> material(
          constitutive_parameters);

        timer.enter_subsection("Initialize symbolic CL");
        Magnetoviscoelastic_Constitutive_Law_SD<dim> material_sd(
          constitutive_parameters, optimizer_type, optimization_flags);
        timer.leave_subsection();

        run_rheological_experiment(experimental_parameters,
                                   material,
                                   material_sd,
                                   timer,
                                   experimental_parameters.output_filename_rd);

        std::cout << "... all calculations are correct!" << std::endl;
      }
    }

  } // namespace CoupledConstitutiveLaws

} // namespace Step71


// @sect3{The main() function}

// The main function only calls the driver functions for the two sets of
// examples that are to be executed.
int main(int argc, char *argv[])
{
  Step71::SimpleExample::run();
  Step71::CoupledConstitutiveLaws::run(argc, argv);

  return 0;
}
