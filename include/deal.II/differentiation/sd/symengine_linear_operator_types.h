// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_differentiation_sd_symengine_linear_operator_types_h
#define dealii_differentiation_sd_symengine_linear_operator_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
// Low level
#  include <symengine/basic.h>
#  include <symengine/symengine_exception.h>
#  include <symengine/symengine_rcp.h>

// Number types
#  include <symengine/complex_double.h>
#  include <symengine/integer.h>
#  include <symengine/number.h>
#  include <symengine/rational.h>
#  include <symengine/real_double.h>

// Number operations
#  include <symengine/add.h>
#  include <symengine/mul.h>

// Optimisation
#  include <symengine/lambda_double.h>
#  include <symengine/visitor.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/numbers.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_la_vector_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <deal.II/lac/linear_operator_tools.h>

#  include <algorithm>
#  include <complex>
#  include <iterator>
#  include <type_traits>
#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;


    // Forward declarations
    template <typename Range, typename Domain, typename Payload>
    class LinearOperatorWrapper;


    namespace internal
    {
      template <typename Range, typename Domain, typename Payload>
      static SE::RCP<const SE::Basic>
      make_symengine_rcp(const LinearOperator<Range, Domain, Payload> &op,
                         const std::string &                           name);

      template <typename Range, typename Domain, typename Payload>
      const LinearOperator<Range, Domain, Payload> &
      evaluate_symengine_number(
        const LinearOperatorWrapper<Range, Domain, Payload> &op);

      template <typename Range, typename Domain, typename Payload>
      bool
      is_a_LinearOperatorWrapper(const SE::RCP<const SE::Basic> &value);
    } // namespace internal



    /**
     * A number class to hold sparse matrices wrapped as LinearOperators.
     *
     * @author Jean-Paul Pelteret, 2018
     */
    template <typename Range_, typename Domain_, typename Payload_>
    class LinearOperatorWrapper : public SE::NumberWrapper, public SE::Evaluate
    {
    public:
      typedef Range_                                    Range;
      typedef Domain_                                   Domain;
      typedef Payload_                                  Payload;
      typedef LinearOperator<Range_, Domain_, Payload_> LinearOperatorType;


      /**
       * Class constructor.
       *
       * This is marked as explicit so that there are no accidental conversions
       * of compatible numbers to this class type.
       */
      explicit LinearOperatorWrapper(const LinearOperatorType &op,
                                     const std::string &       name);

      /**
       * Class destructor.
       */
      virtual ~LinearOperatorWrapper();

      /**
       * @name Value retrieval
       */
      //@{

      /**
       * Returns a copy of the @p value of the underlying linear operator.
       */
      inline LinearOperatorType
      get_value()
      {
        return this->value;
      }

      /**
       * Returns the underlying linear operator by reference.
       */
      inline const LinearOperatorType &
      get_value() const
      {
        return this->value;
      }

      /**
       * Returns the name associated with the underlying linear operator.
       */
      inline const std::string &
      get_name() const
      {
        return this->name;
      }

      //@}

      /**
       * @name Internally used methods overloaded from the SymEngine classes
       */
      //@{

      // ==== METHODS FROM Basic ====

      /**
       * Return the size of the hash for this class.
       *
       * This overloads the same function from the SymEngine::Basic class.
       */
      virtual SE::hash_t
      __hash__() const;

      /**
       * An equality comparator that returns <tt>true</tt> if the @p other class
       * instance takes the same value as this one.
       *
       * This overloads the same function from the SymEngine::Basic class.
       */
      virtual bool
      __eq__(const SE::Basic &other) const;

      /**
       * A comparator that returns <tt>0</tt> if the @p other class
       * instance takes the same value as this one, <tt>-1</tt> if its value
       * is less than that of the other class, otherwise <tt>+1</tt>.
       *
       * This overloads the same function from the SymEngine::Basic class.
       */
      virtual int
      compare(const SE::Basic &other) const;

      // ==== METHODS FROM Number ====

      /**
       * Nominally returns <tt>true</tt> if the @p value stored in this class
       * instance is an exact number. Since auto-differentiable number are
       * floating point values, this function always returns <tt>false</tt>.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual bool
      is_exact() const;

      /**
       * Returns <tt>true</tt> if the @p value stored in this class
       * instance is positive.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual bool
      is_positive() const;

      /**
       * Returns <tt>true</tt> if the @p value stored in this class
       * instance is negative.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual bool
      is_negative() const;

      /**
       * Returns <tt>true</tt> if the @p value stored in this class
       * instance is precisely zero.
       *
       * This overloads the same function from the SymEngine::Number class.
       *
       * SymEngine uses this result to filer out and simplify some operations.
       * Since auto-differentiable numbers may be zero-valued but have non-zero
       * partial derivatives, this function is set to always return false
       * thereby correctly keeping track of the sensitivities for subsequent
       * operations.
       */
      virtual bool
      is_zero() const;

      /**
       * Nominally returns <tt>true</tt> if the @p value stored in this class
       * instance is precisely equal to <tt>+1</tt>. Since auto-differentiable
       * numbers are floating point values, this function always returns
       * <tt>false</tt>.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual bool
      is_one() const;

      /**
       * Nominally returns <tt>true</tt> if the @p value stored in this class
       * instance is precisely equal to <tt>-1</tt>. Since auto-differentiable
       * number are floating point values, this function always returns
       * <tt>false</tt>.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual bool
      is_minus_one() const;

      /**
       * Returns <tt>true</tt> if the @p value stored in this class
       * instance of a complex type.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual bool
      is_complex() const;

      /**
       * Returns an object that can numerically evaluate an instance of
       * this class.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::Evaluate &
      get_eval() const;

      // ==== METHODS FROM NumberWrapper ====

      /**
       * Returns a string that expresses the value of the stored
       * auto-differentiable number.
       *
       * This overloads the same function from the SymEngine::NumberWrapper
       * class.
       */
      virtual std::string
      __str__() const;

      /**
       * Returns an SymEngine pointer to a number representing the numerical
       * value of the underlying auto-differentiable number. Since this does
       * not make much sense in this application, calling this function results
       * in an error being thrown.
       *
       * This overloads the same function from the SymEngine::NumberWrapper
       * class.
       */
      virtual SE::RCP<const SE::Number>
      eval(long bits) const;

      //@}

      /**
       * @name Basic mathematical operations
       */
      //@{

      /**
       * Computes the addition of another SymEngine::Number to this class
       * instance's @p value .
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      add(const SE::Number &other) const;

      /**
       * Computes the subtraction of another SymEngine::Number from this class
       * instance's @p value .
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      sub(const SE::Number &other) const;

      /**
       * Computes the subtraction of this class instance's @p value from that of
       * another SymEngine::Number (i.e. "reciprocal subtraction").
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      rsub(const SE::Number &other) const;

      /**
       * Computes the multiplication of another SymEngine::Number with this
       * class
       * instance's @p value .
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      mul(const SE::Number &other) const;

      /**
       * Computes the division of this class instance's @p value by that of
       * another SymEngine::Number.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      div(const SE::Number &other) const;

      /**
       * Computes the division of another SymEngine::Number by this class
       * instance's @p value  (i.e. "reciprocal division").
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      rdiv(const SE::Number &other) const;

      //@}

      /**
       * @name Power functions
       */
      //@{

      /**
       * Computes the result of this class instance's @p value raised to the
       * power of another SymEngine::Number.
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      pow(const SE::Number &other) const;

      /**
       * Computes the result of another SymEngine::Number raised to the power
       * of this class instance's @p value (i.e. the "reciprocal power").
       *
       * This overloads the same function from the SymEngine::Number class.
       */
      virtual SE::RCP<const SE::Number>
      rpow(const SE::Number &other) const;

      /**
       * Computes the result of the natural number raised to the power
       * of the ADNumberWrapper @p x  (i.e. the exponential function of @p x ).
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      exp(const SE::Basic &x) const;

      /**
       * Computes the natural number logarithm of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      log(const SE::Basic &x) const;

      //@}

      /**
       * @name Trignometric functions
       */
      //@{

      /**
       * Computes the sine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      sin(const SE::Basic &x) const;

      /**
       * Computes the cosine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      cos(const SE::Basic &x) const;

      /**
       * Computes the tangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      tan(const SE::Basic &x) const;

      /**
       * Computes the cosecant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      csc(const SE::Basic &x) const;

      /**
       * Computes the secant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      sec(const SE::Basic &x) const;

      /**
       * Computes the cotangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      cot(const SE::Basic &x) const;

      /**
       * Computes the inverse sine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      asin(const SE::Basic &x) const;

      /**
       * Computes the inverse cosine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      acos(const SE::Basic &x) const;

      /**
       * Computes the inverse tangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      atan(const SE::Basic &x) const;

      /**
       * Computes the inverse cosecant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      acsc(const SE::Basic &x) const;

      /**
       * Computes the inverse secant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      asec(const SE::Basic &x) const;

      /**
       * Computes the inverse cotangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      acot(const SE::Basic &x) const;

      //@}

      /**
       * @name Hyperbolic trignometric functions
       */
      //@{

      /**
       * Computes the hyperbolic sine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      sinh(const SE::Basic &x) const;

      /**
       * Computes the hyperbolic cosine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      cosh(const SE::Basic &x) const;

      /**
       * Computes the hyperbolic tangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      tanh(const SE::Basic &x) const;

      /**
       * Computes the hyperbolic cosecant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      csch(const SE::Basic &x) const;

      /**
       * Computes the hyperbolic secant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      sech(const SE::Basic &x) const;

      /**
       * Computes the hyperbolic cotangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      coth(const SE::Basic &x) const;

      /**
       * Computes the inverse hyperbolic sine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      asinh(const SE::Basic &x) const;

      /**
       * Computes the inverse hyperbolic cosine of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      acosh(const SE::Basic &x) const;

      /**
       * Computes the inverse hyperbolic tangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      atanh(const SE::Basic &x) const;

      /**
       * Computes the inverse hyperbolic cosecant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      acsch(const SE::Basic &x) const;

      /**
       * Computes the inverse hyperbolic secant of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      asech(const SE::Basic &x) const;

      /**
       * Computes the inverse hyperbolic cotangent of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      acoth(const SE::Basic &x) const;

      //@}

      /**
       * @name Other functions
       */
      //@{


      /**
       * Computes the absolute value of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      abs(const SE::Basic &x) const;

      /**
       * Computes the rounded-down value of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      floor(const Basic &) const;

      /**
       * Computes the rounded-up value of the ADNumberWrapper @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      ceiling(const Basic &) const;

      /**
       * Computes the error function with the ADNumberWrapper argument @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      erf(const SE::Basic &) const;

      /**
       * Computes the complementary error function with the ADNumberWrapper
       * argument @p x .
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const Basic>
      erfc(const SE::Basic &) const;

      /**
       * Computes the gamma function of the ADNumberWrapper @p x .
       * This function is not implemented, and will throw an error if called.
       *
       * This overloads the same function from the SymEngine::Evaluate class.
       */
      virtual SE::RCP<const SE::Basic>
      gamma(const SE::Basic &x) const;

      //@}

    private:
      /**
       * The value of the LinearOperator that is being stored.
       */
      LinearOperatorType value;

      /**
       * A symbolic name to define the linear operation that this
       * symbol represents
       */
      std::string name;

    }; // class LinearOperatorWrapper


    /* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN


    namespace internal
    {
      template <typename Range, typename Domain, typename Payload>
      inline SE::RCP<const SE::Basic>
      make_symengine_rcp(const LinearOperator<Range, Domain, Payload> &op,
                         const std::string &                           name)
      {
        return SE::make_rcp<
          const LinearOperatorWrapper<Range, Domain, Payload>>(op, name);
      }


      template <typename Range, typename Domain, typename Payload>
      const LinearOperatorWrapper<Range, Domain, Payload> &
      evaluate_symengine_number(
        const LinearOperatorWrapper<Range, Domain, Payload> &op)
      {
        return op.get_value();
      }



      template <typename Range, typename Domain, typename Payload>
      bool
      is_a_LinearOperatorWrapper(const SE::Basic &value)
      {
        return dynamic_cast<
          const LinearOperatorWrapper<Range, Domain, Payload> *const>(&value);
      }


      template <typename Range, typename Domain, typename Payload>
      struct LinearOperatorEvaluator
      {
        static LinearOperator<Range, Domain, Payload>
        value(const SE::RCP<const SE::Basic> &value)
        {
          if (internal::is_a_LinearOperatorWrapper<Range, Domain, Payload>(
                *value))
            {
              return SE::down_cast<
                       const LinearOperator<Range, Domain, Payload> &>(*value)
                .get_value();
            }
          else
            {
              AssertThrow(false, ExcNotImplemented());
              return LinearOperator<Range, Domain, Payload>();
            }
        }
      };

    } // namespace internal


    /**
     * An error declaration stating that the object
     * is not an LinearOperatorWrapper type.
     */
    DeclExceptionMsg(
      ExcNotLinearOperatorWrapper,
      "Object not of type LinearOperatorWrapper<Range,Domain,Payload>.");

    /**
     * An error declaration stating that the object
     * is not an LinearOperatorWrapper type.
     */
    DeclException1(
      ExcOpNotDefinedByLinearOperatorWrapper,
      const SE::Basic *,
      "Operation is not defined by the type LinearOperatorWrapper<Range,Domain,Payload>: " +
        arg1->__str__());


    template <typename Range, typename Domain, typename Payload>
    LinearOperatorWrapper<Range, Domain, Payload>::LinearOperatorWrapper(
      const LinearOperator<Range, Domain, Payload> &op,
      const std::string &                           name)
    {
      this->value = op;
      this->name  = name;
    }


    template <typename Range, typename Domain, typename Payload>
    LinearOperatorWrapper<Range, Domain, Payload>::~LinearOperatorWrapper()
    {}


    template <typename Range, typename Domain, typename Payload>
    inline SE::hash_t
    LinearOperatorWrapper<Range, Domain, Payload>::__hash__() const
    {
      return std::hash<std::string>{}(__str__());
    }


    template <typename Range, typename Domain, typename Payload>
    bool
    LinearOperatorWrapper<Range, Domain, Payload>::__eq__(
      const SE::Basic &other) const
    {
      if (internal::is_a_LinearOperatorWrapper<Range, Domain, Payload>(other))
        {
          // With linear operators, we cannot make direct comparison of
          // any values, so instead we assume that we can do something
          // sensible with their memory addresses.
          const LinearOperatorWrapper<Range, Domain, Payload> &op_other =
            SE::down_cast<
              const LinearOperatorWrapper<Range, Domain, Payload> &>(other);
          return &(this->value) == &(op_other.value);
        }
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    int
    LinearOperatorWrapper<Range, Domain, Payload>::compare(
      const SE::Basic &other) const
    {
      Assert((internal::is_a_LinearOperatorWrapper<Range, Domain, Payload>(
               other)),
             ExcNotLinearOperatorWrapper());

      // With linear operators, we cannot make direct comparison of
      // any values, so instead we assume that we can do something
      // sensible with their memory addresses.
      const LinearOperatorWrapper<Range, Domain, Payload> &op_other =
        SE::down_cast<const LinearOperatorWrapper<Range, Domain, Payload> &>(
          other);
      if (&(this->value) == &(op_other.value))
        return 0;
      return (&(this->value) < &(op_other.value)) ? -1 : 1;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_exact() const
    {
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_positive() const
    {
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_negative() const
    {
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_zero() const
    {
      return this->value.is_null_operator;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_one() const
    {
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_minus_one() const
    {
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    inline bool
    LinearOperatorWrapper<Range, Domain, Payload>::is_complex() const
    {
      return false;
    }


    template <typename Range, typename Domain, typename Payload>
    SE::Evaluate &
    LinearOperatorWrapper<Range, Domain, Payload>::get_eval() const
    {
      // This class derives from SE::Evaluate, so why shouldn't it
      // be able to evaluate itself?
      return const_cast<LinearOperatorWrapper<Range, Domain, Payload> &>(*this);
    }


    template <typename Range, typename Domain, typename Payload>
    std::string
    LinearOperatorWrapper<Range, Domain, Payload>::__str__() const
    {
      std::stringstream ss;
      ss << "LO[";
      ss << this->name;
      ss << "]";
      return ss.str();
    }



    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::eval(long /*bits*/) const
    {
      // Used in the EvalXYZVisitor classes
      // Don't know why we'd need this when we already
      // have the get_eval() virtual function
      AssertThrow(false, ExcNotImplemented());

      // Needed for SymEngine::StrPrinter
      return SE::make_rcp<const LinearOperatorWrapper<Range, Domain, Payload>>(
        get_value(), get_name());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::add(
      const SE::Number &other) const
    {
      if (internal::is_a_LinearOperatorWrapper<Range, Domain, Payload>(other))
        {
          const LinearOperatorWrapper<Range, Domain, Payload> &other_op =
            SE::down_cast<
              const LinearOperatorWrapper<Range, Domain, Payload> &>(other);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            this->value + other_op.value,
            "(" + this->name + "+" + other_op.name + ")");
        }
      // Due to the way SymEngine works, we may be adding some constants when
      // expanding expressions. The following conditions take these into
      // consideration and return some reasonable result.
      else if (SE::is_a<SE::Integer>(other) && other.is_zero())
        {
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(this->value,
                                                                 this->name);
        }

      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return other.add(*this);
    }

    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::sub(
      const SE::Number &other) const
    {
      if (internal::is_a_LinearOperatorWrapper<Range, Domain, Payload>(other))
        {
          const LinearOperatorWrapper<Range, Domain, Payload> &other_op =
            SE::down_cast<
              const LinearOperatorWrapper<Range, Domain, Payload> &>(other);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            this->value - other_op.value,
            "(" + this->name + "-" + other_op.name + ")");
        }
      // Due to the way SymEngine works, we may be adding some constants when
      // expanding expressions. The following conditions take these into
      // consideration and return some reasonable result.
      else if (SE::is_a<SE::Integer>(other) && other.is_zero())
        {
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(this->value,
                                                                 this->name);
        }

      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return other.rsub(*this);
    }

    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::rsub(
      const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::mul(
      const SE::Number &other) const
    {
      if (internal::is_a_LinearOperatorWrapper<Range, Domain, Payload>(other))
        {
          const LinearOperatorWrapper<Range, Domain, Payload> &other_op =
            SE::down_cast<
              const LinearOperatorWrapper<Range, Domain, Payload> &>(other);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            this->value * other_op.value, this->name + "*" + other_op.name);
        }
      else if (other.is_one())
        {
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(this->value,
                                                                 this->name);
        }
      else if (internal::is_a_LAVectorWrapper<Domain>(other))
        {
          const LAVectorWrapper<Domain> &vec_other =
            SE::down_cast<const LAVectorWrapper<Domain> &>(other);
          // TODO: Make this work for PackagedOperations (if this makes sense).
          //       Currently we force evaluation to get the right return type
          //          const Range result = this->value * vec_other.get_value();
          const PackagedOperation<Range> result =
            this->value * vec_other.get_value();
          return SE::make_rcp<const LAVectorWrapper<Range>>(
            result, this->name + "*" + vec_other.get_name());
        }
      else if (SE::is_a<SE::RealDouble>(other))
        {
          const SE::RealDouble &other_nmbr =
            SE::down_cast<const SE::RealDouble &>(other);
          const double value = internal::evaluate_symengine_number(other_nmbr);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            this->value * value,
            this->name + "*" + dealii::Utilities::to_string(value));
        }
      else if (SE::is_a<SE::Complex>(other))
        {
          AssertThrow(false, ExcNotImplemented());
          return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());

          //          const SE::Complex &other_nmbr = SE::down_cast<const
          //          SE::Complex &>(other); return
          //          internal::make_symengine_rcp(this->value *
          //          ADComplexType(other_op));
        }
      else if (SE::is_a<SE::Integer>(other))
        {
          const SE::Integer &other_nmbr =
            SE::down_cast<const SE::Integer &>(other);
          const int value = internal::evaluate_symengine_number(other_nmbr);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            this->value * value,
            this->name + "*" + dealii::Utilities::to_string(value));
        }
      else if (SE::is_a<SE::Rational>(other))
        {
          const SE::Rational &other_nmbr =
            SE::down_cast<const SE::Rational &>(other);
          const double value = internal::evaluate_symengine_number(other_nmbr);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            this->value * value,
            this->name + "*" + dealii::Utilities::to_string(value));
        }

      // We can guarentee that other number types don't know anything about
      // this object, so we can't do anything more...
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return other.mul(*this);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::div(
      const SE::Number &other) const
    {
      if (other.is_one())
        {
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(this->value,
                                                                 this->name);
        }

      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::rdiv(
      const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    namespace internal
    {
      struct LOFixedPower
      {
        template <typename LinearOperatorType>
        static LinearOperatorType
        value(const LinearOperatorType &value, const int exponent)
        {
          Assert(exponent > 0, ExcInternalError());
          if (exponent == 1)
            return value;
          else if (exponent == 2)
            return value * value;
          else if (exponent == 3)
            return value * value * value;
          else if (exponent == 4)
            return value * value * value * value;
          else
            {
              if (exponent % 2 == 0)
                {
                  const int                bisected_exponent = exponent / 2;
                  const LinearOperatorType half_value =
                    LOFixedPower::value(value, bisected_exponent);
                  return half_value * half_value;
                }
              else
                {
                  const int bisected_exponent = (exponent - 1) / 2;
                  const LinearOperatorType nearly_half_value =
                    LOFixedPower::value(value, bisected_exponent);
                  return value * nearly_half_value * nearly_half_value;
                }
            }

          AssertThrow(false, ExcInternalError());
          return null_operator(value);
        }
      };
    } // namespace internal


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::pow(
      const SE::Number &other) const
    {
      if (SE::is_a<SE::Integer>(other))
        {
          const SE::Integer &other_nmbr =
            SE::down_cast<const SE::Integer &>(other);
          const int exponent = internal::evaluate_symengine_number(other_nmbr);
          return SE::make_rcp<
            const LinearOperatorWrapper<Range, Domain, Payload>>(
            internal::LOFixedPower::value(this->value, exponent),
            this->name + "^" + dealii::Utilities::to_string(exponent));
        }

      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return other.rpow(*this);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::rpow(
      const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::exp(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::log(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sin(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::cos(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }

    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::tan(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::cot(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::csc(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sec(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asin(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acos(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::atan(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acsc(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asec(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acot(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sinh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::cosh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::tanh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::csch(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sech(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::coth(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asinh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acosh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::atanh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acsch(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asech(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acoth(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::abs(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::floor(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::ceiling(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::erf(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::erfc(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::gamma(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE


#endif
