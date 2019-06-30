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

#  include <deal.II/differentiation/sd/symengine_la_vector_types.h>
#  include <deal.II/differentiation/sd/symengine_number_utilities.h>

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


    //    namespace internal
    //    {
    //      /**
    //       * A base class for vistors that evaluate symbolic expressions
    //       constructed
    //       * from auto-differentiable numbers.
    //       *
    //       * We don't normally use this class, but in some special instances
    //       it must
    //       * be called to evaluate a symbolic expression that, for whatever
    //       reason,
    //       * has not been fully resolved.
    //       *
    //       * @author Jean-Paul Pelteret, Isuru Fernando, 2017
    //       */
    //      template <typename ADNumberType>
    //      class ADNumberWrapperVisitor_Base;
    //
    //
    //      /**
    //       * A vistor class that evaluates symbolic expressions constructed
    //       from
    //       * auto-differentiable numbers.
    //       *
    //       * We don't normally use this class, but in some special instances
    //       it must
    //       * be called to evaluate a symbolic expression that, for whatever
    //       reason,
    //       * has not been fully resolved.
    //       *
    //       * @author Jean-Paul Pelteret, Isuru Fernando, 2017
    //       */
    //      template<typename ADNumberType, typename = void>
    //      class ADNumberWrapperVisitor;
    //
    //
    //      /**
    //       * A class that implements common subexpression elimination
    //       * for lambda function visitor classes. Since SymEngine has
    //       * their own implementation for real number types, we'll
    //       * restrict this classes use to auto-differentiable number
    //       * types.
    //       *
    //       * @author Jean-Paul Pelteret, Isuru Fernando, 2017
    //       */
    //      template<typename ADNumberType, typename FunctionType, typename =
    //      void> class ADNumberWrapperCSELambdaVisitor;
    //
    //
    //      /**
    //       * A class that converts dictionary-based operations performed
    //       * on ADNumberWrappers to the equivalent actions as performed
    //       * using lambda functions. This should thus be quicker in terms
    //       * of evaluation time than ordinary dictionary substitution.
    //       *
    //       * @author Jean-Paul Pelteret, Isuru Fernando, 2017
    //       */
    //      template<typename ADNumberType, typename = void>
    //      class ADNumberWrapperLambdaVisitor_Base;
    //
    //
    //      /**
    //       * A class derived from ADNumberWrapperLambdaVisitor_Base that
    //       * separately adds extensions valid only for real or complex
    //       * numbers.
    //       *
    //       * @author Jean-Paul Pelteret, Isuru Fernando, 2017
    //       */
    //      template<typename ADNumberType, typename = void>
    //      class ADNumberWrapperLambdaVisitor;
    //
    //    }


    /* --------------------------- inline and template functions
     * ------------------------- */


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


      //      /**
      //       * Evaluate a SymEngine number to a real auto-differentiable
      //       number
      //       */
      //      template<typename ReturnType>
      //      struct Evaluator<ReturnType,
      //        typename std::enable_if<
      //        AD::is_real_valued_ad_number<ReturnType>::value == true
      //        >::type>
      //      {
      //
      //        /**
      //         * This function is OK for ADOL-C number types
      //         * but not Sacado::Fad::DFad.
      //         * With Trilinos's copy constructor, returning by copy results
      //         in a
      //         * Trilinos Sacado::Fad::DFad<double> with the data
      //         * "10 [ 2 2 ]" gets copied to "10 [ 2 2 ] [ ]"
      //         * This means that its transformed from a number with
      //         * value "10" (a double) and sensitivities "2 ; 2" (both
      //         doubles) to one with
      //         * "value" "10 [2;2]" and no sensitivities. Thats not good :-/
      //         * See test/symengine/sd_ad-function_verification_01b.cc
      //         * and test/symengine/sd_ad-function_verification_02b.cc
      //         */
      //        static ReturnType
      //        evaluate (const SE::RCP<const SE::Basic> &value)
      //        {
      //          if (internal::is_an_ADNumberWrapper<ReturnType>(value))
      //            {
      //// AssertThrow(Differentiation::AD::is_sacado_number<ReturnType>::value
      ///== false, /                          ExcMessage("This evaluate()
      /// function is not safe for returning Trilinos types. " / "Please use the
      /// alternative version that takes in the return value " / "as a function
      /// argument."));
      ////
      ////              return SE::down_cast<const ADNumberWrapper<ReturnType>
      ///&>(*value).get_value();
      //
      //              if
      //              (Differentiation::AD::is_sacado_number<ReturnType>::value
      //              == false)
      //                return SE::down_cast<const ADNumberWrapper<ReturnType>
      //                &>(*value).get_value();
      //              else
      //                {
      //                  // Sacado doesn't know how to initialise its own
      //                  (nested)
      //                  // data types in a safe way via copy constructors.
      //                  However,
      //                  // all is good if we first perform a deep copy and
      //                  then
      //                  // return the result. Hooray for efficiency :-/
      //                  ReturnType out (0.0);
      //                  evaluate(out,value);
      //                  return out;
      //                }
      //            }
      //          else if (SE::is_a_Number(*value))
      //            {
      //              // Its possible that we've evaluated a term down to a
      //              constant return ReturnType(Evaluator<typename
      //              AD::ADNumberTraits<ReturnType>::scalar_type>::value(value));
      //            }
      //          else
      //            {
      //              // Deal with some special cases wherein there exist
      //              remaining
      //              // operations that need to be processed
      //              try
      //                {
      //                  ADNumberWrapperVisitor<ReturnType> visitor;
      //                  return visitor.apply(*value);
      //                }
      //              catch (...)
      //                {
      //                  std::ostringstream sstr;
      //                  sstr << *value;
      //                  AssertThrow(false,
      //                              ExcMessage("The underlying SymEngine
      //                              number type is either not an
      //                              ADNumberWrapper "
      //                                         "or it is a special case for
      //                                         ADNumberWrapper could not be
      //                                         evaluated. \n" "Value: " +
      //                                         sstr.str()));
      //                  return ReturnType(0.0);
      //                }
      //            }
      //        }
      //
      //
      //        /**
      //         * This function is OK for ADOL-C number types
      //         * as well as Sacado::Fad::DFad. This is the
      //         * solution to the problems seen in
      //         * See test/symengine/sd_ad-function_verification_01b.cc
      //         * and test/symengine/sd_ad-function_verification_02b.cc
      //         */
      //        static void
      //        evaluate (ReturnType &out,
      //                  const SE::RCP<const SE::Basic> &value)
      //        {
      //          if (internal::is_an_ADNumberWrapper<ReturnType>(value))
      //            {
      //              out = SE::down_cast<const ADNumberWrapper<ReturnType>
      //              &>(*value).get_value();
      //            }
      //          else if (SE::is_a_Number(*value))
      //            {
      //              // Its possible that we've evaluated a term down to a
      //              constant out = ReturnType(Evaluator<typename
      //              AD::ADNumberTraits<ReturnType>::scalar_type>::value(value));
      //            }
      //          else
      //            {
      //              // Deal with some special cases wherein there exist
      //              remaining
      //              // operations that need to be processed
      //              try
      //                {
      //                  ADNumberWrapperVisitor<ReturnType> visitor;
      //                  out = visitor.apply(*value);
      //                }
      //              catch (...)
      //                {
      //                  std::ostringstream sstr;
      //                  sstr << *value;
      //                  AssertThrow(false,
      //                              ExcMessage("The underlying SymEngine
      //                              number type is either not an
      //                              ADNumberWrapper "
      //                                         "or it is a special case for
      //                                         ADNumberWrapper could not be
      //                                         evaluated. \n" "Value: " +
      //                                         sstr.str()));
      //                  out =
      //                  dealii::internal::NumberType<ReturnType>::value(0.0);
      //                }
      //            }
      //        }
      //      };
      //
      //
      //      /**
      //       * Evaluate a SymEngine number to a complex auto-differentiable
      //       number
      //       */
      //      template<typename ReturnType>
      //      struct Evaluator<ReturnType,
      //        typename std::enable_if<
      //        AD::is_complex_valued_ad_number<ReturnType>::value == true
      //        >::type>
      //      {
      //        static ReturnType
      //        evaluate (const SE::RCP<const SE::Basic> &value)
      //        {
      //          AssertThrow(false, ExcNotImplemented());
      //
      //          return ReturnType(0.0);
      //        }
      //      };
      //
      //
      //
      //      template<typename ReturnType>
      //      struct ADNumber_to_Evaluator<ReturnType,
      //        typename std::enable_if<
      //        std::is_arithmetic<ReturnType>::value
      //        >::type>
      //      {
      //        static ReturnType
      //        evaluate (const SE::RCP<const SE::Basic> &value)
      //        {
      //          // Note: Extra caution should be observed when expecting this
      //          part of the code
      //          //       to be evaluated. This function can flatten an
      //          ADNumberType to,
      //          //       for example, a double that can then be recast to the
      //          ADNumberType without
      //          //       any notice. This may cause you lots of headaches if
      //          its used incorrectly
      //          //       during taping, so be careful!
      //          // Extra note: We cannot use
      //          SE::is_a<ADNumberWrapper<ADNumberType> here to differentiate
      //          //             between the different AD types. In fact, we
      //          cannot use any
      //          //             SE::is_a<NumberWrapper derivative> because this
      //          simply checks type codes
      //          //             (all of which are the same for NumberWrapper
      //          derivatives), not whether
      //          //             they are indeed the correct type. So we'll have
      //          to do this the "slow way"
      //          //             instead of
      //          //             return
      //          AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////#ifdef DEAL_II_WITH_ADOLC
      ////          try
      ////            {
      ////              typedef
      /// AD::NumberTraits<double,Differentiation::AD::NumberTypes::adolc_taped>
      /// ADNumberTraits; /              typedef typename
      /// ADNumberTraits::ad_type ad_type; /              const SE::RCP<const
      /// ADNumberWrapper<ad_type> > ptr_ad_value = SE::rcp_dynamic_cast<const
      /// ADNumberWrapper<ad_type>
      ///>(value); /              return
      /// AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////            }
      ////          catch (const std::runtime_error &error)
      ////            {
      ////              // Nope, that wasn't a wrapped taped ADOL-C number...
      ////              (void)error;
      ////            }
      ////
      ////          try
      ////            {
      ////              typedef
      /// AD::NumberTraits<double,Differentiation::AD::NumberTypes::adolc_tapeless>
      /// ADNumberTraits; /              typedef typename
      /// ADNumberTraits::ad_type ad_type; /              const SE::RCP<const
      /// ADNumberWrapper<ad_type> > ptr_ad_value = SE::rcp_dynamic_cast<const
      /// ADNumberWrapper<ad_type>
      ///>(value); /              return
      /// AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////            }
      ////          catch (const std::runtime_error &error)
      ////            {
      ////              // Nope, that wasn't a wrapped tapeless ADOL-C number...
      ////              (void)error;
      ////            }
      ////#endif
      ////
      ////#ifdef DEAL_II_WITH_TRILINOS
      ////          try
      ////            {
      ////              typedef
      /// AD::NumberTraits<double,Differentiation::AD::NumberTypes::sacado_dfad>
      /// ADNumberTraits; /              typedef typename
      /// ADNumberTraits::ad_type ad_type; /              const SE::RCP<const
      /// ADNumberWrapper<ad_type> > ptr_ad_value = SE::rcp_dynamic_cast<const
      /// ADNumberWrapper<ad_type>
      ///>(value); /              return
      /// AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////            }
      ////          catch (const std::runtime_error &error)
      ////            {
      ////              // Nope, that wasn't a wrapped Sacado::Fad::DFad
      /// number... /              (void)error; /            }
      ////
      ////          try
      ////            {
      ////              typedef
      /// AD::NumberTraits<double,Differentiation::AD::NumberTypes::sacado_dfad_dfad>
      /// ADNumberTraits; /              typedef typename
      /// ADNumberTraits::ad_type ad_type; /              const SE::RCP<const
      /// ADNumberWrapper<ad_type> > ptr_ad_value = SE::rcp_dynamic_cast<const
      /// ADNumberWrapper<ad_type>
      ///>(value); /              return
      /// AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////            }
      ////          catch (const std::runtime_error &error)
      ////            {
      ////              // Nope, that wasn't a wrapped nested Sacado::Fad::DFad
      /// number... /              (void)error; /            }
      ////
      ////          try
      ////            {
      ////              typedef
      /// AD::NumberTraits<double,Differentiation::AD::NumberTypes::sacado_rad>
      /// ADNumberTraits; /              typedef typename
      /// ADNumberTraits::ad_type ad_type; /              const SE::RCP<const
      /// ADNumberWrapper<ad_type> > ptr_ad_value = SE::rcp_dynamic_cast<const
      /// ADNumberWrapper<ad_type>
      ///>(value); /              return
      /// AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////            }
      ////          catch (const std::runtime_error &error)
      ////            {
      ////              // Nope, that wasn't a wrapped Sacado::Rad number...
      ////              (void)error;
      ////            }
      ////
      ////          try
      ////            {
      ////              typedef
      /// AD::NumberTraits<double,Differentiation::AD::NumberTypes::sacado_rad_dfad>
      /// ADNumberTraits; /              typedef typename
      /// ADNumberTraits::ad_type ad_type; /              const SE::RCP<const
      /// ADNumberWrapper<ad_type> > ptr_ad_value = SE::rcp_dynamic_cast<const
      /// ADNumberWrapper<ad_type>
      ///>(value); /              return
      /// AD::internal::NumberType<ReturnType>::value(Evaluator<ad_type>::value(value));
      ////            }
      ////          catch (const std::runtime_error &error)
      ////            {
      ////              // Nope, that wasn't a wrapped nested
      /// Sacado::Rad<Sacado::Fad::DFad> number... /              (void)error; /
      ///}
      ////#endif
      //
      //#define AD_Evaluation_Branches(ADNumberTypeCode) \
//  { \
//    if (scalar_type_code == ScalarNumberTypes::is_float) \
//      { \
//        typedef float ScalarType; \
//        return
      //        AD::internal::NumberType<ReturnType>::value(evaluate<ScalarType,ADNumberTypeCode>(value));
      //        \
//      } \
//    else if (scalar_type_code == ScalarNumberTypes::is_double) \
//      { \
//        typedef double ScalarType; \
//        return
      //        AD::internal::NumberType<ReturnType>::value(evaluate<ScalarType,ADNumberTypeCode>(value));
      //        \
//      } \
//    else \
//      { \
//        AssertThrow(false, \
//                    ExcMessage("Evaluator not implemented for this
      //                    scalar type and"                       \
//                               "AD number type code "
      //                               #ADNumberTypeCode)); \
//      } \
//  }
      //
      //          if (const ADNumberWrapperBase *const p_base =
      //          dynamic_cast<const ADNumberWrapperBase *const>(&*value))
      //            {
      //              const enum ScalarNumberTypes &scalar_type_code =
      //              p_base->scalar_type_code; const enum AD::NumberTypes
      //              &ad_type_code = p_base->ad_type_code;
      //
      //#ifdef DEAL_II_WITH_ADOLC
      //              if (ad_type_code == AD::NumberTypes::adolc_taped)
      //                {
      //                  AD_Evaluation_Branches(AD::NumberTypes::adolc_taped);
      //                }
      //              else if (ad_type_code == AD::NumberTypes::adolc_tapeless)
      //                {
      //                  AD_Evaluation_Branches(AD::NumberTypes::adolc_tapeless);
      //                }
      //#endif
      //
      //#ifdef DEAL_II_WITH_TRILINOS
      //              if (ad_type_code == AD::NumberTypes::sacado_dfad)
      //                {
      //                  AD_Evaluation_Branches(AD::NumberTypes::sacado_dfad);
      //                }
      //              else if (ad_type_code ==
      //              AD::NumberTypes::sacado_dfad_dfad)
      //                {
      //                  AD_Evaluation_Branches(AD::NumberTypes::sacado_dfad_dfad);
      //                }
      //              else if (ad_type_code == AD::NumberTypes::sacado_rad)
      //                {
      //                  AD_Evaluation_Branches(AD::NumberTypes::sacado_rad);
      //                }
      //              else if (ad_type_code == AD::NumberTypes::sacado_rad_dfad)
      //                {
      //                  AD_Evaluation_Branches(AD::NumberTypes::sacado_rad_dfad);
      //                }
      //#endif
      //
      //              AssertThrow(false, ExcMessage("Evaluator not implemented
      //              for this AD number type"));
      //            }
      //
      //          AssertThrow(false, ExcMessage("Unable to convert this wrapped
      //          AD number to a float")); return ReturnType(0.0);
      //        }
      //
      //#undef AD_Evaluation_Branches
      //
      //      private:
      //
      //        template<typename ScalarType, enum AD::NumberTypes
      //        ADNumberTypeCode> static typename
      //        AD::NumberTraits<ScalarType,ADNumberTypeCode>::ad_type evaluate
      //        (const SE::RCP<const SE::Basic> &value)
      //        {
      //          typedef AD::NumberTraits<ScalarType,ADNumberTypeCode>
      //          ADNumberTraits; typedef typename ADNumberTraits::ad_type
      //          ad_type; Assert(is_an_ADNumberWrapper<ad_type>(value),
      //          ExcInternalError());
      //
      //          try
      //            {
      //              const SE::RCP<const ADNumberWrapper<ad_type> >
      //              ptr_ad_value = SE::rcp_dynamic_cast<const
      //              ADNumberWrapper<ad_type> >(value); return
      //              Evaluator<ad_type>::value(value);
      //            }
      //          catch (const std::runtime_error &error)
      //            {
      //              AssertThrow(false, ExcMessage("Something went wrong during
      //              AD number evaluation")); return ad_type(0.0);
      //            }
      //        }
      //      };


      //      /* --------------------------- ADNumberWrapperVisitor_Base
      //      ------------------------- */
      //
      //
      //      template <typename ADNumberType>
      //      class ADNumberWrapperVisitor_Base : public
      //      SE::BaseVisitor<ADNumberWrapperVisitor_Base<ADNumberType> >
      //      {
      //      protected:
      //
      //        // The 'result' variable is assigned into at the very end of
      //        each visit()
      //        // methods below. The only place where these methods are called
      //        from is the
      //        // line 'b.accept(*this)' in apply() and the 'result' is
      //        immediately
      //        // returned. Thus no corruption can happen and apply() can be
      //        safely called
      //        // recursively.
      //        ADNumberType result;
      //
      //      public:
      //
      //        ADNumberType
      //        apply(const SE::Basic &b);
      //
      //        void
      //        bvisit (const SE::Constant &x);
      //
      //        void
      //        bvisit(const SE::Integer &x);
      //
      //        void
      //        bvisit(const SE::Rational &x);
      //
      //        void
      //        bvisit(const SE::RealDouble &x);
      //
      //        void
      //        bvisit(const SE::NumberWrapper &x);
      //
      //        void
      //        bvisit (const ADNumberWrapper<ADNumberType> &x);
      //
      //        void
      //        bvisit(const SE::Symbol &);
      //
      //        void
      //        bvisit(const SE::Basic &);
      //
      //        void
      //        bvisit(const SE::Add &op);
      //
      //        void
      //        bvisit(const SE::Mul &op);
      //
      //        void
      //        bvisit(const SE::Pow &op);
      //
      //        void
      //        bvisit(const SE::Log &op);
      //
      //        void
      //        bvisit(const SE::Sin &op);
      //
      //        void
      //        bvisit(const SE::Cos &op);
      //
      //        void
      //        bvisit(const SE::Tan &op);
      //
      //        void
      //        bvisit(const SE::Csc &op);
      //
      //        void
      //        bvisit(const SE::Sec &op);
      //
      //        void
      //        bvisit(const SE::Cot &op);
      //
      //        void
      //        bvisit(const SE::ASin &op);
      //
      //        void
      //        bvisit(const SE::ACos &op);
      //
      //        void
      //        bvisit(const SE::ATan &op);
      //
      //        void
      //        bvisit(const SE::ASec &op);
      //
      //        void
      //        bvisit(const SE::ACsc &op);
      //
      //        void
      //        bvisit(const SE::ACot &op);
      //
      //        void
      //        bvisit(const SE::Sinh &op);
      //
      //        void
      //        bvisit(const SE::Cosh &op);
      //
      //        void
      //        bvisit(const SE::Tanh &op);
      //
      //        void
      //        bvisit(const SE::Csch &op);
      //
      //        void
      //        bvisit(const SE::Sech &op);
      //
      //        void
      //        bvisit(const SE::Coth &op);
      //
      //        void
      //        bvisit(const SE::ASinh &op);
      //
      //        void
      //        bvisit(const SE::ACosh &op);
      //
      //        void
      //        bvisit(const SE::ATanh &op);
      //
      //        void
      //        bvisit(const SE::ACsch &op);
      //
      //        void
      //        bvisit(const SE::ASech &op);
      //
      //        void
      //        bvisit(const SE::ACoth &op);
      //
      //        void
      //        bvisit(const SE::Abs &op);
      //
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperVisitor
      //      ------------------------- */
      //
      //
      //      template<typename ADNumberType>
      //      class ADNumberWrapperVisitor<ADNumberType, typename
      //      std::enable_if<
      //        AD::is_real_valued_ad_number<ADNumberType>::value>::type>
      //        : public SE::BaseVisitor<
      //        ADNumberWrapperVisitor<ADNumberType, typename std::enable_if<
      //        AD::is_real_valued_ad_number<ADNumberType>::value>::type>,
      //        ADNumberWrapperVisitor_Base<ADNumberType> >
      //      {
      //        typedef typename AD::ADNumberTraits<ADNumberType>::real_type
      //        ADRealType; typedef typename
      //        AD::ADNumberTraits<ADNumberType>::complex_type ADComplexType;
      //
      //        static_assert(AD::is_ad_number<ADNumberType>::value == true,
      //                      "Only auto-differentiable numbers can be used with
      //                      this class."
      //                     );
      //        static_assert(boost::is_complex<ADComplexType>::value == true,
      //                      "Expected a complex number type from
      //                      AD::NumberTraits."
      //                     );
      //
      //      public:
      //
      //        // Classes not implemented are
      //        // Subs, Gamma, LogGamma, UpperGamma, LowerGamma, Dirichlet_eta,
      //        // Zeta, LeviCivita, KroneckerDelta, LambertW
      //        // Derivative, Complex, ComplexDouble, ComplexMPC
      //
      //        using ADNumberWrapperVisitor_Base<ADNumberType>::bvisit;
      //        using ADNumberWrapperVisitor_Base<ADNumberType>::apply;
      //        using ADNumberWrapperVisitor_Base<ADNumberType>::result;
      //
      //        void
      //        bvisit (const SE::ATan2 &op);
      //
      //        void
      //        bvisit (const SE::Erf &op);
      //
      //        void
      //        bvisit (const SE::Erfc &op);
      //
      //        void
      //        bvisit (const SE::Max &op);
      //
      //        void
      //        bvisit (const SE::Min &op);
      //      };
      //
      //      template<typename ADNumberType>
      //      class ADNumberWrapperVisitor<ADNumberType, typename
      //      std::enable_if<
      //        AD::is_complex_valued_ad_number<ADNumberType>::value>::type>
      //        : public SE::BaseVisitor<
      //        ADNumberWrapperVisitor<ADNumberType, typename std::enable_if<
      //        AD::is_complex_valued_ad_number<ADNumberType>::value>::type>,
      //        ADNumberWrapperVisitor_Base<ADNumberType> >
      //      {
      //        typedef typename AD::ADNumberTraits<ADNumberType>::real_type
      //        ADRealType; typedef typename
      //        AD::ADNumberTraits<ADNumberType>::complex_type ADComplexType;
      //
      //        static_assert(AD::is_ad_number<ADNumberType>::value == true,
      //                      "Only auto-differentiable numbers can be used with
      //                      this class."
      //                     );
      //        static_assert(boost::is_complex<ADComplexType>::value == true,
      //                      "Expected a complex number type from
      //                      AD::NumberTraits."
      //                     );
      //
      //      public:
      //        // Classes not implemented are
      //        // ATan2, Erf, Erfc, Max, Min,
      //        // Subs, UpperGamma, LowerGamma, Dirichlet_eta, Zeta
      //        // LeviCivita, KroneckerDelta, FunctionSymbol, LambertW
      //        // Derivative, ComplexMPC
      //
      //        using ADNumberWrapperVisitor_Base<ADNumberType>::bvisit;
      //        using ADNumberWrapperVisitor_Base<ADNumberType>::apply;
      //        using ADNumberWrapperVisitor_Base<ADNumberType>::result;
      //
      //        void
      //        bvisit(const SE::Complex &x);
      //
      //        void
      //        bvisit(const SE::ComplexDouble &x);
      //      };
      //
      //
      //      /* --------------------------- UseStackedLambdaVisitor
      //      ------------------------- */
      //
      //
      //      /**
      //       * Define a helper class that decides whether we'll use a
      //       * memory pool to record/store the
      //       */
      //      template<typename ADNumberType, typename = void>
      //      struct UseStackedLambdaVisitor : std::false_type
      //      {
      //        static_assert(AD::is_tapeless_ad_number<ADNumberType>::value,
      //                      "Should not be using stacked Lambda visitor!");
      //        static_assert(!std::is_same<ADNumberType,adouble>::value,
      //                      "Should be using stacked Lambda visitor!");
      //      };
      //
      //
      //      /**
      //       * Taped AD numbers are the primary motivation for storing
      //       * data in a shared memory pool
      //       */
      //      template<typename ADNumberType>
      //      struct UseStackedLambdaVisitor<ADNumberType, typename
      //      std::enable_if<
      //        AD::is_taped_ad_number<ADNumberType>::value
      //        >::type> : std::true_type
      //      {};
      //
      //
      //      /**
      //       * Trilinos numbers appear to have some difficulty in
      //       * using the alternate implementation of the
      //       * LambdaVisitor or CSELambdaVisitor.
      //       * If we don't make this specialisation then the tests
      //       * symengine/sd_ad-cse_01_[6,8,10,12]
      //       * symengine/step-44-sd_ad-helper_res_lin_0[3-6]a_1
      //       * symengine/step-44-sd_ad-helper_res_lin_0[3-6]a_2
      //       * will fail. Presumably this has to do with the degedation
      //       * of data integrity when using copy constructors
      //       * (the alternate implementation uses deep copying to
      //       * transfer data between data entries).
      //       */
      ////      template<typename ADNumberType>
      ////      struct UseStackedLambdaVisitor<ADNumberType, typename
      /// std::enable_if< /        AD::is_sacado_number<ADNumberType>::value /
      ///>::type> : std::true_type /      {};
      //
      //
      //      /* ---------------------------
      //      ADNumberWrapperCSELambdaVisitor_Base ------------------------- */
      //
      //
      //      /* --------------------------- ADNumberWrapperCSELambdaVisitor
      //      (stacked) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for taped auto-differentiable types
      //       */
      //      template<typename ADNumberType, typename FunctionType>
      //      class ADNumberWrapperCSELambdaVisitor<ADNumberType, FunctionType,
      //      typename std::enable_if<
      //        internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //      {
      //        /**
      //         * Intermediate functions
      //         */
      //        std::vector<FunctionType>  intermediate_functions;
      //
      //        /**
      //         * Values  reduced expressions
      //         */
      //        std::vector<ADNumberType>  intermediate_results;
      //
      //        /**
      //         * A lookup table that provides the index in the
      //         * @p intermediate_results corresponding to a
      //         * tested symbol.
      //         */
      //        std::map<SE::RCP<const SE::Basic>, unsigned int,
      //        SE::RCPBasicKeyLess> intermediate_functions_map;
      //
      //      public:
      //
      //        /**
      //         * Perform common subexpression elimination
      //         *
      //         * Here we build the reduced expressions for the @p dependent_functions
      //         * as well as a list of intermediate, repeating symbolic
      //         expressions that
      //         * are extracted @p dependent_functions. This operation leads to the
      //         * elimination of repeated expressions, so they only have to be
      //         evaluated
      //         * once.
      //         */
      //        void
      //        init(const SE::vec_basic &dependent_functions,
      //             std::vector<FunctionType>                      &results,
      //             std::function<FunctionType(const SE::Basic &)> &apply)
      //        {
      //          intermediate_functions.clear();
      //          intermediate_results.clear();
      //          intermediate_functions_map.clear();
      //
      //          SE::vec_pair  intermediate_symbols_exprs;
      //          SE::vec_basic reduced_exprs;
      //
      //          // After the next call, the data stored in replacements is
      //          structured
      //          // as follows:
      //          //
      //          // replacements[i] := [f, f(x)]
      //          // replacements[i].first  = intermediate function label "f"
      //          // replacements[i].second = intermediate function definition
      //          "f(x)"
      //          //
      //          // It is to be evaluated top down (i.e. index 0 to
      //          replacements.size()),
      //          // with the results going back into the substitution map for
      //          // the next levels. So for each "i", "x" are the superset of
      //          the input
      //          // values and the previously evaluated [f_0(x), f_1(x), ...,
      //          f_{i-1}(x)].
      //          //
      //          // The final result is a set of reduced expressions
      //          // that must be computed after the replacement
      //          // values have been computed.
      //          SE::cse(intermediate_symbols_exprs, reduced_exprs,
      //          dependent_functions);
      //
      //          for (auto &replacement : intermediate_symbols_exprs)
      //            {
      //              const auto intermediate_result =
      //              apply(*(replacement.second));
      //              // Store the replacement symbol values in a dictionary for
      //              // faster lookup for initialization
      //              intermediate_functions_map[replacement.first]
      //                = intermediate_functions.size();
      //              // Store it in a vector for faster use in call
      //              intermediate_functions.push_back(intermediate_result);
      //            }
      //
      //          intermediate_results.resize(intermediate_functions.size());
      //
      //          // Generate functions for all the reduced expressions
      //          // owned by the calling class / function. This call to apply()
      //          // may in turn indirectly call our specialisation for bvisit,
      //          // which expects a non-empty intermediate function map.
      //          AssertThrow(reduced_exprs.size() ==
      //          dependent_functions.size(),
      //                      ExcDimensionMismatch(reduced_exprs.size(),
      //                      dependent_functions.size()));
      //          for (unsigned i = 0; i < dependent_functions.size(); i++)
      //            results.push_back(apply(*reduced_exprs[i]));
      //
      //          // We don't need the intermediate function map anymore,
      //          // so it can be safely cleared
      //          intermediate_functions_map.clear();
      //        }
      //
      //        template<typename ReturnValueType>
      //        void
      //        call(const ADNumberType *substitution_values)
      //        {
      //          Assert(n_intermediate_expressions() > 0, ExcInternalError());
      //          Assert(intermediate_results.size() ==
      //          intermediate_functions.size(),
      //                 ExcDimensionMismatch(intermediate_results.size(),intermediate_functions.size()));
      //
      //          ReturnValueType out (0.0, false);
      //          for (unsigned i = 0; i < intermediate_functions.size(); ++i)
      //            {
      //              intermediate_functions[i](out, substitution_values);
      //              intermediate_results[i] = out.value;
      //              out.set(0.0,false);
      //            }
      //        }
      //
      //
      //        // Stacked type
      //        template<typename ReturnValueType>
      //        bool
      //        bvisit(const SE::Symbol &symb,
      //               FunctionType     &func)
      //        {
      //          Assert(intermediate_functions_map.size() > 0,
      //          ExcNotInitialized());
      //
      //          const auto it =
      //          intermediate_functions_map.find(symb.rcp_from_this()); if (it
      //          != intermediate_functions_map.end())
      //            {
      //              const unsigned index = it->second;
      //              func = [=](ReturnValueType &out, const ADNumberType *)
      //              {
      //                out.set(intermediate_results[index],true);
      //              };
      //
      //              return true; // Found
      //            }
      //
      //          return false; // Not found
      //        }
      //
      //        template <typename Stream>
      //        void
      //        print (Stream &stream) const
      //        {
      //          stream
      //              << "Common subexpression elimination: \n"
      //              << "  Number of intermediate reduced expressions: "
      //              << n_intermediate_expressions()
      //              << std::endl;
      //        }
      //
      //        /**
      //         * Return a flag stating whether we've performed CSE or not.
      //         */
      //        bool
      //        executed() const
      //        {
      //          return n_intermediate_expressions() > 0;
      //        }
      //
      //        /**
      //         * The number of intermediate expressions that must
      //         * be evaluated as part of the collection of common
      //         * subexpressions.
      //         */
      //        unsigned int
      //        n_intermediate_expressions() const
      //        {
      //          return intermediate_functions.size();
      //        }
      //      };
      //
      //
      //
      //      /* --------------------------- ADNumberWrapperCSELambdaVisitor
      //      (non-stacked) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for tapeless auto-differentiable types
      //       */
      //      template<typename ADNumberType, typename FunctionType>
      //      class ADNumberWrapperCSELambdaVisitor<ADNumberType, FunctionType,
      //      typename std::enable_if<
      //        !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //      {
      //        /**
      //         * Intermediate functions
      //         */
      //        std::vector<FunctionType>  intermediate_functions;
      //
      //        /**
      //         * Values  reduced expressions
      //         */
      //        std::vector<ADNumberType>  intermediate_results;
      //
      //        /**
      //         * A lookup table that provides the index in the
      //         * @p intermediate_results corresponding to a
      //         * tested symbol.
      //         */
      //        std::map<SE::RCP<const SE::Basic>, unsigned int,
      //        SE::RCPBasicKeyLess> intermediate_functions_map;
      //
      //      public:
      //
      //        /**
      //         * Perform common subexpression elimination
      //         *
      //         * Here we build the reduced expressions for the @p dependent_functions
      //         * as well as a list of intermediate, repeating symbolic
      //         expressions that
      //         * are extracted @p dependent_functions. This operation leads to the
      //         * elimination of repeated expressions, so they only have to be
      //         evaluated
      //         * once.
      //         */
      //        void
      //        init(const SE::vec_basic &dependent_functions,
      //             std::vector<FunctionType>                      &results,
      //             std::function<FunctionType(const SE::Basic &)> &apply)
      //        {
      //          intermediate_functions.clear();
      //          intermediate_results.clear();
      //          intermediate_functions_map.clear();
      //
      //          SE::vec_pair  intermediate_symbols_exprs;
      //          SE::vec_basic reduced_exprs;
      //
      //          // After the next call, the data stored in replacements is
      //          structured
      //          // as follows:
      //          //
      //          // replacements[i] := [f, f(x)]
      //          // replacements[i].first  = intermediate function label "f"
      //          // replacements[i].second = intermediate function definition
      //          "f(x)"
      //          //
      //          // It is to be evaluated top down (i.e. index 0 to
      //          replacements.size()),
      //          // with the results going back into the substitution map for
      //          // the next levels. So for each "i", "x" are the superset of
      //          the input
      //          // values and the previously evaluated [f_0(x), f_1(x), ...,
      //          f_{i-1}(x)].
      //          //
      //          // The final result is a set of reduced expressions
      //          // that must be computed after the replacement
      //          // values have been computed.
      //          SE::cse(intermediate_symbols_exprs, reduced_exprs,
      //          dependent_functions);
      //
      //          for (auto &replacement : intermediate_symbols_exprs)
      //            {
      //              const auto intermediate_result =
      //              apply(*(replacement.second));
      //              // Store the replacement symbol values in a dictionary for
      //              // faster lookup for initialization
      //              intermediate_functions_map[replacement.first]
      //                = intermediate_functions.size();
      //              // Store it in a vector for faster use in call
      //              intermediate_functions.push_back(intermediate_result);
      //            }
      //
      //          intermediate_results.resize(intermediate_functions.size());
      //
      //          // Generate functions for all the reduced expressions
      //          // owned by the calling class / function. This call to apply()
      //          // mayin turn indirectly call our specialisation for bvisit,
      //          // which expects a non-empty intermediate function map.
      //          AssertThrow(reduced_exprs.size() ==
      //          dependent_functions.size(),
      //                      ExcDimensionMismatch(reduced_exprs.size(),
      //                      dependent_functions.size()));
      //          for (unsigned i = 0; i < dependent_functions.size(); i++)
      //            results.push_back(apply(*reduced_exprs[i]));
      //
      //          // We don't need the intermediate function map anymore,
      //          // so it can be safely cleared
      //          intermediate_functions_map.clear();
      //        }
      //
      //
      //        void
      //        call(const ADNumberType *substitution_values)
      //        {
      //          Assert(n_intermediate_expressions() > 0, ExcInternalError());
      //          Assert(intermediate_results.size() ==
      //          intermediate_functions.size(),
      //                 ExcDimensionMismatch(intermediate_results.size(),intermediate_functions.size()));
      //
      //          for (unsigned i = 0; i < intermediate_functions.size(); ++i)
      //            intermediate_results[i] =
      //            intermediate_functions[i](substitution_values);
      //        }
      //
      //
      //        bool
      //        bvisit(const SE::Symbol &symb,
      //               FunctionType     &func)
      //        {
      //          Assert(intermediate_functions_map.size() > 0,
      //          ExcNotInitialized());
      //
      //          const auto it =
      //          intermediate_functions_map.find(symb.rcp_from_this()); if (it
      //          != intermediate_functions_map.end())
      //            {
      //              const unsigned index = it->second;
      //              func = [=](const ADNumberType *)
      //              {
      //                return intermediate_results[index];
      //              };
      //
      //              return true; // Found
      //            }
      //
      //          return false; // Not found
      //        }
      //
      //        template <typename Stream>
      //        void
      //        print (Stream &stream) const
      //        {
      //          stream
      //              << "Common subexpression elimination: \n"
      //              << "  Number of intermediate reduced expressions: "
      //              << n_intermediate_expressions()
      //              << std::endl;
      //        }
      //
      //        /**
      //         * Return a flag stating whether we've performed CSE or not.
      //         */
      //        bool
      //        executed() const
      //        {
      //          return n_intermediate_expressions() > 0;
      //        }
      //
      //        /**
      //         * The number of intermediate expressions that must
      //         * be evaluated as part of the collection of common
      //         * subexpressions.
      //         */
      //        unsigned int
      //        n_intermediate_expressions() const
      //        {
      //          return intermediate_functions.size();
      //        }
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperLambdaVisitor_Base
      //      (stacked) ------------------------- */
      //
      //
      //
      //      // This annoying little struct is necessary to circumvent
      //      // some current deficiencies in ADOL-C.
      //      // The power function (and perhaps others) does not admit
      //      // an arbitrary range for the base and exponent if both of
      //      // the arguments are ADOL-C doubles.
      //      // However, our lambda function must also have a preset
      //      // return type, which is always an ADOL-C number.
      //      // At the same time, none of the stored SymEngine numbers
      //      // within the SymEngine ops can tell us if its referencing
      //      // the result of an ADOL-C operation or not.
      //      // To satisfy this conflict, use this struct to keep tabs
      //      // on the state of the result (i.e. whether it really stores
      //      // a sensitive ADOL-C number or not) and modify some
      //      // operations based on its contents.
      //      //
      //      // This deals with the problem that, when using lambda
      //      // functions, the return type is fixed. This means that we
      //      // always cast ints, doubles etc. to ADNumbers. For AD math
      //      // functions, like pow, there are limitations when all arguments
      //      // are ADNumbers. So in this patch I track the source of the
      //      // result of operations and mark whether the value really is
      //      // an ADNumber or just a normal number saved as an ADNumber.
      //      // This allows one to manually switch the math functions that
      //      // are called, and prevent a range restriction for the evaluated
      //      // numbers.
      //      template <typename ADNumberType>
      //      struct ReturnValue
      //      {
      //        ReturnValue (const ADNumberType &value,
      //                     const bool         &has_sensitivities)
      //          : value (value),
      //            has_sensitivities (has_sensitivities)
      //        {}
      //
      //        // By default, we will assume that we're storing
      //        // a value that tracks some independent data
      //        ReturnValue (const ADNumberType &value)
      //          : ReturnValue (value,true)
      //        {}
      //
      //        template<typename T>
      //        void
      //        set (const T    &value_,
      //             const bool &has_sensitivities_)
      //        {
      //          value = value_;
      //          has_sensitivities = has_sensitivities_;
      //        }
      //
      //        void
      //        reset()
      //        {
      //          value =
      //          dealii::internal::NumberType<ADNumberType>::value(0.0);
      //          has_sensitivities = false;
      //        }
      //
      //        /**
      //         * For optimum performance, this number should be generated
      //         * from some shared memory pool.
      //         */
      //        ADNumberType value;
      //        bool         has_sensitivities;
      //      };
      //
      //      template <typename ADNumberType>
      //      bool
      //      has_sensitivities (const ReturnValue<ADNumberType> &lhs,
      //                         const ReturnValue<ADNumberType> &rhs)
      //      {
      //        return lhs.has_sensitivities | rhs.has_sensitivities;
      //      }
      //
      //      template<typename ADNumberType,typename... Args>
      //      bool
      //      has_sensitivities (const ReturnValue<ADNumberType> &arg1,
      //                         const ReturnValue<ADNumberType> &arg2,
      //                         const Args &...                  args)
      //      {
      //        return has_sensitivities(args..., arg1) | arg2.has_sensitivities;
      //      }
      //
      //
      //      /**
      //       * Class definition for taped auto-differentiable types
      //       */
      //      template <typename ADNumberType>
      //      class ADNumberWrapperLambdaVisitor_Base<ADNumberType, typename
      //      std::enable_if<
      //        internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //        : public
      //        SE::BaseVisitor<ADNumberWrapperLambdaVisitor_Base<ADNumberType,
      //        typename std::enable_if<
      //        internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type> >
      //      {
      //      protected:
      //
      //        typedef ReturnValue<ADNumberType> ReturnValueType;
      //
      //        /**
      //         * A typedef for the evaluation function signature
      //         */
      //        typedef std::function<void(ReturnValueType &out, const
      //        ADNumberType *values)> FunctionType;
      //
      //        /**
      //         * The set of symbols that represent the independent variables
      //         */
      //        SE::vec_basic symbols;
      //
      //        /**
      //         * The set of constructed lambda functions that represent
      //         * evaluation of the dependent functions
      //         */
      //        std::vector<FunctionType> functions;
      //
      //        /**
      //         * An intermediate value used to build the vector of dependent
      //         functions.
      //         *
      //         * The @p func variable is assigned into at the very end of each visit()
      //         * methods below. The only place where these methods are called
      //         from is the
      //         * line 'b.accept(*this)' in apply() and the @p func is immediately
      //         * returned. Thus no corruption can happen and apply() can be
      //         safely called
      //         * recursively.
      //         */
      //        FunctionType func;
      //
      //        /**
      //         * A helper class that will be used to perform common
      //         * subexpression elimination
      //         */
      //        ADNumberWrapperCSELambdaVisitor<ADNumberType,FunctionType> cse;
      //
      //        /**
      //         * Type definition for a shared memory pool
      //         */
      //        typedef
      //        dealii::internal::SharedPool<ReturnValueType,dealii::internal::Return_to_Pool<ReturnValueType>
      //        > ValuePoolType;
      //
      //        /**
      //         * A storage stack that will be used to store intermediate
      //         results
      //         */
      //        ValuePoolType value_pool;
      //
      //      public:
      //
      //        ADNumberWrapperLambdaVisitor_Base ();
      //
      //        virtual ~ADNumberWrapperLambdaVisitor_Base ();
      //
      //        void
      //        init (const SE::vec_basic &symbs,
      //              const SE::Basic     &func,
      //              const bool           use_cse = false);
      //
      //        void
      //        init(const SE::vec_basic &symbs,
      //             const SE::vec_basic &funcs,
      //             const bool           use_cse = false);
      //
      //        FunctionType
      //        apply (const SE::Basic &x);
      //
      //        ADNumberType
      //        call (const std::vector<ADNumberType> &values);
      //
      //        void
      //        call(ADNumberType       *outputs,
      //             const ADNumberType *values);
      //
      //        template <typename Stream>
      //        void
      //        print (Stream     &stream,
      //               const bool  print_independent_symbols = false,
      //               const bool  print_dependent_functions = false,
      //               const bool  print_cse_reductions      = false) const
      //        {
      //          // Check to see if CSE has been performed
      //          if (print_cse_reductions && cse.executed())
      //            cse.print(stream);
      //        }
      //
      //        void
      //        bvisit (const SE::Integer &x);
      //
      //        void
      //        bvisit (const SE::Rational &x);
      //
      //        void
      //        bvisit (const SE::RealDouble &x);
      //
      //        void
      //        bvisit (const SE::Constant &x);
      //
      //        void
      //        bvisit (const SE::Symbol &symb);
      //
      //        void
      //        bvisit (const ADNumberWrapper<ADNumberType> &x);
      //
      //        void
      //        bvisit (const SE::Basic &x);
      //
      //        void
      //        bvisit (const SE::Add &op);
      //
      //        void
      //        bvisit (const SE::Mul &op);
      //
      //        void
      //        bvisit (const SE::Pow &op);
      //
      //        void
      //        bvisit (const SE::Log &op);
      //
      //        void
      //        bvisit (const SE::Sin &op);
      //
      //        void
      //        bvisit (const SE::Cos &op);
      //
      //        void
      //        bvisit (const SE::Tan &op);
      //
      //        void
      //        bvisit (const SE::Csc &op);
      //
      //        void
      //        bvisit (const SE::Sec &op);
      //
      //        void
      //        bvisit (const SE::Cot &op);
      //
      //        void
      //        bvisit (const SE::ASin &op);
      //
      //        void
      //        bvisit (const SE::ACos &op);
      //
      //        void
      //        bvisit (const SE::ATan &op);
      //
      //        void
      //        bvisit (const SE::ACsc &op);
      //
      //        void
      //        bvisit (const SE::ASec &op);
      //
      //        void
      //        bvisit (const SE::ACot &op);
      //
      //        void
      //        bvisit (const SE::Sinh &op);
      //
      //        void
      //        bvisit (const SE::Cosh &op);
      //
      //        void
      //        bvisit (const SE::Tanh &op);
      //
      //        void
      //        bvisit (const SE::Csch &op);
      //
      //        void
      //        bvisit (const SE::Sech &op);
      //
      //        void
      //        bvisit (const SE::Coth &op);
      //
      //        void
      //        bvisit (const SE::ASinh &op);
      //
      //        void
      //        bvisit (const SE::ACosh &op);
      //
      //        void
      //        bvisit (const SE::ATanh &op);
      //
      //        void
      //        bvisit (const SE::ACsch &op);
      //
      //        void
      //        bvisit (const SE::ASech &op);
      //
      //        void
      //        bvisit (const SE::ACoth &op);
      //
      //        void
      //        bvisit (const SE::Abs &op);
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperLambdaVisitor_Base
      //      (non-stacked) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for tapeless auto-differentiable types
      //       */
      //      template <typename ADNumberType>
      //      class ADNumberWrapperLambdaVisitor_Base<ADNumberType, typename
      //      std::enable_if<
      //        !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //        : public
      //        SE::BaseVisitor<ADNumberWrapperLambdaVisitor_Base<ADNumberType,
      //        typename std::enable_if<
      //        !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type> >
      //      {
      //      protected:
      //
      //        // The 'func' variable is assigned into at the very end of each
      //        visit()
      //        // methods below. The only place where these methods are called
      //        from is the
      //        // line 'b.accept(*this)' in apply() and the 'func' is
      //        immediately
      //        // returned. Thus no corruption can happen and apply() can be
      //        safely called
      //        // recursively.
      //
      //        /**
      //         * A typedef for the evaluation function signature
      //         */
      //        typedef std::function<ADNumberType(const ADNumberType *values)>
      //        FunctionType;
      //
      //        /**
      //         * The set of symbols that represent the independent variables
      //         */
      //        SE::vec_basic symbols;
      //
      //        /**
      //         * The set of constructed lambda functions that represent
      //         * evaluation of the dependent functions
      //         */
      //        std::vector<FunctionType> functions;
      //
      //        /**
      //         * An intermediate value used to build the vector of dependent
      //         functions.
      //         *
      //         * The @p func variable is assigned into at the very end of each visit()
      //         * methods below. The only place where these methods are called
      //         from is the
      //         * line 'b.accept(*this)' in apply() and the @p func is immediately
      //         * returned. Thus no corruption can happen and apply() can be
      //         safely called
      //         * recursively.
      //         */
      //        FunctionType func;
      //
      //        /**
      //         * A helper class that will be used to perform common
      //         * subexpression elimination
      //         */
      //        ADNumberWrapperCSELambdaVisitor<ADNumberType,FunctionType> cse;
      //
      //      public:
      //
      //        virtual ~ADNumberWrapperLambdaVisitor_Base () {}
      //
      //        void
      //        init (const SE::vec_basic &symbs,
      //              const SE::Basic     &func,
      //              const bool           use_cse = false);
      //
      //        void
      //        init(const SE::vec_basic &symbs,
      //             const SE::vec_basic &funcs,
      //             const bool           use_cse = false);
      //
      //        FunctionType
      //        apply (const SE::Basic &x);
      //
      //        ADNumberType
      //        call (const std::vector<ADNumberType> &values);
      //
      //        void
      //        call(ADNumberType       *outputs,
      //             const ADNumberType *values);
      //
      //        template <typename Stream>
      //        void
      //        print (Stream     &stream,
      //               const bool  print_independent_symbols = false,
      //               const bool  print_dependent_functions = false,
      //               const bool  print_cse_reductions      = false) const
      //        {
      //          // Check to see if CSE has been performed
      //          if (print_cse_reductions && cse.executed())
      //            cse.print(stream);
      //        }
      //
      //        void
      //        bvisit (const SE::Integer &x);
      //
      //        void
      //        bvisit (const SE::Rational &x);
      //
      //        void
      //        bvisit (const SE::RealDouble &x);
      //
      //        void
      //        bvisit (const SE::Constant &x);
      //
      //        void
      //        bvisit (const SE::Symbol &symb);
      //
      //        void
      //        bvisit (const ADNumberWrapper<ADNumberType> &x);
      //
      //        void
      //        bvisit (const SE::Basic &x);
      //
      //        void
      //        bvisit (const SE::Add &op);
      //
      //        void
      //        bvisit (const SE::Mul &op);
      //
      //        void
      //        bvisit (const SE::Pow &op);
      //
      //        void
      //        bvisit (const SE::Log &op);
      //
      //        void
      //        bvisit (const SE::Sin &op);
      //
      //        void
      //        bvisit (const SE::Cos &op);
      //
      //        void
      //        bvisit (const SE::Tan &op);
      //
      //        void
      //        bvisit (const SE::Csc &op);
      //
      //        void
      //        bvisit (const SE::Sec &op);
      //
      //        void
      //        bvisit (const SE::Cot &op);
      //
      //        void
      //        bvisit (const SE::ASin &op);
      //
      //        void
      //        bvisit (const SE::ACos &op);
      //
      //        void
      //        bvisit (const SE::ATan &op);
      //
      //        void
      //        bvisit (const SE::ACsc &op);
      //
      //        void
      //        bvisit (const SE::ASec &op);
      //
      //        void
      //        bvisit (const SE::ACot &op);
      //
      //        void
      //        bvisit (const SE::Sinh &op);
      //
      //        void
      //        bvisit (const SE::Cosh &op);
      //
      //        void
      //        bvisit (const SE::Tanh &op);
      //
      //        void
      //        bvisit (const SE::Csch &op);
      //
      //        void
      //        bvisit (const SE::Sech &op);
      //
      //        void
      //        bvisit (const SE::Coth &op);
      //
      //        void
      //        bvisit (const SE::ASinh &op);
      //
      //        void
      //        bvisit (const SE::ACosh &op);
      //
      //        void
      //        bvisit (const SE::ATanh &op);
      //
      //        void
      //        bvisit (const SE::ACsch &op);
      //
      //        void
      //        bvisit (const SE::ASech &op);
      //
      //        void
      //        bvisit (const SE::ACoth &op);
      //
      //        void
      //        bvisit (const SE::Abs &op);
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperLambdaVisitor
      //      (taped, real) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for taped real auto-differentiable types
      //       */
      //      template<typename ADNumberType>
      //      class ADNumberWrapperLambdaVisitor<ADNumberType,
      //        typename std::enable_if<
      //        AD::is_real_valued_ad_number<ADNumberType>::value &&
      //        internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //        : public
      //        SE::BaseVisitor<ADNumberWrapperLambdaVisitor<ADNumberType,
      //          typename std::enable_if<
      //          AD::is_real_valued_ad_number<ADNumberType>::value &&
      //          internal::UseStackedLambdaVisitor<ADNumberType>::value
      //          >::type>,
      //          ADNumberWrapperLambdaVisitor_Base<ADNumberType> >
      //      {
      //        typedef typename AD::ADNumberTraits<ADNumberType>::real_type
      //        ADRealType; typedef typename
      //        AD::ADNumberTraits<ADNumberType>::complex_type ADComplexType;
      //
      //        static_assert(AD::is_ad_number<ADNumberType>::value == true,
      //                      "Only auto-differentiable numbers can be used with
      //                      this class."
      //                     );
      //        static_assert(boost::is_complex<ADComplexType>::value == true,
      //                      "Expected a complex number type from
      //                      AD::NumberTraits."
      //                     );
      //
      //      protected:
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::FunctionType;
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::ReturnValueType;
      ////        using
      /// ADNumberWrapperLambdaVisitor_Base<ADNumberType>::has_sensitivities;
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::ValuePoolType;
      //
      //      public:
      //        // Classes not implemented are
      //        // Subs, UpperGamma, LowerGamma, Dirichlet_eta, Zeta
      //        // LeviCivita, KroneckerDelta, FunctionSymbol, LambertW
      //        // Derivative, Complex, ComplexDouble, ComplexMPC
      //
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::bvisit;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::init;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::call;
      //
      //        virtual ~ADNumberWrapperLambdaVisitor ();
      //
      //        void
      //        bvisit (const SE::ATan2 &op);
      //
      //        void
      //        bvisit (const SE::Erf &op);
      //
      //        void
      //        bvisit (const SE::Erfc &op);
      //
      //        void
      //        bvisit (const SE::Max &op);
      //
      //        void
      //        bvisit (const SE::Min &op);
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperLambdaVisitor
      //      (non-stacked, real) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for tapeless real auto-differentiable types
      //       */
      //      template<typename ADNumberType>
      //      class ADNumberWrapperLambdaVisitor<ADNumberType,
      //        typename std::enable_if<
      //        AD::is_real_valued_ad_number<ADNumberType>::value &&
      //        !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //        : public
      //        SE::BaseVisitor<ADNumberWrapperLambdaVisitor<ADNumberType,
      //          typename std::enable_if<
      //          AD::is_real_valued_ad_number<ADNumberType>::value &&
      //          !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //          >::type>,
      //          ADNumberWrapperLambdaVisitor_Base<ADNumberType> >
      //      {
      //        typedef typename AD::ADNumberTraits<ADNumberType>::real_type
      //        ADRealType; typedef typename
      //        AD::ADNumberTraits<ADNumberType>::complex_type ADComplexType;
      //
      //        static_assert(AD::is_ad_number<ADNumberType>::value == true,
      //                      "Only auto-differentiable numbers can be used with
      //                      this class."
      //                     );
      //        static_assert(boost::is_complex<ADComplexType>::value == true,
      //                      "Expected a complex number type from
      //                      AD::NumberTraits."
      //                     );
      //
      //      protected:
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::FunctionType;
      //
      //      public:
      //        // Classes not implemented are
      //        // Subs, UpperGamma, LowerGamma, Dirichlet_eta, Zeta
      //        // LeviCivita, KroneckerDelta, FunctionSymbol, LambertW
      //        // Derivative, Complex, ComplexDouble, ComplexMPC
      //
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::bvisit;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::init;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::call;
      //
      //        virtual ~ADNumberWrapperLambdaVisitor ();
      //
      //        void
      //        bvisit (const SE::ATan2 &op);
      //
      //        void
      //        bvisit (const SE::Erf &op);
      //
      //        void
      //        bvisit (const SE::Erfc &op);
      //
      //        void
      //        bvisit (const SE::Max &op);
      //
      //        void
      //        bvisit (const SE::Min &op);
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperLambdaVisitor
      //      (taped, complex) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for taped complex auto-differentiable types
      //       */
      //      template<typename ADNumberType>
      //      class ADNumberWrapperLambdaVisitor<ADNumberType,
      //        typename std::enable_if<
      //        AD::is_complex_valued_ad_number<ADNumberType>::value &&
      //        internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //        : public SE::BaseVisitor<
      //        ADNumberWrapperLambdaVisitor<ADNumberType,
      //        typename std::enable_if<
      //        AD::is_complex_valued_ad_number<ADNumberType>::value &&
      //        internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>,
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType> >
      //      {
      //        typedef typename AD::ADNumberTraits<ADNumberType>::real_type
      //        ADRealType; typedef typename
      //        AD::ADNumberTraits<ADNumberType>::complex_type ADComplexType;
      //
      //        static_assert(AD::is_ad_number<ADNumberType>::value == true,
      //                      "Only auto-differentiable numbers can be used with
      //                      this class."
      //                     );
      //        static_assert(boost::is_complex<ADComplexType>::value == true,
      //                      "Expected a complex number type from
      //                      AD::NumberTraits."
      //                     );
      //
      //      protected:
      ////    typedef typename
      /// ADNumberWrapperLambdaVisitor_Base<ADNumberType>::ReturnValueType
      /// ReturnValueType;
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::FunctionType;
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::ReturnValueType;
      ////        using
      /// ADNumberWrapperLambdaVisitor_Base<ADNumberType>::has_sensitivities;
      //
      //      public:
      //        // Classes not implemented are
      //        // ATan2, Erf, Erfc, Max, Min,
      //        // Subs, UpperGamma, LowerGamma, Dirichlet_eta, Zeta
      //        // LeviCivita, KroneckerDelta, FunctionSymbol, LambertW
      //        // Derivative, ComplexMPC
      //
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::bvisit;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::init;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::call;
      //
      //        virtual ~ADNumberWrapperLambdaVisitor ();
      //
      //        void
      //        bvisit (const SE::Complex &x);
      //
      //        void
      //        bvisit (const SE::ComplexDouble &x);
      //      };
      //
      //
      //      /* --------------------------- ADNumberWrapperLambdaVisitor
      //      (non-stacked, complex) ------------------------- */
      //
      //
      //      /**
      //       * Class definition for tapeless complex auto-differentiable types
      //       */
      //      template<typename ADNumberType>
      //      class ADNumberWrapperLambdaVisitor<ADNumberType,
      //        typename std::enable_if<
      //        AD::is_complex_valued_ad_number<ADNumberType>::value &&
      //        !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>
      //        : public SE::BaseVisitor<
      //        ADNumberWrapperLambdaVisitor<ADNumberType,
      //        typename std::enable_if<
      //        AD::is_complex_valued_ad_number<ADNumberType>::value &&
      //        !internal::UseStackedLambdaVisitor<ADNumberType>::value
      //        >::type>,
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType> >
      //      {
      //        typedef typename AD::ADNumberTraits<ADNumberType>::real_type
      //        ADRealType; typedef typename
      //        AD::ADNumberTraits<ADNumberType>::complex_type ADComplexType;
      //
      //        static_assert(AD::is_ad_number<ADNumberType>::value == true,
      //                      "Only auto-differentiable numbers can be used with
      //                      this class."
      //                     );
      //        static_assert(boost::is_complex<ADComplexType>::value == true,
      //                      "Expected a complex number type from
      //                      AD::NumberTraits."
      //                     );
      //
      //      protected:
      //        using typename
      //        ADNumberWrapperLambdaVisitor_Base<ADNumberType>::FunctionType;
      //
      //      public:
      //        // Classes not implemented are
      //        // ATan2, Erf, Erfc, Max, Min,
      //        // Subs, UpperGamma, LowerGamma, Dirichlet_eta, Zeta
      //        // LeviCivita, KroneckerDelta, FunctionSymbol, LambertW
      //        // Derivative, ComplexMPC
      //
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::bvisit;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::init;
      //        using ADNumberWrapperLambdaVisitor_Base<ADNumberType>::call;
      //
      //        virtual ~ADNumberWrapperLambdaVisitor ();
      //
      //        void
      //        bvisit (const SE::Complex &x);
      //
      //        void
      //        bvisit (const SE::ComplexDouble &x);
      //      };


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
      return internal::make_symengine_rcp(0.0);
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
          return internal::make_symengine_rcp(0.0);

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
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Number>
    LinearOperatorWrapper<Range, Domain, Payload>::rdiv(
      const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&other));
      return internal::make_symengine_rcp(0.0);
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
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::exp(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::log(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sin(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::cos(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }

    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::tan(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::cot(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::csc(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sec(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asin(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acos(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::atan(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acsc(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asec(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acot(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sinh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::cosh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::tanh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::csch(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::sech(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::coth(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asinh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acosh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::atanh(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acsch(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::asech(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::acoth(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::abs(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::floor(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::ceiling(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::erf(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::erfc(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


    template <typename Range, typename Domain, typename Payload>
    SE::RCP<const SE::Basic>
    LinearOperatorWrapper<Range, Domain, Payload>::gamma(
      const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLinearOperatorWrapper(&x));
      return internal::make_symengine_rcp(0.0);
    }


#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE


#endif
