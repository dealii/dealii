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

#ifndef dealii_differentiation_sd_symengine_la_vector_types_h
#define dealii_differentiation_sd_symengine_la_vector_types_h

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
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/la_vector.h>
#  include <deal.II/lac/linear_operator.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector.h>

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


    namespace internal
    {
      namespace
      {
        // Something to help us work out if an template parameter
        // refers to a LA::Vector. Since we only support LinearOperators
        // as sparse matrix wrappers within this framework, we can simply
        // check to see if the input type is compatible with the
        // initialization check associated with this class.
        //
        // https://stackoverflow.com/a/44229779
        template <class T, std::size_t = sizeof(T)>
        std::true_type
        vector_template_valid_impl(T *);

        std::false_type
        vector_template_valid_impl(...) __attribute__((unused));

        template <class V, typename T = void>
        struct SDIsVector : std::false_type
        {};

        template <class V, typename T = void>
        struct SDIsBlockVector : std::false_type
        {};

        // Need to avoid instantiating Vector<ADNumberType>, which is blocked by
        // a static assert
        template <class T>
        struct SDIsVector<T,
                          typename std::enable_if<
                            Differentiation::AD::is_ad_number<T>::value>::type>
          : std::false_type
        {};

        template <class T>
        struct SDIsBlockVector<
          T,
          typename std::enable_if<
            Differentiation::AD::is_ad_number<T>::value>::type>
          : std::false_type
        {};

        // Now that we've filtered out accidental instantiations for AD numbers,
        // we can go ahead and use the built-in check for a block vector. Again,
        // we must be careful not to instantiate the Vector classes for AD
        // types. We do this by adding the value check only after we have
        // enabled / disabled this specialization based on the potentially
        // invalid input type.
        template <class V>
        struct SDIsVector<V,
                          typename std::enable_if<
                            !Differentiation::AD::is_ad_number<V>::value>::type>
        {
          static const bool value = IsBlockVector<BlockVector<V>>::value;
        };

        template <class V>
        struct SDIsBlockVector<
          V,
          typename std::enable_if<
            !Differentiation::AD::is_ad_number<V>::value>::type>
        {
          static const bool value = IsBlockVector<V>::value;
        };

        template <class T>
        using vector_template_valid = decltype(
          typename std::enable_if<!std::is_arithmetic<T>::value>::type(),
          typename std::enable_if<SDIsVector<T>::value ||
                                  SDIsBlockVector<T>::value>::type(),
          vector_template_valid_impl(
            std::declval<dealii::internal::LinearOperatorImplementation::
                           ReinitHelper<T> *>()));
      } // namespace

      template <
        typename VectorType,
        typename = typename std::enable_if<
          internal::vector_template_valid<VectorType>::value == true>::type>
      static SE::RCP<const SE::Basic>
      make_symengine_rcp(const VectorType &vec, const std::string &name);

      template <typename VectorType>
      const VectorType &
      evaluate_symengine_number(const VectorType &vec);

      template <typename VectorType>
      bool
      is_a_LAVectorWrapper(const SE::RCP<const SE::Basic> &value);
    } // namespace internal



    /**
     * A number class to hold sparse matrices wrapped as LinearOperators.
     *
     * @author Jean-Paul Pelteret, 2018
     */
    template <typename VectorType_>
    class LAVectorWrapper : public SE::NumberWrapper, public SE::Evaluate
    {
      static_assert(
        internal::vector_template_valid<VectorType_>::value == true,
        "Only linear algebra vectors compatible with LinearOperators "
        "can be used with this class.");

    public:
      typedef VectorType_ VectorType;


      /**
       * Class constructor.
       *
       * This is marked as explicit so that there are no accidental conversions
       * of compatible numbers to this class type.
       */
      explicit LAVectorWrapper(const VectorType &vec, const std::string &name);

      /**
       * Class destructor.
       */
      virtual ~LAVectorWrapper();

      /**
       * @name Value retrieval
       */
      //@{

      /**
       * Returns a copy of the @p value of the underlying vector.
       */
      inline VectorType
      get_value()
      {
        return this->value;
      }

      /**
       * Returns the underlying vector by reference.
       */
      inline const VectorType &
      get_value() const
      {
        return this->value;
      }

      /**
       * Returns the name associated with the underlying vector.
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
       * The value of the linear algebra vector that is being stored.
       */
      VectorType value;

      /**
       * A symbolic name to define the vector that this
       * symbol represents
       */
      std::string name;

    }; // class LAVectorWrapper


    /* ------------------- inline and template functions  ----------------- */


#  ifndef DOXYGEN


    namespace internal
    {
      template <typename VectorType, typename>
      inline SE::RCP<const SE::Basic>
      make_symengine_rcp(const VectorType &vec, const std::string &name)
      {
        return SE::make_rcp<const LAVectorWrapper<VectorType>>(vec, name);
      }


      template <typename VectorType>
      const VectorType &
      evaluate_symengine_number(const VectorType &vec)
      {
        return vec.get_value();
      }



      template <typename VectorType>
      bool
      is_a_LAVectorWrapper(const SE::Basic &value)
      {
        return dynamic_cast<const LAVectorWrapper<VectorType> *const>(&value);
      }



      /**
       * Evaluate a SymEngine number to deal.II Vector type supported by
       * the LinearOperator class.
       */
      template <typename ReturnType>
      struct Evaluator<
        ReturnType,
        typename std::enable_if<
          internal::vector_template_valid<ReturnType>::value == true>::type>
      {
        static ReturnType
        value(const SE::RCP<const SE::Basic> &value)
        {
          if (internal::is_a_LAVectorWrapper<ReturnType>(*value))
            {
              return SE::down_cast<const LAVectorWrapper<ReturnType> &>(*value)
                .get_value();
            }
          else
            {
              AssertThrow(false, ExcNotImplemented());
              return ReturnType();
            }
        }
      };

    } // namespace internal


    /**
     * An error declaration stating that the object
     * is not an LAVectorWrapper type.
     */
    DeclExceptionMsg(ExcNotLAVectorWrapper,
                     "Object not of type LAVectorWrapper<VectorType>.");

    /**
     * An error declaration stating that the object
     * is not an LAVectorWrapper type.
     */
    DeclException1(
      ExcOpNotDefinedByLAVectorWrapper,
      const SE::Basic *,
      "Operation is not defined by the type LAVectorWrapper<VectorType>: " +
        arg1->__str__());


    template <typename VectorType>
    LAVectorWrapper<VectorType>::LAVectorWrapper(const VectorType & vec,
                                                 const std::string &name)
    {
      this->value = vec;
      this->name  = name;
    }


    template <typename VectorType>
    LAVectorWrapper<VectorType>::~LAVectorWrapper()
    {}


    template <typename VectorType>
    inline SE::hash_t
    LAVectorWrapper<VectorType>::__hash__() const
    {
      return std::hash<std::string>{}(__str__());
    }


    template <typename VectorType>
    bool
    LAVectorWrapper<VectorType>::__eq__(const SE::Basic &other) const
    {
      if (internal::is_a_LAVectorWrapper<VectorType>(other))
        {
          // With linear operators, we cannot make direct comparison of
          // any values, so instead we assume that we can do something
          // sensible with their memory addresses.
          const LAVectorWrapper<VectorType> &op_other =
            SE::down_cast<const LAVectorWrapper<VectorType> &>(other);
          return &(this->value) == &(op_other.value);
        }
      return false;
    }


    template <typename VectorType>
    int
    LAVectorWrapper<VectorType>::compare(const SE::Basic &other) const
    {
      Assert((internal::is_a_LAVectorWrapper<VectorType>(other)),
             ExcNotLAVectorWrapper());

      // With vectors, we could check the L2-norm
      // and sort them in this manner, but this is
      // an expensive operation. So we just sort by
      // memory address instead.
      const LAVectorWrapper<VectorType> &op_other =
        SE::down_cast<const LAVectorWrapper<VectorType> &>(other);
      if (&(this->value) == &(op_other.value))
        return 0;
      return (&(this->value) < &(op_other.value)) ? -1 : 1;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_exact() const
    {
      return false;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_positive() const
    {
      return false;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_negative() const
    {
      return false;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_zero() const
    {
      return this->value.l2_norm() == 0.0;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_one() const
    {
      return false;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_minus_one() const
    {
      return false;
    }


    template <typename VectorType>
    inline bool
    LAVectorWrapper<VectorType>::is_complex() const
    {
      return false;
    }


    template <typename VectorType>
    SE::Evaluate &
    LAVectorWrapper<VectorType>::get_eval() const
    {
      // This class derives from SE::Evaluate, so why shouldn't it
      // be able to evaluate itself?
      return const_cast<LAVectorWrapper<VectorType> &>(*this);
    }


    template <typename VectorType>
    std::string
    LAVectorWrapper<VectorType>::__str__() const
    {
      std::stringstream ss;
      ss << "Vec[";
      ss << this->name;
      ss << "]";
      return ss.str();
    }



    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::eval(long /*bits*/) const
    {
      // Used in the EvalXYZVisitor classes
      // Don't know why we'd need this when we already
      // have the get_eval() virtual function
      AssertThrow(false, ExcNotImplemented());

      // Needed for SymEngine::StrPrinter
      return SE::make_rcp<const LAVectorWrapper<VectorType>>(get_value(),
                                                             get_name());
    }


    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::add(const SE::Number &other) const
    {
      if (internal::is_a_LAVectorWrapper<VectorType>(other))
        {
          const LAVectorWrapper<VectorType> &other_op =
            SE::down_cast<const LAVectorWrapper<VectorType> &>(other);
          VectorType result(this->value);
          result += other_op.value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, "(" + this->name + "+" + other_op.name + ")");
        }
      // Due to the way SymEngine works, we may be adding some constants when
      // expanding expressions. The following conditions take these into
      // consideration and return some reasonable result.
      else if (SE::is_a<SE::Integer>(other) && other.is_zero())
        {
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(this->value,
                                                                 this->name);
        }

      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return other.add(*this);
    }

    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::sub(const SE::Number &other) const
    {
      if (internal::is_a_LAVectorWrapper<VectorType>(other))
        {
          const LAVectorWrapper<VectorType> &other_op =
            SE::down_cast<const LAVectorWrapper<VectorType> &>(other);
          VectorType result(this->value);
          result -= other_op.value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, "(" + this->name + "-" + other_op.name + ")");
        }
      // Due to the way SymEngine works, we may be adding some constants when
      // expanding expressions. The following conditions take these into
      // consideration and return some reasonable result.
      else if (SE::is_a<SE::Integer>(other) && other.is_zero())
        {
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(this->value,
                                                                 this->name);
        }

      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return other.rsub(*this);
    }

    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::rsub(const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::mul(const SE::Number &other) const
    {
      if (internal::is_a_LAVectorWrapper<VectorType>(other))
        {
          const LAVectorWrapper<VectorType> &other_op =
            SE::down_cast<const LAVectorWrapper<VectorType> &>(other);
          // Vector inner product returns a number
          return internal::make_symengine_rcp(this->value * other_op.value);
        }
      else if (other.is_one())
        {
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(this->value,
                                                                 this->name);
        }
      else if (SE::is_a<SE::RealDouble>(other))
        {
          const SE::RealDouble &other_nmbr =
            SE::down_cast<const SE::RealDouble &>(other);
          const double value = internal::evaluate_symengine_number(other_nmbr);
          VectorType   result(this->value);
          result *= value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, this->name + "*" + dealii::Utilities::to_string(value));
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
          const signed long int value =
            internal::evaluate_symengine_number(other_nmbr);
          VectorType result(this->value);
          result *= value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, this->name + "*" + dealii::Utilities::to_string(value));
        }
      else if (SE::is_a<SE::Rational>(other))
        {
          const SE::Rational &other_nmbr =
            SE::down_cast<const SE::Rational &>(other);
          const double value = internal::evaluate_symengine_number(other_nmbr);
          VectorType   result(this->value);
          result *= value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, this->name + "*" + dealii::Utilities::to_string(value));
        }

      //        // We can guarentee that other number types don't know anything
      //        about
      //        // this object, so we can't do anything more...
      //        AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));

      // At this point we would typically throw an error, but in the case that
      // the "other" parameter is a LinearOperatorWrapper then we'd actually
      // want it to perform the multiplication operation on this vector. These
      // situations arise because SymEngine is unaware that these operations are
      // not commutative, and so it happily rearranges the order of the
      // variables thinking that this produces a valid operation. See, for
      // example, symengine/symengine_linear_operator_wrapper_03.cc We construct
      // an operation symb_lo_A_inv*(symb_la_u+symb_la_v) which it internally
      // rearranges to (symb_la_u+symb_la_v)*symb_lo_A_inv and, therefore, the
      // RHS object is a LAVectorWrapper, i.e. this object.
      return other.mul(*this);
    }


    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::div(const SE::Number &other) const
    {
      if (internal::is_a_LAVectorWrapper<VectorType>(other))
        {
          AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
          return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
        }
      else if (other.is_one())
        {
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(this->value,
                                                                 this->name);
        }
      else if (SE::is_a<SE::RealDouble>(other))
        {
          const SE::RealDouble &other_nmbr =
            SE::down_cast<const SE::RealDouble &>(other);
          const double value = internal::evaluate_symengine_number(other_nmbr);
          VectorType   result(this->value);
          result /= value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, this->name + "/" + dealii::Utilities::to_string(value));
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
          const signed long int value =
            internal::evaluate_symengine_number(other_nmbr);
          VectorType result(this->value);
          result /= value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, this->name + "/" + dealii::Utilities::to_string(value));
        }
      else if (SE::is_a<SE::Rational>(other))
        {
          const SE::Rational &other_nmbr =
            SE::down_cast<const SE::Rational &>(other);
          const double value = internal::evaluate_symengine_number(other_nmbr);
          VectorType   result(this->value);
          result /= value;
          return SE::make_rcp<const LAVectorWrapper<VectorType>>(
            result, this->name + "/" + dealii::Utilities::to_string(value));
        }

      // We can guarentee that other number types don't know anything about
      // this object, so we can't do anything more...
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return other.mul(*this);
    }


    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::rdiv(const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::pow(const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Number>
    LAVectorWrapper<VectorType>::rpow(const SE::Number &other) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&other));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::exp(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::log(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::sin(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::cos(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }

    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::tan(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::cot(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::csc(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::sec(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::asin(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::acos(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::atan(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::acsc(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::asec(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::acot(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::sinh(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::cosh(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::tanh(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::csch(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::sech(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::coth(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::asinh(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::acosh(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::atanh(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::acsch(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::asech(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::acoth(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::abs(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::floor(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::ceiling(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::erf(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::erfc(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


    template <typename VectorType>
    SE::RCP<const SE::Basic>
    LAVectorWrapper<VectorType>::gamma(const SE::Basic &x) const
    {
      AssertThrow(false, ExcOpNotDefinedByLAVectorWrapper(&x));
      return SE::rcp_dynamic_cast<const SE::Number>(SD::Expression(0.0).get_RCP());
    }


#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE


#endif
