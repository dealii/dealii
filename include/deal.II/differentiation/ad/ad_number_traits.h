// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_differentiation_ad_ad_number_traits_h
#define dealii_differentiation_ad_ad_number_traits_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <deal.II/differentiation/ad/ad_number_types.h>

#include <boost/type_traits.hpp>

#include <complex>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace AD
  {
    /**
     * A number traits class to help describe some characteristic
     * information about auto-differentiable numbers.
     *
     * @tparam ScalarType A real or complex floating point number.
     * @tparam ADNumberTypeCode An enumeration specifying the type
     *         code for the supported auto-differentiable counterpart
     *         to the given @p ScalarType.
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename ScalarType,
              enum NumberTypes ADNumberTypeCode,
              typename T = void>
    struct NumberTraits;



    /**
     * A number traits class to help describe some characteristic
     * information about auto-differentiable numbers. This class is only
     * specialized (and, therefore, only made constructible) if the
     * @p ADNumberType template parameter is indeed a supported
     * auto-differentiable number.
     *
     * @tparam ADNumberType A type corresponding to a supported
     *         auto-differentiable number.
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename ADNumberType, typename T = void>
    struct ADNumberTraits;


    /**
     * This namespace defines the classes that help provide a unified interface
     * to all auto-differentiation numbers.
     */
    namespace internal
    {
      // The following three classes, namely ADNumberInfoFromEnum, Marking, and
      // ExtractData, are those that need to be implemented for each new
      // auto-differentiable number type. This information is then used by
      // NumberTraits and ADNumberTraits to provide a uniform interface, as used
      // by our drivers, to the underlying number types.


      /**
       * A struct that defines some fundamental information about a
       * auto-differentiable number based on the @p ScalarType and the
       * AD-enumeration selected by @p ADNumberTypeCode. This information is
       * used in other convenience classes and templated functions to
       * automatically determine information about the auto-differentiable
       * number that has been selected to wrap the @p ScalarType.
       *
       * The specializations of this class have to implement the following
       * member data and type definitions:
       * @code
       *
       *   // State whether the auto-differentiable number uses taping or not.
       *   static const bool             is_taped;
       *   // The real-type for the auto-differentiable number
       *   using real_type = <ADNumberType>;
       *   // The type of number returned when taking the first derivative of the @p real_type.
       *   using derivative_type = <Scalar/ADNumberType>;
       *   // The number of derivative levels computable from the @p real_type.
       *   static const unsigned int     n_supported_derivative_levels;
       *
       * @endcode
       *
       * @author Jean-Paul Pelteret, 2017
       */
      template <typename ScalarType,
                enum NumberTypes ADNumberTypeCode,
                typename = void>
      struct ADNumberInfoFromEnum;


      /**
       * A struct to assist with the marking of AD numbers that represent
       * independent and dependent variables.
       *
       * The specializations of this class have to implement the following
       * member functions:
       * @code
       *
       *  // Initialize the state of an independent variable.
       *  static void
       *  independent_variable(const ScalarType   &in,
       *                       const unsigned int  index,
       *                       const unsigned int  n_independent_variables,
       *                       ADNumberType       &out)
       *
       *  // Initialize the state of a dependent variable.
       *  static void
       *  dependent_variable(ADNumberType       &out,
       *                     const ADNumberType &func);
       *
       * @endcode
       * where @p ADNumberType is the auto-differentiable number type and
       * @p ScalarType is its floating point counterpart.
       *
       * @tparam ADNumberType A type corresponding to a supported
       *         auto-differentiable number.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       *
       * @author Jean-Paul Pelteret, 2017
       */
      template <typename ADNumberType, typename T = void>
      struct Marking;


      /**
       * A struct to help extract certain information associated with
       * auto-differentiable numbers.
       *
       * The specializations of this class have to implement the following
       * member data and type definitions:
       * @code
       *
       *   // Extract the real scalar value.
       *   static scalar_type
       *   value (const ADNumberType &x);
       *
       *   // Extract the number of directional derivatives.
       *   static unsigned int
       *   n_directional_derivatives (const ADNumberType &x);
       *
       *   // Extract the directional derivative in the specified @p direction.
       *   static derivative_type
       *   directional_derivative (const ADNumberType &x,
       *                           const unsigned int  direction);
       *
       * @endcode
       * where @p ADNumberType is the auto-differentiable number type.
       *
       * @tparam ADNumberType A type corresponding to a supported
       *         auto-differentiable number.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       *
       * @author Jean-Paul Pelteret, 2017
       */
      template <typename ADNumberType, typename T = void>
      struct ExtractData;


      /**
       * A struct that checks that the data expected to be stored in a
       * specialization of the ADNumberInfoFromEnum struct has been supplied. By
       * default it is assumed that the input type does not satisfy the
       * necessary conditions to construct this class.
       *
       * @tparam ADNumberTrait A class that examined whether it contains the necessary
       *         information to satisfy the requirements for being an internally
       *         supported auto-differentiable number.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       *
       * @author Jean-Paul Pelteret, 2017
       */
      template <typename ADNumberTrait, typename T = void>
      struct HasRequiredADInfo;


      /**
       * Provide a convenience function to assist in the casting of some
       * number types to other number types. On top of the standard class
       * definition given in @p base/numbers.h , this extension allows the
       * conversion of automatic-differentiation numbers to generic floats.
       *
       * This is necessary because ADOL-C doesn't provide a convenient
       * way to convert from an @p ADNumberType to floats (@p T) other than
       * the real-type equivalent that its associated with. For Sacado, and
       * likely other AD number types,  the floating point value stored in
       * an @p ADNumberType must be extracted through some function that is
       * specific to each type of AD number. This requires some specialist
       * intervention to get at this data.
       *
       * @author Jean-Paul Pelteret, 2017
       */
      template <typename T>
      struct NumberType;


      /**
       * A small struct to remove the @p std::complex wrapper
       * around a number.
       *
       * @author Jean-Paul Pelteret, 2017
       */
      template <typename Number>
      struct RemoveComplexWrapper;

    } // namespace internal


    /**
     * A struct to indicate whether a given @p NumberType is a supported
     * auto-differentiable number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename NumberType>
    struct is_ad_number;


    /**
     * A struct to indicate whether a given @p NumberType is a taped
     * auto-differentiable number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename NumberType, typename = void>
    struct is_taped_ad_number;


    /**
     * A struct to indicate whether a given @p NumberType is a tapeless
     * auto-differentiable number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename NumberType, typename = void>
    struct is_tapeless_ad_number;


    /**
     * A struct to indicate whether a given @p NumberType is a real-valued
     * auto-differentiable number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename NumberType, typename = void>
    struct is_real_valued_ad_number;


    /**
     * A struct to indicate whether a given @p NumberType is a complex-valued
     * auto-differentiable number or not. By default, numbers are not
     * considered to have the necessary characteristics to fulfill this
     * condition.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename NumberType, typename = void>
    struct is_complex_valued_ad_number;

  } // namespace AD
} // namespace Differentiation


/* ----------- inline and template functions and specializations ----------- */


#ifndef DOXYGEN


namespace Differentiation
{
  namespace AD
  {
    namespace internal
    {
      template <typename ADNumberTrait, typename>
      struct HasRequiredADInfo : std::false_type
      {};


      /**
       * Specialization to detect whether the input AD number
       * is internally supported or not. In particular, we
       * check to see that it is not an floating point type, that it
       * has been assigned a type_code and has other basic
       * characteristics necessary for the internal interface of
       * the AD drivers.
       *
       * The implementation of this struct follows this suggestion:
       * https://stackoverflow.com/a/16000226
       */
      template <typename ADNumberTrait>
      struct HasRequiredADInfo<
        ADNumberTrait,
        decltype((void)ADNumberTrait::type_code,
                 (void)ADNumberTrait::is_taped,
                 (void)std::declval<typename ADNumberTrait::real_type>(),
                 (void)std::declval<typename ADNumberTrait::derivative_type>(),
                 void())>
        : std::conditional<
            std::is_floating_point<typename ADNumberTrait::real_type>::value,
            std::false_type,
            std::true_type>::type
      {};


      /**
       * A specialization for the information struct for arithmetic / floating
       * point numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::none,
        typename std::enable_if<
          std::is_floating_point<ScalarType>::value>::type>
      {
        static const bool is_taped                              = false;
        using real_type                                         = ScalarType;
        using derivative_type                                   = ScalarType;
        static const unsigned int n_supported_derivative_levels = 0;
      };


      /**
       * A dummy specialization for floating point numbers. This is helpful
       * for nested auto-differentiable numbers, where a recursive marking
       * mechanism can be employed (e.g. Sacado types).
       */
      template <typename ScalarType>
      struct Marking<ScalarType,
                     typename std::enable_if<
                       std::is_floating_point<ScalarType>::value>::type>
      {
        /**
         * Initialize the state of an independent variable.
         *
         * With a nested marking approach it is sometimes
         * necessary to initialise the value of an intermediate
         * value that may be a floating point number.
         */
        template <typename ADNumberType>
        static void
        independent_variable(const ScalarType &in,
                             const unsigned int,
                             const unsigned int,
                             ADNumberType &out)
        {
          out = in;
        }

        /*
         * Initialize the state of a dependent variable.
         */
        template <typename ADNumberType>
        static void
        dependent_variable(ADNumberType &, const ScalarType &)
        {
          AssertThrow(
            false,
            ExcMessage(
              "Floating point numbers cannot be marked as dependent variables."));
        }
      };


      /**
       * A specialization of the marking strategy for complex numbers.
       */
      template <typename ADNumberType>
      struct Marking<
        ADNumberType,
        typename std::enable_if<boost::is_complex<ADNumberType>::value>::type>
      {
        /*
         * Initialize the state of an independent variable.
         */
        template <typename ScalarType>
        static void
        independent_variable(const ScalarType &in,
                             const unsigned int,
                             const unsigned int,
                             ADNumberType &out)
        {
          AssertThrow(
            false,
            ExcMessage(
              "Marking for complex numbers has not yet been implemented."));
          out = in;
        }

        /*
         * Initialize the state of a dependent variable.
         */
        template <typename ScalarType>
        static void
        dependent_variable(ADNumberType &, const ScalarType &)
        {
          AssertThrow(
            false,
            ExcMessage(
              "Marking for complex numbers has not yet been implemented."));
        }
      };

    } // namespace internal


    template <typename NumberType, typename>
    struct is_taped_ad_number : std::false_type
    {};


    template <typename NumberType, typename>
    struct is_tapeless_ad_number : std::false_type
    {};


    template <typename NumberType, typename>
    struct is_real_valued_ad_number : std::false_type
    {};


    template <typename NumberType, typename>
    struct is_complex_valued_ad_number : std::false_type
    {};


    /**
     * We use the specialization of the HasRequiredADInfo struct
     * to ensure that only internally supported numbers are
     * considered AD numbers.
     */
    template <typename NumberType>
    struct is_ad_number
      : internal::HasRequiredADInfo<
          ADNumberTraits<typename std::decay<NumberType>::type>>
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (real or complex) taped auto-differentiable number.
     */
    template <typename NumberType>
    struct is_taped_ad_number<
      NumberType,
      typename std::enable_if<
        ADNumberTraits<typename std::decay<NumberType>::type>::is_taped>::type>
      : std::true_type
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (real or complex) tapeless auto-differentiable number.
     */
    template <typename NumberType>
    struct is_tapeless_ad_number<
      NumberType,
      typename std::enable_if<ADNumberTraits<
        typename std::decay<NumberType>::type>::is_tapeless>::type>
      : std::true_type
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (taped or tapeless) real-valued auto-differentiable
     * number.
     */
    template <typename NumberType>
    struct is_real_valued_ad_number<
      NumberType,
      typename std::enable_if<ADNumberTraits<
        typename std::decay<NumberType>::type>::is_real_valued>::type>
      : std::true_type
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (taped or tapeless) complex-valued auto-differentiable
     * number.
     */
    template <typename NumberType>
    struct is_complex_valued_ad_number<
      NumberType,
      typename std::enable_if<ADNumberTraits<
        typename std::decay<NumberType>::type>::is_complex_valued>::type>
      : std::true_type
    {};


    namespace internal
    {
      /**
       * Specialization of the selection struct which sets the type as being
       * one that is the value type of the underlying complex number.
       */
      template <typename Number>
      struct RemoveComplexWrapper
      {
        using type = Number;
      };


      /**
       * Specialization of the selection struct which sets the value type
       * to that resulting from the recursive removal of the complex
       * number wrapper.
       */
      template <typename Number>
      struct RemoveComplexWrapper<std::complex<Number>>
      {
        using type = typename RemoveComplexWrapper<Number>::type;
      };


      /**
       * A dummy specialization for floating point numbers. This is helpful
       * for nested auto-differentiable numbers, where a recursive marking
       * mechanism can be employed (e.g. Sacado types).
       */
      template <typename NumberType>
      struct ExtractData<NumberType,
                         typename std::enable_if<
                           std::is_floating_point<NumberType>::value>::type>
      {
        /**
         * Extract the floating point value.
         */
        static const NumberType &
        value(const NumberType &x)
        {
          return x;
        }


        /**
         * Extract the number of directional derivatives.
         */
        static unsigned int
        n_directional_derivatives(const NumberType &)
        {
          return 0;
        }


        /**
         * Extract the directional derivative in the specified @p direction.
         */
        static NumberType
        directional_derivative(const NumberType &, const unsigned int)
        {
          return 0.0;
        }
      };



      /**
       * A struct specialization to help extract certain information associated
       * with complex auto-differentiable numbers.
       */
      template <typename ADNumberType>
      struct ExtractData<std::complex<ADNumberType>>
      {
        static_assert(Differentiation::AD::is_ad_number<ADNumberType>::value,
                      "Expected an auto-differentiable number.");


        /**
         * Extract the floating point value.
         */
        static std::complex<typename ADNumberTraits<ADNumberType>::scalar_type>
        value(const std::complex<ADNumberType> &x)
        {
          return std::complex<
            typename ADNumberTraits<ADNumberType>::scalar_type>(
            ExtractData<ADNumberType>::value(x.real()),
            ExtractData<ADNumberType>::value(x.imag()));
        }


        /**
         * Extract the number of directional derivatives.
         */
        static unsigned int
        n_directional_derivatives(const std::complex<ADNumberType> &x)
        {
          return ExtractData<ADNumberType>::n_directional_derivatives(x.real());
        }


        /**
         * Extract the directional derivative in the specified @p direction.
         */
        static std::complex<
          typename ADNumberTraits<ADNumberType>::derivative_type>
        directional_derivative(const std::complex<ADNumberType> &x,
                               const unsigned int                direction)
        {
          return std::complex<
            typename ADNumberTraits<ADNumberType>::derivative_type>(
            ExtractData<ADNumberType>::directional_derivative(x.real(),
                                                              direction),
            ExtractData<ADNumberType>::directional_derivative(x.imag(),
                                                              direction));
        }
      };


      template <typename T>
      struct NumberType
      {
        /**
         * Standard number conversion
         */
        template <typename F>
        static auto
        value(const F &f,
              typename std::enable_if<!is_ad_number<F>::value>::type * =
                nullptr) -> decltype(dealii::internal::NumberType<T>::value(f))
        {
          // We call the other function defined in the numbers
          // header to take care of all of the usual cases.
          return dealii::internal::NumberType<T>::value(f);
        }

        /**
         * Conversion from an AD number to a scalar number.
         * We wish to ensure that @p T is a scalar type because,
         * in general, the conversion between AD numbers is either
         * non-trivial or an invalid operation.
         */
        template <typename F>
        static T
        value(const F &f,
              typename std::enable_if<is_ad_number<F>::value &&
                                      std::is_floating_point<T>::value>::type
                * = nullptr)
        {
          // We recursively call this function in case the AD number is a
          // nested one. The recursion ends when the extracted value is
          // a floating point number.
          return NumberType<T>::value(ExtractData<F>::value(f));
        }

        /**
         * Conversion from an AD number another AD number.
         *
         * Since this is the most generic case we'll assume that
         * the return type is constructible from the input type.
         */
        template <typename F>
        static T
        value(const F &f,
              typename std::enable_if<is_ad_number<F>::value &&
                                      is_ad_number<T>::value>::type * = nullptr)
        {
          return T(f);
        }
      };

      template <typename T>
      struct NumberType<std::complex<T>>
      {
        /**
         * Standard complex number conversion
         */
        template <typename F>
        static auto
        value(
          const F &f,
          typename std::enable_if<!is_ad_number<F>::value>::type * = nullptr)
          -> decltype(dealii::internal::NumberType<std::complex<T>>::value(f))
        {
          // We call the other function defined in the numbers
          // header to take care of all of the usual cases.
          return dealii::internal::NumberType<std::complex<T>>::value(f);
        }


        /**
         * Conversion from a complex AD number to another
         * complex number templated on a scalar number.
         */
        template <typename F>
        static std::complex<T>
        value(const F &f,
              typename std::enable_if<is_ad_number<F>::value &&
                                      std::is_floating_point<T>::value>::type
                * = nullptr)
        {
          // We recursively call this function in case the AD number is a
          // nested one. The recursion ends when the extracted value is
          // a floating point number.
          return std::complex<T>(
            NumberType<T>::value(ExtractData<F>::value(f)));
        }

        template <typename F>
        static std::complex<T>
        value(const std::complex<F> &f)
        {
          // Deal with the two parts of the input complex
          // number individually.
          return std::complex<T>(NumberType<T>::value(f.real()),
                                 NumberType<T>::value(f.imag()));
        }
      };

    } // namespace internal



    /**
     * Specialization of the general NumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     *
     * For each new @p ADNumberTypeCode enumeration, a data structure called
     * internal::ADNumberInfoFromEnum is to be defined, from which we require
     * the following basic information determine all of the required information
     * and type traits for our helper classes:
     *   - An alias called @p scalar_type, which defines the scalar or
     *     floating-point valued counterpart to the auto-differentiable number
     * type. This can be real or complex valued.
     *   - An alias called @p derivative_type, which defines the
     *     number type for directional derivatives.
     *   - A boolean called @p is_taped, which defines whether the auto-differentiable
     *     number is a taped number or not.
     *
     * This specific version specializes the generic case where the template
     * @p ScalarType represents a floating or complex number that is templated
     * on a floating point type.
     */
    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    struct NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>
    {
      /**
       * The type of taping used
       */
      static constexpr enum NumberTypes type_code = ADNumberTypeCode;

      // The clang compiler does not seem to like these
      // variables being defined as constant expressions
      // (the tests <adolc|sacado>/ad_number_traits_02 will
      // fail with linking errors). However, GCC complains
      // about the use of non-constant expressions in
      // std::conditional.
#  ifdef __clang__

      /**
       * A flag to indicate whether the number is of
       * the taped variety or not
       */
      static const bool is_taped;


      /**
       * A flag to indicate whether the number is of
       * the tapeless variety or not
       */
      static const bool is_tapeless;


      /**
       * A flag to indicate whether the number represents
       * a real value
       */
      static const bool is_real_valued;


      /**
       * A flag to indicate whether the number represents
       * a complex value
       */
      static const bool is_complex_valued;


      /**
       * The number of directional derivatives that can be
       * taken with this auto-differentiable number
       */
      static const unsigned int n_supported_derivative_levels;

#  else

      /**
       * A flag to indicate whether the number is of
       * the taped variety or not
       */
      static constexpr bool is_taped = internal::ADNumberInfoFromEnum<
        typename internal::RemoveComplexWrapper<ScalarType>::type,
        ADNumberTypeCode>::is_taped;


      /**
       * A flag to indicate whether the number is of
       * the tapeless variety or not
       */
      static constexpr bool is_tapeless =
        !(NumberTraits<ScalarType, ADNumberTypeCode>::is_taped);


      /**
       * A flag to indicate whether the number represents
       * a real value
       */
      static constexpr bool is_real_valued =
        (!boost::is_complex<ScalarType>::value);


      /**
       * A flag to indicate whether the number represents
       * a complex value
       */
      static constexpr bool is_complex_valued =
        !(NumberTraits<ScalarType, ADNumberTypeCode>::is_real_valued);


      /**
       * The number of directional derivatives that can be
       * taken with this auto-differentiable number
       */
      static constexpr unsigned int n_supported_derivative_levels =
        internal::ADNumberInfoFromEnum<
          typename internal::RemoveComplexWrapper<ScalarType>::type,
          ADNumberTypeCode>::n_supported_derivative_levels;

#  endif


      /**
       * Underlying floating point value type.
       * This could real-valued or complex-valued.
       */
      using scalar_type = ScalarType;


      /**
       * Type for real numbers
       */
      using real_type = typename internal::ADNumberInfoFromEnum<
        typename internal::RemoveComplexWrapper<ScalarType>::type,
        ADNumberTypeCode>::real_type;


      /**
       * Type for complex numbers
       */
      using complex_type = std::complex<real_type>;


      /**
       * The actual auto-differentiable number type
       */
      using ad_type = typename std::
        conditional<is_real_valued, real_type, complex_type>::type;

      /**
       * The actual auto-differentiable number directional derivative type
       */
      using derivative_type = typename std::conditional<
        is_real_valued,
        typename internal::ADNumberInfoFromEnum<
          typename internal::RemoveComplexWrapper<ScalarType>::type,
          ADNumberTypeCode>::derivative_type,
        std::complex<typename internal::ADNumberInfoFromEnum<
          typename internal::RemoveComplexWrapper<ScalarType>::type,
          ADNumberTypeCode>::derivative_type>>::type;


      /**
       * Extract the value of an auto-differentiable number
       */
      static scalar_type get_scalar_value(const ad_type &x)
      {
        // Some tricky conversion cases to consider here:
        // - Nested AD numbers
        // - std::complex<double> --> std::complex<float>
        //   e.g. when ScalarType = float and ADNumberTypeCode = adolc_taped
        // Therefore, we use the internal casting mechanism
        // provided by the internal::NumberType struct.
        return internal::NumberType<scalar_type>::value(
          internal::ExtractData<ad_type>::value(x));
      }


      /**
       * Extract the derivative value of an auto-differentiable number
       */
      static derivative_type get_directional_derivative(
        const ad_type &x, const unsigned int direction)
      {
        return internal::ExtractData<ad_type>::directional_derivative(
          x, direction);
      }


      /**
       * Extract the number of directional derivatives value tracked by
       * an auto-differentiable number
       */
      static unsigned int n_directional_derivatives(const ad_type &x)
      {
        return internal::ExtractData<ad_type>::n_directional_derivatives(x);
      }


      static_assert((is_real_valued == true ?
                       std::is_same<ad_type, real_type>::value :
                       std::is_same<ad_type, complex_type>::value),
                    "Incorrect template type selected for ad_type");

      static_assert((is_complex_valued == true ?
                       boost::is_complex<scalar_type>::value :
                       true),
                    "Expected a complex float_type");

      static_assert((is_complex_valued == true ?
                       boost::is_complex<ad_type>::value :
                       true),
                    "Expected a complex ad_type");
    };

#  ifdef __clang__

    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    const bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_taped =
      internal::ADNumberInfoFromEnum<
        typename internal::RemoveComplexWrapper<ScalarType>::type,
        ADNumberTypeCode>::is_taped;


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    const bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_tapeless =
      !(NumberTraits<ScalarType, ADNumberTypeCode>::is_taped);


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    const bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_real_valued =
      (!boost::is_complex<ScalarType>::value);


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    const bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_complex_valued =
      !(NumberTraits<ScalarType, ADNumberTypeCode>::is_real_valued);


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    const unsigned int NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::n_supported_derivative_levels =
      internal::ADNumberInfoFromEnum<
        typename internal::RemoveComplexWrapper<ScalarType>::type,
        ADNumberTypeCode>::n_supported_derivative_levels;

#  else

    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    constexpr bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_taped;


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    constexpr bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_tapeless;


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    constexpr bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_real_valued;


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    constexpr bool NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_complex_valued;


    template <typename ScalarType, enum NumberTypes ADNumberTypeCode>
    constexpr unsigned int NumberTraits<
      ScalarType,
      ADNumberTypeCode,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::n_supported_derivative_levels;

#  endif



    /**
     * A dummy specialization for floating point numbers. This is necessary to
     * deal with the special case of the ADNumberTypeCode that represents a
     * floating point number.
     *
     * This specific version specializes the generic case where the template
     * @p ScalarType represents a floating point type, or complex number that
     * is templated on a floating point type.
     */
    template <typename ScalarType>
    struct NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>
    {
      /**
       * The internal number type code.
       */
      static constexpr enum NumberTypes type_code = NumberTypes::none;

      // The clang compiler does not seem to like these
      // variables being defined as constant expressions
      // (the tests <adolc|sacado>/ad_number_traits_02 will
      // fail with linking errors). However, GCC complains
      // about the use of non-constant expressions in
      // std::conditional.
#  ifdef __clang__

      /**
       * A flag to indicate whether the number is of
       * the taped variety or not.
       */
      static const bool is_taped;


      /**
       * A flag to indicate whether the number is of
       * the tapeless variety or not.
       */
      static const bool is_tapeless;


      /**
       * A flag to indicate whether the number represents
       * a real value.
       */
      static const bool is_real_valued;


      /**
       * A flag to indicate whether the number represents
       * a complex value.
       */
      static const bool is_complex_valued;


      /**
       * The number of directional derivatives that can be
       * taken with this (floating point) number.
       */
      static const unsigned int n_supported_derivative_levels;

#  else

      /**
       * A flag to indicate whether the number is of
       * the taped variety or not.
       */
      static constexpr bool is_taped = false;


      /**
       * A flag to indicate whether the number is of
       * the tapeless variety or not.
       */
      static constexpr bool is_tapeless = false;


      /**
       * A flag to indicate whether the number represents
       * a real value.
       */
      static constexpr bool is_real_valued =
        (!boost::is_complex<ScalarType>::value);


      /**
       * A flag to indicate whether the number represents
       * a complex value.
       */
      static constexpr bool is_complex_valued = !is_real_valued;


      /**
       * The number of directional derivatives that can be
       * taken with this (floating point) number.
       */
      static constexpr unsigned int n_supported_derivative_levels = 0;

#  endif


      /**
       * Underlying floating point value type.
       * This could real-valued or complex-valued.
       */
      using scalar_type = ScalarType;


      /**
       * Type for real numbers.
       */
      using real_type =
        typename dealii::numbers::NumberTraits<scalar_type>::real_type;


      /**
       * Type for complex numbers.
       */
      using complex_type = std::complex<real_type>;


      /**
       * The actual auto-differentiable number type.
       */
      using ad_type = ScalarType;

      /**
       * The actual auto-differentiable number directional derivative type.
       */
      using derivative_type = ScalarType;


      /**
       * Extract the value of an auto-differentiable number.
       */
      static scalar_type
      get_scalar_value(const ad_type &x)
      {
        return x;
      }


      /**
       * Extract the derivative value of an auto-differentiable number.
       */
      static derivative_type
      get_directional_derivative(const ad_type &, const unsigned int)
      {
        Assert(
          false,
          ExcMessage(
            "Floating point/arithmetic numbers have no directional derivatives."));
        return derivative_type();
      }


      /**
       * Extract the number of directional derivatives value tracked by
       * an auto-differentiable number.
       */
      static unsigned int
      n_directional_derivatives(const ad_type &)
      {
        Assert(
          false,
          ExcMessage(
            "Floating point/arithmetic numbers have no directional derivatives."));
        return 0;
      }
    };

#  ifdef __clang__

    template <typename ScalarType>
    const bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_taped = false;


    template <typename ScalarType>
    const bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_tapeless = false;


    template <typename ScalarType>
    const bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_real_valued =
      (!boost::is_complex<ScalarType>::value);


    template <typename ScalarType>
    const bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_complex_valued =
      !(NumberTraits<ScalarType, NumberTypes::none>::is_real_valued);


    template <typename ScalarType>
    const unsigned int NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::n_supported_derivative_levels =
      0;

#  else

    template <typename ScalarType>
    constexpr bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_taped;


    template <typename ScalarType>
    constexpr bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_tapeless;


    template <typename ScalarType>
    constexpr bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_real_valued;


    template <typename ScalarType>
    constexpr bool NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::is_complex_valued;


    template <typename ScalarType>
    constexpr unsigned int NumberTraits<
      ScalarType,
      NumberTypes::none,
      typename std::enable_if<
        std::is_floating_point<ScalarType>::value ||
        (boost::is_complex<ScalarType>::value &&
         std::is_floating_point<typename internal::RemoveComplexWrapper<
           ScalarType>::type>::value)>::type>::n_supported_derivative_levels;

#  endif


    /**
     * A dummy specialization for floating point numbers.
     *
     * This is necessary to deal with two situations:
     * 1. when the dummy "none" number type is to be used, and
     * 2. the tricky case of higher-order derivative extraction from AD numbers.
     *
     * @note Sacado nests the directional derivatives within Sacado numbers, while
     * ADOL-C does not. So, starting off from a scalar function f(x), the number
     * type resulting from the computation of a directional derivative f'(x) of
     * that function is a floating point number for ADOL-C number types and
     * Sacado::FAD<double>, but that of a Sacado::FAD<Sacado::FAD<double>> is a
     * Sacado::FAD<double>.
     *
     * For this reason, when trying to extract higher-order derivatives of a
     * number that does not support them the input value to this function may be
     * a scalar type.
     */
    template <typename ScalarType>
    struct ADNumberTraits<
      ScalarType,
      typename std::enable_if<std::is_floating_point<ScalarType>::value>::type>
      : NumberTraits<ScalarType, NumberTypes::none>
    {};

    /**
     * A dummy specialization for complex floating point numbers.
     */
    template <typename ComplexScalarType>
    struct ADNumberTraits<
      ComplexScalarType,
      typename std::enable_if<
        boost::is_complex<ComplexScalarType>::value &&
        std::is_floating_point<typename ComplexScalarType::value_type>::value>::
        type> : NumberTraits<ComplexScalarType, NumberTypes::none>
    {};

  } // namespace AD
} // namespace Differentiation

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
