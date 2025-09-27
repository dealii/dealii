// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_adolc_number_types_h
#define dealii_differentiation_ad_adolc_number_types_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/ad_number_types.h>

#ifdef DEAL_II_WITH_ADOLC
#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double
#  include <adolc/internal/adolc_settings.h>
#  include <adolc/internal/adubfunc.h> // Taped double math functions
#endif

#include <complex>
#include <limits>
#include <type_traits>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace AD
  {
    /**
     * A struct to indicate whether a given @p NumberType is an
     * ADOL-C number or not. By default, numbers are not considered to
     * have the necessary characteristics to fulfill this condition.
     */
    template <typename NumberType, typename = void>
    struct is_adolc_number : std::false_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a taped
     * ADOL-C number or not. By default, numbers are not considered to
     * have the necessary characteristics to fulfill this condition.
     */
    template <typename NumberType, typename = void>
    struct is_adolc_taped_number : std::false_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a tapeless
     * ADOL-C number or not. By default, numbers are not considered to
     * have the necessary characteristics to fulfill this condition.
     */
    template <typename NumberType, typename = void>
    struct is_adolc_tapeless_number : std::false_type
    {};
  } // namespace AD
} // namespace Differentiation



#ifdef DEAL_II_WITH_ADOLC

/**
 * An exception which states that a function has been disabled due to the
 * configuration of ADOL-C with the advanced branching feature enabled.
 *
 * @ingroup Exceptions
 */
DeclExceptionMsg(ExcADOLCAdvancedBranching,
                 "This function has not yet been implemented for taped ADOL-C "
                 "numbers when the advanced branching feature is activated.");


/* ----------- inline and template functions and specializations ----------- */


#  ifndef DOXYGEN


namespace Differentiation
{
  namespace AD
  {
    namespace internal
    {
      /**
       * A specialization for the information struct for taped ADOL-C
       * numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::adolc_taped,
        std::enable_if_t<std::is_floating_point_v<ScalarType>>>
      {
        static const bool is_taped = true;
        using real_type            = adouble;
        using derivative_type      = double;
        static const unsigned int n_supported_derivative_levels =
          std::numeric_limits<unsigned int>::max();
      };


      /**
       * A specialization for the information struct for tapeless ADOL-C
       * numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::adolc_tapeless,
        std::enable_if_t<std::is_floating_point_v<ScalarType>>>
      {
        static const bool is_taped                              = false;
        using real_type                                         = adtl::adouble;
        using derivative_type                                   = double;
        static const unsigned int n_supported_derivative_levels = 1;
      };


      template <typename ADNumberType>
      struct Marking<
        ADNumberType,
        std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                           NumberTypes::adolc_taped &&
                         ADNumberTraits<ADNumberType>::is_real_valued>>
      {
        using scalar_type = typename ADNumberTraits<ADNumberType>::scalar_type;

        /*
         * Initialize the state of an independent variable.
         */
        static void
        independent_variable(const scalar_type &in,
                             const unsigned int,
                             const unsigned int,
                             ADNumberType &out)
        {
          out <<= in;
        }

        /*
         * Initialize the state of a dependent variable.
         *
         * @note The second argument must be writable, so we
         * simply pass a copy instead of a non-constant reference.
         */
        static void
        dependent_variable(ADNumberType &out, ADNumberType func)
        {
          // Store the value only (strip it of all sensitivities)
          out = ADNumberTraits<ADNumberType>::get_scalar_value(func);
          // Mark as a dependent variable
          scalar_type tmp;
          func >>= tmp;
        }
      };

      template <typename ADNumberType>
      struct Marking<
        ADNumberType,
        std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                           NumberTypes::adolc_tapeless &&
                         ADNumberTraits<ADNumberType>::is_real_valued>>
      {
        using scalar_type = typename ADNumberTraits<ADNumberType>::scalar_type;

        /*
         * Initialize the state of an independent variable.
         */
        static void
        independent_variable(const scalar_type &in,
                             const unsigned int index,
                             const unsigned int,
                             ADNumberType &out)
        {
          // It is important that the tapeless variables have their values set
          // before defining their directional derivative index
          out = in;

          // Violating this condition when will result in an ADOL-C internal
          // error. We could rather always throw here in order to provide a
          // less cryptic message.
          AssertThrow(index < adtl::getNumDir(),
                      ExcMessage(
                        "The index number of the independent variable being "
                        "marked is greater than the number of independent "
                        "variables that have been declared."));
          out.setADValue(index, 1 /*seed value for first derivative*/);
        }

        /*
         * Initialize the state of a dependent variable.
         */
        static void
        dependent_variable(ADNumberType &out, const ADNumberType &func)
        {
          // Simply transfer value with sensitivities
          out = 0.0;
          out = func;
        }
      };


      /**
       * A struct to help extract certain information associated with
       * taped ADOL-C auto-differentiable numbers.
       */
      template <>
      struct ExtractData<adouble>
      {
        /**
         * Extract the real value.
         */
        static double
        value(const adouble &x)
        {
          return x.getValue();
        }


        /**
         * Extract the number of directional derivatives.
         *
         * This information is not available for taped numbers, and this
         * function exists for aesthetic/compatibility reasons only.
         */
        static unsigned int
        n_directional_derivatives(const adouble &)
        {
          return 0;
        }


        /**
         * Extract the directional derivative in the specified @p direction.
         *
         * It is in fact not possible to perform this operation in this manner,
         * for taped numbers, and this function exists for
         * aesthetic/compatibility reasons only.
         */
        static double
        directional_derivative(const adouble &, const unsigned int)
        {
          AssertThrow(false,
                      ExcMessage(
                        "The derivative values for taped ADOL-C numbers must be"
                        " computed through the ::gradient function."));
          return 0.0;
        }
      };


      /**
       * A struct to help extract certain information associated with
       * tapeless ADOL-C auto-differentiable numbers.
       */
      template <>
      struct ExtractData<adtl::adouble>
      {
        /**
         * Extract the floating point value.
         */
        static double
        value(const adtl::adouble &x)
        {
          return x.getValue();
        }


        /**
         * Extract the number of directional derivatives.
         */
        static unsigned int
        n_directional_derivatives(const adtl::adouble &)
        {
          // This is a global function call...
          return adtl::getNumDir();
        }


        /**
         * Extract the directional derivative in the specified @p direction.
         */
        static double
        directional_derivative(const adtl::adouble &x,
                               const unsigned int   direction)
        {
          Assert(
            direction < n_directional_derivatives(x),
            ExcMessage(
              "Requested directional derivative is greater than the number "
              "registered by ADOL-C."));
          return x.getADValue(direction);
        }
      };

    } // namespace internal



    /**
     * Specialization of the general AdolCWrappers::ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * taped ADOL-C (real) double.
     *
     * @note In this case the number traits are the same as those for a taped double.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<ADNumberType, adouble>>>
      : NumberTraits<double, NumberTypes::adolc_taped>
    {
      static_assert(std::is_same_v<ad_type, adouble>,
                    "Incorrect template type selected for taped ad_type");
      static_assert(is_taped == true, "Incorrect setting for taping");
    };



    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * taped ADOL-C complex double.
     *
     * @note In this case the number traits are the same as those for a taped complex
     * double.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<ADNumberType, std::complex<adouble>>>>
      : NumberTraits<std::complex<double>, NumberTypes::adolc_taped>
    {
      static_assert(std::is_same_v<ad_type, std::complex<adouble>>,
                    "Incorrect template type selected for taped ad_type");
      static_assert(is_taped == true, "Incorrect setting for taping");
    };



    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * tapeless ADOL-C (real) double.
     *
     * @note In this case the number traits are the same as those for a tapeless double.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<ADNumberType, adtl::adouble>>>
      : NumberTraits<double, NumberTypes::adolc_tapeless>
    {
      static_assert(std::is_same_v<ad_type, adtl::adouble>,
                    "Incorrect template type selected for tapeless ad_type");
      static_assert(is_tapeless == true, "Incorrect setting for taping");
    };



    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * tapeless ADOL-C complex double.
     *
     * @note In this case the number traits are the same as those for a tapeless
     * complex double.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<
        std::is_same_v<ADNumberType, std::complex<adtl::adouble>>>>
      : NumberTraits<std::complex<double>, NumberTypes::adolc_tapeless>
    {
      static_assert(std::is_same_v<ad_type, std::complex<adtl::adouble>>,
                    "Incorrect template type selected for tapeless ad_type");
      static_assert(is_tapeless == true, "Incorrect setting for taping");
    };



    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) taped ADOL-C number type.
     */
    template <>
    struct NumberTraits<adouble, NumberTypes::adolc_taped>
      : NumberTraits<typename ADNumberTraits<adouble>::scalar_type,
                     NumberTypes::adolc_taped>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) taped ADOL-C complex number type.
     */
    template <>
    struct NumberTraits<std::complex<adouble>, NumberTypes::adolc_taped>
      : NumberTraits<
          typename ADNumberTraits<std::complex<adouble>>::scalar_type,
          NumberTypes::adolc_taped>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) tapeless ADOL-C number type.
     */
    template <>
    struct NumberTraits<adtl::adouble, NumberTypes::adolc_tapeless>
      : NumberTraits<typename ADNumberTraits<adtl::adouble>::scalar_type,
                     NumberTypes::adolc_tapeless>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) tapeless ADOL-C complex number type.
     */
    template <>
    struct NumberTraits<std::complex<adtl::adouble>,
                        NumberTypes::adolc_tapeless>
      : NumberTraits<
          typename ADNumberTraits<std::complex<adtl::adouble>>::scalar_type,
          NumberTypes::adolc_tapeless>
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (real or complex) taped ADOL-C number.
     */
    template <typename NumberType>
    struct is_adolc_taped_number<
      NumberType,
      std::enable_if_t<ADNumberTraits<std::decay_t<NumberType>>::type_code ==
                       NumberTypes::adolc_taped>> : std::true_type
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (real or complex) tapeless ADOL-C number.
     */
    template <typename NumberType>
    struct is_adolc_tapeless_number<
      NumberType,
      std::enable_if_t<ADNumberTraits<std::decay_t<NumberType>>::type_code ==
                       NumberTypes::adolc_tapeless>> : std::true_type
    {};


    /**
     * Specialization of the struct for the case when the input template
     * parameter is a (real or complex; taped or tapeless) ADOL-C number.
     */
    template <typename NumberType>
    struct is_adolc_number<
      NumberType,
      std::enable_if_t<is_adolc_taped_number<NumberType>::value ||
                       is_adolc_tapeless_number<NumberType>::value>>
      : std::true_type
    {};

  } // namespace AD
} // namespace Differentiation


#  endif // DOXYGEN


#endif // DEAL_II_WITH_ADOLC

DEAL_II_NAMESPACE_CLOSE


#endif
