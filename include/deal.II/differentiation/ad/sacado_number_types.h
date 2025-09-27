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

#ifndef dealii_differentiation_ad_sacado_number_types_h
#define dealii_differentiation_ad_sacado_number_types_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/ad_number_types.h>

#ifdef DEAL_II_TRILINOS_WITH_SACADO
// It appears that some versions of Trilinos do not directly or indirectly
// include all the headers for all forward and reverse Sacado AD types
// in Sacado.hpp, so we also directly include these here as a precaution:
#  include <Sacado.hpp>
#  include <Sacado_Fad_DFad.hpp>
#  include <Sacado_trad.hpp>

#  include <complex>
#endif

#include <type_traits>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace AD
  {
    /**
     * A struct to indicate whether a given @p NumberType is a
     * Sacado number or not. By default, numbers are not considered to
     * have the necessary characteristics to fulfill this condition.
     */
    template <typename NumberType, typename = void>
    struct is_sacado_number : std::false_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a supported Sacado::Fad
     * number or not. By default, numbers are not considered to have the
     * necessary characteristics to fulfill this condition.
     */
    template <typename NumberType, typename = void>
    struct is_sacado_dfad_number : std::false_type
    {};


    /**
     * A struct to indicate whether a given @p NumberType is a supported Sacado::Rad
     * number or not. By default, numbers are not considered to have the
     * necessary characteristics to fulfill this condition.
     */
    template <typename NumberType, typename = void>
    struct is_sacado_rad_number : std::false_type
    {};

  } // namespace AD
} // namespace Differentiation



#ifdef DEAL_II_TRILINOS_WITH_SACADO


namespace Differentiation
{
  namespace AD
  {
    namespace internal
    {
      /**
       * A struct that provides a uniform interface to critical
       * implementation details of a @p SacadoNumber. It defines
       * various number types, and records how many levels of
       * differentiation the number is able to support.
       */
      template <typename SacadoNumber, typename = void>
      struct SacadoNumberInfo;

    } // namespace internal



  } // namespace AD
} // namespace Differentiation


/* ----------- inline and template functions and specializations ----------- */


#  ifndef DOXYGEN

namespace Differentiation
{
  namespace AD
  {
    namespace internal
    {
      // The documentation on Sacado numbers is pretty sparse and/or hard to
      // navigate. As a point of reference, see
      // https://trilinos.org/docs/dev/packages/sacado/doc/html/classSacado_1_1Fad_1_1SimpleFad.html
      // for semi-applicable documentation for the Sacado::Fad::Dfad class.
      // and the examples in
      // https://github.com/trilinos/Trilinos/tree/master/packages/sacado/example
      //
      // If one dares to venture there, the relevant files for the classes
      // supported here are:
      //
      // Forward-mode auto-differentiable types:
      // https://github.com/trilinos/Trilinos/blob/master/packages/sacado/src/sacado_dfad_DFad.hpp
      // https://github.com/trilinos/Trilinos/blob/master/packages/sacado/src/sacado_dfad_GeneralFad.hpp
      //
      // Reverse-mode auto-differentiable types:
      // https://github.com/trilinos/Trilinos/blob/master/packages/sacado/src/Sacado_trad.hpp


      /**
       * Specialization for Sacado::Fad numbers
       */
      template <typename SacadoNumber>
      struct SacadoNumberInfo<
        SacadoNumber,
        std::enable_if_t<
          std::is_same_v<SacadoNumber,
                         Sacado::Fad::DFad<typename SacadoNumber::value_type>>>>
      {
        using ad_type         = SacadoNumber;
        using scalar_type     = typename ad_type::scalar_type;
        using value_type      = typename ad_type::value_type;
        using derivative_type = typename ad_type::value_type;

        static const unsigned int n_supported_derivative_levels =
          1 + SacadoNumberInfo<derivative_type>::n_supported_derivative_levels;
      };


      /**
       * Specialization for Sacado::Rad numbers
       */
      template <typename SacadoNumber>
      struct SacadoNumberInfo<
        SacadoNumber,
        std::enable_if_t<std::is_same_v<
          SacadoNumber,
          Sacado::Rad::ADvar<typename SacadoNumber::value_type>>>>
      {
        using ad_type         = SacadoNumber;
        using scalar_type     = typename ad_type::ADVari::scalar_type;
        using value_type      = typename ad_type::ADVari::value_type;
        using derivative_type = typename ad_type::ADVari::value_type;

        static const unsigned int n_supported_derivative_levels =
          1 + SacadoNumberInfo<derivative_type>::n_supported_derivative_levels;
      };


      /**
       * Specialization for floating point numbers.
       *
       * This is required as a termination point for the recursive
       * templates used in the above specializations.
       */
      template <typename Number>
      struct SacadoNumberInfo<
        Number,
        std::enable_if_t<std::is_arithmetic_v<std::decay_t<Number>>>>
      {
        static const unsigned int n_supported_derivative_levels = 0;
      };


      /**
       * A specialization for the information struct for Sacado dynamic forward
       * auto-differentiable numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::sacado_dfad,
        std::enable_if_t<std::is_floating_point_v<ScalarType>>>
      {
        static const bool is_taped = false;
        using real_type            = Sacado::Fad::DFad<ScalarType>;
        using derivative_type =
          typename SacadoNumberInfo<real_type>::derivative_type;
        static const unsigned int n_supported_derivative_levels =
          SacadoNumberInfo<real_type>::n_supported_derivative_levels;
      };


      /**
       * A specialization for the information struct for nested Sacado dynamic
       * forward auto-differentiable numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::sacado_dfad_dfad,
        std::enable_if_t<std::is_floating_point_v<ScalarType>>>
      {
        static const bool is_taped = false;
        using real_type = Sacado::Fad::DFad<Sacado::Fad::DFad<ScalarType>>;
        using derivative_type =
          typename SacadoNumberInfo<real_type>::derivative_type;
        static const unsigned int n_supported_derivative_levels =
          SacadoNumberInfo<real_type>::n_supported_derivative_levels;
      };


      /**
       * A specialization for the information struct for Sacado dynamic reverse
       * auto-differentiable numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::sacado_rad,
        std::enable_if_t<std::is_floating_point_v<ScalarType>>>
      {
        static const bool is_taped = false;
        using real_type            = Sacado::Rad::ADvar<ScalarType>;
        using derivative_type =
          typename SacadoNumberInfo<real_type>::derivative_type;
        static const unsigned int n_supported_derivative_levels =
          SacadoNumberInfo<real_type>::n_supported_derivative_levels;
      };


      /**
       * A specialization for the information struct for Sacado dynamic nested
       * reverse-forward auto-differentiable numbers.
       */
      template <typename ScalarType>
      struct ADNumberInfoFromEnum<
        ScalarType,
        Differentiation::AD::NumberTypes::sacado_rad_dfad,
        std::enable_if_t<std::is_floating_point_v<ScalarType>>>
      {
        static const bool is_taped = false;
        using real_type = Sacado::Rad::ADvar<Sacado::Fad::DFad<ScalarType>>;
        using derivative_type =
          typename SacadoNumberInfo<real_type>::derivative_type;
        static const unsigned int n_supported_derivative_levels =
          SacadoNumberInfo<real_type>::n_supported_derivative_levels;
      };


      /**
       * Specialization of the marking strategy for Sacado::Fad::DFad
       * auto-differentiable numbers
       */
      template <typename NumberType>
      struct Marking<Sacado::Fad::DFad<NumberType>>
      {
        using ad_type =
          typename SacadoNumberInfo<Sacado::Fad::DFad<NumberType>>::ad_type;
        using derivative_type = typename SacadoNumberInfo<
          Sacado::Fad::DFad<NumberType>>::derivative_type;
        using scalar_type =
          typename SacadoNumberInfo<Sacado::Fad::DFad<NumberType>>::scalar_type;

        /*
         * Initialize the state of an independent variable.
         */
        static void
        independent_variable(const scalar_type &in,
                             const unsigned int index,
                             const unsigned int n_independent_variables,
                             ad_type           &out)
        {
          // It is required that we first initialize the outer number before
          // any of the nested ones.
          out = ad_type(n_independent_variables, index, in);

          // Initialize potential nested directional derivatives
          Marking<derivative_type>::independent_variable(
            in, index, n_independent_variables, out.val());
        }

        /*
         * Initialize the state of a dependent variable.
         */
        static void
        dependent_variable(ad_type &out, const ad_type &func)
        {
          out = func;
        }
      };


      /**
       * Specialization of the marking strategy for Sacado::Rad::ADvar
       * auto-differentiable numbers.
       */
      template <typename NumberType>
      struct Marking<Sacado::Rad::ADvar<NumberType>>
      {
        using ad_type =
          typename SacadoNumberInfo<Sacado::Rad::ADvar<NumberType>>::ad_type;
        using derivative_type = typename SacadoNumberInfo<
          Sacado::Rad::ADvar<NumberType>>::derivative_type;
        using scalar_type = typename SacadoNumberInfo<
          Sacado::Rad::ADvar<NumberType>>::scalar_type;

        /*
         * Initialize the state of an independent variable.
         */
        static void
        independent_variable(const scalar_type &in,
                             const unsigned int index,
                             const unsigned int n_independent_variables,
                             ad_type           &out)
        {
          // For Sacado::Rad::ADvar numbers, we have to initialize the
          // ADNumber with an already fully-configured value. This means
          // that if this nests another ADNumber then the nested number
          // must already be set up and ready for use.

          // Initialize potential nested directional derivatives
          derivative_type derivative_initializer;
          Marking<derivative_type>::independent_variable(
            in, index, n_independent_variables, derivative_initializer);

          // Initialize the outer ad_type
          out = derivative_initializer;
        }

        /*
         * Initialize the state of a dependent variable.
         */
        static void
        dependent_variable(ad_type &out, const ad_type &func)
        {
          out = func;
        }
      };


      /**
       * A struct to help extract certain information associated with
       * Sacado dynamic reverse auto-differentiable numbers. The @p NumberType
       * can be either a floating point number or another Sacado type.
       */
      template <typename NumberType>
      struct ExtractData<Sacado::Fad::DFad<NumberType>>
      {
        using derivative_type = typename SacadoNumberInfo<
          Sacado::Fad::DFad<NumberType>>::derivative_type;
        using scalar_type =
          typename SacadoNumberInfo<Sacado::Fad::DFad<NumberType>>::scalar_type;
        using value_type =
          typename SacadoNumberInfo<Sacado::Fad::DFad<NumberType>>::value_type;

        /**
         * Extract the real scalar value.
         */
        static scalar_type
        value(const Sacado::Fad::DFad<NumberType> &x)
        {
          return ExtractData<value_type>::value(x.val());
        }


        /**
         * Extract the number of directional derivatives.
         */
        static unsigned int
        n_directional_derivatives(const Sacado::Fad::DFad<NumberType> &x)
        {
          return x.size();
        }


        /**
         * Extract the directional derivative in the specified @p direction.
         */
        static derivative_type
        directional_derivative(const Sacado::Fad::DFad<NumberType> &x,
                               const unsigned int                   direction)
        {
          if (x.hasFastAccess())
            return x.fastAccessDx(direction);
          else
            return x.dx(direction);
        }
      };


      /**
       * A struct to help extract certain information associated with
       * Sacado dynamic reverse auto-differentiable numbers. The @p NumberType
       * can be either a floating point number or another Sacado type.
       */
      template <typename NumberType>
      struct ExtractData<Sacado::Rad::ADvar<NumberType>>
      {
        using derivative_type = typename SacadoNumberInfo<
          Sacado::Rad::ADvar<NumberType>>::derivative_type;
        using scalar_type = typename SacadoNumberInfo<
          Sacado::Rad::ADvar<NumberType>>::scalar_type;
        using value_type =
          typename SacadoNumberInfo<Sacado::Rad::ADvar<NumberType>>::value_type;

        /**
         * Extract the real scalar value.
         */
        static scalar_type
        value(const Sacado::Rad::ADvar<NumberType> &x)
        {
          return ExtractData<value_type>::value(x.val());
        }


        /**
         * Extract the number of directional derivatives.
         */
        static unsigned int
        n_directional_derivatives(const Sacado::Rad::ADvar<NumberType> &)
        {
          // There are as many directional derivatives as there are
          // independent variables, but each independent variable can
          // only return one directional derivative.
          return 1;
        }


        /**
         * Extract the directional derivative in the specified @p direction.
         *
         * @note For reverse-mode AD, @p x should represent an independent
         * variable with respect to which one wishes to take the derivative
         * @p df/dx of a dependent function @p f(x).
         */
        static derivative_type
        directional_derivative(const Sacado::Rad::ADvar<NumberType> &x,
                               const unsigned int)
        {
          return x.adj();
        }
      };

    } // namespace internal


    /* -------------- NumberTypes::sacado_dfad -------------- */


    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * Sacado::Fad::DFad number templated on a floating point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<
        std::is_same_v<ADNumberType,
                       Sacado::Fad::DFad<typename ADNumberType::scalar_type>>>>
      : NumberTraits<typename ADNumberType::scalar_type,
                     NumberTypes::sacado_dfad>
    {};


    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * complex Sacado::Fad::DFad number templated on a floating point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<
        std::is_same_v<ADNumberType,
                       std::complex<Sacado::Fad::DFad<
                         typename ADNumberType::value_type::scalar_type>>>>>
      : NumberTraits<
          std::complex<typename ADNumberType::value_type::scalar_type>,
          NumberTypes::sacado_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Fad::DFad<float>, NumberTypes::sacado_dfad>
      : NumberTraits<
          typename ADNumberTraits<Sacado::Fad::DFad<float>>::scalar_type,
          NumberTypes::sacado_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) complex Sacado number type.
     */
    template <>
    struct NumberTraits<std::complex<Sacado::Fad::DFad<float>>,
                        NumberTypes::sacado_dfad>
      : NumberTraits<typename ADNumberTraits<
                       std::complex<Sacado::Fad::DFad<float>>>::scalar_type,
                     NumberTypes::sacado_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Fad::DFad<double>, NumberTypes::sacado_dfad>
      : NumberTraits<
          typename ADNumberTraits<Sacado::Fad::DFad<double>>::scalar_type,
          NumberTypes::sacado_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) complex Sacado number type.
     */
    template <>
    struct NumberTraits<std::complex<Sacado::Fad::DFad<double>>,
                        NumberTypes::sacado_dfad>
      : NumberTraits<typename ADNumberTraits<
                       std::complex<Sacado::Fad::DFad<double>>>::scalar_type,
                     NumberTypes::sacado_dfad>
    {};


    /* -------------- NumberTypes::sacado_rad -------------- */


    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * Sacado::Rad::ADvar number templated on a floating point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<
        ADNumberType,
        Sacado::Rad::ADvar<typename ADNumberType::ADVari::scalar_type>>>>
      : NumberTraits<typename ADNumberType::ADVari::scalar_type,
                     NumberTypes::sacado_rad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Rad::ADvar<float>, NumberTypes::sacado_rad>
      : NumberTraits<
          typename ADNumberTraits<Sacado::Rad::ADvar<float>>::scalar_type,
          NumberTypes::sacado_rad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Rad::ADvar<double>, NumberTypes::sacado_rad>
      : NumberTraits<
          typename ADNumberTraits<Sacado::Rad::ADvar<double>>::scalar_type,
          NumberTypes::sacado_rad>
    {};


#    ifdef DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD


    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * complex Sacado::Rad::ADvar number templated on a float type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<
        ADNumberType,
        std::complex<Sacado::Rad::ADvar<
          typename ADNumberType::value_type::ADVari::scalar_type>>>>>
      : NumberTraits<
          std::complex<typename ADNumberType::value_type::ADVari::scalar_type>,
          NumberTypes::sacado_rad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) complex Sacado number type.
     */
    template <>
    struct NumberTraits<std::complex<Sacado::Rad::ADvar<float>>,
                        NumberTypes::sacado_rad>
      : NumberTraits<typename ADNumberTraits<
                       std::complex<Sacado::Rad::ADvar<float>>>::scalar_type,
                     NumberTypes::sacado_rad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) complex Sacado number type.
     */
    template <>
    struct NumberTraits<std::complex<Sacado::Rad::ADvar<double>>,
                        NumberTypes::sacado_rad>
      : NumberTraits<typename ADNumberTraits<
                       std::complex<Sacado::Rad::ADvar<double>>>::scalar_type,
                     NumberTypes::sacado_rad>
    {};


#    endif


    /* -------------- NumberTypes::sacado_dfad_dfad -------------- */

    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * once nested Sacado::Fad::DFad number, with the inner number templated on
     * a floating point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<ADNumberType,
                                      Sacado::Fad::DFad<Sacado::Fad::DFad<
                                        typename ADNumberType::scalar_type>>>>>
      : NumberTraits<typename ADNumberType::scalar_type,
                     NumberTypes::sacado_dfad_dfad>
    {};


    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * complex nested Sacado::Fad::DFad number templated on a floating
     * point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<
        std::is_same_v<ADNumberType,
                       std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<
                         typename ADNumberType::value_type::scalar_type>>>>>>
      : NumberTraits<
          std::complex<typename ADNumberType::value_type::scalar_type>,
          NumberTypes::sacado_dfad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>,
                        NumberTypes::sacado_dfad_dfad>
      : NumberTraits<typename ADNumberTraits<Sacado::Fad::DFad<
                       Sacado::Fad::DFad<float>>>::scalar_type,
                     NumberTypes::sacado_dfad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested complex Sacado number type.
     */
    template <>
    struct NumberTraits<
      std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>>,
      NumberTypes::sacado_dfad_dfad>
      : NumberTraits<typename ADNumberTraits<std::complex<Sacado::Fad::DFad<
                       Sacado::Fad::DFad<float>>>>::scalar_type,
                     NumberTypes::sacado_dfad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>,
                        NumberTypes::sacado_dfad_dfad>
      : NumberTraits<typename ADNumberTraits<Sacado::Fad::DFad<
                       Sacado::Fad::DFad<double>>>::scalar_type,
                     NumberTypes::sacado_dfad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested complex Sacado number type.
     */
    template <>
    struct NumberTraits<
      std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>,
      NumberTypes::sacado_dfad_dfad>
      : NumberTraits<typename ADNumberTraits<std::complex<Sacado::Fad::DFad<
                       Sacado::Fad::DFad<double>>>>::scalar_type,
                     NumberTypes::sacado_dfad_dfad>
    {};


    /* -------------- NumberTypes::sacado_rad_dfad -------------- */

    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * once nested Sacado::Fad::DFad number, with the inner number templated on
     * a floating point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<ADNumberType,
                          std::enable_if_t<std::is_same_v<
                            ADNumberType,
                            Sacado::Rad::ADvar<Sacado::Fad::DFad<
                              typename ADNumberType::ADVari::scalar_type>>>>>
      : NumberTraits<typename ADNumberType::ADVari::scalar_type,
                     NumberTypes::sacado_rad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>,
                        NumberTypes::sacado_rad_dfad>
      : NumberTraits<typename ADNumberTraits<Sacado::Rad::ADvar<
                       Sacado::Fad::DFad<float>>>::scalar_type,
                     NumberTypes::sacado_rad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested Sacado number type.
     */
    template <>
    struct NumberTraits<Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>,
                        NumberTypes::sacado_rad_dfad>
      : NumberTraits<typename ADNumberTraits<Sacado::Rad::ADvar<
                       Sacado::Fad::DFad<double>>>::scalar_type,
                     NumberTypes::sacado_rad_dfad>
    {};


#    ifdef DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD


    /**
     * Specialization of the general ADNumberTraits class that
     * provides relevant information for auto-differentiable numbers.
     * This specialization is for the case where @p ADNumberType is an
     * complex nested Sacado::Fad::DFad number templated on a floating
     * point type.
     */
    template <typename ADNumberType>
    struct ADNumberTraits<
      ADNumberType,
      std::enable_if_t<std::is_same_v<
        ADNumberType,
        std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<
          typename ADNumberType::value_type::ADVari::scalar_type>>>>>>
      : NumberTraits<
          std::complex<typename ADNumberType::value_type::ADVari::scalar_type>,
          NumberTypes::sacado_rad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested complex Sacado number type.
     */
    template <>
    struct NumberTraits<
      std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>>,
      NumberTypes::sacado_rad_dfad>
      : NumberTraits<typename ADNumberTraits<std::complex<Sacado::Rad::ADvar<
                       Sacado::Fad::DFad<float>>>>::scalar_type,
                     NumberTypes::sacado_rad_dfad>
    {};


    /**
     * Specialization of the NumberTraits struct for
     * the (otherwise disabled) nested complex Sacado number type.
     */
    template <>
    struct NumberTraits<
      std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>>,
      NumberTypes::sacado_rad_dfad>
      : NumberTraits<typename ADNumberTraits<std::complex<Sacado::Rad::ADvar<
                       Sacado::Fad::DFad<double>>>>::scalar_type,
                     NumberTypes::sacado_rad_dfad>
    {};


#    endif


    /* -------------- Additional type traits -------------- */


    template <typename NumberType>
    struct is_sacado_dfad_number<
      NumberType,
      std::enable_if_t<ADNumberTraits<std::decay_t<NumberType>>::type_code ==
                         NumberTypes::sacado_dfad ||
                       ADNumberTraits<std::decay_t<NumberType>>::type_code ==
                         NumberTypes::sacado_dfad_dfad>> : std::true_type
    {};


    template <typename NumberType>
    struct is_sacado_dfad_number<
      NumberType,
      std::enable_if_t<
        std::is_same_v<NumberType,
                       Sacado::Fad::Expr<typename NumberType::value_type>>>>
      : std::true_type
    {};


    template <typename NumberType>
    struct is_sacado_rad_number<
      NumberType,
      std::enable_if_t<ADNumberTraits<std::decay_t<NumberType>>::type_code ==
                         NumberTypes::sacado_rad ||
                       ADNumberTraits<std::decay_t<NumberType>>::type_code ==
                         NumberTypes::sacado_rad_dfad>> : std::true_type
    {};


    template <typename NumberType>
    struct is_sacado_rad_number<
      NumberType,
      std::enable_if_t<
        std::is_same_v<NumberType,
                       Sacado::Rad::ADvari<Sacado::Fad::DFad<
                         typename NumberType::ADVari::scalar_type>>>>>
      : std::true_type
    {};


    template <typename NumberType>
    struct is_sacado_number<
      NumberType,
      std::enable_if_t<is_sacado_dfad_number<NumberType>::value ||
                       is_sacado_rad_number<NumberType>::value>>
      : std::true_type
    {};

  } // namespace AD
} // namespace Differentiation


#  endif // DOXYGEN


#endif // DEAL_II_TRILINOS_WITH_SACADO


DEAL_II_NAMESPACE_CLOSE

#endif
