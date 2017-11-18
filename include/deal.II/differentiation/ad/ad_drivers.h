// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_differentiation_ad_ad_drivers_h
#define dealii_differentiation_ad_ad_drivers_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_ADOLC

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#endif // DEAL_II_WITH_ADOLC

#include <vector>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace AD
  {

    /**
     * @addtogroup Exceptions
     */
    //@{

    /**
     * Exception denoting that a class requires some specialization
     * in order to be used.
     */
    DeclExceptionMsg (ExcRequiresADNumberSpecialization,
                      "This function is called in a class that is expected to be specialized "
                      "for auto-differentiable numbers.");

    /**
     * Exception denoting that Adol-C is a required feature.
     */
    DeclExceptionMsg (ExcRequiresAdolC,
                      "This function is only available if deal.II is compiled with ADOL-C.");

    /**
     * This exception is raised whenever the an auto-differentiable number does not
     * support the required number of derivative operations
     *
     * The first parameter to the constructor is the number of derivative operations
     * that it provides, and the second is the minimum number that are required.
     * Both parameters are of type <tt>int</tt>.
     */
    DeclException2 (ExcSupportedDerivativeLevels,
                    std::size_t, std::size_t,
                    << "The number of derivative levels that this auto-differentiable number supports is "
                    << arg1 << ", but it is required that it supports at least " << arg2 << " levels.");

    //@}

    /**
     * A driver class for taped auto-differentiable numbers.
     *
     * It is intended that this class be specialized for the valid
     * combinations of auto-differentiable numbers and output scalar
     * number types.
     *
     * @tparam ADNumberType A type corresponding to a supported
     *         auto-differentiable number.
     * @tparam ScalarType A real or complex floating point number.
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template<typename ADNumberType, typename ScalarType, typename T = void>
    struct TapedDrivers
    {
      // This dummy class definition safely supports compilation
      // against tapeless numbers or taped number types that have
      // not yet been implemented.

      /**
       * @name Taping
       */
      //@{

      /**
       * Enable the recording mode for a given tape.
       *
       * @param[in] tape_index The index of the tape to be written
       * @param[in] keep_independent_values Determines whether the numerical
       *            values of all independent variables are recorded in the
       *            tape buffer.
       */
      static void
      enable_taping(const unsigned int &tape_index,
                    const bool         &keep_independent_values);

      /**
       * Disable the recording mode for a given tape.
       *
       * @param[in] active_tape_index The index of the (currently active) tape
       *            to be finalized and potentially written to file.
       * @param[in] write_tapes_to_file A flag that specified whether the tape
       *            should be written to file or kept in memory.
       */
      static void
      disable_taping(const unsigned int &active_tape_index,
                     const bool         &write_tapes_to_file);

      /**
       * Prints the statistics regarding the usage of the tapes.
       *
       * @param[in] stream The output stream to which the values are to be written.
       * @param[in] tape_index The index of the tape to get the statistics of.
       */
      template<typename Stream>
      static void
      print_tape_stats(Stream             &stream,
                       const unsigned int &tape_index);

      //@}

      /**
       * @name Drivers for scalar functions (one dependent variable)
       */
      //@{

      /**
       * Computes the value of the scalar field.
       *
       * @param[in] active_tape_index The index of the tape on which the dependent
       *            function is recorded.
       * @param[in] n_independent_variables The number of independent variables
       *            whose sensitivities were tracked.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       *
       * @return The scalar values of the function.
       */
      static ScalarType
      value (const unsigned int            &active_tape_index,
             const unsigned int            &n_independent_variables,
             const std::vector<ScalarType> &independent_variables);

      /**
       * Computes the gradient of the scalar field with respect to all independent
       * variables.
       *
       * @param[in] active_tape_index The index of the tape on which the dependent
       *            function is recorded.
       * @param[in] n_independent_variables The number of independent variables
       *            whose sensitivities were tracked.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] gradient The scalar values of the dependent function's gradients.
       */
      static void
      gradient (const unsigned int            &active_tape_index,
                const unsigned int            &n_independent_variables,
                const std::vector<ScalarType> &independent_variables,
                Vector<ScalarType>            &gradient);

      /**
       * Computes the hessian of the scalar field with respect to all independent
       * variables.
       *
       * @param[in] active_tape_index The index of the tape on which the dependent
       *            function is recorded.
       * @param[in] n_independent_variables The number of independent variables
       *            whose sensitivities were tracked.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] hessian The scalar values of the dependent function's hessian.
       */
      static void
      hessian (const unsigned int            &active_tape_index,
               const unsigned int            &n_independent_variables,
               const std::vector<ScalarType> &independent_variables,
               FullMatrix<ScalarType>        &hessian);

      //@}

      /**
       * @name Drivers for vector functions (multiple dependent variables)
       */
      //@{

      /**
       * Computes the values of the vector field.
       *
       * @param[in] active_tape_index The index of the tape on which the dependent
       *            function is recorded.
       * @param[in] n_dependent_variables The number of dependent variables.
       * @param[in] n_independent_variables The number of independent variables
       *            whose sensitivities were tracked.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] values The scalar values of the dependent functions.
       */
      static void
      values (const unsigned int            &active_tape_index,
              const unsigned int            &n_dependent_variables,
              const unsigned int            &n_independent_variables,
              const std::vector<ScalarType> &independent_variables,
              Vector<ScalarType>            &values);

      /**
       * Computes the gradient of the vector field.
       *
       * @param[in] active_tape_index The index of the tape on which the dependent
       *            function is recorded.
       * @param[in] n_dependent_variables The number of dependent variables.
       * @param[in] n_independent_variables The number of independent variables
       *            whose sensitivities were tracked.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] jacobian The scalar values of the dependent function's jacobian.
       */
      static void
      jacobian (const unsigned int            &active_tape_index,
                const unsigned int            &n_dependent_variables,
                const unsigned int            &n_independent_variables,
                const std::vector<ScalarType> &independent_variables,
                FullMatrix<ScalarType>        &jacobian);

      //@}

    };



    /**
     * A prototype driver class for tapeless auto-differentiable numbers.
     *
     * It is intended that this class be specialized for the valid
     * combinations of auto-differentiable numbers and output scalar
     * number types.
     *
     * @tparam ADNumberType A type corresponding to a supported
     *         auto-differentiable number.
     * @tparam ScalarType A real or complex floating point number.
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template<typename ADNumberType, typename ScalarType, typename T = void>
    struct TapelessDrivers
    {
      // This dummy class definition safely supports compilation
      // against taped numbers or tapeless number types that have
      // not yet been implemented.

      /**
       * @name Configuration
       */
      //@{

      /**
       * In the event that the tapeless mode requires <i>a priori</i> knowledge
       * of how many directional derivatives might need to be computed, this function
       * informs the auto-differentiable library of what this number is.
       *
       * @param[in] n_independent_variables The number of independent variables
       *            that will be used in the entire duration of the
       *            simulation.
       *
       * @warning For Adol-C tapeless numbers, the value given to @p n_independent_variables
       *          should be the <b>maximum</b> number of independent variables that will be
       *          used in the entire duration of the simulation. This is important in the
       *          context of, for example, hp-FEM and for multiple constitutive models with
       *          a different number of fields from which a linearization must be computed.
       */
      static void
      initialize (const unsigned int &n_independent_variables);

      //@}

      /**
       * @name Drivers for scalar functions
       */
      //@{

      /**
       * Computes the value of the scalar field.
       *
       * @param[in] dependent_variables The dependent variables whose values are to
       *            be extracted.
       *
       * @return The scalar values of the function.
       */
      static ScalarType
      value (const std::vector<ADNumberType> &dependent_variables);

      /**
       * Computes the gradient of the scalar field with respect to all independent
       * variables.
       *
       * @param[in] independent_variables The independent variables whose sensitivities
       *            were tracked.
       * @param[in] dependent_variables The (single) dependent variable whose gradients
       *            are to be extracted.
       * @param[out] gradient The scalar values of the dependent function's gradients.
       */
      static void
      gradient (const std::vector<ADNumberType> &independent_variables,
                const std::vector<ADNumberType> &dependent_variables,
                Vector<ScalarType>              &gradient);

      /**
       * Computes the hessian of the scalar field with respect to all independent
       * variables.
       *
       * @param[in] independent_variables The independent variables whose sensitivities
       *            were tracked.
       * @param[in] dependent_variables The (single) dependent variable whose hessians
       *            are to be extracted.
       * @param[out] hessian The scalar values of the function's hessian.
       */
      static void
      hessian (const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType>          &hessian);

      //@}

      /**
       * @name Drivers for vector functions
       */
      //@{

      /**
       * Computes the values of the vector field.
       *
       * @param[in] dependent_variables The dependent variables  whose hessians
       *            are to be extracted.
       * @param[out] values The scalar values of the dependent functions.
       */
      static void
      values (const std::vector<ADNumberType> &dependent_variables,
              Vector<ScalarType>              &values);

      /**
       * Computes the gradient of the vector field.
       *
       * @param[in] independent_variables The independent variables whose sensitivities
       *            were tracked.
       * @param[in] dependent_variables The dependent variables whose jacobian
       *            are to be extracted.
       * @param[out] jacobian The scalar values of the function's jacobian.
       */
      static void
      jacobian (const std::vector<ADNumberType> &independent_variables,
                const std::vector<ADNumberType> &dependent_variables,
                FullMatrix<ScalarType>          &jacobian);

      //@}

    };

  }
}



/* --------------------------- inline and template functions ------------------------- */


#ifndef DOXYGEN

namespace Differentiation
{
  namespace AD
  {

    // -------------   TapedDrivers   -------------


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::enable_taping(
      const unsigned int &,
      const bool &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::disable_taping(
      const unsigned int &,
      const bool &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    template<typename Stream>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::print_tape_stats(
      Stream             &stream,
      const unsigned int &tape_index)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    ScalarType
    TapedDrivers<ADNumberType,ScalarType,T>::value (
      const unsigned int &,
      const unsigned int &,
      const std::vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return ScalarType(0.0);
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::gradient (
      const unsigned int &,
      const unsigned int &,
      const std::vector<ScalarType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::hessian (
      const unsigned int &,
      const unsigned int &,
      const std::vector<ScalarType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::values (
      const unsigned int &,
      const unsigned int &,
      const unsigned int &,
      const std::vector<ScalarType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType,ScalarType,T>::jacobian (
      const unsigned int &,
      const unsigned int &,
      const unsigned int &,
      const std::vector<ScalarType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    // Specialization for taped Adol-C auto-differentiable numbers.
    //
    // Note: In the case of Adol-C taped numbers, the associated scalar
    // type is always expected to be a double. So we need to make a further
    // specialization when ScalarType is a float.
    template<typename ADNumberType>
    struct TapedDrivers<ADNumberType,double,typename std::enable_if<
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::adolc_taped
      >::type>
    {
      typedef double scalar_type;

      // === Taping ===

      static void
      enable_taping(const unsigned int &tape_index,
                    const bool         &keep)
      {
        trace_on(tape_index,keep);
      }

      static void
      disable_taping(const unsigned int &active_tape_index,
                     const bool         &write_tapes_to_file)
      {
#ifdef DEAL_II_WITH_ADOLC
        if (write_tapes_to_file)
          {
            trace_off(active_tape_index); // Slow
            std::vector<std::size_t> counts (STAT_SIZE);
            ::tapestats(active_tape_index, counts.data());
          }
        else
          trace_off(); // Fast(er)
#else
        AssertThrow(false, ExcRequiresAdolC());
#endif
      }

      template<typename Stream>
      static void
      print_tape_stats(Stream             &stream,
                       const unsigned int &tape_index)
      {
        // See Adol-C manual section 2.1
        // and adolc/taping.h
        std::vector<std::size_t> counts (STAT_SIZE);
        ::tapestats(tape_index, counts.data());
        Assert(counts.size() >= 18, ExcInternalError());
        stream
            << "Tape index: " << tape_index << "\n"
            << "Number of independent variables: " << counts[0] << "\n"
            << "Number of dependent variables:   " << counts[1] << "\n"
            << "Max number of live, active variables: " << counts[2] << "\n"
            << "Size of taylor stack (number of overwrites): " << counts[3] << "\n"
            << "Operations buffer size: " << counts[4] << "\n"
            << "Total number of recorded operations: " << counts[5] << "\n"
            << "Operations file written or not: " << counts[6] << "\n"
            << "Overall number of locations: " << counts[7] << "\n"
            << "Locations file written or not: " << counts[8] << "\n"
            << "Overall number of values: " << counts[9] << "\n"
            << "Values file written or not: " << counts[10] << "\n"
            << "Locations buffer size: " << counts[11] << "\n"
            << "Values buffer size: " << counts[12] << "\n"
            << "Taylor buffer size: " << counts[13] << "\n"
            << "Number of eq_*_prod for sparsity pattern: " << counts[14] << "\n"
            << "Use of 'min_op', deferred to 'abs_op' for piecewise calculations: " << counts[15] << "\n"
            << "Number of 'abs' calls that can switch branch: " << counts[16] << "\n"
            << "Number of parameters (doubles) interchangable without retaping: " << counts[17] << "\n"
            << std::flush;
      }

      // === Scalar drivers ===

      static scalar_type
      value (const unsigned int             &active_tape_index,
             const unsigned int             &n_independent_variables,
             const std::vector<scalar_type> &independent_variables)
      {

        scalar_type *f = new scalar_type();

#ifdef DEAL_II_WITH_ADOLC
        ::function(active_tape_index,
                   1, // Only one dependent variable
                   n_independent_variables,
                   const_cast<scalar_type *>(independent_variables.data()),
                   f);
#else
        AssertThrow(false, ExcRequiresAdolC());
#endif

        const scalar_type value = f[0];

        // Cleanup :-/
        delete f;
        f = nullptr;

        return value;
      }

      static void
      gradient (const unsigned int             &active_tape_index,
                const unsigned int             &n_independent_variables,
                const std::vector<scalar_type> &independent_variables,
                Vector<scalar_type>            &gradient)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,1));
        Assert(gradient.size() == independent_variables.size(),
               ExcDimensionMismatch(gradient.size(),independent_variables.size()));

        scalar_type *g = new scalar_type[n_independent_variables];

#ifdef DEAL_II_WITH_ADOLC
        ::gradient(active_tape_index,
                   n_independent_variables,
                   const_cast<scalar_type *>(independent_variables.data()),
                   g);
#else
        AssertThrow(false, ExcRequiresAdolC());
#endif

        for (unsigned int i=0; i<n_independent_variables; ++i)
          gradient[i] = g[i];

        // Cleanup :-/
        delete[] g;
        g = nullptr;
      }

      static void
      hessian (const unsigned int             &active_tape_index,
               const unsigned int             &n_independent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type>        &hessian)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,2));
        Assert(hessian.m() == independent_variables.size(),
               ExcDimensionMismatch(hessian.m(),independent_variables.size()));
        Assert(hessian.n() == independent_variables.size(),
               ExcDimensionMismatch(hessian.n(),independent_variables.size()));

        scalar_type **H = new scalar_type*[n_independent_variables];
        for (unsigned int i=0; i<n_independent_variables; ++i)
          H[i] = new scalar_type[i+1]; // Symmetry

#ifdef DEAL_II_WITH_ADOLC
        ::hessian(active_tape_index,
                  n_independent_variables,
                  const_cast<scalar_type *>(independent_variables.data()),
                  H);
#else
        AssertThrow(false, ExcRequiresAdolC());
#endif

        for (unsigned int i=0; i<n_independent_variables; i++)
          for (unsigned int j=0; j<i+1; j++)
            {
              hessian[i][j] = H[i][j];
              if (i != j)
                hessian[j][i] = H[i][j]; // Symmetry
            }

        // Cleanup :-/
        for (unsigned int i=0; i<n_independent_variables; i++)
          delete[] H[i];
        delete[] H;
        H = nullptr;
      }

      // === Vector drivers ===

      static void
      values (const unsigned int             &active_tape_index,
              const unsigned int             &n_dependent_variables,
              const unsigned int             &n_independent_variables,
              const std::vector<scalar_type> &independent_variables,
              Vector<scalar_type>            &values)
      {
        Assert(values.size() == n_dependent_variables,
               ExcDimensionMismatch(values.size(),n_dependent_variables));

        scalar_type *f = new scalar_type[n_dependent_variables];

#ifdef DEAL_II_WITH_ADOLC
        ::function(active_tape_index,
                   n_dependent_variables,
                   n_independent_variables,
                   const_cast<scalar_type *>(independent_variables.data()),
                   f);
#else
        AssertThrow(false, ExcRequiresAdolC());
#endif

        for (unsigned int i=0; i<n_dependent_variables; i++)
          values[i] = f[i];

        // Cleanup :-/
        delete[] f;
        f = nullptr;
      }

      static void
      jacobian (const unsigned int             &active_tape_index,
                const unsigned int             &n_dependent_variables,
                const unsigned int             &n_independent_variables,
                const std::vector<scalar_type> &independent_variables,
                FullMatrix<scalar_type>        &jacobian)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,1));
        Assert(jacobian.m() == n_dependent_variables,
               ExcDimensionMismatch(jacobian.m(),n_dependent_variables));
        Assert(jacobian.n() == independent_variables.size(),
               ExcDimensionMismatch(jacobian.n(),independent_variables.size()));

        scalar_type **J = new scalar_type*[n_dependent_variables];
        for (unsigned int i=0; i<n_dependent_variables; ++i)
          J[i] = new scalar_type[n_independent_variables];

#ifdef DEAL_II_WITH_ADOLC
        ::jacobian(active_tape_index,
                   n_dependent_variables,
                   n_independent_variables,
                   independent_variables.data(),
                   J);
#else
        AssertThrow(false, ExcRequiresAdolC());
#endif

        for (unsigned int i=0; i<n_dependent_variables; i++)
          for (unsigned int j=0; j<n_independent_variables; j++)
            jacobian[i][j] = J[i][j];

        // Cleanup :-/
        for (unsigned int i=0; i<n_dependent_variables; i++)
          delete[] J[i];
        delete[] J;
        J = nullptr;
      }
    };


    // Specialization for Adol-C taped numbers. It is expected that the
    // scalar return type for this class is a float.
    //
    // Note: Adol-C only has drivers for doubles, and so floats are
    // not intrinsically supported. This wrapper struct works around
    // the issue when necessary.
    template<typename ADNumberType>
    struct TapedDrivers<ADNumberType,float,typename std::enable_if<
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::adolc_taped
      >::type>
    {
      typedef float scalar_type;

      static std::vector<double>
      vector_float_to_double (const std::vector<float> &in)
      {
        std::vector<double> out (in.size());
        std::copy(in.begin(), in.end(), out.begin());
        return out;
      }

      // === Taping ===

      static void
      enable_taping(const unsigned int &tape_index,
                    const bool         &keep)
      {
        TapedDrivers<ADNumberType,double>::enable_taping(tape_index,keep);
      }

      static void
      disable_taping(const unsigned int &active_tape_index,
                     const bool         &write_tapes_to_file)
      {
        TapedDrivers<ADNumberType,double>::disable_taping(active_tape_index,write_tapes_to_file);
      }

      template<typename Stream>
      static void
      print_tape_stats(Stream             &stream,
                       const unsigned int &tape_index)
      {
        TapedDrivers<ADNumberType,double>::print_tape_stats(stream,tape_index);
      }

      // === Scalar drivers ===

      static scalar_type
      value (const unsigned int             &active_tape_index,
             const unsigned int             &n_independent_variables,
             const std::vector<scalar_type> &independent_variables)
      {

        return TapedDrivers<ADNumberType,double>::value(
                 active_tape_index,
                 n_independent_variables,
                 vector_float_to_double(independent_variables));
      }

      static void
      gradient (const unsigned int             &active_tape_index,
                const unsigned int             &n_independent_variables,
                const std::vector<scalar_type> &independent_variables,
                Vector<scalar_type>            &gradient)
      {
        Vector<double> gradient_double (gradient.size());
        TapedDrivers<ADNumberType,double>::gradient(
          active_tape_index,
          n_independent_variables,
          vector_float_to_double(independent_variables),
          gradient_double);
        gradient = gradient_double;
      }

      static void
      hessian (const unsigned int             &active_tape_index,
               const unsigned int             &n_independent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type>        &hessian)
      {
        FullMatrix<double> hessian_double (hessian.m(), hessian.n());
        TapedDrivers<ADNumberType,double>::hessian(
          active_tape_index,
          n_independent_variables,
          vector_float_to_double(independent_variables),
          hessian_double);
        hessian = hessian_double;
      }

      // === Vector drivers ===

      static void
      values (const unsigned int             &active_tape_index,
              const unsigned int             &n_dependent_variables,
              const unsigned int             &n_independent_variables,
              const std::vector<scalar_type> &independent_variables,
              Vector<scalar_type>            &values)
      {
        Vector<double> values_double (values.size());
        TapedDrivers<ADNumberType,double>::values(
          active_tape_index,
          n_dependent_variables,
          n_independent_variables,
          vector_float_to_double(independent_variables),
          values_double);
        values = values_double;
      }

      static void
      jacobian (const unsigned int             &active_tape_index,
                const unsigned int             &n_dependent_variables,
                const unsigned int             &n_independent_variables,
                const std::vector<scalar_type> &independent_variables,
                FullMatrix<scalar_type>        &jacobian)
      {
        FullMatrix<double> jacobian_double (jacobian.m(), jacobian.n());
        TapedDrivers<ADNumberType,double>::jacobian(
          active_tape_index,
          n_dependent_variables,
          n_independent_variables,
          vector_float_to_double(independent_variables),
          jacobian_double);
        jacobian = jacobian_double;
      }
    };


    // -------------   TapelessDrivers   -------------


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType,ScalarType,T>::initialize (
      const unsigned int &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    ScalarType
    TapelessDrivers<ADNumberType,ScalarType,T>::value (
      const std::vector<ADNumberType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return ScalarType(0.0);
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType,ScalarType,T>::gradient (
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType,ScalarType,T>::hessian (
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType,ScalarType,T>::values (
      const std::vector<ADNumberType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template<typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType,ScalarType,T>::jacobian (
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    namespace internal
    {

      // A dummy function to define the active dependent variable when using
      // reverse-mode AD.
      template<typename ADNumberType>
      static typename std::enable_if<
      !is_sacado_rad_number<ADNumberType>::value
      >::type
      reverse_mode_dependent_variable_activation (ADNumberType &)
      {

      }

#ifdef DEAL_II_WITH_TRILINOS


      // Define the active dependent variable when using reverse-mode AD.
      //
      // If there are multiple dependent variables then it is necessary to
      // inform the independent variables, from which the adjoints are computed,
      // which dependent variable they are computing the gradients with respect
      // to. This function broadcasts this information.
      template<typename ADNumberType>
      static typename std::enable_if<
      is_sacado_rad_number<ADNumberType>::value
      >::type
      reverse_mode_dependent_variable_activation (ADNumberType &dependent_variable)
      {
        // Compute all gradients (adjoints) for this
        // reverse-mode Sacado dependent variable.
        // For reverse-mode Sacado numbers it is necessary to broadcast to
        // all independent variables that it is time to compute gradients.
        // For one dependent variable one would just need to all
        // ADNumberType::Gradcomp(), but since we have a more
        // generic implementation for vectors of dependent variables
        // (vector mode) we default to the complex case.
        ADNumberType::Outvar_Gradcomp(dependent_variable);
      }

#endif


      // A dummy function to enable vector mode for tapeless
      // auto-differentiable numbers.
      template<typename ADNumberType>
      static typename std::enable_if<
      !(is_adolc_number<ADNumberType>::value &&
        is_tapeless_ad_number<ADNumberType>::value)
      >::type
      configure_tapeless_mode (const unsigned int)
      {

      }

#ifdef DEAL_II_WITH_ADOLC


      // Enable vector mode for Adol-C tapeless numbers.
      //
      // This function checks to see if its legal to increase the maximum
      // number of directional derivatives to be considered during calculations.
      // If not then it throws an error.
      template<typename ADNumberType>
      static typename std::enable_if<
      is_adolc_number<ADNumberType>::value &&
      is_tapeless_ad_number<ADNumberType>::value
      >::type
      configure_tapeless_mode (const unsigned int n_directional_derivatives)
      {
        // See Adol-C manual section 7.1
        //
        // NOTE: It is critical that this is done for tapeless mode BEFORE
        // any adtl::adouble are created. If this is not done, then we see
        // this scary warning:
        //
        // "
        // ADOL-C Warning: Tapeless: Setting numDir could change memory
        // allocation of derivatives in existing adoubles and may lead to
        // erroneous results or memory corruption
        // "
        //
        // So we use this dummy function to configure this setting before
        // we create and initialize our class data
        const std::size_t n_live_variables = adtl::refcounter::getNumLiveVar();
        if (n_live_variables == 0)
          {
            adtl::setNumDir(n_directional_derivatives);
          }
        else
          {
            // So there are some live active variables floating around. Here we
            // check if we ask to increase the number of number of computable
            // directional derivatives. If this really is necessary then its
            // absolutely vital that there exist no live variables before doing
            // so.
            const std::size_t n_set_directional_derivatives = adtl::getNumDir();
            if (n_directional_derivatives > n_set_directional_derivatives)
              AssertThrow(n_live_variables == 0,
                          ExcMessage("There are currently " +
                                     Utilities::to_string(n_live_variables) + " live "
                                     "adtl::adouble variables in existence. They currently "
                                     "assume " +
                                     Utilities::to_string(n_set_directional_derivatives) + " directional derivatives "
                                     "but you wish to increase this to " +
                                     Utilities::to_string(n_directional_derivatives) + ". \n"
                                     "To safely change (or more specifically in this case, "
                                     "increase) the number of directional derivatives, there "
                                     "must be no tapeless doubles in local/global scope."));
          }
      }

#endif

    }


    // Specialization for auto-differentiable numbers that use
    // reverse mode to compute the first derivatives (and, if supported,
    // forward mode for the second).
    template<typename ADNumberType, typename ScalarType>
    struct TapelessDrivers<ADNumberType,ScalarType,typename std::enable_if<
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_rad ||
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_rad_dfad
      >::type>
    {

      // === Configuration ===

      static void
      initialize (const unsigned int &n_independent_variables)
      {
        internal::configure_tapeless_mode<ADNumberType>(n_independent_variables);
      }

      // === Scalar drivers ===

      static ScalarType
      value (const std::vector<ADNumberType> &dependent_variables)
      {
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(),1));
        return ADNumberTraits<ADNumberType>::get_scalar_value(dependent_variables[0]);
      }

      static void
      gradient (const std::vector<ADNumberType> &independent_variables,
                const std::vector<ADNumberType> &dependent_variables,
                Vector<ScalarType>              &gradient)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,1));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(),1));
        Assert(gradient.size() == independent_variables.size(),
               ExcDimensionMismatch(gradient.size(),independent_variables.size()));

        // In reverse mode, the gradients are computed from the
        // independent variables (i.e. the adjoint)
        internal::reverse_mode_dependent_variable_activation(const_cast<ADNumberType &>(dependent_variables[0]));
        const std::size_t n_independent_variables = independent_variables.size();
        for (unsigned int i=0; i<n_independent_variables; i++)
          gradient[i] = internal::NumberType<ScalarType>::value(
                          ADNumberTraits<ADNumberType>::get_directional_derivative(
                            independent_variables[i],
                            0 /*This number doesn't really matter*/));
      }

      static void
      hessian (const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType>          &hessian)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,2));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(),1));
        Assert(hessian.m() == independent_variables.size(),
               ExcDimensionMismatch(hessian.m(),independent_variables.size()));
        Assert(hessian.n() == independent_variables.size(),
               ExcDimensionMismatch(hessian.n(),independent_variables.size()));

        // In reverse mode, the gradients are computed from the
        // independent variables (i.e. the adjoint)
        internal::reverse_mode_dependent_variable_activation(const_cast<ADNumberType &>(dependent_variables[0]));
        const std::size_t n_independent_variables = independent_variables.size();
        for (unsigned int i=0; i<n_independent_variables; i++)
          {
            typedef typename ADNumberTraits<ADNumberType>::derivative_type derivative_type;
            const derivative_type gradient_i
              = ADNumberTraits<ADNumberType>::get_directional_derivative(independent_variables[i], i);

            for (unsigned int j=0; j <= i; ++j) // Symmetry
              {
                // Extract higher-order directional derivatives. Depending on the AD number type,
                // the result may be another AD number or a floating point value.
                const ScalarType hessian_ij
                  = internal::NumberType<ScalarType>::value(
                      ADNumberTraits<derivative_type>::get_directional_derivative(gradient_i, j));
                hessian[i][j] = hessian_ij;
                if (i != j)
                  hessian[j][i] = hessian_ij;  // Symmetry
              }
          }
      }

      // === Vector drivers ===

      static void
      values (const std::vector<ADNumberType> &dependent_variables,
              Vector<ScalarType>              &values)
      {
        Assert(values.size() == dependent_variables.size(),
               ExcDimensionMismatch(values.size(),dependent_variables.size()));

        const std::size_t n_dependent_variables = dependent_variables.size();
        for (unsigned int i=0; i<n_dependent_variables; i++)
          values[i] = ADNumberTraits<ADNumberType>::get_scalar_value(dependent_variables[i]);
      }

      static void
      jacobian (const std::vector<ADNumberType> &independent_variables,
                const std::vector<ADNumberType> &dependent_variables,
                FullMatrix<ScalarType>          &jacobian)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,1));
        Assert(jacobian.m() == dependent_variables.size(),
               ExcDimensionMismatch(jacobian.m(),dependent_variables.size()));
        Assert(jacobian.n() == independent_variables.size(),
               ExcDimensionMismatch(jacobian.n(),independent_variables.size()));

        const std::size_t n_independent_variables = independent_variables.size();
        const std::size_t n_dependent_variables = dependent_variables.size();

        // In reverse mode, the gradients are computed from the
        // independent variables (i.e. the adjoint).
        // For a demonstration of why this accumulation process is
        // required, see the unit tests
        // sacado/basic_01b.cc and sacado/basic_02b.cc
        // Here we also take into consideration the derivative type:
        // The Sacado number may be of the nested variety, in which
        // case the effect of the accumulation process on the
        // sensitivities of the nested number need to be accounted for.
        typedef typename ADNumberTraits<ADNumberType>::derivative_type AccumulationType;
        std::vector<AccumulationType> rad_accumulation (
          n_independent_variables,
          dealii::internal::NumberType<AccumulationType>::value(0.0));

        for (unsigned int i=0; i<n_dependent_variables; i++)
          {
            internal::reverse_mode_dependent_variable_activation(
              const_cast<ADNumberType &>(dependent_variables[i]));
            for (unsigned int j=0; j<n_independent_variables; j++)
              {
                const AccumulationType df_i_dx_j
                  = ADNumberTraits<ADNumberType>::get_directional_derivative(
                      independent_variables[j], i /*This number doesn't really matter*/)
                    - rad_accumulation[j];
                jacobian[i][j] = internal::NumberType<ScalarType>::value(df_i_dx_j);
                rad_accumulation[j] += df_i_dx_j;
              }
          }
      }

    };


    // Specialization for auto-differentiable numbers that use
    // forward mode to compute the first (and, if supported, second) derivatives.
    template<typename ADNumberType, typename ScalarType>
    struct TapelessDrivers<ADNumberType,ScalarType,typename std::enable_if<
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::adolc_tapeless ||
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_dfad ||
      ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_dfad_dfad
      >::type>
    {

      // === Configuration ===

      static void
      initialize (const unsigned int &n_independent_variables)
      {
        internal::configure_tapeless_mode<ADNumberType>(n_independent_variables);
      }

      // === Scalar drivers ===

      static ScalarType
      value (const std::vector<ADNumberType> &dependent_variables)
      {
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(),1));
        return ADNumberTraits<ADNumberType>::get_scalar_value(dependent_variables[0]);
      }

      static void
      gradient (const std::vector<ADNumberType> &independent_variables,
                const std::vector<ADNumberType> &dependent_variables,
                Vector<ScalarType>              &gradient)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,1));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(),1));
        Assert(gradient.size() == independent_variables.size(),
               ExcDimensionMismatch(gradient.size(),independent_variables.size()));

        // In forward mode, the gradients are computed from the
        // dependent variables
        const std::size_t n_independent_variables = independent_variables.size();
        for (unsigned int i=0; i<n_independent_variables; i++)
          gradient[i] = internal::NumberType<ScalarType>::value(
                          ADNumberTraits<ADNumberType>::get_directional_derivative(
                            dependent_variables[0], i));
      }

      static void
      hessian (const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType>          &hessian)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,2));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(),1));
        Assert(hessian.m() == independent_variables.size(),
               ExcDimensionMismatch(hessian.m(),independent_variables.size()));
        Assert(hessian.n() == independent_variables.size(),
               ExcDimensionMismatch(hessian.n(),independent_variables.size()));

        // In forward mode, the gradients are computed from the
        // dependent variables
        const std::size_t n_independent_variables = independent_variables.size();
        for (unsigned int i=0; i<n_independent_variables; i++)
          {
            typedef typename ADNumberTraits<ADNumberType>::derivative_type derivative_type;
            const derivative_type gradient_i
              = ADNumberTraits<ADNumberType>::get_directional_derivative(dependent_variables[0], i);

            for (unsigned int j=0; j <= i; ++j) // Symmetry
              {
                // Extract higher-order directional derivatives. Depending on the AD number type,
                // the result may be another AD number or a floating point value.
                const ScalarType hessian_ij
                  = internal::NumberType<ScalarType>::value(
                      ADNumberTraits<derivative_type>::get_directional_derivative(gradient_i, j));
                hessian[i][j] = hessian_ij;
                if (i != j)
                  hessian[j][i] = hessian_ij;  // Symmetry
              }
          }
      }

      // === Vector drivers ===

      static void
      values (const std::vector<ADNumberType> &dependent_variables,
              Vector<ScalarType>              &values)
      {
        Assert(values.size() == dependent_variables.size(),
               ExcDimensionMismatch(values.size(),dependent_variables.size()));

        const std::size_t n_dependent_variables = dependent_variables.size();
        for (unsigned int i=0; i<n_dependent_variables; i++)
          values[i] = ADNumberTraits<ADNumberType>::get_scalar_value(dependent_variables[i]);
      }

      static void
      jacobian (const std::vector<ADNumberType> &independent_variables,
                const std::vector<ADNumberType> &dependent_variables,
                FullMatrix<ScalarType>          &jacobian)
      {
        Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
               ExcSupportedDerivativeLevels(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,1));
        Assert(jacobian.m() == dependent_variables.size(),
               ExcDimensionMismatch(jacobian.m(),dependent_variables.size()));
        Assert(jacobian.n() == independent_variables.size(),
               ExcDimensionMismatch(jacobian.n(),independent_variables.size()));

        const std::size_t n_independent_variables = independent_variables.size();
        const std::size_t n_dependent_variables = dependent_variables.size();

        // In forward mode, the gradients are computed from the
        // dependent variables
        for (unsigned int i=0; i<n_dependent_variables; i++)
          for (unsigned int j=0; j<n_independent_variables; j++)
            jacobian[i][j] = internal::NumberType<ScalarType>::value(
                               ADNumberTraits<ADNumberType>::get_directional_derivative(dependent_variables[i], j));
      }

    };


  } // namespace AD
} // namespace Differentiation


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE


#endif
