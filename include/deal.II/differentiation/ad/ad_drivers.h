// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_ADOLC

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <adolc/adolc_fatalerror.h>
#  include <adolc/drivers/drivers.h>
#  include <adolc/internal/usrparms.h>
#  include <adolc/taping.h>
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
    DeclExceptionMsg(
      ExcRequiresADNumberSpecialization,
      "This function is called in a class that is expected to be specialized "
      "for auto-differentiable numbers.");

    /**
     * Exception denoting that ADOL-C is a required feature.
     */
    DeclExceptionMsg(
      ExcRequiresADOLC,
      "This function is only available if deal.II is compiled with ADOL-C.");

    /**
     * This exception is raised whenever the an auto-differentiable number does
     * not support the required number of derivative operations
     *
     * The first parameter to the constructor is the number of derivative
     * operations that it provides, and the second is the minimum number that
     * are required. Both parameters are of type <tt>int</tt>.
     */
    DeclException2(
      ExcSupportedDerivativeLevels,
      std::size_t,
      std::size_t,
      << "The number of derivative levels that this auto-differentiable number type supports is "
      << arg1
      << ", but to perform the intended operation the number must support at least "
      << arg2 << " levels.");

    //@}


    /**
     * A struct defining a collection of types used within the context of
     * auto-differentiable numbers.
     */
    template <typename ADNumberType, typename T = void>
    struct Types
    {
      /**
       * Typedef for tape indices.
       */
      using tape_index = unsigned int;

      /**
       * Typedef for tape buffer sizes.
       */
      using tape_buffer_sizes = unsigned int;
    }; // struct Types


    /**
     * A struct defining a collection of special numbers used within the
     * context of auto-differentiable numbers.
     */
    template <typename ADNumberType, typename T = void>
    struct Numbers
    {
      /**
       * A tape index that is unusable and can be used to invalidate recording
       * operations.
       */
      // Note: We stipulate in the documentation for the helper classes that the
      // valid tape index range is (invalid_tape_index,max_tape_index).
      static const typename Types<ADNumberType>::tape_index invalid_tape_index =
        0u;

      /**
       * The maximum number of tapes that can be written on one process.
       */
      static const typename Types<ADNumberType>::tape_index max_tape_index =
        std::numeric_limits<typename Types<ADNumberType>::tape_index>::max();
    }; // struct Numbers


    /**
     * A driver class for taped auto-differentiable numbers.
     *
     * It is intended that this class be specialized for the valid
     * combinations of auto-differentiable numbers and output scalar
     * number types.
     *
     * @tparam ADNumberType A type corresponding to a supported
     *         auto-differentiable number.
     * @tparam ScalarType A real or complex floating point number type
     *         that is the scalar value type used for input to, and output
     *         from, operations performed with auto-differentiable numbers.
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename ADNumberType, typename ScalarType, typename T = void>
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
       * Return whether or not this class is tracking calculations performed
       * with its marked independent variables.
       */
      bool
      is_recording() const;

      /**
       * Return the tape number which is currently activated for recording or
       * reading.
       */
      typename Types<ADNumberType>::tape_index
      active_tape_index() const;

      /**
       * Return whether or not a tape number has already been used
       * or registered.
       */
      bool
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      /**
       * Return whether or not the numerical values of all independent
       * variables are recorded in the tape buffer.
       */
      bool
      keep_independent_values() const;

      /**
       * Set the buffer sizes for the next active tape.
       *
       * This function must be called before start_recording_operations() for
       * it to have any influence on the memory allocated to the next recorded
       * tape.
       *
       * @note This function only has an effect when using ADOL-C numbers. As
       * stated by the ADOL-C manual, it may be desirable to create a file
       * ".adolcrc" in the program run directory and set the buffer size
       * therein.
       * Alternatively, this function can be used to override the settings for
       * any given tape, or can be used in the event that no ".adolcrc" file
       * exists.
       * The default value for each buffer is set at 64 megabytes, a
       * heuristically chosen value thought to be appropriate for use within
       * the context of finite element analysis when considering coupled
       * problems with multiple vector-valued fields discretised by higher
       * order shape functions, as well as complex constitutive laws.
       *
       * @param[in] obufsize ADOL-C operations buffer size
       * @param[in] lbufsize ADOL-C locations buffer size
       * @param[in] vbufsize ADOL-C value buffer size
       * @param[in] tbufsize ADOL-C Taylor buffer size
       */
      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes obufsize = 64 *
                                                                         1024 *
                                                                         1024,
        const typename Types<ADNumberType>::tape_buffer_sizes lbufsize = 64 *
                                                                         1024 *
                                                                         1024,
        const typename Types<ADNumberType>::tape_buffer_sizes vbufsize = 64 *
                                                                         1024 *
                                                                         1024,
        const typename Types<ADNumberType>::tape_buffer_sizes tbufsize = 64 *
                                                                         1024 *
                                                                         1024);

      /**
       * Enable the recording mode for a given tape.
       *
       * @param[in] tape_index The index of the tape to be written.
       * @param[in] keep_independent_values Determines whether the numerical
       *            values of all independent variables are recorded in the
       *            tape buffer.
       */
      void
      start_taping(const typename Types<ADNumberType>::tape_index tape_index,
                   const bool keep_independent_values);

      /**
       * Disable the recording mode for a given tape.
       *
       * @param[in] active_tape_index The index of the (currently active) tape
       *            to be finalized and potentially written to file.
       * @param[in] write_tapes_to_file A flag that specified whether the tape
       *            should be written to file or kept in memory.
       */
      void
      stop_taping(
        const typename Types<ADNumberType>::tape_index active_tape_index,
        const bool                                     write_tapes_to_file);

      /**
       * Select a tape to record to or read from.
       *
       * This function activates a tape, but depending on whether @p read_mode
       * is set, the tape is either taken as previously written to (and put
       * into read-only mode), or cleared for (re-)taping.
       *
       * @param[in] tape_index The index of the tape to be written to/read
       *            from.
       *
       * @note The chosen tape index must be greater than
       * Numbers<ADNumberType>::invalid_tape_index and less than
       * Numbers<ADNumberType>::max_tape_index.
       */
      void
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index);

      /**
       * Reset the state of the class.
       *
       * @note This also resets the active tape number to an invalid number, and
       * deactivates the recording mode for taped variables.
       */
      void
      reset(const bool clear_registered_tapes);

      /**
       * Print the status of all queriable data. Exactly what is printed and
       * its format depends on the @p ad_type, as is determined by the
       * @p ADNumberTypeCode template parameter.
       *
       * @param[in] stream The output stream to which the values are to be
       * written.
       */
      void
      print(std::ostream &stream) const;

      /**
       * Print the statistics regarding the usage of the tapes.
       *
       * @param[in] tape_index The index of the tape to get the statistics of.
       * @param[out] stream The output stream to which the values are to be
       *            written.
       */
      static void
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream &                                 stream);

      //@}

      /**
       * @name Drivers for scalar functions (one dependent variable)
       */
      //@{

      /**
       * Compute the value of the scalar field.
       *
       * @param[in] active_tape_index The index of the tape on which the
       *            dependent function is recorded.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       *
       * @return The scalar value of the function.
       */
      static ScalarType
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<ScalarType> &independent_variables);

      /**
       * Compute the gradient of the scalar field with respect to all
       * independent variables.
       *
       * @param[in] active_tape_index The index of the tape on which the
       *            dependent function is recorded.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] gradient The values of the dependent function's
       *             gradients. It is expected that this vector be of the
       *             correct size (with length
       *             <code>n_independent_variables</code>).
       */
      static void
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<ScalarType> &independent_variables,
               Vector<ScalarType> &           gradient);

      /**
       * Compute the Hessian of the scalar field with respect to all
       * independent variables.
       *
       * @param[in] active_tape_index The index of the tape on which the
       *            dependent function is recorded.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] hessian The values of the dependent function's
       *             Hessian. It is expected that this matrix be of the correct
       *             size (with dimensions
       *             <code>n_independent_variables</code>$\times$<code>n_independent_variables</code>).
       */
      static void
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<ScalarType> &independent_variables,
              FullMatrix<ScalarType> &       hessian);

      //@}

      /**
       * @name Drivers for vector functions (multiple dependent variables)
       */
      //@{

      /**
       * Compute the values of the vector field.
       *
       * @param[in] active_tape_index The index of the tape on which the
       *            dependent function is recorded.
       * @param[in] n_dependent_variables The number of dependent variables.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] values The component values of the dependent functions.
       *             It is expected that this vector be of the correct size
       *             (with length <code>n_dependent_variables</code>).
       */
      static void
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int             n_dependent_variables,
             const std::vector<ScalarType> &independent_variables,
             Vector<ScalarType> &           values);

      /**
       * Compute the Jacobian of the vector field.
       *
       * The Jacobian of a vector field is in essence the gradient of each
       * dependent variable with respect to all independent variables.
       * This operation is therefore analogous to the gradient() operation
       * performed on a collection of scalar valued fields.
       *
       * @param[in] active_tape_index The index of the tape on which the
       *            dependent function is recorded.
       * @param[in] n_dependent_variables The number of dependent variables.
       * @param[in] independent_variables The scalar values of the independent
       *            variables whose sensitivities were tracked.
       * @param[out] jacobian The component values of the dependent functions'
       *             Jacobian. It is expected that this matrix be of the correct
       *             size (with dimensions
       *             <code>n_dependent_variables</code>$\times$<code>n_independent_variables</code>).
       */
      static void
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int             n_dependent_variables,
               const std::vector<ScalarType> &independent_variables,
               FullMatrix<ScalarType> &       jacobian);

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
     * @tparam ScalarType A real or complex floating point number type
     *         that is the scalar value type used for input to, and output
     *         from, operations performed with auto-differentiable numbers.
     * @tparam T An arbitrary type resulting from the application of
     *         the SFINAE idiom to selectively specialize this class.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <typename ADNumberType, typename ScalarType, typename T = void>
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
       * of how many directional derivatives might need to be computed, this
       * function informs the auto-differention library of what this number
       * is.
       *
       * @param[in] n_independent_variables The number of independent variables
       *            that will be used for the entire duration of the
       *            simulation.
       *
       * @warning For ADOL-C tapeless numbers, the value given to
       * @p n_independent_variables should be the <b>maximum</b> number of
       * independent variables that will be used for the entire duration of
       * the simulation. This is important in the context of, for example,
       * hp-FEM and for multiple constitutive models with a different number of
       * fields from which a linearization must be computed.
       */
      static void
      initialize_global_environment(const unsigned int n_independent_variables);

      //@}

      /**
       * Operation status
       */
      //@{

      /**
       * Set a flag that states that we can safely mark dependent variables
       * within the current phase of operations.
       */
      void
      allow_dependent_variable_marking();

      /**
       * Set a flag that states that we cannot safely mark dependent variables
       * within the current phase of operations.
       */
      void
      prevent_dependent_variable_marking();

      /**
       * Query a flag as to whether or not dependent variables can be marked
       * within the current phase of operations.
       */
      bool
      is_dependent_variable_marking_allowed() const;

      //@}

      /**
       * @name Drivers for scalar functions
       */
      //@{

      /**
       * Compute the value of the scalar field.
       *
       * @param[in] dependent_variables The dependent variables whose values are
       *            to be extracted.
       *
       * @return The scalar value of the function.
       */
      static ScalarType
      value(const std::vector<ADNumberType> &dependent_variables);

      /**
       * Compute the gradient of the scalar field with respect to all
       * independent variables.
       *
       * @param[in] independent_variables The independent variables whose
       *            sensitivities were tracked.
       * @param[in] dependent_variables The (single) dependent variable whose
       *            gradients are to be extracted.
       * @param[out] gradient The values of the dependent function's
       *             gradients. It is expected that this vector be of the
       *             correct size (with length
       *             <code>n_independent_variables</code>).
       */
      static void
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType> &             gradient);

      /**
       * Compute the Hessian of the scalar field with respect to all
       * independent variables.
       *
       * @param[in] independent_variables The independent variables whose
       *            sensitivities were tracked.
       * @param[in] dependent_variables The (single) dependent variable whose
       *            Hessians are to be extracted.
       * @param[out] hessian The values of the dependent function's
       *             Hessian. It is expected that this matrix be of the correct
       *             size (with dimensions
       *             <code>n_independent_variables</code>$\times$<code>n_independent_variables</code>).
       */
      static void
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType> &         hessian);

      //@}

      /**
       * @name Drivers for vector functions
       */
      //@{

      /**
       * Compute the values of the vector field.
       *
       * @param[in] dependent_variables The dependent variables whose Hessians
       *            are to be extracted.
       * @param[out] values The component values of the dependent functions.
       *             It is expected that this vector be of the correct size
       *             (with length <code>n_dependent_variables</code>).
       */
      static void
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType> &             values);

      /**
       * Compute the Jacobian of the vector field.
       *
       * The Jacobian of a vector field is in essence the gradient of each
       * dependent variable with respect to all independent variables.
       * This operation is therefore analogous to the gradient() operation
       * performed on a collection of scalar valued fields.
       *
       * @param[in] independent_variables The independent variables whose
       *            sensitivities were tracked.
       * @param[in] dependent_variables The dependent variables whose Jacobian
       *            are to be extracted.
       * @param[out] jacobian The component values of the dependent functions'
       *             Jacobian. It is expected that this matrix be of the correct
       *             size (with dimensions
       *             <code>n_dependent_variables</code>$\times$<code>n_independent_variables</code>).
       */
      static void
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType> &         jacobian);

      //@}
    };

  } // namespace AD
} // namespace Differentiation



/* --------------------- inline and template functions --------------------- */


#ifndef DOXYGEN

namespace Differentiation
{
  namespace AD
  {
    // -------------   TapedDrivers   -------------


    template <typename ADNumberType, typename ScalarType, typename T>
    bool
    TapedDrivers<ADNumberType, ScalarType, T>::is_recording() const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return false;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    typename Types<ADNumberType>::tape_index
    TapedDrivers<ADNumberType, ScalarType, T>::active_tape_index() const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return Numbers<ADNumberType>::invalid_tape_index;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    bool
    TapedDrivers<ADNumberType, ScalarType, T>::is_registered_tape(
      const typename Types<ADNumberType>::tape_index) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return false;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    bool
    TapedDrivers<ADNumberType, ScalarType, T>::keep_independent_values() const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return false;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::set_tape_buffer_sizes(
      const typename Types<ADNumberType>::tape_buffer_sizes,
      const typename Types<ADNumberType>::tape_buffer_sizes,
      const typename Types<ADNumberType>::tape_buffer_sizes,
      const typename Types<ADNumberType>::tape_buffer_sizes)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::start_taping(
      const typename Types<ADNumberType>::tape_index,
      const bool)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::stop_taping(
      const typename Types<ADNumberType>::tape_index,
      const bool)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::activate_tape(
      const typename Types<ADNumberType>::tape_index)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::reset(const bool)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::print(std::ostream &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::print_tape_stats(
      const typename Types<ADNumberType>::tape_index,
      std::ostream &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    ScalarType
    TapedDrivers<ADNumberType, ScalarType, T>::value(
      const typename Types<ADNumberType>::tape_index,
      const std::vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return ScalarType(0.0);
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::gradient(
      const typename Types<ADNumberType>::tape_index,
      const std::vector<ScalarType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::hessian(
      const typename Types<ADNumberType>::tape_index,
      const std::vector<ScalarType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::values(
      const typename Types<ADNumberType>::tape_index,
      const unsigned int,
      const std::vector<ScalarType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::jacobian(
      const typename Types<ADNumberType>::tape_index,
      const unsigned int,
      const std::vector<ScalarType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


#  ifdef DEAL_II_WITH_ADOLC

    /**
     * Specialization for taped ADOL-C auto-differentiable numbers.
     */
    template <typename ADNumberType>
    struct Types<
      ADNumberType,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_taped>::type>
    {
      /**
       * Typedef for tape indices. ADOL-C uses short integers, so
       * we restrict ourselves to similar types.
       */
      using tape_index = unsigned short;

      /**
       * Typedef for tape buffer sizes.
       */
      using tape_buffer_sizes = unsigned int;
    }; // struct Types


    /**
     * Specialization for taped ADOL-C auto-differentiable numbers.
     */
    template <typename ADNumberType>
    struct Numbers<
      ADNumberType,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_taped>::type>
    {
      /**
       * A tape index that is unusable and can be used to invalidate recording
       * operations.
       *
       * @note ADOL-C doesn't allow us to record to this reserved tape (i.e. can't
       * write it to file), so we can safely use it as an invalidation case. In
       * general, we want the user to be able to record to a tape if they'd
       * like.
       */
      static const typename Types<ADNumberType>::tape_index invalid_tape_index =
        0;

      /**
       * The maximum number of tapes that can be written on one process.
       */
      // Note: This value is a limitation of ADOL-C, and not something that we
      // have control over. See test adolc/helper_tape_index_01.cc for
      // verification that we cannot use or exceed this value. This value is
      // defined as TBUFNUM; see
      // https://gitlab.com/adol-c/adol-c/blob/master/ADOL-C/include/adolc/internal/usrparms.h#L34
      static const typename Types<ADNumberType>::tape_index max_tape_index =
        TBUFNUM;
    }; // struct Numbers


    /**
     * Specialization for taped ADOL-C auto-differentiable numbers.
     *
     * Note: In the case of ADOL-C taped numbers, the associated scalar
     * type is always expected to be a double. So we need to make a further
     * specialization when ScalarType is a float.
     */
    template <typename ADNumberType>
    struct TapedDrivers<
      ADNumberType,
      double,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_taped>::type>
    {
      using scalar_type = double;

      TapedDrivers()
        : active_tape(Numbers<ADNumberType>::invalid_tape_index)
        , keep_values(true)
        , is_recording_flag(false)
        , use_stored_taped_buffer_sizes(false)
        , obufsize(0u)
        , lbufsize(0u)
        , vbufsize(0u)
        , tbufsize(0u)
      {}

      // === Taping ===

      bool
      is_recording() const
      {
        return is_recording_flag;
      }

      typename Types<ADNumberType>::tape_index
      active_tape_index() const
      {
        return active_tape;
      }

      bool
      keep_independent_values() const
      {
        return keep_values;
      }

      bool
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const
      {
        // See https://gitlab.com/adol-c/adol-c/issues/11#note_108341333
        try
          {
            std::vector<std::size_t> counts(STAT_SIZE);
            ::tapestats(tape_index, counts.data());
            return true;
          }
        catch (const ::FatalError &exc)
          {
            return false;
          }
      }

      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes in_obufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes in_lbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes in_vbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes in_tbufsize)
      {
        // When valid for the chosen AD number type, these values will be used
        // the next time start_recording_operations() is called.
        obufsize                      = in_obufsize;
        lbufsize                      = in_lbufsize;
        vbufsize                      = in_vbufsize;
        tbufsize                      = in_tbufsize;
        use_stored_taped_buffer_sizes = true;
      }

      void
      start_taping(const typename Types<ADNumberType>::tape_index tape_index,
                   const bool keep_independent_values)
      {
        if (use_stored_taped_buffer_sizes)
          trace_on(tape_index,
                   keep_independent_values,
                   obufsize,
                   lbufsize,
                   vbufsize,
                   tbufsize);
        else
          trace_on(tape_index, keep_independent_values);

        // Set some other flags to their indicated / required values
        keep_values       = keep_independent_values;
        is_recording_flag = true;
      }

      void
      stop_taping(
        const typename Types<ADNumberType>::tape_index active_tape_index,
        const bool                                     write_tapes_to_file)
      {
        if (write_tapes_to_file)
          trace_off(active_tape_index); // Slow
        else
          trace_off(); // Fast(er)

        // Now that we've turned tracing off, we've definitely
        // stopped all tape recording.
        is_recording_flag = false;

        // If the keep_values flag is set, then we expect the user to use this
        // tape immediately after recording it. There is therefore no need to
        // invalidate it. However, there is now also no way to double-check
        // that the newly recorded tape is indeed the active tape.
        if (keep_independent_values() == false)
          active_tape = Numbers<ADNumberType>::invalid_tape_index;
      }

      void
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index)
      {
        active_tape = tape_index;
      }

      void
      reset(const bool clear_registered_tapes)
      {
        active_tape       = Numbers<ADNumberType>::invalid_tape_index;
        is_recording_flag = false;
        if (clear_registered_tapes)
          {
            const std::vector<typename Types<ADNumberType>::tape_index>
              registered_tape_indices = get_registered_tape_indices();
            for (const auto &tape_index : registered_tape_indices)
              removeTape(tape_index, TapeRemovalType::ADOLC_REMOVE_COMPLETELY);
          }
      }

      void
      print(std::ostream &stream) const
      {
        const std::vector<typename Types<ADNumberType>::tape_index>
          registered_tape_indices = get_registered_tape_indices();
        stream << "Registered tapes: ";
        auto it_registered_tape = registered_tape_indices.begin();
        for (unsigned int i = 0; i < registered_tape_indices.size();
             ++i, ++it_registered_tape)
          stream << *it_registered_tape
                 << (i < (registered_tape_indices.size() - 1) ? "," : "");
        stream << "\n";

        stream << "Keep values? " << keep_independent_values() << "\n";
        stream << "Use stored tape buffer sizes? "
               << use_stored_taped_buffer_sizes << "\n";
      }

      static void
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream &                                 stream)
      {
        // See ADOL-C manual section 2.1
        // and adolc/taping.h
        std::vector<std::size_t> counts(STAT_SIZE);
        ::tapestats(tape_index, counts.data());
        Assert(counts.size() >= 18, ExcInternalError());
        stream
          << "Tape index: " << tape_index << "\n"
          << "Number of independent variables: " << counts[0] << "\n"
          << "Number of dependent variables:   " << counts[1] << "\n"
          << "Max number of live, active variables: " << counts[2] << "\n"
          << "Size of taylor stack (number of overwrites): " << counts[3]
          << "\n"
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
          << "Use of 'min_op', deferred to 'abs_op' for piecewise calculations: "
          << counts[15] << "\n"
          << "Number of 'abs' calls that can switch branch: " << counts[16]
          << "\n"
          << "Number of parameters (doubles) interchangeable without retaping: "
          << counts[17] << "\n"
          << std::flush;
      }


      // === Scalar drivers ===

      static scalar_type
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<scalar_type> &independent_variables)
      {
        scalar_type value = 0.0;

        ::function(active_tape_index,
                   1, // Only one dependent variable
                   independent_variables.size(),
                   const_cast<double *>(independent_variables.data()),
                   &value);

        return value;
      }

      static void
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<scalar_type> &independent_variables,
               Vector<scalar_type> &           gradient)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            1));
        Assert(gradient.size() == independent_variables.size(),
               ExcDimensionMismatch(gradient.size(),
                                    independent_variables.size()));

        // Note: ADOL-C's ::gradient function expects a *double as the last
        // parameter. Here we take advantage of the fact that the data in the
        // Vector class is aligned (e.g. stored as an Array)
        ::gradient(active_tape_index,
                   independent_variables.size(),
                   const_cast<scalar_type *>(independent_variables.data()),
                   &gradient[0]);
      }

      static void
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<scalar_type> &independent_variables,
              FullMatrix<scalar_type> &       hessian)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            2));
        Assert(hessian.m() == independent_variables.size(),
               ExcDimensionMismatch(hessian.m(), independent_variables.size()));
        Assert(hessian.n() == independent_variables.size(),
               ExcDimensionMismatch(hessian.n(), independent_variables.size()));

        const unsigned int n_independent_variables =
          independent_variables.size();
        std::vector<scalar_type *> H(n_independent_variables);
        for (unsigned int i = 0; i < n_independent_variables; ++i)
          H[i] = &(hessian[i][0]);

        ::hessian(active_tape_index,
                  n_independent_variables,
                  const_cast<scalar_type *>(independent_variables.data()),
                  H.data());

        // ADOL-C builds only the lower-triangular part of the
        // symmetric Hessian, so we should copy the relevant
        // entries into the upper triangular part.
        for (unsigned int i = 0; i < n_independent_variables; i++)
          for (unsigned int j = 0; j < i; j++)
            hessian[j][i] = hessian[i][j]; // Symmetry
      }

      // === Vector drivers ===

      static void
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int              n_dependent_variables,
             const std::vector<scalar_type> &independent_variables,
             Vector<scalar_type> &           values)
      {
        Assert(values.size() == n_dependent_variables,
               ExcDimensionMismatch(values.size(), n_dependent_variables));

        // Note: ADOL-C's ::function function expects a *double as the last
        // parameter. Here we take advantage of the fact that the data in the
        // Vector class is aligned (e.g. stored as an Array)
        ::function(active_tape_index,
                   n_dependent_variables,
                   independent_variables.size(),
                   const_cast<scalar_type *>(independent_variables.data()),
                   &values[0]);
      }

      static void
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int              n_dependent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type> &       jacobian)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            1));
        Assert(jacobian.m() == n_dependent_variables,
               ExcDimensionMismatch(jacobian.m(), n_dependent_variables));
        Assert(jacobian.n() == independent_variables.size(),
               ExcDimensionMismatch(jacobian.n(),
                                    independent_variables.size()));

        std::vector<scalar_type *> J(n_dependent_variables);
        for (unsigned int i = 0; i < n_dependent_variables; ++i)
          J[i] = &(jacobian[i][0]);

        ::jacobian(active_tape_index,
                   n_dependent_variables,
                   independent_variables.size(),
                   independent_variables.data(),
                   J.data());
      }

    protected:
      /**
       * Index of the tape that is currently in use. It is this tape that will
       * be recorded to or read from when performing computations using "taped"
       * auto-differentiable numbers.
       */
      typename Types<ADNumberType>::tape_index active_tape;

      /**
       * Mark whether we're going to inform taped data structures to retain
       * the coefficients ("Taylors" in ADOL-C nomenclature) stored on the
       * tape so that they can be evaluated again at a later stage.
       */
      bool keep_values;

      /**
       * Mark whether we're currently recording a tape. Dependent on the state
       * of this flag, only a restricted set of operations are allowable.
       */
      bool is_recording_flag;

      /**
       * A flag indicating that we should preferentially use the user-defined
       * taped buffer sizes as opposed to either the default values selected
       * by the AD library (or, in the case of ADOL-C, defined in an
       * ".adolcrc" file).
       */
      bool use_stored_taped_buffer_sizes;

      /**
       * ADOL-C operations buffer size.
       */
      typename Types<ADNumberType>::tape_buffer_sizes obufsize;

      /**
       * ADOL-C locations buffer size.
       */
      typename Types<ADNumberType>::tape_buffer_sizes lbufsize;

      /**
       * ADOL-C value buffer size.
       */
      typename Types<ADNumberType>::tape_buffer_sizes vbufsize;

      /**
       * ADOL-C Taylor buffer size.
       */
      typename Types<ADNumberType>::tape_buffer_sizes tbufsize;

      /**
       * Return a list of registered tape indices
       */
      std::vector<typename Types<ADNumberType>::tape_index>
      get_registered_tape_indices() const
      {
        // We've chosen to use unsigned shorts for the tape
        // index type (a safety precaution) so we need to
        // perform a conversion betwwen ADOL-C's native tape
        // index type and that chosen by us.
        std::vector<short> registered_tape_indices_s;
        cachedTraceTags(registered_tape_indices_s);

        std::vector<typename Types<ADNumberType>::tape_index>
          registered_tape_indices(registered_tape_indices_s.size());
        std::copy(registered_tape_indices_s.begin(),
                  registered_tape_indices_s.end(),
                  registered_tape_indices.begin());
        return registered_tape_indices;
      }
    };

#  else // DEAL_II_WITH_ADOLC

    /**
     * Specialization for taped ADOL-C auto-differentiable numbers.
     *
     * Although we could revert to the default definition for the
     * unspecialized TapedDrivers class, we add this specialization
     * to provide a more descriptive error message if any of its
     * static member functions are called.
     */
    template <typename ADNumberType>
    struct TapedDrivers<
      ADNumberType,
      double,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_taped>::type>
    {
      using scalar_type = double;

      // === Taping ===

      bool
      is_recording() const
      {
        AssertThrow(false, ExcRequiresADOLC());
        return false;
      }

      typename Types<ADNumberType>::tape_index
      active_tape_index() const
      {
        AssertThrow(false, ExcRequiresADOLC());
        return Numbers<ADNumberType>::invalid_tape_index;
      }

      bool
      keep_independent_values() const
      {
        AssertThrow(false, ExcRequiresADOLC());
        return false;
      }

      bool
      is_registered_tape(const typename Types<ADNumberType>::tape_index) const
      {
        AssertThrow(false, ExcRequiresADOLC());
        return false;
      }

      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      void
      start_taping(const typename Types<ADNumberType>::tape_index, const bool)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      void
      stop_taping(const typename Types<ADNumberType>::tape_index, const bool)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      void
      activate_tape(const typename Types<ADNumberType>::tape_index)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      void
      reset(const bool)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      void
      print(std::ostream &) const
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      static void
      print_tape_stats(const typename Types<ADNumberType>::tape_index,
                       std::ostream &)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }


      // === Scalar drivers ===

      static scalar_type
      value(const typename Types<ADNumberType>::tape_index,
            const std::vector<scalar_type> &)
      {
        AssertThrow(false, ExcRequiresADOLC());
        return 0.0;
      }

      static void
      gradient(const typename Types<ADNumberType>::tape_index,
               const std::vector<scalar_type> &,
               Vector<scalar_type> &)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      static void
      hessian(const typename Types<ADNumberType>::tape_index,
              const std::vector<scalar_type> &,
              FullMatrix<scalar_type> &)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      // === Vector drivers ===

      static void
      values(const typename Types<ADNumberType>::tape_index,
             const unsigned int,
             const std::vector<scalar_type> &,
             Vector<scalar_type> &)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

      static void
      jacobian(const typename Types<ADNumberType>::tape_index,
               const unsigned int,
               const std::vector<scalar_type> &,
               FullMatrix<scalar_type> &)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }
    };

#  endif // DEAL_II_WITH_ADOLC


    /**
     * Specialization for ADOL-C taped numbers. It is expected that the
     * scalar return type for this class is a float.
     *
     * @note ADOL-C only has drivers for doubles, and so floats are
     * not intrinsically supported. This wrapper struct works around
     * the issue when necessary.
     */
    template <typename ADNumberType>
    struct TapedDrivers<
      ADNumberType,
      float,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_taped>::type>
    {
      using scalar_type = float;

      // === Taping ===

      bool
      is_recording() const
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        return taped_driver.is_recording();
      }

      typename Types<ADNumberType>::tape_index
      active_tape_index() const
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        return taped_driver.active_tape_index();
      }

      bool
      keep_independent_values() const
      {
        return taped_driver.keep_independent_values();
      }

      bool
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        return taped_driver.is_registered_tape(tape_index);
      }

      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes obufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes lbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes vbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes tbufsize)
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        taped_driver.set_tape_buffer_sizes(obufsize,
                                           lbufsize,
                                           vbufsize,
                                           tbufsize);
      }

      void
      start_taping(const typename Types<ADNumberType>::tape_index tape_index,
                   const bool keep_independent_values)
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        taped_driver.start_taping(tape_index, keep_independent_values);
      }

      void
      stop_taping(
        const typename Types<ADNumberType>::tape_index active_tape_index,
        const bool                                     write_tapes_to_file)
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        taped_driver.stop_taping(active_tape_index, write_tapes_to_file);
      }

      void
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index)
      {
        taped_driver.activate_tape(tape_index);
      }

      void
      reset(const bool clear_registered_tapes)
      {
        taped_driver.reset(clear_registered_tapes);
      }

      void
      print(std::ostream &stream) const
      {
        taped_driver.print(stream);
      }

      static void
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream &                                 stream)
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        TapedDrivers<ADNumberType, double>::print_tape_stats(tape_index,
                                                             stream);
      }

      // === Scalar drivers ===

      static scalar_type
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<scalar_type> &independent_variables)
      {
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        return TapedDrivers<ADNumberType, double>::value(
          active_tape_index, vector_float_to_double(independent_variables));
      }

      static void
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<scalar_type> &independent_variables,
               Vector<scalar_type> &           gradient)
      {
        Vector<double> gradient_double(gradient.size());
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        TapedDrivers<ADNumberType, double>::gradient(active_tape_index,
                                                     vector_float_to_double(
                                                       independent_variables),
                                                     gradient_double);
        gradient = gradient_double;
      }

      static void
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<scalar_type> &independent_variables,
              FullMatrix<scalar_type> &       hessian)
      {
        FullMatrix<double> hessian_double(hessian.m(), hessian.n());
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        TapedDrivers<ADNumberType, double>::hessian(active_tape_index,
                                                    vector_float_to_double(
                                                      independent_variables),
                                                    hessian_double);
        hessian = hessian_double;
      }

      // === Vector drivers ===

      static void
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int              n_dependent_variables,
             const std::vector<scalar_type> &independent_variables,
             Vector<scalar_type> &           values)
      {
        Vector<double> values_double(values.size());
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        TapedDrivers<ADNumberType, double>::values(active_tape_index,
                                                   n_dependent_variables,
                                                   vector_float_to_double(
                                                     independent_variables),
                                                   values_double);
        values = values_double;
      }

      static void
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int              n_dependent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type> &       jacobian)
      {
        FullMatrix<double> jacobian_double(jacobian.m(), jacobian.n());
        // ADOL-C only supports 'double', not 'float', so we can forward to
        // the 'double' implementation of this function
        TapedDrivers<ADNumberType, double>::jacobian(active_tape_index,
                                                     n_dependent_variables,
                                                     vector_float_to_double(
                                                       independent_variables),
                                                     jacobian_double);
        jacobian = jacobian_double;
      }

    private:
      /**
       * Copy a vector of floats into a vector of doubles
       */
      static std::vector<double>
      vector_float_to_double(const std::vector<float> &in)
      {
        std::vector<double> out(in.size());
        std::copy(in.begin(), in.end(), out.begin());
        return out;
      }

      /**
       * The object that actually takes care of the taping
       */
      TapedDrivers<ADNumberType, double> taped_driver;
    };


    // -------------   TapelessDrivers   -------------


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::initialize_global_environment(
      const unsigned int)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }

    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::
      allow_dependent_variable_marking()
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }

    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::
      prevent_dependent_variable_marking()
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }

    template <typename ADNumberType, typename ScalarType, typename T>
    bool
    TapelessDrivers<ADNumberType, ScalarType, T>::
      is_dependent_variable_marking_allowed() const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return false;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    ScalarType
    TapelessDrivers<ADNumberType, ScalarType, T>::value(
      const std::vector<ADNumberType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return ScalarType(0.0);
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::gradient(
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::hessian(
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::values(
      const std::vector<ADNumberType> &,
      Vector<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::jacobian(
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      FullMatrix<ScalarType> &)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    namespace internal
    {
      /**
       * A dummy function to define the active dependent variable when using
       * reverse-mode AD.
       */
      template <typename ADNumberType>
      typename std::enable_if<!(ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::sacado_rad ||
                                ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::sacado_rad_dfad)>::type
      reverse_mode_dependent_variable_activation(ADNumberType &)
      {}

#  ifdef DEAL_II_TRILINOS_WITH_SACADO


      /**
       * Define the active dependent variable when using reverse-mode AD.
       *
       * If there are multiple dependent variables then it is necessary to
       * inform the independent variables, from which the adjoints are computed,
       * which dependent variable they are computing the gradients with respect
       * to. This function broadcasts this information.
       */
      template <typename ADNumberType>
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::sacado_rad ||
                              ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::sacado_rad_dfad>::type
      reverse_mode_dependent_variable_activation(
        ADNumberType &dependent_variable)
      {
        // Compute all gradients (adjoints) for this
        // reverse-mode Sacado dependent variable.
        // For reverse-mode Sacado numbers it is necessary to broadcast to
        // all independent variables that it is time to compute gradients.
        // For one dependent variable one would just need to call
        // ADNumberType::Gradcomp(), but since we have a more
        // generic implementation for vectors of dependent variables
        // (vector mode) we default to the complex case.
        ADNumberType::Outvar_Gradcomp(dependent_variable);
      }

#  endif


      /**
       * A dummy function to enable vector mode for tapeless
       * auto-differentiable numbers.
       */
      template <typename ADNumberType>
      typename std::enable_if<!(ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::adolc_tapeless)>::type
      configure_tapeless_mode(const unsigned int)
      {}


#  ifdef DEAL_II_WITH_ADOLC


      /**
       * Enable vector mode for ADOL-C tapeless numbers.
       *
       * This function checks to see if its legal to increase the maximum
       * number of directional derivatives to be considered during calculations.
       * If not then it throws an error.
       */
      template <typename ADNumberType>
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_tapeless>::type
      configure_tapeless_mode(const unsigned int n_directional_derivatives)
      {
#    ifdef DEAL_II_ADOLC_WITH_TAPELESS_REFCOUNTING
        // See ADOL-C manual section 7.1
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
            // check if we ask to increase the number of computable
            // directional derivatives. If this really is necessary then it's
            // absolutely vital that there exist no live variables before doing
            // so.
            const std::size_t n_set_directional_derivatives = adtl::getNumDir();
            if (n_directional_derivatives > n_set_directional_derivatives)
              AssertThrow(
                n_live_variables == 0,
                ExcMessage(
                  "There are currently " +
                  Utilities::to_string(n_live_variables) +
                  " live "
                  "adtl::adouble variables in existence. They currently "
                  "assume " +
                  Utilities::to_string(n_set_directional_derivatives) +
                  " directional derivatives "
                  "but you wish to increase this to " +
                  Utilities::to_string(n_directional_derivatives) +
                  ". \n"
                  "To safely change (or more specifically in this case, "
                  "increase) the number of directional derivatives, there "
                  "must be no tapeless doubles in local/global scope."));
          }
#    else
        // If ADOL-C is not configured with tapeless number reference counting
        // then there is no way that we can guarantee that the following call
        // is safe. No comment... :-/
        adtl::setNumDir(n_directional_derivatives);
#    endif
      }

#  else // DEAL_II_WITH_ADOLC

      template <typename ADNumberType>
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_tapeless>::type
      configure_tapeless_mode(const unsigned int /*n_directional_derivatives*/)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

#  endif

    } // namespace internal


    /**
     * Specialization for auto-differentiable numbers that use
     * reverse mode to compute the first derivatives (and, if supported,
     * forward mode for the second).
     */
    template <typename ADNumberType, typename ScalarType>
    struct TapelessDrivers<
      ADNumberType,
      ScalarType,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::sacado_rad ||
                              ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::sacado_rad_dfad>::type>
    {
      TapelessDrivers()
        : dependent_variable_marking_safe(false)
      {}

      // === Configuration ===

      static void
      initialize_global_environment(const unsigned int n_independent_variables)
      {
        internal::configure_tapeless_mode<ADNumberType>(
          n_independent_variables);
      }

      // === Operation status ===
      void
      allow_dependent_variable_marking()
      {
        dependent_variable_marking_safe = true;
      }

      void
      prevent_dependent_variable_marking()
      {
        dependent_variable_marking_safe = false;
      }

      bool
      is_dependent_variable_marking_allowed() const
      {
        return dependent_variable_marking_safe;
      }

      // === Scalar drivers ===

      static ScalarType
      value(const std::vector<ADNumberType> &dependent_variables)
      {
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(), 1));
        return ADNumberTraits<ADNumberType>::get_scalar_value(
          dependent_variables[0]);
      }

      static void
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType> &             gradient)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            1));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(), 1));
        Assert(gradient.size() == independent_variables.size(),
               ExcDimensionMismatch(gradient.size(),
                                    independent_variables.size()));

        // In reverse mode, the gradients are computed from the
        // independent variables (i.e. the adjoint)
        internal::reverse_mode_dependent_variable_activation(
          const_cast<ADNumberType &>(dependent_variables[0]));
        const std::size_t n_independent_variables =
          independent_variables.size();
        for (unsigned int i = 0; i < n_independent_variables; i++)
          gradient[i] = internal::NumberType<ScalarType>::value(
            ADNumberTraits<ADNumberType>::get_directional_derivative(
              independent_variables[i],
              0 /*This number doesn't really matter*/));
      }

      static void
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType> &         hessian)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            2));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(), 1));
        Assert(hessian.m() == independent_variables.size(),
               ExcDimensionMismatch(hessian.m(), independent_variables.size()));
        Assert(hessian.n() == independent_variables.size(),
               ExcDimensionMismatch(hessian.n(), independent_variables.size()));

        // In reverse mode, the gradients are computed from the
        // independent variables (i.e. the adjoint)
        internal::reverse_mode_dependent_variable_activation(
          const_cast<ADNumberType &>(dependent_variables[0]));
        const std::size_t n_independent_variables =
          independent_variables.size();
        for (unsigned int i = 0; i < n_independent_variables; i++)
          {
            using derivative_type =
              typename ADNumberTraits<ADNumberType>::derivative_type;
            const derivative_type gradient_i =
              ADNumberTraits<ADNumberType>::get_directional_derivative(
                independent_variables[i], i);

            for (unsigned int j = 0; j <= i; ++j) // Symmetry
              {
                // Extract higher-order directional derivatives. Depending on
                // the AD number type, the result may be another AD number or a
                // floating point value.
                const ScalarType hessian_ij =
                  internal::NumberType<ScalarType>::value(
                    ADNumberTraits<derivative_type>::get_directional_derivative(
                      gradient_i, j));
                hessian[i][j] = hessian_ij;
                if (i != j)
                  hessian[j][i] = hessian_ij; // Symmetry
              }
          }
      }

      // === Vector drivers ===

      static void
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType> &             values)
      {
        Assert(values.size() == dependent_variables.size(),
               ExcDimensionMismatch(values.size(), dependent_variables.size()));

        const std::size_t n_dependent_variables = dependent_variables.size();
        for (unsigned int i = 0; i < n_dependent_variables; i++)
          values[i] = ADNumberTraits<ADNumberType>::get_scalar_value(
            dependent_variables[i]);
      }

      static void
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType> &         jacobian)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            1));
        Assert(jacobian.m() == dependent_variables.size(),
               ExcDimensionMismatch(jacobian.m(), dependent_variables.size()));
        Assert(jacobian.n() == independent_variables.size(),
               ExcDimensionMismatch(jacobian.n(),
                                    independent_variables.size()));

        const std::size_t n_independent_variables =
          independent_variables.size();
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
        using accumulation_type =
          typename ADNumberTraits<ADNumberType>::derivative_type;
        std::vector<accumulation_type> rad_accumulation(
          n_independent_variables,
          dealii::internal::NumberType<accumulation_type>::value(0.0));

        for (unsigned int i = 0; i < n_dependent_variables; i++)
          {
            internal::reverse_mode_dependent_variable_activation(
              const_cast<ADNumberType &>(dependent_variables[i]));
            for (unsigned int j = 0; j < n_independent_variables; j++)
              {
                const accumulation_type df_i_dx_j =
                  ADNumberTraits<ADNumberType>::get_directional_derivative(
                    independent_variables[j],
                    i /*This number doesn't really matter*/) -
                  rad_accumulation[j];
                jacobian[i][j] =
                  internal::NumberType<ScalarType>::value(df_i_dx_j);
                rad_accumulation[j] += df_i_dx_j;
              }
          }
      }

    private:
      /**
       * A flag that states whether or not dependent variables can be marked
       * within the current phase of operations.
       */
      bool dependent_variable_marking_safe;
    };


    /**
     * Specialization for auto-differentiable numbers that use
     * forward mode to compute the first (and, if supported, second)
     * derivatives.
     */
    template <typename ADNumberType, typename ScalarType>
    struct TapelessDrivers<
      ADNumberType,
      ScalarType,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::adolc_tapeless ||
                              ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::sacado_dfad ||
                              ADNumberTraits<ADNumberType>::type_code ==
                                NumberTypes::sacado_dfad_dfad>::type>
    {
      TapelessDrivers()
        : dependent_variable_marking_safe(false)
      {}

      // === Configuration ===

      static void
      initialize_global_environment(const unsigned int n_independent_variables)
      {
        internal::configure_tapeless_mode<ADNumberType>(
          n_independent_variables);
      }

      // === Operation status ===
      void
      allow_dependent_variable_marking()
      {
        dependent_variable_marking_safe = true;
      }

      void
      prevent_dependent_variable_marking()
      {
        dependent_variable_marking_safe = false;
      }

      bool
      is_dependent_variable_marking_allowed() const
      {
        return dependent_variable_marking_safe;
      }

      // === Scalar drivers ===

      static ScalarType
      value(const std::vector<ADNumberType> &dependent_variables)
      {
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(), 1));
        return ADNumberTraits<ADNumberType>::get_scalar_value(
          dependent_variables[0]);
      }

      static void
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType> &             gradient)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            1));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(), 1));
        Assert(gradient.size() == independent_variables.size(),
               ExcDimensionMismatch(gradient.size(),
                                    independent_variables.size()));

        // In forward mode, the gradients are computed from the
        // dependent variables
        const std::size_t n_independent_variables =
          independent_variables.size();
        for (unsigned int i = 0; i < n_independent_variables; i++)
          gradient[i] = internal::NumberType<ScalarType>::value(
            ADNumberTraits<ADNumberType>::get_directional_derivative(
              dependent_variables[0], i));
      }

      static void
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType> &         hessian)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            2));
        Assert(dependent_variables.size() == 1,
               ExcDimensionMismatch(dependent_variables.size(), 1));
        Assert(hessian.m() == independent_variables.size(),
               ExcDimensionMismatch(hessian.m(), independent_variables.size()));
        Assert(hessian.n() == independent_variables.size(),
               ExcDimensionMismatch(hessian.n(), independent_variables.size()));

        // In forward mode, the gradients are computed from the
        // dependent variables
        const std::size_t n_independent_variables =
          independent_variables.size();
        for (unsigned int i = 0; i < n_independent_variables; i++)
          {
            using derivative_type =
              typename ADNumberTraits<ADNumberType>::derivative_type;
            const derivative_type gradient_i =
              ADNumberTraits<ADNumberType>::get_directional_derivative(
                dependent_variables[0], i);

            for (unsigned int j = 0; j <= i; ++j) // Symmetry
              {
                // Extract higher-order directional derivatives. Depending on
                // the AD number type, the result may be another AD number or a
                // floating point value.
                const ScalarType hessian_ij =
                  internal::NumberType<ScalarType>::value(
                    ADNumberTraits<derivative_type>::get_directional_derivative(
                      gradient_i, j));
                hessian[i][j] = hessian_ij;
                if (i != j)
                  hessian[j][i] = hessian_ij; // Symmetry
              }
          }
      }

      // === Vector drivers ===

      static void
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType> &             values)
      {
        Assert(values.size() == dependent_variables.size(),
               ExcDimensionMismatch(values.size(), dependent_variables.size()));

        const std::size_t n_dependent_variables = dependent_variables.size();
        for (unsigned int i = 0; i < n_dependent_variables; i++)
          values[i] = ADNumberTraits<ADNumberType>::get_scalar_value(
            dependent_variables[i]);
      }

      static void
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType> &         jacobian)
      {
        Assert(
          AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1,
          ExcSupportedDerivativeLevels(
            AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
            1));
        Assert(jacobian.m() == dependent_variables.size(),
               ExcDimensionMismatch(jacobian.m(), dependent_variables.size()));
        Assert(jacobian.n() == independent_variables.size(),
               ExcDimensionMismatch(jacobian.n(),
                                    independent_variables.size()));

        const std::size_t n_independent_variables =
          independent_variables.size();
        const std::size_t n_dependent_variables = dependent_variables.size();

        // In forward mode, the gradients are computed from the
        // dependent variables
        for (unsigned int i = 0; i < n_dependent_variables; i++)
          for (unsigned int j = 0; j < n_independent_variables; j++)
            jacobian[i][j] = internal::NumberType<ScalarType>::value(
              ADNumberTraits<ADNumberType>::get_directional_derivative(
                dependent_variables[i], j));
      }

    private:
      /**
       * A flag that states whether or not dependent variables can be marked
       * within the current phase of operations.
       */
      bool dependent_variable_marking_safe;
    };

  } // namespace AD
} // namespace Differentiation


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE


#endif
