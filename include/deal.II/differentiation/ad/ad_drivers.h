// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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

#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/ad_number_types.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_ADOLC

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <adolc/internal/usrparms.h>
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
     * A prototype driver class for taped auto-differentiable numbers.
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
       * Return a list of registered tape indices.
       */
      std::vector<typename Types<ADNumberType>::tape_index>
      get_registered_tape_indices() const;

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
       * Return a flag that, when <code>true</code>, indicates that the retaping
       * of the dependent function for the chosen @p tape_index is necessary for
       * a reliable computation to be performed.
       * This may be necessary if a sign comparison within branched operations
       * yields different results to those computed at the original tape
       * evaluation point.
       *
       * This issue, known as "branch switching", can be clarified by means of
       * a trivial, contrived example:
       * @code
       * ADNumberType func (ADNumberType x, ADNumberType y, ADNumberType z)
       * {
       *   if (x < y)
       *     return y;
       *   else
       *     return x*z;
       * }
       * @endcode
       * During taping, the conditional statement may be either <tt>true</tt> or
       * <tt>false</tt>, and the result (with its sensitivities) returned by
       * this function.
       * The AD library doesn't just record the parse tree of the operations
       * applied on the branch chosen at the time to taping, but also checks
       * that the condition continues to be satisfied. For some other evaluation
       * of the tape (i.e. for some different inputs @p x and @p y), the other
       * branch of the conditional check may be chosen. The result of following
       * this code path has not been recorded on the tape, and therefore cannot
       * be evaluated. In such a case, the underlying AD library will be able to
       * tell you that it is necessary to re-record the tape at the new
       * evaluation point in order to resolve the new code branch. This function
       * can be used to find out whether this is so.
       *
       * @note The chosen tape index must be greater than
       * Numbers<ADNumberType>::invalid_tape_index and less than
       * Numbers<ADNumberType>::max_tape_index.
       */
      bool
      requires_retaping(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      /**
       * Return a flag that, when <code>true</code>, indicates that the retaping
       * of the dependent function is necessary for a reliable computation to be
       * performed on the currently active tape.
       * This may be necessary if a sign comparison within branched operations
       * yields different results to those computed at the original tape
       * evaluation point.
       *
       * This issue, known as "branch switching", can be clarified by means of
       * a trivial, contrived example:
       * @code
       * ADNumberType func (ADNumberType x, ADNumberType y, ADNumberType z)
       * {
       *   if (x < y)
       *     return y;
       *   else
       *     return x*z;
       * }
       * @endcode
       * During taping, the conditional statement may be either <tt>true</tt> or
       * <tt>false</tt>, and the result (with its sensitivities) returned by
       * this function.
       * The AD library doesn't just record the parse tree of the operations
       * applied on the branch chosen at the time to taping, but also checks
       * that the condition continues to be satisfied. For some other evaluation
       * of the tape (i.e. for some different inputs @p x and @p y), the other
       * branch of the conditional check may be chosen. The result of following
       * this code path has not been recorded on the tape, and therefore cannot
       * be evaluated. In such a case, the underlying AD library will be able to
       * tell you that it is necessary to re-record the tape at the new
       * evaluation point in order to resolve the new code branch. This function
       * can be used to find out whether this is so.
       */
      bool
      last_action_requires_retaping() const;

      /**
       * Completely erases the tape with the given @p tape_index.
       */
      void
      remove_tape(const typename Types<ADNumberType>::tape_index tape_index);

      /**
       * Reset the state of the class.
       *
       * @note This also resets the active tape number to an invalid number, and
       * deactivates the recording mode for taped variables.
       */
      void
      reset(const bool clear_registered_tapes);

      /**
       * Print the status of all queryable data. Exactly what is printed and
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
      void
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream &                                 stream) const;

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
      ScalarType
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<ScalarType> &independent_variables) const;

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
      void
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<ScalarType> &independent_variables,
               Vector<ScalarType> &           gradient) const;

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
      void
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<ScalarType> &independent_variables,
              FullMatrix<ScalarType> &       hessian) const;

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
      void
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int             n_dependent_variables,
             const std::vector<ScalarType> &independent_variables,
             Vector<ScalarType> &           values) const;

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
      void
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int             n_dependent_variables,
               const std::vector<ScalarType> &independent_variables,
               FullMatrix<ScalarType> &       jacobian) const;

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
      ScalarType
      value(const std::vector<ADNumberType> &dependent_variables) const;

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
      void
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType> &             gradient) const;

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
      void
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType> &         hessian) const;

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
      void
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType> &             values) const;

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
      void
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType> &         jacobian) const;

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
#    ifdef __clang__
      static const typename Types<ADNumberType>::tape_index max_tape_index =
        TBUFNUM;
#    else
      // For some reason, the test adolc/helper_tape_index_01 indicates that
      // ADOL-C does not reliably perform correct computations for the full
      // range of tape indices when GCC is the compiler. So we limit this number
      // according to the results of the test.
      static const typename Types<ADNumberType>::tape_index max_tape_index =
        TBUFNUM - 2;
#    endif
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

      /**
       * Constructor
       */
      TapedDrivers();


      /**
       * @name Taping
       */
      //@{

      bool
      is_recording() const;

      typename Types<ADNumberType>::tape_index
      active_tape_index() const;

      bool
      keep_independent_values() const;

      bool
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes in_obufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes in_lbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes in_vbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes in_tbufsize);

      void
      start_taping(const typename Types<ADNumberType>::tape_index tape_index,
                   const bool keep_independent_values);

      void
      stop_taping(
        const typename Types<ADNumberType>::tape_index active_tape_index,
        const bool                                     write_tapes_to_file);

      std::vector<typename Types<ADNumberType>::tape_index>
      get_registered_tape_indices() const;

      void
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index);

      bool
      requires_retaping(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      bool
      last_action_requires_retaping() const;

      void
      remove_tape(const typename Types<ADNumberType>::tape_index tape_index);

      void
      reset(const bool clear_registered_tapes);

      void
      print(std::ostream &stream) const;

      void
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream &                                 stream) const;

      //@}

      /**
       * @name Drivers for scalar functions (one dependent variable)
       */
      //@{

      scalar_type
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<scalar_type> &independent_variables) const;

      void
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<scalar_type> &independent_variables,
               Vector<scalar_type> &           gradient) const;

      void
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<scalar_type> &independent_variables,
              FullMatrix<scalar_type> &       hessian) const;

      //@}

      /**
       * @name Drivers for vector functions (multiple dependent variables)
       */
      //@{

      void
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int              n_dependent_variables,
             const std::vector<scalar_type> &independent_variables,
             Vector<scalar_type> &           values) const;

      void
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int              n_dependent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type> &       jacobian) const;

      //@}

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
       * The status of the last function or derivative evaluation performed on
       * the selected tape. As quoted from the ADOL-C manual, this can take
       * on one of six different values with the following interpretation:
       *
       * - <b>+3:</b> The function is locally analytic.
       * - <b>+2:</b> The function is locally analytic but the sparsity
       *     structure (compared to the situation at the taping point) may have
       *     changed, e.g. while at taping arguments <code>fmax(a,b)</code>
       *     returned <code>a</code>, we get <code>b</code> at the argument
       *     currently used.
       * - <b>+1:</b> At least one of the functions <code>fmin</code>,
       *     <code>fmax</code> or <code>fabs</code> is evaluated at a tie or
       *     zero, respectively. Hence, the function to be differentiated is
       *     Lipschitz-continuous but possibly non-differentiable.
       * - <b>0:</b>  Some arithmetic comparison involving adoubles yields a
       *     tie. Hence, the function to be differentiated may be discontinuous.
       * - <b>-1:</b> An <code>adouble</code> comparison yields different
       *     results from the evaluation point at which the tape was generated.
       * - <b>-2:</b> The argument of a user-defined quadrature has changed from
       *     the evaluation point at which the tape was generated.
       *
       * When the @p status variable takes a negative value, retaping of the dependent
       * function is necessary for a reliable computation to be performed.
       * This status can be queried through the requires_retaping() and
       * last_action_requires_retaping() functions.
       */
      mutable std::map<typename Types<ADNumberType>::tape_index, int> status;

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
    };

#  else

    /**
     * Specialization for taped ADOL-C auto-differentiable numbers.
     *
     * Although we could revert to the default definition for the
     * unspecialized TapedDrivers class, we add this specialization
     * to provide a more descriptive error message if any of its
     * member functions are called.
     */
    template <typename ADNumberType>
    struct TapedDrivers<
      ADNumberType,
      double,
      typename std::enable_if<ADNumberTraits<ADNumberType>::type_code ==
                              NumberTypes::adolc_taped>::type>
    {
      using scalar_type = double;

      /**
       * @name Taping
       */
      //@{

      bool
      is_recording() const;

      typename Types<ADNumberType>::tape_index
      active_tape_index() const;

      bool
      keep_independent_values() const;

      bool
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes);

      void
      start_taping(const typename Types<ADNumberType>::tape_index, const bool);

      void
      stop_taping(const typename Types<ADNumberType>::tape_index, const bool);

      std::vector<typename Types<ADNumberType>::tape_index>
      get_registered_tape_indices() const;

      void
      activate_tape(const typename Types<ADNumberType>::tape_index);

      bool
      requires_retaping(const typename Types<ADNumberType>::tape_index) const;

      bool
      last_action_requires_retaping() const;

      void
      remove_tape(const typename Types<ADNumberType>::tape_index);

      void
      reset(const bool);

      void
      print(std::ostream &stream) const;

      void
      print_tape_stats(const typename Types<ADNumberType>::tape_index,
                       std::ostream &) const;

      //@}

      /**
       * @name Drivers for scalar functions (one dependent variable)
       */
      //@{

      scalar_type
      value(const typename Types<ADNumberType>::tape_index,
            const std::vector<scalar_type> &) const;

      void
      gradient(const typename Types<ADNumberType>::tape_index,
               const std::vector<scalar_type> &,
               Vector<scalar_type> &) const;

      void
      hessian(const typename Types<ADNumberType>::tape_index,
              const std::vector<scalar_type> &,
              FullMatrix<scalar_type> &) const;

      //@}

      /**
       * @name Drivers for vector functions (multiple dependent variables)
       */
      //@{

      void
      values(const typename Types<ADNumberType>::tape_index,
             const unsigned int,
             const std::vector<scalar_type> &,
             Vector<scalar_type> &) const;

      void
      jacobian(const typename Types<ADNumberType>::tape_index,
               const unsigned int,
               const std::vector<scalar_type> &,
               FullMatrix<scalar_type> &) const;

      //@}
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

      /**
       * @name Taping
       */
      //@{

      bool
      is_recording() const;

      typename Types<ADNumberType>::tape_index
      active_tape_index() const;

      bool
      keep_independent_values() const;

      bool
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      void
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes obufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes lbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes vbufsize,
        const typename Types<ADNumberType>::tape_buffer_sizes tbufsize);

      void
      start_taping(const typename Types<ADNumberType>::tape_index tape_index,
                   const bool keep_independent_values);

      void
      stop_taping(
        const typename Types<ADNumberType>::tape_index active_tape_index,
        const bool                                     write_tapes_to_file);

      std::vector<typename Types<ADNumberType>::tape_index>
      get_registered_tape_indices() const;

      void
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index);

      bool
      requires_retaping(
        const typename Types<ADNumberType>::tape_index tape_index) const;

      bool
      last_action_requires_retaping() const;

      void
      remove_tape(const typename Types<ADNumberType>::tape_index tape_index);

      void
      reset(const bool clear_registered_tapes);

      void
      print(std::ostream &stream) const;

      void
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream &                                 stream) const;

      //@}

      /**
       * @name Drivers for scalar functions (one dependent variable)
       */
      //@{

      scalar_type
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<scalar_type> &independent_variables) const;

      void
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<scalar_type> &independent_variables,
               Vector<scalar_type> &           gradient) const;

      void
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<scalar_type> &independent_variables,
              FullMatrix<scalar_type> &       hessian) const;

      //@}

      /**
       * @name Drivers for vector functions (multiple dependent variables)
       */
      //@{

      void
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int              n_dependent_variables,
             const std::vector<scalar_type> &independent_variables,
             Vector<scalar_type> &           values) const;

      void
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int              n_dependent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type> &       jacobian) const;

      //@}

    private:
      /**
       * Copy a vector of floats into a vector of doubles
       */
      std::vector<double>
      vector_float_to_double(const std::vector<float> &in) const;

      /**
       * The object that actually takes care of the taping
       */
      TapedDrivers<ADNumberType, double> taped_driver;
    };


    // -------------   TapelessDrivers   -------------


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
      /**
       * Constructor
       */
      TapelessDrivers();

      /**
       * @name Configuration
       */
      //@{

      static void
      initialize_global_environment(const unsigned int n_independent_variables);

      //@}

      /**
       * Operation status
       */
      //@{

      void
      allow_dependent_variable_marking();

      void
      prevent_dependent_variable_marking();

      bool
      is_dependent_variable_marking_allowed() const;

      //@}

      /**
       * @name Drivers for scalar functions
       */
      //@{

      ScalarType
      value(const std::vector<ADNumberType> &dependent_variables) const;

      void
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType> &             gradient) const;

      void
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType> &         hessian) const;

      //@}

      /**
       * @name Drivers for vector functions
       */
      //@{

      void
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType> &             values) const;

      void
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType> &         jacobian) const;

      //@}

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
      /**
       * Constructor
       */
      TapelessDrivers();

      /**
       * @name Configuration
       */
      //@{

      static void
      initialize_global_environment(const unsigned int n_independent_variables);

      //@}

      /**
       * Operation status
       */
      //@{

      void
      allow_dependent_variable_marking();

      void
      prevent_dependent_variable_marking();

      bool
      is_dependent_variable_marking_allowed() const;

      //@}

      /**
       * @name Drivers for scalar functions
       */
      //@{

      ScalarType
      value(const std::vector<ADNumberType> &dependent_variables) const;

      void
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType> &             gradient) const;

      void
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType> &         hessian) const;

      //@}

      /**
       * @name Drivers for vector functions
       */
      //@{

      void
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType> &             values) const;

      void
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType> &         jacobian) const;

      //@}

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
