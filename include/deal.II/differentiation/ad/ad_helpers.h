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

#ifndef dealii_differentiation_ad_ad_helpers_h
#define dealii_differentiation_ad_ad_helpers_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)

#  include <deal.II/base/numbers.h>
#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/differentiation/ad/ad_drivers.h>
#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/adolc_number_types.h>
#  include <deal.II/differentiation/ad/adolc_product_types.h>
#  include <deal.II/differentiation/ad/sacado_number_types.h>
#  include <deal.II/differentiation/ad/sacado_product_types.h>

#  include <deal.II/fe/fe_values_extractors.h>

#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <algorithm>
#  include <iostream>
#  include <iterator>
#  include <numeric>
#  include <set>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace AD
  {
    /**
     * A base helper class that facilitates the evaluation of the derivative(s)
     * of a number of user-defined dependent variables $\mathbf{f}(\mathbf{X})$
     * with respect to a set of independent variables $\mathbf{X}$, that is
     * $\dfrac{d^{i} \mathbf{f}(\mathbf{X})}{d \mathbf{X}^{i}}$.
     *
     * This class is templated on the floating point type @p scalar_type of the
     * number that we'd like to differentiate, as well as an enumeration
     * indicating the @p ADNumberTypeCode . The @p ADNumberTypeCode dictates
     * which auto-differentiation library is to be used, and what the nature of
     * the underlying auto-differentiable number is. Refer to the
     * @ref auto_symb_diff module for more details in this regard.
     *
     * For all of the classes derived from this base class, there are two
     * possible ways that the code in which they are used can be structured.
     * The one approach is effectively a subset of the other, and which might
     * be necessary to use depends on the nature of the chosen
     * auto-differentiable number.
     *
     * When "tapeless" numbers are employed, the most simple code structure
     * would be of the following form:
     *
     * @code
     *   // Initialize AD helper
     *   ADHelperType<tapeless_AD_typecode> ad_helper (...);
     *
     *   // Register independent variables
     *   ad_helper.register_independent_variable(...);
     *
     *   // Extract the sensitive equivalent of the independent variables.
     *   // They are the auto-differentiable counterparts to the values
     *   // used as arguments to the register_independent_variable() function.
     *   // The operations conducted with these AD numbers will be tracked.
     *   const auto ad_independent_variables
     *     = ad_helper.get_sensitive_variables(...);
     *
     *   // Use the sensitive variables to compute the dependent variables.
     *   const auto ad_dependent_variables = func(sensitive_variables);
     *
     *   // Register the dependent variables with the helper class
     *   ad_helper.register_dependent_variables(ad_dependent_variables);
     *
     *   // Compute derivatives of the dependent variables
     *   const auto derivatives = ad_helper.compute_gradients();
     * @endcode
     *
     * Note that since the specialized classes interpret the independent
     * variables in different ways, above represents only an outline of the
     * steps taken to compute derivatives. More specific examples are outlined
     * in the individual classes that specialize this base class.
     *
     * When "taped" numbers are to be used, the above code should be wrapped by
     * a few more lines of code to manage the taping procedure:
     *
     * @code
     *   // Initialize AD helper
     *   ADHelperType<taped_or_tapeless_AD_typecode> ad_helper (...);
     *
     *   // An optional call to set the amount of memory to be allocated to
     *   // storing taped data
     *   ad_helper.set_tape_buffer_sizes();
     *
     *    // Select a tape number to record to
     *   const typename Types<ad_type>::tape_index  tape_index = ...;
     *
     *   // Indicate that we are about to start tracing the operations for
     *   // function evaluation on the tape. If this tape has already been used
     *   // (i.e. the operations are already recorded) then we (optionally)
     *   // load the tape and reuse this data.
     *   const bool is_recording
     *     = ad_helper.start_recording_operations(tape_index);
     *   if (is_recording == true)
     *   {
     *     // This is the "recording" phase of the operations.
     *     // In this block one places the majority of the operations described
     *     // in the previous code block. The set of operations that are
     *     // conducted here therefore includes the following steps:
     *     // - Register independent variables
     *     // - Extract the sensitive equivalent of the independent variables
     *     // - Use the sensitive variables to compute the dependent variables
     *     // - Register the dependent variables with the helper class
     *
     *     // Indicate that have completed tracing the operations onto the tape.
     *     ad_helper.stop_recording_operations(false); // write_tapes_to_file
     *   }
     *   else
     *   {
     *     // This is the "tape reuse" phase of the operations.
     *     // Here we will leverage the already traced operations that reside
     *     // on a tape, and simply revaluate the tape at a different point
     *     // to get the function values and their derivatives.
     *
     *     // Load the existing tape to be reused
     *     ad_helper.activate_recorded_tape(tape_no);
     *
     *     // Set the new values of the independent variables where the recorded
     *     // dependent functions are to be evaluated (and differentiated
     *     // around).
     *     ad_helper.set_independent_variable(...);
     *   }
     *
     *   // Compute derivatives of the dependent variables
     *   const auto derivatives = ad_helper.compute_gradients();
     * @endcode
     *
     * The second approach outlined here is more general than the first, and
     * will work equally well for both taped and tapeless auto-differentiable
     * numbers.
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @todo Make this class thread safe for Sacado number and ADOL-C tapeless
     * numbers (if supported).
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class HelperBase
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename AD::NumberTraits<ScalarType, ADNumberTypeCode>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename AD::NumberTraits<ScalarType, ADNumberTypeCode>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{f}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e., the dimension of the
       * image space.
       */
      HelperBase(const unsigned int n_independent_variables,
                 const unsigned int n_dependent_variables);

      /**
       * Destructor
       */
      virtual ~HelperBase() = default;

      //@}

      /**
       * @name Interrogation of internal information
       */
      //@{

      /**
       * Return the number of independent variables that this object expects to
       * work with. This is the dimension of the domain space.
       */
      std::size_t
      n_independent_variables() const;

      /**
       * Return the number of dependent variables that this object expects to
       * operate on. This is the dimension of the image space.
       */
      std::size_t
      n_dependent_variables() const;

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
       * Print the values currently assigned to the independent variables.
       *
       * @param[in] stream The output stream to which the values are to be
       * written.
       */
      void
      print_values(std::ostream &stream) const;

      /**
       * Print the statistics regarding the usage of the tapes.
       *
       * @param[in] tape_index The index of the tape to get the statistics of.
       * @param[out] stream The output stream to which the values are to be
       * written.
       *
       * @note This function only produces meaningful output when @p ad_type
       * is a taped auto-differentiable number.
       */
      void
      print_tape_stats(const typename Types<ad_type>::tape_index tape_index,
                       std::ostream &                            stream) const;

      //@}

      /**
       * @name Operations specific to tapeless mode
       */
      //@{

      /**
       * Pre-specify the number of @p independent_variables to be used in
       * tapeless mode.
       *
       * Although this function is called internally in the HelperBase
       * constructor, there may be occasions when ADOL-C tapeless numbers
       * (<tt>adtl::adoubles</tt>) are created before an instance of this class
       * is created. This function therefore allows one to declare at the
       * earliest possible instance how many directional derivatives will be
       * considered in tapeless mode.
       *
       * @warning With @p ensure_persistent_setting set to <tt>true</tt> when
       * the @p ad_type is an ADOL-C tapeless number, calling this function
       * leaves the set number of directional derivatives in a persistent state.
       * It will therefore not be possible to further modify the number of
       * directional derivatives to be tracked by <tt>adtl::adoubles</tt>'s
       * during course of the program's execution.
       */
      static void
      configure_tapeless_mode(const unsigned int n_independent_variables,
                              const bool ensure_persistent_setting = true);

      //@}

      /**
       * @name Operations specific to taped mode: Recording tapes
       */
      //@{

      /**
       * Reset the state of the helper class.
       *
       * When an instance of an HelperBase is stored as a class member object
       * with the intention to reuse its instance, it may be necessary to reset
       * the state of the object before use. This is because, internally, there
       * is error checking performed to ensure that the correct
       * auto-differentiable data is being tracked and used only when
       * appropriate. This function clears all member data and, therefore,
       * allows the state of all internal flags to be safely reset to their
       * initial state.
       *
       * In the rare case that the number of independent or dependent variables
       * has changed, this can also reconfigured by passing in the appropriate
       * arguments to the function.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{f}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e., the dimension of the
       * image space.
       * @param[in] clear_registered_tapes A flag that indicates the that
       * list of @p registered_tapes must be cleared.
       * If set to <tt>true</tt> then the data structure that tracks which
       * tapes have been recorded is cleared as well. It is then expected that
       * any preexisting tapes be re-recorded.
       *
       * @note This also resets the active tape number to an invalid number, and
       * deactivates the recording mode for taped variables.
       */
      virtual void
      reset(const unsigned int n_independent_variables =
              dealii::numbers::invalid_unsigned_int,
            const unsigned int n_dependent_variables =
              dealii::numbers::invalid_unsigned_int,
            const bool clear_registered_tapes = true);

      /**
       * Return whether or not this class is tracking calculations performed
       * with its marked independent variables.
       */
      bool
      is_recording() const;

      /**
       * Return the tape index which is currently activated for recording or
       * reading.
       */
      typename Types<ad_type>::tape_index
      active_tape_index() const;

      /**
       * Return whether or not a tape number has already been used
       * or registered.
       */
      bool
      is_registered_tape(
        const typename Types<ad_type>::tape_index tape_index) const;

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
        const typename Types<ad_type>::tape_buffer_sizes obufsize = 64 * 1024 *
                                                                    1024,
        const typename Types<ad_type>::tape_buffer_sizes lbufsize = 64 * 1024 *
                                                                    1024,
        const typename Types<ad_type>::tape_buffer_sizes vbufsize = 64 * 1024 *
                                                                    1024,
        const typename Types<ad_type>::tape_buffer_sizes tbufsize = 64 * 1024 *
                                                                    1024);

      /**
       * Enable recording mode for a given tape. The use of this function is
       * mandatory if the auto-differentiable number is a taped type. However,
       * for the purpose of developing generic code, it can also be safely
       * called for tapeless auto-differentiable numbers.
       *
       * The operations that take place between this function call and that
       * of stop_recording_operations() are recorded to the tape and can
       * be replayed and reevaluated as necessary.
       *
       * The typical set of operations to be performed during this "recording"
       * phase (between the calls to start_recording_operations() and
       * stop_recording_operations() ) are:
       *   - Definition of some independent variables via
       *     register_independent_variable() or
       *     register_independent_variables(). These define the branch of
       *     operations tracked by the tape. If the @p keep flag is set to
       *     <tt>true</tt> then these represent precisely the point about which
       *     the function derivatives are to be computed. If the @p keep flag is
       *     set to <tt>false</tt> then these only represent dummy values, and
       *     the point at which the function derivatives are to be computed must
       *     be set by calling set_independent_variables() again.
       *   - Extraction of a set of independent variables of auto-differentiable
       *     type using get_sensitive_variables(). These are then tracked during
       *     later computations.
       *   - Defining the dependent variables via register_dependent_variable()
       *     or register_dependent_variables(). These are the functions that
       *     will be differentiated with respect to the independent variables.
       *
       * @param[in] tape_index The index of the tape to be written
       * @param[in] overwrite_tape Express whether tapes are allowed to be
       * overwritten. If <tt>true</tt> then any existing tape with a given
       * @p tape_index will be destroyed and a new tape traced over it.
       * @param[in] keep_independent_values Determines whether the numerical
       * values of all independent variables are recorded in the tape buffer.
       * If true, then the tape can be immediately used to perform computations
       * after recording is complete.
       *
       * @note During the recording phase, no value(), gradient(), hessian(),
       * or jacobian() operations can be performed.
       *
       * @note The chosen tape index must be greater than
       * Numbers<ad_type>::invalid_tape_index and less than
       * Numbers<ad_type>::max_tape_index.
       */
      bool
      start_recording_operations(
        const typename Types<ad_type>::tape_index tape_index,
        const bool                                overwrite_tape = false,
        const bool keep_independent_values                       = true);

      /**
       * Disable recording mode for a given tape. The use of this function is
       * mandatory if the auto-differentiable number is a taped type. However,
       * for the purpose of developing generic code, it can also be safely
       * called for tapeless auto-differentiable numbers.
       *
       * @note After this function call, the tape is considered ready for use and
       * operations such as value(), gradient() or hessian() can be executed.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      stop_recording_operations(const bool write_tapes_to_file = false);

      /**
       * Select a pre-recorded tape to read from.
       *
       * @param[in] tape_index The index of the tape to be read from.
       *
       * @note The chosen tape index must be greater than
       * Numbers<ad_type>::invalid_tape_index and less than
       * Numbers<ad_type>::max_tape_index.
       */
      void
      activate_recorded_tape(
        const typename Types<ad_type>::tape_index tape_index);

      /**
       * Return a flag that, when <code>true</code>, indicates that the
       * retaping of the dependent function is necessary for a reliable
       * computation to be performed on a tape with the given @p tape_index.
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
       * For the output of this function to be meaningful, it must be called
       * after activate_recorded_tape() is called and the new evaluation
       * point for the tape (i.e. values of the independent variables) have
       * been set and subsequently used (i.e. in the determination of the values
       * or derivatives of the dependent variables).
       */
      bool
      recorded_tape_requires_retaping(
        const typename Types<ad_type>::tape_index tape_index) const;

      /**
       * Return a flag that, when <code>true</code>, indicates that the
       * retaping of the dependent function is necessary for a reliable
       * computation to be performed on the currently active tape.
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
       * For the output of this function to be meaningful, it must be called
       * after activate_recorded_tape() is called and the new evaluation
       * point for the tape (i.e. values of the independent variables) have
       * been set and subsequently used (i.e. in the determination of the values
       * or derivatives of the dependent variables).
       */
      bool
      active_tape_requires_retaping() const;

      /**
       * Clears and removes the currently active tape.
       *
       * This is typically only necessary when branch switching is detected on
       * the original tape at evaluation point. This state can be checked using
       * the active_tape_requires_retaping() function.
       */
      void
      clear_active_tape();

      //@}

    protected:
      /**
       * @name Drivers and taping
       */
      //@{

      /**
       * An object used to help manage stored tapes.
       *
       * In the event that the @p ad_type is a tapeless AD type, then the object
       * constructed here is, effectively, a dummy one.
       */
      TapedDrivers<ad_type, scalar_type> taped_driver;

      /**
       * An object used to help manage tapeless data structures.
       *
       * In the event that the @p ad_type is a taped AD type, then the object
       * constructed here is, effectively, a dummy one.
       */
      TapelessDrivers<ad_type, scalar_type> tapeless_driver;

      /**
       * Select a tape to record to or read from.
       *
       * This function activates a tape, but depending on whether @p read_mode
       * is set, the tape is either taken as previously written to (and put
       * into read-only mode), or cleared for (re-)taping.
       *
       * @param[in] tape_index The index of the tape to be written to/read
       *            from.
       * @param[in] read_mode A flag that marks whether or not we expect to
       *            read data from a preexisting tape.
       *
       * @note The chosen tape index must be greater than
       * Numbers<ad_type>::invalid_tape_index and less than
       * Numbers<ad_type>::max_tape_index.
       */
      void
      activate_tape(const typename Types<ad_type>::tape_index tape_index,
                    const bool                                read_mode);

      //@}

      /**
       * @name Independent variables
       */
      //@{

      /**
       * A set of independent variables $\mathbf{X}$ that differentiation will
       * be performed with respect to.
       *
       * The gradients and Hessians of dependent variables will be computed
       * at these finite values.
       */
      mutable std::vector<scalar_type> independent_variable_values;

      /**
       * A set of sensitive independent variables $\mathbf{X}$ that
       * differentiation will be performed with respect to.
       *
       * The gradients and Hessians of dependent variables will be computed
       * using these configured AD numbers. Note that only reverse-mode AD
       * requires that the sensitive independent variables be stored.
       */
      mutable std::vector<ad_type> independent_variables;

      /**
       * A list of registered independent variables that have been manipulated
       * for a given set of operations.
       */
      std::vector<bool> registered_independent_variable_values;

      /**
       * A list of registered independent variables that have been extracted and
       * their sensitivities marked.
       */
      mutable std::vector<bool> registered_marked_independent_variables;

      /**
       * Reset the boolean vector @p registered_independent_variable_values that
       * indicates which independent variables we've been manipulating for
       * the current set of operations.
       */
      void
      reset_registered_independent_variables();

      /**
       * Set the actual value of the independent variable $X_{i}$.
       *
       * @param[in] index The index in the vector of independent variables.
       * @param[in] value The value to set the index'd independent variable to.
       */
      void
      set_sensitivity_value(const unsigned int index, const scalar_type &value);

      /**
       * Initialize an independent variable $X_{i}$ such that subsequent
       * operations performed with it are tracked.
       *
       * @note Care must be taken to mark each independent variable only once.
       *
       * @note The order in which the independent variables are marked defines the
       * order of all future internal operations. They must be manipulated in
       * the same order as that in which they are first marked. If not then,
       * for example, ADOL-C won't throw an error, but rather it might complain
       * nonsensically during later computations or produce garbage results.
       *
       * @param[in] index The index in the vector of independent variables.
       * @param[out] out An auto-differentiable number that is ready for use in
       * computations. The operations that are performed with it are recorded on
       * the tape and will be considered in the computation of dependent
       * variable values.
       */
      void
      mark_independent_variable(const unsigned int index, ad_type &out) const;

      /**
       * Finalize the state of the independent variables before use.
       *
       * This step and the storage of the independent variables is done
       * separately because some derived classes may offer the capability
       * to add independent variables in a staggered manner. This function
       * is to be triggered when these values are considered finalized
       * and we can safely initialize the sensitive equivalents of those
       * values.
       */
      void
      finalize_sensitive_independent_variables() const;

      /**
       * Initialize an independent variable $X_{i}$.
       *
       * @param[out] out An auto-differentiable number that is ready for use in
       * standard computations. The operations that are performed with it are
       * not recorded on the tape, and so should only be used when not in
       * recording mode.
       * @param[in] index The index in the vector of independent variables.
       */
      void
      initialize_non_sensitive_independent_variable(const unsigned int index,
                                                    ad_type &out) const;

      /**
       * The number of independent variables that have been manipulated within a
       * set of operations.
       */
      unsigned int
      n_registered_independent_variables() const;

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * The set of dependent variables $\mathbf{f}(\mathbf{X})$ of which the
       * derivatives with respect to $\mathbf{X}$ will be computed.
       *
       * The gradients and Hessians of these dependent variables will be
       * computed at the values $\mathbf{X}$ set with the
       * set_sensitivity_value() function.
       *
       * @note These are stored as an @p ad_type so that we can use them to
       * compute function values and directional derivatives in the case that
       * tapeless numbers are used
       */
      std::vector<ad_type> dependent_variables;

      /**
       * A list of registered dependent variables.
       */
      std::vector<bool> registered_marked_dependent_variables;

      /**
       * Reset the boolean vector @p registered_marked_dependent_variables that
       * indicates which independent variables have been manipulated by the
       * current set of operations. All entries in the vector are set to the
       * value of the @p flag.
       */
      void
      reset_registered_dependent_variables(const bool flag = false);

      /**
       * The number of dependent variables that have been registered.
       */
      unsigned int
      n_registered_dependent_variables() const;

      /**
       * Register the definition of the index'th dependent variable
       * $f(\mathbf{X})$.
       *
       * @param[in] index The index of the entry in the global list of dependent
       * variables that this function belongs to.
       * @param[in] func The recorded function that defines a dependent
       * variable.
       *
       * @note Each dependent variable must only be registered once.
       */
      void
      register_dependent_variable(const unsigned int index,
                                  const ad_type &    func);

      //@}

    }; // class HelperBase



    /**
     * A general helper class that facilitates the evaluation of a vector of
     * functions, as well as its first derivatives (their Jacobian).
     * This class would typically be used to compute the linearization of a
     * set of local nonlinear equations, but can also be used as the basis of
     * the linearization of the residual vector defined on the level of a finite
     * element (for example, in order to compute the Jacobian matrix necessary
     * in Newton-type solvers for nonlinear problems).
     *
     * @note When using the cell-level taped AD methods in 3d and/or with higher
     * order elements, it is incredibly easy to exceed the tape buffer size.
     * The reason for this is two-fold:
     *   1. there are many independent variables (the local
     *      degrees-of-freedom) to take the derivatives with respect to, and
     *   2. the expressions for the dependent variables (each being a component
     *      of the residual vector) in terms of all of the independent variables
     *      are lengthy, especially when non-trivial constitutive laws are
     *      considered.
     * These buffer variables dictate the amount of memory allocated to a tape
     * before it is written to file (at a significant performance loss).
     * Therefore for ADOL-C taped AD numbers, it may be desirable to
     * create a file ".adolcrc" in the program run directory and set the buffer
     * size therein (as is suggested by the ADOL-C manual). For example, the
     * following settings increase the default buffer size by 128 times:
     * @code
     * "OBUFSIZE" "67108864"
     * "LBUFSIZE" "67108864"
     * "VBUFSIZE" "67108864"
     * "TBUFSIZE" "67108864"
     * @endcode
     * Note that the quotation marks are mandatory.
     * An alternative approach that allows for run-time decision making is to
     * use the HelperBase::set_tape_buffer_sizes() function before starting
     * taping (as done via the HelperBase::start_recording_operations()
     * function).
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class CellLevelBase : public HelperBase<ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{f}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e., the dimension of the
       * image space.
       */
      CellLevelBase(const unsigned int n_independent_variables,
                    const unsigned int n_dependent_variables);

      /**
       * Destructor
       */
      virtual ~CellLevelBase() = default;

      //@}

      /**
       * @name Independent variables
       */
      //@{

      /**
       * Register the complete set of independent variables $\mathbf{X}$ that
       * represent the local degree of freedom values.
       *
       * @param[in] dof_values A vector field associated with local
       * degree of freedom values on the current finite element. These define
       * the values of all independent variables. When considering taped AD
       * numbers with branching functions, to avoid potential issues with branch
       * switching it may be a good idea to choose these values close or equal
       * to those that will be later evaluated and linearized around.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      register_dof_values(const std::vector<scalar_type> &dof_values);

      /**
       * Register the complete set of independent variables $\mathbf{X}$ that
       * represent the local degree of freedom values.
       *
       * @param[in] values A global field from which the values of all
       * independent variables will be extracted. This typically will be the
       * solution vector around which point a residual vector is to be
       * computed and around which linearization is to occur.
       * When considering taped AD numbers with branching functions, to avoid
       * potential issues with branch switching it may be a good idea to choose
       * these values close or equal to those that will be later evaluated and
       * linearized around.
       * @param[in] local_dof_indices A vector of degree of freedom indices from
       * which to extract the local degree of freedom values. This would
       * typically obtained by calling <code>cell->get_dof_indices()</code>.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      template <typename VectorType>
      void
      register_dof_values(
        const VectorType &                                  values,
        const std::vector<dealii::types::global_dof_index> &local_dof_indices);

      /**
       * Return the complete set of degree of freedom values as represented by
       * auto-differentiable numbers. These are the independent
       * variables $\mathbf{X}$ about which the solution is linearized.
       *
       * This function indicates to the AD library that implements the
       * auto-differentiable number type that operations performed on these
       * numbers are to be tracked so they are considered "sensitive"
       * variables. This is, therefore, the set of variables with which one
       * would then perform computations, and based on which one can then
       * extract both the value of the function and its derivatives with the
       * member functions below. The values of the components of the returned
       * object are initialized to the values set with
       * register_independent_variable().
       *
       * @return An array of auto-differentiable type numbers representing the
       * local degree of freedom values.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      const std::vector<ad_type> &
      get_sensitive_dof_values() const;

      //@}

      /**
       * @name Operations specific to taped mode: Reusing tapes
       */
      //@{

      /**
       * Set the values for the independent variables $\mathbf{X}$, i.e., the
       * linearization point.
       *
       * @param[in] dof_values A vector field associated with local
       * degree of freedom values on the current finite element. These define
       * the values of all independent variables.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note If the @p keep_independent_values flag has been set when
       * HelperBase::start_recording_operations() is called then the tape is
       * immediately usable after creation, and the values of the independent
       * variables set by register_dof_values() are those at which the function
       * is to be evaluated. In this case, a separate call to this function is
       * not strictly necessary.
       */
      void
      set_dof_values(const std::vector<scalar_type> &dof_values);

      /**
       * Set the values for the independent variables $\mathbf{X}$, i.e., the
       * linearization point.
       *
       * @param[in] values A vector field from which the values of all
       * independent variables is to be extracted.
       * @param[in] local_dof_indices A vector of degree of freedom indices from
       * which to extract the local degree of freedom values. This would
       * typically obtained by calling <code>cell->get_dof_indices()</code>.
       *
       * @note If the @p keep_independent_values flag has been set when
       * HelperBase::start_recording_operations() is called then the tape is
       * immediately usable after creation, and the values of the independent
       * variables set by register_dof_values() are those at which the function
       * is to be evaluated. In this case, a separate call to this function is
       * not strictly necessary.
       */
      template <typename VectorType>
      void
      set_dof_values(
        const VectorType &                                  values,
        const std::vector<dealii::types::global_dof_index> &local_dof_indices);

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Compute the value of the residual vector field
       * $\mathbf{r}(\mathbf{X})$.
       *
       * @param[out] residual A Vector object with the value for each component
       * of the vector field evaluated at the point defined by the independent
       * variable values.
       *
       * @note The size of the @p residual vector is determined by the derived
       * classes, as it depends on the order of the dependent variable(s)
       * derivative(s) that it represents. Code examples that show how to use
       * this interface will be provided in the documentation of the derived
       * classes.
       */
      virtual void
      compute_residual(Vector<scalar_type> &residual) const = 0;

      /**
       * Compute the gradient (first derivative) of the residual vector field
       * with respect to all independent variables, i.e.
       * @f[
       *   \frac{\partial\mathbf{r}(\mathbf{X})}{\partial\mathbf{X}}
       * @f]
       *
       * @param[out] linearization A FullMatrix with the gradient of each
       * component of the vector field evaluated at the point defined by the
       * independent variable values.
       *
       * @note The dimensions of the @p linearization matrix is determined by
       * the derived classes, as it depends on the order of the dependent
       * variable(s) derivative(s) that it represents. Code examples that show
       * how to use this interface will be provided in the documentation of
       * the derived classes.
       */
      virtual void
      compute_linearization(FullMatrix<scalar_type> &linearization) const = 0;

      //@}

    }; // class CellLevelBase



    /**
     * A helper class that facilitates the implementation of a generic
     * (incremental) variational formulation from which the computation of the
     * residual vector, as well as its linearization, is automated. This class
     * would typically be used to derive the residual vector and tangent matrix
     * (defined on the level of a cell), or a linearized system of
     * equations, starting from a scalar energy functional.
     *
     * An example of its usage in the case of a residual and tangent
     * computations might be as follows (in this case we'll compute the
     * linearization of a finite-strain hyperelastic solid from a stored/strain
     * energy density function):
     *
     * @code
     *   // Existing data structures:
     *   Vector<double> solution (...); // Or another vector type
     *   std::vector<types::global_dof_index> local_dof_indices (...);
     *   const FEValuesExtractors::Vector u_fe (...);
     *   FEValues<dim> fe_values (...);
     *   const unsigned int n_q_points (...);
     *   FullMatrix<double> cell_matrix (...);
     *   Vector<double> cell_rhs (...);
     *
     *   // Assembly loop:
     *   for (auto &cell : ...)
     *   {
     *     cell->get_dof_indices(local_dof_indices);
     *     const unsigned int n_independent_variables =
     *       local_dof_indices.size();
     *
     *     // Create some aliases for the AD helper.
     *     // In the example, the AD_typecode used for the template argument can
     *     // be refer to either a taped or tapeless type.
     *     using ADHelper = AD::EnergyFunctional<...>;
     *     using ADNumberType = typename ADHelper::ad_type;
     *
     *     // Create and initialize an instance of the helper class.
     *     ADHelper ad_helper(n_independent_variables);
     *
     *     // Initialize the local data structures for assembly.
     *     // This is also taken care of by the ADHelper, so this step could
     *     // be skipped.
     *     cell_rhs.reinit(n_independent_variables);
     *     cell_matrix.reinit(n_independent_variables,n_independent_variables);
     *
     *     // An optional call to set the amount of memory to be allocated to
     *     // storing taped data.
     *     // If using a taped AD number then we would likely want to increase
     *     // the buffer size from the default values as the expression for each
     *     // residual component will likely be lengthy, and therefore memory
     *     // intensive.
     *     ad_helper.set_tape_buffer_sizes(...);
     *
     *     // If using a taped AD number, then at this point we would initiate
     *     // taping of the expression for the energy for this FE type and
     *     // material combination:
     *
     *     // Select a tape number to record to
     *     const typename Types<ad_type>::tape_index tape_index = ...;
     *
     *     // Indicate that we are about to start tracing the operations for
     *     // function evaluation on the tape. If this tape has already been
     *     // used (i.e. the operations are already recorded) then we
     *     // (optionally) load the tape and reuse this data.
     *     const bool is_recording
     *       = ad_helper.start_recording_operations(tape_index);
     *
     *     // The steps that follow in the recording phase are required for
     *     // tapeless methods as well.
     *     if (is_recording == true)
     *     {
     *       // This is the "recording" phase of the operations.
     *       // First, we set the values for all DoFs.
     *       ad_helper.register_dof_values(solution, local_dof_indices);
     *
     *       // Then we get the complete set of degree of freedom values as
     *       // represented by auto-differentiable numbers. The operations
     *       // performed with these variables are tracked by the AD library
     *       // from this point until stop_recording_operations() is called.
     *       const std::vector<ADNumberType> dof_values_ad
     *         = ad_helper.get_sensitive_dof_values();
     *
     *       // Then we do some problem specific tasks, the first being to
     *       // compute all values, gradients, etc. based on sensitive AD DoF
     *       // values. Here we are fetching the displacement gradients at each
     *       // quadrature point.
     *       std::vector<Tensor<2, dim, ADNumberType>> Grad_u(
     *         n_q_points, Tensor<2, dim, ADNumberType>());
     *       fe_values[u_fe].get_function_gradients_from_local_dof_values(
     *         dof_values_ad, Grad_u);
     *
     *       // This variable stores the cell total energy.
     *       // IMPORTANT: Note that it is hand-initialized with a value of
     *       // zero. This is a highly recommended practise, as some AD numbers
     *       // appear not to safely initialize their internal data structures.
     *       ADNumberType energy_ad = ADNumberType(0.0);
     *
     *       // Compute the cell total energy = (internal + external) energies
     *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
     *       {
     *         // Calculate the deformation gradient at this quadrature point
     *         const Tensor<2, dim, ADNumberType> F =
     *           unit_symmetric_tensor<dim>() + Grad_u[q_point];
     *         Assert(numbers::value_is_greater_than(determinant(F), 0.0),
     *                ExcMessage("Negative determinant of the deformation "
     *                           "gradient detected!"));
     *
     *         // Add contribution of the internal energy:
     *         // Integrate the stored energy density function with the current
     *         // solution.
     *         energy_ad += get_Psi(F) * fe_values.JxW(q_point);
     *       }
     *
     *       // Add contribution from external energy:
     *       // Loop over faces and accumulate external energy into cell
     *       // total energy.
     *       for (unsigned int face : ...)
     *         if (cell->face(face)->at_boundary())
     *           energy_ad += ...
     *
     *       // Register the definition of the total cell energy
     *       ad_helper.register_energy_functional(energy_ad);
     *
     *       // Indicate that we have completed tracing the operations onto
     *       // the tape.
     *       ad_helper.stop_recording_operations(false); // write_tapes_to_file
     *     }
     *     else
     *     {
     *       // This is the "tape reuse" phase of the operations.
     *       // Here we will leverage the already traced operations that reside
     *       // on a tape, and simply re-evaluate the tape at a different point
     *       // to get the function values and their derivatives.
     *
     *       // Load the existing tape to be reused
     *       ad_helper.activate_recorded_tape(tape_index);
     *
     *       // Set the new values of the independent variables where the
     *       // recorded dependent functions are to be evaluated (and
     *       // differentiated around).
     *       ad_helper.set_dof_values(solution, local_dof_indices);
     *     }
     *
     *     // Compute the residual values and their Jacobian at the
     *     // evaluation point
     *     ad_helper.compute_residual(cell_rhs);
     *     cell_rhs *= -1.0; // RHS = - residual
     *     ad_helper.compute_linearization(cell_matrix);
     *   }
     * @endcode
     *
     * In most use cases, and in particular in the code example shown above,
     * the number of independent variables equals the number of
     * <code>dofs_per_cell</code> for the used finite element.
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class EnergyFunctional : public CellLevelBase<ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\Psi(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       *
       * @note There is only one dependent variable associated with the total
       * energy attributed to the local finite element. That is to say, this
       * class assumes that the (local) right hand side and matrix contribution
       * is computed from the first and second derivatives of a scalar
       * function $\Psi(\mathbf{X})$.
       */
      EnergyFunctional(const unsigned int n_independent_variables);

      /**
       * Destructor
       */
      virtual ~EnergyFunctional() = default;

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Register the definition of the total cell energy
       * $\Psi(\mathbf{X})$.
       *
       * @param[in] energy A recorded function that defines the total cell
       * energy. This represents the single dependent variable from which both
       * the residual and its linearization are to be computed.
       *
       * @note For this class that expects only a single scalar dependent
       * variable, this function must only be called once per tape.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      register_energy_functional(const ad_type &energy);

      /**
       * Evaluation of the total scalar energy functional for a chosen set of
       * degree of freedom values, i.e.
       * @f[
       *   \Psi(\mathbf{X}) \vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are obtained by calling
       * CellLevelBase::set_dof_values().
       *
       * @return The value of the energy functional at the evaluation point
       * corresponding to a chosen set of local degree of freedom values.
       */
      scalar_type
      compute_energy() const;

      /**
       * Evaluation of the residual for a chosen set of degree of freedom
       * values. Underlying this is the computation of the gradient (first
       * derivative) of the scalar function $\Psi$ with respect to all
       * independent variables, i.e.
       * @f[
       *   \mathbf{r}(\mathbf{X}) =
       * \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{X}}
       * \Big\vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are obtained by calling
       * CellLevelBase::set_dof_values().
       *
       * @param[out] residual A Vector object, for which the value for each
       * entry represents the residual value for the corresponding local
       * degree of freedom. The output @p residual vector has a length
       * corresponding to @p n_independent_variables.
       */
      void
      compute_residual(Vector<scalar_type> &residual) const override;

      /**
       * Compute the linearization of the residual vector around a chosen set
       * of degree of freedom values. Underlying this is the computation of the
       * Hessian (second derivative) of the scalar function $\Psi$ with respect
       * to all independent variables, i.e.
       * @f[
       *   \frac{\partial\mathbf{r}(\mathbf{X})}{\partial\mathbf{X}}
       *     =
       * \frac{\partial^{2}\Psi(\mathbf{X})}{\partial\mathbf{X}
       * \otimes \partial\mathbf{X}} \Big\vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are obtained by calling
       * CellLevelBase::set_dof_values().
       *
       * @param[out] linearization A FullMatrix representing the linearization
       * of the residual vector. The output @p linearization matrix has
       * dimensions corresponding to
       * <code>n_independent_variables</code>$\times$<code>n_independent_variables</code>.
       */
      virtual void
      compute_linearization(
        FullMatrix<scalar_type> &linearization) const override;

      //@}

    }; // class EnergyFunctional



    /**
     * A helper class that facilitates the evaluation and automated
     * linearization of a vector of functions that represents a residual vector
     * (as computed from some corresponding local degree of freedom values).
     * This class would typically be used to compute the linearization of a
     * residual vector defined on the level of a cell, or for local
     * nonlinear equations.
     *
     * An example of its usage in the case of a residual linearization
     * might be as follows (in this case we'll compute the
     * linearization of a finite-strain magneto-elastic solid from the residual,
     * as constructed from the Piola-Kirchoff stress and magnetic induction
     * assuming the magnetic scalar potential formulation):
     *
     * @code
     *   // Existing data structures:
     *   Vector<double> solution (...); // Or another vector type
     *   std::vector<types::global_dof_index> local_dof_indices (...);
     *   const FEValuesExtractors::Vector u_fe (...);
     *   const FEValuesExtractors::Scalar msp_fe (...);
     *   const unsigned int u_block (...);
     *   const unsigned int msp_block (...);
     *   FEValues<dim> fe_values (...);
     *   const unsigned int n_q_points (...);
     *   FullMatrix<double> cell_matrix (...);
     *   Vector<double> cell_rhs (...);
     *
     *   // Assembly loop:
     *   for (auto &cell : ...)
     *   {
     *     cell->get_dof_indices(local_dof_indices);
     *     const unsigned int n_independent_variables
     *       = local_dof_indices.size();
     *     const unsigned int n_dependent_variables
     *       = local_dof_indices.size();
     *
     *     // Create some aliases for the AD helper.
     *     // In this example, we strictly assume that we're using tapeless
     *     // AD types, and so the AD_typecode used in the template argument
     *     // must refer to one of these types. See the example for the
     *     // EnergyFunctional for details on how to extend
     *     // support to taped AD numbers.
     *     using ADHelper = AD::ResidualLinearization<...>;
     *     using ADNumberType = typename ADHelper::ad_type;
     *
     *     // Create and initialize an instance of the helper class.
     *     ADHelper ad_helper(n_independent_variables,n_dependent_variables);
     *
     *     // Initialize the local data structures for assembly.
     *     // This is also taken care of by the ADHelper, so this step could
     *     // be skipped.
     *     cell_rhs.reinit(n_dependent_variables);
     *     cell_matrix.reinit(n_independent_variables,n_dependent_variables);
     *
     *     // This next code block corresponds to the "recording" phase, where
     *     // the operations performed with the AD numbers are tracked and
     *     // differentiation is performed.
     *     {
     *       // First, we set the values for all DoFs.
     *       ad_helper.register_dof_values(solution, local_dof_indices);
     *
     *       // Then we get the complete set of degree of freedom values as
     *       // represented by auto-differentiable numbers. The operations
     *       // performed with these variables are tracked by the AD library
     *       // from this point until stop_recording_operations() is called.
     *       const std::vector<ADNumberType> dof_values_ad
     *         = ad_helper.get_sensitive_dof_values();
     *
     *       // Then we do some problem specific tasks, the first being to
     *       // compute all values, gradients, etc. based on sensitive AD DoF
     *       // values. Here we are fetching the displacement gradients at each
     *       // quadrature point, as well as the gradients of the magnetic
     *       // scalar potential field.
     *       std::vector<Tensor<2, dim, ADNumberType>> Grad_u(
     *         n_q_points, Tensor<2, dim, ADNumberType>());
     *       std::vector<Tensor<1, dim, ADNumberType>> Grad_msp(
     *         n_q_points, Tensor<1, dim, ADNumberType>());
     *       fe_values[u_fe].get_function_gradients_from_local_dof_values(
     *         dof_values_ad, Grad_u);
     *       fe_values[msp_fe].get_function_gradients_from_local_dof_values(
     *         dof_values_ad, Grad_msp);
     *
     *       // This variable stores the cell residual vector contributions.
     *       // IMPORTANT: Note that each entry is hand-initialized with a value
     *       // of zero. This is a highly recommended practise, as some AD
     *       // numbers appear not to safely initialize their internal data
     *       // structures.
     *       std::vector<ADNumberType> residual_ad (
     *         n_dependent_variables, ADNumberType(0.0));
     *
     *       // Compute the cell total residual
     *       //   = (internal + external) contributions
     *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
     *       {
     *         // Calculate the deformation gradient and magnetic field at this
     *         // quadrature point
     *         const Tensor<2, dim, ADNumberType> F =
     *           unit_symmetric_tensor<dim>() + Grad_u[q_point];
     *         const Tensor<1, dim, ADNumberType> H = -Grad_msp[q_point];
     *         Assert(numbers::value_is_greater_than(determinant(F), 0.0),
     *                ExcMessage("Negative determinant of the deformation "
     *                           "gradient detected!"));
     *
     *         // Extract some configuration dependent variables from our
     *         // nonlinear constitutive law for the current quadrature point.
     *         // In this way they only have to be computed once per quadrature
     *         // point.
     *         const SymmetricTensor<2,dim,ad_type> S = get_S(F,H);
     *         const Tensor<1,dim,ad_type>          B = get_B(F,H);
     *
     *         // Define some position-dependent aliases, to make the assembly
     *         // process easier to follow.
     *         const double JxW = fe_values.JxW(q_point);
     *
     *         // Add contribution of the internal forces:
     *         // Note that we assemble the residual contribution directly
     *         // as it is this vector that is to be automatically linearized.
     *         for (unsigned int I = 0; I < n_dofs_per_cell; ++I)
     *         {
     *           // Determine the component and block associated with
     *           // the I'th local degree of freedom.
     *           const unsigned int block_I =
     *             fe.system_to_base_index(I).first.first;
     *
     *           if (block_I == u_block) // u-terms
     *           {
     *             // Variation of the Green-Lagrange strain tensor
     *             // associated with the I'th vector-valued basis function.
     *             const SymmetricTensor<2,dim,double> dE_I
     *               = symmetrize(transpose(F)
     *               * fe_values[u_fe].gradient(I,q_point));
     *
     *             residual_ad[I] += (dE_I*S) * JxW;
     *           }
     *           else if (block_I == msp_block)
     *           {
     *             // Variation of the magnetic field vector associated with
     *             // the I'th scalar-valued basis function
     *             const Tensor<1,dim,double> dH_I
     *               = -fe_values[msp_fe].gradient(I, q_point);
     *
     *             residual_ad[I] -= (dH_I*B) * JxW;
     *           }
     *         }
     *       }
     *
     *       // Add contribution from external sources. If these contributions
     *       // are also solution dependent then they will also be consistently
     *       // linearized.
     *       // Loop over faces and accumulate external contributions into the
     *       // cell total residual.
     *       for (unsigned int face : ...)
     *         if (cell->face(face)->at_boundary())
     *           residual_ad[I] += ...
     *
     *       // Register the definition of the cell residual
     *       ad_helper.register_residual_vector(residual_ad);
     *     }
     *
     *     // Compute the residual values and their Jacobian at the
     *     // evaluation point
     *     ad_helper.compute_residual(cell_rhs);
     *     cell_rhs *= -1.0; // RHS = - residual
     *     ad_helper.compute_linearization(cell_matrix);
     *   }
     * @endcode
     *
     * In most use cases, and in particular in the code example shown above,
     * both the number of independent and dependent variables equals the
     * number of <code>dofs_per_cell</code> for the used finite element.
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class ResidualLinearization
      : public CellLevelBase<ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{r}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{r}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{r}$, i.e., the dimension of the
       * image space.
       */
      ResidualLinearization(const unsigned int n_independent_variables,
                            const unsigned int n_dependent_variables);

      /**
       * Destructor
       */
      virtual ~ResidualLinearization() = default;

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Register the definition of the cell residual vector
       * $\mathbf{r}(\mathbf{X})$.
       *
       * @param[in] residual A vector of recorded functions that defines the
       * residual. The components of this vector represents the dependent
       * variables.
       *
       * @note For this class that expects only vector fields of dependent
       * variables, this function must only be called once per tape.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      register_residual_vector(const std::vector<ad_type> &residual);

      /**
       * Evaluation of the residual for a chosen set of degree of freedom
       * values. This corresponds to the computation of the residual vector,
       * i.e.
       * @f[
       *   \mathbf{r}(\mathbf{X}) \vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are obtained by calling
       * CellLevelBase::set_dof_values().
       *
       * @param[out] residual A Vector object, for which the value for each
       * entry represents the residual value for the corresponding local
       * degree of freedom. The output @p residual vector has a length
       * corresponding to @p n_dependent_variables.
       */
      virtual void
      compute_residual(Vector<scalar_type> &residual) const override;

      /**
       * Compute the linearization of the residual vector around a chosen set
       * of degree of freedom values. Underlying this is the computation of the
       * gradient (first derivative) of the residual vector $\mathbf{r}$ with
       * respect to all independent variables, i.e.
       * @f[
       *   \frac{\partial\mathbf{r}(\mathbf{X})}{\partial\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are obtained by calling
       * CellLevelBase::set_dof_values().
       *
       * @param[out] linearization A FullMatrix representing the linearization
       * of the residual vector. The output @p linearization matrix has
       * dimensions corresponding to
       * <code>n_dependent_variables</code>$\times$<code>n_independent_variables</code>.
       */
      virtual void
      compute_linearization(
        FullMatrix<scalar_type> &linearization) const override;

      //@}

    }; // class ResidualLinearization



    namespace internal
    {
      /**
       * A helper struct that assists with the extraction of data associated
       * with fields that are defined by FEExtractors.
       */
      template <int dim, typename ExtractorType>
      struct Extractor;


      /**
       * A helper struct that assists with the extraction of data associated
       * with fields that are defined by FEExtractors.
       * This particular specialization is for scalar fields.
       */
      template <int dim>
      struct Extractor<dim, FEValuesExtractors::Scalar>
      {
        /**
         * The number of components of the field.
         */
        static const unsigned int n_components = 1;

        /**
         * The tensor rank of the field.
         */
        static const unsigned int rank = 0;

        /**
         * The tensor type associated with this field.
         */
        template <typename NumberType>
        using tensor_type = Tensor<rank, dim, NumberType>;

        static_assert(
          n_components == tensor_type<double>::n_independent_components,
          "The number of components doesn't match that of the corresponding tensor type.");
        static_assert(
          rank == tensor_type<double>::rank,
          "The rank doesn't match that of the corresponding tensor type.");

        /**
         * The value type associated with this field.
         */
        // Note: FEValuesViews::Scalar::tensor_type is double, so we can't
        // use it (FEValuesViews) in this context.
        // In fact, sadly, all FEValuesViews objects expect doubles as value
        // types.
        template <typename NumberType>
        using value_type = NumberType;

        /**
         * The gradient type associated with this field.
         */
        template <typename NumberType>
        using gradient_type = Tensor<rank + 1, dim, NumberType>; // NumberType;

        /**
         * Return the first global component of this field.
         */
        static inline unsigned int
        first_component(const FEValuesExtractors::Scalar &extractor)
        {
          return extractor.component;
        }

        /**
         * Return a flag that indicates if the input @p unrolled_index
         * corresponds to a symmetric component of the field.
         *
         * For a scalar field, the single component is defined
         * to not have a symmetric counterpart.
         */
        static bool
        symmetric_component(const unsigned int unrolled_index)
        {
          (void)unrolled_index;
          return false;
        }

        /**
         * Return the local unrolled component corresponding to
         * @p column_offset entry of the @p table_indices.
         *
         * For a scalar field, the local component is always
         * equal to zero.
         */
        template <typename IndexType = unsigned int, int rank_in>
        static IndexType
        local_component(const TableIndices<rank_in> &table_indices,
                        const unsigned int           column_offset)
        {
          Assert(column_offset <= rank_in, ExcInternalError());
          (void)table_indices;
          (void)column_offset;
          return 0;
        }
      };


      /**
       * A helper struct that assists with the extraction of data associated
       * with fields that are defined by FEExtractors.
       * This particular specialization is for vector fields.
       */
      template <int dim>
      struct Extractor<dim, FEValuesExtractors::Vector>
      {
        /**
         * The number of components of the field.
         */
        static const unsigned int n_components = dim;

        /**
         * The tensor rank of the field.
         */
        static const unsigned int rank = 1;

        /**
         * The tensor type associated with this field.
         */
        template <typename NumberType>
        using tensor_type = Tensor<rank, dim, NumberType>;

        static_assert(
          n_components == tensor_type<double>::n_independent_components,
          "The number of components doesn't match that of the corresponding tensor type.");
        static_assert(
          rank == tensor_type<double>::rank,
          "The rank doesn't match that of the corresponding tensor type.");

        /**
         * The value type associated with this field.
         */
        template <typename NumberType>
        using value_type = tensor_type<NumberType>;

        /**
         * The gradient type associated with this field.
         */
        template <typename NumberType>
        using gradient_type = Tensor<rank + 1, dim, NumberType>;

        /**
         * Return the first global component of this field.
         */
        static inline unsigned int
        first_component(const FEValuesExtractors::Vector &extractor)
        {
          return extractor.first_vector_component;
        }

        /**
         *
         * Return a flag that indicates if the input @p unrolled_index
         * corresponds to a symmetric component of the field.
         *
         * For a vector field, the none of the vector components
         * have a symmetric counterpart.
         */
        static bool
        symmetric_component(const unsigned int unrolled_index)
        {
          (void)unrolled_index;
          return false;
        }

        /**
         * Return the table index corresponding to
         * @p column_offset entry of the input @p table_indices.
         */
        template <int rank_in>
        static TableIndices<rank>
        table_index_view(const TableIndices<rank_in> &table_indices,
                         const unsigned int           column_offset)
        {
          Assert(0 + column_offset < rank_in, ExcInternalError());
          return TableIndices<rank>(table_indices[column_offset]);
        }

        /**
         * Return the local unrolled component corresponding to
         * @p column_offset entry of the @p table_indices.
         *
         * This function computes and returns a local component
         * associated with the extractor's @p tensor_type from a
         * set of @p table_indices that are generally associated
         * with a tensor of equal or greater order. In effect, it
         * creates a view of a selected number of indices of the
         * input table, and interprets that subtable's indices as
         * the local index to be returned. Since the @p table_indices
         * may be of size greater than the extractor's @p rank,
         * the @p column_offset specifies the first index of the
         * input table to create the view from.
         */
        template <typename IndexType = unsigned int, int rank_in>
        static IndexType
        local_component(const TableIndices<rank_in> &table_indices,
                        const unsigned int           column_offset)
        {
          static_assert(
            rank_in >= rank,
            "Cannot extract more table indices than the input table has!");
          using TensorType = tensor_type<double>;
          return TensorType::component_to_unrolled_index(
            table_index_view(table_indices, column_offset));
        }
      };


      /**
       * A helper struct that assists with the extraction of data associated
       * with fields that are defined by FEExtractors.
       * This particular specialization is for rank-1 tensor fields.
       */
      template <int dim>
      struct Extractor<dim, FEValuesExtractors::Tensor<1>>
      {
        /**
         * The number of components of the field.
         */
        static const unsigned int n_components =
          Tensor<1, dim>::n_independent_components;

        /**
         * The tensor rank of the field.
         */
        static const unsigned int rank = 1;

        /**
         * The tensor type associated with this field.
         */
        template <typename NumberType>
        using tensor_type = Tensor<rank, dim, NumberType>;

        /**
         * The value type associated with this field.
         */
        template <typename NumberType>
        using value_type = tensor_type<NumberType>;

        /**
         * The gradient type associated with this field.
         */
        template <typename NumberType>
        using gradient_type = Tensor<rank + 1, dim, NumberType>;

        /**
         * Return the first global component of this field.
         */
        static inline unsigned int
        first_component(const FEValuesExtractors::Tensor<1> &extractor)
        {
          return extractor.first_tensor_component;
        }

        /**
         * Return a flag that indicates if the input @p unrolled_index
         * corresponds to a symmetric component of the field.
         *
         * For a vector field, the none of the vector components
         * have a symmetric counterpart.
         */
        static bool
        symmetric_component(const unsigned int unrolled_index)
        {
          (void)unrolled_index;
          return false;
        }

        /**
         * Return the table index corresponding to
         * @p column_offset entry of the input @p table_indices.
         */
        template <int rank_in>
        static TableIndices<rank>
        table_index_view(const TableIndices<rank_in> &table_indices,
                         const unsigned int           column_offset)
        {
          Assert(column_offset < rank_in, ExcInternalError());
          return TableIndices<rank>(table_indices[column_offset]);
        }

        /**
         * Return the local unrolled component corresponding to
         * a subset of table indices from the input @p table_indices.
         *
         * This function computes and returns a local component
         * associated with the extractor's @p tensor_type from a
         * set of @p table_indices that are generally associated
         * with a tensor of equal or greater order. In effect, it
         * creates a view of a selected number of indices of the
         * input table, and interprets that subtable's indices as
         * the local index to be returned. Since the @p table_indices
         * may be of size greater than the extractor's @p rank,
         * the @p column_offset specifies the first index of the
         * input table to create the view from.
         */
        template <typename IndexType = unsigned int, int rank_in>
        static IndexType
        local_component(const TableIndices<rank_in> &table_indices,
                        const unsigned int           column_offset)
        {
          static_assert(
            rank_in >= rank,
            "Cannot extract more table indices than the input table has!");
          using TensorType = tensor_type<double>;
          return TensorType::component_to_unrolled_index(
            table_index_view(table_indices, column_offset));
        }
      };


      /**
       * A helper struct that assists with the extraction of data associated
       * with fields that are defined by FEExtractors.
       * This particular specialization is for rank-2 tensor fields.
       */
      template <int dim>
      struct Extractor<dim, FEValuesExtractors::Tensor<2>>
      {
        /**
         * The number of components of the field.
         */
        static const unsigned int n_components =
          Tensor<2, dim>::n_independent_components;

        /**
         * The tensor rank of the field.
         */
        static const unsigned int rank = Tensor<2, dim>::rank;

        /**
         * The tensor type associated with this field.
         */
        template <typename NumberType>
        using tensor_type = Tensor<rank, dim, NumberType>;

        /**
         * The value type associated with this field.
         */
        template <typename NumberType>
        using value_type = tensor_type<NumberType>;

        /**
         * The gradient type associated with this field.
         */
        template <typename NumberType>
        using gradient_type = Tensor<rank + 1, dim, NumberType>;

        /**
         * Return the first global component of this field.
         */
        static inline unsigned int
        first_component(const FEValuesExtractors::Tensor<2> &extractor)
        {
          return extractor.first_tensor_component;
        }

        /**
         * Return a flag that indicates if the input @p unrolled_index
         * corresponds to a symmetric component of the field.
         *
         * For a rank-2 tensor field, the none of the tensor
         * components have a symmetric counterpart.
         */
        static bool
        symmetric_component(const unsigned int unrolled_index)
        {
          (void)unrolled_index;
          return false;
        }

        /**
         * Return the table indices corresponding to
         * @p column_offset entry of the input @p table_indices.
         */
        template <int rank_in>
        static TableIndices<rank>
        table_index_view(const TableIndices<rank_in> &table_indices,
                         const unsigned int           column_offset)
        {
          Assert(column_offset < rank_in, ExcInternalError());
          Assert(column_offset + 1 < rank_in, ExcInternalError());
          return TableIndices<rank>(table_indices[column_offset],
                                    table_indices[column_offset + 1]);
        }

        /**
         * Return the local unrolled component corresponding to
         * @p column_offset entry of the @p table_indices.
         *
         * This function computes and returns a local component
         * associated with the extractor's @p tensor_type from a
         * set of @p table_indices that are generally associated
         * with a tensor of equal or greater order. In effect, it
         * creates a view of a selected number of indices of the
         * input table, and interprets that subtable's indices as
         * the local index to be returned. Since the @p table_indices
         * may be of size greater than the extractor's @p rank,
         * the @p column_offset specifies the first index of the
         * input table to create the view from.
         */
        template <typename IndexType = unsigned int, int rank_in>
        static IndexType
        local_component(const TableIndices<rank_in> &table_indices,
                        const unsigned int           column_offset)
        {
          static_assert(
            rank_in >= rank,
            "Cannot extract more table indices than the input table has!");
          using TensorType = tensor_type<double>;
          return TensorType::component_to_unrolled_index(
            table_index_view(table_indices, column_offset));
        }
      };


      /**
       * A helper struct that assists with the extraction of data associated
       * with fields that are defined by FEExtractors.
       * This particular specialization is for rank-2 symmetric tensor fields.
       */
      template <int dim>
      struct Extractor<dim, FEValuesExtractors::SymmetricTensor<2>>
      {
        /**
         * The number of components of the field.
         */
        static const unsigned int n_components =
          SymmetricTensor<2, dim>::n_independent_components;

        /**
         * The tensor rank of the field.
         */
        static const unsigned int rank = SymmetricTensor<2, dim>::rank;

        /**
         * The tensor type associated with this field.
         */
        template <typename NumberType>
        using tensor_type = SymmetricTensor<rank, dim, NumberType>;

        /**
         * The value type associated with this field.
         */
        template <typename NumberType>
        using value_type = tensor_type<NumberType>;

        /**
         * The gradient type associated with this field.
         */
        template <typename NumberType>
        using gradient_type = Tensor<rank + 1, dim, NumberType>;

        /**
         * Return the first global component of this field.
         */
        static inline unsigned int
        first_component(const FEValuesExtractors::SymmetricTensor<2> &extractor)
        {
          return extractor.first_tensor_component;
        }

        /**
         * Return a flag that indicates if the input @p unrolled_index
         * corresponds to a symmetric component of the field.
         *
         * For a rank-2 symmetric tensor field, each of the
         * off-diagonal components have a symmetric counterpart,
         * while the diagonal components do not.
         */
        static bool
        symmetric_component(const unsigned int unrolled_index)
        {
          const TableIndices<2> table_indices =
            tensor_type<double>::unrolled_to_component_indices(unrolled_index);
          return table_indices[0] != table_indices[1];
        }

        /**
         * Return the table indices corresponding to
         * @p column_offset entry of the input @p table_indices.
         */
        template <int rank_in>
        static TableIndices<rank>
        table_index_view(const TableIndices<rank_in> &table_indices,
                         const unsigned int           column_offset)
        {
          Assert(column_offset < rank_in, ExcInternalError());
          Assert(column_offset + 1 < rank_in, ExcInternalError());
          return TableIndices<rank>(table_indices[column_offset],
                                    table_indices[column_offset + 1]);
        }

        /**
         * Return the local unrolled component corresponding to
         * @p column_offset entry of the @p table_indices.
         *
         * This function computes and returns a local component
         * associated with the extractor's @p tensor_type from a
         * set of @p table_indices that are generally associated
         * with a tensor of equal or greater order. In effect, it
         * creates a view of a selected number of indices of the
         * input table, and interprets that subtable's indices as
         * the local index to be returned. Since the @p table_indices
         * may be of size greater than the extractor's @p rank,
         * the @p column_offset specifies the first index of the
         * input table to create the view from.
         */
        template <typename IndexType = unsigned int, int rank_in>
        static IndexType
        local_component(const TableIndices<rank_in> &table_indices,
                        const unsigned int           column_offset)
        {
          static_assert(
            rank_in >= rank,
            "Cannot extract more table indices than the input table has!");
          using TensorType = tensor_type<double>;
          return TensorType::component_to_unrolled_index(
            table_index_view(table_indices, column_offset));
        }
      };


      /**
       * A helper struct that defines the return type of gradient (first
       * derivative) calculations of scalar fields with respect to a field
       * defined by the @p ExtractorType template parameter.
       */
      template <int dim, typename NumberType, typename ExtractorType>
      struct ScalarFieldGradient
      {
        /**
         * The type associated with computing the gradient of a scalar
         * field with respect to the given @p ExtractorType.
         */
        using type =
          typename Extractor<dim,
                             ExtractorType>::template tensor_type<NumberType>;
      };


      /**
       * An intermediate helper struct that defines the return type of Hessian
       * (second derivative) calculations of scalar fields with respect to
       * fields defined by the two extractor-type template parameters.
       * The first, @p ExtractorType_Row, defines the field that the first
       * derivatives are taken with respect to while the second,
       * @p ExtractorType_Col, defines the field that the second derivatives
       * are taken with respect to.
       */
      template <typename ExtractorType_Row, typename ExtractorType_Col>
      struct HessianType
      {
        /**
         * The type associated with computing the gradient of a scalar
         * field with respect to the given @p ExtractorType_Row
         * followed by the @p ExtractorType_Col.
         *
         * @note We set the return type for
         * HessianType<FEExtractor::Vector,FEExtractor::Vector>
         * as a normal Tensor. This is because if one has two vector components,
         * the coupling tensor (i.e. Hessian component<FE::V_1,FE::V_2>) is in
         * general not symmetric.
         */
        template <int rank, int dim, typename NumberType>
        using type = Tensor<rank, dim, NumberType>;
      };


      /**
       * An intermediate helper struct that defines the return type of Hessian
       * (second derivative) calculations of scalar fields with respect to
       * fields defined by the two extractor-type template parameters. This
       * particular specialization is for taking the first derivative with
       * respect to a symmetric tensor field, and the second with respect to a
       * scalar field.
       */
      template <>
      struct HessianType<FEValuesExtractors::SymmetricTensor<2>,
                         FEValuesExtractors::Scalar>
      {
        /**
         * The type associated with computing the gradient of a scalar
         * field with respect to the given
         * <code>ExtractorType_Row =
         * FEValuesExtractors::SymmetricTensor<2></code> followed by the
         * <code>ExtractorType_Col = FEValuesExtractors::Scalar</code>.
         */
        template <int /*rank*/, int dim, typename NumberType>
        using type = SymmetricTensor<2 /*rank*/, dim, NumberType>;
      };


      /**
       * An intermediate helper struct that defines the return type of Hessian
       * (second derivative) calculations of scalar fields with respect to
       * fields defined by the two extractor-type template parameters. This
       * particular specialization is for taking the first derivative with
       * respect to a scalar field, and the second with respect to a symmetric
       * tensor field.
       */
      template <>
      struct HessianType<FEValuesExtractors::Scalar,
                         FEValuesExtractors::SymmetricTensor<2>>
      {
        /**
         * The type associated with computing the gradient of a scalar
         * field with respect to the given
         * <code>ExtractorType_Row = FEValuesExtractors::Scalar</code>
         * followed by the
         * <code>ExtractorType_Col =
         * FEValuesExtractors::SymmetricTensor<2></code>.
         */
        template <int /*rank*/, int dim, typename NumberType>
        using type = SymmetricTensor<2 /*rank*/, dim, NumberType>;
      };


      /**
       * An intermediate helper struct that defines the return type of Hessian
       * (second derivative) calculations of scalar fields with respect to
       * fields defined by the two extractor-type template parameters. This
       * particular specialization is for taking both the first and second
       * derivatives with respect to symmetric tensor fields.
       */
      template <>
      struct HessianType<FEValuesExtractors::SymmetricTensor<2>,
                         FEValuesExtractors::SymmetricTensor<2>>
      {
        /**
         * The type associated with computing the gradient of a scalar
         * field with respect to the given
         * <code>ExtractorType_Row =
         * FEValuesExtractors::SymmetricTensor<2></code> followed by the
         * <code>ExtractorType_Col =
         * FEValuesExtractors::SymmetricTensor<2></code>.
         */
        template <int /*rank*/, int dim, typename NumberType>
        using type = SymmetricTensor<4 /*rank*/, dim, NumberType>;
      };


      /**
       * A helper struct that defines the final return type of Hessian
       * (second derivative) calculations of scalar fields with respect to
       * fields defined by the two extractor-type
       * template parameters. The first, @p ExtractorType_Row, defines the field
       * that the first derivatives are taken with respect to while the second,
       * @p ExtractorType_Col, defines the field that the second derivatives
       * are taken with respect to.
       */
      template <int dim,
                typename NumberType,
                typename ExtractorType_Row,
                typename ExtractorType_Col>
      struct ScalarFieldHessian
      {
        /**
         * The tensor rank of the resulting derivative computation.
         */
        static const int rank = Extractor<dim, ExtractorType_Row>::rank +
                                Extractor<dim, ExtractorType_Col>::rank;

        /**
         * The type associated with computing the Hessian of a scalar
         * field with first respect to the field defined by the
         * @p ExtractorType_Row and then with respect to the field defined by
         * the @p ExtractorType_Col.
         */
        using type =
          typename HessianType<ExtractorType_Row, ExtractorType_Col>::
            template type<rank, dim, NumberType>;
      };


      /**
       * A helper struct that defines the return type of value computations
       * of vector fields the @p ExtractorType_Field template parameter.
       */
      template <int dim, typename NumberType, typename ExtractorType_Field>
      using VectorFieldValue =
        ScalarFieldGradient<dim, NumberType, ExtractorType_Field>;


      /**
       * A helper struct that defines the final return type of Jacobian
       * (first derivative) calculations of vector fields with respect to
       * fields as defined by the two extractor-type template parameters.
       * The first, @p ExtractorType_Field, defines the field from which
       * the initial field values are computed while the second,
       * @p ExtractorType_Derivative, defines the field that the derivatives
       * are taken with respect to.
       */
      template <int dim,
                typename NumberType,
                typename ExtractorType_Field,
                typename ExtractorType_Derivative>
      using VectorFieldJacobian = ScalarFieldHessian<dim,
                                                     NumberType,
                                                     ExtractorType_Field,
                                                     ExtractorType_Derivative>;


      /**
       * Return a global view of the field component indices that correspond to
       * the input @p extractor. For this general function the
       * @p ignore_symmetries flag has no effect.
       */
      template <int dim,
                typename IndexType = unsigned int,
                typename ExtractorType>
      std::vector<IndexType>
      extract_field_component_indices(const ExtractorType &extractor,
                                      const bool ignore_symmetries = true)
      {
        (void)ignore_symmetries;
        const IndexType n_components =
          internal::Extractor<dim, ExtractorType>::n_components;
        const IndexType comp_first =
          internal::Extractor<dim, ExtractorType>::first_component(extractor);
        std::vector<IndexType> indices(n_components);
        std::iota(indices.begin(), indices.end(), comp_first);
        return indices;
      }


      /**
       * Return a global view of the field component indices that correspond to
       * the input FEValuesExtractors::SymmetricTensor @p extractor_symm_tensor.
       * If the @p ignore_symmetries is set <code>true</code>, then all
       * component of the tensor are considered to be independent. If set to
       * <code>false</code>, then the set of returned indices will contain
       * duplicate entries for components that are symmetric.
       */
      template <int dim, typename IndexType = unsigned int>
      std::vector<IndexType>
      extract_field_component_indices(
        const FEValuesExtractors::SymmetricTensor<2> &extractor_symm_tensor,
        const bool                                    ignore_symmetries = true)
      {
        using ExtractorType = FEValuesExtractors::SymmetricTensor<2>;
        const IndexType n_components =
          internal::Extractor<dim, ExtractorType>::n_components;
        if (ignore_symmetries == true)
          {
            const IndexType comp_first =
              internal::Extractor<dim, ExtractorType>::first_component(
                extractor_symm_tensor);
            std::vector<IndexType> indices(n_components);
            std::iota(indices.begin(), indices.end(), comp_first);
            return indices;
          }
        else
          {
            // First get all of the indices of the non-symmetric tensor
            const FEValuesExtractors::Tensor<2> extractor_tensor(
              extractor_symm_tensor.first_tensor_component);
            std::vector<IndexType> indices =
              extract_field_component_indices<dim>(extractor_tensor, true);

            // Then we overwrite any illegal entries with the equivalent indices
            // from the symmetric tensor
            for (unsigned int i = 0; i < indices.size(); ++i)
              {
                if (indices[i] >= n_components)
                  {
                    const TableIndices<2> ti_tensor =
                      Tensor<2, dim>::unrolled_to_component_indices(indices[i]);
                    const IndexType sti_new_index =
                      SymmetricTensor<2, dim>::component_to_unrolled_index(
                        ti_tensor);
                    indices[i] = sti_new_index;
                  }
              }

            return indices;
          }
      }


      /**
       * Set the unrolled component given by @p index in the generic tensor
       * @p t to the given @p value.
       */
      template <typename TensorType, typename NumberType>
      inline void
      set_tensor_entry(TensorType &       t,
                       const unsigned int unrolled_index,
                       const NumberType & value)
      {
        // Where possible, set values using TableIndices
        AssertIndexRange(unrolled_index, t.n_independent_components);
        t[TensorType::unrolled_to_component_indices(unrolled_index)] = value;
      }


      /**
       * Set the unrolled component given by @p index in the rank-0 tensor
       * @p t to the given @p value.
       */
      template <int dim, typename NumberType>
      inline void set_tensor_entry(Tensor<0, dim, NumberType> &t,
                                   const unsigned int          unrolled_index,
                                   const NumberType &          value)
      {
        AssertIndexRange(unrolled_index, 1);
        (void)unrolled_index;
        t = value;
      }


      /**
       * Set the value of @p t to the given @p value.
       * This function exists to provide compatibility with similar functions
       * that exist for use with the tensor classes.
       */
      template <typename NumberType>
      inline void
      set_tensor_entry(NumberType &       t,
                       const unsigned int unrolled_index,
                       const NumberType & value)
      {
        AssertIndexRange(unrolled_index, 1);
        (void)unrolled_index;
        t = value;
      }


      /**
       * Set the unrolled component given by the @p index_row and
       * the @p index_col in the fourth-order symmetric tensor
       * @p t to the given @p value.
       */
      template <int dim, typename NumberType>
      inline void set_tensor_entry(SymmetricTensor<4, dim, NumberType> &t,
                                   const unsigned int unrolled_index_row,
                                   const unsigned int unrolled_index_col,
                                   const NumberType & value)
      {
        // Fourth order symmetric tensors require a specialized interface
        // to extract values.
        using SubTensorType = SymmetricTensor<2, dim, NumberType>;
        AssertIndexRange(unrolled_index_row,
                         SubTensorType::n_independent_components);
        AssertIndexRange(unrolled_index_col,
                         SubTensorType::n_independent_components);
        const TableIndices<2> indices_row =
          SubTensorType::unrolled_to_component_indices(unrolled_index_row);
        const TableIndices<2> indices_col =
          SubTensorType::unrolled_to_component_indices(unrolled_index_col);
        t[indices_row[0]][indices_row[1]][indices_col[0]][indices_col[1]] =
          value;
      }


      /**
       * Return the value of the @p index'th unrolled component of the
       * generic tensor @p t.
       */
      template <int rank,
                int dim,
                typename NumberType,
                template <int, int, typename> class TensorType>
      inline NumberType
      get_tensor_entry(const TensorType<rank, dim, NumberType> &t,
                       const unsigned int                       unrolled_index)
      {
        // Where possible, get values using TableIndices
        AssertIndexRange(unrolled_index, t.n_independent_components);
        return t[TensorType<rank, dim, NumberType>::
                   unrolled_to_component_indices(unrolled_index)];
      }


      /**
       * Return the value of the @p index'th unrolled component of the
       * rank-0 tensor @p t.
       */
      template <int dim,
                typename NumberType,
                template <int, int, typename> class TensorType>
      inline NumberType
      get_tensor_entry(const TensorType<0, dim, NumberType> &t,
                       const unsigned int                    unrolled_index)
      {
        AssertIndexRange(unrolled_index, 1);
        (void)unrolled_index;
        return t;
      }


      /**
       * Return the value of @p t.
       * This function exists to provide compatibility with similar functions
       * that exist for use with the tensor classes.
       */
      template <typename NumberType>
      inline const NumberType &
      get_tensor_entry(const NumberType &t, const unsigned int unrolled_index)
      {
        AssertIndexRange(unrolled_index, 1);
        (void)unrolled_index;
        return t;
      }


      /**
       * Return a reference to the entry stored in the @p index'th unrolled
       * component of the generic tensor @p t.
       */
      template <int rank,
                int dim,
                typename NumberType,
                template <int, int, typename> class TensorType>
      inline NumberType &
      get_tensor_entry(TensorType<rank, dim, NumberType> &t,
                       const unsigned int                 unrolled_index)
      {
        // Where possible, get values using TableIndices
        AssertIndexRange(unrolled_index, t.n_independent_components);
        return t[TensorType<rank, dim, NumberType>::
                   unrolled_to_component_indices(unrolled_index)];
      }


      /**
       * Return a reference to the entry stored in the @p index'th unrolled
       * component of the rank-0 tensor @p t.
       */
      template <int dim,
                typename NumberType,
                template <int, int, typename> class TensorType>
      NumberType &get_tensor_entry(TensorType<0, dim, NumberType> &t,
                                   const unsigned int              index)
      {
        AssertIndexRange(index, 1);
        (void)index;
        return t;
      }


      /**
       * Return a reference to  @p t.
       * This function exists to provide compatibility with similar functions
       * that exist for use with the tensor classes.
       */
      template <typename NumberType>
      inline NumberType &
      get_tensor_entry(NumberType &t, const unsigned int index)
      {
        AssertIndexRange(index, 1);
        (void)index;
        return t;
      }

    } // namespace internal



    /**
     * A base helper class that facilitates the evaluation of point-wise defined
     * functions. This is the point-wise counterpart of the
     * CellLevelBase class, and was conceived for computations at a
     * continuum point, or quadrature point, rather than for finite-element
     * level calculations. That being said, the interface to this and the
     * derived classes are sufficiently generic that the dependent function(s)
     * and their argument(s), that are the independent variables, can be
     * interpreted in any manner that the user may choose.
     *
     * As it offers a field-based interface, this class would
     * typically be used to compute the derivatives of a constitutive law
     * defined at a quadrature point; however, it may also be used in other
     * contexts, such as to compute the linearization of a set of local
     * nonlinear equations.
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class PointLevelFunctionsBase
      : public HelperBase<ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the dimension of the associated input and output
       * tensor types.
       */
      static const unsigned int dimension = dim;

      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{f}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e., the dimension of the
       * image space.
       */
      PointLevelFunctionsBase(const unsigned int n_independent_variables,
                              const unsigned int n_dependent_variables);

      /**
       * Destructor
       */
      virtual ~PointLevelFunctionsBase() = default;

      //@}

      /**
       * @name Independent variables
       */
      //@{

      /**
       * @copydoc HelperBase::reset()
       */
      virtual void
      reset(const unsigned int n_independent_variables =
              dealii::numbers::invalid_unsigned_int,
            const unsigned int n_dependent_variables =
              dealii::numbers::invalid_unsigned_int,
            const bool clear_registered_tapes = true) override;

      /**
       * Register the complete set of independent variables $\mathbf{X}$.
       *
       * @param[in] values A field that defines the values of all independent
       * variables. When considering taped AD numbers with branching functions,
       * to avoid potential issues with branch switching it may be a good idea
       * to choose these values close or equal to those that will be later
       * evaluated and differentiated around.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      register_independent_variables(const std::vector<scalar_type> &values);

      /**
       * Register the subset of independent variables
       * $\mathbf{A} \subset \mathbf{X}$.
       *
       * @param[in] value A field that defines a number of independent
       * variables. When considering taped AD numbers with branching functions,
       * to avoid potential issues with branch switching it may be a good idea
       * to choose these values close or equal to those that will be later
       * evaluated and differentiated around.
       * @param[in] extractor An extractor associated with the input field
       * variables. This effectively defines which components of the global set
       * of independent variables this field is associated with.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note The input extractor must correspond to the input @p ValueType.
       * So, for example, if a value is a rank-1 tensor
       * (i.e. of type Tensor<1,dim,scalar_type>), then the extractor must
       * be an FEValuesExtractors::Vector or FEValuesExtractors::Tensor<1>.
       *
       * @note This function may be repeatedly used until a call to
       * finalize_sensitive_independent_variables() or
       * get_sensitive_variables() is made.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      template <typename ValueType, typename ExtractorType>
      void
      register_independent_variable(const ValueType &    value,
                                    const ExtractorType &extractor);

      /**
       * Return the complete set of independent variables as represented by
       * auto-differentiable numbers. These are the independent
       * variables $\mathbf{X}$ at which the dependent values are evaluated
       * and differentiated.
       *
       * This function indicates to the AD library that implements the
       * auto-differentiable number type that operations performed on these
       * numbers are to be tracked so they are considered "sensitive"
       * variables. This is, therefore, the set of variables with which one
       * would then perform computations, and based on which one can then
       * extract both the value of the function and its derivatives with the
       * member functions below. The values of the components of the returned
       * object are initialized to the values set with
       * register_independent_variable().
       *
       * @return An array of auto-differentiable type numbers.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      const std::vector<ad_type> &
      get_sensitive_variables() const;

      /*
       * Extract a subset of the independent variables as represented by
       * auto-differentiable numbers. These are the independent
       * variables $\mathbf{A} \subset \mathbf{X}$ at which the dependent values
       * are evaluated and differentiated.
       *
       * This function indicates to the AD library that implements the
       * auto-differentiable number type that operations performed on these
       * numbers are to be tracked so they are considered "sensitive"
       * variables. This is, therefore, the set of variables with which one
       * would then perform computations, and based on which one can then
       * extract both the value of the function and its derivatives with the
       * member functions below. The values of the components of the returned
       * object are initialized to the values set with
       * register_independent_variable().
       *
       * @param[in] extractor An extractor associated with the input field
       * variables. This effectively defines which components of the global set
       * of independent variables this field is associated with.
       * @return An object of auto-differentiable type numbers. The return type is
       * based on the input extractor, and will be either a scalar,
       * Tensor<1,dim>, Tensor<2,dim>, or SymmetricTensor<2,dim>.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      template <typename ExtractorType>
      typename internal::Extractor<dim,
                                   ExtractorType>::template tensor_type<ad_type>
      get_sensitive_variables(const ExtractorType &extractor) const;

      //@}

      /**
       * @name Operations specific to taped mode: Reusing tapes
       */
      //@{

      /**
       * Set the values for the independent variables $\mathbf{X}$.
       *
       * @param[in] values A vector that defines the values of all
       * independent variables.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note If the @p keep_independent_values flag has been set when
       * HelperBase::start_recording_operations() is called then the tape is
       * immediately usable after creation, and the values of the independent
       * variables set by register_independent_variables() are those at which
       * the function is to be evaluated. In this case, a separate call to this
       * function is not strictly necessary.
       */
      void
      set_independent_variables(const std::vector<scalar_type> &values);

      /**
       * Set the values for a subset of independent variables
       * $\mathbf{A} \subset \mathbf{X}$.
       *
       * @param[in] value A field that defines the values of a number of
       * independent variables.
       * @param[in] extractor An extractor associated with the input field
       * variables. This effectively defines which components of the global set
       * of independent variables this field is associated with.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note The input extractor must correspond to the input @p ValueType.
       * So, for example, if a value is a rank-1 tensor
       * (i.e. of type Tensor<1,dim,scalar_type>), then the extractor must
       * be an FEValuesExtractors::Vector or FEValuesExtractors::Tensor<1>.
       *
       * @note If the @p keep_independent_values flag has been set when
       * HelperBase::start_recording_operations() is called then the tape is
       * immediately usable after creation, and the values of the independent
       * variables set by register_independent_variable() are those at which the
       * function is to be evaluated. In this case, a separate call to this
       * function is not strictly necessary.
       */
      template <typename ValueType, typename ExtractorType>
      void
      set_independent_variable(const ValueType &    value,
                               const ExtractorType &extractor);

      //@}

    protected:
      /**
       * @name Independent variables
       */
      //@{

      /**
       * Set the actual value of the independent variable $X_{i}$.
       *
       * @param[in] index The index in the vector of independent variables.
       * @param[in] symmetric_component Mark whether this index relates to a
       * component of a field that has a symmetric counterpart
       * (e.g. if @p index represents an off-diagonal entry in a symmetric
       * tensor).
       * @param[in] value The value to set the index'd independent variable to.
       */
      void
      set_sensitivity_value(const unsigned int index,
                            const bool         symmetric_component,
                            const scalar_type &value);

      /**
       * Return whether the @p index'th independent variables is one for which
       * we must take into account symmetry when extracting their gradient or
       * Hessian values.
       */
      bool
      is_symmetric_independent_variable(const unsigned int index) const;

      /**
       * Return the number of independent variables that have been marked as
       * being components of a symmetric field.
       */
      unsigned int
      n_symmetric_independent_variables() const;

      //@}

    private:
      /**
       * @name Independent variables
       */
      //@{

      /**
       * The independent variables for which we must take into account symmetry
       * when extracting their gradient or Hessian values.
       */
      std::vector<bool> symmetric_independent_variables;

      //@}

    }; // class PointLevelFunctionsBase



    /**
     * A helper class that facilitates the evaluation of a scalar function,
     * its first derivatives (gradient), and its second derivatives (Hessian).
     * This class would typically be used to compute the first and second
     * derivatives of a <b>stored energy function</b> defined at a quadrature
     * point. It can also be used to compute derivatives of any other scalar
     * field so long as all its dependencies on the independent variables are
     * explicit (that is to say, no independent variables may have some implicit
     * dependence on one another).
     *
     * An example of its usage in the case of a multi-field constitutive law
     * might be as follows:
     * @code
     *   // Define some extractors that will help us set independent variables
     *   // and later get the computed values related to the dependent
     *   // variables. Each of these extractors is related to the gradient of a
     *   // component of the solution field (in this case, displacement and
     *   // magnetic scalar potential). Here "C" is the right Cauchy-Green
     *   // tensor and "H" is the magnetic field.
     *   const FEValuesExtractors::SymmetricTensor<2> C_dofs (0);
     *   const FEValuesExtractors::Vector             H_dofs
     *     (dealii::SymmetricTensor<2,dim>::n_independent_components);
     *   const unsigned int n_independent_variables =
     *     SymmetricTensor<2,dim>::n_independent_components +
     *     Tensor<1,dim>::n_independent_components;
     *
     *   // Define the helper that we will use in the AD computations for our
     *   // scalar energy function. Note that we expect it to return values of
     *   // type double.
     *   ScalarFunction<dim,...> ad_helper (n_independent_variables);
     *   using ADNumberType = typename ADHelper::ad_type;
     *
     *   // Compute the fields that provide the independent values.
     *   // When the tape is being replayed, these should be set to something
     *   // meaningful.
     *   const Tensor<1,dim> H = ...;
     *   const SymmetricTensor<2,dim> C = ...;
     *
     *   // If using a taped AD number, then at this point we would initiate
     *   // taping of the expression for the material stored energy function
     *   // for this particular set of material parameters:
     *
     *   // Select a tape number to record to
     *   const typename Types<ADNumberType>::tape_index tape_index = ...;
     *
     *   // Indicate that we are about to start tracing the operations for
     *   // function evaluation on the tape. If this tape has already been
     *   // used (i.e. the operations are already recorded) then we
     *   // (optionally) load the tape and reuse this data.
     *   const bool is_recording
     *     = ad_helper.start_recording_operations(tape_index);
     *
     *   // The steps that follow in the recording phase are required for
     *   // tapeless methods as well.
     *   if (is_recording == true)
     *   {
     *     // This is the "recording" phase of the operations.
     *
     *     // First, we set the values for all fields.
     *     // These could happily be set to anything, unless the function will
     *     // be evaluated along a branch not otherwise traversed during later
     *     // use. For this reason, in this example instead of using some dummy
     *     // values, we'll actually map out the function at the same point
     *     // around which we'll later linearize it.
     *     ad_helper.register_independent_variable(H, H_dofs);
     *     ad_helper.register_independent_variable(C, C_dofs);
     *
     *     // NOTE: We have to extract the sensitivities in the order we wish to
     *     // introduce them. So this means we have to do it by logical order
     *     // of the extractors that we've created.
     *     const SymmetricTensor<2,dim,ADNumberType> C_AD =
     *       ad_helper.get_sensitive_variables(C_dofs);
     *     const Tensor<1,dim,ADNumberType>          H_AD =
     *       ad_helper.get_sensitive_variables(H_dofs);
     *
     *     // Here we define the material stored energy function.
     *     // This example is sufficiently complex to warrant the use of AD to,
     *     // at the very least, verify an unassisted implementation.
     *     const double mu_e = 10;          // Shear modulus
     *     const double lambda_e = 15;      // Lam&eacute; parameter
     *     const double mu_0 = 4*M_PI*1e-7; // Magnetic permeability constant
     *     const double mu_r = 5;           // Relative magnetic permeability
     *
     *     const ADNumberType J = std::sqrt(determinant(C_AD));
     *     const SymmetricTensor<2,dim,ADNumberType> C_inv_AD = invert(C_AD);
     *     const ADNumberType psi =
     *       0.5*mu_e*(1.0+std::tanh((H_AD*H_AD)/100.0))*
     *         (trace(C_AD) - dim - 2*std::log(J)) +
     *       lambda_e*std::log(J)*std::log(J) -
     *       0.5*mu_0*mu_r*J*H_AD*C_inv_AD*H_AD;
     *
     *     // Register the definition of the total stored energy
     *     ad_helper.register_dependent_variable(psi_CH);
     *
     *     // Indicate that we have completed tracing the operations onto
     *     // the tape.
     *     ad_helper.stop_recording_operations(false); // write_tapes_to_file
     *   }
     *   else
     *   {
     *     // This is the "tape reuse" phase of the operations.
     *     // Here we will leverage the already traced operations that reside
     *     // on a tape, and simply re-evaluate the tape at a different point
     *     // to get the function values and their derivatives.
     *
     *     // Load the existing tape to be reused
     *     ad_helper.activate_recorded_tape(tape_index);
     *
     *     // Set the new values of the independent variables where the
     *     // recorded dependent functions are to be evaluated (and
     *     // differentiated around).
     *     ad_helper.set_independent_variable(C, C_dofs);
     *     ad_helper.set_independent_variable(H, H_dofs);
     *   }
     *
     *   // Play the tape and store the output function value, its gradient and
     *   // linearization. These are expensive to compute, so we'll do this once
     *   // and extract the desired values from these intermediate outputs.
     *   Vector<double> Dpsi (ad_helper.n_dependent_variables());
     *   FullMatrix<double> D2psi (ad_helper.n_dependent_variables(),
     *                             ad_helper.n_independent_variables());
     *   const double psi = ad_helper.compute_value();
     *   ad_helper.compute_gradient(Dpsi);
     *   ad_helper.compute_hessian(D2psi);
     *
     *   // Extract the desired components of the gradient vector and Hessian
     *   // matrix. In this example, we use them to compute the Piola-Kirchhoff
     *   // stress tensor and its associated tangent, defined by thermodynamic
     *   // arguments as S = 2*dpsi/dC and HH = 2*dS/dC...
     *   const SymmetricTensor<2,dim> S =
     *     2.0*ad_helper.extract_gradient_component(Dpsi,C_dofs);
     *   const SymmetricTensor<4,dim> HH =
     *     4.0*ad_helper.extract_hessian_component(D2psi,C_dofs,C_dofs);
     *
     *   // ... the magnetic induction and its associated tangent defined
     *   // as B = -dpsi/dH and BB = dB/dH...
     *   const Tensor<1,dim> B =
     *     -ad_helper.extract_gradient_component(Dpsi,H_dofs);
     *   const SymmetricTensor<2,dim> BB =
     *     -symmetrize(ad_helper.extract_hessian_component(D2psi,H_dofs,H_dofs));
     *
     *   // ... and finally the magnetoelastic coupling tangent, defined
     *   // as PP = -dS/dH = -d/dH(2*dpsi/dC). Here the order of the extractor
     *   // arguments is especially important, as it dictates the order in which
     *   // the directional derivatives are taken.
     *   const Tensor<3,dim,double> PP =
     *     -2.0*ad_helper.extract_hessian_component(D2psi,C_dofs,H_dofs)
     * @endcode
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class ScalarFunction
      : public PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{f}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       */
      ScalarFunction(const unsigned int n_independent_variables);

      /**
       * Destructor.
       */
      virtual ~ScalarFunction() = default;

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Register the definition of the scalar field $\Psi(\mathbf{X})$.
       *
       * @param[in] func The recorded function that defines a dependent
       * variable.
       *
       * @note For this class that expects only one dependent variable, this
       * function must only be called once per tape.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      register_dependent_variable(const ad_type &func);

      /**
       * Compute the value of the scalar field $\Psi(\mathbf{X})$ using the
       * tape as opposed to executing the source code.
       *
       * @return A scalar object with the value for the scalar field evaluated
       * at the point defined by the independent variable values.
       */
      scalar_type
      compute_value() const;

      /**
       * Compute the gradient (first derivative) of the scalar field with
       * respect to all independent variables, i.e.
       * @f[
       *   \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{X}}
       * @f]
       *
       * @param[out] gradient A Vector with the values for the scalar field
       * gradient (first derivatives) evaluated at the point defined by the
       * independent variable values. The output @p gradient vector has a length
       * corresponding to @p n_independent_variables.
       */
      void
      compute_gradient(Vector<scalar_type> &gradient) const;

      /**
       * Compute the Hessian (second derivative)  of the scalar field with
       * respect to all independent variables, i.e.
       * @f[
       *   \frac{\partial^{2}\Psi(\mathbf{X})}{\partial\mathbf{X} \otimes
       * \partial\mathbf{X}}
       * @f]
       *
       * @param[out] hessian A FullMatrix with the values for the scalar field
       * Hessian (second derivatives) evaluated at the point defined by the
       * independent variable values. The output @p hessian matrix has
       * dimensions corresponding to
       * <code>n_independent_variables</code>$\times$<code>n_independent_variables</code>.
       */
      void
      compute_hessian(FullMatrix<scalar_type> &hessian) const;

      /**
       * Extract the function gradient for a subset of independent variables
       * $\mathbf{A} \subset \mathbf{X}$, i.e.
       * @f[
       *   \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}}
       * @f]
       *
       * @param[in] gradient The gradient of the scalar function with respect to
       * all independent variables, i.e., that returned by compute_gradient().
       * @param[in] extractor_row An extractor associated with the input field
       * variables. This effectively defines which components of the global set
       * of independent variables this field is associated with.
       *
       * @return A Tensor or SymmetricTensor with its rank and symmetries
       * determined by the @p extractor_row.
       * This corresponds to subsetting a whole set of rows of the
       * gradient vector, scaling those entries to take account of tensor
       * symmetries, and then reshaping the (sub-)vector so obtained into a
       * tensor, the final result.
       * For example, if
       * @p extractor_row is a FEValuesExtractors::Vector and
       * @p extractor_col is a FEValuesExtractors::Tensor,
       * then the returned object is a Tensor of rank 3, with its first
       * index associated with the field corresponding to the row extractor and
       * the second and third indices associated with the field corresponding to
       * the column extractor.
       * Similarly, if
       * @p extractor_row is a FEValuesExtractors::SymmetricTensor and
       * @p extractor_col is a FEValuesExtractors::SymmetricTensor,
       * then the returned object is a SymmetricTensor of rank 4, with its first
       * two indices associated with the field corresponding to the row
       * extractor and the last two indices associated with the field
       * corresponding to the column extractor.
       */
      template <typename ExtractorType_Row>
      static typename internal::
        ScalarFieldGradient<dim, scalar_type, ExtractorType_Row>::type
        extract_gradient_component(const Vector<scalar_type> &gradient,
                                   const ExtractorType_Row &  extractor_row);

      /**
       * Extract the function Hessian for a subset of independent variables
       * $\mathbf{A},\mathbf{B} \subset \mathbf{X}$, i.e.
       * @f[
       *   \frac{}{\partial\mathbf{B}} \left[
       * \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}} \right] =
       * \frac{\partial^{2}\Psi(\mathbf{X})}{\partial\mathbf{B} \otimes
       * \partial\mathbf{A}}
       * @f]
       *
       * @param[in] hessian The Hessian of the scalar function with respect to
       * all independent variables, i.e., that returned by compute_hessian().
       * @param[in] extractor_row An extractor associated with the input field
       * variables for which the first index of the Hessian is extracted.
       * @param[in] extractor_col An extractor associated with the input field
       * variables for which the second index of the Hessian is extracted.
       *
       * @return A Tensor or SymmetricTensor with its rank and symmetries
       * determined by the @p extractor_row and @p extractor_col .
       * This corresponds to subsetting a whole set of rows and columns of the
       * Hessian matrix, scaling those entries to take account of tensor
       * symmetries, and then reshaping the (sub-)matrix so obtained into a
       * tensor, the final result.
       * For example, if
       * @p extractor_row is a FEValuesExtractors::Vector and
       * @p extractor_col is a FEValuesExtractors::Tensor,
       * then the returned object is a Tensor of rank 3, with its first
       * index associated with the field corresponding to the row extractor and
       * the second and third indices associated with the field corresponding to
       * the column extractor.
       * Similarly, if
       * @p extractor_row is a FEValuesExtractors::SymmetricTensor and
       * @p extractor_col is a FEValuesExtractors::SymmetricTensor,
       * then the returned object is a SymmetricTensor of rank 4, with its first
       * two indices associated with the field corresponding to the row
       * extractor and the last two indices associated with the field
       * corresponding to the column extractor.
       */
      template <typename ExtractorType_Row, typename ExtractorType_Col>
      static typename internal::ScalarFieldHessian<dim,
                                                   scalar_type,
                                                   ExtractorType_Row,
                                                   ExtractorType_Col>::type
      extract_hessian_component(const FullMatrix<scalar_type> &hessian,
                                const ExtractorType_Row &      extractor_row,
                                const ExtractorType_Col &      extractor_col);

      /**
       * Extract the function Hessian for a subset of independent variables
       * $\mathbf{A},\mathbf{B} \subset \mathbf{X}$, i.e.
       * @f[
       *   \frac{}{\partial\mathbf{B}} \left[
       * \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}} \right]
       * @f]
       *
       * This function is a specialization of the above for rank-0 tensors
       * (scalars). This corresponds to extracting a single entry of the
       * Hessian matrix because both extractors imply selection of just a
       * single row or column of the matrix.
       */
      static Tensor<0, dim, scalar_type>
      extract_hessian_component(
        const FullMatrix<scalar_type> &   hessian,
        const FEValuesExtractors::Scalar &extractor_row,
        const FEValuesExtractors::Scalar &extractor_col);

      /**
       * Extract the function Hessian for a subset of independent variables
       * $\mathbf{A},\mathbf{B} \subset \mathbf{X}$, i.e.
       * @f[
       *   \frac{}{\partial\mathbf{B}} \left[
       * \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}} \right]
       * @f]
       *
       * This function is a specialization of the above for rank-4 symmetric
       * tensors.
       */
      static SymmetricTensor<4, dim, scalar_type>
      extract_hessian_component(
        const FullMatrix<scalar_type> &               hessian,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_row,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_col);

      //@}

    }; // class ScalarFunction



    /**
     * A helper class that facilitates the evaluation of a vector of functions,
     * typically one that represents a collection of coupled, multi-dimensional
     * fields. This class would typically be used to compute the linearization
     * of a set of kinetic field variables defined at the quadrature point
     * level.
     *
     * An example of its usage in the case of linearizing the kinetic variables
     * derived from a multi-field constitutive law might be as follows:
     * @code
     *   // Define some extractors that will help us set independent variables
     *   // and later get the computed values related to the dependent
     *   // variables. Each of these extractors is related to the gradient of a
     *   // component of the solution field (in this case, displacement and
     *   // magnetic scalar potential). Here "C" is the right Cauchy-Green
     *   // tensor and "H" is the magnetic field.
     *   const FEValuesExtractors::SymmetricTensor<2> C_dofs (0);
     *   const FEValuesExtractors::Vector             H_dofs
     *     (dealii::SymmetricTensor<2,dim>::n_independent_components);
     *   const unsigned int n_independent_variables =
     *     SymmetricTensor<2,dim>::n_independent_components +
     *     Tensor<1,dim>::n_independent_components;
     *
     *   // Declare how many dependent variables we expect to compute.
     *   // In this case, we will be computing a stress field (a symmetric
     *   // rank-2 tensor) and the magnetic induction (a vector field).
     *   // At the same time we define some additional extractors associated
     *   // with these kinetic fields. In general, these need not be of the same
     *   // layout as the independent variables.
     *   const FEValuesExtractors::SymmetricTensor<2> S_dofs (0);
     *   const FEValuesExtractors::Vector             B_dofs
     *     (dealii::SymmetricTensor<2,dim>::n_independent_components);
     *   const unsigned int n_dependent_variables =
     *     SymmetricTensor<2,dim>::n_independent_components +
     *     Tensor<1,dim>::n_independent_components;
     *
     *   // Define the helper that we will use in the AD computations for our
     *   // scalar energy function. Note that we expect it to return values of
     *   // type double.
     *   VectorFunction<dim,double> ad_helper (n_independent_variables,
     *                                                 n_dependent_variables);
     *   using ADNumberType = typename ADHelper::ad_type;
     *
     *   // Compute the fields that provide the independent values.
     *   // When the tape is being replayed, these should be set to something
     *   // meaningful.
     *   const Tensor<1,dim> H = ...;
     *   const SymmetricTensor<2,dim> C = ...;
     *
     *   // If using a taped AD number, then at this point we would initiate
     *   // taping of the expression for the material stored energy function
     *   // for this particular set of material parameters:
     *
     *   // Select a tape number to record to
     *   const typename Types<ADNumberType>::tape_index tape_index = ...;
     *
     *   // Indicate that we are about to start tracing the operations for
     *   // function evaluation on the tape. If this tape has already been
     *   // used (i.e. the operations are already recorded) then we
     *   // (optionally) load the tape and reuse this data.
     *   const bool is_recording
     *     = ad_helper.start_recording_operations(tape_index);
     *
     *   // The steps that follow in the recording phase are required for
     *   // tapeless methods as well.
     *   if (is_recording == true)
     *   {
     *     // This is the "recording" phase of the operations.
     *
     *     // First, we set the values for all fields.
     *     // These could happily be set to anything, unless the function will
     *     // be evaluated along a branch not otherwise traversed during later
     *     // use. For this reason, in this example instead of using some dummy
     *     // values, we'll actually map out the function at the same point
     *     // around which we'll later linearize it.
     *     ad_helper.register_independent_variable(H, H_dofs);
     *     ad_helper.register_independent_variable(C, C_dofs);
     *
     *     // NOTE: We have to extract the sensitivities in the order we wish to
     *     // introduce them. So this means we have to do it by logical order
     *     // of the extractors that we've created.
     *     const SymmetricTensor<2,dim,ADNumberType> C_AD =
     *       ad_helper.get_sensitive_variables(C_dofs);
     *     const Tensor<1,dim,ADNumberType>          H_AD =
     *       ad_helper.get_sensitive_variables(H_dofs);
     *
     *     // Here we define the stress and magnetic induction in terms
     *     // of the independent values C_AD and H_AD.
     *     const SymmetricTensor<2, dim, ad_type> S_AD = ...;
     *     const Tensor<1, dim, ad_type>          B_AD = ...;
     *
     *     // Register the definition of the kinetic fields. The second
     *     // argument to the function provides a non-overlapping ordering
     *     // of the
     *     ad_helper.register_dependent_variable(S_AD, S_dofs);
     *     ad_helper.register_dependent_variable(B_AD, B_dofs);
     *
     *     // Indicate that we have completed tracing the operations onto
     *     // the tape.
     *     ad_helper.stop_recording_operations(false); // write_tapes_to_file
     *   }
     *   else
     *   {
     *     // This is the "tape reuse" phase of the operations.
     *     // Here we will leverage the already traced operations that reside
     *     // on a tape, and simply re-evaluate the tape at a different point
     *     // to get the function values and their derivatives.
     *
     *     // Load the existing tape to be reused
     *     ad_helper.activate_recorded_tape(tape_index);
     *
     *     // Set the new values of the independent variables where the
     *     // recorded dependent functions are to be evaluated (and
     *     // differentiated around).
     *     ad_helper.set_independent_variable(C, C_dofs);
     *     ad_helper.set_independent_variable(H, H_dofs);
     *   }
     *
     *   // Play the tape and store the output function value, its gradient and
     *   // linearization. These are expensive to compute, so we'll do this once
     *   // and extract the desired values from these intermediate outputs.
     *   Vector<double> values (ad_helper.n_dependent_variables());
     *   FullMatrix<double> jacobian (ad_helper.n_dependent_variables(),
     *                                ad_helper.n_independent_variables());
     *   ad_helper.compute_values(values);
     *   ad_helper.compute_jacobian(jacobian);
     *
     *   // Extract the desired components of the value vector and its Jacobian
     *   // matrix. In this example, we use them to compute the Piola-Kirchhoff
     *   // stress tensor S and its associated tangent, defined
     *   // as HH = 2*dS/dC...
     *   const SymmetricTensor<2,dim> S =
     *     ad_helper.extract_value_component(values,S_dofs);
     *   const SymmetricTensor<4,dim> HH =
     *     2.0*ad_helper.extract_jacobian_component(jacobian,S_dofs,C_dofs);
     *
     *   // ... the magnetic induction B and its associated tangent defined
     *   // as BB = dB/dH...
     *   const Tensor<1,dim> B =
     *     ad_helper.extract_value_component(values,H_dofs);
     *   const SymmetricTensor<2,dim> BB =
     *     symmetrize(ad_helper.extract_jacobian_component(jacobian,B_dofs,H_dofs));
     *
     *   // ... and finally the magnetoelastic coupling tangent, defined
     *   // as PP = -dS/dH. Here the order of the extractor arguments is
     *   // especially important, as it dictates the field that is being
     *   // differentiated, and which the directional derivatives are being
     *   // computed.
     *   const Tensor<3,dim,double> PP =
     *     -ad_helper.extract_jacobian_component(jacobian,S_dofs,H_dofs)
     * @endcode
     *
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class VectorFunction
      : public PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
       * The constructor for the class.
       *
       * @param[in] n_independent_variables The number of independent variables
       * that will be used in the definition of the functions that it is
       * desired to compute the sensitivities of. In the computation of
       * $\mathbf{f}(\mathbf{X})$, this will be the number of inputs
       * $\mathbf{X}$, i.e., the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e., the dimension of the
       * image space.
       */
      VectorFunction(const unsigned int n_independent_variables,
                     const unsigned int n_dependent_variables);

      /**
       * Destructor.
       */
      virtual ~VectorFunction() = default;

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Register the definition of the vector field
       * $\boldsymbol{\Psi}(\mathbf{X})$.
       *
       * @param[in] funcs A vector of recorded functions that defines the
       * dependent variables.
       *
       * @note For this class that expects only vector field of dependent
       * variables, this function must only be called once per tape.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      void
      register_dependent_variables(const std::vector<ad_type> &funcs);

      /**
       * Register the definition of the vector field
       * $\hat{\mathbf{g}}(\mathbf{X}) \subset \boldsymbol{\Psi}(\mathbf{X})$
       * that may represent a subset of the dependent variables.
       *
       * @param[in] funcs The recorded functions that define a set of dependent
       * variables.
       * @param[in] extractor An extractor associated with the input field
       * variables. This effectively defines which components of the global set
       * of dependent variables this field is associated with.
       *
       * @note The input extractor must correspond to the input @p ValueType.
       * So, for example, if a value is a rank-1 tensor
       * (i.e. of type Tensor<1,dim,scalar_type>), then the extractor must
       * be an FEValuesExtractors::Vector or FEValuesExtractors::Tensor<1>.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      template <typename ValueType, typename ExtractorType>
      void
      register_dependent_variable(const ValueType &    funcs,
                                  const ExtractorType &extractor);

      /**
       * Compute the value of the vector field $\boldsymbol{\Psi}(\mathbf{X})$.
       *
       * @param[out] values A Vector object with the value for each component
       * of the vector field evaluated at the point defined by the independent
       * variable values. The output @p values vector has a length
       * corresponding to @p n_dependent_variables.
       */
      void
      compute_values(Vector<scalar_type> &values) const;

      /**
       * Compute the Jacobian (first derivative) of the vector field with
       * respect to all independent variables, i.e.
       * @f[
       *   \mathbf{J}(\boldsymbol{\Psi})
       *      = \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{X}}
       * @f]
       *
       * @param[out] jacobian A FullMatrix with the gradient of each component
       * of the vector field evaluated at the point defined by the independent
       * variable values. The output @p jacobian matrix has
       * dimensions corresponding to
       * <code>n_dependent_variables</code>$\times$<code>n_independent_variables</code>.
       */
      void
      compute_jacobian(FullMatrix<scalar_type> &jacobian) const;


      /**
       * Extract the set of functions' values for a subset of dependent
       * variables
       * $\mathbf{g} \subset \boldsymbol{\Psi}(\mathbf{X})$.
       *
       * @param[in] values A Vector object with the value for each component of
       * the vector field evaluated at the point defined by the independent
       * variable values.
       * @param[in] extractor_row An extractor associated with the input field
       * variables. This effectively defines which components of the global set
       * of dependent variables this field is associated with.
       */
      template <typename ExtractorType_Row>
      static typename internal::
        VectorFieldValue<dim, scalar_type, ExtractorType_Row>::type
        extract_value_component(const Vector<scalar_type> &values,
                                const ExtractorType_Row &  extractor_row);

      /**
       * Extract the Jacobian of the subset of dependent functions
       * $\mathbf{g} \subset \boldsymbol{\Psi}(\mathbf{X})$
       * for a subset of independent variables
       * $\mathbf{A} \subset \mathbf{X}$, i.e.
       * @f[
       *   \mathbf{J}(\mathbf{g})
       *      = \frac{\partial\mathbf{g}(\mathbf{X})}{\partial\mathbf{A}}
       * @f]
       * The first index of the Jacobian matrix $\mathbf{J}(\mathbf{g})$
       * relates to the dependent variables, while the second index relates
       * to the independent variables.
       *
       * @param[in] jacobian The Jacobian of the vector function with respect to
       * all independent variables, i.e., that returned by compute_jacobian().
       * @param[in] extractor_row An extractor associated with the input field
       * variables for which the first index of the Jacobian is extracted.
       * This effectively defines the correspondence between components of the
       * global set of dependent variables and the field (representing a
       * subset of dependent functions) associated with the extractor.
       * @param[in] extractor_col An extractor associated with the input field
       * variables for which the second index of the Jacobian is extracted.
       * This effectively defines the correspondence between components of the
       * global set of independent variables and the field (representing a
       * subset of independent variables) associated with the extractor.
       *
       * @return A Tensor or SymmetricTensor with its rank and symmetries
       * determined by the @p extractor_row and @p extractor_col .
       * This corresponds to subsetting a whole set of rows and columns of the
       * Jacobian matrix, scaling those entries to take account of tensor
       * symmetries, and then reshaping the (sub-)matrix so obtained into a
       * tensor, the final result.
       * For example, if
       * @p extractor_row is a FEValuesExtractors::Vector and
       * @p extractor_col is a FEValuesExtractors::Tensor,
       * then the returned object is a Tensor of rank 3, with its first
       * index associated with the field corresponding to the row extractor and
       * the second and third indices associated with the field corresponding to
       * the column extractor.
       * Similarly, if
       * @p extractor_row is a FEValuesExtractors::SymmetricTensor and
       * @p extractor_col is a FEValuesExtractors::SymmetricTensor,
       * then the returned object is a SymmetricTensor of rank 4, with its first
       * two indices associated with the field corresponding to the row
       * extractor and the last two indices associated with the field
       * corresponding to the column extractor.
       */
      template <typename ExtractorType_Row, typename ExtractorType_Col>
      static typename internal::VectorFieldJacobian<dim,
                                                    scalar_type,
                                                    ExtractorType_Row,
                                                    ExtractorType_Col>::type
      extract_jacobian_component(const FullMatrix<scalar_type> &jacobian,
                                 const ExtractorType_Row &      extractor_row,
                                 const ExtractorType_Col &      extractor_col);

      /**
       * Extract the Jacobian of the subset of dependent functions
       * $\mathbf{g} \subset \boldsymbol{\Psi}(\mathbf{X})$
       * for a subset of independent variables
       * $\mathbf{A} \subset \mathbf{X}$, i.e.
       * @f[
       *   \mathbf{J}(\mathbf{g})
       *      = \frac{\partial\mathbf{g}(\mathbf{X})}{\partial\mathbf{A}}
       * @f]
       *
       * This function is a specialization of the above for rank-0 tensors
       * (scalars). This corresponds to extracting a single entry of the
       * Jacobian matrix because both extractors imply selection of just a
       * single row or column of the matrix.
       */
      static Tensor<0, dim, scalar_type>
      extract_jacobian_component(
        const FullMatrix<scalar_type> &   jacobian,
        const FEValuesExtractors::Scalar &extractor_row,
        const FEValuesExtractors::Scalar &extractor_col);

      /**
       * Extract the Jacobian of the subset of dependent functions
       * $\mathbf{g} \subset \boldsymbol{\Psi}(\mathbf{X})$
       * for a subset of independent variables
       * $\mathbf{A} \subset \mathbf{X}$, i.e.
       * @f[
       *   \mathbf{J}(\mathbf{g})
       *      = \frac{\partial\mathbf{g}(\mathbf{X})}{\partial\mathbf{A}}
       * @f]
       *
       * This function is a specialization of the above for rank-4 symmetric
       * tensors.
       */
      static SymmetricTensor<4, dim, scalar_type>
      extract_jacobian_component(
        const FullMatrix<scalar_type> &               jacobian,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_row,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_col);

      //@}

    }; // class VectorFunction


  } // namespace AD
} // namespace Differentiation


/* ----------------- inline and template functions ----------------- */


#  ifndef DOXYGEN

namespace Differentiation
{
  namespace AD
  {
    /* ----------------- CellLevelBase ----------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template <typename VectorType>
    void
    CellLevelBase<ADNumberTypeCode, ScalarType>::register_dof_values(
      const VectorType &                                  values,
      const std::vector<dealii::types::global_dof_index> &local_dof_indices)
    {
      // This is actually the same thing the set_dof_values() function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
      Assert(
        local_dof_indices.size() == this->n_independent_variables(),
        ExcMessage(
          "Degree of freedom index vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        {
          Assert(this->registered_independent_variable_values[i] == false,
                 ExcMessage("Independent variables already registered."));
        }
      set_dof_values(values, local_dof_indices);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template <typename VectorType>
    void
    CellLevelBase<ADNumberTypeCode, ScalarType>::set_dof_values(
      const VectorType &                                  values,
      const std::vector<dealii::types::global_dof_index> &local_dof_indices)
    {
      Assert(local_dof_indices.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        HelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
          i, values[local_dof_indices[i]]);
    }



    /* ----------------- PointLevelFunctionsBase ----------------- */



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ValueType, typename ExtractorType>
    void
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      register_independent_variable(const ValueType &    value,
                                    const ExtractorType &extractor)
    {
      // This is actually the same thing as the set_independent_variable
      // function, in the sense that we simply populate our array of independent
      // values with a meaningful number. However, in this case we need to
      // double check that we're not registering these variables twice
#    ifdef DEBUG
      const std::vector<unsigned int> index_set(
        internal::extract_field_component_indices<dim>(extractor));
      for (const unsigned int index : index_set)
        {
          Assert(
            this->registered_independent_variable_values[index] == false,
            ExcMessage(
              "Overlapping indices for independent variables. "
              "One or more indices associated with the field that "
              "is being registered as an independent variable have "
              "already been associated with another field. This suggests "
              "that the component offsets used to construct their counterpart "
              "extractors are incompatible with one another. Make sure that "
              "the first component for each extractor properly takes into "
              "account the dimensionality of the preceding fields."));
        }
#    endif
      set_independent_variable(value, extractor);
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ValueType, typename ExtractorType>
    void
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      set_independent_variable(const ValueType &    value,
                               const ExtractorType &extractor)
    {
      const std::vector<unsigned int> index_set(
        internal::extract_field_component_indices<dim>(extractor));
      for (unsigned int i = 0; i < index_set.size(); ++i)
        {
          set_sensitivity_value(
            index_set[i],
            internal::Extractor<dim, ExtractorType>::symmetric_component(i),
            internal::get_tensor_entry(value, i));
        }
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ExtractorType>
    typename internal::Extractor<dim, ExtractorType>::template tensor_type<
      typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type>
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      get_sensitive_variables(const ExtractorType &extractor) const
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // If necessary, finalize the internally stored vector of
      // AD numbers that represents the independent variables
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size() ==
               this->n_independent_variables(),
             ExcDimensionMismatch(this->independent_variables.size(),
                                  this->n_independent_variables()));

      const std::vector<unsigned int> index_set(
        internal::extract_field_component_indices<dim>(extractor));
      typename internal::Extractor<dim,
                                   ExtractorType>::template tensor_type<ad_type>
        out;

      for (unsigned int i = 0; i < index_set.size(); ++i)
        {
          const unsigned int index = index_set[i];
          Assert(index < this->n_independent_variables(), ExcInternalError());
          Assert(this->registered_independent_variable_values[index] == true,
                 ExcInternalError());
          internal::get_tensor_entry(out, i) =
            this->independent_variables[index];
        }

      return out;
    }



    /* ----------------- ScalarFunction ----------------- */



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ExtractorType_Row>
    typename internal::ScalarFieldGradient<
      dim,
      typename ScalarFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type,
      ExtractorType_Row>::type
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_gradient_component(const Vector<scalar_type> &gradient,
                                 const ExtractorType_Row &  extractor_row)
    {
      // NOTE: The order of components must be consistently defined throughout
      // this class.
      typename internal::
        ScalarFieldGradient<dim, scalar_type, ExtractorType_Row>::type out;

      // Get indexsets for the subblock from which we wish to extract the
      // gradient values
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(extractor_row));
      Assert(out.n_independent_components == row_index_set.size(),
             ExcMessage("Not all tensor components have been extracted!"));
      for (unsigned int r = 0; r < row_index_set.size(); ++r)
        internal::set_tensor_entry(out, r, gradient[row_index_set[r]]);

      return out;
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ExtractorType_Row, typename ExtractorType_Col>
    typename internal::ScalarFieldHessian<
      dim,
      typename ScalarFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type,
      ExtractorType_Row,
      ExtractorType_Col>::type
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_hessian_component(const FullMatrix<scalar_type> &hessian,
                                const ExtractorType_Row &      extractor_row,
                                const ExtractorType_Col &      extractor_col)
    {
      using InternalHessian      = internal::ScalarFieldHessian<dim,
                                                           scalar_type,
                                                           ExtractorType_Row,
                                                           ExtractorType_Col>;
      using InternalExtractorRow = internal::Extractor<dim, ExtractorType_Row>;
      using InternalExtractorCol = internal::Extractor<dim, ExtractorType_Col>;
      using HessianType          = typename InternalHessian::type;

      // NOTE: The order of components must be consistently defined throughout
      // this class.
      HessianType out;

      // Get indexsets for the subblocks from which we wish to extract the
      // Hessian values
      // NOTE: Here we have to do some clever accounting when the
      // one extractor is a symmetric Tensor and the other is not, e.g.
      // <SymmTensor,Vector>. In this scenario the return type is a
      // non-symmetric Tensor<3,dim> but we have to fetch information from a
      // SymmTensor row/column that has too few entries to fill the output
      // tensor. So we must duplicate the relevant entries in the row/column
      // indexset to fetch off-diagonal components that are Otherwise
      // non-existent in a SymmTensor.
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(
          extractor_row, false /*ignore_symmetries*/));
      const std::vector<unsigned int> col_index_set(
        internal::extract_field_component_indices<dim>(
          extractor_col, false /*ignore_symmetries*/));

      for (unsigned int index = 0;
           index < HessianType::n_independent_components;
           ++index)
        {
          const TableIndices<HessianType::rank> ti_out =
            HessianType::unrolled_to_component_indices(index);
          const unsigned int r =
            InternalExtractorRow::local_component(ti_out, 0);
          const unsigned int c =
            InternalExtractorCol::local_component(ti_out,
                                                  InternalExtractorRow::rank);

          internal::set_tensor_entry(
            out, index, hessian[row_index_set[r]][col_index_set[c]]);
        }

      return out;
    }



    /* ----------------- VectorFunction ----------------- */



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ValueType, typename ExtractorType>
    void
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::
      register_dependent_variable(const ValueType &    funcs,
                                  const ExtractorType &extractor)
    {
      const std::vector<unsigned int> index_set(
        internal::extract_field_component_indices<dim>(extractor));
      for (unsigned int i = 0; i < index_set.size(); ++i)
        {
          Assert(this->registered_marked_dependent_variables[index_set[i]] ==
                   false,
                 ExcMessage("Overlapping indices for dependent variables."));
          HelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
            index_set[i], internal::get_tensor_entry(funcs, i));
        }
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ExtractorType_Row>
    typename internal::VectorFieldValue<
      dim,
      typename VectorFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type,
      ExtractorType_Row>::type
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::extract_value_component(
      const Vector<scalar_type> &values,
      const ExtractorType_Row &  extractor_row)
    {
      // NOTE: The order of components must be consistently defined throughout
      // this class.
      typename internal::VectorFieldValue<dim, scalar_type, ExtractorType_Row>::
        type out;

      // Get indexsets for the subblock from which we wish to extract the
      // gradient values
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(extractor_row));
      Assert(out.n_independent_components == row_index_set.size(),
             ExcMessage("Not all tensor components have been extracted!"));
      for (unsigned int r = 0; r < row_index_set.size(); ++r)
        internal::set_tensor_entry(out, r, values[row_index_set[r]]);

      return out;
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    template <typename ExtractorType_Row, typename ExtractorType_Col>
    typename internal::VectorFieldJacobian<
      dim,
      typename VectorFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type,
      ExtractorType_Row,
      ExtractorType_Col>::type
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_jacobian_component(const FullMatrix<scalar_type> &jacobian,
                                 const ExtractorType_Row &      extractor_row,
                                 const ExtractorType_Col &      extractor_col)
    {
      using InternalJacobian     = internal::VectorFieldJacobian<dim,
                                                             scalar_type,
                                                             ExtractorType_Row,
                                                             ExtractorType_Col>;
      using InternalExtractorRow = internal::Extractor<dim, ExtractorType_Row>;
      using InternalExtractorCol = internal::Extractor<dim, ExtractorType_Col>;
      using JacobianType         = typename InternalJacobian::type;

      // NOTE: The order of components must be consistently defined throughout
      // this class.
      JacobianType out;

      // Get indexsets for the subblocks from which we wish to extract the
      // Hessian values.
      // NOTE: Here we have to do some clever accounting when the
      // one extractor is a symmetric Tensor and the other is not, e.g.
      // <SymmTensor,Vector>. In this scenario the return type is a
      // non-symmetric Tensor<3,dim> but we have to fetch information from a
      // SymmTensor row/column that has too few entries to fill the output
      // tensor. So we must duplicate the relevant entries in the row/column
      // indexset to fetch off-diagonal components that are Otherwise
      // non-existent in a SymmTensor.
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(
          extractor_row, false /*ignore_symmetries*/));
      const std::vector<unsigned int> col_index_set(
        internal::extract_field_component_indices<dim>(
          extractor_col, false /*ignore_symmetries*/));

      for (unsigned int index = 0;
           index < JacobianType::n_independent_components;
           ++index)
        {
          const TableIndices<JacobianType::rank> ti_out =
            JacobianType::unrolled_to_component_indices(index);
          const unsigned int r =
            InternalExtractorRow::local_component(ti_out, 0);
          const unsigned int c =
            InternalExtractorCol::local_component(ti_out,
                                                  InternalExtractorRow::rank);

          internal::set_tensor_entry(
            out, index, jacobian[row_index_set[r]][col_index_set[c]]);
        }

      return out;
    }


  } // namespace AD
} // namespace Differentiation


#  endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)

#endif // dealii_differentiation_ad_ad_helpers_h
