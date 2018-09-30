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

#  include <deal.II/differentiation/ad/ad_drivers.h>
#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/adolc_number_types.h>
#  include <deal.II/differentiation/ad/adolc_product_types.h>
#  include <deal.II/differentiation/ad/sacado_number_types.h>
#  include <deal.II/differentiation/ad/sacado_product_types.h>

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
     *   const types::tape_index tape_index = ...;
     *
     *   // Indicate that we are about to start tracing the operations for
     *   // function evaluation on the tape. If this tape has already been used
     *   // (i.e. the operations are already recorded) then we (optionally) load
     *   // the tape and reuse this data.
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
    class ADHelperBase
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
       * $\mathbf{X}$, i.e. the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e. the dimension of the
       * image space.
       */
      ADHelperBase(const unsigned int n_independent_variables,
                   const unsigned int n_dependent_variables);

      /**
       * Destructor
       */
      virtual ~ADHelperBase();

      //@}

      /**
       * @name Interrogation of internal information
       */
      //@{

      /**
       * Returns the number of independent variables that this object expects to
       * work with. This is the dimension of the domain space.
       */
      std::size_t
      n_independent_variables() const;

      /**
       * Returns the number of dependent variables that this object expects to
       * operate on. This is the dimension of the image space.
       */
      std::size_t
      n_dependent_variables() const;

      /**
       * Prints the status of all queryable data. Exactly what is printed and
       * its format depends on the @p ad_type, as is determined by the
       * @p ADNumberTypeCode template parameter.
       *
       * @param[in] stream The output stream to which the values are to be
       * written.
       */
      void
      print(std::ostream &stream) const;

      /**
       * Prints the values currently assigned to the independent variables.
       *
       * @param[in] stream The output stream to which the values are to be
       * written.
       */
      void
      print_values(std::ostream &stream) const;

      /**
       * Prints the statistics regarding the usage of the tapes.
       *
       * @param[in] tape_index The index of the tape to get the statistics of.
       * @param[out] stream The output stream to which the values are to be
       * written.
       *
       * @note This function only produces meaningful output when @p ad_type
       * is a taped auto-differentiable number.
       */
      void
      print_tape_stats(const types::tape_index tape_index,
                       std::ostream &          stream) const;

      //@}

      /**
       * @name Operations specific to tapeless mode
       */
      //@{

      /**
       * Pre-specify the number of @p independent_variables to be used in
       * tapeless mode.
       *
       * Although this function is called internally in the ADHelperBase
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
       * When an instance of an ADHelperBase is stored as a class member object
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
       * $\mathbf{X}$, i.e. the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e. the dimension of the
       * image space.
       * @param[in] clear_registered_tapes A flag that indicates the that
       * list of @p registered_tapes must be cleared.
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
       * Returns whether or not this class is tracking calculations performed
       * with its marked independent variables.
       */
      bool
      is_recording() const;

      /**
       * Returns the tape number which is currently activated for recording or
       * reading.
       */
      types::tape_index
      active_tape() const;

      /**
       * Returns whether or not a tape number has already been used
       * or registered.
       */
      bool
      is_registered_tape(const types::tape_index tape_index) const;

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
      set_tape_buffer_sizes(const types::tape_buffer_sizes obufsize = 67108864,
                            const types::tape_buffer_sizes lbufsize = 67108864,
                            const types::tape_buffer_sizes vbufsize = 67108864,
                            const types::tape_buffer_sizes tbufsize = 67108864);

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
       * @param[in] keep_values Determines whether the numerical values of all
       * independent variables are recorded in the tape buffer. If true, then
       * the tape can be immediately used to perform computations after
       * recording is complete.
       *
       * @note During the recording phase, no value(), gradient(), hessian(),
       * or jacobian() operations can be performed.
       *
       * @note The chosen tape index must be greater than
       * numbers::invalid_tape_index and less than numbers::max_tape_index.
       */
      bool
      start_recording_operations(const types::tape_index tape_index,
                                 const bool              overwrite_tape = false,
                                 const bool              keep_values    = true);

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
       * numbers::invalid_tape_index and less than numbers::max_tape_index.
       */
      void
      activate_recorded_tape(const types::tape_index tape_index);

      //@}

    protected:
      /**
       * @name Taping
       */
      //@{

      /**
       * Index of the tape that is currently in use. It is this tape that will
       * be recorded to or read from when performing computations using "taped"
       * auto-differentiable numbers.
       */
      types::tape_index active_tape_index;

      /**
       * A collection of tapes that have been recorded to on this process.
       *
       * It is important to keep track of this so that one doesn't accidentally
       * record over a tape (unless specifically instructed to) and that one
       * doesn't try to use a tape that doesn't exist.
       */
      static std::set<types::tape_index> registered_tapes;

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
      types::tape_buffer_sizes obufsize;

      /**
       * ADOL-C locations buffer size.
       */
      types::tape_buffer_sizes lbufsize;

      /**
       * ADOL-C value buffer size.
       */
      types::tape_buffer_sizes vbufsize;

      /**
       * ADOL-C Taylor buffer size.
       */
      types::tape_buffer_sizes tbufsize;

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
       * numbers::invalid_tape_index and less than numbers::max_tape_index.
       */
      void
      activate_tape(const types::tape_index tape_index, const bool read_mode);

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

    private:
      /**
       * @name Miscellaneous
       */
      //@{

      /**
       * A counter keeping track of the number of helpers in existence.
       *
       * This is only important information for when we use taped number types.
       * As the tapes are stored in some global register, they exist independent
       * of these helpers. However, it is assumed that when all helpers go
       * out of scope then the tapes can be written over.
       */
      static unsigned int n_helpers;

      //@}

    }; // class ADHelperBase



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
     *   1. there are are many independent variables (the local
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
     * use the ADHelperBase::set_tape_buffer_sizes() function before starting
     * taping (as done via the ADHelperBase::start_recording_operations()
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
    class ADHelperCellLevelBase
      : public ADHelperBase<ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename ADHelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename ADHelperBase<ADNumberTypeCode, ScalarType>::ad_type;

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
       * $\mathbf{X}$, i.e. the dimension of the domain space.
       * @param[in] n_dependent_variables The number of scalar functions to be
       * defined that will have a sensitivity to the given independent
       * variables. In the computation of $\mathbf{f}(\mathbf{X})$, this will
       * be the number of outputs $\mathbf{f}$, i.e. the dimension of the
       * image space.
       */
      ADHelperCellLevelBase(const unsigned int n_independent_variables,
                            const unsigned int n_dependent_variables);

      /**
       * Destructor
       */
      virtual ~ADHelperCellLevelBase() = default;

      //@}

      /**
       * @name Independent variables
       */
      //@{

      /**
       * Register the complete set of independent variables $\mathbf{X}$ that
       * represent the local degree-of-freedom values.
       *
       * @param[in] dof_values A vector field associated with local
       * degree-of-freedom values on the current finite element. These define
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
       * represent the local degree-of-freedom values.
       *
       * @param[in] values A global field from which the values of all
       * independent variables will be extracted. This typically will be the
       * solution vector around which point a residual vector is to be
       * computed and around which linearization is to occur.
       * When considering taped AD numbers with branching functions, to avoid
       * potential issues with branch switching it may be a good idea to choose
       * these values close or equal to those that will be later evaluated and
       * linearized around.
       * @param[in] local_dof_indices A vector of degree-of-freedom indices from
       * which to extract the local degree-of-freedom values. This would
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
       * Returns the complete set of degree-of-freedom values as represented by
       * auto-differentiable numbers. These are the independent
       * variables $\mathbf{X}$ about which the solution is linearized.
       *
       * It is indicated to the AD library that operations performed with these
       * numbers are to be tracked, so they are considered "sensitive"
       * variables. This is, therefore, the set of variables with which one
       * would then perform computations, and based on which one can then
       * extract both the value of the function and its derivatives with the
       * member functions below. The values of the components of the returned
       * object are initialized to the values set with
       * register_independent_variable().
       *
       * @return An array of auto-differentiable type numbers representing the
       * local degree-of-freedom values.
       *
       * @note For taped AD numbers, this operation is only valid in recording mode.
       */
      const std::vector<ad_type> &
      get_sensitive_dof_values();

      //@}

      /**
       * @name Post-processing
       */
      //@{

      /*
       * Returns the complete set of degree-of-freedom values of
       * auto-differentiable number type. These store the same scalar values as
       * the independent variables $\mathbf{X}$ about which the solution is
       * linearized.
       *
       * Operations performed with these numbers are not tracked by the AD,
       * libraries so they are considered "non-sensitive" variables.
       * The values of the components of the returned object are initialized to
       * the values set with register_dof_values().
       *
       * @return An array of auto-differentiable type numbers representing the
       * local degree-of-freedom values.
       *
       * @note This function is not typically used within the context of automatic
       * differentation computations, but can make performing substutitions in
       * other post-processing computations more straight forward.
       *
       * @note For taped AD numbers, this operation is only valid outside recording mode.
       */
      std::vector<ad_type>
      get_non_sensitive_dof_values() const;

      //@}

      /**
       * @name Operations specific to taped mode: Reusing tapes
       */
      //@{

      /**
       * Set the values for the independent variables $\mathbf{X}$, i.e. the
       * linearization point.
       *
       * @param[in] dof_values A vector field associated with local
       * degree-of-freedom values on the current finite element. These define
       * the values of all independent variables.
       *
       * @note The input value type must correspond to this class's @p scalar_type.
       * Depending on the selected @p ADNumberTypeCode, this may or may not
       * correspond with the @p ScalarType prescribed as a template argument.
       *
       * @note If the keep flag has been set when
       * ADHelperBase::start_recording_operations() is called then the tape is
       * immediately usable after creation, and the values of the independent
       * variables set by register_dof_values() are those at which the function
       * is to be evaluated. In this case, a separate call to this function is
       * not strictly necessary.
       */
      void
      set_dof_values(const std::vector<scalar_type> &dof_values);

      /**
       * Set the values for the independent variables $\mathbf{X}$, i.e. the
       * linearization point.
       *
       * @param[in] values A vector field from which the values of all
       * independent variables is to be extracted.
       * @param[in] local_dof_indices A vector of degree-of-freedom indices from
       * which to extract the local degree-of-freedom values. This would
       * typically obtained by calling <code>cell->get_dof_indices()</code>.
       *
       * @note If the keep flag has been set when
       * ADHelperBase::start_recording_operations() is called then the tape is
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
       * Computes the value of the residual vector field
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
       * Computes the gradient (first derivative) of the residual vector field
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

    }; // class ADHelperCellLevelBase



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
     *   for (auto cell & : ...)
     *   {
     *     cell->get_dof_indices(local_dof_indices);
     *     const unsigned int n_independent_variables =
     *       local_dof_indices.size();
     *
     *     // Create some aliases for the AD helper.
     *     // In the example, the AD_typecode used for the template argument can
     *     // be refer to either a taped or tapeless type.
     *     using ADHelper = AD::ADHelperEnergyFunctional<...>;
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
     *     const types::tape_index tape_index = ...;
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
     *       // Then we get the complete set of degree-of-freedom values as
     *       // represented by auto-differentiable numbers. The operations
     *       // performed with these variables are tracked by the AD library
     *       // from this point until stop_recording_operations() is called.
     *       const std::vector<ADNumberType> dof_values_ad
     *         = ad_helper.get_sensitive_dof_values();
     *
     *       // Then we do some problem specific tasks, the first being to
     *       // compute all values, gradients etc. based on sensitive AD DoF
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
     *                ExcMessage("Negative jacobian detected!"));
     *
     *         // Add contribution of the internal energy:
     *         // Integrate the stored energy density function with the current
     *         // solution.
     *         energy_ad += get_Psi(F) * fe_values.JxW(q_point);
     *       }
     *
     *       // Add contribution from external energy:
     *       // Loop over faces and accumulate external energy into cell
     *       // total energy
     *       // energy_ad += ...
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
     * @warning ADOL-C does not support the standard threading models used by
     * deal.II, so this class should @b not be embedded within a multithreaded
     * function when using ADOL-C number types. It is, however, suitable for use
     * in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016, 2017, 2018
     */
    template <enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType = double>
    class ADHelperEnergyFunctional
      : public ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>
    {
    public:
      /**
       * Type definition for the floating point number type that is used in,
       * and results from, all computations.
       */
      using scalar_type =
        typename ADHelperBase<ADNumberTypeCode, ScalarType>::scalar_type;

      /**
       * Type definition for the auto-differentiation number type that is used
       * in all computations.
       */
      using ad_type =
        typename ADHelperBase<ADNumberTypeCode, ScalarType>::ad_type;

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
       * $\mathbf{X}$, i.e. the dimension of the domain space.
       *
       * @note There is only one dependent variable associated with the total
       * energy attributed to the local finite element. That is to say, this
       * class assumes that the (local) right hand side and matrix contribution
       * is computed from the first and second derivatives of a scalar
       * function $\Psi(\mathbf{X})$.
       */
      ADHelperEnergyFunctional(const unsigned int n_independent_variables);

      /**
       * Destructor
       */
      virtual ~ADHelperEnergyFunctional() = default;

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
       * degree-of-freedom values, i.e.
       * @f[
       *   \Psi(\mathbf{X}) \vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are by calling
       * ADHelperCellLevelBase::set_dof_values().
       *
       * @return The value of the energy functional at the evaluation point
       * corresponding to a chosen set of local degree-of freedom values.
       */
      scalar_type
      compute_energy() const;

      /**
       * Evaluation of the residual for a chosen set of degree-of-freedom
       * values. Underlying this is the computation of the gradient (first
       * derivative) of the scalar function $\Psi$ with respect to all
       * independent variables, i.e.
       * @f[
       *   \mathbf{r}(\mathbf{X}) =
       * \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{X}}
       * \Big\vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are by calling
       * ADHelperCellLevelBase::set_dof_values().
       *
       * @param[out] residual A Vector object, for which the value for each
       * entry represents the residual value for the corresponding local
       * degree-of freedom. The output @p residual vector has a length
       * corresponding to @p n_independent_variables.
       */
      void
      compute_residual(Vector<scalar_type> &residual) const override;

      /**
       * Computes the linearization of the residual vector around a chosen set
       * of degree-of-freedom values. Underlying this is the computation of the
       * Hessian (second derivative) of the scalar function $\Psi$ with respect
       * to all independent variables, i.e.
       * @f[
       *   \frac{\partial\mathbf{r}(\mathbf{X})}{\partial\mathbf{X}}
       *     =
       * \frac{\partial^{2}\Psi(\mathbf{X})}{\partial\mathbf{X}
       * \otimes \partial\mathbf{X}} \Big\vert_{\mathbf{X}}
       * @f]
       *
       * The values at the evaluation point $\mathbf{X}$ are by calling
       * ADHelperCellLevelBase::set_dof_values().
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

    }; // class ADHelperEnergyFunctional


  } // namespace AD
} // namespace Differentiation


/* ----------------- inline and template functions ----------------- */


#  ifndef DOXYGEN

namespace Differentiation
{
  namespace AD
  {
    /* ----------------- ADHelperCellLevelBase ----------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template <typename VectorType>
    void
    ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::register_dof_values(
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
          "Degree-of-freedom index vector size does not match number of independent variables"));
#    ifdef DEBUG
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        {
          Assert(this->registered_independent_variable_values[i] == false,
                 ExcMessage("Independent variables already registered."));
        }
#    endif
      set_dof_values(values, local_dof_indices);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template <typename VectorType>
    void
    ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::set_dof_values(
      const VectorType &                                  values,
      const std::vector<dealii::types::global_dof_index> &local_dof_indices)
    {
      Assert(local_dof_indices.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        ADHelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
          i, values[local_dof_indices[i]]);
    }


  } // namespace AD
} // namespace Differentiation


#  endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)

#endif // dealii__adolc_helpers_h
