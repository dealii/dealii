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

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)

#  include <deal.II/differentiation/ad/ad_drivers.h>
#  include <deal.II/differentiation/ad/ad_helpers.h>

#  include <type_traits>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace AD
  {
    /* -------------------------- HelperBase -------------------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    HelperBase<ADNumberTypeCode, ScalarType>::HelperBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      : independent_variable_values(
          n_independent_variables,
          dealii::internal::NumberType<scalar_type>::value(0.0))
      , registered_independent_variable_values(n_independent_variables, false)
      , registered_marked_independent_variables(n_independent_variables, false)
      , registered_marked_dependent_variables(n_dependent_variables, false)
    {
      // We have enabled the compilation of this class for arithmetic
      // types (i.e. ADNumberTypeCode == NumberTypes::none), but we
      // can't actually do anything with them. Lets not advance any further
      // and seemingly allow any operations that will not give any
      // sensible results.
      Assert(ADNumberTypeCode != NumberTypes::none,
             ExcMessage(
               "Floating point/arithmetic numbers have no derivatives."));
      Assert(
        AD::ADNumberTraits<ad_type>::n_supported_derivative_levels >= 1,
        ExcMessage(
          "The AD number type does not support the calculation of any derivatives."));

      // Tapeless mode must be configured before any active live
      // variables are created.
      if (AD::is_tapeless_ad_number<ad_type>::value)
        {
          configure_tapeless_mode(n_independent_variables,
                                  false /*ensure_persistent_setting*/);
        }

      // For safety, we ensure that the entries in this vector *really* are
      // initialized correctly by sending in the constructed zero-valued
      // initializer.
      dependent_variables.resize(n_dependent_variables,
                                 dealii::internal::NumberType<ad_type>::value(
                                   0.0));
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode,
               ScalarType>::reset_registered_independent_variables()
    {
      std::fill(registered_independent_variable_values.begin(),
                registered_independent_variable_values.end(),
                false);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::
      reset_registered_dependent_variables(const bool flag)
    {
      std::fill(registered_marked_dependent_variables.begin(),
                registered_marked_dependent_variables.end(),
                flag);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
      const unsigned int index,
      const scalar_type &value)
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        {
          // A dummy call in case the user does not encapsulate a set of
          // calls to tapeless ADHelpers with the initial and final calls to
          // [start,stop]_recording_operations.
          if (this->is_recording() == false)
            start_recording_operations(1 /*tape index*/);

          Assert(this->is_recording() == true,
                 ExcMessage(
                   "Cannot change the value of an independent variable "
                   "of the tapeless variety while this class is not set "
                   "in recording operations."));
        }
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(
        index < n_independent_variables(),
        ExcMessage(
          "Trying to set the value of a non-existent independent variable."));

      independent_variable_values[index]            = value;
      registered_independent_variable_values[index] = true;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::mark_independent_variable(
      const unsigned int index,
      ad_type &          out) const
    {
      Assert(index < n_independent_variables(), ExcInternalError());
      Assert(registered_independent_variable_values[index] == true,
             ExcInternalError());

      if (index > 0)
        {
          Assert(
            registered_marked_independent_variables[index - 1] == true,
            ExcMessage(
              "Need to extract sensitivities in the order they're created."));
        }

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(active_tape_index() != Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(is_recording() == true,
                 ExcMessage(
                   "The marking of independent variables is only valid "
                   "during recording."));
        }

      internal::Marking<ad_type>::independent_variable(
        independent_variable_values[index],
        index,
        this->n_independent_variables(),
        out);
      registered_marked_independent_variables[index] = true;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode,
               ScalarType>::finalize_sensitive_independent_variables() const
    {
      // Double check that we've actually registered all DoFs
      Assert(n_registered_independent_variables() == n_independent_variables(),
             ExcMessage("Not all values of sensitivities have been recorded!"));

      // This should happen only once
      if (this->independent_variables.size() == 0)
        {
          this->independent_variables.resize(
            this->n_independent_variables(),
            dealii::internal::NumberType<ad_type>::value(0.0));

          // Indicate the sensitivity that each entry represents
          for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
            this->mark_independent_variable(i, this->independent_variables[i]);
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::
      initialize_non_sensitive_independent_variable(const unsigned int index,
                                                    ad_type &out) const
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(active_tape_index() != Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(is_recording() == false,
             ExcMessage(
               "The initialization of non-sensitive independent variables is "
               "only valid outside of recording operations."));

      Assert(index < n_independent_variables(), ExcInternalError());
      Assert(registered_independent_variable_values[index] == true,
             ExcInternalError());

      out = independent_variable_values[index];
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int
    HelperBase<ADNumberTypeCode,
               ScalarType>::n_registered_independent_variables() const
    {
      return std::count(registered_independent_variable_values.begin(),
                        registered_independent_variable_values.end(),
                        true);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::size_t
    HelperBase<ADNumberTypeCode, ScalarType>::n_independent_variables() const
    {
      return independent_variable_values.size();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int
    HelperBase<ADNumberTypeCode, ScalarType>::n_registered_dependent_variables()
      const
    {
      return std::count(registered_marked_dependent_variables.begin(),
                        registered_marked_dependent_variables.end(),
                        true);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::size_t
    HelperBase<ADNumberTypeCode, ScalarType>::n_dependent_variables() const
    {
      return dependent_variables.size();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    HelperBase<ADNumberTypeCode, ScalarType>::is_recording() const
    {
      if (AD::is_taped_ad_number<ad_type>::value)
        return taped_driver.is_recording();
      else
        return tapeless_driver.is_dependent_variable_marking_allowed();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    typename Types<
      typename HelperBase<ADNumberTypeCode, ScalarType>::ad_type>::tape_index
    HelperBase<ADNumberTypeCode, ScalarType>::active_tape_index() const
    {
      if (AD::is_taped_ad_number<ad_type>::value)
        return taped_driver.active_tape_index();
      else
        return 1;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    HelperBase<ADNumberTypeCode, ScalarType>::is_registered_tape(
      const typename Types<ad_type>::tape_index tape_index) const
    {
      if (AD::is_taped_ad_number<ad_type>::value)
        return taped_driver.is_registered_tape(tape_index);
      else
        return true;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::print(std::ostream &stream) const
    {
      // Store stream flags
      const std::ios_base::fmtflags stream_flags(stream.flags());
      // Set stream to print booleans as "true"/"false"
      stream.setf(std::ios_base::boolalpha);

      stream << "Active tape index: " << active_tape_index() << "\n";
      stream << "Recording? " << is_recording() << "\n";
      stream << std::flush;

      if (AD::is_taped_ad_number<ad_type>::value)
        taped_driver.print(stream);

      stream << "Registered independent variables: "
             << "\n";
      for (unsigned int i = 0; i < n_independent_variables(); i++)
        stream << registered_independent_variable_values[i]
               << (i < (n_independent_variables() - 1) ? "," : "");
      stream << std::endl;

      stream << "Independent variable values: "
             << "\n";
      print_values(stream);

      stream << "Registered marked independent variables: "
             << "\n";
      for (unsigned int i = 0; i < n_independent_variables(); i++)
        stream << registered_marked_independent_variables[i]
               << (i < (n_independent_variables() - 1) ? "," : "")
               << std::flush;
      stream << std::endl;

      stream << "Dependent variable values: "
             << "\n";
      for (unsigned int i = 0; i < n_dependent_variables(); i++)
        stream << dependent_variables[i]
               << (i < (n_dependent_variables() - 1) ? "," : "");
      stream << std::endl;

      stream << "Registered dependent variables: "
             << "\n";
      for (unsigned int i = 0; i < n_dependent_variables(); i++)
        stream << registered_marked_dependent_variables[i]
               << (i < (n_dependent_variables() - 1) ? "," : "");
      stream << std::endl;

      // Restore stream flags
      stream.flags(stream_flags);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::print_values(
      std::ostream &stream) const
    {
      for (unsigned int i = 0; i < n_independent_variables(); i++)
        stream << independent_variable_values[i]
               << (i < (n_independent_variables() - 1) ? "," : "")
               << std::flush;
      stream << std::endl;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::print_tape_stats(
      const typename Types<ad_type>::tape_index tape_index,
      std::ostream &                            stream) const
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        return;

      Assert(is_registered_tape(tape_index),
             ExcMessage("Tape number not registered"));

      this->taped_driver.print_tape_stats(tape_index, stream);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::reset(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables,
      const bool         clear_registered_tapes)
    {
      const unsigned int new_n_independent_variables =
        (n_independent_variables != dealii::numbers::invalid_unsigned_int ?
           n_independent_variables :
           this->n_independent_variables());
      const unsigned int new_n_dependent_variables =
        (n_dependent_variables != dealii::numbers::invalid_unsigned_int ?
           n_dependent_variables :
           this->n_dependent_variables());

      // Here we clear our vectors of AD data with a reallocation of memory
      // (i.e. we *nuke* them entirely). Why we do this differs for each
      // AD type:
      // - ADOL-C taped mode must have their tapes fully cleared of marked data
      //   before the tapes can be overwritten.
      // - ADOL-C tapeless mode must be configured for before any active live
      //   variables are created.
      // - Reverse-mode Sacado numbers to must be destroyed to reset their
      //   accumulations.
      // - Forward-mode Sacado numbers have no specific requirements, but it
      //   doesn't really hurt to perform this operation anyway.
      {
        std::vector<ad_type>().swap(independent_variables);
        std::vector<ad_type>().swap(dependent_variables);
      }

      // Tapeless mode must be configured before any active live
      // variables are created.
      if (AD::is_tapeless_ad_number<ad_type>::value)
        {
          configure_tapeless_mode(new_n_independent_variables,
                                  false /*ensure_persistent_setting*/);
        }
      if (AD::is_taped_ad_number<ad_type>::value)
        taped_driver.reset(clear_registered_tapes);

      independent_variable_values = std::vector<scalar_type>(
        new_n_independent_variables,
        dealii::internal::NumberType<scalar_type>::value(0.0));
      registered_independent_variable_values =
        std::vector<bool>(new_n_independent_variables, false);
      registered_marked_independent_variables =
        std::vector<bool>(new_n_independent_variables, false);
      dependent_variables =
        std::vector<ad_type>(new_n_dependent_variables,
                             dealii::internal::NumberType<ad_type>::value(0.0));
      registered_marked_dependent_variables =
        std::vector<bool>(new_n_dependent_variables, false);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::configure_tapeless_mode(
      const unsigned int n_independent_variables,
      const bool         ensure_persistent_setting)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        return;

      // Try to safely initialize the global environment
      TapelessDrivers<ad_type, scalar_type>::initialize_global_environment(
        n_independent_variables);

      if (ensure_persistent_setting == true)
        if (is_adolc_tapeless_number<ad_type>::value)
          {
            // In order to ensure that the settings remain for the entire
            // duration of the simulation, we create a global live variable
            // that doesn't go out of scope.
            static ad_type num = 0.0;
            (void)num;
          }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::activate_recorded_tape(
      const typename Types<ad_type>::tape_index tape_index)
    {
      activate_tape(tape_index, true /*read_mode*/);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    HelperBase<ADNumberTypeCode, ScalarType>::recorded_tape_requires_retaping(
      const typename Types<ad_type>::tape_index tape_index) const
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        return false;

      return taped_driver.requires_retaping(tape_index);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    HelperBase<ADNumberTypeCode, ScalarType>::active_tape_requires_retaping()
      const
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        return false;

      return taped_driver.last_action_requires_retaping();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::clear_active_tape()
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        return;

      taped_driver.remove_tape(taped_driver.active_tape_index());
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::activate_tape(
      const typename Types<ad_type>::tape_index tape_index,
      const bool                                read_mode)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(tape_index != Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(tape_index < Numbers<ad_type>::max_tape_index,
                 ExcMessage("Tape index exceeds maximum allowable value"));
          taped_driver.activate_tape(tape_index);
          reset_registered_independent_variables();

          // A tape may have been defined by a different ADHelper, so in this
          // case we ignore the fact that any dependent variables within the
          // current data structure have not been marked as dependents
          if (read_mode == true)
            {
              Assert(is_registered_tape(tape_index),
                     ExcMessage("Tape number not registered"));
              reset_registered_dependent_variables(true);
              Assert(n_registered_dependent_variables() ==
                       n_dependent_variables(),
                     ExcMessage("Not all dependent variables have been set!"));
            }
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::set_tape_buffer_sizes(
      const typename Types<ad_type>::tape_buffer_sizes obufsize,
      const typename Types<ad_type>::tape_buffer_sizes lbufsize,
      const typename Types<ad_type>::tape_buffer_sizes vbufsize,
      const typename Types<ad_type>::tape_buffer_sizes tbufsize)
    {
      // When valid for the chosen AD number type, these values will be used the
      // next time start_recording_operations() is called.
      if (ADNumberTraits<ad_type>::is_taped == true)
        taped_driver.set_tape_buffer_sizes(obufsize,
                                           lbufsize,
                                           vbufsize,
                                           tbufsize);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    HelperBase<ADNumberTypeCode, ScalarType>::start_recording_operations(
      const typename Types<ad_type>::tape_index tape_index,
      const bool                                overwrite_tape,
      const bool                                keep_independent_values)
    {
      // Define this here for clarity when this flag is used later.
      const bool read_mode = false;

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          if (overwrite_tape != true)
            {
              Assert(is_recording() == false,
                     ExcMessage("Already recording..."));
            }

          // Check conditions to enable tracing
          if (is_registered_tape(tape_index) == false || overwrite_tape == true)
            {
              // Setup the data structures for this class in the
              // appropriate manner
              activate_tape(tape_index, read_mode);

              // Start taping
              taped_driver.start_taping(active_tape_index(),
                                        keep_independent_values);

              // Clear the flags that state which independent and
              // dependent variables have been registered
              reset_registered_independent_variables();
              reset_registered_dependent_variables();
            }
          else
            {
              Assert(is_recording() == false,
                     ExcMessage(
                       "Tape recording is unexpectedly still enabled."));

              // Now we activate the pre-recorded tape so that its immediately
              // available for use
              activate_recorded_tape(tape_index);
            }
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());

          // Set the flag that states that we can safely mark dependent
          // variables within this current phase of operations
          tapeless_driver.allow_dependent_variable_marking();

          // Dummy call to ensure that the intuitively correct
          // value for the active tape (whether "valid" or not)
          // is always returned to the user.
          activate_tape(tape_index, read_mode);
        }

      return is_recording();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::stop_recording_operations(
      const bool write_tapes_to_file)
    {
      Assert(is_recording() == true, ExcMessage("Not currently recording..."));

      // Double check that we've actually registered all DoFs
      Assert(n_registered_independent_variables() == n_independent_variables(),
             ExcMessage("Not all values of sensitivities have been recorded!"));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          // Stop tracing
          taped_driver.stop_taping(active_tape_index(), write_tapes_to_file);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          // Double check that we've actually registered dependent variables
          Assert(n_registered_dependent_variables() == n_dependent_variables(),
                 ExcMessage("Not all dependent variables have been set!"));

          // By changing this flag, we ensure that the we can no longer
          // legally alter the values of the dependent variables using
          // set_dependent_variable(). This is important because the value of
          // the tapeless independent variables are set and finalized when
          // mark_dependent_variable() is called. So we cannot allow this to
          // be done when not in the "recording" phase.
          tapeless_driver.prevent_dependent_variable_marking();
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    HelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
      const unsigned int index,
      const ad_type &    func)
    {
      Assert(index < n_dependent_variables(), ExcMessage("Index out of range"));
      Assert(registered_marked_dependent_variables[index] == false,
             ExcMessage(
               "This dependent variable has already been registered."));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(active_tape_index() != Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(is_recording() == true,
                 ExcMessage(
                   "Must be recording when registering dependent variables."));
        }

      // Register the given dependent variable
      internal::Marking<ad_type>::dependent_variable(dependent_variables[index],
                                                     func);
      registered_marked_dependent_variables[index] = true;
    }



    /* -------------------- CellLevelBase -------------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    CellLevelBase<ADNumberTypeCode, ScalarType>::CellLevelBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      : HelperBase<ADNumberTypeCode, ScalarType>(n_independent_variables,
                                                 n_dependent_variables)
    {}



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    CellLevelBase<ADNumberTypeCode, ScalarType>::register_dof_values(
      const std::vector<scalar_type> &dof_values)
    {
      // This is actually the same thing the set_independent_variable function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
      Assert(dof_values.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        {
          Assert(this->registered_independent_variable_values[i] == false,
                 ExcMessage("Independent variable value already registered."));
        }
      set_dof_values(dof_values);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    const std::vector<
      typename CellLevelBase<ADNumberTypeCode, ScalarType>::ad_type> &
    CellLevelBase<ADNumberTypeCode, ScalarType>::get_sensitive_dof_values()
      const
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // If necessary, initialize the internally stored vector of
      // AD numbers that represents the independent variables
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size() ==
               this->n_independent_variables(),
             ExcDimensionMismatch(this->independent_variables.size(),
                                  this->n_independent_variables()));

      return this->independent_variables;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    CellLevelBase<ADNumberTypeCode, ScalarType>::set_dof_values(
      const std::vector<scalar_type> &values)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        HelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
          i, values[i]);
    }



    /* ------------------ EnergyFunctional ------------------ */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    EnergyFunctional<ADNumberTypeCode, ScalarType>::EnergyFunctional(
      const unsigned int n_independent_variables)
      : CellLevelBase<ADNumberTypeCode, ScalarType>(n_independent_variables, 1)
    {}



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    EnergyFunctional<ADNumberTypeCode, ScalarType>::register_energy_functional(
      const ad_type &energy)
    {
      Assert(this->n_dependent_variables() == 1, ExcInternalError());
      HelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
        0, energy);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    typename EnergyFunctional<ADNumberTypeCode, ScalarType>::scalar_type
    EnergyFunctional<ADNumberTypeCode, ScalarType>::compute_energy() const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(
        this->n_dependent_variables() == 1,
        ExcMessage(
          "The EnergyFunctional class expects there to be only one dependent variable."));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute value while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          return this->taped_driver.value(this->active_tape_index(),
                                          this->independent_variable_values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          return this->tapeless_driver.value(this->dependent_variables);
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    EnergyFunctional<ADNumberTypeCode, ScalarType>::compute_residual(
      Vector<scalar_type> &gradient) const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(
        this->n_dependent_variables() == 1,
        ExcMessage(
          "The EnergyFunctional class expects there to be only one dependent variable."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::gradient().
      if (gradient.size() != this->n_independent_variables())
        gradient.reinit(this->n_independent_variables(),
                        true /*omit_zeroing_entries*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute gradient while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.gradient(this->active_tape_index(),
                                      this->independent_variable_values,
                                      gradient);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          this->tapeless_driver.gradient(this->independent_variables,
                                         this->dependent_variables,
                                         gradient);
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    EnergyFunctional<ADNumberTypeCode, ScalarType>::compute_linearization(
      FullMatrix<scalar_type> &hessian) const
    {
      Assert(AD::ADNumberTraits<ad_type>::n_supported_derivative_levels >= 2,
             ExcMessage(
               "Cannot computed function Hessian: AD number type does "
               "not support the calculation of second order derivatives."));

      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false))
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(
        this->n_dependent_variables() == 1,
        ExcMessage(
          "The EnergyFunctional class expects there to be only one dependent variable."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::hessian().
      if (hessian.m() != this->n_independent_variables() ||
          hessian.n() != this->n_independent_variables())
        hessian.reinit({this->n_independent_variables(),
                        this->n_independent_variables()},
                       true /*omit_default_initialization*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute hessian while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.hessian(this->active_tape_index(),
                                     this->independent_variable_values,
                                     hessian);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          this->tapeless_driver.hessian(this->independent_variables,
                                        this->dependent_variables,
                                        hessian);
        }
    }


    /* ------------------- ResidualLinearization ------------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ResidualLinearization<ADNumberTypeCode, ScalarType>::ResidualLinearization(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      : CellLevelBase<ADNumberTypeCode, ScalarType>(n_independent_variables,
                                                    n_dependent_variables)
    {}



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ResidualLinearization<ADNumberTypeCode, ScalarType>::
      register_residual_vector(const std::vector<ad_type> &residual)
    {
      Assert(residual.size() == this->n_dependent_variables(),
             ExcMessage(
               "Vector size does not match number of dependent variables"));
      for (unsigned int i = 0; i < this->n_dependent_variables(); ++i)
        HelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
          i, residual[i]);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ResidualLinearization<ADNumberTypeCode, ScalarType>::compute_residual(
      Vector<scalar_type> &values) const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::values().
      if (values.size() != this->n_dependent_variables())
        values.reinit(this->n_dependent_variables(),
                      true /*omit_zeroing_entries*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute values while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.values(this->active_tape_index(),
                                    this->n_dependent_variables(),
                                    this->independent_variable_values,
                                    values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          this->tapeless_driver.values(this->dependent_variables, values);
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ResidualLinearization<ADNumberTypeCode, ScalarType>::compute_linearization(
      FullMatrix<scalar_type> &jacobian) const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::jacobian().
      if (jacobian.m() != this->n_dependent_variables() ||
          jacobian.n() != this->n_independent_variables())
        jacobian.reinit({this->n_dependent_variables(),
                         this->n_independent_variables()},
                        true /*omit_default_initialization*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute hessian while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.jacobian(this->active_tape_index(),
                                      this->n_dependent_variables(),
                                      this->independent_variable_values,
                                      jacobian);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          this->tapeless_driver.jacobian(this->independent_variables,
                                         this->dependent_variables,
                                         jacobian);
        }
    }



    /* ----------------- PointLevelFunctionsBase  ----------------- */



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      PointLevelFunctionsBase(const unsigned int n_independent_variables,
                              const unsigned int n_dependent_variables)
      : HelperBase<ADNumberTypeCode, ScalarType>(n_independent_variables,
                                                 n_dependent_variables)
      , symmetric_independent_variables(n_independent_variables, false)
    {}



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::reset(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables,
      const bool         clear_registered_tapes)
    {
      HelperBase<ADNumberTypeCode, ScalarType>::reset(n_independent_variables,
                                                      n_dependent_variables,
                                                      clear_registered_tapes);

      const unsigned int new_n_independent_variables =
        (n_independent_variables != dealii::numbers::invalid_unsigned_int ?
           n_independent_variables :
           this->n_independent_variables());
      symmetric_independent_variables =
        std::vector<bool>(new_n_independent_variables, false);
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    bool
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      is_symmetric_independent_variable(const unsigned int index) const
    {
      Assert(index < symmetric_independent_variables.size(),
             ExcInternalError());
      return symmetric_independent_variables[index];
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    unsigned int
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      n_symmetric_independent_variables() const
    {
      return std::count(symmetric_independent_variables.begin(),
                        symmetric_independent_variables.end(),
                        true);
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      register_independent_variables(const std::vector<scalar_type> &values)
    {
      // This is actually the same thing the set_independent_variable function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        {
          Assert(this->registered_independent_variable_values[i] == false,
                 ExcMessage("Independent variable value already registered."));
        }
      set_independent_variables(values);
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    const std::vector<typename PointLevelFunctionsBase<dim,
                                                       ADNumberTypeCode,
                                                       ScalarType>::ad_type> &
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      get_sensitive_variables() const
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // Just in case the user has not done so, we repeat the call to
      // initialize the internally stored vector of AD numbers that
      // represents the independent variables.
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size() ==
               this->n_independent_variables(),
             ExcDimensionMismatch(this->independent_variables.size(),
                                  this->n_independent_variables()));

      return this->independent_variables;
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      set_sensitivity_value(const unsigned int index,
                            const bool         symmetric_component,
                            const scalar_type &value)
    {
      HelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(index,
                                                                      value);
      Assert(
        index < this->n_independent_variables(),
        ExcMessage(
          "Trying to set the symmetry flag of a non-existent independent variable."));
      Assert(index < symmetric_independent_variables.size(),
             ExcInternalError());
      symmetric_independent_variables[index] = symmetric_component;
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>::
      set_independent_variables(const std::vector<scalar_type> &values)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        HelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
          i, values[i]);
    }



    /* -------------------- ScalarFunction -------------------- */



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::ScalarFunction(
      const unsigned int n_independent_variables)
      : PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>(
          n_independent_variables,
          1)
    {}



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::
      register_dependent_variable(const ad_type &func)
    {
      Assert(this->n_dependent_variables() == 1, ExcInternalError());
      HelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
        0, func);
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    typename ScalarFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::compute_value() const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(
        this->n_dependent_variables() == 1,
        ExcMessage(
          "The ScalarFunction class expects there to be only one dependent variable."));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute values while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          return this->taped_driver.value(this->active_tape_index(),
                                          this->independent_variable_values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          return this->tapeless_driver.value(this->dependent_variables);
        }
    }


    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::compute_gradient(
      Vector<scalar_type> &gradient) const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(
        this->n_dependent_variables() == 1,
        ExcMessage(
          "The ScalarFunction class expects there to be only one dependent variable."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::gradient().
      if (gradient.size() != this->n_independent_variables())
        gradient.reinit(this->n_independent_variables(),
                        true /*omit_zeroing_entries*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute gradient while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.gradient(this->active_tape_index(),
                                      this->independent_variable_values,
                                      gradient);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          this->tapeless_driver.gradient(this->independent_variables,
                                         this->dependent_variables,
                                         gradient);
        }

      // Account for symmetries of tensor components
      for (unsigned int i = 0; i < this->n_independent_variables(); i++)
        {
          if (this->is_symmetric_independent_variable(i) == true)
            gradient[i] *= 0.5;
        }
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::compute_hessian(
      FullMatrix<scalar_type> &hessian) const
    {
      Assert(AD::ADNumberTraits<ad_type>::n_supported_derivative_levels >= 2,
             ExcMessage(
               "Cannot computed function Hessian: AD number type does "
               "not support the calculation of second order derivatives."));

      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false))
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(
        this->n_dependent_variables() == 1,
        ExcMessage(
          "The ScalarFunction class expects there to be only one dependent variable."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::hessian().
      if (hessian.m() != this->n_independent_variables() ||
          hessian.n() != this->n_independent_variables())
        hessian.reinit({this->n_independent_variables(),
                        this->n_independent_variables()},
                       true /*omit_default_initialization*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute Hessian while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.hessian(this->active_tape_index(),
                                     this->independent_variable_values,
                                     hessian);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          this->tapeless_driver.hessian(this->independent_variables,
                                        this->dependent_variables,
                                        hessian);
        }

      // Account for symmetries of tensor components
      for (unsigned int i = 0; i < this->n_independent_variables(); i++)
        for (unsigned int j = 0; j < i + 1; j++)
          {
            if (this->is_symmetric_independent_variable(i) == true &&
                this->is_symmetric_independent_variable(j) == true)
              {
                hessian[i][j] *= 0.25;
                if (i != j)
                  hessian[j][i] *= 0.25;
              }
            else if ((this->is_symmetric_independent_variable(i) == true &&
                      this->is_symmetric_independent_variable(j) == false) ||
                     (this->is_symmetric_independent_variable(j) == true &&
                      this->is_symmetric_independent_variable(i) == false))
              {
                hessian[i][j] *= 0.5;
                if (i != j)
                  hessian[j][i] *= 0.5;
              }
          }
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    Tensor<
      0,
      dim,
      typename ScalarFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type>
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_hessian_component(const FullMatrix<scalar_type> &   hessian,
                                const FEValuesExtractors::Scalar &extractor_row,
                                const FEValuesExtractors::Scalar &extractor_col)
    {
      // NOTE: It is necessary to make special provision for the case when the
      // HessianType is scalar. Unfortunately Tensor<0,dim> does not provide
      // the function unrolled_to_component_indices!
      // NOTE: The order of components must be consistently defined throughout
      // this class.
      Tensor<0, dim, scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the
      // matrix values
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set(
        internal::extract_field_component_indices<dim>(extractor_col));
      Assert(row_index_set.size() == 1, ExcInternalError());
      Assert(col_index_set.size() == 1, ExcInternalError());

      internal::set_tensor_entry(out,
                                 0,
                                 hessian[row_index_set[0]][col_index_set[0]]);

      return out;
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    SymmetricTensor<
      4,
      dim,
      typename ScalarFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type>
    ScalarFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_hessian_component(
        const FullMatrix<scalar_type> &               hessian,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_row,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_col)
    {
      // NOTE: The order of components must be consistently defined throughout
      // this class.
      // NOTE: We require a specialisation for rank-4 symmetric tensors because
      // they do not define their rank, and setting data using TableIndices is
      // somewhat specialised as well.
      SymmetricTensor<4, dim, scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the
      // matrix values
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set(
        internal::extract_field_component_indices<dim>(extractor_col));

      for (unsigned int r = 0; r < row_index_set.size(); ++r)
        for (unsigned int c = 0; c < col_index_set.size(); ++c)
          {
            internal::set_tensor_entry(
              out, r, c, hessian[row_index_set[r]][col_index_set[c]]);
          }

      return out;
    }



    /* -------------------- VectorFunction -------------------- */



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::VectorFunction(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      : PointLevelFunctionsBase<dim, ADNumberTypeCode, ScalarType>(
          n_independent_variables,
          n_dependent_variables)
    {}



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::
      register_dependent_variables(const std::vector<ad_type> &funcs)
    {
      Assert(funcs.size() == this->n_dependent_variables(),
             ExcMessage(
               "Vector size does not match number of dependent variables"));
      for (unsigned int i = 0; i < this->n_dependent_variables(); ++i)
        HelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
          i, funcs[i]);
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::compute_values(
      Vector<scalar_type> &values) const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::values().
      if (values.size() != this->n_dependent_variables())
        values.reinit(this->n_dependent_variables(),
                      true /*omit_zeroing_entries*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute values while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.values(this->active_tape_index(),
                                    this->n_dependent_variables(),
                                    this->independent_variable_values,
                                    values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          this->tapeless_driver.values(this->dependent_variables, values);
        }
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    void
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::compute_jacobian(
      FullMatrix<scalar_type> &jacobian) const
    {
      if ((ADNumberTraits<ad_type>::is_taped == true &&
           this->taped_driver.keep_independent_values() == false) ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(
            this->n_registered_independent_variables() ==
              this->n_independent_variables(),
            ExcMessage(
              "Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_registered_dependent_variables() ==
               this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately in the succeeding call to
      // Drivers::jacobian().
      if (jacobian.m() != this->n_dependent_variables() ||
          jacobian.n() != this->n_independent_variables())
        jacobian.reinit({this->n_dependent_variables(),
                         this->n_independent_variables()},
                        true /*omit_default_initialization*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape_index() !=
                   Numbers<ad_type>::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute Jacobian while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          this->taped_driver.jacobian(this->active_tape_index(),
                                      this->n_dependent_variables(),
                                      this->independent_variable_values,
                                      jacobian);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variables.size(),
                                      this->n_independent_variables()));

          this->tapeless_driver.jacobian(this->independent_variables,
                                         this->dependent_variables,
                                         jacobian);
        }

      for (unsigned int j = 0; j < this->n_independent_variables(); j++)
        {
          // Because we perform just a single differentiation
          // operation with respect to the "column" variables,
          // we only need to consider them for symmetry conditions.
          if (this->is_symmetric_independent_variable(j) == true)
            for (unsigned int i = 0; i < this->n_dependent_variables(); i++)
              jacobian[i][j] *= 0.5;
        }
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    Tensor<
      0,
      dim,
      typename VectorFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type>
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_jacobian_component(
        const FullMatrix<scalar_type> &   jacobian,
        const FEValuesExtractors::Scalar &extractor_row,
        const FEValuesExtractors::Scalar &extractor_col)
    {
      // NOTE: It is necessary to make special provision for the case when the
      // HessianType is scalar. Unfortunately Tensor<0,dim> does not provide
      // the function unrolled_to_component_indices!
      // NOTE: The order of components must be consistently defined throughout
      // this class.
      Tensor<0, dim, scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the
      // matrix values
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set(
        internal::extract_field_component_indices<dim>(extractor_col));
      Assert(row_index_set.size() == 1, ExcInternalError());
      Assert(col_index_set.size() == 1, ExcInternalError());

      internal::set_tensor_entry(out,
                                 0,
                                 jacobian[row_index_set[0]][col_index_set[0]]);

      return out;
    }



    template <int                  dim,
              enum AD::NumberTypes ADNumberTypeCode,
              typename ScalarType>
    SymmetricTensor<
      4,
      dim,
      typename VectorFunction<dim, ADNumberTypeCode, ScalarType>::scalar_type>
    VectorFunction<dim, ADNumberTypeCode, ScalarType>::
      extract_jacobian_component(
        const FullMatrix<scalar_type> &               jacobian,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_row,
        const FEValuesExtractors::SymmetricTensor<2> &extractor_col)
    {
      // NOTE: The order of components must be consistently defined throughout
      // this class.
      // NOTE: We require a specialisation for rank-4 symmetric tensors because
      // they do not define their rank, and setting data using TableIndices is
      // somewhat specialised as well.
      SymmetricTensor<4, dim, scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the
      // matrix values
      const std::vector<unsigned int> row_index_set(
        internal::extract_field_component_indices<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set(
        internal::extract_field_component_indices<dim>(extractor_col));

      for (unsigned int r = 0; r < row_index_set.size(); ++r)
        for (unsigned int c = 0; c < col_index_set.size(); ++c)
          {
            internal::set_tensor_entry(
              out, r, c, jacobian[row_index_set[r]][col_index_set[c]]);
          }

      return out;
    }


  } // namespace AD
} // namespace Differentiation


/* --- Explicit instantiations --- */
#  include "ad_helpers.inst"

#  ifdef DEAL_II_WITH_ADOLC
#    include "ad_helpers.inst1"
#  endif
#  ifdef DEAL_II_TRILINOS_WITH_SACADO
#    include "ad_helpers.inst2"
#  endif


DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)
