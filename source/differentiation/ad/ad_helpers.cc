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
    /* -------------------------- ADHelperBase -------------------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::set<types::tape_index>
      ADHelperBase<ADNumberTypeCode, ScalarType>::registered_tapes;



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int ADHelperBase<ADNumberTypeCode, ScalarType>::n_helpers = 0;



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperBase<ADNumberTypeCode, ScalarType>::ADHelperBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      : active_tape_index(numbers::invalid_tape_index)
      , keep_values(true)
      , is_recording_flag(false)
      , use_stored_taped_buffer_sizes(false)
      , obufsize(0u)
      , lbufsize(0u)
      , vbufsize(0u)
      , tbufsize(0u)
      , independent_variable_values(
          n_independent_variables,
          dealii::internal::NumberType<scalar_type>::value(0.0))
      , registered_independent_variable_values(n_independent_variables, false)
      , registered_marked_independent_variables(n_independent_variables, false)
      , registered_marked_dependent_variables(n_dependent_variables, false)
    {
      ++n_helpers;

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
    ADHelperBase<ADNumberTypeCode, ScalarType>::~ADHelperBase()
    {
      --n_helpers;

      // Clear any static data when there are no more helpers in scope.
      // This means that we've effectively marked any of these tapes to
      // be safely overwritten.
      if (n_helpers == 0)
        registered_tapes.clear();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode,
                 ScalarType>::reset_registered_independent_variables()
    {
      for (typename std::vector<bool>::iterator it =
             registered_independent_variable_values.begin();
           it != registered_independent_variable_values.end();
           ++it)
        *it = false;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::
      reset_registered_dependent_variables(const bool flag)
    {
      for (typename std::vector<bool>::iterator it =
             registered_marked_dependent_variables.begin();
           it != registered_marked_dependent_variables.end();
           ++it)
        *it = flag;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
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
                   "of the tapeless variety while this class is not set"
                   "in recording operations."));
        }
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
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
    ADHelperBase<ADNumberTypeCode, ScalarType>::mark_independent_variable(
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
          Assert(active_tape() != numbers::invalid_tape_index,
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
    ADHelperBase<ADNumberTypeCode,
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
    ADHelperBase<ADNumberTypeCode, ScalarType>::
      initialize_non_sensitive_independent_variable(const unsigned int index,
                                                    ad_type &out) const
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(active_tape() != numbers::invalid_tape_index,
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
    ADHelperBase<ADNumberTypeCode,
                 ScalarType>::n_registered_independent_variables() const
    {
      return std::count(registered_independent_variable_values.begin(),
                        registered_independent_variable_values.end(),
                        true);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::size_t
    ADHelperBase<ADNumberTypeCode, ScalarType>::n_independent_variables() const
    {
      return independent_variable_values.size();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int
    ADHelperBase<ADNumberTypeCode,
                 ScalarType>::n_registered_dependent_variables() const
    {
      return std::count(registered_marked_dependent_variables.begin(),
                        registered_marked_dependent_variables.end(),
                        true);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::size_t
    ADHelperBase<ADNumberTypeCode, ScalarType>::n_dependent_variables() const
    {
      return dependent_variables.size();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    ADHelperBase<ADNumberTypeCode, ScalarType>::is_recording() const
    {
      return is_recording_flag;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    types::tape_index
    ADHelperBase<ADNumberTypeCode, ScalarType>::active_tape() const
    {
      return active_tape_index;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    ADHelperBase<ADNumberTypeCode, ScalarType>::is_registered_tape(
      const types::tape_index tape_index) const
    {
      return registered_tapes.find(tape_index) != registered_tapes.end();
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::print(
      std::ostream &stream) const
    {
      // Store stream flags
      const std::ios_base::fmtflags stream_flags(stream.flags());
      // Set stream to print booleans as "true"/"false"
      stream.setf(std::ios_base::boolalpha);

      stream << "Active tape index: " << active_tape_index << "\n";
      stream << "Registered tapes: ";
      auto it_registered_tape = registered_tapes.begin();
      for (unsigned int i = 0; i < registered_tapes.size();
           ++i, ++it_registered_tape)
        stream << *it_registered_tape
               << (i < (registered_tapes.size() - 1) ? "," : "");
      stream << "\n";
      stream << "Recording? " << is_recording() << "\n";
      stream << "Keep values? " << keep_values << "\n";
      stream << "Use stored tape buffer sizes? "
             << use_stored_taped_buffer_sizes << "\n";
      stream << std::flush;

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
    ADHelperBase<ADNumberTypeCode, ScalarType>::print_values(
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
    ADHelperBase<ADNumberTypeCode, ScalarType>::print_tape_stats(
      const types::tape_index tape_index,
      std::ostream &          stream) const
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        return;

      Assert(is_registered_tape(tape_index),
             ExcMessage("Tape number not registered"));

      TapedDrivers<ad_type, scalar_type>::print_tape_stats(tape_index, stream);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::reset(
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
          configure_tapeless_mode(n_independent_variables,
                                  false /*ensure_persistent_setting*/);
        }

      active_tape_index = numbers::invalid_tape_index;
      is_recording_flag = false;
      if (clear_registered_tapes)
        registered_tapes.clear();

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
    ADHelperBase<ADNumberTypeCode, ScalarType>::configure_tapeless_mode(
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
          }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::activate_recorded_tape(
      const types::tape_index tape_index)
    {
      activate_tape(tape_index, true /*read_mode*/);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::activate_tape(
      const types::tape_index tape_index,
      const bool              read_mode)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(tape_index != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(tape_index < numbers::max_tape_index,
                 ExcMessage("Tape index exceeds maximum allowable value"));
          Assert(is_registered_tape(tape_index),
                 ExcMessage("Tape number not registered"));
          active_tape_index = tape_index;
          reset_registered_independent_variables();

          // A tape may have been defined by a different ADHelper, so in this
          // case we ignore the fact that any dependent variables within the
          // current data structure have not been marked as dependents
          if (read_mode == true)
            {
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
          active_tape_index = tape_index; // This is, in effect, a dummy value
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::set_tape_buffer_sizes(
      const types::tape_buffer_sizes in_obufsize,
      const types::tape_buffer_sizes in_lbufsize,
      const types::tape_buffer_sizes in_vbufsize,
      const types::tape_buffer_sizes in_tbufsize)
    {
      // When valid for the chosen AD number type, these values will be used the
      // next time start_recording_operations() is called.
      obufsize                      = in_obufsize;
      lbufsize                      = in_lbufsize;
      vbufsize                      = in_vbufsize;
      tbufsize                      = in_tbufsize;
      use_stored_taped_buffer_sizes = true;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    ADHelperBase<ADNumberTypeCode, ScalarType>::start_recording_operations(
      const types::tape_index tape_index,
      const bool              overwrite_tape,
      const bool              keep)
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
          if (is_registered_tape(tape_index) == false || overwrite_tape == true)
            {
              // Register the tape...
              registered_tapes.insert(tape_index);
              // ... and setup the data structures for this class in the
              // appropriate manner
              activate_tape(tape_index, read_mode);

              // Start taping
              if (use_stored_taped_buffer_sizes)
                {
                  Assert(obufsize > 0u,
                         ExcMessage("Buffer value not initialized."));
                  Assert(lbufsize > 0u,
                         ExcMessage("Buffer value not initialized."));
                  Assert(vbufsize > 0u,
                         ExcMessage("Buffer value not initialized."));
                  Assert(tbufsize > 0u,
                         ExcMessage("Buffer value not initialized."));
                  TapedDrivers<ad_type, scalar_type>::start_taping(
                    active_tape(),
                    keep,
                    obufsize,
                    lbufsize,
                    vbufsize,
                    tbufsize);

                  // Reset this, as we don't assume that the same tape buffer
                  // sizes are used for every tape
                  use_stored_taped_buffer_sizes = false;
                  obufsize                      = 0u;
                  lbufsize                      = 0u;
                  vbufsize                      = 0u;
                  tbufsize                      = 0u;
                }
              else
                {
                  TapedDrivers<ad_type, scalar_type>::start_taping(
                    active_tape(), keep);
                }

              // Clear the flags that state which independent and
              // dependent variables have been registered
              reset_registered_independent_variables();
              reset_registered_dependent_variables();

              // Set some other flags to their indicated / required values
              keep_values       = keep;
              is_recording_flag = true;
            }
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          is_recording_flag = true;

          // Dummy call to ensure that the intuitively correct
          // value for the active tape (whether "valid" or not)
          // is always returned to the user.
          activate_tape(tape_index, read_mode);
        }

      return is_recording_flag;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::stop_recording_operations(
      const bool write_tapes_to_file)
    {
      Assert(is_recording() == true, ExcMessage("Not currently recording..."));

      // Double check that we've actually registered all DoFs
      Assert(n_registered_independent_variables() == n_independent_variables(),
             ExcMessage("Not all values of sensitivities have been recorded!"));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          // Stop taping
          TapedDrivers<ad_type, scalar_type>::stop_taping(active_tape(),
                                                          write_tapes_to_file);

          // Now that we've turned tracing off, we've definitely
          // stopped all tape recording.
          is_recording_flag = false;

          // If the keep_values flag is set, then we expect the user to use this
          // tape immediately after recording it. There is therefore no need to
          // invalidate it. However, there is now also no way to double-check
          // that the newly recorded tape is indeed the active tape.
          if (keep_values == false)
            active_tape_index = numbers::invalid_tape_index;
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
          // set_independent_variable(). This is important because the value of
          // the tapeless independent variables are set and finalized when
          // mark_independent_variable() is called. So we cannot allow this to
          // be done when not in the "recording" phase
          is_recording_flag = false;
          active_tape_index = numbers::invalid_tape_index;
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
      const unsigned int index,
      const ad_type &    func)
    {
      Assert(index < n_dependent_variables(), ExcMessage("Index out of range"));
      Assert(registered_marked_dependent_variables[index] == false,
             ExcMessage(
               "This dependent variable has already been registered."));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(active_tape() != numbers::invalid_tape_index,
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



    /* -------------------- ADHelperCellLevelBase -------------------- */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::ADHelperCellLevelBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      : ADHelperBase<ADNumberTypeCode, ScalarType>(n_independent_variables,
                                                   n_dependent_variables)
    {}



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::register_dof_values(
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
      typename ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::ad_type> &
    ADHelperCellLevelBase<ADNumberTypeCode,
                          ScalarType>::get_sensitive_dof_values()
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // If necessary, initialize the internally stored vector of
      // AD numbers that represents the independent variables
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size() ==
               this->n_independent_variables(),
             ExcInternalError());

      return this->independent_variables;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::vector<
      typename ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::ad_type>
    ADHelperCellLevelBase<ADNumberTypeCode,
                          ScalarType>::get_non_sensitive_dof_values() const
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      std::vector<ad_type> out(this->n_independent_variables(),
                               dealii::internal::NumberType<ad_type>::value(
                                 0.0));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        this->initialize_non_sensitive_independent_variable(i, out[i]);

      return out;
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>::set_dof_values(
      const std::vector<scalar_type> &values)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage(
               "Vector size does not match number of independent variables"));
      for (unsigned int i = 0; i < this->n_independent_variables(); ++i)
        ADHelperBase<ADNumberTypeCode, ScalarType>::set_sensitivity_value(
          i, values[i]);
    }



    /* ------------------ ADHelperEnergyFunctional ------------------ */



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperEnergyFunctional<ADNumberTypeCode, ScalarType>::
      ADHelperEnergyFunctional(const unsigned int n_independent_variables)
      : ADHelperCellLevelBase<ADNumberTypeCode, ScalarType>(
          n_independent_variables,
          1)
    {}



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperEnergyFunctional<ADNumberTypeCode, ScalarType>::
      register_energy_functional(const ad_type &energy)
    {
      Assert(this->n_dependent_variables() == 1, ExcInternalError());
      ADHelperBase<ADNumberTypeCode, ScalarType>::register_dependent_variable(
        0, energy);
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    typename ADHelperEnergyFunctional<ADNumberTypeCode, ScalarType>::scalar_type
    ADHelperEnergyFunctional<ADNumberTypeCode, ScalarType>::compute_energy()
      const
    {
      if (this->keep_values == false ||
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
          "The ADHelperEnergyFunctional class expects there to be only one dependent variable."));

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute value while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          return TapedDrivers<ad_type, scalar_type>::value(
            this->active_tape(), this->independent_variable_values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcInternalError());

          return TapelessDrivers<ad_type, scalar_type>::value(
            this->dependent_variables);
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperEnergyFunctional<ADNumberTypeCode, ScalarType>::compute_residual(
      Vector<scalar_type> &gradient) const
    {
      if (this->keep_values == false ||
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
          "The ADHelperEnergyFunctional class expects there to be only one dependent variable."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately.
      if (gradient.size() != this->n_independent_variables())
        gradient.reinit(this->n_independent_variables(),
                        true /*omit_zeroing_entries*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute gradient while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          TapedDrivers<ad_type, scalar_type>::gradient(
            this->active_tape(), this->independent_variable_values, gradient);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcInternalError());

          TapelessDrivers<ad_type, scalar_type>::gradient(
            this->independent_variables, this->dependent_variables, gradient);
        }
    }



    template <enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperEnergyFunctional<ADNumberTypeCode, ScalarType>::
      compute_linearization(FullMatrix<scalar_type> &hessian) const
    {
      Assert(AD::ADNumberTraits<ad_type>::n_supported_derivative_levels >= 2,
             ExcMessage(
               "Cannot computed function Hessian: AD number type does"
               "not support the calculation of second order derivatives."));

      if (this->keep_values == false)
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
          "The ADHelperEnergyFunctional class expects there to be only one dependent variable."));

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately.
      if (hessian.m() != this->n_independent_variables() &&
          hessian.n() != this->n_independent_variables())
        hessian.reinit({this->n_independent_variables(),
                        this->n_independent_variables()},
                       true /*omit_default_initialization*/);

      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape() != numbers::invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording() == false,
                 ExcMessage(
                   "Cannot compute hessian while tape is being recorded."));
          Assert(this->independent_variable_values.size() ==
                   this->n_independent_variables(),
                 ExcDimensionMismatch(this->independent_variable_values.size(),
                                      this->n_independent_variables()));

          TapedDrivers<ad_type, scalar_type>::hessian(
            this->active_tape(), this->independent_variable_values, hessian);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true,
                 ExcInternalError());
          Assert(this->independent_variables.size() ==
                   this->n_independent_variables(),
                 ExcInternalError());
          TapelessDrivers<ad_type, scalar_type>::hessian(
            this->independent_variables, this->dependent_variables, hessian);
        }
    }


  } // namespace AD
} // namespace Differentiation


/* --- Explicit instantiations --- */
#  ifdef DEAL_II_WITH_ADOLC
#    include "ad_helpers.inst1"
#  endif
#  ifdef DEAL_II_TRILINOS_WITH_SACADO
#    include "ad_helpers.inst2"
#  endif


DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)
