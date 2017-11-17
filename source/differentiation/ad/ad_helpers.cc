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

#if defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_WITH_TRILINOS)

#include <deal.II/differentiation/ad/ad_helpers.h>

#include <type_traits>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace AD
  {

// -------------------------- ADHelperBase ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ADHelperBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      :
      active_tape_index(internal::configure_adtl_and_return_tape_index<ad_type>(n_independent_variables, invalid_tape_index)),
      keep_values(true),
      is_recording(false),
      independent_variable_values(n_independent_variables,
                                  dealii::internal::NumberType<scalar_type>::value(0.0)),
      touched_independent_variables(n_independent_variables,false),
      touched_marked_independent_variables(n_independent_variables,false),
      dependent_variables(n_dependent_variables,
                          dealii::internal::NumberType<ad_type>::value(0.0)),
      touched_dependent_variables(n_dependent_variables,false)
    {
      // TODO: Tapeless: Create some sort of shared static mutex to ensure
      // that only one ADHelperBase can exist in scope at a time?
      // Bad things will happen in the case that adtl::setNumDir is called with
      // different values during the course of tracking some tapeless variables.
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::~ADHelperBase()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::reset_touched_independent_variables ()
    {
      for (typename std::vector<bool>::iterator
           it = touched_independent_variables.begin();
           it != touched_independent_variables.end(); ++it)
        *it = false;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::set_sensitivity_value (
      const scalar_type    &value,
      const unsigned int  index)
    {
      Assert(index < n_independent_variables(),
             ExcMessage("Trying to set the value of a non-existent independent variable."));
      independent_variable_values[index]           = value;
      touched_independent_variables[index]   = true;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::mark_independent_variable (
      ad_type       &out,
      const unsigned int  index) const
    {
      Assert(index < n_independent_variables(), ExcInternalError());
      Assert(touched_independent_variables[index] == true, ExcInternalError());

      if (index > 0)
        {
          Assert(touched_marked_independent_variables[index-1] == true,
                 ExcMessage("Need to extract sensitivities in the order they're created."));
        }

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(active_tape()!=invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(is_recording==true,
                 ExcMessage("Only valid during recording"));
        }

      internal::Marking<ad_type>::independent_variable(
        independent_variable_values[index], index,
        this->n_independent_variables(), out);
      touched_marked_independent_variables[index] = true;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::finalize_sensitive_independent_variables () const
    {
      // Double check that we've actually touched all DoFs
      Assert(n_touched_independent_variables() == n_independent_variables(),
             ExcMessage("Not all values of sensitivities have been recorded!"));

      // This should happen only once
      if (this->independent_variables.size() == 0)
        {
          this->independent_variables = std::vector<ad_type> (this->n_independent_variables(), dealii::internal::NumberType<ad_type>::value(0.0));

          // Indicate the sensitivity that each entry represents
          for (unsigned int i=0; i<this->n_independent_variables(); ++i)
            this->mark_independent_variable(this->independent_variables[i], i);
        }
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::get_independent_variable (
      ad_type       &out,
      const unsigned int  index) const
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(active_tape()!=invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(is_recording==false,
             ExcMessage("Only valid outside of recording"));

      Assert(index < n_independent_variables(), ExcInternalError());
      Assert(touched_independent_variables[index] == true, ExcInternalError());

      out = independent_variable_values[index];
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::n_touched_independent_variables () const
    {
      return std::accumulate(touched_independent_variables.begin(),
                             touched_independent_variables.end(), 0u);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::size_t
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::n_independent_variables() const
    {
      return independent_variable_values.size();
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::n_touched_dependent_variables () const
    {
      return std::accumulate(touched_dependent_variables.begin(),
                             touched_dependent_variables.end(), 0u);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::size_t
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::n_dependent_variables() const
    {
      return dependent_variables.size();
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    int
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::active_tape() const
    {
      return active_tape_index;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::print_values(std::ostream &stream) const
    {
      for (unsigned int i=0; i<n_independent_variables(); i++)
        stream
            << independent_variable_values[i]
            << (i<(n_independent_variables()-1) ? "," : "")
            << std::flush;
      stream << std::endl;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode,typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::print_tape_stats(
      std::ostream       &stream,
      const unsigned int &tape_index) const
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        return;

      Assert(registered_tapes.find(tape_index) != registered_tapes.end(),
             ExcMessage("Tape number not registered"));

#ifdef DEAL_II_WITH_ADOLC
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
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
              << "Number of parameters (doubles) interchangable without retaping: " << counts[17] << "\n";
        }
#endif
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::reset (const unsigned int n_independent_variables,
                                                          const unsigned int n_dependent_variables)
    {
      active_tape_index = internal::configure_adtl_and_return_tape_index<ad_type>(n_independent_variables, invalid_tape_index);
      registered_tapes.clear();
      is_recording = false;

      const unsigned int new_n_independent_variables =
        (n_independent_variables != 0 ? n_independent_variables : this->n_independent_variables());

      independent_variable_values =  std::vector<scalar_type>(new_n_independent_variables,
                                                              dealii::internal::NumberType<scalar_type>::value(0.0));
      touched_independent_variables = std::vector<bool>(new_n_independent_variables,false);
      touched_marked_independent_variables = std::vector<bool>(new_n_independent_variables,false);

      const unsigned int new_n_dependent_variables =
        (n_dependent_variables != 0 ? n_dependent_variables : this->n_dependent_variables());

      dependent_variables = std::vector<ad_type>(new_n_dependent_variables,
                                                 dealii::internal::NumberType<ad_type>::value(0.0));
      touched_dependent_variables = std::vector<bool>(new_n_dependent_variables,false);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::configure_tapeless_mode (const unsigned int &n_independent_variables)
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_tapeless)
        {
          internal::configure_adtl<ad_type>(n_independent_variables);
          // In order to ensure that the settings remain for the entire duration of the simulation,
          // we create a global live variable that doesn't go out of scope.
          static ad_type num = 0.0;
        }
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::activate_tape(const unsigned int &tape_index)
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(tape_index!=invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(tape_index<max_tape_index,
                 ExcMessage("Tape index exceeds maximum allowable value"));
          Assert(registered_tapes.find(tape_index) != registered_tapes.end(),
                 ExcMessage("Tape number not registered"));
          active_tape_index = tape_index;
          reset_touched_independent_variables();
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true, ExcInternalError());
          active_tape_index = 1; // Some dummy value
        }
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::enable_record_sensitivities(
      const unsigned int &tape_index,
      const bool          overwrite_tape,
      const bool          keep)
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          if (overwrite_tape != true)
            {
              Assert(is_recording == false,
                     ExcMessage("Already recording..."));
            }
          if (registered_tapes.find(tape_index) == registered_tapes.end() ||
              overwrite_tape == true)
            {
              registered_tapes.insert(tape_index);
              activate_tape(tape_index);
              trace_on(active_tape(),keep);
              reset_touched_independent_variables();
              keep_values = keep;
              is_recording = true;
            }
        }
      else
        {
          // Tapeless mode
          is_recording = true;
          activate_tape(tape_index);
        }

      return is_recording;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::disable_record_sensitivities(
      const bool write_tapes_to_file)
    {
      Assert(is_recording == true,
             ExcMessage("Not currently recording..."));

      // Double check that we've actually touched all DoFs
      Assert(n_touched_independent_variables() == n_independent_variables(),
             ExcMessage("Not all values of sensitivities have been recorded!"));

      if (ADNumberTraits<ad_type>::is_tapeless == true)
        {
          // Double check that we've actually touched dependent variables
          Assert(n_touched_dependent_variables() == n_dependent_variables(),
                 ExcMessage("Not all dependent variables have been set!"));

          // By changing this flag, we ensure that the we can no longer
          // alter the values of the dependent variables using
          // set_independent_variable(). This is important because the value of
          // the tapeless independent variables are set and finalized when
          // mark_independent_variable() is called. So we cannot allow this to
          // be done when not in the "recording" phase
          is_recording = false;
          active_tape_index = invalid_tape_index;
          return;
        }
      else
        {
#ifdef DEAL_II_WITH_ADOLC
          if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
            {
              if (write_tapes_to_file)
                {
                  trace_off(active_tape()); // Slow

                  std::vector<std::size_t> counts (STAT_SIZE);
                  ::tapestats(active_tape(), counts.data());
                }
              else
                trace_off(); // Fast(er)
            }
#endif

          // Now that we've turned tracing off, we've definitely
          // stopped all tape recording.
          is_recording = false;

          // If the keep_values flag is set, then we expect the user to use this tape
          // immediately after recording it. There is therefore no need to invalidate
          // it. However, there is now also no way to double-check that the newly
          // recorded tape is indeed the active tape.
          if (keep_values == false)
            active_tape_index = invalid_tape_index;
        }
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperBase<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable (
      const ad_type &func,
      const unsigned int  index)
    {
      Assert(index<n_dependent_variables(),
             ExcMessage("Index out of range"));
      Assert(touched_dependent_variables[index] == false,
             ExcMessage("This dependent variable has already been registered."));

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(active_tape()!=invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(is_recording == true,
                 ExcMessage("Must be recording when registering dependent variables."));
        }

      // See Adol-C manual section 1.4
      // Note: Even though it appears that we're doing nothing in particular here,
      // the following steps are in fact being recorded on the Adol-C tape.
      internal::Marking<ad_type>::dependent_variable(dependent_variables[index], func);
      touched_dependent_variables[index] = true;
    }



// -------------------------- ADHelperPointLevelFunctionsBase ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::ADHelperPointLevelFunctionsBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      :
      ADHelperBase<dim,ADNumberTypeCode,ScalarType>(n_independent_variables, n_dependent_variables),
      symmetric_independent_variables(n_independent_variables,false)
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::~ADHelperPointLevelFunctionsBase()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::reset (
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
    {
      ADHelperBase<dim,ADNumberTypeCode,ScalarType>::reset(n_independent_variables,n_dependent_variables);

      const unsigned int new_n_independent_variables =
        (n_independent_variables != 0 ? n_independent_variables : this->n_independent_variables());
      symmetric_independent_variables = std::vector<bool>(new_n_independent_variables,false);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    bool
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::is_symmetric_independent_variable (const unsigned int index) const
    {
      Assert(index<symmetric_independent_variables.size(),
             ExcInternalError());
      return symmetric_independent_variables[index];
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    unsigned int
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::n_symmetric_independent_variables () const
    {
      return std::accumulate(symmetric_independent_variables.begin(),
                             symmetric_independent_variables.end(), 0u);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::register_independent_variables (const std::vector<scalar_type> &values)
    {
      // This is actually the same thing the set_independent_variable function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage("Vector size does not match number of independent variables"));
#ifdef DEBUG
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        {
          Assert(this->touched_independent_variables[i] == false,
                 ExcMessage("Independent variables already registered."));
        }
#endif
      set_independent_variables(values);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    const std::vector<typename ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::ad_type> &
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::get_sensitive_variables ()
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // If necessary, initialize the internally stored vector of
      // AD numbers that represents the independent variables
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size()==this->n_independent_variables(),
             ExcInternalError());

      return this->independent_variables;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::set_sensitivity_value (
      const scalar_type    &value,
      const unsigned int  index,
      const bool          symmetric_dof)
    {
      ADHelperBase<dim,ADNumberTypeCode,ScalarType>::set_sensitivity_value(value,index);
      Assert(index < this->n_independent_variables(),
             ExcMessage("Trying to set the symmetry flag of a non-existent independent variable."));
      Assert(index < symmetric_independent_variables.size(),
             ExcInternalError());
      symmetric_independent_variables[index] = symmetric_dof;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::set_independent_variables (const std::vector<scalar_type> &values)
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage("Vector size does not match number of independent variables"));
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        ADHelperBase<dim,ADNumberTypeCode,ScalarType>::set_sensitivity_value(values[i], i);
    }



    // -------------------------- Internal functions ----------------------



    namespace internal
    {
      namespace
      {
        /**
         * Adol-C only has drivers for doubles, and so floats are not intrinsically
         * supported. This wrapper struct works around the issue when necessary.
         */
        template<typename ScalarType>
        struct AdolCTapedDrivers;

        template<>
        struct AdolCTapedDrivers<double>
        {
          typedef double scalar_type;

          // === Scalar drivers ===

          static scalar_type
          value (const unsigned int             &active_tape,
                 const unsigned int             &n_independent_variables,
                 const std::vector<scalar_type> &independent_variables)
          {
            double *f = new double();
            ::function(active_tape,
                       1, // Only one dependent variable
                       n_independent_variables,
                       const_cast<double *>(independent_variables.data()),
                       f);

            const double value = f[0];

            // Cleanup :-/
            delete f;
            f = nullptr;

            return value;
          }

          static void
          gradient (Vector<scalar_type>            &gradient,
                    const unsigned int             &active_tape,
                    const unsigned int             &n_independent_variables,
                    const std::vector<scalar_type> &independent_variables)
          {
            Assert(gradient.size() == n_independent_variables,
                   ExcMessage("The length of the gradient vector must be equal "
                              "to the number of independent variables."));

            scalar_type *g = new scalar_type[n_independent_variables]; // Use smart pointers or std_cxx11::array?
            ::gradient(active_tape,
                       n_independent_variables,
                       const_cast<scalar_type *>(independent_variables.data()),
                       g);

            for (unsigned int i=0; i<n_independent_variables; ++i)
              gradient[i] = g[i];

            // Cleanup :-/
            delete[] g;
            g = nullptr;
          }

          static void
          hessian (FullMatrix<scalar_type>        &hessian,
                   const unsigned int             &active_tape,
                   const unsigned int             &n_independent_variables,
                   const std::vector<scalar_type> &independent_variables)
          {
            Assert(hessian.m() == n_independent_variables,
                   ExcMessage("The hessian row length must be equal "
                              "to the number of independent variables."));

            Assert(hessian.n() == n_independent_variables,
                   ExcMessage("The hessian column length must be equal "
                              "to the number of independent variables."));

            scalar_type **H = new scalar_type*[n_independent_variables]; // Use smart pointers or std_cxx11::array?
            for (unsigned int i=0; i<n_independent_variables; ++i)
              H[i] = new scalar_type[i+1]; // Symmetry

            ::hessian(active_tape,
                      n_independent_variables,
                      const_cast<scalar_type *>(independent_variables.data()),
                      H);

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
          values (Vector<scalar_type>            &values,
                  const unsigned int             &active_tape,
                  const unsigned int             &n_dependent_variables,
                  const unsigned int             &n_independent_variables,
                  const std::vector<scalar_type> &independent_variables)
          {
            Assert(values.size() == n_dependent_variables,
                   ExcMessage("The length of the dependent function vector must "
                              " be equal to the number of dependent variables."));

            scalar_type *f = new scalar_type[n_dependent_variables]; // Use smart pointers or std_cxx11::array?
            ::function(active_tape,
                       n_dependent_variables,
                       n_independent_variables,
                       const_cast<scalar_type *>(independent_variables.data()),
                       f);

            for (unsigned int i=0; i<n_dependent_variables; i++)
              values[i] = f[i];

            // Cleanup :-/
            delete[] f;
            f = nullptr;
          }

          static void
          jacobian (FullMatrix<scalar_type>        &jacobian,
                    const unsigned int             &active_tape,
                    const unsigned int             &n_dependent_variables,
                    const unsigned int             &n_independent_variables,
                    const std::vector<scalar_type> &independent_variables)
          {
            scalar_type **J = new scalar_type*[n_dependent_variables]; // Use smart pointers or std_cxx11::array?
            for (unsigned int i=0; i<n_dependent_variables; ++i)
              J[i] = new scalar_type[n_independent_variables];

            ::jacobian(active_tape,
                       n_dependent_variables,
                       n_independent_variables,
                       independent_variables.data(),
                       J);

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


        template<>
        struct AdolCTapedDrivers<float>
        {
          typedef float scalar_type;

          static std::vector<double>
          vector_float_to_double (const std::vector<float> &in)
          {
            std::vector<double> out (in.size());
            std::copy(in.begin(), in.end(), out.begin());
            return out;
          }

          // === Scalar drivers ===

          static scalar_type
          value (const unsigned int             &active_tape,
                 const unsigned int             &n_independent_variables,
                 const std::vector<scalar_type> &independent_variables)
          {

            return AdolCTapedDrivers<double>::value(
                     active_tape,
                     n_independent_variables,
                     vector_float_to_double(independent_variables));
          }

          static void
          gradient (Vector<scalar_type>            &gradient,
                    const unsigned int             &active_tape,
                    const unsigned int             &n_independent_variables,
                    const std::vector<scalar_type> &independent_variables)
          {
            Vector<double> gradient_double (gradient.size());
            AdolCTapedDrivers<double>::gradient(
              gradient_double,
              active_tape,
              n_independent_variables,
              vector_float_to_double(independent_variables));
            gradient = gradient_double;
          }

          static void
          hessian (FullMatrix<scalar_type>        &hessian,
                   const unsigned int             &active_tape,
                   const unsigned int             &n_independent_variables,
                   const std::vector<scalar_type> &independent_variables)
          {
            FullMatrix<double> hessian_double (hessian.m(), hessian.n());
            AdolCTapedDrivers<double>::hessian(
              hessian_double,
              active_tape,
              n_independent_variables,
              vector_float_to_double(independent_variables));
            hessian = hessian_double;
          }

          // === Vector drivers ===

          static void
          values (Vector<scalar_type>            &values,
                  const unsigned int             &active_tape,
                  const unsigned int             &n_dependent_variables,
                  const unsigned int             &n_independent_variables,
                  const std::vector<scalar_type> &independent_variables)
          {
            Vector<double> values_double (values.size());
            AdolCTapedDrivers<double>::values(
              values_double,
              active_tape,
              n_dependent_variables,
              n_independent_variables,
              vector_float_to_double(independent_variables));
            values = values_double;
          }

          static void
          jacobian (FullMatrix<scalar_type>        &jacobian,
                    const unsigned int             &active_tape,
                    const unsigned int             &n_dependent_variables,
                    const unsigned int             &n_independent_variables,
                    const std::vector<scalar_type> &independent_variables)
          {
            FullMatrix<double> jacobian_double (jacobian.m(), jacobian.n());
            AdolCTapedDrivers<double>::jacobian(
              jacobian_double,
              active_tape,
              n_dependent_variables,
              n_independent_variables,
              vector_float_to_double(independent_variables));
            jacobian = jacobian_double;
          }
        };

      }
    }



    // -------------------------- ADHelperScalarFunction ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::ADHelperScalarFunction(const unsigned int n_independent_variables)
      :
      ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>(n_independent_variables, 1)
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::~ADHelperScalarFunction()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable (const ad_type &func)
    {
      Assert(this->n_dependent_variables() == 1, ExcInternalError());
      ADHelperBase<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable(func,0);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::compute_value() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(this->n_dependent_variables() == 1,
             ExcMessage("Only valid for one dependent variable."));

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute values while tape is being recorded."));

          // This is a work-around for the fact that Adol-C does not perform
          // calculations with floats (only doubles)
          return internal::AdolCTapedDrivers<scalar_type>::value(
                   this->active_tape(),
                   this->n_independent_variables(),
                   this->independent_variable_values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true, ExcInternalError());
          return ADNumberTraits<ad_type>::get_scalar_value(this->dependent_variables[0]);
        }
    }


    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    Vector<typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::compute_gradient() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(this->n_dependent_variables() == 1,
             ExcMessage("Only valid for one dependent variable."));

      Vector<scalar_type> gradient (this->n_independent_variables());

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute gradient while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::gradient(
            gradient,
            this->active_tape(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else if (ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad ||
               ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad_dfad)
        {
          Assert(this->independent_variables.size() == this->n_independent_variables(), ExcInternalError());
          // In reverse mode, the gradients are computed from the
          // independent variables (i.e. the adjoint)
          internal::reverse_mode_dependent_variable_activation(const_cast<ad_type &>(this->dependent_variables[0]));
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            gradient[i] = internal::NumberType<scalar_type>::value(ADNumberTraits<ad_type>::get_directional_derivative(
                                                                     this->independent_variables[i],
                                                                     0 /*This number doesn't really matter*/));
        }
      else
        {
          Assert((ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_tapeless ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad_dfad),
                 ExcMessage("An unexpected AD type has fallen through to the default case."));
          // In forward mode, the gradients are computed from the
          // dependent variables
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            gradient[i] = internal::NumberType<scalar_type>::value(ADNumberTraits<ad_type>::get_directional_derivative(
                                                                     this->dependent_variables[0], i));
        }

      // Account for symmetries of tensor components
      for (unsigned int i=0; i<this->n_independent_variables(); i++)
        {
          if (this->is_symmetric_independent_variable(i) == true)
            gradient[i] *= 0.5;
        }

      return gradient;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    FullMatrix<typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::compute_hessian() const
    {
      AssertThrow(AD::ADNumberTraits<ad_type>::n_supported_derivative_levels >= 2,
                  ExcMessage("Cannot computed function Hessian: AD number type does not support the calculation of second order derivatives."));

      if (this->keep_values == false)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(this->n_dependent_variables() == 1,
             ExcMessage("Only valid for one dependent variable."));

      FullMatrix<scalar_type> hessian (this->n_independent_variables(),
                                       this->n_independent_variables());

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute hessian while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::hessian(
            hessian,
            this->active_tape(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else if (ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad_dfad)
        {
          Assert(this->independent_variables.size() == this->n_independent_variables(), ExcInternalError());
          // In reverse mode, the gradients are computed from the
          // independent variables (i.e. the adjoint)
          internal::reverse_mode_dependent_variable_activation(const_cast<ad_type &>(this->dependent_variables[0]));
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            {
              typedef typename ADNumberTraits<ad_type>::derivative_type derivative_type;
              const derivative_type gradient_i
                = ADNumberTraits<ad_type>::get_directional_derivative(this->independent_variables[i], i);

              for (unsigned int j=0; j <= i; ++j) // Symmetry
                {
                  // Extract higher-order directional derivatives. Depending on the AD number type,
                  // the result may be another AD number or a floating point value.
                  const scalar_type hessian_ij = internal::NumberType<scalar_type>::value(
                                                   ADNumberTraits<derivative_type>::get_directional_derivative(gradient_i, j));
                  hessian[i][j] = hessian_ij;
                  if (i != j)
                    hessian[j][i] = hessian_ij;  // Symmetry
                }
            }
        }
      else
        {
          Assert((ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad_dfad),
                 ExcMessage("An unexpected AD type has fallen through to the default case."));
          // In forward mode, the gradients are computed from the
          // dependent variables
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            {
              typedef typename ADNumberTraits<ad_type>::derivative_type derivative_type;
              const derivative_type gradient_i
                = ADNumberTraits<ad_type>::get_directional_derivative(this->dependent_variables[0], i);

              for (unsigned int j=0; j <= i; ++j) // Symmetry
                {
                  const scalar_type hessian_ij = internal::NumberType<scalar_type>::value(
                                                   ADNumberTraits<derivative_type>::get_directional_derivative(gradient_i, j));
                  hessian[i][j] = hessian_ij;
                  if (i != j)
                    hessian[j][i] = hessian_ij;  // Symmetry
                }
            }
        }

      // Account for symmetries of tensor components
      for (unsigned int i=0; i<this->n_independent_variables(); i++)
        for (unsigned int j=0; j<i+1; j++)
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

      return hessian;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    Tensor<0,dim,typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::extract_hessian_component(
      const FullMatrix<scalar_type>      &hessian,
      const FEValuesExtractors::Scalar &extractor_row,
      const FEValuesExtractors::Scalar &extractor_col) const
    {
      // NOTE: It is necessary to make special provision for the case when the HessianType
      //       is scalar. Unfortunately Tensor<0,dim> does not provide the function
      //       unrolled_to_component_indices!
      // NOTE: The order of components must be consistently defined throughout this class.
      Tensor<0,dim,scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the matrix values
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set (internal::extract_index_set<dim>(extractor_col));
      Assert(row_index_set.size() == 1, ExcInternalError());
      Assert(col_index_set.size() == 1, ExcInternalError());

      internal::set_tensor_entry(out, 0,
                                 hessian[row_index_set[0]][col_index_set[0]]);

      return out;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    SymmetricTensor<4,dim,typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::extract_hessian_component(
      const FullMatrix<scalar_type>                   &hessian,
      const FEValuesExtractors::SymmetricTensor<2>  &extractor_row,
      const FEValuesExtractors::SymmetricTensor<2>  &extractor_col) const
    {
      // NOTE: The order of components must be consistently defined throughout this class.
      // NOTE: We require a specialisation for rank-4 symmetric tensors because they
      //       do not define their rank, and setting data using TableIndices is somewhat
      //       specialised as well.
      SymmetricTensor<4,dim,scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the matrix values
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set (internal::extract_index_set<dim>(extractor_col));

      for (unsigned int r=0; r<row_index_set.size(); ++r)
        for (unsigned int c=0; c<col_index_set.size(); ++c)
          {
            internal::set_tensor_entry(out, r, c,
                                       hessian[row_index_set[r]][col_index_set[c]]);
          }

      return out;
    }



// -------------------------- ADHelperVectorFunction ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::ADHelperVectorFunction(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      :
      ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>(n_independent_variables,
          n_dependent_variables)
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::~ADHelperVectorFunction()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::register_dependent_variables (const std::vector<ad_type> &funcs)
    {
      Assert(funcs.size() == this->n_dependent_variables(),
             ExcMessage("Vector size does not match number of dependent variables"));
      for (unsigned int i=0; i<this->n_dependent_variables(); ++i)
        ADHelperBase<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable(funcs[i],i);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    Vector<typename ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::compute_values() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Vector<scalar_type> values (this->n_dependent_variables());
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute values while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::values(
            values,
            this->active_tape(),
            this->n_dependent_variables(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true, ExcInternalError());
          for (unsigned int i=0; i<this->n_dependent_variables(); i++)
            values[i] = ADNumberTraits<ad_type>::get_scalar_value(this->dependent_variables[i]);
        }

      return values;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    FullMatrix<typename ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::compute_jacobian() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      FullMatrix<scalar_type> jacobian (this->n_dependent_variables(),
                                        this->n_independent_variables());
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute jacobian while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::jacobian(
            jacobian,
            this->active_tape(),
            this->n_dependent_variables(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else if (ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad ||
               ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad_dfad)
        {
          Assert(this->independent_variables.size() == this->n_independent_variables(), ExcInternalError());
          // In reverse mode, the gradients are computed from the
          // independent variables (i.e. the adjoint).
          // For a demonstration of why this accumulation process is
          // required, see the unit tests
          // sacado/basic_01b.cc and sacado/basic_02b.cc
          // Here we also take into consideration the derivative type:
          // The Sacado number may be of the nested variety, in which
          // case the effect of the accumulation process on the
          // sensitivities of the nested number need to be accounted for.
          typedef typename ADNumberTraits<ad_type>::derivative_type AccumulationType;
          std::vector<AccumulationType> rad_accumulation (
            this->n_independent_variables(),
            dealii::internal::NumberType<AccumulationType>::value(0.0));
          for (unsigned int i=0; i<this->n_dependent_variables(); i++)
            {
              internal::reverse_mode_dependent_variable_activation(const_cast<ad_type &>(this->dependent_variables[i]));
              for (unsigned int j=0; j<this->n_independent_variables(); j++)
                {
                  const AccumulationType df_i_dx_j
                    = ADNumberTraits<ad_type>::get_directional_derivative(
                        this->independent_variables[j], i /*This number doesn't really matter*/)
                      - rad_accumulation[j];
                  jacobian[i][j] = internal::NumberType<scalar_type>::value(df_i_dx_j);
                  rad_accumulation[j] += df_i_dx_j;
                }
            }
        }
      else
        {
          Assert((ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_tapeless ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad_dfad),
                 ExcMessage("An unexpected AD type has fallen through to the default case."));
          // In forward mode, the gradients are computed from the
          // dependent variables
          for (unsigned int i=0; i<this->n_dependent_variables(); i++)
            for (unsigned int j=0; j<this->n_independent_variables(); j++)
              jacobian[i][j] = internal::NumberType<scalar_type>::value(ADNumberTraits<ad_type>::get_directional_derivative(
                                                                          this->dependent_variables[i], j));
        }

      for (unsigned int i=0; i<this->n_dependent_variables(); i++)
        for (unsigned int j=0; j<this->n_independent_variables(); j++)
          // Because we perform just a single differentiation
          // operation with respect to the "column" variables,
          // we only need to consider them for symmetry conditions.
          if (this->is_symmetric_independent_variable(j) == true)
            jacobian[i][j] *= 0.5;

      return jacobian;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    Tensor<0,dim,typename ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::extract_jacobian_component(
      const FullMatrix<scalar_type>      &jacobian,
      const FEValuesExtractors::Scalar &extractor_row,
      const FEValuesExtractors::Scalar &extractor_col) const
    {
      // NOTE: It is necessary to make special provision for the case when the HessianType
      //       is scalar. Unfortunately Tensor<0,dim> does not provide the function
      //       unrolled_to_component_indices!
      // NOTE: The order of components must be consistently defined throughout this class.
      Tensor<0,dim,scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the matrix values
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set (internal::extract_index_set<dim>(extractor_col));
      Assert(row_index_set.size() == 1, ExcInternalError());
      Assert(col_index_set.size() == 1, ExcInternalError());

      internal::set_tensor_entry(out, 0,
                                 jacobian[row_index_set[0]][col_index_set[0]]);

      return out;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    SymmetricTensor<4,dim,typename ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::extract_jacobian_component(
      const FullMatrix<scalar_type>                   &jacobian,
      const FEValuesExtractors::SymmetricTensor<2>  &extractor_row,
      const FEValuesExtractors::SymmetricTensor<2>  &extractor_col) const
    {
      // NOTE: The order of components must be consistently defined throughout this class.
      // NOTE: We require a specialisation for rank-4 symmetric tensors because they
      //       do not define their rank, and setting data using TableIndices is somewhat
      //       specialised as well.
      SymmetricTensor<4,dim,scalar_type> out;

      // Get indexsets for the subblocks from which we wish to extract the matrix values
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row));
      const std::vector<unsigned int> col_index_set (internal::extract_index_set<dim>(extractor_col));

      for (unsigned int r=0; r<row_index_set.size(); ++r)
        for (unsigned int c=0; c<col_index_set.size(); ++c)
          {
            internal::set_tensor_entry(out, r, c,
                                       jacobian[row_index_set[r]][col_index_set[c]]);
          }

      return out;
    }



// -------------------------- ADHelperCellLevelBase ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::ADHelperCellLevelBase(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      :
      ADHelperBase<dim,ADNumberTypeCode,ScalarType>(n_independent_variables,
                                                    n_dependent_variables)
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::~ADHelperCellLevelBase()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::register_dof_values (const std::vector<scalar_type> &dof_values)
    {
      // This is actually the same thing the set_independent_variable function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
      Assert(dof_values.size() == this->n_independent_variables(),
             ExcMessage("Vector size does not match number of independent variables"));
#ifdef DEBUG
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        {
          Assert(this->touched_independent_variables[i] == false,
                 ExcMessage("Independent variables already registered."));
        }
#endif
      set_dof_values(dof_values);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    const std::vector<typename ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::ad_type> &
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::get_sensitive_dof_values ()
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // If necessary, initialize the internally stored vector of
      // AD numbers that represents the independent variables
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size()==this->n_independent_variables(),
             ExcInternalError());

      return this->independent_variables;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    std::vector<typename ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::ad_type>
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::get_non_sensitive_dof_values ()
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      std::vector<ad_type> out (this->n_independent_variables(), dealii::internal::NumberType<ad_type>::value(0.0));
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        this->get_independent_variable(out[i], i);

      return out;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::set_dof_values (const std::vector<scalar_type> &values)
    {
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(values.size() == this->n_independent_variables(),
             ExcMessage("Vector size does not match number of independent variables"));
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        ADHelperBase<dim,ADNumberTypeCode,ScalarType>::set_sensitivity_value(values[i], i);
    }



    // -------------------------- ADHelperVariationalFormulation ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::ADHelperVariationalFormulation(
      const unsigned int n_independent_variables)
      :
      ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>(n_independent_variables,1)
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::~ADHelperVariationalFormulation()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::register_energy_functional (const ad_type &energy)
    {
      Assert(this->n_dependent_variables() == 1, ExcInternalError());
      ADHelperBase<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable(energy,0);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    Vector<typename ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::compute_residual() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(this->n_dependent_variables() == 1,
             ExcMessage("Only valid for one dependent variable."));

      Vector<scalar_type> gradient (this->n_independent_variables());

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute gradient while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::gradient(
            gradient,
            this->active_tape(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else if (ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad ||
               ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad_dfad)
        {
          Assert(this->independent_variables.size() == this->n_independent_variables(), ExcInternalError());
          // In reverse mode, the gradients are computed from the
          // independent variables (i.e. the adjoint)
          internal::reverse_mode_dependent_variable_activation(const_cast<ad_type &>(this->dependent_variables[0]));
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            gradient[i] = internal::NumberType<scalar_type>::value(ADNumberTraits<ad_type>::get_directional_derivative(
                                                                     this->independent_variables[i],
                                                                     0 /*This number doesn't really matter*/));
        }
      else
        {

          Assert((ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_tapeless ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad_dfad),
                 ExcMessage("An unexpected AD type has fallen through to the default case."));
          // In forward mode, the gradients are computed from the
          // dependent variables
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            gradient[i] = internal::NumberType<scalar_type>::value(ADNumberTraits<ad_type>::get_directional_derivative(
                                                                     this->dependent_variables[0], i));
        }

      return gradient;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    FullMatrix<typename ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperVariationalFormulation<dim,ADNumberTypeCode,ScalarType>::compute_linearization() const
    {
      AssertThrow(AD::ADNumberTraits<ad_type>::n_supported_derivative_levels >= 2,
                  ExcMessage("Cannot computed functional linearization: AD number type does not support the calculation of second order derivatives."));

      if (this->keep_values == false)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Assert(this->n_dependent_variables() == 1,
             ExcMessage("Only valid for one dependent variable."));

      FullMatrix<scalar_type> hessian (this->n_independent_variables(),
                                       this->n_independent_variables());

      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute hessian while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::hessian(
            hessian,
            this->active_tape(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else if (ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad_dfad)
        {
          Assert(this->independent_variables.size() == this->n_independent_variables(), ExcInternalError());
          // In reverse mode, the gradients are computed from the
          // independent variables (i.e. the adjoint)
          internal::reverse_mode_dependent_variable_activation(const_cast<ad_type &>(this->dependent_variables[0]));
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            {
              typedef typename ADNumberTraits<ad_type>::derivative_type derivative_type;
              const derivative_type gradient_i
                = ADNumberTraits<ad_type>::get_directional_derivative(this->independent_variables[i], i);

              for (unsigned int j=0; j <= i; ++j) // Symmetry
                {
                  const scalar_type hessian_ij = internal::NumberType<scalar_type>::value(
                                                   ADNumberTraits<derivative_type>::get_directional_derivative(gradient_i, j));
                  hessian[i][j] = hessian_ij;
                  if (i != j)
                    hessian[j][i] = hessian_ij;  // Symmetry
                }
            }
        }
      else
        {
          Assert((ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad_dfad),
                 ExcMessage("An unexpected AD type has fallen through to the default case."));
          // In forward mode, the gradients are computed from the
          // dependent variables
          for (unsigned int i=0; i<this->n_independent_variables(); i++)
            {
              typedef typename ADNumberTraits<ad_type>::derivative_type derivative_type;
              const derivative_type gradient_i
                = ADNumberTraits<ad_type>::get_directional_derivative(this->dependent_variables[0], i);

              for (unsigned int j=0; j <= i; ++j) // Symmetry
                {
                  const scalar_type hessian_ij = internal::NumberType<scalar_type>::value(
                                                   ADNumberTraits<derivative_type>::get_directional_derivative(gradient_i, j));
                  hessian[i][j] = hessian_ij;
                  if (i != j)
                    hessian[j][i] = hessian_ij;  // Symmetry
                }
            }
        }

      return hessian;
    }


// -------------------------- ADHelperResidualLinearisation ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::ADHelperResidualLinearisation(
      const unsigned int n_independent_variables,
      const unsigned int n_dependent_variables)
      :
      ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>(n_independent_variables,
                                                             n_dependent_variables)
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::~ADHelperResidualLinearisation()
    { }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    void
    ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::register_residual_vector (const std::vector<ad_type> &residual)
    {
      Assert(residual.size() == this->n_dependent_variables(),
             ExcMessage("Vector size does not match number of dependent variables"));
      for (unsigned int i=0; i<this->n_dependent_variables(); ++i)
        ADHelperBase<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable(residual[i],i);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    Vector<typename ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::compute_residual() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      Vector<scalar_type> values (this->n_dependent_variables());
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute values while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::values(
            values,
            this->active_tape(),
            this->n_dependent_variables(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else
        {
          Assert(ADNumberTraits<ad_type>::is_tapeless == true, ExcInternalError());
          for (unsigned int i=0; i<this->n_dependent_variables(); i++)
            values[i] = ADNumberTraits<ad_type>::get_scalar_value(this->dependent_variables[i]);
        }

      return values;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    FullMatrix<typename ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::scalar_type>
    ADHelperResidualLinearisation<dim,ADNumberTypeCode,ScalarType>::compute_linearization() const
    {
      if (this->keep_values == false ||
          ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->n_touched_independent_variables() == this->n_independent_variables(),
                 ExcMessage("Not all values of sensitivities have been registered or subsequently set!"));
        }
      Assert(this->n_touched_dependent_variables() == this->n_dependent_variables(),
             ExcMessage("Not all dependent variables have been registered."));

      FullMatrix<scalar_type> jacobian (this->n_dependent_variables(),
                                        this->n_independent_variables());
      if (ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_taped)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
          Assert(this->is_recording == false,
                 ExcMessage("Cannot compute hessian while tape is being recorded."));

          internal::AdolCTapedDrivers<scalar_type>::jacobian(
            jacobian,
            this->active_tape(),
            this->n_dependent_variables(),
            this->n_independent_variables(),
            this->independent_variable_values);
        }
      else if (ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad ||
               ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_rad_dfad)
        {
          Assert(this->independent_variables.size() == this->n_independent_variables(), ExcInternalError());
          // In reverse mode, the gradients are computed from the
          // independent variables (i.e. the adjoint).
          // For a demonstration of why this accumulation process is
          // required, see the unit tests
          // sacado/basic_01b.cc and sacado/basic_02b.cc
          // Here we also take into consideration the derivative type:
          // The Sacado number may be of the nested variety, in which
          // case the effect of the accumulation process on the
          // sensitivities of the nested number need to be accounted for.
          typedef typename ADNumberTraits<ad_type>::derivative_type AccumulationType;
          std::vector<AccumulationType> rad_accumulation (
            this->n_independent_variables(),
            dealii::internal::NumberType<AccumulationType>::value(0.0));
          for (unsigned int i=0; i<this->n_dependent_variables(); i++)
            {
              internal::reverse_mode_dependent_variable_activation(const_cast<ad_type &>(this->dependent_variables[i]));
              for (unsigned int j=0; j<this->n_independent_variables(); j++)
                {
                  const AccumulationType df_i_dx_j
                    = ADNumberTraits<ad_type>::get_directional_derivative(
                        this->independent_variables[j], i /*This number doesn't really matter*/)
                      - rad_accumulation[j];
                  jacobian[i][j] = internal::NumberType<scalar_type>::value(df_i_dx_j);
                  rad_accumulation[j] += df_i_dx_j;
                }
            }
        }
      else
        {
          Assert((ADNumberTraits<ad_type>::type_code == NumberTypes::adolc_tapeless ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad ||
                  ADNumberTraits<ad_type>::type_code == NumberTypes::sacado_dfad_dfad),
                 ExcMessage("An unexpected AD type has fallen through to the default case."));
          // In forward mode, the gradients are computed from the
          // dependent variables
          for (unsigned int i=0; i<this->n_dependent_variables(); i++)
            for (unsigned int j=0; j<this->n_independent_variables(); j++)
              jacobian[i][j] = internal::NumberType<scalar_type>::value(ADNumberTraits<ad_type>::get_directional_derivative(
                                                                          this->dependent_variables[i], j));
        }

      return jacobian;
    }


  } // namespace AD
} // namespace Differentiation


/* --- Explicit instantiations --- */
#ifdef DEAL_II_WITH_ADOLC
#include "ad_helpers.inst1"
#endif
#ifdef DEAL_II_WITH_TRILINOS
#include "ad_helpers.inst2"
#endif


DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_WITH_TRILINOS)
