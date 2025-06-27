// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/types.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/differentiation/ad/ad_drivers.h>
#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/ad_number_types.h>
#  include <deal.II/differentiation/ad/adolc_number_types.h>
#  include <deal.II/differentiation/ad/sacado_number_types.h>

#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/vector.h>

#  ifdef DEAL_II_WITH_ADOLC
#    include <adolc/adolc_fatalerror.h>
#    include <adolc/drivers/drivers.h>
#    include <adolc/taping.h>
#  endif // DEAL_II_WITH_ADOLC

#  include <vector>


DEAL_II_NAMESPACE_OPEN


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
    std::vector<typename Types<ADNumberType>::tape_index>
    TapedDrivers<ADNumberType, ScalarType, T>::get_registered_tape_indices()
      const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return std::vector<typename Types<ADNumberType>::tape_index>();
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::activate_tape(
      const typename Types<ADNumberType>::tape_index)
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    bool
    TapedDrivers<ADNumberType, ScalarType, T>::requires_retaping(
      const typename Types<ADNumberType>::tape_index) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return false;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    bool
    TapedDrivers<ADNumberType, ScalarType, T>::last_action_requires_retaping()
      const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return false;
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::remove_tape(
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
      std::ostream &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    ScalarType
    TapedDrivers<ADNumberType, ScalarType, T>::value(
      const typename Types<ADNumberType>::tape_index,
      const std::vector<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return ScalarType(0.0);
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::gradient(
      const typename Types<ADNumberType>::tape_index,
      const std::vector<ScalarType> &,
      Vector<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::hessian(
      const typename Types<ADNumberType>::tape_index,
      const std::vector<ScalarType> &,
      FullMatrix<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::values(
      const typename Types<ADNumberType>::tape_index,
      const unsigned int,
      const std::vector<ScalarType> &,
      Vector<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapedDrivers<ADNumberType, ScalarType, T>::jacobian(
      const typename Types<ADNumberType>::tape_index,
      const unsigned int,
      const std::vector<ScalarType> &,
      FullMatrix<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }



#  ifdef DEAL_II_WITH_ADOLC

#    ifndef DOXYGEN
    // Specialization for taped ADOL-C auto-differentiable numbers.

    template <typename ADNumberType>
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::TapedDrivers()
      : active_tape(Numbers<ADNumberType>::invalid_tape_index)
      , keep_values(true)
      , is_recording_flag(false)
      , use_stored_taped_buffer_sizes(false)
      , obufsize(0u)
      , lbufsize(0u)
      , vbufsize(0u)
      , tbufsize(0u)
    {}
#    endif


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::is_recording()
      const
    {
      return is_recording_flag;
    }


    template <typename ADNumberType>
    typename Types<ADNumberType>::tape_index
    TapedDrivers<
      ADNumberType,
      double,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_taped>>::active_tape_index() const
    {
      return active_tape;
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      keep_independent_values() const
    {
      return keep_values;
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const
    {
      // Sigh... This is a mess :-/
      // The most succinct way to get this piece of information, would be to
      // use the getTapeInfos() function, but would come at the expense of
      // creating an inactive tape data within ADOL-C's global store. For
      // hints as to why this is the way it is, see the getTapeInfos()
      // function in
      // https://gitlab.com/adol-c/adol-c/blob/master/ADOL-C/src/tape_handling.cpp
      // An alternative solution would be to manually access the tape data;
      // see the removeTape() function in
      // https://gitlab.com/adol-c/adol-c/blob/master/ADOL-C/src/tape_handling.cpp
      // with the consideration of the #defines in
      // https://gitlab.com/adol-c/adol-c/blob/master/ADOL-C/src/taping_p.h
      // as to how this would be performed.
      // Doing things "manually" (the second way) without creating the
      // additional data object would be a lot more work...
      //
      // Both of the above solutions would be possible IF ADOL-C exposed
      // this data object or a method to access it to the outside world.
      // But they don't :-(
      // Instead, what we'll have to do is get the statistics for this tape
      // and make our own determination as to whether or not this tape exists.
      // This effectively executes the first solution, with even more
      // overhead! If the tape is in existence, we can assume that it should a
      // non-zero number of dependent and independent variables. Those are
      // stored in the zeroth and first entries of the statistics vector.
      //
      // But, oh wait... Surprise! This will trigger an error if the tape
      // doesn't exist at all! So lets first check their tape info cache to
      // see if the tape REALLY exists (i.e. has been touched, even if nothing
      // has been written to it) before trying to access it. It'll only take
      // an O(n) search, but at this point who really cares about efficiency?
      //
      // It has been suggested in
      // https://gitlab.com/adol-c/adol-c/issues/11
      // that a simply try-catch block around ::tapestats is the solution that
      // we want here. Unfortunately this results in unwanted pollution of
      // the terminal, of the form
      // ADOL-C error: reading integer tape number 4!
      //               >>> File or directory not found! <<<
      // , every time a query is made about a non-existent tape.
      // So either way we have to guard that check with something more
      // conservative so that we don't output useless messages for our users.
      const std::vector<typename Types<ADNumberType>::tape_index>
                 registered_tape_indices = get_registered_tape_indices();
      const auto it = std::find(registered_tape_indices.begin(),
                                registered_tape_indices.end(),
                                tape_index);
      if (it == registered_tape_indices.end())
        return false;

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


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
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


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
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


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
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


    template <typename ADNumberType>
    std::vector<typename Types<ADNumberType>::tape_index>
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      get_registered_tape_indices() const
    {
      // We've chosen to use unsigned shorts for the tape
      // index type (a safety precaution) so we need to
      // perform a conversion between ADOL-C's native tape
      // index type and that chosen by us.
      std::vector<short> registered_tape_indices_s;
      cachedTraceTags(registered_tape_indices_s);

      return std::vector<typename Types<ADNumberType>::tape_index>(
        registered_tape_indices_s.begin(), registered_tape_indices_s.end());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index)
    {
      active_tape = tape_index;
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      requires_retaping(
        const typename Types<ADNumberType>::tape_index tape_index) const
    {
      Assert(
        is_registered_tape(tape_index) == true,
        ExcMessage(
          "Cannot ask for the status of a tape that has not yet been recorded and used."));

      const auto it_status_tape = status.find(tape_index);

      // This tape's status has not been found in the map. This could be
      // because a non-existent tape_index has been used as an argument, or
      // because the tape exists but has not been used (in a way that
      // initiated a status update). For example, on can create a tape in one
      // section of code and then query the status of all existing tapes in a
      // completely different section of code that knows nothing about the
      // first tape. There is no prerequisite that the first tape is ever
      // used, and it can therefore not have a status. So, in this case
      // there's not much we can do other than to report that the tape does
      // not require retaping.
      if (it_status_tape == status.end())
        return false;

      const auto status_tape = it_status_tape->second;

      // See ADOL-C manual section 1.7 and comments in last paragraph of
      // section 3.1
      Assert(
        status_tape < 4 && status_tape >= -2,
        ExcMessage(
          "The tape status is not within the range specified within the ADOL-C documentation."));
      return (status_tape < 0);
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      last_action_requires_retaping() const
    {
      return requires_retaping(active_tape);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      remove_tape(const typename Types<ADNumberType>::tape_index tape_index)
    {
      Assert(is_registered_tape(tape_index),
             ExcMessage(
               "This tape does not exist, and therefore cannot be cleared."));
      removeTape(tape_index, TapeRemovalType::ADOLC_REMOVE_COMPLETELY);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      reset(const bool clear_registered_tapes)
    {
      active_tape       = Numbers<ADNumberType>::invalid_tape_index;
      is_recording_flag = false;
      status.clear();
      if (clear_registered_tapes)
        {
          const std::vector<typename Types<ADNumberType>::tape_index>
            registered_tape_indices = get_registered_tape_indices();
          for (const auto &tape_index : registered_tape_indices)
            remove_tape(tape_index);
        }
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::print(std::ostream
                                                                      &stream)
      const
    {
      const std::vector<typename Types<ADNumberType>::tape_index>
        registered_tape_indices = get_registered_tape_indices();
      stream << "Registered tapes and their status: ";
      auto it_registered_tape = registered_tape_indices.begin();
      for (unsigned int i = 0; i < registered_tape_indices.size();
           ++i, ++it_registered_tape)
        {
          const auto tape_index     = *it_registered_tape;
          const auto it_status_tape = status.find(tape_index);
          Assert(it_status_tape != status.end(),
                 ExcMessage(
                   "This tape's status has not been found in the map."));
          const auto status_tape = it_status_tape->second;

          stream << tape_index << "->" << status_tape
                 << (i < (registered_tape_indices.size() - 1) ? "," : "");
        }
      stream << '\n';

      stream << "Keep values? " << keep_independent_values() << '\n';
      stream << "Use stored tape buffer sizes? "
             << use_stored_taped_buffer_sizes << '\n';
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream                                  &stream) const
    {
      // See ADOL-C manual section 2.1
      // and adolc/taping.h
      std::vector<std::size_t> counts(STAT_SIZE);
      ::tapestats(tape_index, counts.data());
      Assert(counts.size() >= 18, ExcInternalError());
      stream
        << "Tape index: " << tape_index << '\n'
        << "Number of independent variables: " << counts[0] << '\n'
        << "Number of dependent variables:   " << counts[1] << '\n'
        << "Max number of live, active variables: " << counts[2] << '\n'
        << "Size of taylor stack (number of overwrites): " << counts[3] << '\n'
        << "Operations buffer size: " << counts[4] << '\n'
        << "Total number of recorded operations: " << counts[5] << '\n'
        << "Operations file written or not: " << counts[6] << '\n'
        << "Overall number of locations: " << counts[7] << '\n'
        << "Locations file written or not: " << counts[8] << '\n'
        << "Overall number of values: " << counts[9] << '\n'
        << "Values file written or not: " << counts[10] << '\n'
        << "Locations buffer size: " << counts[11] << '\n'
        << "Values buffer size: " << counts[12] << '\n'
        << "Taylor buffer size: " << counts[13] << '\n'
        << "Number of eq_*_prod for sparsity pattern: " << counts[14] << '\n'
        << "Use of 'min_op', deferred to 'abs_op' for piecewise calculations: "
        << counts[15] << '\n'
        << "Number of 'abs' calls that can switch branch: " << counts[16]
        << '\n'
        << "Number of parameters (doubles) interchangeable without retaping: "
        << counts[17] << '\n'
        << std::flush;
    }


    template <typename ADNumberType>
    typename TapedDrivers<
      ADNumberType,
      double,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_taped>>::scalar_type
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<scalar_type> &independent_variables) const
    {
      Assert(is_registered_tape(active_tape_index),
             ExcMessage("This tape has not yet been recorded."));

      scalar_type value = 0.0;

      status[active_tape_index] =
        ::function(active_tape_index,
                   1, // Only one dependent variable
                   independent_variables.size(),
                   const_cast<double *>(independent_variables.data()),
                   &value);

      return value;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<scalar_type> &independent_variables,
               Vector<scalar_type>            &gradient) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               1,
             ExcSupportedDerivativeLevels(
               AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
               1));
      Assert(gradient.size() == independent_variables.size(),
             ExcDimensionMismatch(gradient.size(),
                                  independent_variables.size()));
      Assert(is_registered_tape(active_tape_index),
             ExcMessage("This tape has not yet been recorded."));

      // Note: ADOL-C's ::gradient function expects a *double as the last
      // parameter. Here we take advantage of the fact that the data in the
      // Vector class is aligned (e.g. stored as an Array)
      status[active_tape_index] =
        ::gradient(active_tape_index,
                   independent_variables.size(),
                   const_cast<scalar_type *>(independent_variables.data()),
                   gradient.data());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<scalar_type> &independent_variables,
              FullMatrix<scalar_type>        &hessian) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               2,
             ExcSupportedDerivativeLevels(
               AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
               2));
      Assert(hessian.m() == independent_variables.size(),
             ExcDimensionMismatch(hessian.m(), independent_variables.size()));
      Assert(hessian.n() == independent_variables.size(),
             ExcDimensionMismatch(hessian.n(), independent_variables.size()));
      Assert(is_registered_tape(active_tape_index),
             ExcMessage("This tape has not yet been recorded."));

      const unsigned int n_independent_variables = independent_variables.size();
      std::vector<scalar_type *> H(n_independent_variables);
      for (unsigned int i = 0; i < n_independent_variables; ++i)
        H[i] = &hessian[i][0];

      status[active_tape_index] =
        ::hessian(active_tape_index,
                  n_independent_variables,
                  const_cast<scalar_type *>(independent_variables.data()),
                  H.data());

      // ADOL-C builds only the lower-triangular part of the
      // symmetric Hessian, so we should copy the relevant
      // entries into the upper triangular part.
      for (unsigned int i = 0; i < n_independent_variables; ++i)
        for (unsigned int j = 0; j < i; ++j)
          hessian[j][i] = hessian[i][j]; // Symmetry
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int              n_dependent_variables,
             const std::vector<scalar_type> &independent_variables,
             Vector<scalar_type>            &values) const
    {
      Assert(values.size() == n_dependent_variables,
             ExcDimensionMismatch(values.size(), n_dependent_variables));
      Assert(is_registered_tape(active_tape_index),
             ExcMessage("This tape has not yet been recorded."));

      // Note: ADOL-C's ::function function expects a *double as the last
      // parameter. Here we take advantage of the fact that the data in the
      // Vector class is aligned (e.g. stored as an Array)
      status[active_tape_index] =
        ::function(active_tape_index,
                   n_dependent_variables,
                   independent_variables.size(),
                   const_cast<scalar_type *>(independent_variables.data()),
                   values.data());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int              n_dependent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type>        &jacobian) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               1,
             ExcSupportedDerivativeLevels(
               AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
               1));
      Assert(jacobian.m() == n_dependent_variables,
             ExcDimensionMismatch(jacobian.m(), n_dependent_variables));
      Assert(jacobian.n() == independent_variables.size(),
             ExcDimensionMismatch(jacobian.n(), independent_variables.size()));
      Assert(is_registered_tape(active_tape_index),
             ExcMessage("This tape has not yet been recorded."));

      std::vector<scalar_type *> J(n_dependent_variables);
      for (unsigned int i = 0; i < n_dependent_variables; ++i)
        J[i] = &jacobian[i][0];

      status[active_tape_index] = ::jacobian(active_tape_index,
                                             n_dependent_variables,
                                             independent_variables.size(),
                                             independent_variables.data(),
                                             J.data());
    }

#  else

    // Specialization for taped ADOL-C auto-differentiable numbers.

    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::is_recording()
      const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return false;
    }


    template <typename ADNumberType>
    typename Types<ADNumberType>::tape_index
    TapedDrivers<
      ADNumberType,
      double,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_taped>>::active_tape_index() const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return Numbers<ADNumberType>::invalid_tape_index;
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      keep_independent_values() const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return false;
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      is_registered_tape(const typename Types<ADNumberType>::tape_index) const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return false;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      set_tape_buffer_sizes(
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes,
        const typename Types<ADNumberType>::tape_buffer_sizes)
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      start_taping(const typename Types<ADNumberType>::tape_index, const bool)
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      stop_taping(const typename Types<ADNumberType>::tape_index, const bool)
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    std::vector<typename Types<ADNumberType>::tape_index>
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      get_registered_tape_indices() const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return std::vector<typename Types<ADNumberType>::tape_index>();
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      activate_tape(const typename Types<ADNumberType>::tape_index)
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      requires_retaping(const typename Types<ADNumberType>::tape_index) const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return false;
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      last_action_requires_retaping() const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return false;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      remove_tape(const typename Types<ADNumberType>::tape_index)
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::reset(const bool)
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::print(std::ostream
                                                                      &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      print_tape_stats(const typename Types<ADNumberType>::tape_index,
                       std::ostream &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    typename TapedDrivers<
      ADNumberType,
      double,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_taped>>::scalar_type
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      value(const typename Types<ADNumberType>::tape_index,
            const std::vector<scalar_type> &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
      return 0.0;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      gradient(const typename Types<ADNumberType>::tape_index,
               const std::vector<scalar_type> &,
               Vector<scalar_type> &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      hessian(const typename Types<ADNumberType>::tape_index,
              const std::vector<scalar_type> &,
              FullMatrix<scalar_type> &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      values(const typename Types<ADNumberType>::tape_index,
             const unsigned int,
             const std::vector<scalar_type> &,
             Vector<scalar_type> &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 double,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      jacobian(const typename Types<ADNumberType>::tape_index,
               const unsigned int,
               const std::vector<scalar_type> &,
               FullMatrix<scalar_type> &) const
    {
      AssertThrow(false, ExcRequiresADOLC());
    }

#  endif // DEAL_II_WITH_ADOLC


    // Specialization for ADOL-C taped numbers. It is expected that the
    // scalar return type for this class is a float.

    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::is_recording()
      const
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      return taped_driver.is_recording();
    }


    template <typename ADNumberType>
    typename Types<ADNumberType>::tape_index
    TapedDrivers<
      ADNumberType,
      float,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_taped>>::active_tape_index() const
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      return taped_driver.active_tape_index();
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      keep_independent_values() const
    {
      return taped_driver.keep_independent_values();
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      is_registered_tape(
        const typename Types<ADNumberType>::tape_index tape_index) const
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      return taped_driver.is_registered_tape(tape_index);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
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


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      start_taping(const typename Types<ADNumberType>::tape_index tape_index,
                   const bool keep_independent_values)
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.start_taping(tape_index, keep_independent_values);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      stop_taping(
        const typename Types<ADNumberType>::tape_index active_tape_index,
        const bool                                     write_tapes_to_file)
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.stop_taping(active_tape_index, write_tapes_to_file);
    }


    template <typename ADNumberType>
    std::vector<typename Types<ADNumberType>::tape_index>
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      get_registered_tape_indices() const
    {
      return taped_driver.get_registered_tape_indices();
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      activate_tape(const typename Types<ADNumberType>::tape_index tape_index)
    {
      taped_driver.activate_tape(tape_index);
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      requires_retaping(
        const typename Types<ADNumberType>::tape_index tape_index) const
    {
      return taped_driver.requires_retaping(tape_index);
    }


    template <typename ADNumberType>
    bool
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      last_action_requires_retaping() const
    {
      return taped_driver.last_action_requires_retaping();
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      remove_tape(const typename Types<ADNumberType>::tape_index tape_index)
    {
      taped_driver.remove_tape(tape_index);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      reset(const bool clear_registered_tapes)
    {
      taped_driver.reset(clear_registered_tapes);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::print(std::ostream
                                                                      &stream)
      const
    {
      taped_driver.print(stream);
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      print_tape_stats(
        const typename Types<ADNumberType>::tape_index tape_index,
        std::ostream                                  &stream) const
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.print_tape_stats(tape_index, stream);
    }


    template <typename ADNumberType>
    typename TapedDrivers<
      ADNumberType,
      float,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_taped>>::scalar_type
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      value(const typename Types<ADNumberType>::tape_index active_tape_index,
            const std::vector<scalar_type> &independent_variables) const
    {
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      return taped_driver.value(active_tape_index,
                                vector_float_to_double(independent_variables));
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      gradient(const typename Types<ADNumberType>::tape_index active_tape_index,
               const std::vector<scalar_type> &independent_variables,
               Vector<scalar_type>            &gradient) const
    {
      Vector<double> gradient_double(gradient.size());
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.gradient(active_tape_index,
                            vector_float_to_double(independent_variables),
                            gradient_double);
      gradient = gradient_double;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      hessian(const typename Types<ADNumberType>::tape_index active_tape_index,
              const std::vector<scalar_type> &independent_variables,
              FullMatrix<scalar_type>        &hessian) const
    {
      FullMatrix<double> hessian_double(hessian.m(), hessian.n());
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.hessian(active_tape_index,
                           vector_float_to_double(independent_variables),
                           hessian_double);
      hessian = hessian_double;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      values(const typename Types<ADNumberType>::tape_index active_tape_index,
             const unsigned int              n_dependent_variables,
             const std::vector<scalar_type> &independent_variables,
             Vector<scalar_type>            &values) const
    {
      Vector<double> values_double(values.size());
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.values(active_tape_index,
                          n_dependent_variables,
                          vector_float_to_double(independent_variables),
                          values_double);
      values = values_double;
    }


    template <typename ADNumberType>
    void
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      jacobian(const typename Types<ADNumberType>::tape_index active_tape_index,
               const unsigned int              n_dependent_variables,
               const std::vector<scalar_type> &independent_variables,
               FullMatrix<scalar_type>        &jacobian) const
    {
      FullMatrix<double> jacobian_double(jacobian.m(), jacobian.n());
      // ADOL-C only supports 'double', not 'float', so we can forward to
      // the 'double' implementation of this function
      taped_driver.jacobian(active_tape_index,
                            n_dependent_variables,
                            vector_float_to_double(independent_variables),
                            jacobian_double);
      jacobian = jacobian_double;
    }


#  ifndef DOXYGEN
    template <typename ADNumberType>
    std::vector<double>
    TapedDrivers<ADNumberType,
                 float,
                 std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                  NumberTypes::adolc_taped>>::
      vector_float_to_double(const std::vector<float> &in) const
    {
      std::vector<double> out(in.size());
      std::copy(in.begin(), in.end(), out.begin());
      return out;
    }
#  endif


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
      const std::vector<ADNumberType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
      return ScalarType(0.0);
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::gradient(
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      Vector<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::hessian(
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      FullMatrix<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::values(
      const std::vector<ADNumberType> &,
      Vector<ScalarType> &) const
    {
      AssertThrow(false, ExcRequiresADNumberSpecialization());
    }


    template <typename ADNumberType, typename ScalarType, typename T>
    void
    TapelessDrivers<ADNumberType, ScalarType, T>::jacobian(
      const std::vector<ADNumberType> &,
      const std::vector<ADNumberType> &,
      FullMatrix<ScalarType> &) const
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
      std::enable_if_t<!(ADNumberTraits<ADNumberType>::type_code ==
                           NumberTypes::sacado_rad ||
                         ADNumberTraits<ADNumberType>::type_code ==
                           NumberTypes::sacado_rad_dfad)>
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
      std::enable_if_t<
        ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_rad ||
        ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_rad_dfad>
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
      std::enable_if_t<!(ADNumberTraits<ADNumberType>::type_code ==
                         NumberTypes::adolc_tapeless)>
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
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_tapeless>
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
                  "There are currently " + std::to_string(n_live_variables) +
                  " live "
                  "adtl::adouble variables in existence. They currently "
                  "assume " +
                  std::to_string(n_set_directional_derivatives) +
                  " directional derivatives "
                  "but you wish to increase this to " +
                  std::to_string(n_directional_derivatives) +
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
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                       NumberTypes::adolc_tapeless>
      configure_tapeless_mode(const unsigned int /*n_directional_derivatives*/)
      {
        AssertThrow(false, ExcRequiresADOLC());
      }

#  endif

    } // namespace internal



    // Specialization for auto-differentiable numbers that use
    // reverse mode to compute the first derivatives (and, if supported,
    // forward mode for the second).

#  ifndef DOXYGEN
    template <typename ADNumberType, typename ScalarType>
    TapelessDrivers<
      ADNumberType,
      ScalarType,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                         NumberTypes::sacado_rad ||
                       ADNumberTraits<ADNumberType>::type_code ==
                         NumberTypes::sacado_rad_dfad>>::TapelessDrivers()
      : dependent_variable_marking_safe(false)
    {}
#  endif


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      initialize_global_environment(const unsigned int n_independent_variables)
    {
      internal::configure_tapeless_mode<ADNumberType>(n_independent_variables);
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<
      ADNumberType,
      ScalarType,
      std::enable_if_t<
        ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_rad ||
        ADNumberTraits<ADNumberType>::type_code ==
          NumberTypes::sacado_rad_dfad>>::allow_dependent_variable_marking()
    {
      dependent_variable_marking_safe = true;
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<
      ADNumberType,
      ScalarType,
      std::enable_if_t<
        ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_rad ||
        ADNumberTraits<ADNumberType>::type_code ==
          NumberTypes::sacado_rad_dfad>>::prevent_dependent_variable_marking()
    {
      dependent_variable_marking_safe = false;
    }


    template <typename ADNumberType, typename ScalarType>
    bool
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      is_dependent_variable_marking_allowed() const
    {
      return dependent_variable_marking_safe;
    }


    template <typename ADNumberType, typename ScalarType>
    ScalarType
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      value(const std::vector<ADNumberType> &dependent_variables) const
    {
      Assert(dependent_variables.size() == 1,
             ExcDimensionMismatch(dependent_variables.size(), 1));
      return ADNumberTraits<ADNumberType>::get_scalar_value(
        dependent_variables[0]);
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType>              &gradient) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               1,
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
      const std::size_t n_independent_variables = independent_variables.size();
      for (unsigned int i = 0; i < n_independent_variables; ++i)
        gradient[i] = internal::NumberType<ScalarType>::value(
          ADNumberTraits<ADNumberType>::get_directional_derivative(
            independent_variables[i], 0 /*This number doesn't really matter*/));
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType>          &hessian) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               2,
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
      const std::size_t n_independent_variables = independent_variables.size();
      for (unsigned int i = 0; i < n_independent_variables; ++i)
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


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType>              &values) const
    {
      Assert(values.size() == dependent_variables.size(),
             ExcDimensionMismatch(values.size(), dependent_variables.size()));

      const std::size_t n_dependent_variables = dependent_variables.size();
      for (unsigned int i = 0; i < n_dependent_variables; ++i)
        values[i] = ADNumberTraits<ADNumberType>::get_scalar_value(
          dependent_variables[i]);
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_rad_dfad>>::
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType>          &jacobian) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               1,
             ExcSupportedDerivativeLevels(
               AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
               1));
      Assert(jacobian.m() == dependent_variables.size(),
             ExcDimensionMismatch(jacobian.m(), dependent_variables.size()));
      Assert(jacobian.n() == independent_variables.size(),
             ExcDimensionMismatch(jacobian.n(), independent_variables.size()));

      const std::size_t n_independent_variables = independent_variables.size();
      const std::size_t n_dependent_variables   = dependent_variables.size();

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

      for (unsigned int i = 0; i < n_dependent_variables; ++i)
        {
          internal::reverse_mode_dependent_variable_activation(
            const_cast<ADNumberType &>(dependent_variables[i]));
          for (unsigned int j = 0; j < n_independent_variables; ++j)
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



    // Specialization for auto-differentiable numbers that use
    // forward mode to compute the first (and, if supported, second)
    // derivatives.

#  ifndef DOXYGEN
    template <typename ADNumberType, typename ScalarType>
    TapelessDrivers<
      ADNumberType,
      ScalarType,
      std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                         NumberTypes::adolc_tapeless ||
                       ADNumberTraits<ADNumberType>::type_code ==
                         NumberTypes::sacado_dfad ||
                       ADNumberTraits<ADNumberType>::type_code ==
                         NumberTypes::sacado_dfad_dfad>>::TapelessDrivers()
      : dependent_variable_marking_safe(false)
    {}
#  endif


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      initialize_global_environment(const unsigned int n_independent_variables)
    {
      internal::configure_tapeless_mode<ADNumberType>(n_independent_variables);
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<
      ADNumberType,
      ScalarType,
      std::enable_if_t<
        ADNumberTraits<ADNumberType>::type_code ==
          NumberTypes::adolc_tapeless ||
        ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_dfad ||
        ADNumberTraits<ADNumberType>::type_code ==
          NumberTypes::sacado_dfad_dfad>>::allow_dependent_variable_marking()
    {
      dependent_variable_marking_safe = true;
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<
      ADNumberType,
      ScalarType,
      std::enable_if_t<
        ADNumberTraits<ADNumberType>::type_code ==
          NumberTypes::adolc_tapeless ||
        ADNumberTraits<ADNumberType>::type_code == NumberTypes::sacado_dfad ||
        ADNumberTraits<ADNumberType>::type_code ==
          NumberTypes::sacado_dfad_dfad>>::prevent_dependent_variable_marking()
    {
      dependent_variable_marking_safe = false;
    }


    template <typename ADNumberType, typename ScalarType>
    bool
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      is_dependent_variable_marking_allowed() const
    {
      return dependent_variable_marking_safe;
    }


    template <typename ADNumberType, typename ScalarType>
    ScalarType
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      value(const std::vector<ADNumberType> &dependent_variables) const
    {
      Assert(dependent_variables.size() == 1,
             ExcDimensionMismatch(dependent_variables.size(), 1));
      return ADNumberTraits<ADNumberType>::get_scalar_value(
        dependent_variables[0]);
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      gradient(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               Vector<ScalarType>              &gradient) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               1,
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
      const std::size_t n_independent_variables = independent_variables.size();
      for (unsigned int i = 0; i < n_independent_variables; ++i)
        gradient[i] = internal::NumberType<ScalarType>::value(
          ADNumberTraits<ADNumberType>::get_directional_derivative(
            dependent_variables[0], i));
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      hessian(const std::vector<ADNumberType> &independent_variables,
              const std::vector<ADNumberType> &dependent_variables,
              FullMatrix<ScalarType>          &hessian) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               2,
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
      const std::size_t n_independent_variables = independent_variables.size();
      for (unsigned int i = 0; i < n_independent_variables; ++i)
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


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      values(const std::vector<ADNumberType> &dependent_variables,
             Vector<ScalarType>              &values) const
    {
      Assert(values.size() == dependent_variables.size(),
             ExcDimensionMismatch(values.size(), dependent_variables.size()));

      const std::size_t n_dependent_variables = dependent_variables.size();
      for (unsigned int i = 0; i < n_dependent_variables; ++i)
        values[i] = ADNumberTraits<ADNumberType>::get_scalar_value(
          dependent_variables[i]);
    }


    template <typename ADNumberType, typename ScalarType>
    void
    TapelessDrivers<ADNumberType,
                    ScalarType,
                    std::enable_if_t<ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::adolc_tapeless ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad ||
                                     ADNumberTraits<ADNumberType>::type_code ==
                                       NumberTypes::sacado_dfad_dfad>>::
      jacobian(const std::vector<ADNumberType> &independent_variables,
               const std::vector<ADNumberType> &dependent_variables,
               FullMatrix<ScalarType>          &jacobian) const
    {
      Assert(AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >=
               1,
             ExcSupportedDerivativeLevels(
               AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels,
               1));
      Assert(jacobian.m() == dependent_variables.size(),
             ExcDimensionMismatch(jacobian.m(), dependent_variables.size()));
      Assert(jacobian.n() == independent_variables.size(),
             ExcDimensionMismatch(jacobian.n(), independent_variables.size()));

      const std::size_t n_independent_variables = independent_variables.size();
      const std::size_t n_dependent_variables   = dependent_variables.size();

      // In forward mode, the gradients are computed from the
      // dependent variables
      for (unsigned int i = 0; i < n_dependent_variables; ++i)
        for (unsigned int j = 0; j < n_independent_variables; ++j)
          jacobian[i][j] = internal::NumberType<ScalarType>::value(
            ADNumberTraits<ADNumberType>::get_directional_derivative(
              dependent_variables[i], j));
    }


  } // namespace AD
} // namespace Differentiation


/* --- Explicit instantiations --- */
// We don't build the .inst files if deal.II isn't configured with the
// external dependencies, but doxygen doesn't know that and tries to
// find that file anyway for parsing -- which then of course it fails
// on. So exclude the following from doxygen consideration.
#  ifndef DOXYGEN
#    include "differentiation/ad/ad_drivers.inst"
#    ifdef DEAL_II_WITH_ADOLC
#      include "differentiation/ad/ad_drivers.inst1"
#    endif
#    ifdef DEAL_II_TRILINOS_WITH_SACADO
#      include "differentiation/ad/ad_drivers.inst2"
#    endif
#  endif


DEAL_II_NAMESPACE_CLOSE


#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_TRILINOS_WITH_SACADO)
