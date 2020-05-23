
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// Test that some cases that we expect to fail do in fact do so:
// Registering more than one dependent variable

#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"

namespace AD = dealii::Differentiation::AD;

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
bool
test_register_multiple_dep_vars()
{
  typedef AD::ScalarFunction<dim, ad_type_code, number_t> ADHelper;
  typedef typename ADHelper::ad_type                      ADNumberType;

  std::cout << "*** Test variables: Scalar + Scalar (coupled), "
            << (AD::ADNumberTraits<ADNumberType>::is_taped == true ? "Taped" :
                                                                     "Tapeless")
            << std::endl;

  try
    {
      const FEValuesExtractors::Scalar s1_dof(0);
      const FEValuesExtractors::Scalar s2_dof(1);
      const unsigned int               n_AD_components = 2;
      ADHelper                         ad_helper(n_AD_components);
      ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                         // default values

      double s1 = 3.1;
      double s2 = 5.9;

      // Configure tape
      const int  tape_no = 1;
      const bool is_recording =
        ad_helper.start_recording_operations(tape_no /*material_id*/,
                                             true /*overwrite_tape*/,
                                             true /*keep*/);

      ad_helper.register_independent_variable(s1, s1_dof);
      ad_helper.register_independent_variable(s2, s2_dof);

      const ADNumberType s1_ad = ad_helper.get_sensitive_variables(s1_dof);
      const ADNumberType s2_ad = ad_helper.get_sensitive_variables(s2_dof);

      const ADNumberType psi(s1_ad * s2_ad);

      // NOTE: Should only be calling this function once for this helper class
      ad_helper.register_dependent_variable(psi);
      ad_helper.register_dependent_variable(psi); // Wrong!

      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Unexpected success!" << std::endl;
      return true;
    }
  catch (...)
    {
      std::cout << "Expected failure." << std::endl;
      return false;
    }
}

int
main()
{
  initlog();

  deal_II_exceptions::disable_abort_on_exception();
#if DEBUG
  // Asserts should be triggered
  const bool expected_result = false;
#else
  // User beware: Asserts ignored
  const bool expected_result = true;
#endif

  const unsigned int dim = 2;
  AssertThrow(
    (test_register_multiple_dep_vars<dim,
                                     double,
                                     AD::NumberTypes::adolc_taped>() ==
     expected_result),
    ExcInternalError());
  AssertThrow(
    (test_register_multiple_dep_vars<dim,
                                     double,
                                     AD::NumberTypes::adolc_tapeless>() ==
     expected_result),
    ExcInternalError());

  deallog << "OK" << std::endl;
}
