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


// Test that some cases that we expect to fail do in fact do so:
// Overlapping indices of independent variables

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
test_overlapping_ind_var_components()
{
  using ADHelper     = AD::ScalarFunction<dim, ad_type_code, number_t>;
  using ADNumberType = typename ADHelper::ad_type;

  std::cout << "*** Test variables: Vector + Scalar (coupled), "
            << (AD::ADNumberTraits<ADNumberType>::is_taped == true ? "Taped" :
                                                                     "Tapeless")
            << std::endl;

  try
    {
      // NOTE: Second extractor overlaps first ont.
      const FEValuesExtractors::Vector v1_dof(0);
      const FEValuesExtractors::Scalar s2_dof(1);
      const unsigned int               n_AD_components =
        Tensor<1, dim>::n_independent_components + 1;
      ADHelper ad_helper(n_AD_components);
      ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                         // default values

      Tensor<1, dim> v1;
      double         s2 = 5.9;

      // Configure tape
      const int  tape_no = 1;
      const bool is_recording =
        ad_helper.start_recording_operations(tape_no /*material_id*/,
                                             true /*overwrite_tape*/,
                                             true /*keep*/);

      ad_helper.register_independent_variable(v1, v1_dof);
      ad_helper.register_independent_variable(s2, s2_dof);

      const Tensor<1, dim, ADNumberType> v1_ad =
        ad_helper.get_sensitive_variables(v1_dof);
      const ADNumberType s2_ad = ad_helper.get_sensitive_variables(s2_dof);

      const ADNumberType psi((v1_ad * v1_ad) * s2_ad);

      ad_helper.register_dependent_variable(psi);
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
  bool expected_result;
  if constexpr (running_in_debug_mode())
    {
      // Asserts should be triggered
      expected_result = false;
    }
  else
    {
      // User beware: Asserts ignored
      expected_result = true;
    }

  const unsigned int dim = 2;
  AssertThrow(
    (test_overlapping_ind_var_components<dim,
                                         double,
                                         AD::NumberTypes::adolc_taped>() ==
     expected_result),
    ExcInternalError());
  AssertThrow(
    (test_overlapping_ind_var_components<dim,
                                         double,
                                         AD::NumberTypes::adolc_tapeless>() ==
     expected_result),
    ExcInternalError());

  deallog << "OK" << std::endl;
}
