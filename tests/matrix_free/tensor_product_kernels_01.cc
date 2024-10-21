// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests the evaluate_general variant of the apply_matrix_vector_product
// function where either one, two or three vectors are used at the same time

#include <deal.II/matrix_free/tensor_product_kernels.h>

#include "../tests.h"

template <bool transpose, bool add>
void
test()
{
  const unsigned int n_columns = 20;
  const unsigned int n_rows    = 20;
  double             matrix[n_rows * n_columns];
  for (unsigned int i = 0; i < n_columns * n_rows; ++i)
    matrix[i] = i;



  double in0[n_columns];
  double in1[n_columns];
  double in2[n_columns];
  for (unsigned int i = 0; i < n_columns; ++i)
    {
      in0[i] = i;
      in1[i] = 2 * i;
      in2[i] = 3 * i;
    }


  // Reference solutions
  double out0_0[n_rows] = {0};
  double out0_1[n_rows] = {1};
  double out0_2[n_rows] = {2};

  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false>(matrix, in0, out0_0, n_rows, n_columns, 1, 1);
  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false>(matrix, in1, out0_1, n_rows, n_columns, 1, 1);
  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false>(matrix, in2, out0_2, n_rows, n_columns, 1, 1);

  // Compute 2 vectors at once
  double out1_0[n_rows] = {0};
  double out1_1[n_rows] = {1};
  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false>(
    matrix, in0, in1, out1_0, out1_1, n_rows, n_columns, 1, 1);

  // Compute 3 vectors at once
  double out2_0[n_rows] = {0};
  double out2_1[n_rows] = {1};
  double out2_2[n_rows] = {2};
  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false>(
    matrix, in0, in1, in2, out2_0, out2_1, out2_2, n_rows, n_columns, 1, 1);


  // Check if outputs are the same
  deallog << "Check results" << std::endl;
  for (unsigned int i = 0; i < n_rows; ++i)
    {
      if ((out0_0[i] != out1_0[i]) || (out0_1[i] != out1_1[i]))
        deallog << "Error 2 components at entry " << i << std::endl;
      else if ((out0_0[i] != out2_0[i]) || (out0_1[i] != out2_1[i]) ||
               (out0_2[i] != out2_2[i]))
        deallog << "Error 3 components at entry " << i << std::endl;
      else
        deallog << "Correct line " << i << std::endl;
    }
  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(9);
  test<false, false>();
  test<true, false>();
  test<false, true>();
  test<true, true>();
}
