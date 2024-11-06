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
  const unsigned int n_columns_raw = 20;
  const unsigned int n_rows_raw    = 105;

  const unsigned int n_columns = transpose ? n_rows_raw : n_columns_raw;
  const unsigned int n_rows    = transpose ? n_columns_raw : n_rows_raw;

  double matrix[n_rows * n_columns];
  for (unsigned int i = 0; i < n_columns * n_rows; ++i)
    matrix[i] = 0.1 * (i + 1);

  double in[n_columns];
  double in_batched_2[2 * n_columns];
  double in_batched[3 * n_columns];
  for (unsigned int i = 0; i < n_columns; ++i)
    {
      in[i]                         = i + 1;
      in_batched[i]                 = i + 1;
      in_batched[n_columns + i]     = i + 1;
      in_batched_2[i]               = i + 1;
      in_batched_2[n_columns + i]   = i + 1;
      in_batched[2 * n_columns + i] = i + 1;
    }


  // Reference solutions
  double out[n_rows];
  double out_batched_2[2 * n_rows];
  double out_batched[3 * n_rows];

  for (unsigned int i = 0; i < n_rows; ++i)
    {
      out[i]                      = i + 1;
      out_batched[i]              = i + 1;
      out_batched[n_rows + i]     = i + 1;
      out_batched_2[i]            = i + 1;
      out_batched_2[n_rows + i]   = i + 1;
      out_batched[2 * n_rows + i] = i + 1;
    }

  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*n_rows*/ n_rows_raw,
    /*n_columns*/ n_columns_raw,
    /*strid_in*/ 1,
    /*strid_out*/ 1,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    double,
    double>(matrix, in, out);

  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false,
    double,
    double,
    /*n_components*/ 3>(
    matrix, in_batched, out_batched, n_rows_raw, n_columns_raw, 1, 1);

  dealii::internal::apply_matrix_vector_product<
    dealii::internal::EvaluatorVariant::evaluate_general,
    dealii::internal::EvaluatorQuantity::value,
    /*transpose_matrix*/ transpose,
    /*add*/ add,
    /*consider_strides*/ false,
    double,
    double,
    /*n_components*/ 2>(
    matrix, in_batched_2, out_batched_2, n_rows_raw, n_columns_raw, 1, 1);


  // Check if outputs are the same
  deallog << "Check results for transpose = " << transpose
          << " and add = " << add << std::endl;
  for (unsigned int i = 0; i < n_rows; ++i)
    {
      if ((out[i] != out_batched[i]) || (out[i] != out_batched[n_rows + i]) ||
          (out[i] != out_batched[2 * n_rows + i]))
        deallog << "Error at entry " << i << std::endl;
      else
        deallog << "Correct line " << i << std::endl;
    }
  for (unsigned int i = 0; i < n_rows; ++i)
    {
      if ((out[i] != out_batched_2[i]) || (out[i] != out_batched_2[n_rows + i]))
        deallog << "Error at entry " << i << std::endl;
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
