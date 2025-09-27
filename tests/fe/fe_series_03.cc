// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test FESeries::process_coefficients()

#include <deal.II/fe/fe_series.h>

#include <iostream>

#include "../tests.h"


std::pair<bool, unsigned int>
pred_ind(const TableIndices<2> &ind)
{
  return std::make_pair(true, ind[0] + ind[1]);
}

void
test2d(const VectorTools::NormType norm)
{
  const unsigned int dim = 2;
  const unsigned int N   = 4;
  Table<dim, double> coefficients(4, 4);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      coefficients(i, j) = i * N + j;

  std::pair<std::vector<unsigned int>, std::vector<double>> res =
    FESeries::process_coefficients<2, double>(coefficients, pred_ind, norm);

  for (unsigned int i = 0; i < res.first.size(); ++i)
    deallog << res.first[i] << " : " << res.second[i] << std::endl;
}



int
main()
{
  initlog();

  deallog << "L2_norm" << std::endl;
  test2d(VectorTools::L2_norm);
  deallog << "L1_norm" << std::endl;
  test2d(VectorTools::L1_norm);
  deallog << "Linfty_norm" << std::endl;
  test2d(VectorTools::Linfty_norm);
  deallog << "mean" << std::endl;
  test2d(VectorTools::mean);
}
