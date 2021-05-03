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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// same as fe_series_03, but with ignoring some coefficients by predicate.

#include <deal.II/fe/fe_series.h>

#include <iostream>

#include "../tests.h"


std::pair<bool, unsigned int>
pred_ind(const TableIndices<2> &ind)
{
  const unsigned int val = ind[0] + ind[1];
  if (val < 4)
    return std::make_pair(true, val);
  else
    return std::make_pair(false, val);
}

void
test2d(const VectorTools::NormType norm)
{
  const unsigned int dim = 2;
  const unsigned int N   = 4;
  Table<dim, double> coefficients(4, 4);
  for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
      coefficients(i, j) = i * N + j;

  std::pair<std::vector<unsigned int>, std::vector<double>> res =
    FESeries::process_coefficients<2, double>(coefficients, pred_ind, norm);

  for (unsigned int i = 0; i < res.first.size(); i++)
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
