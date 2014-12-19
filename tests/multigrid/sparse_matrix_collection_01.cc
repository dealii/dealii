// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/multigrid/sparse_matrix_collection.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;



int main()
{
  initlog();

  mg::SparseMatrixCollection<float> smc;
  smc.resize(0, 5);

  deallog << "matrix      " << smc.matrix.min_level() << '-' << smc.matrix.max_level() << std::endl;
  deallog << "matrix_in   " << smc.matrix_in.min_level() << '-' << smc.matrix_in.max_level() << std::endl;
  deallog << "matrix_out  " << smc.matrix_out.min_level() << '-' << smc.matrix_out.max_level() << std::endl;
  deallog << "matrix_up   " << smc.matrix_up.min_level() << '-' << smc.matrix_up.max_level() << std::endl;
  deallog << "matrix_down " << smc.matrix_down.min_level() << '-' << smc.matrix_down.max_level() << std::endl;
}
