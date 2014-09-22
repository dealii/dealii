// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



// check FullMatrix::swap_col and FullMatrix::swap_row for nonsymmetric
// matrices; there used to be a bug in the implementation


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "output";


template <typename number>
void
check ()
{
  FullMatrix<number> m (5, 7);
  fill_matrix (m);
  print_matrix (m);

  m.swap_col (2, 4);
  print_matrix (m);

  m.swap_row (2, 4);
  print_matrix (m);
}

