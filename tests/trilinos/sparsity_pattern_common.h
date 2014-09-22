// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// check

#include "../tests.h"
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include "../lac/testmatrix.h"

const unsigned int N = 15;


// reinitialize sparsity patterns for 5-point star
void do_reinit (TrilinosWrappers::SparsityPattern &sp)
{
  sp.reinit((N-1)*(N-1), (N-1)*(N-1));
}



void build_sparsity (TrilinosWrappers::SparsityPattern &sparsity_pattern)
{
  // generate usual 5-point sparsity pattern
  do_reinit (sparsity_pattern);
  FDMatrix(N,N).five_point_structure (sparsity_pattern);
  sparsity_pattern.compress ();

  deallog << sparsity_pattern.n_rows() << " "
          << sparsity_pattern.n_cols() << " "
          << sparsity_pattern.bandwidth() << " "
          << sparsity_pattern.n_nonzero_elements()
          << std::endl;
}



void row_length ()
{
  TrilinosWrappers::SparsityPattern sparsity_pattern;
  build_sparsity (sparsity_pattern);

  for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
    deallog << sparsity_pattern.row_length(i) << std::endl;

  deallog << "OK" << std::endl;
}



void print_gnuplot ()
{
  TrilinosWrappers::SparsityPattern sparsity_pattern;
  build_sparsity (sparsity_pattern);

  sparsity_pattern.print_gnuplot(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



void print ()
{
  TrilinosWrappers::SparsityPattern sparsity_pattern;
  build_sparsity (sparsity_pattern);

  sparsity_pattern.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



