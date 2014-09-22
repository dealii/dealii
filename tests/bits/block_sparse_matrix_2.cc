// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


// BlockSparseMatrix::clear used to forget to reset all sizes to zero

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>



int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria,0,1);
  tria.refine_global (1);

  FE_Q_Hierarchical<2> fe (1);
  DoFHandler<2> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.reinit (2,2);
  sparsity_pattern.collect_sizes ();

  sparsity_pattern.block(0,0).reinit (dof_handler.n_dofs(),
                                      dof_handler.n_dofs(),
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(0,1).reinit (dof_handler.n_dofs(),
                                      1,
                                      1);
  sparsity_pattern.block(1,0).reinit (1,
                                      dof_handler.n_dofs(),
                                      dof_handler.n_dofs());
  sparsity_pattern.block(1,1).reinit (1,
                                      1,
                                      1);
  sparsity_pattern.collect_sizes ();

  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern.block(0,0));

  for (unsigned int j=0; j<dof_handler.n_dofs (); ++j)
    {
      sparsity_pattern.block(0,1).add (j,0);
      sparsity_pattern.block(1,0).add (0,j);
    }
  sparsity_pattern.block(1,1).add (0,0);
  sparsity_pattern.compress ();

  BlockSparseMatrix<double> B;
  B.reinit (sparsity_pattern);

  // check some sizes
  Assert (B.m() == B.block(0,0).m() + B.block(1,1).m(),
          ExcInternalError());
  Assert (B.n() == B.block(0,0).n() + B.block(1,1).n(),
          ExcInternalError());
  Assert (B.n_block_rows() == 2, ExcInternalError());
  Assert (B.n_block_cols() == 2, ExcInternalError());

  // then clear, and check again
  B.clear ();
  Assert (B.m() == 0, ExcInternalError());
  Assert (B.n() == 0, ExcInternalError());
  Assert (B.n_block_rows() == 0, ExcInternalError());
  Assert (B.n_block_cols() == 0, ExcInternalError());

  deallog << "OK" << std::endl;

  return 0;
}
