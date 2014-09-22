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


// test by Brian: check some of the scaling operations on matrices

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>



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

  QGauss<2> qr (2);
  FEValues<2> fe_values (dof_handler.get_fe(),
                         qr, UpdateFlags(update_values    |
                                         update_gradients |
                                         update_q_points  |
                                         update_JxW_values));

  MatrixTools::create_laplace_matrix (dof_handler, qr, B.block(0,0));
  B.block(1,1).add (0,0,1.);

  B.print_formatted (deallog.get_file_stream (),3,false);

  B.block(0,0) *= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);

  B *= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);

  B /= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);

  B.block(0,0) /= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);

  B = 0;

  return 0;
}
