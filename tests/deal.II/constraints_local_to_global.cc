// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



// this function tests the correctness of the implementation of
// ConstraintMatrix::distribute_local_to_global for FullMatrix by comparing
// the results with a sparse matrix. As a test case, we use a square mesh that
// is refined once globally and then the first cell is refined adaptively.

#include "../tests.h"

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <fstream>
#include <iostream>
#include <complex>

std::ofstream logfile("output");

template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.begin()->face(0)->set_boundary_indicator(1);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim> fe (1);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  VectorTools::interpolate_boundary_values (dof, 1, ZeroFunction<dim>(),
                                            constraints);
  constraints.close();

  SparsityPattern sparsity;
  {
    CompressedSimpleSparsityPattern csp (dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern (dof, csp, constraints, false);
    sparsity.copy_from (csp);
  }
  SparseMatrix<double> sparse (sparsity);
  FullMatrix<double> full (dof.n_dofs(), dof.n_dofs());

  FullMatrix<double> local_mat (fe.dofs_per_cell, fe.dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);

  // loop over cells, fill local matrix with
  // random values, insert both into sparse and
  // full matrix. Make some random entries equal
  // to zero
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof.begin_active(), endc = dof.end();
  unsigned int counter = 0;
  for ( ; cell != endc; ++cell)
    {
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j, ++counter)
          if (counter % 42 == 0)
            local_mat(i,j) = 0;
          else
            local_mat (i,j) = (double)Testing::rand() / RAND_MAX;
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (local_mat, local_dof_indices,
                                              sparse);
      constraints.distribute_local_to_global (local_mat, local_dof_indices,
                                              full);
    }

  // now check that the entries are indeed the
  // same by copying the sparse matrix into a
  // full matrix and checking the Frobenius norm
  // of the difference matrix
  FullMatrix<double> ref;
  ref.copy_from (sparse);
  full.add (-1., ref);
  deallog << "Difference between full and sparse matrix: "
          << full.frobenius_norm() << std::endl;
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-14);

  test<2>();
}

