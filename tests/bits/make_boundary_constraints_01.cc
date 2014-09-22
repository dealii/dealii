// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// check DoFTools::make_zero_boundary_constraints by comparing
// apply_boundary_values with a zero function to
// make_zero_boundary_constraints

#include "../tests.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

  // refine the mesh in a random way
  triangulation.refine_global (4-dim);
  for (unsigned int i=0; i<11-2*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
      for (unsigned int index=0; cell != triangulation.end(); ++cell, ++index)
        if (index % (3*dim) == 0)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement ();
    }
  deallog << "Number of cells: "
          << triangulation.n_active_cells() << std::endl;

  // assign quadrature, set up a
  // DoFHandler, and distribute dofs
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs : "
          << dof_handler.n_dofs() << std::endl;

  // then set up a sparsity pattern
  // and two matrices and vectors on
  // top of it.
  SparsityPattern sparsity (dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity);
  sparsity.compress ();
  SparseMatrix<double> A(sparsity), B(sparsity);
  Vector<double> a1 (dof_handler.n_dofs());

  // initialize object denoting zero
  // boundary values and boundary
  // constraints
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ConstantFunction<dim>(1.),
                                            boundary_values);
  ConstraintMatrix constraints;
  constraints.clear();
  DoFTools::make_zero_boundary_constraints (dof_handler,
                                            constraints);
  constraints.close();

  // then fill two matrices by
  // setting up bogus matrix entries
  // and for A applying constraints
  // right away and for B applying
  // constraints later on
  std::vector<types::global_dof_index> local_dofs (fe.dofs_per_cell);
  FullMatrix<double> local_matrix (fe.dofs_per_cell, fe.dofs_per_cell);
  Vector<double> local_vector (fe.dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      cell->get_dof_indices (local_dofs);
      local_matrix = 0;

      // create local matrices
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
            local_matrix(i,j) = (i+1.)*(j+1.)*(local_dofs[i]+1.)*(local_dofs[j]+1.);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        local_vector(i) = (i+1.)*(local_dofs[i]+1.);

      // for matrix A and vector a1 apply boundary
      // values
      MatrixTools::local_apply_boundary_values (boundary_values, local_dofs,
                                                local_matrix, local_vector,
                                                true);
      cell->distribute_local_to_global (local_matrix, A);
      cell->distribute_local_to_global (local_vector, a1);

      // for matrix B cast constraint
      // matrix
      constraints.distribute_local_to_global (local_matrix,
                                              local_dofs,
                                              B);
    }

  // here comes the check: compare
  // the l1_norm of matrices A and B,
  // their difference should be zero
  deallog << "|A| = " << A.l1_norm() << std::endl;
  deallog << "|B| = " << B.l1_norm() << std::endl;
  Assert (A.l1_norm() - B.l1_norm() == 0,
          ExcInternalError());
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test<1> ();
      test<2> ();
      test<3> ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
