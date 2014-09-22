// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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



// check MatrixTools::local_apply_boundary_values with elimination of
// columns

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
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;

  // set up a DoFHandler and compute hanging
  // node constraints
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  // then set up a sparsity pattern and two
  // matrices on top of it. similar for two
  // vectors
  SparsityPattern sparsity (dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity);
  sparsity.compress ();
  SparseMatrix<double> A(sparsity), B(sparsity);
  Vector<double> b1 (dof_handler.n_dofs());
  Vector<double> b2 (dof_handler.n_dofs());

  // then fill the two matrices and vectors
  // by setting up bogus matrix entries and
  // (1) writing them into the matrix and
  // applying boundary values later on, or
  // (2) applying them right away
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ConstantFunction<dim>(1.),
                                            boundary_values);

  std::vector<types::global_dof_index> local_dofs (fe.dofs_per_cell);
  FullMatrix<double> local_matrix (fe.dofs_per_cell, fe.dofs_per_cell);
  Vector<double> local_vector (fe.dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      cell->get_dof_indices (local_dofs);
      local_matrix = 0;
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          local_matrix(i,j) = (i+1.)*(j+1.)*(local_dofs[i]+1.)*(local_dofs[j]+1.);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        local_vector(i) = (i+1.)*(local_dofs[i]+1.);

      // copy local to global by ourselves
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          A.add (local_dofs[i], local_dofs[j], local_matrix(i,j));
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        b1(local_dofs[i]) += local_vector(i);

      // or let other functions do that after
      // removing boundary values
      MatrixTools::local_apply_boundary_values (boundary_values, local_dofs,
                                                local_matrix, local_vector,
                                                true);
      cell->distribute_local_to_global (local_matrix, B);
      cell->distribute_local_to_global (local_vector, b2);
    }

  // for A, remove boundary values only now.
  Vector<double> x (dof_handler.n_dofs());
  MatrixTools::apply_boundary_values (boundary_values, A, x, b1, true);

  // now comes the check: we subtract B from
  // A, and make sure that the result is zero
  A.add (-1., B);
  deallog << "|A|=" << A.frobenius_norm() << std::endl;
  deallog << "|B|=" << B.frobenius_norm() << std::endl;
  Assert (A.frobenius_norm() < 1e-12*B.frobenius_norm(),
          ExcInternalError());

  // similar for b1 and b2
  b1 -= b2;
  deallog << "|b1|=" << b1.l2_norm() << std::endl;
  deallog << "|b2|=" << b2.l2_norm() << std::endl;
  Assert (b1.l2_norm() < 1e-12*b2.l2_norm(),
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
