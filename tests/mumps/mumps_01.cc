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


// test the mumps sparse direct solver on a mass matrix
// - check several ways to solve

#include "../tests.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>


void solve_and_check(const SparseMatrix<double> & M,
		     const Vector<double> & rhs,
		     const Vector<double> & solution)
{
  {
    SparseDirectMUMPS solver;
    solver.initialize (M);
    Vector<double> dst(rhs.size());
    solver.vmult(dst, rhs);
    dst -= solution;
    Assert(dst.l2_norm() < 1e-9, ExcInternalError());
  }
  {
    SparseDirectMUMPS solver;
    solver.initialize (M, rhs);
    Vector<double> dst(rhs.size());
    solver.solve(dst);
    dst -= solution;
    Assert(dst.l2_norm() < 1e-9, ExcInternalError());
  }
}



template <int dim>
void test ()
{
  deallog << dim << 'd' << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria,0,1);
  tria.refine_global (1);

  // destroy the uniformity of the matrix by
  // refining one cell
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global(8-2*dim);

  FE_Q<dim> fe (1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  deallog << "Number of dofs = " << dof_handler.n_dofs() << std::endl;

  SparsityPattern sparsity_pattern;
  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();

  SparseMatrix<double> B;
  B.reinit (sparsity_pattern);

  QGauss<dim> qr (2);
  MatrixTools::create_mass_matrix (dof_handler, qr, B);

  // compute a decomposition of the matrix
  SparseDirectMUMPS Binv;
  Binv.initialize (B);

  // for a number of different solution
  // vectors, make up a matching rhs vector
  // and check what the solver finds
  for (unsigned int i=0; i<3; ++i)
    {
      Vector<double> solution (dof_handler.n_dofs());
      Vector<double> x (dof_handler.n_dofs());
      Vector<double> b (dof_handler.n_dofs());

      for (unsigned int j=0; j<dof_handler.n_dofs(); ++j)
        solution(j) = j+j*(i+1)*(i+1);

      B.vmult (b, solution);

      Binv.vmult (x,b);

      x -= solution;
      deallog << "relative norm distance = "
              << x.l2_norm() / solution.l2_norm()
              << std::endl;
      deallog << "absolute norms = "
              << x.l2_norm() << ' ' << solution.l2_norm()
              << std::endl;
      Assert (x.l2_norm() / solution.l2_norm() < 1e-8,
              ExcInternalError());

      
				       //also check solving in different ways:
      solve_and_check(B, b, solution);
    }
}


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  initlog();
  
  test<1> ();
  test<2> ();
  test<3> ();
}
