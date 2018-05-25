// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// Create a mass matrix for the RT_Bubbles element and try to invert it.
// Same as the rt_bubbles_5 test except we use a library function to
// create the mass matrix

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_rt_bubbles.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/matrix_tools.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 5


std::ofstream logfile("output");

template <int dim>
void
test(const unsigned int degree)
{
  FE_RT_Bubbles<dim> fe_rt_bubbles(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe_rt_bubbles);

  QTrapez<1>         q_trapez;
  const unsigned int div = 4;
  QIterated<dim>     q(q_trapez, div);

  const unsigned int dofs_per_cell = fe_rt_bubbles.dofs_per_cell;
  SparsityPattern    sp(dofs_per_cell, dofs_per_cell, dofs_per_cell);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
      sp.add(i, j);
  sp.compress();
  SparseMatrix<double> mass_matrix(sp);

  MatrixTools::create_mass_matrix(dof, q, mass_matrix);

  mass_matrix.print_formatted(logfile, 3, false, 0, "0", 1);

  SolverControl           solver_control(3 * dofs_per_cell, 1e-8);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              solver(solver_control, vector_memory);

  Vector<double> tmp1(dofs_per_cell), tmp2(dofs_per_cell);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    tmp1(i) = 1. * Testing::rand() / RAND_MAX;

  deallog << "solving degree = " << degree << std::endl;
  check_solver_within_range(
    solver.solve(mass_matrix, tmp2, tmp1, PreconditionIdentity()),
    solver_control.last_step(),
    3,
    45);
}



int
main()
{
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  for (unsigned int i = 1; i < 4; ++i)
    {
      test<2>(i);
      test<3>(i);
    }

  return 0;
}
