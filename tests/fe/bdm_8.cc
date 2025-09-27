// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// build a mass matrix for the BDM element and try to invert it. we had trouble
// with this at one time

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 5


template <int dim>
void
test(const unsigned int degree)
{
  FE_BDM<dim>        fe_rt(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  DoFHandler<dim>                         dof(tr);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(fe_rt);

  QTrapezoid<1>      q_trapez;
  const unsigned int div = 4;
  QIterated<dim>     q(q_trapez, div);
  FEValues<dim>      fe(fe_rt, q, update_values | update_JxW_values);
  fe.reinit(c);

  const unsigned int dofs_per_cell = fe_rt.dofs_per_cell;
  FullMatrix<double> mass_matrix(dofs_per_cell, dofs_per_cell);

  Assert(fe.get_fe().n_components() == dim, ExcInternalError());


  for (unsigned int q_point = 0; q_point < q.size(); ++q_point)
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        for (unsigned int d = 0; d < dim; ++d)
          mass_matrix(i, j) +=
            (fe.shape_value_component(i, q_point, d) *
             fe.shape_value_component(j, q_point, d) * fe.JxW(q_point));

  mass_matrix.print_formatted(deallog.get_file_stream(), 3, false, 0, "0", 1);

  SolverControl           solver_control(2 * dofs_per_cell, 1e-8);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              solver(solver_control, vector_memory);

  Vector<double> tmp1(dofs_per_cell), tmp2(dofs_per_cell);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    tmp1(i) = random_value<double>();

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
  initlog();
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;

  for (unsigned int i = 1; i < 4; ++i)
    {
      test<2>(i);
      if (i < 3)
        test<3>(i);
    }

  return 0;
}
