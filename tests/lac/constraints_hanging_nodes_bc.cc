// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this function tests the correctness of the implementation of mixing hanging
// node constraints and Dirichlet boundary constraints. This is trivial in 2D
// because then hanging nodes never appear on the boundary, but needs to be
// checked in 3D. To do this, compare a constraint matrix that we fill
// manually with one that we get from using the library functions
// make_hanging_node_constraints and interpolate_boundary_values.

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  // refine some of the cells.
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (unsigned int counter = 0; cell != endc; ++cell, ++counter)
    if (counter % 5 == 0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  cell = tria.begin_active();
  endc = tria.end();
  for (unsigned int counter = 0; cell != endc; ++cell, ++counter)
    if (counter % 8 == 0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> correct_constraints, library_constraints;

  DoFTools::make_hanging_node_constraints(dof, correct_constraints);
  library_constraints.merge(correct_constraints);

  {
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(
      dof, 0, Functions::ConstantFunction<dim>(1.), boundary_values);
    std::map<types::global_dof_index, double>::const_iterator boundary_value =
      boundary_values.begin();
    for (; boundary_value != boundary_values.end(); ++boundary_value)
      {
        if (!correct_constraints.is_constrained(boundary_value->first))
          {
            correct_constraints.constrain_dof_to_zero(boundary_value->first);
            correct_constraints.set_inhomogeneity(boundary_value->first,
                                                  boundary_value->second);
          }
      }
  }
  correct_constraints.close();

  deallog << "Number of DoFs: " << dof.n_dofs() << std::endl
          << "Number of hanging nodes: " << library_constraints.n_constraints()
          << std::endl
          << "Total number of constraints: "
          << correct_constraints.n_constraints() << std::endl;

  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ConstantFunction<dim>(1.),
                                           library_constraints);
  library_constraints.close();

  // the two constraint matrices should look the
  // same, so go through them and check
  deallog << "Check that both constraint matrices are identical... ";
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      AssertThrow(correct_constraints.is_constrained(i) ==
                    library_constraints.is_constrained(i),
                  ExcInternalError());
      using constraint_format =
        const std::vector<std::pair<types::global_dof_index, double>> &;
      if (correct_constraints.is_constrained(i))
        {
          constraint_format correct =
            *correct_constraints.get_constraint_entries(i);
          constraint_format library =
            *library_constraints.get_constraint_entries(i);
          AssertThrow(correct.size() == library.size(), ExcInternalError());
          for (unsigned int q = 0; q < correct.size(); ++q)
            {
              AssertThrow(correct[q].first == library[q].first,
                          ExcInternalError());
              AssertThrow(std::fabs(correct[q].second - library[q].second) <
                            1e-14,
                          ExcInternalError());
            }
          AssertThrow(std::fabs(correct_constraints.get_inhomogeneity(i) -
                                library_constraints.get_inhomogeneity(i)) <
                        1e-14,
                      ExcInternalError());
        }
    }
  deallog << "OK." << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  {
    deallog.push("3d");
    test<3>();
    deallog.pop();
  }
}
