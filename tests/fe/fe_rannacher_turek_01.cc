// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Interfaces being tested
#include <deal.II/base/polynomials_rannacher_turek.h>

#include <deal.II/fe/fe_rannacher_turek.h>
// Interfaces needed for testing
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Known results for the Rannacher Turek element. All output should be zero.

void
test_known_values()
{
  PolynomialsRannacherTurek<2> pols;
  Point<2>                     p(0.5, 0.5);
  for (unsigned int i = 0; i < 4; ++i)
    {
      deallog << pols.compute_value(i, p) - 0.25 << std::endl;
    }
}

void
test_n_dofs()
{
  FE_RannacherTurek<2> fe_ratu;

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);

  DoFHandler<2> dofh(tria);
  dofh.distribute_dofs(fe_ratu);

  deallog << dofh.n_dofs() - 12 << std::endl;
}

void
test_nodal_matrix()
{
  FE_RannacherTurek<2> fe;

  FullMatrix<double> N(4, 4);
  // compute_node_matrix does not check for number of components but
  // simply assumes that there are dim components. Thus we do it ourselves.
  const unsigned int           n_dofs = fe.dofs_per_cell;
  const std::vector<Point<2>> &points = fe.get_generalized_support_points();
  std::vector<Vector<double>>  values(points.size(), Vector<double>(1));
  std::vector<double>          local_dofs(n_dofs);

  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      for (unsigned int k = 0; k < values.size(); ++k)
        values[k][0] = fe.shape_value(i, points[k]);
      fe.convert_generalized_support_point_values_to_dof_values(values,
                                                                local_dofs);

      for (unsigned int j = 0; j < n_dofs; ++j)
        N(j, i) = local_dofs[j];
    }

  for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        {
          if (i == j)
            {
              deallog << N(i, j) - 1.0;
            }
          else
            {
              deallog << N(i, j) - 0.0;
            }
          if (j + 1 < 4)
            deallog << ' ';
          else
            deallog << std::endl;
        }
    }
}

void
test_interpolation()
{
  Triangulation<2> tr;
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(2);

  FE_RannacherTurek<2> fe;
  const unsigned int   n_dofs = fe.dofs_per_cell;

  DoFHandler<2> dofh(tr);
  dofh.distribute_dofs(fe);

  Vector<double> input_vector(dofh.n_dofs());
  for (unsigned int i = 0; i < input_vector.size(); ++i)
    {
      input_vector[i] = double(i);
    }

  Quadrature<2> quadrature(fe.get_generalized_support_points());
  FEValues<2>   fev(fe, quadrature, update_values | update_JxW_values);

  using cell_it = DoFHandler<2>::cell_iterator;
  cell_it cell  = dofh.begin_active();
  for (; cell != dofh.end(); ++cell)
    {
      fev.reinit(cell);

      std::vector<Vector<double>> values(quadrature.size(), Vector<double>(1));
      fev.get_function_values(input_vector, values);

      std::vector<double> interpolated_local_dofs(n_dofs);
      fe.convert_generalized_support_point_values_to_dof_values(
        values, interpolated_local_dofs);

      Vector<double> local_dofs(n_dofs);
      cell->get_dof_values(input_vector, local_dofs);

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          deallog << local_dofs[j] - interpolated_local_dofs[j] << ' ';
        }
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  test_known_values();
  test_n_dofs();
  test_nodal_matrix();
  test_interpolation();

  return 0;
}
