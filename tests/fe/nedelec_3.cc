// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 2



template <int dim>
inline void
check(const unsigned int p)
{
  deallog << dim << "-D" << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  // work on a once-refined grid
  tr.refine_global(1);

  FE_Nedelec<dim> fe_ned(p);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe_ned);

  // generate a function on the
  // coarse grid (with shape function
  // values equal to the index of the
  // line)
  Vector<double> coarse_grid_values(fe_ned.dofs_per_cell);
  for (unsigned int i = 0; i < coarse_grid_values.size(); ++i)
    coarse_grid_values(i) = i;

  // then transfer this function
  // from the coarse grid to the
  // child cells.
  Vector<double> values(dof.n_dofs());
  for (typename DoFHandler<dim>::active_cell_iterator c = dof.begin_active();
       c != dof.end();
       ++c)
    {
      Vector<double> tmp(fe_ned.dofs_per_cell);
      fe_ned.get_prolongation_matrix(c->index()).vmult(tmp, coarse_grid_values);
      c->set_dof_values(tmp, values);
    };


  // then output these values at the
  // quadrature points of all cells
  // of the finer grid
  QTrapezoid<dim>             quadrature;
  std::vector<Vector<double>> shape_values(quadrature.size(),
                                           Vector<double>(dim));
  FEValues<dim>               fe(fe_ned,
                   quadrature,
                   update_values | update_quadrature_points);

  for (typename DoFHandler<dim>::active_cell_iterator c = dof.begin_active();
       c != dof.end();
       ++c)
    {
      deallog << "  CELL " << c << std::endl;
      fe.reinit(c);
      fe.get_function_values(values, shape_values);

      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          deallog << ", xq=" << fe.quadrature_point(q) << ", f=[";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << (d == 0 ? "" : " ") << shape_values[q](d);

          deallog << ']' << std::endl;
        };

      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog << "Degree 0: " << std::endl;
  check<2>(0);
  check<3>(0);
  deallog << "Degree 1: " << std::endl;
  check<2>(1);
  check<3>(1);
  return 0;
}
