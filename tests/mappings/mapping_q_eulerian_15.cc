// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests that a MappingQEulerian mapping returns
// identical values to a MappingQ if the shift vector
// is zero.

#include "../tests.h"

// all include files you need here

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <string>

template <int dim, int spacedim>
void
test(unsigned int degree)
{
  deallog << "Dim: " << dim << ". Spacedim: " << spacedim
          << ". Degree: " << degree << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_shell(tria, Point<spacedim>(), 0.5, 1);
  FESystem<dim, spacedim> fe(FE_Q<dim, spacedim>(degree), dim);

  DoFHandler<dim, spacedim> shift_dh(tria);

  shift_dh.distribute_dofs(fe);

  Vector<double> shift(shift_dh.n_dofs());

  QGauss<dim>                                     quad(degree + 1);
  MappingQEulerian<dim, Vector<double>, spacedim> mapping_q_eulerian(degree,
                                                                     shift_dh,
                                                                     shift);

  MappingQ<dim, spacedim> mapping_q(degree);

  typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();

  UpdateFlags update_flags = UpdateFlags(update_quadrature_points);

  FEValues<dim> fe_values_q(mapping_q, fe, quad, update_flags);
  FEValues<dim> fe_values_q_eulerian(mapping_q_eulerian,
                                     fe,
                                     quad,
                                     update_flags);

  for (; cell != endc; ++cell)
    {
      fe_values_q.reinit(cell);
      fe_values_q_eulerian.reinit(cell);

      deallog << cell << std::endl;
      for (unsigned int q = 0; q < quad.size(); ++q)
        {
          deallog << "MappingQ: " << fe_values_q.quadrature_point(q)
                  << ". MappingQEulerian: "
                  << fe_values_q_eulerian.quadrature_point(q)
                  << ". Normalized difference: "
                  << (fe_values_q.quadrature_point(q) -
                      fe_values_q_eulerian.quadrature_point(q)) /
                       fe_values_q.quadrature_point(q).norm()
                  << std::endl;
        }
    }
}

int
main()
{
  initlog();

  test<2, 2>(1);
  test<2, 2>(2);
  test<2, 2>(3);
}
