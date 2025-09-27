// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test GridTools::minimal_cell_diameter and GridTools::maximal_cell_diameter
// with a mapping

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 2


template <int dim, int spacedim>
void
test(const unsigned int degree)
{
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;
  deallog << "degree = " << degree << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim, spacedim>     fe(degree);
  FESystem<dim, spacedim> fe_sys(fe, spacedim);

  DoFHandler<dim, spacedim> dof_sys(tria);
  dof_sys.distribute_dofs(fe_sys);

  Vector<double> euler;
  euler.reinit(dof_sys.n_dofs());
  const ComponentMask mask(spacedim, true);

  VectorTools::get_position_vector(dof_sys, euler, mask);

  // By using MappingQEulerian with the position vector, we are making
  // everything bigger by a factor 2.
  MappingQEulerian<dim, Vector<double>, spacedim> map_fe(degree,
                                                         dof_sys,
                                                         euler);

  deallog << "Min diameter        : " << GridTools::minimal_cell_diameter(tria)
          << std::endl
          << "Max diameter        : " << GridTools::maximal_cell_diameter(tria)
          << std::endl
          << "Min mapped diameter : "
          << GridTools::minimal_cell_diameter(tria, map_fe) << std::endl
          << "Max mapped diameter : "
          << GridTools::maximal_cell_diameter(tria, map_fe) << std::endl;
}

int
main()
{
  initlog();

  for (unsigned int d = 1; d < 4; ++d)
    {
      test<2, 2>(d);
      test<2, 3>(d);
      test<3, 3>(d);
    }
}
