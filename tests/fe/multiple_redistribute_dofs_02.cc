// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test multiple redistribute dofs.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/fe_field_function.h>

#include <iostream>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  deallog << "Testing for dim = " << dim << ", spacedim = " << spacedim
          << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  hp::FECollection<dim, spacedim> fe_collection(FE_Q<dim, spacedim>(1));

  DoFHandler<dim, spacedim> dh(tria);
  dh.distribute_dofs(fe_collection);

  // This stores a pointer to the FE in dh.
  hp::FEValues<dim, spacedim> fe_v(dh.get_fe_collection(),
                                   hp::QCollection<dim>(QGauss<dim>(1)),
                                   update_values);
  fe_v.reinit(dh.begin_active());

  tria.refine_global(1);
  dh.distribute_dofs(fe_collection);

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
