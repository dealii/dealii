//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that calling hp::DoFhandler::distribute_dofs again with the same
// hp::FECollection doesn't recreate the copy.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog();

  constexpr int dim      = 2;
  constexpr int spacedim = 2;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  FE_Q<dim, spacedim> fe(1);

  DoFHandler<dim, spacedim> dh(tria);
  dh.distribute_dofs(fe);

  ObserverPointer<const FiniteElement<dim, spacedim>> fe_p(&dh.get_fe());

  dh.distribute_dofs(*fe_p);

  deallog << "OK" << std::endl;
}
