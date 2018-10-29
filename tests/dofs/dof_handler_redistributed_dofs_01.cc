//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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

  SmartPointer<const FiniteElement<dim, spacedim>> fe_p(&dh.get_fe());

  dh.distribute_dofs(*fe_p);

  deallog << "OK" << std::endl;
}
