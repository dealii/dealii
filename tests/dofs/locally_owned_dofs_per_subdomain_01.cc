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


// Test that the structures reporting DoF indexing are sized correctly for
// the hp::DoFHandler. In particular, we're interested in the case when one
// subdomain has no DoFs associated with it.

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <fstream>

#include "../tests.h"



int
main()
{
  initlog();

  const int dim = 2;

  const unsigned int n_subdomains = 2;

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Nothing<dim>());

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);
  GridTools::partition_triangulation(n_subdomains, tria);

  DoFHandler<dim> hp_dof_handler(tria);
  for (auto &cell : hp_dof_handler.active_cell_iterators())
    {
      if (cell == hp_dof_handler.begin_active())
        cell->set_active_fe_index(0); // FE_Q
      else
        cell->set_active_fe_index(1); // FE_Nothing
    }
  hp_dof_handler.distribute_dofs(fe_collection);

  Assert(DoFTools::locally_owned_dofs_per_subdomain(hp_dof_handler).size() ==
           n_subdomains,
         ExcDimensionMismatch(
           DoFTools::locally_owned_dofs_per_subdomain(hp_dof_handler).size(),
           n_subdomains));
  Assert(DoFTools::locally_relevant_dofs_per_subdomain(hp_dof_handler).size() ==
           n_subdomains,
         ExcDimensionMismatch(
           DoFTools::locally_relevant_dofs_per_subdomain(hp_dof_handler).size(),
           n_subdomains));

  deallog << "OK" << std::endl;
}
