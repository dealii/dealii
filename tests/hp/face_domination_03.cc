// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Make sure both sides of a face know which of the two elements dominates


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
class MixedFECollection
{
public:
  MixedFECollection();
  ~MixedFECollection();

  void
  run();

  Triangulation<dim> triangulation;

  FESystem<dim> continuous_fe;
  FESystem<dim> discontinuous_fe;

  DoFHandler<dim>       dof_handler;
  hp::FECollection<dim> fe_collection;
};

template <int dim>
MixedFECollection<dim>::MixedFECollection()
  : continuous_fe(FE_Q<dim>(1), dim)
  , discontinuous_fe(FE_DGQ<dim>(1), dim)
  , dof_handler(triangulation)
{}

template <int dim>
MixedFECollection<dim>::~MixedFECollection()
{
  dof_handler.clear();
}

template <int dim>
void
MixedFECollection<dim>::run()
{
  // add two a CG and a DG finite element object to fe_collection
  fe_collection.push_back(continuous_fe);
  fe_collection.push_back(discontinuous_fe);

  // produce a simple grid with 16 cells
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  // looping over all cells and assigning the FE_DG object to the first cell
  // that comes up -- works. looping over all cells and assigning the FE_DG
  // object to the interior cells -- doesn't work.
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      //          if (cell == dof_handler.begin_active()) // this works
      if (cell->center().norm() < 0.5) // this doesn't
        {
          cell->set_active_fe_index(cell->active_fe_index() + 1);
        }
      deallog << cell->center() << "\t" << cell->active_fe_index() << std::endl;
    }

  dof_handler.distribute_dofs(fe_collection);

  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells() << std::endl
          << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl
          << "   Number of vertices: " << triangulation.n_vertices()
          << std::endl;
}



int
main()
{
  initlog();

  MixedFECollection<2> mixed_FECollection_problem;
  mixed_FECollection_problem.run();
}
