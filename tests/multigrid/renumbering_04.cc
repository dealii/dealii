// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
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


// Test DoFRenumbering::block_wise(dh, level)

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <algorithm>

#include "../tests.h"

using namespace std;

template <int dim>
void
check(FiniteElement<dim> &fe)
{
  deallog << std::endl << "**** " << fe.get_name() << std::endl;

  Triangulation<dim> tria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);
  dh.distribute_mg_dofs();

  deallog << "** before:" << std::endl;
  {
    unsigned int                         n_dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dofs(n_dofs_per_cell);
    for (unsigned int level = 0; level < tria.n_levels(); ++level)
      {
        deallog << "* level " << level << std::endl;
        for (typename DoFHandler<dim>::cell_iterator cell = dh.begin(level);
             cell != dh.end(level);
             ++cell)
          {
            deallog << "cell " << cell->id() << ":" << std::endl;
            cell->get_mg_dof_indices(local_dofs);
            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              deallog << local_dofs[i] << " ";
            deallog << std::endl;
          }
      }
  }

  for (unsigned int level = 0; level < tria.n_levels(); ++level)
    DoFRenumbering::block_wise(dh, level);

  deallog << std::endl << "** after:" << std::endl;
  {
    unsigned int                         n_dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dofs(n_dofs_per_cell);
    for (unsigned int level = 0; level < tria.n_levels(); ++level)
      {
        deallog << "* level " << level << std::endl;
        for (typename DoFHandler<dim>::cell_iterator cell = dh.begin(level);
             cell != dh.end(level);
             ++cell)
          {
            deallog << "cell " << cell->id() << ":" << std::endl;
            cell->get_mg_dof_indices(local_dofs);
            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              deallog << local_dofs[i] << " ";
            deallog << std::endl;
          }
      }
  }
}



int
main()
{
  initlog(__FILE__);

  {
    FESystem<1> fe(FE_Q<1>(1), 1, FE_DGQ<1>(0), 1);
    check<1>(fe);
  }


  {
    FESystem<2> fe(FE_RaviartThomas<2>(0), 1, FE_DGQ<2>(0), 1);
    check<2>(fe);
  }
  {
    FESystem<2> fe(FE_Q<2>(1), 1, FE_DGQ<2>(0), 1);
    check<2>(fe);
  }

  {
    FESystem<3> fe(FE_Q<3>(1), 1, FE_DGQ<3>(0), 1);
    check<3>(fe);
  }

  deallog << "OK" << endl;
}
