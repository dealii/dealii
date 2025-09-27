// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that we can do things like cell->face() in 1d as well. here:
// test cell->face(0)->get_dof_indices()
// compared to _06, we now test for an hp-DoFHandler

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int spacedim>
void
print_dofs(const typename DoFHandler<1, spacedim>::face_iterator &i,
           const unsigned int                                     fe_index,
           const unsigned int                                     n)
{
  std::vector<types::global_dof_index> dof_indices(n);
  i->get_dof_indices(dof_indices, fe_index);
  for (unsigned int i = 0; i < n; ++i)
    deallog << dof_indices[i] << ' ';
  deallog << std::endl;
}



template <int spacedim>
void
print_dofs(const typename DoFHandler<1, spacedim>::cell_iterator &i,
           const unsigned int                                     n)
{
  std::vector<types::global_dof_index> dof_indices(n);
  i->get_dof_indices(dof_indices);
  for (unsigned int i = 0; i < n; ++i)
    deallog << dof_indices[i] << ' ';
  deallog << std::endl;
}



template <int spacedim>
void
test()
{
  Triangulation<1, spacedim> tria;
  GridGenerator::hyper_cube(tria);

  FESystem<1, spacedim> fe1(FE_Q<1, spacedim>(2), 1, FE_Q<1, spacedim>(1), 1);
  FESystem<1, spacedim> fe2(FE_Q<1, spacedim>(3), 1, FE_Q<1, spacedim>(2), 1);
  hp::FECollection<1, spacedim> fe_collection;
  fe_collection.push_back(fe1);
  fe_collection.push_back(fe2);

  DoFHandler<1, spacedim> dof_handler(tria);
  dof_handler.begin_active()->set_active_fe_index(0);
  dof_handler.distribute_dofs(fe_collection);

  deallog << "Coarse mesh:" << std::endl;
  print_dofs<spacedim>(dof_handler.begin_active()->face(0),
                       0,
                       fe1.dofs_per_face);
  print_dofs<spacedim>(dof_handler.begin_active()->face(1),
                       0,
                       fe1.dofs_per_face);

  tria.refine_global(2);
  {
    unsigned int index = 0;
    for (typename DoFHandler<1, spacedim>::active_cell_iterator cell =
           dof_handler.begin_active();
         cell != dof_handler.end();
         ++cell, index = (index + 1) % fe_collection.size())
      cell->set_active_fe_index(index);
  }
  dof_handler.distribute_dofs(fe_collection);

  for (typename DoFHandler<1, spacedim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      deallog << "Cell: " << cell
              << ", active_fe_index=" << cell->active_fe_index() << std::endl;

      print_dofs<spacedim>(cell, cell->get_fe().dofs_per_cell);
      print_dofs<spacedim>(cell->face(0),
                           cell->active_fe_index(),
                           cell->get_fe().dofs_per_face);
      print_dofs<spacedim>(cell->face(1),
                           cell->active_fe_index(),
                           cell->get_fe().dofs_per_face);
    }
}



int
main()
{
  initlog();

  test<1>();
  test<2>();

  return 0;
}
