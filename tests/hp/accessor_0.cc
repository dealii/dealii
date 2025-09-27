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

// Ensure that the newly added instantiations for DoFAccessor<0, ...>
// (dof_index and get_fe) work correctly for an hp::DoFHandler.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::boolalpha;


  Triangulation<1> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  hp::FECollection<1> fe_collection;
  DoFHandler<1>       dof_handler(triangulation);
  fe_collection.push_back(FE_Q<1>(2));
  fe_collection.push_back(FE_Q<1>(4));
  fe_collection.push_back(FE_Q<1>(6));

  const unsigned int n_fe_indices = 3;
  {
    typename DoFHandler<1>::active_cell_iterator cell =
      dof_handler.begin_active();
    dof_handler.begin_active()->set_active_fe_index(1);
    ++cell; // go to cell 1
    ++cell; // go to cell 2
    ++cell; // go to cell 3
    cell->set_active_fe_index(2);
  }
  dof_handler.distribute_dofs(fe_collection);

  // dof_index
  // get_fe

  std::vector<types::global_dof_index> dof_indices;

  typename DoFHandler<1>::active_cell_iterator cell =
                                                 dof_handler.begin_active(),
                                               endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      deallog << "===================================" << std::endl;
      deallog << "cell center: " << cell->center() << std::endl;

      dof_indices.resize(fe_collection[cell->active_fe_index()].dofs_per_cell);
      cell->get_dof_indices(dof_indices);
      deallog << "cell dofs: ";
      for (unsigned int dof_n = 0; dof_n < dof_indices.size(); ++dof_n)
        {
          deallog << dof_indices[dof_n];
          if (dof_n != dof_indices.size() - 1)
            {
              deallog << ", ";
            }
        }
      deallog << std::endl;

      // see if we have a neighbor on the right. If so, the common vertex
      // should be associated with two FE indices.
      const typename DoFHandler<1>::active_cell_iterator neighbor =
        cell->neighbor(1);
      if (neighbor != dof_handler.end())
        {
          const unsigned int current_index  = cell->active_fe_index();
          const unsigned int neighbor_index = neighbor->active_fe_index();
          deallog << "dof index (current cell, current index): "
                  << cell->face(1)->dof_index(0, current_index)
                  << " (neighbor cell, current index): "
                  << neighbor->face(0)->dof_index(0, current_index) << std::endl
                  << "dof index (current cell, neighbor index): "
                  << cell->face(1)->dof_index(0, neighbor_index)
                  << " (neighbor cell, neighbor index): "
                  << neighbor->face(0)->dof_index(0, neighbor_index)
                  << std::endl;
        }

      for (unsigned int fe_index = 0; fe_index < n_fe_indices; ++fe_index)
        {
          const bool index_is_active = cell->fe_index_is_active(fe_index);
          deallog << "cell uses fe index " << fe_index << ": "
                  << index_is_active << std::endl;

          for (const unsigned int face_n : GeometryInfo<1>::face_indices())
            {
              AssertThrow(&cell->face(face_n)->get_fe(fe_index) ==
                            &fe_collection[fe_index],
                          ExcMessage("The result of get_fe should always return"
                                     " a known finite element."));

              if (index_is_active)
                {
                  deallog << "vertex dof index: "
                          << cell->face(face_n)->dof_index(0, fe_index)
                          << std::endl;
                }
            }
        }
    }

  deallog << "OK" << std::endl;

  return 0;
}
