// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that filtered iterators can be used with using FEInterfaceValues.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);

  const FE_Q<dim, spacedim> fe(1);
  const QGauss<dim - 1>     quadrature(1);

  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  FEInterfaceValues<dim, spacedim> fe_interface_values(fe,
                                                       quadrature,
                                                       update_default);
  const unsigned int invalid_subface = numbers::invalid_unsigned_int;
  const unsigned int face_index      = 3;

  for (const auto &cell :
       filter_iterators(dof_handler.active_cell_iterators(),
                        IteratorFilters::ActiveFEIndexEqualTo(0)))
    {
      if (cell == dof_handler.begin_active())
        {
          const auto cell_neighbor = cell->neighbor(face_index);

          // This commented out code block identifies that the filtered iterator
          // is of a different type to that returned by cell->neighbor()
          //
          // using C1 = std::decay_t<decltype(cell)>;
          // using C2 = std::decay_t<decltype(cell_neighbor)>;
          // static_assert(std::is_same_v<C1, C2>,
          //               "Cell iterator types are not identical.");

          fe_interface_values.reinit(cell,
                                     face_index,
                                     invalid_subface,
                                     cell->neighbor(face_index),
                                     cell->neighbor_of_neighbor(face_index),
                                     invalid_subface);
        }
    }

  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::ActiveFEIndexEqualTo(0))
    {
      if (cell == dof_handler.begin_active())
        {
          fe_interface_values.reinit(cell,
                                     face_index,
                                     invalid_subface,
                                     cell->neighbor(face_index),
                                     cell->neighbor_of_neighbor(face_index),
                                     invalid_subface);
        }
    }
}


int
main()
{
  initlog();

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  deallog << "OK" << std::endl;
}
