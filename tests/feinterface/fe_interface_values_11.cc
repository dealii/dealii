// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test basic properties of FEInterfaceValues in the hp case.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
inspect_fiv(FEInterfaceValues<dim> &fiv)
{
  deallog << "at_boundary(): " << fiv.at_boundary() << "\n"
          << "n_current_interface_dofs(): " << fiv.n_current_interface_dofs()
          << "\n";

  std::vector<types::global_dof_index> indices =
    fiv.get_interface_dof_indices();
  Assert(indices.size() == fiv.n_current_interface_dofs(), ExcInternalError());

  deallog << "interface_dof_indices: ";
  for (auto i : indices)
    deallog << i << ' ';
  deallog << "\n";


  unsigned int idx = 0;
  for (auto v : indices)
    {
      deallog << "  index " << idx << " global_dof_index:" << v << ":\n";

      const auto pair = fiv.interface_dof_to_dof_indices(idx);
      deallog << "    dof indices: " << static_cast<int>(pair[0]) << " | "
              << static_cast<int>(pair[1]) << "\n";

      ++idx;
    }

  deallog << std::endl;
}



template <int dim>
void
test(const unsigned int p)
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim>          dofh(tria);
  hp::FECollection<dim>    fe_collection;
  hp::QCollection<dim - 1> q_collection;

  fe_collection.push_back(FE_DGQ<dim>(p));
  fe_collection.push_back(FE_DGQ<dim>(p + 1));

  q_collection.push_back(QGauss<dim - 1>(p));
  q_collection.push_back(QGauss<dim - 1>(p + 1));

  // Set different finite elements spaces on the two cells.
  unsigned int fe_index = 0;
  for (const auto &cell : dofh.active_cell_iterators())
    {
      cell->set_active_fe_index(fe_index);
      ++fe_index;
    }

  dofh.distribute_dofs(fe_collection);

  UpdateFlags update_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;

  FEInterfaceValues<dim> fiv(fe_collection, q_collection, update_flags);


  auto cell = dofh.begin();

  deallog << "** interface between cell 0 and 1 **\n";

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    if (!cell->at_boundary(f))
      {
        fiv.reinit(cell,
                   f,
                   numbers::invalid_unsigned_int,
                   cell->neighbor(f),
                   cell->neighbor_of_neighbor(f),
                   numbers::invalid_unsigned_int);

        Assert(fiv.get_fe_face_values(0).get_cell() == cell,
               ExcInternalError());
        Assert(fiv.get_fe_face_values(1).get_cell() == cell->neighbor(f),
               ExcInternalError());
        Assert(fiv.n_current_interface_dofs() ==
                 fe_collection[0].n_dofs_per_cell() +
                   fe_collection[1].n_dofs_per_cell(),
               ExcInternalError());
        Assert(!fiv.at_boundary(), ExcInternalError());

        auto mycell = cell;
        for (unsigned int c = 0; c < 2; ++c)
          {
            std::vector<types::global_dof_index> indices(
              fe_collection[c].n_dofs_per_cell());
            mycell->get_dof_indices(indices);
            deallog << "cell " << c << ": ";
            for (auto i : indices)
              deallog << i << ' ';
            deallog << "\n";
            ++mycell;
          }

        inspect_fiv(fiv);
      }

  deallog << "** boundary interface on cell 0 **\n" << std::endl;

  {
    fiv.reinit(cell, 0);
    Assert(fiv.get_fe_face_values(0).get_cell() == cell, ExcInternalError());
    Assert(fiv.n_current_interface_dofs() == fe_collection[0].n_dofs_per_cell(),
           ExcInternalError());
    Assert(fiv.at_boundary(), ExcInternalError());
    inspect_fiv(fiv);
  }
}



int
main()
{
  initlog();
  for (const unsigned int p : {1, 2})
    test<2>(p);
  for (const unsigned int p : {1})
    test<3>(p);
}
