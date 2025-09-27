// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a version of hp_hanging_node_02 that crashed at the time
// of writing the time

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim>        triangulation;
  hp::FECollection<dim>     fe;
  DoFHandler<dim>           dof_handler(triangulation);
  AffineConstraints<double> hanging_node_constraints;

  FE_Q<dim> fe_1(1), fe_2(2), fe_3(QIterated<1>(QTrapezoid<1>(), 3)),
    fe_4(QIterated<1>(QTrapezoid<1>(), 4));

  fe.push_back(fe_1);
  fe.push_back(fe_2);
  fe.push_back(fe_3);
  fe.push_back(fe_4);



  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);
  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: " << triangulation.n_cells() << std::endl;

  // Now to the p-Method. Assign
  // random active_fe_indices to the
  // different cells.
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    cell->set_active_fe_index(Testing::rand() % fe.size());

  dof_handler.distribute_dofs(fe);

  cell = dof_handler.begin_active();
  for (; cell != endc; ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        deallog << "cell=" << cell << ", face=" << f
                << ", fe_index=" << cell->active_fe_index() << ", dofs=";
        std::vector<types::global_dof_index> dofs(
          fe[cell->active_fe_index()].dofs_per_face);
        cell->face(f)->get_dof_indices(dofs, cell->active_fe_index());
        for (unsigned int i = 0; i < dofs.size(); ++i)
          deallog << dofs[i] << ' ';
        deallog << std::endl;
      }

  DoFTools::make_hanging_node_constraints(dof_handler,
                                          hanging_node_constraints);

  hanging_node_constraints.print(deallog.get_file_stream());

  hanging_node_constraints.close();
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  // this test depends on the right
  // starting value of the random
  // number generator. set it to the
  // same value it has in
  // hp_hanging_nodes_02
  for (unsigned int i = 0; i < 64; ++i)
    Testing::rand();
  test<3>();

  return 0;
}
