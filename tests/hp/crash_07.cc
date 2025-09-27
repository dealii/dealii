// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a distilled version of hp_constraints_q_05 and a couple of other tests
// where DoFs were constrained to other DoFs already constrained

char logname[] = "output";


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"



int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  std::vector<Point<2>> points_glob;
  std::vector<Point<2>> points;

  points_glob.push_back(Point<2>(0.0, 0.0));
  points_glob.push_back(Point<2>(1.0, 0.0));
  points_glob.push_back(Point<2>(1.0, 0.5));
  points_glob.push_back(Point<2>(1.0, 1.0));
  points_glob.push_back(Point<2>(0.6, 0.5));
  points_glob.push_back(Point<2>(0.5, 1.0));
  points_glob.push_back(Point<2>(0.0, 1.0));

  // Prepare cell data
  std::vector<CellData<2>> cells(3);

  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 4;
  cells[0].vertices[3] = 2;
  cells[0].material_id = 0;

  cells[1].vertices[0] = 4;
  cells[1].vertices[1] = 2;
  cells[1].vertices[2] = 5;
  cells[1].vertices[3] = 3;
  cells[1].material_id = 0;

  cells[2].vertices[0] = 0;
  cells[2].vertices[1] = 4;
  cells[2].vertices[2] = 6;
  cells[2].vertices[3] = 5;
  cells[2].material_id = 0;

  Triangulation<2> triangulation;
  triangulation.create_triangulation(points_glob, cells, SubCellData());

  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  hp::FECollection<2> fe;
  fe.push_back(FE_Q<2>(1));
  fe.push_back(FE_Q<2>(2));

  DoFHandler<2> dof_handler(triangulation);

  // distribute fe_indices randomly
  unsigned int cell_no = 0;
  for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell, ++cell_no)
    cell->set_active_fe_index(0);
  (std::next(dof_handler.begin_active()))->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe);

  deallog << "n_dofs=" << dof_handler.n_dofs() << std::endl;

  for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      deallog << "Cell=" << cell << std::endl;
      deallog << "    vertices=" << cell->vertex_index(0) << ' '
              << cell->vertex_index(1) << ' ' << cell->vertex_index(2) << ' '
              << cell->vertex_index(3) << std::endl;
      deallog << "    active_fe_index=" << cell->active_fe_index() << std::endl;

      deallog << "    dofs=";
      std::vector<types::global_dof_index> local_dofs(
        fe[cell->active_fe_index()].dofs_per_cell);
      cell->get_dof_indices(local_dofs);
      for (unsigned int i = 0; i < fe[cell->active_fe_index()].dofs_per_cell;
           ++i)
        deallog << local_dofs[i] << ' ';
      deallog << std::endl;
    }


  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());
}
