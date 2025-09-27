// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// CellDataTransfer test by refinement
// for sequential Triangulation


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/cell_data_transfer.templates.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  // ------ setup ------
  Triangulation<dim, spacedim> tria;

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  // ----- prepare -----
  // set refinement/coarsening flags
  for (auto &cell : tria.active_cell_iterators())
    {
      if (cell->id().to_string() == "0_1:0")
        cell->set_refine_flag();
      else if (((dim == 1) && (cell->parent()->id().to_string() == "1_0:")) ||
               ((dim == 2) && (cell->parent()->id().to_string() == "3_0:")) ||
               ((dim == 3) && (cell->parent()->id().to_string() == "7_0:")))
        cell->set_coarsen_flag();
    }

  // ----- gather -----
  // store parent id of all cells
  std::vector<unsigned int> cell_ids_pre(tria.n_active_cells());
  for (auto &cell : tria.active_cell_iterators())
    {
      const std::string  parent_cellid = cell->parent()->id().to_string();
      const unsigned int parent_coarse_cell_id =
        static_cast<unsigned int>(std::stoul(parent_cellid));
      cell_ids_pre[cell->active_cell_index()] = parent_coarse_cell_id;

      deallog << "cellid=" << cell->id()
              << " parentid=" << cell_ids_pre[cell->active_cell_index()];
      if (cell->coarsen_flag_set())
        deallog << " coarsening";
      else if (cell->refine_flag_set())
        deallog << " refining";
      deallog << std::endl;
    }

  // ----- transfer -----
  CellDataTransfer<dim, spacedim, std::vector<unsigned int>> cell_data_transfer(
    tria);

  cell_data_transfer.prepare_for_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  std::vector<unsigned int> cell_ids_post(tria.n_active_cells());
  cell_data_transfer.unpack(cell_ids_pre, cell_ids_post);

  // ------ verify ------
  // check if all children adopted the correct id
  for (auto &cell : tria.active_cell_iterators())
    deallog << "cellid=" << cell->id()
            << " parentid=" << cell_ids_post[(cell->active_cell_index())]
            << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  initlog();

  deallog.push("1d");
  test<1, 1>();
  deallog.pop();
  deallog.push("2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d");
  test<3, 3>();
  deallog.pop();
}
