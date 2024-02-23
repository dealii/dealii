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

// It checks whether the cell_data_storage, specifically the methods
// initialize() and get_data(), works with two different materials
// types storing two different data structures.
// It makes sure the results are not affected during local refinement and
// and coarsening.

#include <deal.II/base/quadrature_point_data.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"

// history structs for each material
struct MaterialBase
{
  virtual ~MaterialBase()
  {}
};
struct Mat0 : MaterialBase
{
  double x{5.};
};
struct Mat1 : MaterialBase
{
  double y{8.};
};
constexpr unsigned int n_data_points_per_cell = 4;


void
setup_history(Triangulation<2>              &tria,
              CellDataStorage<typename Triangulation<2>::active_cell_iterator,
                              MaterialBase> &storage)
{
  deallog << "initializing history" << std::endl;
  for (auto cell : tria.active_cell_iterators())
    {
      deallog << "initializing cell " << cell->id() << std::endl;
      if (cell->material_id() == 0)
        storage.template initialize<Mat0>(cell, n_data_points_per_cell);
      else if (cell->material_id() == 1)
        storage.template initialize<Mat1>(cell, n_data_points_per_cell);
    }
}

void
read_history(Triangulation<2>              &tria,
             CellDataStorage<typename Triangulation<2>::active_cell_iterator,
                             MaterialBase> &storage)
{
  deallog << "reading history" << std::endl;
  for (auto cell : tria.active_cell_iterators())
    {
      if (cell->material_id() == 0)
        {
          auto data_vec = storage.template get_data<Mat0>(cell);
          deallog << "Cell with material id = 0 contains the data "
                  << data_vec[0]->x << std::endl;
        }
      else if (cell->material_id() == 1)
        {
          auto data_vec = storage.template get_data<Mat1>(cell);
          deallog << "Cell with material id = 1 contains the data "
                  << data_vec[0]->y << std::endl;
        }
    }
}

int
main()
{
  initlog();

  // create a mesh with two cells
  Triangulation<2> tria;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            /*repetitions*/ {2, 1},
                                            /*point1*/ {0., 0.},
                                            /*point2*/ {2., 1.});

  // assign material ids
  auto cell0 = tria.begin_active();
  cell0->set_material_id(0);
  auto cell1 = cell0;
  ++cell1;
  cell1->set_material_id(1);

  // refine one cell
  cell0->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  deallog << "first refinement done" << std::endl;

  // create a history structure and populate it
  CellDataStorage<typename Triangulation<2>::active_cell_iterator, MaterialBase>
    storage;
  setup_history(tria, storage);
  read_history(tria, storage);

  // coarsen the cell and refine the other cell
  for (auto cell : tria.cell_iterators_on_level(1))
    cell->set_coarsen_flag();
  cell1->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  deallog << "second refinement done" << std::endl;

  // initialize history for newly created cells
  setup_history(tria, storage);
  read_history(tria, storage);

  deallog << "OK" << std::endl;

  return 0;
}
