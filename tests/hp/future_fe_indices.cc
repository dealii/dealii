// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2019 by the deal.II authors
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



// check if future FE indices will be set correctly


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Q<dim>(1));

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe_collection);

  // check if indices are initialized correctly
  for (const auto &cell : dh.active_cell_iterators())
    Assert(cell->future_fe_index_set() == false, ExcInternalError());

  // set future index at least once
  auto cell = dh.begin_active();
  cell->set_future_fe_index(1);

  // verify flags
  for (const auto &cell : dh.active_cell_iterators())
    {
      deallog << "cell:" << cell->id().to_string()
              << ", future_fe:" << cell->future_fe_index() << std::endl;
      Assert(&(dh.get_fe(cell->future_fe_index())) == &(cell->get_future_fe()),
             ExcMessage(
               "DoFCellAccessor::get_future_fe() returns the wrong object."));
    }

  // clear all flags and check if all were cleared
  for (const auto &cell : dh.active_cell_iterators())
    {
      cell->clear_future_fe_index();
      Assert(cell->future_fe_index_set() == false, ExcInternalError());
    }
}


int
main()
{
  initlog();

  deallog.push("1D");
  test<1>();
  deallog.pop();

  deallog.push("2D");
  test<2>();
  deallog.pop();

  deallog.push("3D");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
