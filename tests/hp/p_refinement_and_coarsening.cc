// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// Check execution of p-refinement and p-coarsening via
// Triangulation::execute_coarsening_and_refinement()


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

  hp::FECollection<dim> fe;
  for (unsigned int i = 0; i < std::pow(2, dim); ++i)
    fe.push_back(FE_Q<dim>(1));

  DoFHandler<dim> dh(tria);

  // set future_fe_indices
  unsigned int future_feidx = 0;
  for (const auto &cell : dh.active_cell_iterators())
    {
      // check if cell is initialized correctly
      Assert(cell->active_fe_index() == 0, ExcInternalError());

      cell->set_future_fe_index(future_feidx);
      future_feidx = ((future_feidx + 1) < fe.size()) ? future_feidx + 1 : 0;
    }

  dh.distribute_dofs(fe);
  tria.execute_coarsening_and_refinement();

  // check if all flags were cleared and verify fe_indices
  for (const auto &cell : dh.active_cell_iterators())
    {
      Assert(cell->future_fe_index_set() == false, ExcInternalError());

      deallog << "cell:" << cell->id().to_string()
              << ", fe_index:" << cell->active_fe_index() << std::endl;
    }
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
