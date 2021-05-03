// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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



// trigger an error in hp::DoFHandler::create_active_fe_table


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_DGQ<dim>(1));

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe_collection);

  tria.refine_global(1);
  dof_handler.distribute_dofs(fe_collection);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe_collection);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe_collection);
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
