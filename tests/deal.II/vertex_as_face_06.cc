// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// verify that we can do things like cell->face() in 1d as well. here:
// test cell->face(0)->get_dof_indices()


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>


template <int spacedim>
void test ()
{
  Triangulation<1,spacedim> tria;
  GridGenerator::hyper_cube (tria);

  FESystem<1,spacedim> fe(FE_Q<1,spacedim>(2),1,
                          FE_Q<1,spacedim>(1),1);
  DoFHandler<1,spacedim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  std::vector<types::global_dof_index> dof_indices(fe.dofs_per_face);

  deallog << "Coarse mesh:" << std::endl;
  dof_handler.begin_active()->face(0)->get_dof_indices (dof_indices);
  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
    deallog << "Left vertex=" << dof_indices[i] << std::endl;
  dof_handler.begin_active()->face(1)->get_dof_indices (dof_indices);
  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
    deallog << "Right vertex=" << dof_indices[i] << std::endl;

  tria.refine_global (2);
  dof_handler.distribute_dofs (fe);

  for (typename DoFHandler<1,spacedim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      cell->face(0)->get_dof_indices (dof_indices);
      for (unsigned int i=0; i<fe.dofs_per_face; ++i)
        deallog << "Left vertex=" << dof_indices[i] << std::endl;
      cell->face(1)->get_dof_indices (dof_indices);
      for (unsigned int i=0; i<fe.dofs_per_face; ++i)
        deallog << "Right vertex=" << dof_indices[i] << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
