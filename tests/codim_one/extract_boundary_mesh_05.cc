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



#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

using namespace std;


void test ()
{
  const unsigned int dim=2;

  Triangulation<dim-1,dim> boundary_mesh;
  map<Triangulation<dim-1,dim>::cell_iterator,Triangulation<dim,dim>::face_iterator >
  surface_to_volume_mapping;
  Triangulation<dim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);

  surface_to_volume_mapping
    = GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh);

  FE_Q <dim-1,dim>  boundary_fe (1);
  DoFHandler<dim-1,dim> boundary_dh(boundary_mesh);
  boundary_dh.distribute_dofs (boundary_fe);

  deallog << "n_dofs=" << boundary_dh.n_dofs() << std::endl;

  for (DoFHandler<dim-1,dim>::active_cell_iterator
       cell = boundary_dh.begin_active(),
       endc = boundary_dh.end(); cell!=endc; ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
        {
          unsigned int index = cell->vertex_dof_index(v,0);
          deallog << "vertex: " << v
                  << ", global: " << cell->vertex_index(v)
                  << " index: " << index << std::endl;
        }

      deallog << std::endl;
    }
}



int main ()
{

  ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test ();
}
