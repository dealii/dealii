// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
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


// read a 1d 3d file in the VTK format with point data

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"

template <int dim, int spacedim>
void
check_file(const std::string                      name,
           typename GridIn<dim, spacedim>::Format format)
{
  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim>        gridin;
  GridOut                      gridout;

  gridin.attach_triangulation(triangulation);
  gridin.read(name, format);
  deallog << '\t' << triangulation.n_vertices() << '\t'
          << triangulation.n_cells() << std::endl;
  std::map<std::string, Vector<double>> point_data;

  point_data = gridin.get_point_data();
  if (point_data.size() > 0)
    {
      std::map<std::string, Vector<double>>::const_iterator iter_scalar =
        point_data.find("path_distance");
      if (iter_scalar != point_data.end())
        {
          Vector<double> cell_path_distance = iter_scalar->second;
          for (typename Triangulation<dim, spacedim>::active_cell_iterator
                 cell = triangulation.begin_active();
               cell != triangulation.end();
               ++cell)
            {
              unsigned int cell_id = cell->id().get_coarse_cell_id();
              deallog << '\t' << cell_path_distance[cell_id] << std::endl;
            }
        }



      deallog << "Velocity field at vertices" << std::endl;
      std::map<std::string, Vector<double>>::const_iterator iter_vel =
        point_data.find("velocity");
      if (iter_vel != point_data.end())
        {
          Vector<double> vec = iter_vel->second;
          for (typename Triangulation<dim, spacedim>::active_cell_iterator
                 cell = triangulation.begin_active();
               cell != triangulation.end();
               ++cell)
            {
              for (unsigned int v = 0; v < cell->n_vertices(); ++v)
                {
                  // Print data at vertex index
                  for (unsigned int i = 0; i < spacedim; ++i)
                    deallog << vec[cell->vertex_index(v) * spacedim + i] << " ";
                  deallog << std::endl;
                }
            }
        }
    }
}

void
filename_resolution()
{
  check_file<1, 3>(std::string(SOURCE_DIR "/grid_in_vtk_1d_3d_point_data/mesh"),
                   GridIn<1, 3>::vtk);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  filename_resolution();
}
