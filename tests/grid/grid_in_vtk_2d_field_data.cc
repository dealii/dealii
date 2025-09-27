// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2024 by the deal.II authors
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


// read a 2d file in the VTK format with field data

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <string>

#include "../tests.h"

template <int dim>
void
check_file(const std::string name, typename GridIn<dim>::Format format)
{
  Triangulation<dim> triangulation;
  GridIn<dim>        gridin;
  GridOut            gridout;

  gridin.attach_triangulation(triangulation);
  gridin.read(name, format);
  deallog << '\t' << triangulation.n_vertices() << '\t'
          << triangulation.n_cells() << std::endl;
  std::map<std::string, Vector<double>> cell_data;

  cell_data = gridin.get_cell_data();
  if (cell_data.size() > 0)
    {
      std::map<std::string, Vector<double>>::const_iterator iter =
        cell_data.find("Density");
      if (iter != cell_data.end())
        {
          Vector<double> cell_density = iter->second;
          for (typename Triangulation<dim>::active_cell_iterator cell =
                 triangulation.begin_active();
               cell != triangulation.end();
               ++cell)
            {
              unsigned int cell_id = cell->id().get_coarse_cell_id();
              deallog << '\t' << cell_density[cell_id] << std::endl;
            }
        }
    }

  gridout.write_gnuplot(triangulation, deallog.get_file_stream());
}

void
filename_resolution()
{
  check_file<2>(std::string(SOURCE_DIR "/grid_in_vtk_2d_field_data/mesh"),
                GridIn<2>::vtk);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  filename_resolution();
}
