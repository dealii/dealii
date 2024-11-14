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
  std::map<std::string, std::vector<double>> field_data;

  field_data = gridin.get_field_data();
  if (field_data.size() > 0)
    {
      std::map<std::string, std::vector<double>>::const_iterator iter =
        field_data.find("Density");
      if (iter != field_data.end())
        {
          std::vector<double> cell_density = iter->second;
          for (typename Triangulation<dim>::active_cell_iterator cell =
                 triangulation.begin_active();
               cell != triangulation.end();
               ++cell)
            {
              // store cell ID as a string
              std::string text = boost::lexical_cast<std::string>(cell->id());
              // Convert the string containing cell ID to an integer
              boost::char_separator<char>                   sep{"_"};
              boost::tokenizer<boost::char_separator<char>> tokens(text, sep);
              unsigned int                                  cell_id =
                boost::lexical_cast<int>(*(tokens.begin()));
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
