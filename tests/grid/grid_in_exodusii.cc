// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2020 by the deal.II authors
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

// Read in several Exodus meshes and save the results to output files.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <string>

#include "../tests.h"

using namespace dealii;

template <int dim, int spacedim = dim>
void
read_and_print(const std::string &filename,
               const bool         read_as_manifolds = false)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  const typename GridIn<dim, spacedim>::ExodusIIData exodus_data =
    gi.read_exodusii(filename, read_as_manifolds);
  const auto &boundary_id_to_sideset_ids = exodus_data.id_to_sideset_ids;

  for (unsigned int i = 0; i < boundary_id_to_sideset_ids.size(); ++i)
    {
      deallog << "boundary id = " << i << " sideset ids = ";
      for (const int v : boundary_id_to_sideset_ids[i])
        deallog << v << ", ";
      deallog << std::endl;
    }

  for (const auto face : tria.active_face_iterators())
    {
      if (face->at_boundary())
        {
          deallog << "face center = " << face->center();
          if (read_as_manifolds)
            deallog << "   face manifold id = " << face->manifold_id();
          else
            deallog << "   face boundary id = " << face->boundary_id();
          deallog << std::endl;
        }
    }

  deallog << "Number of vertices: " << tria.get_vertices().size() << std::endl;
  deallog << "Number of cells: " << tria.n_cells() << std::endl;

  GridOut go;
  go.write_vtk(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(2);

  deallog << "-----------" << std::endl;
  deallog << "five-square" << std::endl;
  deallog << "-----------" << std::endl;
  read_and_print<2>(SOURCE_DIR "/grids/exodusii/five.e");

  deallog << "-----------------" << std::endl;
  deallog << "five-square in 3D" << std::endl;
  deallog << "-----------------" << std::endl;
  read_and_print<2, 3>(SOURCE_DIR "/grids/exodusii/five-3d.e");

  deallog << "-----------------------------" << std::endl;
  deallog << "Straightforward 3D mixed mesh" << std::endl;
  deallog << "-----------------------------" << std::endl;
  read_and_print<3>(SOURCE_DIR "/grids/exodusii/four.e");

  deallog << "--------------------------" << std::endl;
  deallog << "Higher-order 3D mixed mesh" << std::endl;
  deallog << "--------------------------" << std::endl;
  read_and_print<3>(SOURCE_DIR "/grids/exodusii/four-p.e");

  deallog << "--------------------------------------" << std::endl;
  deallog << "3D mixed mesh with boundary indicators" << std::endl;
  deallog << "--------------------------------------" << std::endl;
  read_and_print<3>(SOURCE_DIR "/grids/exodusii/four-sidesets.e");

  deallog << "--------------------------------------" << std::endl;
  deallog << "3D mixed mesh with manifold indicators" << std::endl;
  deallog << "--------------------------------------" << std::endl;
  read_and_print<3>(SOURCE_DIR "/grids/exodusii/four-sidesets.e", true);


  deallog << "OK" << std::endl;
}
