// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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
#include <deal.II/base/logstream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <unistd.h>


template <int dim, int spacedim>
void write_vtk (const parallel::distributed::Triangulation<dim,spacedim> &tria,
                const char                                *filename)
{
  deallog << "Checksum: "
          << tria.get_checksum ()
          << std::endl;

  tria.write_mesh_vtk (filename);

  // copy the .pvtu and .vtu files
  // into the logstream
  {
    unsigned int c=0;
    std::ifstream in((std::string(filename) + ".pvtu").c_str());
    while (in && ++c<20)
      {
        std::string s;
        std::getline(in, s);
        deallog.get_file_stream() << s << "\n";
      }
  }

  {
    unsigned int c=0;
    std::ifstream in((std::string(filename) + "_0000.vtu").c_str());
    while (in && ++c<20)
      {
        std::string s;
        std::getline(in, s);
        deallog.get_file_stream() << s << "\n";
      }
  }
}


template <int dim>
void assert_tria_equal(const Triangulation<dim> &a, const Triangulation<dim> &b)
{
  Assert (a.n_active_cells() == b.n_active_cells(),
          ExcInternalError());

  std::string file1 = "tmp_grid1";
  std::string file2 = "tmp_grid2";

  std::ofstream out1(file1.c_str());
  GridOut().write_gnuplot (a, out1);
  std::ofstream out2(file2.c_str());
  GridOut().write_gnuplot (b, out2);

  out1.close();
  out2.close();

  //compare the two files
  std::string cmd = std::string("diff -q ")+file1+std::string(" ")+file2;
  Assert (system(cmd.c_str()) == 0, ExcInternalError());

  //and delete them
  std::remove (file1.c_str());
  std::remove (file2.c_str());
}
