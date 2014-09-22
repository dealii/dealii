// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// check GeometryInfo::face_to_cell_vertices

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/geometry_info.h>

#include <fstream>
#include <cstdlib>


template <int dim>
void test ()
{
  deallog << "Checking in " << dim << "d" << std::endl;

  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
      {
        deallog << "Face " << f << ", vertex=" << v << ": ";
        deallog << GeometryInfo<dim>::face_to_cell_vertices(f,v,true)
                << std::endl;
      }

  if (dim == 3)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
        {
          deallog << "Face " << f << ", vertex=" << v
                  << " (reverse orientation): ";
          deallog << GeometryInfo<dim>::face_to_cell_vertices(f,v,false)
                  << std::endl;
        }
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
