// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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


// Create a Triangulation, interpolate its boundary points to a close
// smooth BSpline, and use that as a Boundary Descriptor.

#include "../tests.h"

#include <deal.II/opencascade/utilities.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>

using namespace OpenCASCADE;

void
remove_iges_header(const std::string &in_filename,
                   const std::string &out_filename)
{
  std::ifstream in(in_filename);
  std::ofstream out(out_filename);
  std::string line;
  unsigned int counter = 5;
  while (counter--) std::getline(in,line);
  while (std::getline(in,line))
    out << line << std::endl;
  in.close();
  out.close();
}

template<int spacedim>
void
test()
{
  deallog << "Testing <2," << spacedim << ">" << std::endl;

  Triangulation<2,spacedim> tria1;
  Triangulation<2,spacedim> tria2;
  Triangulation<2,spacedim> tria3;

  GridGenerator::hyper_cube(tria1, 0, 1);
  GridGenerator::hyper_cube(tria2, 2, 3);

  GridGenerator::merge_triangulations(tria1, tria2, tria3);

  tria3.refine_global(1);

  auto v1 = create_curves_from_triangulation_boundary(tria1);
  auto v3 = create_curves_from_triangulation_boundary(tria3);

  if (spacedim == 3)
    {
      Triangulation<2,3> tria4;
      GridGenerator::hyper_sphere(tria4);
      auto v4 = create_curves_from_triangulation_boundary(tria4);

      if (v4.size() != 0)
        deallog << "Not OK!: size is " << v4.size() << std::endl;
      else
        deallog << "OK" << std::endl;
    }

  if (v1.size() != 1)
    deallog << "Not OK: size is " << v1.size() << std::endl;
  else
    deallog << "OK" << std::endl;

  if (v3.size() != 2)
    deallog << "Not OK!: size is " << v3.size() << std::endl;
  else
    deallog << "OK" << std::endl;

  write_IGES(v1[0], "edge1.iges");
  write_IGES(v3[0], "edge30.iges");
  write_IGES(v3[1], "edge31.iges");

  remove_iges_header("edge1.iges", "edge1_noheader.iges");
  remove_iges_header("edge30.iges", "edge30_noheader.iges");
  remove_iges_header("edge31.iges", "edge31_noheader.iges");

  cat_file("edge1_noheader.iges");
  cat_file("edge30_noheader.iges");
  cat_file("edge31_noheader.iges");
}

int
main ()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
