// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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

// computes points in real space starting from some quadrature points on the unit element

#include "../tests.h"

// all include files you need here

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <string>

std::ofstream logfile("output");

template <int dim, int spacedim>
void
test(std::string filename)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename.c_str());
  gi.read_ucd(in);

  GridOut grid_out;
  grid_out.set_flags(GridOutFlags::Ucd(true));
  grid_out.write_ucd(tria, logfile);

  QTrapez<dim>                                                quad;
  MappingQGeneric<dim, spacedim>                              mapping(1);
  typename Triangulation<dim, spacedim>::active_cell_iterator cell
    = tria.begin_active(),
    endc = tria.end();
  Point<spacedim> real;
  Point<dim>      unit;
  double          eps = 1e-10;
  for(; cell != endc; ++cell)
    {
      deallog << cell << std::endl;
      for(unsigned int q = 0; q < quad.size(); ++q)
        {
          real = mapping.transform_unit_to_real_cell(cell, quad.point(q));
          unit = mapping.transform_real_to_unit_cell(cell, real);
          deallog << quad.point(q) << " -> " << real << std::endl;
          if((unit - quad.point(q)).norm() > eps)
            deallog << "Error: " << quad.point(q) << " != " << unit
                    << std::endl;
        }
    }
}

int
main()
{
  deallog.attach(logfile);

  test<1, 2>(SOURCE_DIR "/grids/circle_1.inp");
  test<2, 3>(SOURCE_DIR "/grids/square.inp");
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}
