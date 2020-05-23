// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// test airfoil generator

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
print_triangulation(const Triangulation<dim, dim> &tria)
{
  for (const auto &cell : tria)
    {
      deallog << cell.material_id() << " ";

      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        deallog << (cell.face(f)->at_boundary() ? cell.face(f)->boundary_id() :
                                                  -1)
                << " ";

      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        deallog << cell.vertex(v) << " ";

      deallog << std::endl;
    }
}

template <int dim>
void
check(const std::string type)
{
  deallog.push(type);

  // without periodic boundaries
  deallog.push("regular");
  {
    Triangulation<dim, dim> tria;

    GridGenerator::Airfoil::AdditionalData additional_data;
    additional_data.airfoil_type = type;
    GridGenerator::Airfoil::create_triangulation(tria, additional_data);

    print_triangulation(tria);
  }
  deallog.pop();

  // with periodic boundaries
  deallog.push("periodic");
  {
    std::vector<GridTools::PeriodicFacePair<
      typename Triangulation<dim, dim>::cell_iterator>>
      periodic_faces;

    Triangulation<dim, dim> tria;

    GridGenerator::Airfoil::AdditionalData additional_data;
    additional_data.airfoil_type = type;
    GridGenerator::Airfoil::create_triangulation(tria,
                                                 periodic_faces,
                                                 additional_data);

    print_triangulation(tria);
  }
  deallog.pop();

  deallog.pop();
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<2>("NACA");
  check<2>("Joukowski");
}
