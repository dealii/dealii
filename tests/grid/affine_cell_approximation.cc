// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check GridTools::affine_cell_approximation

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
do_test(const Triangulation<dim, spacedim> &tria)
{
  for (const auto &cell : tria.cell_iterators())
    {
      std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
        vertices;
      for (const unsigned int v : cell->vertex_indices())
        vertices[v] = cell->vertex(v);
      deallog << "Cell " << cell->id() << " with center " << cell->center()
              << ": A = [";
      const auto A =
        GridTools::affine_cell_approximation<dim, spacedim>(vertices).first;
      for (unsigned int d = 0; d < spacedim; ++d)
        {
          for (unsigned int e = 0; e < dim; ++e)
            deallog << A[d][e] << (e < dim - 1 ? " " : "");
          deallog << (d < spacedim - 1 ? ", " : "");
        }
      deallog
        << "] b = "
        << GridTools::affine_cell_approximation<dim, spacedim>(vertices).second
        << std::endl;
    }
}



template <int dim>
void
test1()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  do_test(tria);
}



template <int dim>
void
test2()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  do_test(tria);
}



template <int dim, int spacedim>
void
test3()
{
  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim>        grid_in;
  grid_in.attach_triangulation(triangulation);
  if (dim == 1)
    {
      std::ifstream fname(SOURCE_DIR "/../codim_one/grids/circle_1.inp");
      grid_in.read_ucd(fname);
    }
  else
    {
      std::ifstream fname(SOURCE_DIR "/../codim_one/grids/sphere_0.inp");
      grid_in.read_ucd(fname);
    }

  do_test(triangulation);
}


int
main()
{
  initlog();

  test1<1>();
  test1<2>();
  test1<3>();

  test2<2>();
  test2<3>();

  test3<1, 2>();
  test3<2, 3>();
}
