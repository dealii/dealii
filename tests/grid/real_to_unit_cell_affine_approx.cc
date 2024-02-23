// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check TriaAccessor::real_to_unit_cell_affine_approximation

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim, int spacedim>
void
do_test(const Triangulation<dim, spacedim> &tria, const Point<spacedim> &p)
{
  MappingQ<dim, spacedim> mapping(1);

  for (typename Triangulation<dim, spacedim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    {
      deallog << "Point p=(" << p << ") on cell with id " << cell->id() << ": "
              << cell->real_to_unit_cell_affine_approximation(p) << std::endl;
      Point<dim> mapping_point;
      try
        {
          mapping_point = mapping.transform_real_to_unit_cell(cell, p);
          mapping_point -= cell->real_to_unit_cell_affine_approximation(p);
          deallog << "Distance to mapping point: ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << mapping_point[d] << ' ';
          deallog << std::endl;
        }
      catch (const typename Mapping<dim, spacedim>::ExcTransformationFailed &)
        {
          deallog << "No MappingQ transform possible for this cell and point."
                  << std::endl;
        }
    }
}



template <int dim>
void
test1()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  Point<dim> p;
  for (unsigned int d = 0; d < dim; ++d)
    p[d] = -0.2 + 0.3 * d;

  do_test(tria, p);
}



template <int dim>
void
test2()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  Point<dim> p;
  for (unsigned int d = 0; d < dim; ++d)
    p[d] = -0.4 + 0.5 * d;

  do_test(tria, p);
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

  Point<spacedim> p;
  for (unsigned int d = 0; d < dim; ++d)
    p[d] = -0.4 + 0.5 * d;

  do_test(triangulation, p);
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
