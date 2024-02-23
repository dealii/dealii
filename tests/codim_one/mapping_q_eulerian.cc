// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// computes points in real space starting from some quadrature points on the
// unit element

#include "../tests.h"

// all include files you need here

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <string>

template <int dim, int spacedim>
void
test(std::string filename, unsigned int degree)
{
  Triangulation<dim, spacedim> tria;
  FE_Q<dim, spacedim>          base_fe(degree);
  FESystem<dim, spacedim>      fe(base_fe, spacedim);

  DoFHandler<dim, spacedim> shift_dh(tria);

  GridIn<dim, spacedim> gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);
  gi.read_ucd(in);

  shift_dh.distribute_dofs(fe);

  Vector<double> shift(shift_dh.n_dofs());

  shift.add(+1);

  GridOut grid_out;
  grid_out.set_flags(GridOutFlags::Ucd(true));
  grid_out.write_ucd(tria, deallog.get_file_stream());

  QTrapezoid<dim>                                 quad;
  MappingQEulerian<dim, Vector<double>, spacedim> mapping(degree,
                                                          shift_dh,
                                                          shift);

  typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();
  Point<spacedim> real;
  Point<dim>      unit;
  double          eps = 1e-10;
  for (; cell != endc; ++cell)
    {
      deallog << cell << std::endl;
      for (unsigned int q = 0; q < quad.size(); ++q)
        {
          real = mapping.transform_unit_to_real_cell(cell, quad.point(q));
          unit = mapping.transform_real_to_unit_cell(cell, real);
          deallog << quad.point(q) << " -> " << real << std::endl;
          if ((unit - quad.point(q)).norm() > eps)
            deallog << "Error: " << quad.point(q) << " != " << unit
                    << std::endl;
        }
    }
}

int
main()
{
  initlog();

  test<1, 2>(SOURCE_DIR "/grids/circle_1.inp", 1);
  test<1, 2>(SOURCE_DIR "/grids/circle_1.inp", 2);
  test<2, 3>(SOURCE_DIR "/grids/square.inp", 1);
  test<2, 3>(SOURCE_DIR "/grids/square.inp", 2);
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp", 1);
  // test<2,3>(SOURCE_DIR "/grids/problem.inp",1);
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp", 2);

  return 0;
}
