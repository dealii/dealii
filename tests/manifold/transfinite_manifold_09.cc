// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This test verifies that the transfinite interpolation works on a distorted
// cylinder with some given coordinates

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(9);

  const double       radius = 1.;
  const double       length = 3.;
  constexpr int      dim    = 3;
  Triangulation<dim> tria;

  // position of auxiliary point to achieve an angle of 120 degrees in corner
  // of inner cell
  std::vector<Point<3>> vertices{
    {1.6666666666666663, 0., -0.69035593728849154},
    {2.4033740890369071, 0., 0.046351485081748842},
    {2., 0., -0.82842712474618985},
    {2.7222718241315023, 0., -0.10615530061468759},
    {1.5548344933395739, 0.44817247098067337, -0.64403353438675115},
    {2.2998293798138509, 0.44396471873216636, 0.073847495332288168},
    {1.9264041086678665, 0.53755670408577649, -0.79794270842148385},
    {2.6272400227990982, 0.52932727817817915, -0.10032762350281543}};

  std::vector<CellData<3>> cell_data(1);
  for (unsigned int v = 0; v < 8; ++v)
    cell_data[0].vertices[v] = v;

  SubCellData subcell_data;
  tria.create_triangulation(vertices, cell_data, subcell_data);

  // set cylindrical manifold on face 3
  tria.set_all_manifold_ids(1);
  tria.begin_active()->face(3)->set_all_manifold_ids(2);

  // attach 3 cylindrical manifolds to mesh
  tria.set_manifold(
    2,
    CylindricalManifold<dim>(Point<dim>{std::sqrt(0.5), 0., std::sqrt(0.5)},
                             Point<dim>()));
  TransfiniteInterpolationManifold<dim> transfinite;
  transfinite.initialize(tria);
  tria.set_manifold(1, transfinite);

  deallog << "Cell center " << tria.begin_active()->center(true, true)
          << std::endl;

  tria.refine_global();
  for (auto &cell : tria.active_cell_iterators())
    {
      deallog << "Face centers:" << std::endl;
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        deallog << cell->face(f)->center(true) << std::endl;
      deallog << std::endl;
    }
}
