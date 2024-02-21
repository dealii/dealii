/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2020 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Like convert_hypercube_to_simplex_mesh_01, but also refines the grid once.
 */

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim, int spacedim>
void
create_triangulation(Triangulation<dim, spacedim> &triangulation)
{
  GridGenerator::subdivided_hyper_cube(triangulation, 4);
}

template <int dim>
void
create_triangulation(Triangulation<dim, dim> &triangulation)
{
  GridGenerator::quarter_hyper_ball(triangulation);
}

template <int dim, int spacedim>
void
check_file() // for dim = spaceim
{
  Triangulation<dim, spacedim> in_tria, out_tria;
  create_triangulation(in_tria);

  // make each cell a different material id
  unsigned int m_id = 0;
  for (const auto &cell : in_tria)
    {
      cell.set_material_id(m_id++);
    }

  // set different boundary ids and output
  unsigned int b_id = 0;
  for (const auto &cell : in_tria)
    {
      for (const auto f : cell.face_indices())
        {
          if (cell.face(f)->at_boundary())
            {
              cell.face(f)->set_boundary_id(b_id);
              b_id++;
            }
        }
    }
  in_tria.refine_global(1);

  GridGenerator::convert_hypercube_to_simplex_mesh(in_tria, out_tria);

  // copy manifolds to test global refining
  for (const auto i : in_tria.get_manifold_ids())
    if (i != numbers::flat_manifold_id)
      out_tria.set_manifold(i, in_tria.get_manifold(i));

  // write 2 outputs (total mesh and only surface mesh)
  const auto grid_out = [](const auto &tria,
                           const bool  surface_mesh_only = false) {
    GridOutFlags::Vtk flags;

    if (surface_mesh_only)
      {
        flags.output_cells         = false;
        flags.output_faces         = true;
        flags.output_edges         = false;
        flags.output_only_relevant = false;
      }

    GridOut grid_out;
    grid_out.set_flags(flags);

    grid_out.write_vtk(tria, deallog.get_file_stream());
  };

  grid_out(out_tria);       // total mesh
  grid_out(out_tria, true); // only surface mesh

  deallog << "OK!" << std::endl;
}


int
main()
{
  initlog();
  // TRIANGULAR ELEMENTS
  // dim = spacedim = 2
  deallog.push(
    "2D: conversion triangulation with quad elements to tri elements: ");
  check_file<2, 2>();
  deallog.pop();

  // TETRAHEDRAL ELEMENTS
  // dim = 2, spacedim = 2
  deallog.push(
    "2D: conversion triangulation with quad elements to tri elements: ");
  check_file<2, 3>();
  deallog.pop();

  // TETRAHEDRAL ELEMENTS
  // dim = spacedim = 3
  deallog.push(
    "3D: conversion triangulation with tet elements to hex elements: ");
  check_file<3, 3>();
  deallog.pop();
}
