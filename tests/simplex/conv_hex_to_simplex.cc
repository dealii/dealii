/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Test function GridGenerator::convert_hypercube_to_simplex_mesh() in 2D
 * and 3D (dim = spacedim = 2, 3) on a quarter_hyper_ball() triangulation and
 * for dim = 2 < spacedim = 3 on a subdivided_hyper_cube() triangulation.
 * Therefore, the output is written in .vtk format for entire triangulations
 * as well as surface-only triangulations.
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

using namespace dealii;

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
check_file(unsigned int n_refinements = 0) // for dim = spaceim
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

  GridGenerator::convert_hypercube_to_simplex_mesh(in_tria, out_tria);

  // copy manifolds to test global refining
  for (const auto i : in_tria.get_manifold_ids())
    if (i != numbers::flat_manifold_id)
      out_tria.set_manifold(i, in_tria.get_manifold(i));

  out_tria.refine_global(n_refinements);

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

    // Demonstrate a bug with copying manifold ids more clearly:
    if (dim == 3)
      flags.output_only_relevant = false;

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

  // TETRAHEDRAL ELEMENTS
  // dim = spacedim = 3
  deallog.push(
    "3D: conversion triangulation with tet elements to hex elements + refinement: ");
  check_file<3, 3>(1);
  deallog.pop();
}
