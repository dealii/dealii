// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Step-02 on a simplex mesh. Following incompatible modifications had to be
// made:
//  - The function GridGenerator::hyper_shell() is only working for hypercube
//    meshes. Use this function to create a temporary mesh and convert it to
//    simplex mesh.
//  - Local refinement is replaced by global refinement.


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <fstream>

#include "../tests.h"


void
make_grid(Triangulation<2> &triangulation)
{
  Triangulation<2> triangulation_temp;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation_temp, center, inner_radius, outer_radius, 5);

  GridGenerator::convert_hypercube_to_simplex_mesh(triangulation_temp,
                                                   triangulation);
  for (const auto i : triangulation_temp.get_manifold_ids())
    if (i != numbers::flat_manifold_id)
      triangulation.set_manifold(i, triangulation_temp.get_manifold(i));

  triangulation.refine_global(); // WARNING: no local refinement is performed
}

void
distribute_dofs(DoFHandler<2> &dof_handler)
{
  const FE_SimplexP<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());

  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  sparsity_pattern.print_svg(deallog.get_file_stream());
}

void
renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  sparsity_pattern.print_svg(deallog.get_file_stream());
}

int
main()
{
  initlog();

  Triangulation<2> triangulation;
  make_grid(triangulation);

  DoFHandler<2> dof_handler(triangulation);

  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);
}
