// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include <deal.II/base/exceptions.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <vector>

#include "../tests.h"

template <typename FaceType>
void
print_face(const FaceType &face, const unsigned int level)
{
  AssertThrow(
    face->get_fe(0).n_dofs_per_face() > 0,
    ExcMessage(
      "Please use a finite element with at least 1 DoF for this test"));

  std::vector<types::global_dof_index> mg_dofs(
    face->get_fe(0).n_dofs_per_face());
  face->get_mg_dof_indices(level, mg_dofs);

  // Print face DoFs
  deallog << mg_dofs[0];
  for (unsigned int i_dof = 1; i_dof < mg_dofs.size(); ++i_dof)
    deallog << ", " << mg_dofs[i_dof];
}

template <int dim>
void
check()
{
  using Tria = Triangulation<dim>;
  // Create triangulation
  Tria tria(Tria::limit_level_difference_at_vertices);

  // Create coarse grid
  GridGenerator::hyper_cube(tria, -1, 1, true);

  // Create periodic boundary along x axis
  std::vector<GridTools::PeriodicFacePair<typename Tria::cell_iterator>>
    face_pairs;

  GridTools::collect_periodic_faces(tria, 0, 1, 0, face_pairs);

  tria.add_periodicity(face_pairs);

  // Refine once and then refine the corner again
  tria.refine_global(1);

  tria.begin_active()->set_refine_flag();
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();

  // Assemble DoFs
  DoFHandler<dim> dof_handler(tria);
  FE_Q<dim>       fe_q(1);

  dof_handler.distribute_dofs(fe_q);
  dof_handler.distribute_mg_dofs();

  // Try to assemble MGConstrainedDofs
  MGConstrainedDoFs constrained_dofs;
  constrained_dofs.initialize(dof_handler);

  // Report which DIM is ok
  deallog << "--- DIM " << dim << " OK ---" << std::endl;

  // Debug printouts:
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    {
      deallog << "Level " << level << ": " << std::endl;

      // Print faces that connect to coarser DoF
      deallog << "  Faces neighboring coarser cells:" << std::endl;
      for (const auto level_cell : dof_handler.cell_iterators_on_level(level))
        {
          // Iterate over all faces
          for (const unsigned int i_face : GeometryInfo<dim>::face_indices())
            {
              if ((level_cell->at_boundary(i_face) &&
                   !level_cell->has_periodic_neighbor(i_face)) ||
                  level_cell->neighbor_or_periodic_neighbor(i_face)->level() ==
                    level_cell->level())
                continue;

              // Get face
              const auto face = level_cell->face(i_face);

              deallog << "    ";
              print_face(face, level);

              // Print whether the neighbor is periodic
              if (level_cell->has_periodic_neighbor(i_face))
                deallog << " via periodic boundary" << std::endl;
              else
                deallog << " via internal face" << std::endl;
            }
        }

      // Print hanging node DoFs
      deallog << "  Refinement edge DoFs: ";
      constrained_dofs.get_refinement_edge_indices(level).print(deallog);
      deallog << std::endl;

      // Print faces that have periodic neighbors
      deallog << "  Faces having periodic neighbors of same level:"
              << std::endl;
      for (const auto level_cell : dof_handler.cell_iterators_on_level(level))
        {
          // Iterate over all faces
          for (const unsigned int i_face : GeometryInfo<dim>::face_indices())
            {
              if (!level_cell->has_periodic_neighbor(i_face) ||
                  level_cell->periodic_neighbor(i_face)->level() !=
                    level_cell->level())
                continue;

              typename DoFHandler<dim>::cell_iterator neighbor_cell =
                level_cell->periodic_neighbor(i_face);

              // Each side only once
              if (level_cell->index() > neighbor_cell->index())
                continue;

              unsigned int neighbor_i_face =
                level_cell->periodic_neighbor_face_no(i_face);

              // Get face
              const auto face          = level_cell->face(i_face);
              const auto neighbor_face = neighbor_cell->face(neighbor_i_face);

              deallog << "    ";
              print_face(face, level);
              deallog << " <-> ";
              print_face(neighbor_face, level);

              // Print whether the neighbor is periodic
              if (level_cell->has_periodic_neighbor(i_face))
                deallog << " via periodic boundary" << std::endl;
              else
                deallog << " via internal face" << std::endl;
            }
        }

      // Print constraint matrix
      deallog << "  Constraint matrix:" << std::endl;
      constrained_dofs.get_level_constraint_matrix(level).print(
        deallog.get_file_stream());
    }
}

int
main()
{
  initlog();

  // check<1>();
  check<2>();
  check<3>();
  return 0;
}
