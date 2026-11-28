// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check that DoFAccessor::get_dof_indices() on faces of artificial cells
// returns numbers::invalid_dof_index for all degrees of freedom if any of the
// vertices, the lines, or the face itself is only adjacent to artificial
// cells and does not store degrees of freedom. FE_Q(2) has exactly one degree
// of freedom on every vertex, line, and quad, so the expected indices are
// those of the vertices, then the lines (in 3d), then the face itself.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  std::vector<unsigned int> repetitions(dim, 1);
  repetitions[0] = 8;
  Point<dim> top_right;
  for (unsigned int d = 0; d < dim; ++d)
    top_right[d] = (d == 0 ? 8. : 1.);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<dim>(),
                                            top_right);
  tria.refine_global(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(2));

  unsigned int n_partial         = 0;
  unsigned int n_without_storage = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_artificial())
      for (const auto f : cell->face_indices())
        {
          const auto face = cell->face(f);

          std::vector<types::global_dof_index> expected;
          for (const auto v : face->vertex_indices())
            expected.push_back(face->vertex_dof_index(v, 0));
          if constexpr (dim == 3)
            for (const auto l : face->line_indices())
              expected.push_back(face->line(l)->dof_index(0));
          expected.push_back(face->dof_index(0));

          const auto n_invalid = std::count(expected.begin(),
                                            expected.end(),
                                            numbers::invalid_dof_index);
          if (n_invalid == static_cast<long>(expected.size()))
            ++n_without_storage;
          else if (n_invalid > 0)
            {
              ++n_partial;
              std::fill(expected.begin(),
                        expected.end(),
                        numbers::invalid_dof_index);
            }

          std::vector<types::global_dof_index> dof_indices(
            dof_handler.get_fe().n_dofs_per_face());
          face->get_dof_indices(dof_indices);
          AssertThrow(dof_indices == expected, ExcInternalError());
        }

  deallog << "faces of artificial cells with partial storage: " << n_partial
          << ", without storage: " << n_without_storage << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
