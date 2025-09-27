// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test case with p4est with a generated mesh where the orientation on the
// lines does not match and we would previously run into a bug


#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  constexpr unsigned int                    dim = 2;
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridIn<dim> grid_in;
  grid_in.attach_triangulation(tria);
  std::ifstream in_stream(SOURCE_DIR "/2d_periodicity_01.inp");
  grid_in.read(in_stream, GridIn<dim>::abaqus);

  deallog << "Number of mesh elements: " << tria.n_active_cells() << std::endl;

  const BoundingBox<dim> bounding_box = GridTools::compute_bounding_box(tria);
  for (const auto &cell : tria.cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        for (unsigned int d = 0; d < dim; ++d)
          {
            if (std::abs(face->center()[d] -
                         bounding_box.get_boundary_points().first[d]) < 1e-8)
              face->set_boundary_id(1 + 2 * d);
            else if (std::abs(face->center()[d] -
                              bounding_box.get_boundary_points().second[d]) <
                     1e-8)
              face->set_boundary_id(2 + 2 * d);
          }

  std::vector<
    dealii::GridTools::PeriodicFacePair<Triangulation<dim>::cell_iterator>>
    periodic_face_pairs;

  for (unsigned int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(
      tria, 1 + 2 * d, 2 + 2 * d, d, periodic_face_pairs);

  deallog << "Number of periodic face pairs found: "
          << periodic_face_pairs.size() << std::endl;
  tria.add_periodicity(periodic_face_pairs);

  tria.refine_global();
  deallog << "OK" << std::endl;
}
