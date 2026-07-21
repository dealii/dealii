// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Verify that we correctly handle cells whose material ids are not equal to
// their parent's material id.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include "../grid/tests.h"


template <int dim>
void
test(MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::hyper_cube(basetria);

  basetria.refine_global(2);
  for (auto &cell : basetria.cell_iterators())
    cell->set_material_id(100 * cell->level() + cell->index());

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  for (const auto &cell : tria_pft.cell_iterators())
    deallog << cell << ": " << cell->material_id() << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  const MPI_Comm comm = MPI_COMM_WORLD;

  test<2>(comm);
}
