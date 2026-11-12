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

// Check fixed- and variable-size cell data on unequal partitions, including
// empty ranks. Serialization offsets count cells, not bytes.

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_description.h>

#include <fstream>
#include <vector>

#include "../tests.h"

std::vector<char>
payload(const unsigned int material, const bool variable_size)
{
  return std::vector<char>(variable_size ? material + 1 : 3, material + 7);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize initialization(argc, argv, 1);
  const auto                       communicator = MPI_COMM_WORLD;
  const auto rank    = Utilities::MPI::this_mpi_process(communicator);
  const auto n_ranks = Utilities::MPI::n_mpi_processes(communicator);
  if (rank == 0)
    initlog();

  Triangulation<1> serial;
  GridGenerator::subdivided_hyper_cube(serial, 5);
  for (const auto &cell : serial.active_cell_iterators())
    {
      cell->set_material_id(cell->index() + 1);
      cell->set_subdomain_id(n_ranks == 1 || cell->index() < 3 ? 0 : 1);
    }
  parallel::fullydistributed::Triangulation<1> mesh(communicator);
  mesh.create_triangulation(
    TriangulationDescription::Utilities::create_description_from_triangulation(
      serial, communicator));
  AssertThrow(mesh.n_locally_owned_active_cells() ==
                (n_ranks == 1 ? 5 : (rank == 0 ? 3 : (rank == 1 ? 2 : 0))),
              ExcInternalError());

  for (const bool variable_size : {false, true})
    mesh.register_data_attach(
      [variable_size](const auto &cell, const CellStatus) {
        return payload(cell->material_id(), variable_size);
      },
      variable_size);
  mesh.save("checkpoint");

  // Using the same wrong offset for saving and loading can hide the bug.
  // The file must contain five cell sizes followed by 2+3+4+5+6 payload bytes.
  if (rank == 0)
    {
      std::ifstream file("checkpoint_variable.data",
                         std::ios::binary | std::ios::ate);
      AssertThrow(file.tellg() == 5 * sizeof(unsigned int) + 20,
                  ExcMessage("Checkpoint contains gaps between cell data."));
    }

  parallel::fullydistributed::Triangulation<1> restored(communicator);
  restored.load("checkpoint");
  for (const bool variable_size : {false, true})
    {
      const auto handle = restored.register_data_attach(
        [](const auto &, const CellStatus) { return std::vector<char>(); },
        variable_size);
      unsigned int unpacked = 0;
      restored.notify_ready_to_unpack(
        handle, [&](const auto &cell, const CellStatus, const auto &data) {
          AssertThrow(std::vector<char>(data.begin(), data.end()) ==
                        payload(cell->material_id(), variable_size),
                      ExcMessage("Checkpoint changed cell data."));
          ++unpacked;
        });
      AssertThrow(unpacked == restored.n_locally_owned_active_cells(),
                  ExcInternalError());
      if (rank == 0)
        deallog << (variable_size ? "variable" : "fixed") << " data: OK"
                << std::endl;
    }
}
