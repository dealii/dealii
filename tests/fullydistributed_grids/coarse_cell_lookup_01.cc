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


// Test the mapping from coarse-cell ids to coarse-cell indices of
// parallel::fullydistributed::Triangulation. Rank 1 only owns the first and
// the last coarse cell of a long strip of cells, so its locally relevant
// coarse-cell ids are far apart (exercising the sorted-list lookup), whereas
// rank 0 owns a contiguous block (exercising the dense lookup table). On every
// rank, query every coarse-cell id of the global mesh, including ids that are
// not locally relevant.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <set>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  const unsigned int n_coarse_cells = 64;

  Triangulation<2> basetria;
  GridGenerator::subdivided_hyper_rectangle(basetria,
                                            {n_coarse_cells, 1},
                                            Point<2>(0, 0),
                                            Point<2>(n_coarse_cells, 1));

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(comm);
  for (const auto &cell : basetria.active_cell_iterators())
    if (n_procs > 1 &&
        (cell->index() == 0 || cell->index() == n_coarse_cells - 1))
      cell->set_subdomain_id(1);
    else
      cell->set_subdomain_id(0);

  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  parallel::fullydistributed::Triangulation<2> tria(comm);
  tria.create_triangulation(construction_data);

  std::set<types::coarse_cell_id> relevant_ids;
  for (const auto &cell : tria.cell_iterators_on_level(0))
    relevant_ids.insert(cell->id().get_coarse_cell_id());

  deallog << "locally relevant coarse cells:";
  for (const auto id : relevant_ids)
    deallog << ' ' << id;
  deallog << std::endl;

  unsigned int n_errors = 0;
  for (types::coarse_cell_id id = 0; id < n_coarse_cells; ++id)
    {
      const CellId cell_id(id, {});
      const bool   contained = tria.contains_cell(cell_id);
      if (contained != (relevant_ids.count(id) > 0))
        {
          deallog << "contains_cell() is wrong for coarse cell " << id
                  << std::endl;
          ++n_errors;
        }

      if (contained && tria.create_cell_iterator(cell_id)->id() != cell_id)
        {
          deallog << "create_cell_iterator() is wrong for coarse cell " << id
                  << std::endl;
          ++n_errors;
        }
    }

  // Test an id outside the range
  if (tria.contains_cell(CellId(n_coarse_cells + 10, {})))
    ++n_errors;

  deallog << "errors: " << n_errors << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test(MPI_COMM_WORLD);
}
