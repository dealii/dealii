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


// Tests whether the `reference_cells` vector inside the triangulation is
// updated in case of refinement (Note: Refinement of pyramids introduces Tets.)


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <mpi.h>

#include "../tests.h"

#include "./simplex_grids.h"

template <typename TriaType>
std::string
get_all_ref_cell_types_from_cells(const TriaType    &tria,
                                  const std::string &pre_post)
{
  std::ostringstream oss;
  oss << "Reference cells " << std::setw(4) << std::left << pre_post
      << " refinement from cells:" << std::setw(14) << " ";
  for (const auto &cell : tria.active_cell_iterators())
    oss << cell->reference_cell().to_string() << ", ";
  return oss.str();
}

template <typename TriaType>
std::string
get_ref_cells_in_tria(const TriaType &tria, const std::string &pre_post)
{
  std::ostringstream oss;
  oss << "Ref-Cells " << std::setw(4) << std::left << pre_post
      << " ref. given by tria.get_reference_cells(): ";
  for (const auto &r : tria.get_reference_cells())
    oss << r.to_string() << ", ";
  return oss.str();
}

void
log_gathered_string(const std::string &local_msg, const MPI_Comm &comm)
{
  const unsigned int mpi_process =
    dealii::Utilities::MPI::this_mpi_process(comm);
  const auto process_logs = dealii::Utilities::MPI::all_gather(comm, local_msg);

  for (unsigned int r = 0; r < process_logs.size(); ++r)
    deallog << r << ": " << process_logs[r] << std::endl;
}

template <unsigned int dim>
void
test_serial(const ReferenceCell<dim> ref_cell)
{
  dealii::Triangulation<dim> tria;
  dealii::GridGenerator::reference_cell(tria, ref_cell);

  const auto log_ref_cells = [&tria](const std::string &pre_post) {
    deallog << " Ref-Cells " << std::setw(4) << std::left << pre_post
            << " ref. given by tria.get_reference_cells(): ";
    for (const auto &r : tria.get_reference_cells())
      deallog << r.to_string() << ", ";
    deallog << std::endl;
  };

  const auto log_cells = [&tria](const std::string &pre_post) {
    deallog << " " << get_all_ref_cell_types_from_cells(tria, pre_post)
            << std::endl;
  };

  log_ref_cells("pre");
  log_cells("pre");

  deallog << std::endl;

  tria.refine_global(1);

  log_ref_cells("post");
  log_cells("post");
}

template <unsigned int dim>
void
test_shared(const ReferenceCell<dim> ref_cell)
{
  dealii::parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  dealii::GridGenerator::reference_cell(tria, ref_cell);

  const MPI_Comm comm = tria.get_mpi_communicator();

  log_gathered_string(get_ref_cells_in_tria(tria, "pre"), comm);
  log_gathered_string(get_all_ref_cell_types_from_cells(tria, "pre"), comm);

  deallog << std::endl;

  tria.refine_global(1);

  log_gathered_string(get_ref_cells_in_tria(tria, "post"), comm);
  log_gathered_string(get_all_ref_cell_types_from_cells(tria, "post"), comm);
}

template <unsigned int dim>
void
test_fully_distributed(const ReferenceCell<dim> ref_cell)
{
  dealii::parallel::fullydistributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  const MPI_Comm comm = tria.get_mpi_communicator();

  unsigned int       refinements = 0;
  const unsigned int group_size  = 40;

  const auto serial_grid_generator =
    [&ref_cell, &refinements](dealii::Triangulation<dim, dim> &tria_serial) {
      dealii::GridGenerator::reference_cell(tria_serial, ref_cell);
      // if (refinements > 0)
      tria_serial.refine_global(refinements);
    };

  const auto serial_grid_partitioner =
    [](dealii::Triangulation<dim, dim> &tria_serial,
       const MPI_Comm                   comm_part,
       const unsigned int) {
      dealii::GridTools::partition_triangulation_zorder(
        dealii::Utilities::MPI::n_mpi_processes(comm_part), tria_serial);
    };

  const auto make_description = [&]() {
    return dealii::TriangulationDescription::Utilities::
      create_description_from_triangulation_in_groups<dim, dim>(
        serial_grid_generator,
        serial_grid_partitioner,
        comm,
        group_size,
        dealii::Triangulation<dim>::none,
        dealii::TriangulationDescription::default_setting);
  };

  // Pre-refinement test
  tria.create_triangulation(make_description());

  std::ostringstream msg;
  log_gathered_string(get_ref_cells_in_tria(tria, "pre"), comm);
  log_gathered_string(get_all_ref_cell_types_from_cells(tria, "pre"), comm);

  deallog << std::endl;

  // Post-refinement test
  tria.clear();
  refinements = 1; // implicitly handed to lambda
  tria.create_triangulation(make_description());

  log_gathered_string(get_ref_cells_in_tria(tria, "post"), comm);
  log_gathered_string(get_all_ref_cell_types_from_cells(tria, "post"), comm);
}

enum class TestModes
{
  Serial,
  FullyDist,
  Shared
};

template <TestModes mode = TestModes::Serial, unsigned int dim>
void
run()
{
  auto [name,
        test] = []() -> std::tuple<std::string, void (*)(ReferenceCell<dim>)> {
    switch (mode)
      {
        case TestModes::FullyDist:
          return {"Fully-Distributed Tria", test_fully_distributed<dim>};
        case TestModes::Shared:
          return {"Shared Tria", test_shared<dim>};
        case TestModes::Serial:
        default:
          return {"Serial Tria", test_serial<dim>};
      }
  }();

  const auto cells = ReferenceCells::get_reference_cells_in_dim<dim>();

  for (size_t i = 0; i < cells.size(); ++i)
    {
      deallog << std::endl
              << " ======================================" << std::endl
              << " " << name << " in " << dim << "D" << std::endl
              << std::endl;
      test(cells[i]);
    }
}

int
main(int argc, char **argv)
{
  mpi_initlog();

  dealii::Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  constexpr unsigned int dim = 3; // currently only pyramids (3D) cause issues.

  if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1)
    run<TestModes::Serial, dim>();
  else
    {
      run<TestModes::Shared, dim>();
      // Distributed Triangulations only support isotopically refined hexes so -
      // as of now - no need to test.
      run<TestModes::FullyDist, dim>();
    }
}
