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
void
log_reference_cells(const TriaType &tria, const std::string &prefix)
{
  deallog << prefix;
  for (const auto &r : tria.get_reference_cells())
    deallog << r.to_string() << ", ";
  deallog << std::endl;
}

template <typename TriaType>
std::string
get_active_cell_ref_types(const TriaType &tria)
{
  std::ostringstream oss;
  for (const auto &cell : tria.active_cell_iterators())
    {
      if constexpr (std::is_same_v<
                      TriaType,
                      dealii::parallel::fullydistributed::Triangulation<3>>)
        {
          if (cell->is_artificial())
            continue;
        }
      oss << cell->reference_cell().to_string() << ", ";
    }
  return oss.str();
}

void
log_gathered_string(const std::string &local_msg,
                    const MPI_Comm    &comm,
                    const std::string &header)
{
  const unsigned int mpi_process =
    dealii::Utilities::MPI::this_mpi_process(comm);
  const auto process_logs = dealii::Utilities::MPI::all_gather(comm, local_msg);

  if (mpi_process == 0)
    {
      if (!header.empty())
        deallog << header << std::endl;
      for (unsigned int r = 0; r < process_logs.size(); ++r)
        deallog << r << ": " << process_logs[r] << std::endl;
    }
}

template <unsigned int dim>
void
test_serial(const ReferenceCell<dim> ref_cell)
{
  dealii::Triangulation<dim> tria;
  dealii::GridGenerator::reference_cell(tria, ref_cell);

  log_reference_cells(
    tria, "Ref-Cells pre ref. given by tria.get_reference_cells():  ");
  deallog << "Reference cells pre refinement from cells:" << std::setw(15)
          << " " << get_active_cell_ref_types(tria) << std::endl
          << std::endl;

  tria.refine_global(1);

  log_reference_cells(
    tria, "Ref-Cells post ref. given by tria.get_reference_cells(): ");
  deallog << "Reference cells post refinement from cells:" << std::setw(14)
          << " " << get_active_cell_ref_types(tria) << std::endl;
}

template <unsigned int dim>
void
test_shared(const ReferenceCell<dim> ref_cell)
{
  dealii::parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  dealii::GridGenerator::reference_cell(tria, ref_cell);

  const MPI_Comm comm = tria.get_mpi_communicator();

  std::ostringstream msg;
  msg << "Ref-Cells pre ref. given by tria.get_reference_cells():  ";
  for (const auto &r : tria.get_reference_cells())
    msg << r.to_string() << ", ";
  log_gathered_string(msg.str(), comm, "");

  msg.str("");
  msg << "Reference cells pre refinement from cells:" << std::setw(15) << " "
      << get_active_cell_ref_types(tria);
  log_gathered_string(msg.str(), comm, "");

  tria.refine_global(1);

  msg.str("");
  msg << "Ref-Cells post ref. given by tria.get_reference_cells(): ";
  for (const auto &r : tria.get_reference_cells())
    msg << r.to_string() << ", ";
  log_gathered_string(msg.str(), comm, " ");

  msg.str("");
  msg << "Reference cells post refinement from cells:" << std::setw(14) << " "
      << get_active_cell_ref_types(tria);
  log_gathered_string(msg.str(), comm, "");
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
  msg << "Ref-Cells pre ref. given by tria.get_reference_cells():  ";
  for (const auto &r : tria.get_reference_cells())
    msg << r.to_string() << ", ";
  log_gathered_string(msg.str(), comm, "");

  msg.str("");
  msg << "Reference cells pre refinement from cells:" << std::setw(15) << " "
      << get_active_cell_ref_types(tria);
  log_gathered_string(msg.str(), comm, "");

  // Post-refinement test
  tria.clear();
  refinements = 1; // implicitly handed to lambda
  tria.create_triangulation(make_description());

  msg.str("");
  msg << "Ref-Cells post ref. given by tria.get_reference_cells(): ";
  for (const auto &r : tria.get_reference_cells())
    msg << r.to_string() << ", ";
  log_gathered_string(msg.str(), comm, " ");

  msg.str("");
  msg << "Reference cells post refinement from cells:" << std::setw(14) << " "
      << get_active_cell_ref_types(tria);
  log_gathered_string(msg.str(), comm, "");
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

  const unsigned int mpi_process =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  for (size_t i = 0; i < cells.size(); ++i)
    {
      if (mpi_process == 0)
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
  initlog();

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
