// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test mesh_loop in parallel for GMG level cells

#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/meshworker/mesh_loop.h>

#include "../tests.h"

struct ScratchData
{};
struct CopyData
{
  unsigned     n_cells;
  unsigned int n_own_cells;
  unsigned int n_ghost_cells;

  CopyData()
    : n_cells(0)
    , n_own_cells(0)
    , n_ghost_cells(0)
  {}

  void
  reset()
  {
    n_cells = n_own_cells = n_ghost_cells = 0;
  }
};

using namespace MeshWorker;

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_L(tria);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();


  ScratchData scratch;
  CopyData    copy;

  auto cell = dofh.begin_mg();
  auto endc = dofh.end_mg();

  using Iterator = decltype(cell);

  auto cell_worker = [](const Iterator &cell, ScratchData &s, CopyData &c) {
    deallog << "Cell worker on : " << cell->id()
            << " is_locally_owned_on_level? "
            << cell->is_locally_owned_on_level() << std::endl;
    ++c.n_cells;
    if (cell->is_locally_owned_on_level())
      ++c.n_own_cells;
    else
      ++c.n_ghost_cells;
  };

  auto boundary_worker =
    [](const Iterator &cell, const unsigned int &f, ScratchData &, CopyData &) {
      deallog << "Boundary worker on : " << cell << ", Face : " << f
              << std::endl;
    };

  auto face_worker = [](const Iterator     &cell,
                        const unsigned int &f,
                        const unsigned int &sf,
                        const Iterator     &ncell,
                        const unsigned int &nf,
                        const unsigned int &nsf,
                        ScratchData        &s,
                        CopyData           &c) {
    deallog << "Face worker on : " << cell << ", Neighbor cell : " << ncell
            << ", Face : " << f << ", Neighbor Face : " << nf
            << ", Subface: " << sf << ", Neighbor Subface: " << nsf
            << ", cell->is_locally_owned_on_level(): "
            << cell->is_locally_owned_on_level()
            << ", neighbor->is_locally_owned_on_level(): "
            << ncell->is_locally_owned_on_level() << std::endl;
  };

  CopyData data;
  auto     copier = [&](const CopyData &c) {
    data.n_cells += c.n_cells;
    data.n_own_cells += c.n_own_cells;
    data.n_ghost_cells += c.n_ghost_cells;
  };

  std::function<void(const Iterator &, ScratchData &, CopyData &)>
    empty_cell_worker;
  std::function<
    void(const Iterator &, const unsigned int &, ScratchData &, CopyData &)>
    empty_boundary_worker;

  auto print_summary = [&]() {
    deallog << "n_cells: " << data.n_cells
            << " n_own_cells: " << data.n_own_cells
            << " n_ghost_cells: " << data.n_ghost_cells << std::endl;
  };


  {
    MPILogInitAll log;
    deallog << "OWN CELLS:" << std::endl << std::endl;

    data.reset();
    mesh_loop(
      cell, endc, cell_worker, copier, scratch, copy, assemble_own_cells);
    print_summary();
  }

  {
    MPILogInitAll log;
    deallog << "GHOST and OWN CELLS" << std::endl << std::endl;

    data.reset();
    mesh_loop(cell,
              endc,
              cell_worker,
              copier,
              scratch,
              copy,
              assemble_own_cells | assemble_ghost_cells);
    print_summary();
  }

  {
    MPILogInitAll log;
    deallog << "GHOST CELLS:" << std::endl << std::endl;

    data.reset();
    mesh_loop(
      cell, endc, cell_worker, copier, scratch, copy, assemble_ghost_cells);
    print_summary();
  }

  {
    MPILogInitAll log;
    deallog << "CELLS and FACES" << std::endl << std::endl;

    data.reset();
    mesh_loop(cell,
              endc,
              cell_worker,
              copier,
              scratch,
              copy,
              assemble_own_cells | assemble_own_interior_faces_once |
                assemble_ghost_faces_once,
              empty_boundary_worker,
              face_worker);
    print_summary();
  }

  {
    MPILogInitAll log;
    deallog << "OWN CELLS AND GHOST FACES ONCE" << std::endl << std::endl;

    data.reset();
    mesh_loop(cell,
              endc,
              cell_worker,
              copier,
              scratch,
              copy,
              assemble_own_cells | assemble_ghost_faces_once,
              empty_boundary_worker,
              face_worker);
    print_summary();
  }

  {
    MPILogInitAll log;
    deallog << "GHOST FACES BOTH" << std::endl << std::endl;

    data.reset();
    mesh_loop(cell,
              endc,
              empty_cell_worker,
              copier,
              scratch,
              copy,
              assemble_ghost_faces_both,
              empty_boundary_worker,
              face_worker);
    print_summary();
  }

  {
    MPILogInitAll log;
    deallog << "BOUNDARY FACES" << std::endl << std::endl;

    data.reset();
    mesh_loop(cell,
              endc,
              empty_cell_worker,
              copier,
              scratch,
              copy,
              assemble_boundary_faces,
              boundary_worker);
    print_summary();
  }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  test<2, 2>();
  test<3, 3>();
}
