// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// A test contributed by Praveen C that checks that MeshWorker::loop()
// works correctly when encountering periodic faces while doing
// something on faces.
//
// We make a 3x3 mesh with periodicity in both directions.
//
// Cells are numbered like this.
//
// |-----|-----|-----|
// |  6  |  7  |  8  |
// |-----|-----|-----|
// |  3  |  4  |  5  |
// |-----|-----|-----|
// |  0  |  1  |  2  |
// |-----|-----|-----|
// When you run the code,
//
// ----- Cell no = 0-----
//   Periodic neighbor cell = 2
//   Periodic neighbor cell = 6
// ----- Cell no = 1-----
//   Periodic neighbor cell = 7
// ----- Cell no = 2-----
//   Periodic neighbor cell = 0
//   Periodic neighbor cell = 8
// ----- Cell no = 3-----
//   Periodic neighbor cell = 5
// ----- Cell no = 4-----
// ----- Cell no = 5-----
//   Periodic neighbor cell = 3
// ----- Cell no = 6-----
//   Periodic neighbor cell = 8
//   Periodic neighbor cell = 0
// ----- Cell no = 7-----
//   Periodic neighbor cell = 1
// ----- Cell no = 8-----
//   Periodic neighbor cell = 6
//   Periodic neighbor cell = 2
//
// No of boundary faces = 12
//
// Beginning of mesh_loop:
//
// Interior face: (cell,neigh) = (0,2)
// Interior face: (cell,neigh) = (0,1)
// Interior face: (cell,neigh) = (0,6)
// Interior face: (cell,neigh) = (0,3)
// Interior face: (cell,neigh) = (1,2)
// Interior face: (cell,neigh) = (1,7)
// Interior face: (cell,neigh) = (1,4)
// Interior face: (cell,neigh) = (2,8)
// Interior face: (cell,neigh) = (2,5)
// Interior face: (cell,neigh) = (3,5)
// Interior face: (cell,neigh) = (3,4)
// Interior face: (cell,neigh) = (3,6)
// Interior face: (cell,neigh) = (4,5)
// Interior face: (cell,neigh) = (4,7)
// Interior face: (cell,neigh) = (5,8)
// Interior face: (cell,neigh) = (6,8)
// Interior face: (cell,neigh) = (6,7)
// Interior face: (cell,neigh) = (7,8)
//
// you see that boundary_worker is never called, because there are no
// boundary faces. Faces on the boundary are periodic and they are
// treated as interior faces.
//
// There is a periodic face between cell 0 and cell 2; this face is
// visited only from cell 0, see the output (0,2) is present but (2,0)
// is not present.


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

const int dim = 2;

const std::vector<unsigned int> mesh_size = {3, 3};

typedef parallel::distributed::Triangulation<dim> PTriangulation;

void
test()
{
  PTriangulation triangulation(MPI_COMM_WORLD,
                               Triangulation<dim>::smoothing_on_refinement,
                               parallel::distributed::Triangulation<dim>::
                                 mesh_reconstruction_after_repartitioning);

  const Point<dim> point1(0, 0);
  const Point<dim> point2(1, 1);
  GridGenerator::subdivided_hyper_rectangle(
    triangulation, mesh_size, point1, point2, true);

  typedef typename PTriangulation::cell_iterator Iter;
  std::vector<GridTools::PeriodicFacePair<Iter>> periodicity_vector;

  // Periodic along x
  {
    deallog << "Collect periodic faces along x\n";
    GridTools::collect_periodic_faces(
      triangulation, 0, 1, 0, periodicity_vector);
  }

  // Periodic along y
  {
    deallog << "Collect periodic faces along y\n";
    GridTools::collect_periodic_faces(
      triangulation, 2, 3, 1, periodicity_vector);
  }

  {
    deallog << "Applying periodicity\n";
    triangulation.add_periodicity(periodicity_vector);
  }

  int count = 0;
  for (auto &cell : triangulation.active_cell_iterators())
    {
      cell->set_user_index(count++);
    }

  int no_boundary_faces = 0;
  for (auto &cell : triangulation.active_cell_iterators())
    {
      deallog << "----- Cell no = " << cell->user_index() << "----- \n";
      for (auto f : cell->face_indices())
        if (cell->face(f)->at_boundary())
          {
            ++no_boundary_faces;
            deallog << "\t Periodic neighbor cell = "
                    << cell->periodic_neighbor(f)->user_index() << "\n";
          }
    }

  deallog << "No of boundary faces = " << no_boundary_faces << "\n";

  auto assemble_flags = MeshWorker::assemble_own_interior_faces_once;
  assemble_flags |= MeshWorker::assemble_ghost_faces_once;
  assemble_flags |= MeshWorker::assemble_boundary_faces;

  const auto iterator_range =
    filter_iterators(triangulation.active_cell_iterators(),
                     IteratorFilters::LocallyOwnedCell());

  using ScratchData = MeshWorker::ScratchData<dim, dim>;
  using CopyData    = MeshWorker::CopyData<1, 1, 1>;
  using Iterator    = decltype(triangulation.begin_active());

  const auto  fe = FE_Q<dim>(1);
  ScratchData scratch_data(fe, QGauss<dim>(1), update_values);
  CopyData    copy_data(fe.dofs_per_cell);

  auto face_worker = [](const Iterator    &cell,
                        const unsigned int f,
                        const unsigned int sf,
                        const Iterator    &ncell,
                        const unsigned int nf,
                        const unsigned int nsf,
                        ScratchData       &scratch_data,
                        CopyData          &copy_data) {
    deallog << "Interior face: (cell,neigh) = (" << cell->user_index() << ","
            << ncell->user_index() << ")" << std::endl;
  };

  auto boundary_worker = [](const Iterator    &cell,
                            const unsigned int f,
                            ScratchData       &scratch_data,
                            CopyData          &copy_data) {
    deallog << "Boundary face: cell,f: " << cell->user_index() << " " << f
            << std::endl;
  };

  auto copier = [](const CopyData &cd) {};

  deallog << "Beginning of mesh_loop:\n\n";
  MeshWorker::mesh_loop(iterator_range,
                        nullptr,
                        copier,
                        scratch_data,
                        copy_data,
                        assemble_flags,
                        boundary_worker,
                        face_worker);
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  test();
}
