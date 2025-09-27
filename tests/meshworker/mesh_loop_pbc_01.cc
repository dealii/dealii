// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test that we can use mesh_loop with DG and periodic boundary conditions

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/meshworker/mesh_loop.h>

#include "../tests.h"

template <int dim>
struct ScratchData
{};

struct CopyData
{};

template <int dim>
void
test()
{
  Triangulation<dim> triangulation;

  const FE_DGQ<dim> fe(1);
  DoFHandler<dim>   dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0., 1., true);
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  GridTools::collect_periodic_faces(triangulation, 0, 1, 0, periodicity_vector);
  triangulation.add_periodicity(periodicity_vector);

  triangulation.refine_global(4);
  dof_handler.distribute_dofs(fe);
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator & /*cell*/,
                               ScratchData<dim> & /*scratch_data*/,
                               CopyData & /*copy_data*/) {};

  const auto boundary_worker = [&](const Iterator & /*cell*/,
                                   const unsigned int & /*face_no*/,
                                   ScratchData<dim> & /*scratch_data*/,
                                   CopyData & /*copy_data*/) {};

  const auto face_worker = [&](const Iterator & /*cell*/,
                               const unsigned int & /*f*/,
                               const unsigned int & /*sf*/,
                               const Iterator & /*ncell*/,
                               const unsigned int & /*nf*/,
                               const unsigned int & /*nsf*/,
                               ScratchData<dim> & /*scratch_data*/,
                               CopyData & /*copy_data*/) {};

  const auto copier = [&](const CopyData & /*c*/) {};

  ScratchData<dim> scratch_data;
  CopyData         copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(),
                        dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
  return 0;
}
