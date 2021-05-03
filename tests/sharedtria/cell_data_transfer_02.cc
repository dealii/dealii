// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// CellDataTransfer test by refinement
// for parallel::shared::Triangulation with artificial cells


#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.templates.h>

#include <deal.II/numerics/cell_data_transfer.templates.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  // ------ setup ------
  MPI_Comm                             mpi_communicator(MPI_COMM_WORLD);
  parallel::shared::Triangulation<dim> tria(mpi_communicator,
                                            ::Triangulation<dim>::none,
                                            true);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  // ----- prepare -----
  // set refinement/coarsening flags
  for (auto &cell : tria.active_cell_iterators())
    {
      if (cell->id().to_string() == "0_1:0")
        cell->set_refine_flag();
      else if (((dim == 1) && (cell->parent()->id().to_string() == "1_0:")) ||
               ((dim == 2) && (cell->parent()->id().to_string() == "3_0:")) ||
               ((dim == 3) && (cell->parent()->id().to_string() == "7_0:")))
        cell->set_coarsen_flag();
    }

  // ----- gather -----
  // store parent id of all locally owned cells
  Vector<PetscScalar> cell_ids_pre(tria.n_active_cells());
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const std::string  parent_cellid = cell->parent()->id().to_string();
        const unsigned int parent_coarse_cell_id =
          static_cast<unsigned int>(std::stoul(parent_cellid));
        cell_ids_pre(cell->active_cell_index()) = parent_coarse_cell_id;
      }

  // distribute local vector (as presented in step-18)
  PETScWrappers::MPI::Vector distributed_cell_ids_pre(
    mpi_communicator,
    tria.n_active_cells(),
    tria.n_locally_owned_active_cells());

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      distributed_cell_ids_pre(cell->active_cell_index()) =
        cell_ids_pre(cell->active_cell_index());
  distributed_cell_ids_pre.compress(VectorOperation::insert);

  cell_ids_pre = distributed_cell_ids_pre;

  // output initial situation
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "cellid=" << cell->id() << " parentid="
              << std::real(cell_ids_pre(cell->active_cell_index()));
      if (cell->coarsen_flag_set())
        deallog << " coarsening";
      else if (cell->refine_flag_set())
        deallog << " refining";
      deallog << std::endl;
    }

  // ----- transfer -----
  CellDataTransfer<dim, spacedim, Vector<PetscScalar>> cell_data_transfer(tria);

  cell_data_transfer.prepare_for_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  Vector<PetscScalar> cell_ids_post(tria.n_active_cells());
  cell_data_transfer.unpack(cell_ids_pre, cell_ids_post);

  // ------ verify ------
  // check if all children adopted the correct id
  for (auto &cell : tria.active_cell_iterators())
    deallog << "cellid=" << cell->id() << " parentid="
            << std::real(cell_ids_post(cell->active_cell_index())) << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("1d");
  test<1, 1>();
  deallog.pop();
  deallog.push("2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d");
  test<3, 3>();
  deallog.pop();
}
