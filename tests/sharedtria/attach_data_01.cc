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



// Check the register_data_attach()/notify_ready_to_unpack() mechanism on a
// parallel::shared::Triangulation with artificial cells:
//
//  - The pack and unpack callbacks must only be invoked on cells owned by
//    the present MPI rank (for a to-be-coarsened parent cell / a refined
//    parent cell: on the rank owning the first child). In particular they
//    must never be invoked on artificial cells, on which DoF indices
//    cannot be queried.
//
//  - The attached data must arrive on the new owner of a cell even if the
//    partition changes during execute_coarsening_and_refinement(). We
//    verify this by attaching the local DoF values of a ghosted vector
//    that interpolates an affine function: after refinement and
//    coarsening, the unpacked values have to coincide with the function
//    values at the cell's vertices, on every rank the cell may have been
//    reassigned to.
//
//  - Furthermore, refinement and coarsening flags need only be set on
//    locally owned cells, and calling prepare_coarsening_and_refinement()
//    directly must synchronize them over all ranks before flag smoothing:
//    otherwise coarsening flags of families straddling a partition
//    boundary are deleted and the resulting mesh depends on the
//    partitioning. The cell counts logged by this test therefore have to
//    agree between the mpirun=1 and mpirun=3 output files.

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/cell_status.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"


// An affine function that is reproduced exactly by FE_Q(1) elements on
// every mesh, so that the transferred data can be verified against it:
template <int dim>
double
affine_function(const Point<dim> &p)
{
  double value = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    value += (2. + d) * p[d];
  return value;
}


template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    /*allow_artificial_cells*/ true,
    parallel::shared::Triangulation<dim>::partition_zorder);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(2);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Set up a ghosted vector interpolating the affine function:

  LinearAlgebra::distributed::Vector<double> vector(
    dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler),
    MPI_COMM_WORLD);

  std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell->get_dof_indices(dof_indices);
        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
          if (vector.in_local_range(dof_indices[i]))
            vector(dof_indices[i]) = affine_function<dim>(cell->vertex(i));
      }
  vector.update_ghost_values();

  // Set refinement and coarsening flags on locally owned cells only:

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const auto center = cell->center();

        bool refine = true;
        for (unsigned int d = 0; d < dim; ++d)
          refine = refine && (center[d] < 0.25);

        if (refine)
          cell->set_refine_flag();
        else if (center[0] > 0.5)
          cell->set_coarsen_flag();
      }

  tria.prepare_coarsening_and_refinement();

  // After prepare_coarsening_and_refinement() every rank must hold the
  // complete, partition-independent set of flags:

  unsigned int n_refine_flags  = 0;
  unsigned int n_coarsen_flags = 0;
  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->refine_flag_set())
        ++n_refine_flags;
      if (cell->coarsen_flag_set())
        ++n_coarsen_flags;
    }
  deallog << "flags after prepare: refine=" << n_refine_flags
          << " coarsen=" << n_coarsen_flags << std::endl;
  AssertThrow(Utilities::MPI::max(n_refine_flags, MPI_COMM_WORLD) ==
                Utilities::MPI::min(n_refine_flags, MPI_COMM_WORLD),
              ExcInternalError());
  AssertThrow(Utilities::MPI::max(n_coarsen_flags, MPI_COMM_WORLD) ==
                Utilities::MPI::min(n_coarsen_flags, MPI_COMM_WORLD),
              ExcInternalError());

  // Pack the local DoF values of the ghosted vector on every locally
  // owned cell:

  unsigned int n_packed_persist = 0;
  unsigned int n_packed_refine  = 0;
  unsigned int n_packed_coarsen = 0;

  const auto pack_callback =
    [&](const typename Triangulation<dim>::cell_iterator &cell,
        const CellStatus status) -> std::vector<char> {
    const typename DoFHandler<dim>::cell_iterator dof_cell(&tria,
                                                           cell->level(),
                                                           cell->index(),
                                                           &dof_handler);

    Vector<double> values(fe.dofs_per_cell);
    switch (status)
      {
        case CellStatus::cell_will_persist:
          ++n_packed_persist;
          AssertThrow(cell->is_locally_owned(), ExcInternalError());
          dof_cell->get_dof_values(vector, values);
          break;
        case CellStatus::cell_will_be_refined:
          ++n_packed_refine;
          AssertThrow(cell->is_locally_owned(), ExcInternalError());
          dof_cell->get_dof_values(vector, values);
          break;
        case CellStatus::children_will_be_coarsened:
          ++n_packed_coarsen;
          AssertThrow(cell->child(0)->is_locally_owned(), ExcInternalError());
          dof_cell->get_interpolated_dof_values(vector, values);
          break;
        default:
          AssertThrow(false, ExcInternalError());
      }

    const std::vector<double> buffer(values.begin(), values.end());
    return Utilities::pack(buffer, /*allow_compression=*/false);
  };

  const unsigned int handle =
    tria.register_data_attach(pack_callback,
                              /*returns_variable_size_data=*/false);
  deallog << "handle=" << handle << std::endl;

  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  // Unpack on the new mesh - with a new partition - and verify the
  // received values against the affine function:

  unsigned int n_unpacked_persist = 0;
  unsigned int n_unpacked_refine  = 0;
  unsigned int n_unpacked_coarsen = 0;
  unsigned int n_mismatches       = 0;

  const auto unpack_callback =
    [&](const typename Triangulation<dim>::cell_iterator &cell,
        const CellStatus                                  status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &data_range) {
      switch (status)
        {
          case CellStatus::cell_will_persist:
            ++n_unpacked_persist;
            AssertThrow(cell->is_locally_owned(), ExcInternalError());
            break;
          case CellStatus::cell_will_be_refined:
            // We are handed the (now refined) parent cell:
            ++n_unpacked_refine;
            AssertThrow(cell->child(0)->is_locally_owned(), ExcInternalError());
            break;
          case CellStatus::children_will_be_coarsened:
            ++n_unpacked_coarsen;
            AssertThrow(cell->is_locally_owned(), ExcInternalError());
            break;
          default:
            AssertThrow(false, ExcInternalError());
        }

      const auto values =
        Utilities::unpack<std::vector<double>>(data_range.begin(),
                                               data_range.end(),
                                               /*allow_compression=*/false);
      AssertThrow(values.size() == fe.dofs_per_cell, ExcInternalError());

      // Independent of the status, the packed values are the values of
      // the affine function at the vertices of the cell we are handed:
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        if (std::abs(values[i] - affine_function<dim>(cell->vertex(i))) >
            1.e-10)
          ++n_mismatches;
    };

  tria.notify_ready_to_unpack(handle, unpack_callback);

  // Every cell of the old and new mesh, respectively, must have been
  // visited exactly once over all ranks:

  deallog << "packed:   persist="
          << Utilities::MPI::sum(n_packed_persist, MPI_COMM_WORLD)
          << " refine=" << Utilities::MPI::sum(n_packed_refine, MPI_COMM_WORLD)
          << " coarsen="
          << Utilities::MPI::sum(n_packed_coarsen, MPI_COMM_WORLD) << std::endl;
  deallog << "unpacked: persist="
          << Utilities::MPI::sum(n_unpacked_persist, MPI_COMM_WORLD)
          << " refine="
          << Utilities::MPI::sum(n_unpacked_refine, MPI_COMM_WORLD)
          << " coarsen="
          << Utilities::MPI::sum(n_unpacked_coarsen, MPI_COMM_WORLD)
          << std::endl;
  deallog << "mismatches: " << Utilities::MPI::sum(n_mismatches, MPI_COMM_WORLD)
          << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
