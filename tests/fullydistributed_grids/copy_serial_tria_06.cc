// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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


// Create a serial triangulation, copy it, and create a partitioner
// by matrixfree.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../grid/tests.h"

using namespace dealii;

template <int dim>
void
test(int n_refinements, MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::hyper_cube(basetria);
  basetria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  unsigned int degree = 2;

  // test triangulation
  FE_DGQ<dim>     fe(degree);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);

  // test 1: print the indices of the dofs of each active cell
  {
    std::vector<types::global_dof_index> dofs(Utilities::pow(degree + 1, dim));
    for (auto cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(dofs);

          deallog << std::setw(4) << cell->subdomain_id() << " : ";
          for (const auto dof : dofs)
            deallog << std::setw(4) << dof << ' ';
          deallog << std::endl;
        }
  }

  // test 2: let matrixfree setup all partitioners
  {
    typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags =
      update_gradients | update_JxW_values | update_quadrature_points;
    additional_data.mapping_update_flags_inner_faces =
      update_gradients | update_JxW_values;
    additional_data.mapping_update_flags_boundary_faces =
      update_gradients | update_JxW_values | update_quadrature_points;
    additional_data.hold_all_faces_to_owned_cells = true;
    additional_data.mapping_update_flags_faces_by_cells =
      update_gradients | update_JxW_values | update_quadrature_points;

    MappingQ<dim>             mapping(1);
    QGauss<1>                 quad(degree + 1);
    AffineConstraints<double> constraint;

    MatrixFree<dim, double> matrix_free;
    matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);
  }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const int      n_refinements = 3;
  const MPI_Comm comm          = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(n_refinements, comm);
    deallog.pop();
  }
}
