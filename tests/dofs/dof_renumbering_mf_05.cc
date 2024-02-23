// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check DoFRenumbering::matrix_free_data_locality on a hypercube mesh in
// parallel for an FE_System with a second DoFHandler in MatrixFree.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  deallog << "Test in " << dim << "D with degree " << degree << std::endl;
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, -1., 1.);
  tria.refine_global(5 - dim);
  FESystem<dim>   fe(FE_Q<dim>(degree), dim);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  FESystem<dim>   fe2(FE_DGQ<dim>(degree), dim);
  DoFHandler<dim> dof2(tria);
  dof2.distribute_dofs(fe2);

  using MatrixFreeType = MatrixFree<dim, double, VectorizedArray<double, 1>>;
  typename MatrixFreeType::AdditionalData mf_data;
  mf_data.tasks_parallel_scheme = MatrixFreeType::AdditionalData::none;
  MatrixFreeType            mf;
  AffineConstraints<double> constraints, constraints2;
  mf.reinit(dealii::MappingQ1<dim>(),
            {&dof2, &dof},
            std::vector<const AffineConstraints<double> *>{&constraints2,
                                                           &constraints},
            dealii::QGauss<1>(2),
            mf_data);

  {
    const auto renumber =
      DoFRenumbering::compute_matrix_free_data_locality(dof, mf);

    deallog << "Renumbering no constraints: " << std::endl;
    for (unsigned int i = 0; i < renumber.size(); ++i)
      {
        deallog << renumber[i] << " ";
        if (i % 16 == 15)
          deallog << std::endl;
      }
    deallog << std::endl;
  }

  DoFRenumbering::matrix_free_data_locality(dof, mf);
  std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
  deallog << "New dof indices on cells: " << std::endl;
  for (const auto &cell : dof.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell->get_dof_indices(dof_indices);
        for (const auto i : dof_indices)
          deallog << i << " ";
        deallog << std::endl;
      }
  deallog << std::endl;

  dof.distribute_dofs(fe);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(dim),
                                           constraints);
  constraints.close();

  mf.reinit(dealii::MappingQ1<dim>(),
            {&dof2, &dof},
            std::vector<const AffineConstraints<double> *>{&constraints2,
                                                           &constraints},
            dealii::QGauss<1>(2),
            mf_data);

  {
    const auto renumber =
      DoFRenumbering::compute_matrix_free_data_locality(dof, mf);

    deallog << "Renumbering Dirichlet constraints: " << std::endl;
    for (unsigned int i = 0; i < renumber.size(); ++i)
      {
        deallog << renumber[i] << " ";
        if (i % 16 == 15)
          deallog << std::endl;
      }
    deallog << std::endl;
  }

  DoFRenumbering::matrix_free_data_locality(dof, mf);
  deallog << "New dof indices on cells: " << std::endl;
  for (const auto &cell : dof.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell->get_dof_indices(dof_indices);
        for (const auto i : dof_indices)
          deallog << i << " ";
        deallog << std::endl;
      }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>(2);
  test<2>(4);
  test<3>(2);
  test<3>(3);
}
