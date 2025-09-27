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



// test MatrixFree::reinit when some procs have only FE_Nothing and thus zero
// dofs

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim, int fe_degree>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);

  DoFHandler<dim> dof(tria);
  {
    std::vector<types::fe_index> active_indices(tria.n_active_cells());
    for (unsigned int i = 0; i < tria.n_active_cells() / 2; ++i)
      active_indices[i] = 1;
    dof.set_active_fe_indices(active_indices);

    hp::FECollection<dim> fe{FE_Q<dim>(fe_degree), FE_Nothing<dim>(1)};
    dof.distribute_dofs(fe);
  }

  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing reinit for FE_Nothing + " << dof.get_fe().get_name()
          << std::endl;
  deallog << "# local dofs: " << dof.n_locally_owned_dofs() << std::endl;

  using number = double;
  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.mapping_update_flags =
      (dealii::update_values | dealii::update_gradients |
       dealii::update_JxW_values | dealii::update_quadrature_points);
    mf_data.reinit(MappingQ<dim>(fe_degree),
                   dof,
                   constraints,
                   hp::QCollection<dim>(quad, quad),
                   data);
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
