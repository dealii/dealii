// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// check the setup of ghost faces and elements in case the flag
// MatrixFree::AdditionalData::hold_all_faces_to_owned_cells is set to true

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"



template <int dim, int fe_degree>
void
test()
{
  if (fe_degree > 1)
    return;

  // raise element degree by two to test cubic and quartic shape functions
  // rather than linears and quadratics according to
  // matrix_vector_faces_common.h

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  MappingQ<dim> mapping(dof.get_fe().degree + 1);

  LinearAlgebra::distributed::Vector<double> in, out;

  MatrixFree<dim, double>                          mf_data;
  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);
  data.hold_all_faces_to_owned_cells = true;

  mf_data.reinit(mapping, dof, constraints, quad, data);

  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < in.local_size(); ++i)
    {
      const double entry  = Testing::rand() / (double)RAND_MAX;
      in.local_element(i) = entry;
    }

  MatrixFreeTest<dim,
                 fe_degree,
                 fe_degree + 1,
                 double,
                 LinearAlgebra::distributed::Vector<double>>
    mf(mf_data);
  mf.vmult(out, in);

  deallog << "Norm of result:          " << out.l2_norm() << std::endl;
}
