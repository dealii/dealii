// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test MassOperator initialized on a MatrixFree object with several
// components by comparing to another one that only has a single component

#include "../tests.h"

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

template <int dim>
void
test()
{
  typedef double number;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dof, relevant_set);

  ConstraintMatrix constraints_0(relevant_set), constraints_1(relevant_set);
  VectorTools::interpolate_boundary_values(
    dof, 0, Functions::ZeroFunction<dim>(), constraints_0);
  constraints_0.close();
  constraints_1.close();

  std::shared_ptr<MatrixFree<dim, number>> mf_data_0(
    new MatrixFree<dim, number>());
  std::shared_ptr<MatrixFree<dim, number>> mf_data_1(
    new MatrixFree<dim, number>());
  std::shared_ptr<MatrixFree<dim, number>> mf_data_combined(
    new MatrixFree<dim, number>());
  const QGauss<1>                                  quad(2);
  typename MatrixFree<dim, number>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
  mf_data_0->reinit(dof, constraints_0, quad, data);
  mf_data_1->reinit(dof, constraints_1, quad, data);
  {
    std::vector<const DoFHandler<dim>*>  dof_handlers(2, &dof);
    std::vector<const ConstraintMatrix*> constraint(2);
    constraint[0] = &constraints_0;
    constraint[1] = &constraints_1;
    mf_data_combined->reinit(dof_handlers, constraint, quad, data);
  }

  MatrixFreeOperators::
    MassOperator<dim, 1, 2, 1, LinearAlgebra::distributed::Vector<number>>
      mf_0, mf_1, mf_c0, mf_c1;
  mf_0.initialize(mf_data_0);
  mf_1.initialize(mf_data_1);

  mf_c0.initialize(mf_data_combined, std::vector<unsigned int>(1, 0));
  mf_c1.initialize(mf_data_combined, std::vector<unsigned int>(1, 1));

  LinearAlgebra::distributed::Vector<number> in, out, ref;
  mf_data_0->initialize_dof_vector(in);

  for(unsigned int i = 0; i < in.local_size(); ++i)
    in.local_element(i) = random_value<double>();

  mf_c0.initialize_dof_vector(out);
  mf_c0.initialize_dof_vector(ref);

  mf_0.vmult(ref, in);
  mf_c0.vmult(out, in);
  out -= ref;

  deallog << "Error in component 0: " << out.linfty_norm() << std::endl;

  mf_c1.initialize_dof_vector(out);
  mf_c1.initialize_dof_vector(ref);
  mf_1.vmult(ref, in);
  mf_c1.vmult(out, in);
  out -= ref;

  deallog << "Error in component 1: " << out.linfty_norm() << std::endl;
}

int
main(int argc, char** argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();
  test<2>();
}
