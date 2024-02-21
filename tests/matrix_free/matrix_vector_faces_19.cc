// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// similar to matrix_vector_faces_13 (integrate_scatter, renumbered degrees of
// freedom) but testing a vmult_add operation

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"


template <int dim, int fe_degree_>
void
test()
{
  // raise element degree by two to test cubic and quartic shape functions
  // rather than linears and quadratics according to
  // matrix_vector_faces_common.h

  const unsigned int                        fe_degree = fe_degree_ + 2;
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_mesh(tria);

  if (dim == 2)
    tria.refine_global(1);
  {
    typename Triangulation<dim>::active_cell_iterator cell =
      tria.begin_active();
    typename Triangulation<dim>::active_cell_iterator endc    = tria.end();
    unsigned int                                      counter = 0;
    for (; cell != endc; ++cell, ++counter)
      if (cell->is_locally_owned() && counter % 3 == 0)
        cell->set_refine_flag();
    tria.execute_coarsening_and_refinement();
  }

  FE_DGQHermite<dim> fe(fe_degree);
  DoFHandler<dim>    dof(tria);
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

  LinearAlgebra::distributed::Vector<double> in, out, out_dist;

  MatrixFree<dim, double>                          mf_data;
  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);
  data.initialize_mapping = false;

  mf_data.reinit(mapping, dof, constraints, quad, data);

  std::vector<types::global_dof_index> renumbering;
  mf_data.renumber_dofs(renumbering);
  dof.renumber_dofs(renumbering);

  data.initialize_mapping = true;
  mf_data.reinit(mapping, dof, constraints, quad, data);

  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);
  mf_data.initialize_dof_vector(out_dist);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < in.locally_owned_size(); ++i)
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

  MatrixFreeVariant<dim,
                    fe_degree,
                    fe_degree + 1,
                    double,
                    LinearAlgebra::distributed::Vector<double>>
    mf2(mf_data);
  mf2.vmult_add(out_dist, in);

  out_dist -= out;

  double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << ' ';

  // test again, now doing a matrix-vector product and then an add, which
  // should give twice the reference result
  mf2.vmult(out_dist, in);
  mf2.vmult_add(out_dist, in);
  out_dist -= out;
  out_dist -= out;
  diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << diff_norm << std::endl;
}
