// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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



// similar to matrix_vector_faces_09 but including MPI and a few different
// variants for the loop

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

  mf_data.reinit(mapping, dof, constraints, quad, data);

  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);
  mf_data.initialize_dof_vector(out_dist);

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

  MatrixFreeVariant<dim,
                    fe_degree,
                    fe_degree + 1,
                    double,
                    LinearAlgebra::distributed::Vector<double>>
    mf2(mf_data, false);
  mf2.vmult(out_dist, in);
  mf2.vmult(out_dist, in);

  out_dist -= out;

  const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << std::endl;
}
