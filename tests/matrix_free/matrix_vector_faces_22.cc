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



// extension of matrix_vector_faces_15 (advection, gather_evaluate and
// integrate_scatter) to the vector-valued case in form of multiple components
// with different vector entries but otherwise similar setup.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"



template <int dim, int fe_degree_>
void
test()
{
  // raise element degree by one to test quadratic and cubic shape functions
  // rather than linears and quadratics according to the main function in
  // matrix_vector_faces_common.h

  const unsigned int                        fe_degree = fe_degree_ + 1;
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, -1, 1);

  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_boundary_ids(f);
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  for (unsigned int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(
      tria, 2 * d, 2 * d + 1, d, periodic_faces);
  tria.add_periodicity(periodic_faces);

  tria.refine_global(3 - dim);

  FE_DGQ<dim>               fe(fe_degree);
  FESystem<dim>             fe_system(fe, dim);
  DoFHandler<dim>           dof(tria);
  DoFHandler<dim>           dof_system(tria);
  AffineConstraints<double> constraints;
  constraints.close();

  for (unsigned int test = 0; test < 3; ++test)
    {
      tria.refine_global(1);

      dof.distribute_dofs(fe);
      dof_system.distribute_dofs(fe_system);

      deallog << "Testing " << dof_system.get_fe().get_name();
      deallog << " on " << dof_system.n_dofs() << " DoFs";
      deallog << std::endl;

      // first build a vector-valued system that contains all components in
      // analogy to matrix_vector_faces_14
      LinearAlgebra::distributed::Vector<double> in, out, out_dist;

      MatrixFree<dim, double>                          mf_data;
      const QGauss<1>                                  quad(fe_degree + 2);
      typename MatrixFree<dim, double>::AdditionalData data;
      data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
      data.mapping_update_flags_inner_faces =
        (update_gradients | update_JxW_values);
      data.mapping_update_flags_boundary_faces =
        (update_gradients | update_JxW_values);
      mf_data.reinit(dof_system, constraints, quad, data);

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

      MatrixFreeAdvection<dim,
                          fe_degree,
                          fe_degree + 2,
                          double,
                          LinearAlgebra::distributed::Vector<double>,
                          dim>
        mf2(mf_data, true);

      // now compare the result to a scalar implementation on each of the dim
      // components, using vmult_add in subsequent steps
      mf2.vmult(out, in);
      for (unsigned int d = 0; d < dim; ++d)
        {
          MatrixFreeAdvection<dim,
                              fe_degree,
                              fe_degree + 2,
                              double,
                              LinearAlgebra::distributed::Vector<double>,
                              1>
            mf3(mf_data, true, d);
          if (d == 0)
            mf3.vmult(out_dist, in);
          else
            mf3.vmult_add(out_dist, in);
        }

      out_dist -= out;
      double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
      deallog << "Norm of difference to sum of scalar:     " << diff_norm
              << std::endl;

      if (dim == 3)
        {
          MatrixFreeAdvection<dim,
                              fe_degree,
                              fe_degree + 2,
                              double,
                              LinearAlgebra::distributed::Vector<double>,
                              2>
            mf3(mf_data, true, 0);
          mf3.vmult(out_dist, in);
          MatrixFreeAdvection<dim,
                              fe_degree,
                              fe_degree + 2,
                              double,
                              LinearAlgebra::distributed::Vector<double>,
                              1>
            mf4(mf_data, true, 2);
          mf4.vmult_add(out_dist, in);

          out_dist -= out;
          diff_norm = out_dist.linfty_norm() / out.linfty_norm();
          deallog << "Norm of difference to vector/scalar sum: " << diff_norm
                  << std::endl;
        }
      if (dim == 3)
        {
          MatrixFreeAdvection<dim,
                              fe_degree,
                              fe_degree + 2,
                              double,
                              LinearAlgebra::distributed::Vector<double>,
                              1>
            mf3(mf_data, true, 0);
          mf3.vmult(out_dist, in);
          MatrixFreeAdvection<dim,
                              fe_degree,
                              fe_degree + 2,
                              double,
                              LinearAlgebra::distributed::Vector<double>,
                              2>
            mf4(mf_data, true, 1);
          mf4.vmult_add(out_dist, in);

          out_dist -= out;
          diff_norm = out_dist.linfty_norm() / out.linfty_norm();
          deallog << "Norm of difference to scalar/vector sum: " << diff_norm
                  << std::endl;
        }

      // finally compare to a series of scalar problems
      MatrixFree<dim, double> mf_data_scalar;
      mf_data_scalar.reinit(dof, constraints, quad, data);

      LinearAlgebra::distributed::Vector<double> in_small, out_small, ref_small;
      mf_data_scalar.initialize_dof_vector(in_small);
      mf_data_scalar.initialize_dof_vector(out_small);
      mf_data_scalar.initialize_dof_vector(ref_small);

      MatrixFreeAdvection<dim,
                          fe_degree,
                          fe_degree + 2,
                          double,
                          LinearAlgebra::distributed::Vector<double>,
                          1>
        mf4(mf_data_scalar, true);

      for (unsigned int d = 0; d < dim; ++d)
        {
          std::vector<types::global_dof_index> dof_indices_system(
            fe_system.dofs_per_cell);
          std::vector<types::global_dof_index> dof_indices_scalar(
            fe.dofs_per_cell);
          for (typename DoFHandler<dim>::active_cell_iterator
                 cell_scalar = dof.begin_active(),
                 cell_system = dof_system.begin_active();
               cell_scalar != dof.end();
               ++cell_scalar, ++cell_system)
            if (cell_scalar->is_locally_owned())
              {
                cell_scalar->get_dof_indices(dof_indices_scalar);
                cell_system->get_dof_indices(dof_indices_system);
                for (unsigned int i = 0; i < fe_system.dofs_per_cell; ++i)
                  if (fe_system.system_to_component_index(i).first == d)
                    {
                      in_small(
                        dof_indices_scalar
                          [fe_system.system_to_component_index(i).second]) =
                        in(dof_indices_system[i]);
                      out_small(
                        dof_indices_scalar
                          [fe_system.system_to_component_index(i).second]) =
                        out(dof_indices_system[i]);
                    }
              }

          mf4.vmult(ref_small, in_small);

          out_small -= ref_small;
          diff_norm = out_small.linfty_norm() / out.linfty_norm();
          deallog << "Norm of difference to single scalar:     " << diff_norm
                  << std::endl;
        }
    }
}
