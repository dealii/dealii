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



// similar to matrix_vector_faces_13 (renumbering of degrees of freedom to
// better vectorized access, gather_evaluate and integrate_scatter) but with
// block vectors

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_block_vector.h>

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

  FE_DGQHermite<dim>        fe(fe_degree);
  DoFHandler<dim>           dof(tria);
  AffineConstraints<double> constraints;
  constraints.close();

  for (unsigned int test = 0; test < 2; ++test)
    {
      tria.refine_global(1);

      dof.distribute_dofs(fe);

      deallog << "Testing " << dof.get_fe().get_name();
      deallog << " on " << dof.n_dofs() << " DoFs";
      deallog << std::endl;

      MappingQ<dim> mapping(dof.get_fe().degree + 1);

      LinearAlgebra::distributed::BlockVector<double> in(1);
      LinearAlgebra::distributed::BlockVector<double> out(1);
      LinearAlgebra::distributed::BlockVector<double> out_dist(1);

      MatrixFree<dim, double>                          mf_data;
      const QGauss<1>                                  quad(fe_degree + 1);
      typename MatrixFree<dim, double>::AdditionalData data;
      data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
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

      mf_data.initialize_dof_vector(in.block(0));
      mf_data.initialize_dof_vector(out.block(0));
      mf_data.initialize_dof_vector(out_dist.block(0));

      // Set random seed for reproducibility
      Testing::srand(42);
      for (unsigned int i = 0; i < in.block(0).locally_owned_size(); ++i)
        {
          const double entry           = Testing::rand() / (double)RAND_MAX;
          in.block(0).local_element(i) = entry;
        }

      MatrixFreeTest<dim,
                     fe_degree,
                     fe_degree + 1,
                     double,
                     LinearAlgebra::distributed::BlockVector<double>>
        mf(mf_data);
      mf.vmult(out, in);

      MatrixFreeVariant<dim,
                        fe_degree,
                        fe_degree + 1,
                        double,
                        LinearAlgebra::distributed::BlockVector<double>>
        mf2(mf_data);
      mf2.vmult(out_dist, in);

      out_dist -= out;

      double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
      deallog << "Norm of difference:          " << diff_norm << ' ';

      // test again, now doing matrix-vector product twice
      mf2.vmult(out_dist, in);
      mf2.vmult(out_dist, in);
      out_dist -= out;
      diff_norm = out_dist.linfty_norm() / out.linfty_norm();
      deallog << diff_norm << std::endl;
    }
}
