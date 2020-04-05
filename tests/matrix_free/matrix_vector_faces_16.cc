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



// similar to matrix_vector_faces_13 but using another mesh that triggers a
// more difficult case in terms of the numbers in the mesh

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"


template <int dim, int fe_degree>
void
test()
{
  // only run linears in 3D

  if (dim == 2 || fe_degree > 1)
    return;

  constexpr unsigned int mydim = 3;

  parallel::distributed::Triangulation<mydim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  for (const unsigned int f : GeometryInfo<mydim>::face_indices())
    tria.begin_active()->face(f)->set_all_boundary_ids(f);
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<mydim>::cell_iterator>>
    periodic_faces;
  for (unsigned int d = 0; d < mydim; ++d)
    GridTools::collect_periodic_faces(
      tria, 2 * d, 2 * d + 1, d, periodic_faces);
  tria.add_periodicity(periodic_faces);

  tria.refine_global(4);
  std::cout << "Number of cells: " << tria.n_global_active_cells() << std::endl;

  FE_DGQ<mydim>     fe(1);
  DoFHandler<mydim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  LinearAlgebra::distributed::Vector<double> in, out, in2, out2;

  MatrixFree<mydim, double>                          mf_data;
  const QGauss<1>                                    quad(fe_degree + 1);
  typename MatrixFree<mydim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<mydim, double>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(dof, constraints, quad, data);
  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < in.local_size(); ++i)
    {
      const double entry  = Testing::rand() / (double)RAND_MAX;
      in.local_element(i) = entry;
    }

  MatrixFreeTest<3, 1, 2, double, LinearAlgebra::distributed::Vector<double>>
    mf(mf_data);
  mf.vmult(out, in);

  std::vector<types::global_dof_index> renumbering;
  mf_data.renumber_dofs(renumbering);
  dof.renumber_dofs(renumbering);

  mf_data.reinit(dof, constraints, quad, data);
  mf_data.initialize_dof_vector(in2);
  mf_data.initialize_dof_vector(out2);
  for (unsigned int i = 0; i < in.local_size(); ++i)
    {
      in2(renumbering[i]) = in.local_element(i);
    }

  mf.vmult(out2, in2);

  for (unsigned int i = 0; i < in.local_size(); ++i)
    out2(renumbering[i]) -= out.local_element(i);

  double diff_norm = out2.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << std::endl;
}
