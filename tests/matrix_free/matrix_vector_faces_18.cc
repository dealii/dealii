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



// same as matrix_vector_faces_15 except for a larger mesh

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"


template <int dim, int fe_degree>
void
test()
{
  if (fe_degree > 1)
    return;

  const double       h  = 1;
  const double       Lx = 30 * h;
  const double       Li = 10 * h;
  const double       Ly = 6 * h;
  const double       Lz = 4 * h;
  const double       ER = Ly / (Ly - h);
  Triangulation<dim> left, right, bottom, temp;
  Point<dim> left_one, left_two, right_one, right_two, bottom_one, bottom_two;

  left_one[0] = -Li;
  left_one[1] = h;

  left_two[0] = 0;
  left_two[1] = Ly;

  right_one[0] = 0;
  right_one[1] = h;

  right_two[0] = Lx - Li;
  right_two[1] = Ly;

  bottom_one[0] = 0;
  bottom_one[1] = 0;

  bottom_two[0] = Lx - Li;
  bottom_two[1] = h;

  if (dim == 3)
    {
      left_one[2]   = 0;
      right_one[2]  = 0;
      bottom_one[2] = 0;
      left_two[2]   = Lz;
      right_two[2]  = Lz;
      bottom_two[2] = Lz;
    }

  std::vector<unsigned int> refinements_left(dim, 1);
  std::vector<unsigned int> refinements_right(dim, 1);
  std::vector<unsigned int> refinements_bottom(dim, 1);

  refinements_left[1]   = 5;
  refinements_right[1]  = 5;
  refinements_bottom[1] = 1;
  refinements_left[0]   = 10;
  refinements_right[0]  = 20;
  refinements_bottom[0] = 20;

  if (dim == 3)
    {
      refinements_left[2]   = 4;
      refinements_right[2]  = 4;
      refinements_bottom[2] = 4;
    }
  GridGenerator::subdivided_hyper_rectangle(
    left, refinements_left, left_one, left_two, false);
  GridGenerator::subdivided_hyper_rectangle(
    right, refinements_right, right_one, right_two, false);
  GridGenerator::subdivided_hyper_rectangle(
    bottom, refinements_bottom, bottom_one, bottom_two, false);

  GridGenerator::merge_triangulations(left, right, temp);

  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::merge_triangulations(temp, bottom, triangulation);

  if (dim == 3)
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin();
         cell != triangulation.end();
         ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (cell->face(f)->at_boundary())
          {
            if (std::abs(cell->face(f)->center()[2]) < 1e-12)
              cell->face(f)->set_all_boundary_ids(10);
            else if (std::abs(cell->face(f)->center()[2] - Lz) < 1e-12)
              cell->face(f)->set_all_boundary_ids(11);
            else if (std::abs(cell->face(f)->center()[1] - h) < 1e-12)
              cell->face(f)->set_all_boundary_ids(0);
            else if (std::abs(cell->face(f)->center()[1]) < 1e-12)
              cell->face(f)->set_all_boundary_ids(0);
            else if (std::abs(cell->face(f)->center()[0]) < 1e-12)
              cell->face(f)->set_all_boundary_ids(0);
            else if (std::abs(cell->face(f)->center()[1] - Ly) < 1e-12)
              cell->face(f)->set_all_boundary_ids(3);
            else if (std::abs(cell->face(f)->center()[0] + Li) < 1e-12)
              cell->face(f)->set_all_boundary_ids(1);
            else if (std::abs(cell->face(f)->center()[0] - (Lx - Li)) < 1e-12)
              cell->face(f)->set_all_boundary_ids(2);
          }

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;

  if (dim == 3)
    {
      GridTools::collect_periodic_faces(
        triangulation, 10, 11, 2, periodic_faces);
      triangulation.add_periodicity(periodic_faces);
    }

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof_orig(triangulation);
  DoFHandler<dim> dof(triangulation);
  dof_orig.distribute_dofs(fe);
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

  LinearAlgebra::distributed::Vector<double> in, in_orig, out, out_orig;

  // create MatrixFree with DoFHandler in original numbering
  MatrixFree<dim, double>                          mf_data_orig;
  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);
  data.initialize_mapping = true;

  mf_data_orig.reinit(mapping, dof_orig, constraints, quad, data);
  mf_data_orig.initialize_dof_vector(in_orig);
  mf_data_orig.initialize_dof_vector(out_orig);

  // create MatrixFree with renumbered degrees of freedom
  MatrixFree<dim, double> mf_data;
  data.initialize_mapping = false;

  mf_data.reinit(mapping, dof, constraints, quad, data);

  std::vector<types::global_dof_index> renumbering;
  mf_data.renumber_dofs(renumbering);
  dof.renumber_dofs(renumbering);

  data.initialize_mapping = true;
  mf_data.reinit(mapping, dof, constraints, quad, data);

  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < in_orig.local_size(); ++i)
    {
      const double entry       = Testing::rand() / (double)RAND_MAX;
      in_orig.local_element(i) = entry;
      in(renumbering[i])       = entry;
    }

  MatrixFreeAdvection<dim,
                      fe_degree,
                      fe_degree + 1,
                      double,
                      LinearAlgebra::distributed::Vector<double>>
    mf(mf_data_orig);
  mf.vmult(out_orig, in_orig);

  MatrixFreeAdvection<dim,
                      fe_degree,
                      fe_degree + 1,
                      double,
                      LinearAlgebra::distributed::Vector<double>>
    mf2(mf_data);
  mf2.vmult(out, in);

  for (unsigned int i = 0; i < out.local_size(); ++i)
    out(renumbering[i]) -= out_orig.local_element(i);

  double diff_norm = out.linfty_norm() / out_orig.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << " ";

  // test again, now doing matrix-vector product twice
  mf2.vmult(out, in);
  mf2.vmult(out, in);
  for (unsigned int i = 0; i < out.local_size(); ++i)
    out(renumbering[i]) -= out_orig.local_element(i);
  diff_norm = out.linfty_norm() / out_orig.linfty_norm();
  deallog << diff_norm << std::endl;
}
