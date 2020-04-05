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



// compares the computation of the diagonal using a face integration
// facilities with the alternative reinit(cell_index, face_number) approach

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

#include "create_mesh.h"

template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename number   = double>
class LaplaceOperator : public Subscriptor
{
public:
  typedef number value_type;

  LaplaceOperator(){};

  void
  initialize(const Mapping<dim> &   mapping,
             const DoFHandler<dim> &dof_handler,
             const unsigned int     level = numbers::invalid_unsigned_int)
  {
    const QGauss<1>                                  quad(n_q_points_1d);
    typename MatrixFree<dim, number>::AdditionalData addit_data;
    addit_data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::none;
    addit_data.tasks_block_size = 3;
    addit_data.mg_level         = level;
    addit_data.mapping_update_flags_inner_faces =
      update_JxW_values | update_normal_vectors | update_jacobians;
    addit_data.mapping_update_flags_boundary_faces =
      update_JxW_values | update_normal_vectors | update_jacobians;
    addit_data.mapping_update_flags_faces_by_cells =
      update_JxW_values | update_normal_vectors | update_jacobians;
    AffineConstraints<double> constraints;
    constraints.close();

    data.reinit(mapping, dof_handler, constraints, quad, addit_data);
  }

  void
  compute_diagonal_by_face(
    LinearAlgebra::distributed::Vector<number> &result) const
  {
    int dummy;
    result = 0;
    data.loop(&LaplaceOperator::local_diagonal_dummy,
              &LaplaceOperator::local_diagonal_face,
              &LaplaceOperator::local_diagonal_boundary,
              this,
              result,
              dummy);
  }

  void
  compute_diagonal_by_cell(
    LinearAlgebra::distributed::Vector<number> &result) const
  {
    int dummy;
    result.zero_out_ghosts();
    data.cell_loop(&LaplaceOperator::local_diagonal_by_cell,
                   this,
                   result,
                   dummy);
  }

  void
  initialize_dof_vector(
    LinearAlgebra::distributed::Vector<number> &vector) const
  {
    data.initialize_dof_vector(vector);
  }


private:
  void
  local_diagonal_dummy(const MatrixFree<dim, number> &,
                       LinearAlgebra::distributed::Vector<number> &,
                       const int &,
                       const std::pair<unsigned int, unsigned int> &) const
  {}

  void
  local_diagonal_face(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const int &,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data, true);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi_outer(data,
                                                                         false);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        phi.reinit(face);
        phi_outer.reinit(face);

        VectorizedArray<number> local_diagonal_vector[phi.static_dofs_per_cell];

        // Choose a one-sided definition of sigmaF rather than the actual
        // penalty parameter because the cell-based method cannot read the
        // sigma parameter on the neighbor
        VectorizedArray<number> sigmaF =
          std::abs(
            (phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]) *
          (number)(std::max(fe_degree, 1) * (fe_degree + 1.0)) * 2.0;

        // Compute phi part
        for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
          phi_outer.begin_dof_values()[j] = VectorizedArray<number>();
        phi_outer.evaluate(true, true);
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(true, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                VectorizedArray<number> average_value =
                  (phi.get_value(q) - phi_outer.get_value(q)) * 0.5;
                VectorizedArray<number> average_valgrad =
                  phi.get_normal_derivative(q) +
                  phi_outer.get_normal_derivative(q);
                average_valgrad =
                  average_value * 2. * sigmaF - average_valgrad * 0.5;
                phi.submit_normal_derivative(-average_value, q);
                phi.submit_value(average_valgrad, q);
              }
            phi.integrate(true, true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);

        // Compute phi_outer part
        sigmaF = std::abs((phi.get_normal_vector(0) *
                           phi_outer.inverse_jacobian(0))[dim - 1]) *
                 (number)(std::max(fe_degree, 1) * (fe_degree + 1.0)) * 2.0;
        for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
          phi.begin_dof_values()[j] = VectorizedArray<number>();
        phi.evaluate(true, true);
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi_outer.begin_dof_values()[j] = VectorizedArray<number>();
            phi_outer.begin_dof_values()[i] = 1.;
            phi_outer.evaluate(true, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                VectorizedArray<number> average_value =
                  (phi.get_value(q) - phi_outer.get_value(q)) * 0.5;
                VectorizedArray<number> average_valgrad =
                  phi.get_normal_derivative(q) +
                  phi_outer.get_normal_derivative(q);
                average_valgrad =
                  average_value * 2. * sigmaF - average_valgrad * 0.5;
                phi_outer.submit_normal_derivative(-average_value, q);
                phi_outer.submit_value(-average_valgrad, q);
              }
            phi_outer.integrate(true, true);
            local_diagonal_vector[i] = phi_outer.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi_outer.begin_dof_values()[i] = local_diagonal_vector[i];
        phi_outer.distribute_local_to_global(dst);
      }
  }

  void
  local_diagonal_boundary(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const int &,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        phi.reinit(face);

        VectorizedArray<number> local_diagonal_vector[phi.static_dofs_per_cell];
        VectorizedArray<number> sigmaF =
          std::abs(
            (phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]) *
          (number)(std::max(1, fe_degree) * (fe_degree + 1.0)) * 2.;

        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(true, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                VectorizedArray<number> average_value = phi.get_value(q);
                VectorizedArray<number> average_valgrad =
                  -phi.get_normal_derivative(q);
                average_valgrad += average_value * sigmaF * 2.0;
                phi.submit_normal_derivative(-average_value, q);
                phi.submit_value(average_valgrad, q);
              }
            phi.integrate(true, true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_diagonal_by_cell(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phif(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        VectorizedArray<number>
          local_diagonal_vector[phif.static_dofs_per_cell];
        for (unsigned int i = 0; i < phif.static_dofs_per_cell; ++i)
          local_diagonal_vector[i] = 0.;
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            phif.reinit(cell, face);
            VectorizedArray<number> sigmaF =
              std::abs((phif.get_normal_vector(0) *
                        phif.inverse_jacobian(0))[dim - 1]) *
              (number)(std::max(1, fe_degree) * (fe_degree + 1.0)) * 2.;

            std::array<types::boundary_id, VectorizedArray<number>::size()>
                                    boundary_ids = data.get_faces_by_cells_boundary_id(cell, face);
            VectorizedArray<number> factor_boundary;
            for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
              // interior face
              if (boundary_ids[v] == numbers::invalid_boundary_id)
                factor_boundary[v] = 0.5;
              // Dirichlet boundary
              else
                factor_boundary[v] = 1.0;
            for (unsigned int i = 0; i < phif.dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < phif.dofs_per_cell; ++j)
                  phif.begin_dof_values()[j] = VectorizedArray<number>();
                phif.begin_dof_values()[i] = 1.;
                phif.evaluate(true, true);
                for (unsigned int q = 0; q < phif.n_q_points; ++q)
                  {
                    VectorizedArray<number> average_value =
                      phif.get_value(q) * factor_boundary;
                    VectorizedArray<number> average_valgrad =
                      phif.get_normal_derivative(q) * factor_boundary;
                    average_valgrad =
                      average_value * 2. * sigmaF - average_valgrad;
                    phif.submit_normal_derivative(-average_value, q);
                    phif.submit_value(average_valgrad, q);
                  }
                phif.integrate(true, true);
                local_diagonal_vector[i] += phif.begin_dof_values()[i];
              }
          }
        for (unsigned int i = 0; i < phif.static_dofs_per_cell; ++i)
          phif.begin_dof_values()[i] = local_diagonal_vector[i];
        phif.set_dof_values(dst);
      }
  }

  MatrixFree<dim, number> data;
};



template <int dim, int fe_degree, int n_q_points_1d, typename number>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  create_mesh(tria);
  tria.refine_global(std::max(8 - fe_degree - 2 * dim, 1));

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();
  deallog << "Number of DoFs: " << dof.n_dofs() << std::endl;

  MappingQGeneric<dim>                                   mapping(fe_degree + 1);
  LaplaceOperator<dim, fe_degree, n_q_points_1d, number> fine_matrix;
  fine_matrix.initialize(mapping, dof);

  LinearAlgebra::distributed::Vector<number> res1, res2;
  fine_matrix.initialize_dof_vector(res1);
  fine_matrix.initialize_dof_vector(res2);

  fine_matrix.compute_diagonal_by_face(res1);
  fine_matrix.compute_diagonal_by_cell(res2);
  res2 -= res1;
  deallog << "Error in diagonal on active cells: " << (double)res2.linfty_norm()
          << std::endl;

  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      LaplaceOperator<dim, fe_degree, n_q_points_1d, number> fine_matrix;
      fine_matrix.initialize(mapping, dof, level);

      LinearAlgebra::distributed::Vector<number> res1, res2;
      fine_matrix.initialize_dof_vector(res1);
      fine_matrix.initialize_dof_vector(res2);

      fine_matrix.compute_diagonal_by_face(res1);
      fine_matrix.compute_diagonal_by_cell(res2);
      res2 -= res1;
      deallog << "Error in diagonal on level " << level << ":      "
              << (double)res2.linfty_norm() << std::endl;
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  mpi_initlog();

  {
    deallog.push("2d");
    test<2, 1, 2, double>();
    test<2, 2, 3, double>();
    test<2, 1, 2, float>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1, 2, double>();
    test<3, 2, 3, double>();
    test<3, 1, 2, float>();
    deallog.pop();
  }
}
