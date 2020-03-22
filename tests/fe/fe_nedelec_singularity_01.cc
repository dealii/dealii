// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
//
// By Daniel Garcia-Sanchez, CNRS
//
// Test a maxwell singularity in 3D. Maxwell singularities are common in sharp
// metallic edges such as the Fichera corner. Here we test the elements Nedelec
// and NedelecSZ using the L2 norm and the continuity of the solution.
//
// This test solves the complex valued curl-curl equation in 3D:
//
// curl((1/mu_r)curl(E))
//   -omega^2*epsilon_0*mu_0(epsilon_r-(i*sigma/(omega*epsilon_0)))E
//                                                        = RightHandSide
//
// The manufactured solution is:
//
// E_x = (((x^2 / sqrt(x^2 + y^2)) * (x^2 - (dimension_x / 2)^2) *
//                 (y^2 - (dimension_y / 2)^2) * (z^2 - (dimension_z / 2)^2)) /
//        ((dimension_x / 2)^3 * (dimension_y / 2)^2 * (dimension_z / 2)^2))
// E_y = ( ((y^2 / sqrt(x^2 + y^2)) * (x^2 - (dimension_x / 2)^2) *
//                 (y^2 - (dimension_y / 2)^2) * (z^2 - (dimension_z / 2)^2)) /
//        ((dimension_x / 2)^2 * (dimension_y / 2)^3 * (dimension_z / 2)^2))
// E_z =  10 * (x * (x^2 - (dimension_x / 2)^2) * (y^2 - (dimension_y / 2)^2) /
//                      ((dimension_x / 2)^2 * (dimension_y / 2)^2))
//
// This manufactured solution has a singularity at x = y = 0
//
// The right hand side can be calculated with a symbolic math package such as
// sympy.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



namespace nedelec_singularity
{
  // For the sake of simplicity define the parameters as global variables.
  static const double dimension_x             = 0.04;
  static const double dimension_y             = 0.04;
  static const double dimension_z             = 0.04;
  static const double epsilon_0               = 8.85418782e-12;
  static const double mu_0                    = 1.25663706e-06;
  static const double epsilon_r               = 1;
  static const double mu_r                    = 1;
  static const double sigma                   = 0.0001;
  static const double omega                   = 6e9 * 2 * numbers::PI;
  static unsigned int nb_probe_points         = 100;
  static unsigned int grid_level              = 1;
  static unsigned int coarse_mesh_divisions_z = 3;



  template <int dim>
  class ExactSolution : public Function<dim, std::complex<double>>
  {
  public:
    ExactSolution();
    virtual std::complex<double>
    value(const Point<dim> &p, const unsigned int component) const override;
  };



  template <int dim>
  ExactSolution<dim>::ExactSolution()
    : Function<dim, std::complex<double>>(dim)
  {}



  template <int dim>
  std::complex<double>
  ExactSolution<dim>::value(const Point<dim> & p,
                            const unsigned int component) const
  {
    const double R_x = p[0];
    const double R_y = p[1];
    const double R_z = p[2];

    switch (component)
      {
        case 0:
          return 2 * std::pow(R_x, 2) *
                 (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                 (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) *
                 (4 * std::pow(R_z, 2) - std::pow(dimension_z, 2)) /
                 (std::pow(dimension_x, 3) * std::pow(dimension_y, 2) *
                  std::pow(dimension_z, 2) *
                  std::sqrt(std::pow(R_x, 2) + std::pow(R_y, 2)));
          break;
        case 1:
          return 2 * std::pow(R_x, 2) *
                 (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                 (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) *
                 (4 * std::pow(R_z, 2) - std::pow(dimension_z, 2)) /
                 (std::pow(dimension_x, 3) * std::pow(dimension_y, 2) *
                  std::pow(dimension_z, 2) *
                  std::sqrt(std::pow(R_x, 2) + std::pow(R_y, 2)));
          break;
        case 2:
          return 10 * R_x * (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                 (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) /
                 (std::pow(dimension_x, 2) * std::pow(dimension_y, 2));
          break;
        default:
          Assert(false, ExcNotImplemented());
          return 0;
      }
  }



  template <int dim>
  class RightHandSide : public Function<dim, std::complex<double>>
  {
  public:
    RightHandSide();
    virtual std::complex<double>
    value(const Point<dim> &p, const unsigned int component) const override;
  };



  template <int dim>
  RightHandSide<dim>::RightHandSide()
    : Function<dim, std::complex<double>>(dim)
  {}



  template <int dim>
  std::complex<double>
  RightHandSide<dim>::value(const Point<dim> & p,
                            const unsigned int component) const
  {
    const double R_x = p[0];
    const double R_y = p[1];
    const double R_z = p[2];

    const std::complex<double> I(0, 1);

    switch (component)
      {
        case 0:
          return 2. * R_x *
                 (-R_x * dimension_y * mu_0 * mu_r * omega *
                    std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                    (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                    (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) *
                    (4 * std::pow(R_z, 2) - std::pow(dimension_z, 2)) *
                    (epsilon_0 * epsilon_r * omega - I * sigma) -
                  8 * R_x * dimension_y *
                    std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                    (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                    (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                  (4 * std::pow(R_z, 2) - std::pow(dimension_z, 2)) *
                    (-3 * R_x * std::pow(R_y, 2) * dimension_y *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) -
                     8 * R_x * dimension_y *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) +
                     R_x * dimension_y * (std::pow(R_x, 2) + std::pow(R_y, 2)) *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                       (20 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                     3 * std::pow(R_y, 3) * dimension_x *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                     16 * R_y * dimension_x *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                       (8 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                     2 * R_y * dimension_x *
                       (std::pow(R_x, 2) + std::pow(R_y, 2)) *
                       (std::pow(R_y, 2) * (-16 * std::pow(R_x, 2) +
                                            4 * std::pow(dimension_x, 2)) +
                        std::pow(R_y, 2) * (-16 * std::pow(R_y, 2) +
                                            4 * std::pow(dimension_y, 2)) -
                        (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                          (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2))))) /
                 (std::pow(dimension_x, 3) * std::pow(dimension_y, 3) *
                  std::pow(dimension_z, 2) * mu_r *
                  std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 5.0 / 2.0));
          break;
        case 1:
          return 2. * R_y *
                 (-R_y * dimension_x * mu_0 * mu_r * omega *
                    std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                    (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                    (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) *
                    (4 * std::pow(R_z, 2) - std::pow(dimension_z, 2)) *
                    (epsilon_0 * epsilon_r * omega - I * sigma) -
                  8 * R_y * dimension_x *
                    std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                    (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                    (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                  (4 * std::pow(R_z, 2) - std::pow(dimension_z, 2)) *
                    (3 * std::pow(R_x, 3) * dimension_y *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) -
                     3 * std::pow(R_x, 2) * R_y * dimension_x *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                     16 * R_x * dimension_y *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                       (8 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) +
                     2 * R_x * dimension_y *
                       (std::pow(R_x, 2) + std::pow(R_y, 2)) *
                       (std::pow(R_x, 2) * (-16 * std::pow(R_x, 2) +
                                            4 * std::pow(dimension_x, 2)) +
                        std::pow(R_x, 2) * (-16 * std::pow(R_y, 2) +
                                            4 * std::pow(dimension_y, 2)) -
                        (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                          (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2))) -
                     8 * R_y * dimension_x *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                     R_y * dimension_x * (std::pow(R_x, 2) + std::pow(R_y, 2)) *
                       (20 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)))) /
                 (std::pow(dimension_x, 3) * std::pow(dimension_y, 3) *
                  std::pow(dimension_z, 2) * mu_r *
                  std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 5.0 / 2.0));
          break;
        case 2:
          return 2. *
                 (-5 * R_x * dimension_x * dimension_y *
                    std::pow(dimension_z, 2) * mu_0 * mu_r * omega *
                    std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) *
                    (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                    (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) *
                    (epsilon_0 * epsilon_r * omega - I * sigma) +
                  8 * R_x * dimension_y *
                    (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) *
                    (-std::pow(R_x, 2) * R_z *
                       std::sqrt(std::pow(R_x, 2) + std::pow(R_y, 2)) *
                       (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) +
                     2 * R_z *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2),
                                3.0 / 2.0) *
                       (8 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) -
                     15 * dimension_x * std::pow(dimension_z, 2) *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2)) +
                  8 * dimension_x *
                    (4 * std::pow(R_x, 2) - std::pow(dimension_x, 2)) *
                    (-5 * R_x * dimension_y * std::pow(dimension_z, 2) *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2) -
                     std::pow(R_y, 3) * R_z *
                       std::sqrt(std::pow(R_x, 2) + std::pow(R_y, 2)) *
                       (4 * std::pow(R_y, 2) - std::pow(dimension_y, 2)) +
                     2 * R_y * R_z *
                       std::pow(std::pow(R_x, 2) + std::pow(R_y, 2),
                                3.0 / 2.0) *
                       (8 * std::pow(R_y, 2) - std::pow(dimension_y, 2)))) /
                 (std::pow(dimension_x, 3) * std::pow(dimension_y, 3) *
                  std::pow(dimension_z, 2) * mu_r *
                  std::pow(std::pow(R_x, 2) + std::pow(R_y, 2), 2));
          break;
        default:
          Assert(false, ExcNotImplemented());
          return 0;
      }
  }



  template <int dim, typename FiniteElementT>
  class NedelecSingularity
  {
  public:
    NedelecSingularity();
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    output_results();

    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    const unsigned int                        fe_order;
    const QGauss<dim>                         quadrature_formula;
    FiniteElementT                            fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    AffineConstraints<std::complex<double>>   constraints;
    LinearAlgebraPETSc::MPI::SparseMatrix     system_matrix;
    LinearAlgebraPETSc::MPI::Vector           locally_relevant_solution;
    LinearAlgebraPETSc::MPI::Vector           system_rhs;
  };



  template <int dim, typename FiniteElementT>
  NedelecSingularity<dim, FiniteElementT>::NedelecSingularity()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe_order(1)
    , quadrature_formula(fe_order + 2)
    , fe(fe_order)
    , dof_handler(triangulation)
  {}

  template <int dim, typename FiniteElementT>
  void
  NedelecSingularity<dim, FiniteElementT>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);

    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    const unsigned int first_vector_component = 0;
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      first_vector_component,
      ZeroFunction<dim, std::complex<double>>(dim),
      0,
      constraints);

    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
  }



  template <int dim, typename FiniteElementT>
  void
  NedelecSingularity<dim, FiniteElementT>::assemble_system()
  {
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<std::complex<double>> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<std::complex<double>>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const RightHandSide<dim> right_hand_side;

    std::vector<Vector<std::complex<double>>> rhs_values(
      n_q_points, Vector<std::complex<double>>(dim));

    const FEValuesExtractors::Vector electric_field(0);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            cell_matrix = 0;
            cell_rhs    = 0;

            fe_values.reinit(cell);

            right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                              rhs_values);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                Tensor<1, dim, std::complex<double>> rhs;

                for (unsigned int component = 0; component < dim; ++component)
                  {
                    // Convert vectors to tensors
                    rhs[component] = rhs_values[q][component];
                  }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const Tensor<1, dim> phi_i =
                      fe_values[electric_field].value(i, q);
                    const Tensor<1, dim> curl_phi_i =
                      fe_values[electric_field].curl(i, q);

                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        const Tensor<1, dim> phi_j =
                          fe_values[electric_field].value(j, q);
                        const Tensor<1, dim> curl_phi_j =
                          fe_values[electric_field].curl(j, q);

                        std::complex<double> matrix_sum = 0;

                        matrix_sum +=
                          std::pow(omega, 2) *
                          (-epsilon_0 * mu_0 * epsilon_r * phi_i * phi_j);
                        matrix_sum += omega * std::complex<double>(0, 1) *
                                      mu_0 * sigma * phi_i * phi_j;
                        matrix_sum += (1 / mu_r) * curl_phi_i * curl_phi_j;

                        cell_matrix(i, j) += matrix_sum * fe_values.JxW(q);
                      }

                    cell_rhs(i) += phi_i * rhs * fe_values.JxW(q);
                  }
              }
            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
          }
      }
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }



  template <int dim, typename FiniteElementT>
  void
  NedelecSingularity<dim, FiniteElementT>::solve()
  {
    LinearAlgebraPETSc::MPI::Vector completely_distributed_solution(
      locally_owned_dofs, mpi_communicator);

    SolverControl                    solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);

    constraints.distribute(completely_distributed_solution);
    locally_relevant_solution = completely_distributed_solution;
  }



  template <int dim, typename FiniteElementT>
  void
  NedelecSingularity<dim, FiniteElementT>::output_results()
  {
    {
      const ExactSolution<dim> exact_solution_function;
      Vector<double> difference_per_cell(triangulation.n_active_cells());

      VectorTools::integrate_difference(dof_handler,
                                        locally_relevant_solution,
                                        exact_solution_function,
                                        difference_per_cell,
                                        QGauss<dim>(fe_order + 2),
                                        VectorTools::L2_norm);
      const double L2_error =
        VectorTools::compute_global_error(triangulation,
                                          difference_per_cell,
                                          VectorTools::L2_norm);

      deallog << " L2_error: " << L2_error << std::endl;

      // Check the continuity between between two adjacent elements. Nedelec
      // enforces the continuity only on the tangencial component. Although, if
      // the solution of the PDE is correct, the perpendicular component should
      // also be continuous. An element boundary can be found at
      // x = dimension_x/3.
      const double     delta = dimension_x / 1000.;
      const Point<dim> point_a(dimension_x / 3. - delta, delta, delta);
      const Point<dim> point_b(dimension_x / 3. + delta, delta, delta);
      deallog << " Point_a = " << point_a << std::endl;
      deallog << " Point_b = " << point_b << std::endl;
      Vector<std::complex<double>> solution_a(dim);
      Vector<std::complex<double>> solution_b(dim);
      solution_a = 0;
      solution_b = 0;
      {
        bool point_in_locally_owned_cell;
        auto mapping = StaticMappingQ1<dim>::mapping;
        // find the cell in which this point
        // is, initialize a quadrature rule with
        // it, and then a FEValues object
        const std::pair<typename DoFHandler<dim>::active_cell_iterator,
                        Point<dim>>
          cell_point = GridTools::find_active_cell_around_point(mapping,
                                                                dof_handler,
                                                                point_a);

        point_in_locally_owned_cell = cell_point.first->is_locally_owned();
        if (point_in_locally_owned_cell)
          {
            VectorTools::point_value(dof_handler,
                                     locally_relevant_solution,
                                     point_a,
                                     solution_a);
          }
      }
      {
        bool point_in_locally_owned_cell;
        auto mapping = StaticMappingQ1<dim>::mapping;
        // find the cell in which this point
        // is, initialize a quadrature rule with
        // it, and then a FEValues object
        const std::pair<typename DoFHandler<dim>::active_cell_iterator,
                        Point<dim>>
          cell_point = GridTools::find_active_cell_around_point(mapping,
                                                                dof_handler,
                                                                point_b);

        point_in_locally_owned_cell = cell_point.first->is_locally_owned();
        if (point_in_locally_owned_cell)
          {
            VectorTools::point_value(dof_handler,
                                     locally_relevant_solution,
                                     point_b,
                                     solution_b);
          }
      }
      // Only one process has the solution_a or/and solution_b. This is a simple
      // approach to send solution_a and solution_b to all the processes.
      Utilities::MPI::sum(solution_a, mpi_communicator, solution_a);
      Utilities::MPI::sum(solution_b, mpi_communicator, solution_b);
      deallog << " Solution(point_a) : " << solution_a << std::endl;
      deallog << " Solution(point_b) : " << solution_b << std::endl;
      // Vector does not provide operator-
      deallog << " Solution(point_b) - solution (point_a): "
              << (solution_b -= solution_a) << std::endl;
    }

    {
      std::vector<std::string> solution_names(1, "electric_field_x");
      if (dim >= 2)
        {
          solution_names.emplace_back("electric_field_y");
        }
      if (dim == 3)
        {
          solution_names.emplace_back("electric_field_z");
        }
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(dim, DataComponentInterpretation::component_is_scalar);

      DataOut<dim> data_out;
      data_out.add_data_vector(dof_handler,
                               locally_relevant_solution,
                               solution_names,
                               interpretation);
      Vector<float> subdomain(triangulation.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      const RightHandSide<dim>                  rhs_function;
      const ExactSolution<dim>                  exact_solution_function;
      std::vector<Vector<std::complex<double>>> rhs(
        dim, Vector<std::complex<double>>(triangulation.n_active_cells()));
      std::vector<Vector<std::complex<double>>> exact_solution(
        dim, Vector<std::complex<double>>(triangulation.n_active_cells()));

      // Loop over all the cells
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
                {
                  rhs[dim_idx](cell->active_cell_index()) =
                    rhs_function.value(cell->center(), dim_idx);
                  exact_solution[dim_idx](cell->active_cell_index()) =
                    exact_solution_function.value(cell->center(), dim_idx);
                }
            }
          // And on the cells that we are not interested in, set the respective
          // value in the vector to a random value in order to make sure that if
          // we were somehow wrong about our assumption that these elements
          // would not appear in the output file, that we would find out by
          // looking at the graphical output:
          else
            {
              for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
                {
                  rhs[dim_idx](cell->active_cell_index()) = -1e90;
                }
            }
        }

      for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
        {
          data_out.add_data_vector(rhs[dim_idx],
                                   "rhs_" + std::to_string(dim_idx));
          data_out.add_data_vector(exact_solution[dim_idx],
                                   "exact_solution_" + std::to_string(dim_idx));
        }

      data_out.build_patches(2);

      unsigned int nb_number_positions;
      if (std::is_same<FiniteElementT, FE_Nedelec<dim>>::value)
        {
          data_out.write_vtu_in_parallel("result_nedelec.vtu",
                                         mpi_communicator);
        }
      else if (std::is_same<FiniteElementT, FE_NedelecSZ<dim>>::value)
        {
          data_out.write_vtu_in_parallel("result_nedelec_sz.vtu",
                                         mpi_communicator);
        }
      else
        {
          Assert(false, ExcInternalError());
        }
    }
  }



  template <int dim, typename FiniteElementT>
  void
  NedelecSingularity<dim, FiniteElementT>::run()
  {
    {
      Point<dim> p0;
      p0(0) = -dimension_x / 2;
      p0(1) = -dimension_y / 2;
      p0(2) = -dimension_z / 2;
      Point<dim> p1;
      p1(0) = dimension_x / 2;
      p1(1) = dimension_y / 2;
      p1(2) = dimension_z / 2;
      double smallest_dimension =
        std::min(dimension_z, std::min(dimension_x, dimension_y));
      std::vector<unsigned int> divisions(dim);
      divisions[0] = std::max<unsigned int>(coarse_mesh_divisions_z, 1) *
                     int((p1(0) - p0(0)) / smallest_dimension);
      divisions[1] = std::max<unsigned int>(coarse_mesh_divisions_z, 1) *
                     int((p1(1) - p0(1)) / smallest_dimension);
      divisions[2] = std::max<unsigned int>(coarse_mesh_divisions_z, 1) *
                     int((p1(2) - p0(2)) / smallest_dimension);
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                divisions,
                                                p0,
                                                p1);
    }

    if (grid_level > 0)
      {
        triangulation.refine_global(grid_level);
      }

    setup_system();
    deallog << " Number of active cells :       "
            << triangulation.n_active_cells() << std::endl;
    deallog << " Number of degrees of freedom : " << dof_handler.n_dofs()
            << std::endl;


    assemble_system();
    solve();

    output_results();
  }
} // namespace nedelec_singularity

int
main(int argc, char *argv[])
{
  try
    {
      const int dim = 3;

      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      MPILogInitAll log;

      {
        nedelec_singularity::NedelecSingularity<dim, FE_Nedelec<dim>>
          nedelec_singularity_3d;
        nedelec_singularity_3d.run();
      }

      {
        nedelec_singularity::NedelecSingularity<dim, FE_NedelecSZ<dim>>
          nedelec_singularity_3d;
        nedelec_singularity_3d.run();
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
