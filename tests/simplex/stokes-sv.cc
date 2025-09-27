// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Stokes on a simplex mesh using barycenter refinement to show
// that Scott-Vogelius elements are pointwise divergence free.

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include "deal.II/numerics/vector_tools_mean_value.h"
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

// #define HEX


namespace Step56
{
  using namespace dealii;

#ifdef HEX
  template <int dim>
  using QuadratureT = QGauss<dim>;
#else
  template <int dim>
  using QuadratureT = QGaussSimplex<dim>;
#endif

  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution()
      : Function<dim>(dim + 1)
    {}
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;
  };

  template <>
  double
  Solution<2>::value(const Point<2> &p, const unsigned int component) const
  {
    Assert(component <= 2 + 1, ExcIndexRange(component, 0, 2 + 1));

    using numbers::PI;
    const double x = p(0);
    const double y = p(1);

    if (component == 0)
      return sin(PI * x);
    if (component == 1)
      return -PI * y * cos(PI * x);
    if (component == 2)
      return sin(PI * x) * cos(PI * y);

    return 0;
  }

  template <>
  double
  Solution<3>::value(const Point<3> &p, const unsigned int component) const
  {
    Assert(component <= 3 + 1, ExcIndexRange(component, 0, 3 + 1));

    using numbers::PI;
    const double x = p(0);
    const double y = p(1);
    const double z = p(2);

    if (component == 0)
      return 2.0 * sin(PI * x);
    if (component == 1)
      return -PI * y * cos(PI * x);
    if (component == 2)
      return -PI * z * cos(PI * x);
    if (component == 3)
      return sin(PI * x) * cos(PI * y) * sin(PI * z);

    return 0;
  }

  template <>
  Tensor<1, 2>
  Solution<2>::gradient(const Point<2> &p, const unsigned int component) const
  {
    Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1));

    using numbers::PI;
    const double x = p(0);
    const double y = p(1);

    Tensor<1, 2> return_value;
    if (component == 0)
      {
        return_value[0] = PI * cos(PI * x);
        return_value[1] = 0.0;
      }
    else if (component == 1)
      {
        return_value[0] = y * PI * PI * sin(PI * x);
        return_value[1] = -PI * cos(PI * x);
      }
    else if (component == 2)
      {
        return_value[0] = PI * cos(PI * x) * cos(PI * y);
        return_value[1] = -PI * sin(PI * x) * sin(PI * y);
      }

    return return_value;
  }

  template <>
  Tensor<1, 3>
  Solution<3>::gradient(const Point<3> &p, const unsigned int component) const
  {
    Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1));

    using numbers::PI;
    const double x = p(0);
    const double y = p(1);
    const double z = p(2);

    Tensor<1, 3> return_value;
    if (component == 0)
      {
        return_value[0] = 2 * PI * cos(PI * x);
        return_value[1] = 0.0;
        return_value[2] = 0.0;
      }
    else if (component == 1)
      {
        return_value[0] = y * PI * PI * sin(PI * x);
        return_value[1] = -PI * cos(PI * x);
        return_value[2] = 0.0;
      }
    else if (component == 2)
      {
        return_value[0] = z * PI * PI * sin(PI * x);
        return_value[1] = 0.0;
        return_value[2] = -PI * cos(PI * x);
      }
    else if (component == 3)
      {
        return_value[0] = PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
        return_value[1] = -PI * sin(PI * x) * sin(PI * y) * sin(PI * z);
        return_value[2] = PI * sin(PI * x) * cos(PI * y) * cos(PI * z);
      }

    return return_value;
  }

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim + 1)
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };

  template <>
  double
  RightHandSide<2>::value(const Point<2> &p, const unsigned int component) const
  {
    Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1));

    using numbers::PI;
    double x = p(0);
    double y = p(1);
    if (component == 0)
      return PI * PI * sin(PI * x) + PI * cos(PI * x) * cos(PI * y);
    if (component == 1)
      return -PI * PI * PI * y * cos(PI * x) - PI * sin(PI * y) * sin(PI * x);
    if (component == 2)
      return 0;

    return 0;
  }

  template <>
  double
  RightHandSide<3>::value(const Point<3> &p, const unsigned int component) const
  {
    Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1));

    using numbers::PI;
    double x = p(0);
    double y = p(1);
    double z = p(2);
    if (component == 0)
      return 2 * PI * PI * sin(PI * x) +
             PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
    if (component == 1)
      return -PI * PI * PI * y * cos(PI * x) +
             PI * (-1) * sin(PI * y) * sin(PI * x) * sin(PI * z);
    if (component == 2)
      return -PI * PI * PI * z * cos(PI * x) +
             PI * cos(PI * z) * sin(PI * x) * cos(PI * y);
    if (component == 3)
      return 0;

    return 0;
  }

  template <class PreconditionerAType, class PreconditionerSType>
  class BlockSchurPreconditioner : public EnableObserverPointer
  {
  public:
    BlockSchurPreconditioner(
      const BlockSparseMatrix<double> &system_matrix,
      const SparseMatrix<double>      &schur_complement_matrix,
      const PreconditionerAType       &preconditioner_A,
      const PreconditionerSType       &preconditioner_S);

    void
    vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;

    mutable unsigned int n_iterations_A;
    mutable unsigned int n_iterations_S;

  private:
    const BlockSparseMatrix<double> &system_matrix;
    const SparseMatrix<double>      &schur_complement_matrix;
    const PreconditionerAType       &preconditioner_A;
    const PreconditionerSType       &preconditioner_S;
  };

  template <class PreconditionerAType, class PreconditionerSType>
  BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::
    BlockSchurPreconditioner(
      const BlockSparseMatrix<double> &system_matrix,
      const SparseMatrix<double>      &schur_complement_matrix,
      const PreconditionerAType       &preconditioner_A,
      const PreconditionerSType       &preconditioner_S)
    : n_iterations_A(0)
    , n_iterations_S(0)
    , system_matrix(system_matrix)
    , schur_complement_matrix(schur_complement_matrix)
    , preconditioner_A(preconditioner_A)
    , preconditioner_S(preconditioner_S)
  {}



  template <class PreconditionerAType, class PreconditionerSType>
  void
  BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::vmult(
    BlockVector<double>       &dst,
    const BlockVector<double> &src) const
  {
    Vector<double> utmp(src.block(0));

    {
      n_iterations_S += 1;
      preconditioner_S.vmult(dst.block(1), src.block(1));
      dst.block(1) *= -1.0;
    }

    {
      system_matrix.block(0, 1).vmult(utmp, dst.block(1));
      utmp *= -1.0;
      utmp += src.block(0);
    }

    {
      preconditioner_A.vmult(dst.block(0), utmp);
      n_iterations_A += 1;
    }
  }

  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(const unsigned int pressure_degree);
    void
    run();

  private:
    void
    setup_dofs();
    void
    assemble_system();
    void
    solve();
    void
    compute_errors();
    void
    output_results(const unsigned int refinement_cycle) const;

    const unsigned int pressure_degree;

    Triangulation<dim> triangulation;

#ifdef HEX
    MappingQ1<dim> mapping;
#else
    MappingFE<dim> mapping;
#endif

    FESystem<dim>   velocity_fe;
    FESystem<dim>   fe;
    DoFHandler<dim> dof_handler;
    DoFHandler<dim> velocity_dof_handler;

    AffineConstraints<double> constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    SparseMatrix<double>      pressure_mass_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem(const unsigned int pressure_degree)

    : pressure_degree(pressure_degree)
#ifdef HEX
    , velocity_fe(FE_Q<dim>(pressure_degree + 1), dim)
    , fe(velocity_fe, 1, FE_Q<dim>(pressure_degree), 1)
#else
    , mapping(FE_SimplexP<dim>(1))
    , velocity_fe(FE_SimplexP<dim>(pressure_degree + 1), dim)
    , fe(velocity_fe, 1, FE_SimplexDGP<dim>(pressure_degree), 1)
#endif
    , dof_handler(triangulation)
    , velocity_dof_handler(triangulation)
  {}


  template <int dim>
  void
  StokesProblem<dim>::setup_dofs()
  {
    system_matrix.clear();
    pressure_mass_matrix.clear();

    dof_handler.distribute_dofs(fe);

    std::vector<unsigned int> block_component(2);
    block_component[0] = 0;
    block_component[1] = 1;

    const FEValuesExtractors::Vector velocities(0);

    DoFRenumbering::block_wise(dof_handler);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    {
      constraints.clear();
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Solution<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));

      constraints.close();
    }

    deallog << "\tNumber of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "\tNumber of degrees of freedom: " << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')' << std::endl;

    {
      BlockDynamicSparsityPattern csp(dofs_per_block, dofs_per_block);
      DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
      sparsity_pattern.copy_from(csp);
    }
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dofs_per_block);
    system_rhs.reinit(dofs_per_block);
  }

  template <int dim>
  void
  StokesProblem<dim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const bool assemble_pressure_mass_matrix = true;

    const QuadratureT<dim> quadrature_formula(pressure_degree + 2);

    FEValues<dim> fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const RightHandSide<dim>    right_hand_side;
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1));

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
    std::vector<double>                  div_phi_u(dofs_per_cell);
    std::vector<double>                  phi_p(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        local_matrix = 0;
        local_rhs    = 0;

        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                symgrad_phi_u[k] =
                  fe_values[velocities].symmetric_gradient(k, q);
                div_phi_u[k] = fe_values[velocities].divergence(k, q);
                phi_p[k]     = fe_values[pressure].value(k, q);
              }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j <= i; ++j)
                  {
                    local_matrix(i, j) +=
                      (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                       (assemble_pressure_mass_matrix ? phi_p[i] * phi_p[j] :
                                                        0)) *
                      fe_values.JxW(q);
                  }

                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                local_rhs(i) += fe_values.shape_value(i, q) *
                                rhs_values[q](component_i) * fe_values.JxW(q);
              }
          }

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            local_matrix(i, j) = local_matrix(j, i);

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);
      }

    {
      pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
      pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
      system_matrix.block(1, 1) = 0;
    }
  }

  template <int dim>
  void
  StokesProblem<dim>::solve()
  {
    constraints.set_zero(solution);

    SolverControl solver_control(10000,
                                 1e-10 * system_rhs.l2_norm(),
                                 false,
                                 false);
    unsigned int  n_iterations_A;
    unsigned int  n_iterations_S;

    SolverGMRES<BlockVector<double>> solver(
      solver_control,
      SolverGMRES<BlockVector<double>>::AdditionalData(
        50, true)); // right preconditioning

    {
      SparseILU<double> A_preconditioner;
      A_preconditioner.initialize(system_matrix.block(0, 0));

      SparseILU<double> S_preconditioner;
      S_preconditioner.initialize(pressure_mass_matrix);

      const BlockSchurPreconditioner<SparseILU<double>, SparseILU<double>>
        preconditioner(system_matrix,
                       pressure_mass_matrix,
                       A_preconditioner,
                       S_preconditioner);

      {
        solver.solve(system_matrix, solution, system_rhs, preconditioner);
        n_iterations_A = preconditioner.n_iterations_A;
        n_iterations_S = preconditioner.n_iterations_S;
      }
    }

    constraints.distribute(solution);
  }

  template <int dim>
  void
  StokesProblem<dim>::compute_errors()
  {
    const double mean_pressure =
      VectorTools::compute_mean_value(mapping,
                                      dof_handler,
                                      QuadratureT<dim>(pressure_degree + 2),
                                      solution,
                                      dim);
    VectorTools::add_constant(solution, dof_handler, dim, -mean_pressure);

    const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim + 1);

    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QuadratureT<dim>(pressure_degree + 2),
                                      VectorTools::L2_norm,
                                      &velocity_mask);

    const double Velocity_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QuadratureT<dim>(pressure_degree + 2),
                                      VectorTools::L2_norm,
                                      &pressure_mask);

    const double Pressure_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QuadratureT<dim>(pressure_degree + 2),
                                      VectorTools::H1_norm,
                                      &velocity_mask);

    const double Velocity_H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_norm);

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QuadratureT<dim>(pressure_degree + 2),
                                      VectorTools::Hdiv_seminorm,
                                      &velocity_mask);

    const double Velocity_Hdiv_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::Hdiv_seminorm);
    deallog << std::endl
            << "   Velocity L2 Error: " << Velocity_L2_error << std::endl
            << "   Pressure L2 Error: " << Pressure_L2_error << std::endl
            << "   Velocity H1 Error: " << Velocity_H1_error << std::endl
            << "   Velocity Hdiv Err: " << Velocity_Hdiv_error << std::endl;
  }

  template <int dim>
  void
  StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const
  {
    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.emplace_back("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    data_out.build_patches();

    std::ofstream output("solution-" +
                         Utilities::int_to_string(refinement_cycle, 2) +
                         "-dim-" + Utilities::int_to_string(dim, 2) + ".vtk");
    data_out.write_vtk(output);
  }

  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    // For Scott-Vogelius elements in 3D, we must use a cubic velocity space and
    // so we must use fewer refinements.
    const unsigned int max_refinements_cycles = dim < 3 ? 3 : 2;
    for (unsigned int refinement_cycle = 0;
         refinement_cycle < max_refinements_cycles;
         ++refinement_cycle)
      {
        deallog << "Refinement cycle " << refinement_cycle << std::endl;

        Triangulation<dim> s_tria;
        GridGenerator::subdivided_hyper_cube_with_simplices<dim, dim>(
          s_tria, std::pow(1 + refinement_cycle, 2));
        triangulation.clear();
        GridGenerator::alfeld_split_of_simplex_mesh(s_tria, triangulation);

        deallog << "   Set-up..." << std::endl;
        setup_dofs();

        deallog << "   Assembling..." << std::endl;
        assemble_system();

        deallog << "   Solving..." << std::flush;
        solve();

        compute_errors();

        output_results(refinement_cycle);
      }
  }
} // namespace Step56

int
main()
{
  initlog();
  {
    using namespace Step56;

    const int          degree = 1;
    const int          dim    = 2;
    StokesProblem<dim> flow_problem(degree);

    flow_problem.run();
  }
  {
    using namespace Step56;

    const int          degree = 2; // In 3D, not stable unless degree >= 2.
    const int          dim    = 3;
    StokesProblem<dim> flow_problem(degree);

    flow_problem.run();
  }
  return 0;
}
