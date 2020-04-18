/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2016 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 * Author: Ryan Grove, Clemson University
 *         Timo Heister, Clemson University
 */

// @sect3{Include files}

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/grid/grid_out.h>

// We need to include the following file to do timings:
#include <deal.II/base/timer.h>

// This includes the files necessary for us to use geometric Multigrid
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <iostream>
#include <fstream>

namespace Step56
{
  using namespace dealii;

  // In order to make it easy to switch between the different solvers that are
  // being used, we declare an enum that can be passed as an argument to the
  // constructor of the main class.
  enum class SolverType
  {
    FGMRES_ILU,
    FGMRES_GMG,
    UMFPACK
  };

  // @sect3{Functions for Solution and Righthand side}
  //
  // The class Solution is used to define the boundary conditions and to
  // compute errors of the numerical solution. Note that we need to define the
  // values and gradients in order to compute L2 and H1 errors. Here we
  // decided to separate the implementations for 2d and 3d using template
  // specialization.
  //
  // Note that the first dim components are the velocity components
  // and the last is the pressure.
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution()
      : Function<dim>(dim + 1)
    {}
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };

  template <>
  double Solution<2>::value(const Point<2> &   p,
                            const unsigned int component) const
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
  double Solution<3>::value(const Point<3> &   p,
                            const unsigned int component) const
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

  // Note that for the gradient we need to return a Tensor<1,dim>
  template <>
  Tensor<1, 2> Solution<2>::gradient(const Point<2> &   p,
                                     const unsigned int component) const
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
  Tensor<1, 3> Solution<3>::gradient(const Point<3> &   p,
                                     const unsigned int component) const
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

  // Implementation of $f$. See the introduction for more information.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim + 1)
    {}

    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };

  template <>
  double RightHandSide<2>::value(const Point<2> &   p,
                                 const unsigned int component) const
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
  double RightHandSide<3>::value(const Point<3> &   p,
                                 const unsigned int component) const
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



  // @sect3{ASPECT BlockSchurPreconditioner}

  // In the following, we will implement a preconditioner that expands
  // on the ideas discussed in the Results section of step-22.
  // Specifically, we
  // 1. use an upper block-triangular preconditioner because we want to
  // use right preconditioning.
  // 2. optionally allow using an inner solver for the velocity block instead
  // of a single preconditioner application.
  // 3. do not use InverseMatrix but explicitly call SolverCG.
  // This approach is also used in the ASPECT code
  // (see https://aspect.geodynamics.org) that solves the Stokes equations in
  // the context of simulating convection in the earth mantle, and which
  // has been used to solve problems on many thousands of processors.
  //
  // The bool flag @p do_solve_A in the constructor allows us to either
  // apply the preconditioner for the velocity block once or use an inner
  // iterative solver for a more accurate approximation instead.
  //
  // Notice how we keep track of the sum of the inner iterations
  // (preconditioner applications).
  template <class PreconditionerAType, class PreconditionerSType>
  class BlockSchurPreconditioner : public Subscriptor
  {
  public:
    BlockSchurPreconditioner(
      const BlockSparseMatrix<double> &system_matrix,
      const SparseMatrix<double> &     schur_complement_matrix,
      const PreconditionerAType &      preconditioner_A,
      const PreconditionerSType &      preconditioner_S,
      const bool                       do_solve_A);

    void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;

    mutable unsigned int n_iterations_A;
    mutable unsigned int n_iterations_S;

  private:
    const BlockSparseMatrix<double> &system_matrix;
    const SparseMatrix<double> &     schur_complement_matrix;
    const PreconditionerAType &      preconditioner_A;
    const PreconditionerSType &      preconditioner_S;

    const bool do_solve_A;
  };

  template <class PreconditionerAType, class PreconditionerSType>
  BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::
    BlockSchurPreconditioner(
      const BlockSparseMatrix<double> &system_matrix,
      const SparseMatrix<double> &     schur_complement_matrix,
      const PreconditionerAType &      preconditioner_A,
      const PreconditionerSType &      preconditioner_S,
      const bool                       do_solve_A)
    : n_iterations_A(0)
    , n_iterations_S(0)
    , system_matrix(system_matrix)
    , schur_complement_matrix(schur_complement_matrix)
    , preconditioner_A(preconditioner_A)
    , preconditioner_S(preconditioner_S)
    , do_solve_A(do_solve_A)
  {}



  template <class PreconditionerAType, class PreconditionerSType>
  void
  BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::vmult(
    BlockVector<double> &      dst,
    const BlockVector<double> &src) const
  {
    Vector<double> utmp(src.block(0));

    // First solve with the approximation for S
    {
      SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
      SolverCG<Vector<double>> cg(solver_control);

      dst.block(1) = 0.0;
      cg.solve(schur_complement_matrix,
               dst.block(1),
               src.block(1),
               preconditioner_S);

      n_iterations_S += solver_control.last_step();
      dst.block(1) *= -1.0;
    }

    // Second, apply the top right block (B^T)
    {
      system_matrix.block(0, 1).vmult(utmp, dst.block(1));
      utmp *= -1.0;
      utmp += src.block(0);
    }

    // Finally, either solve with the top left block
    // or just apply one preconditioner sweep
    if (do_solve_A == true)
      {
        SolverControl            solver_control(10000, utmp.l2_norm() * 1e-4);
        SolverCG<Vector<double>> cg(solver_control);

        dst.block(0) = 0.0;
        cg.solve(system_matrix.block(0, 0),
                 dst.block(0),
                 utmp,
                 preconditioner_A);

        n_iterations_A += solver_control.last_step();
      }
    else
      {
        preconditioner_A.vmult(dst.block(0), utmp);
        n_iterations_A += 1;
      }
  }

  // @sect3{The StokesProblem class}
  //
  // This is the main class of the problem.
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(const unsigned int pressure_degree,
                  const SolverType   solver_type);
    void run();

  private:
    void setup_dofs();
    void assemble_system();
    void assemble_multigrid();
    void solve();
    void compute_errors();
    void output_results(const unsigned int refinement_cycle) const;

    const unsigned int pressure_degree;
    const SolverType   solver_type;

    Triangulation<dim> triangulation;
    FESystem<dim>      velocity_fe;
    FESystem<dim>      fe;
    DoFHandler<dim>    dof_handler;
    DoFHandler<dim>    velocity_dof_handler;

    AffineConstraints<double> constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    SparseMatrix<double>      pressure_mass_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

    MGLevelObject<SparsityPattern>      mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double>> mg_matrices;
    MGLevelObject<SparseMatrix<double>> mg_interface_matrices;
    MGConstrainedDoFs                   mg_constrained_dofs;

    TimerOutput computing_timer;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem(const unsigned int pressure_degree,
                                    const SolverType   solver_type)

    : pressure_degree(pressure_degree)
    , solver_type(solver_type)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    ,
    // Finite element for the velocity only:
    velocity_fe(FE_Q<dim>(pressure_degree + 1), dim)
    ,
    // Finite element for the whole system:
    fe(velocity_fe, 1, FE_Q<dim>(pressure_degree), 1)
    , dof_handler(triangulation)
    , velocity_dof_handler(triangulation)
    , computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
  {}



  // @sect4{StokesProblem::setup_dofs}

  // This function sets up the DoFHandler, matrices, vectors, and Multigrid
  // structures (if needed).
  template <int dim>
  void StokesProblem<dim>::setup_dofs()
  {
    TimerOutput::Scope scope(computing_timer, "Setup");

    system_matrix.clear();
    pressure_mass_matrix.clear();

    // The main DoFHandler only needs active DoFs, so we are not calling
    // distribute_mg_dofs() here
    dof_handler.distribute_dofs(fe);

    // This block structure separates the dim velocity components from
    // the pressure component (used for reordering). Note that we have
    // 2 instead of dim+1 blocks like in step-22, because our FESystem
    // is nested and the dim velocity components appear as one block.
    std::vector<unsigned int> block_component(2);
    block_component[0] = 0;
    block_component[1] = 1;

    // Velocities start at component 0:
    const FEValuesExtractors::Vector velocities(0);

    // ILU behaves better if we apply a reordering to reduce fillin. There
    // is no advantage in doing this for the other solvers.
    if (solver_type == SolverType::FGMRES_ILU)
      {
        TimerOutput::Scope ilu_specific(computing_timer, "(ILU specific)");
        DoFRenumbering::Cuthill_McKee(dof_handler);
      }

    // This ensures that all velocities DoFs are enumerated before the
    // pressure unknowns. This allows us to use blocks for vectors and
    // matrices and allows us to get the same DoF numbering for
    // dof_handler and velocity_dof_handler.
    DoFRenumbering::block_wise(dof_handler);

    if (solver_type == SolverType::FGMRES_GMG)
      {
        TimerOutput::Scope multigrid_specific(computing_timer,
                                              "(Multigrid specific)");
        TimerOutput::Scope setup_multigrid(computing_timer,
                                           "Setup - Multigrid");

        // This distributes the active dofs and multigrid dofs for the
        // velocity space in a separate DoFHandler as described in the
        // introduction.
        velocity_dof_handler.distribute_dofs(velocity_fe);
        velocity_dof_handler.distribute_mg_dofs();

        // The following block of code initializes the MGConstrainedDofs
        // (using the boundary conditions for the velocity), and the
        // sparsity patterns and matrices for each level. The resize()
        // function of MGLevelObject<T> will destroy all existing contained
        // objects.
        std::set<types::boundary_id> zero_boundary_ids;
        zero_boundary_ids.insert(0);

        mg_constrained_dofs.clear();
        mg_constrained_dofs.initialize(velocity_dof_handler);
        mg_constrained_dofs.make_zero_boundary_constraints(velocity_dof_handler,
                                                           zero_boundary_ids);
        const unsigned int n_levels = triangulation.n_levels();

        mg_interface_matrices.resize(0, n_levels - 1);
        mg_matrices.resize(0, n_levels - 1);
        mg_sparsity_patterns.resize(0, n_levels - 1);

        for (unsigned int level = 0; level < n_levels; ++level)
          {
            DynamicSparsityPattern csp(velocity_dof_handler.n_dofs(level),
                                       velocity_dof_handler.n_dofs(level));
            MGTools::make_sparsity_pattern(velocity_dof_handler, csp, level);
            mg_sparsity_patterns[level].copy_from(csp);

            mg_matrices[level].reinit(mg_sparsity_patterns[level]);
            mg_interface_matrices[level].reinit(mg_sparsity_patterns[level]);
          }
      }

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    {
      constraints.clear();
      // The following makes use of a component mask for interpolation of the
      // boundary values for the velocity only, which is further explained in
      // the vector valued dealii step-20 tutorial.
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Solution<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));

      // As discussed in the introduction, we need to fix one degree of freedom
      // of the pressure variable to ensure solvability of the problem. We do
      // this here by marking the first pressure dof, which has index n_u as a
      // constrained dof.
      if (solver_type == SolverType::UMFPACK)
        constraints.add_line(n_u);

      constraints.close();
    }

    std::cout << "\tNumber of active cells: " << triangulation.n_active_cells()
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


  // @sect4{StokesProblem::assemble_system}

  // In this function, the system matrix is assembled. We assemble the pressure
  // mass matrix in the (1,1) block (if needed) and move it out of this location
  // at the end of this function.
  template <int dim>
  void StokesProblem<dim>::assemble_system()
  {
    TimerOutput::Scope assemble(computing_timer, "Assemble");
    system_matrix = 0;
    system_rhs    = 0;

    // If true, we will assemble the pressure mass matrix in the (1,1) block:
    const bool assemble_pressure_mass_matrix =
      (solver_type == SolverType::UMFPACK) ? false : true;

    QGauss<dim> quadrature_formula(pressure_degree + 2);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;

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

    if (solver_type != SolverType::UMFPACK)
      {
        pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
        pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
        system_matrix.block(1, 1) = 0;
      }
  }

  // @sect4{StokesProblem::assemble_multigrid}

  // Here, like in step-16, we have a function that assembles the level
  // and interface matrices necessary for the multigrid preconditioner.
  template <int dim>
  void StokesProblem<dim>::assemble_multigrid()
  {
    TimerOutput::Scope multigrid_specific(computing_timer,
                                          "(Multigrid specific)");
    TimerOutput::Scope assemble_multigrid(computing_timer,
                                          "Assemble Multigrid");

    mg_matrices = 0.;

    QGauss<dim> quadrature_formula(pressure_degree + 2);

    FEValues<dim> fe_values(velocity_fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);

    const unsigned int dofs_per_cell = velocity_fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);

    std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);

    std::vector<AffineConstraints<double>> boundary_constraints(
      triangulation.n_levels());
    std::vector<AffineConstraints<double>> boundary_interface_constraints(
      triangulation.n_levels());
    for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
      {
        boundary_constraints[level].add_lines(
          mg_constrained_dofs.get_refinement_edge_indices(level));
        boundary_constraints[level].add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        boundary_constraints[level].close();

        IndexSet idx = mg_constrained_dofs.get_refinement_edge_indices(level) &
                       mg_constrained_dofs.get_boundary_indices(level);

        boundary_interface_constraints[level].add_lines(idx);
        boundary_interface_constraints[level].close();
      }

    // This iterator goes over all cells (not just active)
    for (const auto &cell : velocity_dof_handler.cell_iterators())
      {
        fe_values.reinit(cell);
        cell_matrix = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient(k, q);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j <= i; ++j)
                {
                  cell_matrix(i, j) +=
                    (symgrad_phi_u[i] * symgrad_phi_u[j]) * fe_values.JxW(q);
                }
          }

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            cell_matrix(i, j) = cell_matrix(j, i);

        cell->get_mg_dof_indices(local_dof_indices);

        boundary_constraints[cell->level()].distribute_local_to_global(
          cell_matrix, local_dof_indices, mg_matrices[cell->level()]);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            if (!mg_constrained_dofs.at_refinement_edge(cell->level(),
                                                        local_dof_indices[i]) ||
                mg_constrained_dofs.at_refinement_edge(cell->level(),
                                                       local_dof_indices[j]))
              cell_matrix(i, j) = 0;

        boundary_interface_constraints[cell->level()]
          .distribute_local_to_global(cell_matrix,
                                      local_dof_indices,
                                      mg_interface_matrices[cell->level()]);
      }
  }

  // @sect4{StokesProblem::solve}

  // This function sets up things differently based on if you want to use ILU
  // or GMG as a preconditioner.  Both methods share the same solver (FGMRES)
  // but require a different preconditioner to be initialized. Here we time not
  // only the entire solve function, but we separately time the setup of the
  // preconditioner as well as the solve itself.
  template <int dim>
  void StokesProblem<dim>::solve()
  {
    TimerOutput::Scope solve(computing_timer, "Solve");
    constraints.set_zero(solution);

    if (solver_type == SolverType::UMFPACK)
      {
        computing_timer.enter_subsection("(UMFPACK specific)");
        computing_timer.enter_subsection("Solve - Initialize");

        SparseDirectUMFPACK A_direct;
        A_direct.initialize(system_matrix);

        computing_timer.leave_subsection();
        computing_timer.leave_subsection();

        {
          TimerOutput::Scope solve_backslash(computing_timer,
                                             "Solve - Backslash");
          A_direct.vmult(solution, system_rhs);
        }

        constraints.distribute(solution);
        return;
      }

    // Here we must make sure to solve for the residual with "good enough"
    // accuracy
    SolverControl solver_control(system_matrix.m(),
                                 1e-10 * system_rhs.l2_norm());
    unsigned int  n_iterations_A;
    unsigned int  n_iterations_S;

    // This is used to pass whether or not we want to solve for A inside
    // the preconditioner.  One could change this to false to see if
    // there is still convergence and if so does the program then run
    // faster or slower
    const bool use_expensive = true;

    SolverFGMRES<BlockVector<double>> solver(solver_control);

    if (solver_type == SolverType::FGMRES_ILU)
      {
        computing_timer.enter_subsection("(ILU specific)");
        computing_timer.enter_subsection("Solve - Set-up Preconditioner");

        std::cout << "   Computing preconditioner..." << std::endl
                  << std::flush;

        SparseILU<double> A_preconditioner;
        A_preconditioner.initialize(system_matrix.block(0, 0));

        SparseILU<double> S_preconditioner;
        S_preconditioner.initialize(pressure_mass_matrix);

        const BlockSchurPreconditioner<SparseILU<double>, SparseILU<double>>
          preconditioner(system_matrix,
                         pressure_mass_matrix,
                         A_preconditioner,
                         S_preconditioner,
                         use_expensive);

        computing_timer.leave_subsection();
        computing_timer.leave_subsection();

        {
          TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES");

          solver.solve(system_matrix, solution, system_rhs, preconditioner);
          n_iterations_A = preconditioner.n_iterations_A;
          n_iterations_S = preconditioner.n_iterations_S;
        }
      }
    else
      {
        computing_timer.enter_subsection("(Multigrid specific)");
        computing_timer.enter_subsection("Solve - Set-up Preconditioner");

        // Transfer operators between levels
        MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs);
        mg_transfer.build(velocity_dof_handler);

        // Setup coarse grid solver
        FullMatrix<double> coarse_matrix;
        coarse_matrix.copy_from(mg_matrices[0]);
        MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
        coarse_grid_solver.initialize(coarse_matrix);

        using Smoother = PreconditionSOR<SparseMatrix<double>>;
        mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother;
        mg_smoother.initialize(mg_matrices);
        mg_smoother.set_steps(2);

        // Multigrid, when used as a preconditioner for CG, needs to be a
        // symmetric operator, so the smoother must be symmetric
        mg_smoother.set_symmetric(true);

        mg::Matrix<Vector<double>> mg_matrix(mg_matrices);
        mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices);
        mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices);

        // Now, we are ready to set up the V-cycle operator and the multilevel
        // preconditioner.
        Multigrid<Vector<double>> mg(
          mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother);
        mg.set_edge_matrices(mg_interface_down, mg_interface_up);

        PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
          A_Multigrid(velocity_dof_handler, mg, mg_transfer);

        SparseILU<double> S_preconditioner;
        S_preconditioner.initialize(pressure_mass_matrix,
                                    SparseILU<double>::AdditionalData());

        const BlockSchurPreconditioner<
          PreconditionMG<dim,
                         Vector<double>,
                         MGTransferPrebuilt<Vector<double>>>,
          SparseILU<double>>
          preconditioner(system_matrix,
                         pressure_mass_matrix,
                         A_Multigrid,
                         S_preconditioner,
                         use_expensive);

        computing_timer.leave_subsection();
        computing_timer.leave_subsection();

        {
          TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES");
          solver.solve(system_matrix, solution, system_rhs, preconditioner);
          n_iterations_A = preconditioner.n_iterations_A;
          n_iterations_S = preconditioner.n_iterations_S;
        }
      }

    constraints.distribute(solution);

    std::cout
      << std::endl
      << "\tNumber of FGMRES iterations: " << solver_control.last_step()
      << std::endl
      << "\tTotal number of iterations used for approximation of A inverse: "
      << n_iterations_A << std::endl
      << "\tTotal number of iterations used for approximation of S inverse: "
      << n_iterations_S << std::endl
      << std::endl;
  }


  // @sect4{StokesProblem::process_solution}

  // This function computes the L2 and H1 errors of the solution. For this,
  // we need to make sure the pressure has mean zero.
  template <int dim>
  void StokesProblem<dim>::compute_errors()
  {
    // Compute the mean pressure $\frac{1}{\Omega} \int_{\Omega} p(x) dx $
    // and then subtract it from each pressure coefficient. This will result
    // in a pressure with mean value zero. Here we make use of the fact that
    // the pressure is component $dim$ and that the finite element space
    // is nodal.
    const double mean_pressure = VectorTools::compute_mean_value(
      dof_handler, QGauss<dim>(pressure_degree + 2), solution, dim);
    solution.block(1).add(-mean_pressure);
    std::cout << "   Note: The mean value was adjusted by " << -mean_pressure
              << std::endl;

    const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim + 1);

    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 2),
                                      VectorTools::L2_norm,
                                      &velocity_mask);

    const double Velocity_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 2),
                                      VectorTools::L2_norm,
                                      &pressure_mask);

    const double Pressure_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 2),
                                      VectorTools::H1_norm,
                                      &velocity_mask);

    const double Velocity_H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_norm);

    std::cout << std::endl
              << "   Velocity L2 Error: " << Velocity_L2_error << std::endl
              << "   Pressure L2 Error: " << Pressure_L2_error << std::endl
              << "   Velocity H1 Error: " << Velocity_H1_error << std::endl;
  }


  // @sect4{StokesProblem::output_results}

  // This function generates graphical output like it is done in step-22.
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

    std::ofstream output(
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk");
    data_out.write_vtk(output);
  }



  // @sect4{StokesProblem::run}

  // The last step in the Stokes class is, as usual, the function that
  // generates the initial grid and calls the other functions in the
  // respective order.
  template <int dim>
  void StokesProblem<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(6 - dim);

    if (solver_type == SolverType::FGMRES_ILU)
      std::cout << "Now running with ILU" << std::endl;
    else if (solver_type == SolverType::FGMRES_GMG)
      std::cout << "Now running with Multigrid" << std::endl;
    else
      std::cout << "Now running with UMFPACK" << std::endl;


    for (unsigned int refinement_cycle = 0; refinement_cycle < 3;
         ++refinement_cycle)
      {
        std::cout << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
          triangulation.refine_global(1);

        std::cout << "   Set-up..." << std::endl;
        setup_dofs();

        std::cout << "   Assembling..." << std::endl;
        assemble_system();

        if (solver_type == SolverType::FGMRES_GMG)
          {
            std::cout << "   Assembling Multigrid..." << std::endl;

            assemble_multigrid();
          }

        std::cout << "   Solving..." << std::flush;
        solve();

        compute_errors();

        output_results(refinement_cycle);

        Utilities::System::MemoryStats mem;
        Utilities::System::get_memory_stats(mem);
        std::cout << "   VM Peak: " << mem.VmPeak << std::endl;

        computing_timer.print_summary();
        computing_timer.reset();
      }
  }
} // namespace Step56

// @sect3{The main function}
int main()
{
  try
    {
      using namespace Step56;

      const int degree = 1;
      const int dim    = 3;
      // options for SolverType: UMFPACK FGMRES_ILU FGMRES_GMG
      StokesProblem<dim> flow_problem(degree, SolverType::FGMRES_GMG);

      flow_problem.run();
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
