/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2016 - 2020 by the deal.II authors
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

 *
 * Author: Timo Heister, Clemson University, 2016
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

// The following chunk out code is identical to step-40 and allows
// switching between PETSc and Trilinos:

#include <deal.II/lac/generic_linear_algebra.h>

/* #define FORCE_USE_OF_TRILINOS */

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <cmath>
#include <fstream>
#include <iostream>

namespace Step55
{
  using namespace dealii;

  // @sect3{Linear solvers and preconditioners}

  // We need a few helper classes to represent our solver strategy
  // described in the introduction.

  namespace LinearSolvers
  {
    // This class exposes the action of applying the inverse of a giving
    // matrix via the function InverseMatrix::vmult(). Internally, the
    // inverse is not formed explicitly. Instead, a linear solver with CG
    // is performed. This class extends the InverseMatrix class in step-22
    // with an option to specify a preconditioner, and to allow for different
    // vector types in the vmult function.
    template <class Matrix, class Preconditioner>
    class InverseMatrix : public Subscriptor
    {
    public:
      InverseMatrix(const Matrix &m, const Preconditioner &preconditioner);

      template <typename VectorType>
      void vmult(VectorType &dst, const VectorType &src) const;

    private:
      const SmartPointer<const Matrix> matrix;
      const Preconditioner &           preconditioner;
    };


    template <class Matrix, class Preconditioner>
    InverseMatrix<Matrix, Preconditioner>::InverseMatrix(
      const Matrix &        m,
      const Preconditioner &preconditioner)
      : matrix(&m)
      , preconditioner(preconditioner)
    {}



    template <class Matrix, class Preconditioner>
    template <typename VectorType>
    void
    InverseMatrix<Matrix, Preconditioner>::vmult(VectorType &      dst,
                                                 const VectorType &src) const
    {
      SolverControl solver_control(src.size(), 1e-8 * src.l2_norm());
      SolverCG<LA::MPI::Vector> cg(solver_control);
      dst = 0;

      try
        {
          cg.solve(*matrix, dst, src, preconditioner);
        }
      catch (std::exception &e)
        {
          Assert(false, ExcMessage(e.what()));
        }
    }


    // The class A template class for a simple block diagonal preconditioner
    // for 2x2 matrices.
    template <class PreconditionerA, class PreconditionerS>
    class BlockDiagonalPreconditioner : public Subscriptor
    {
    public:
      BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A,
                                  const PreconditionerS &preconditioner_S);

      void vmult(LA::MPI::BlockVector &      dst,
                 const LA::MPI::BlockVector &src) const;

    private:
      const PreconditionerA &preconditioner_A;
      const PreconditionerS &preconditioner_S;
    };

    template <class PreconditionerA, class PreconditionerS>
    BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>::
      BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A,
                                  const PreconditionerS &preconditioner_S)
      : preconditioner_A(preconditioner_A)
      , preconditioner_S(preconditioner_S)
    {}


    template <class PreconditionerA, class PreconditionerS>
    void BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>::vmult(
      LA::MPI::BlockVector &      dst,
      const LA::MPI::BlockVector &src) const
    {
      preconditioner_A.vmult(dst.block(0), src.block(0));
      preconditioner_S.vmult(dst.block(1), src.block(1));
    }

  } // namespace LinearSolvers

  // @sect3{Problem setup}

  // The following classes represent the right hand side and the exact
  // solution for the test problem.

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim + 1)
    {}

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &  value) const override;
  };


  template <int dim>
  void RightHandSide<dim>::vector_value(const Point<dim> &p,
                                        Vector<double> &  values) const
  {
    const double R_x = p[0];
    const double R_y = p[1];

    const double pi  = numbers::PI;
    const double pi2 = pi * pi;
    values[0] =
      -1.0L / 2.0L * (-2 * sqrt(25.0 + 4 * pi2) + 10.0) *
        exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) -
      0.4 * pi2 * exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) +
      0.1 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 2) *
        exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi);
    values[1] = 0.2 * pi * (-sqrt(25.0 + 4 * pi2) + 5.0) *
                  exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) -
                0.05 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 3) *
                  exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) /
                  pi;
    values[2] = 0;
  }


  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution()
      : Function<dim>(dim + 1)
    {}

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &  value) const override;
  };

  template <int dim>
  void ExactSolution<dim>::vector_value(const Point<dim> &p,
                                        Vector<double> &  values) const
  {
    const double R_x = p[0];
    const double R_y = p[1];

    const double pi  = numbers::PI;
    const double pi2 = pi * pi;
    values[0] =
      -exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 1;
    values[1] = (1.0L / 2.0L) * (-sqrt(25.0 + 4 * pi2) + 5.0) *
                exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) /
                pi;
    values[2] =
      -1.0L / 2.0L * exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) -
      2.0 *
        (-6538034.74494422 +
         0.0134758939981709 * exp(4 * sqrt(25.0 + 4 * pi2))) /
        (-80.0 * exp(3 * sqrt(25.0 + 4 * pi2)) +
         16.0 * sqrt(25.0 + 4 * pi2) * exp(3 * sqrt(25.0 + 4 * pi2))) -
      1634508.68623606 * exp(-3.0 * sqrt(25.0 + 4 * pi2)) /
        (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)) +
      (-0.00673794699908547 * exp(sqrt(25.0 + 4 * pi2)) +
       3269017.37247211 * exp(-3 * sqrt(25.0 + 4 * pi2))) /
        (-8 * sqrt(25.0 + 4 * pi2) + 40.0) +
      0.00336897349954273 * exp(1.0 * sqrt(25.0 + 4 * pi2)) /
        (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2));
  }



  // @sect3{The main program}
  //
  // The main class is very similar to step-40, except that matrices and
  // vectors are now block versions, and we store a std::vector<IndexSet>
  // for owned and relevant DoFs instead of a single IndexSet. We have
  // exactly two IndexSets, one for all velocity unknowns and one for all
  // pressure unknowns.
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(unsigned int velocity_degree);

    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

    unsigned int velocity_degree;
    double       viscosity;
    MPI_Comm     mpi_communicator;

    FESystem<dim>                             fe;
    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;

    AffineConstraints<double> constraints;

    LA::MPI::BlockSparseMatrix system_matrix;
    LA::MPI::BlockSparseMatrix preconditioner_matrix;
    LA::MPI::BlockVector       locally_relevant_solution;
    LA::MPI::BlockVector       system_rhs;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem(unsigned int velocity_degree)
    : velocity_degree(velocity_degree)
    , viscosity(0.1)
    , mpi_communicator(MPI_COMM_WORLD)
    , fe(FE_Q<dim>(velocity_degree), dim, FE_Q<dim>(velocity_degree - 1), 1)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {}


  // The Kovasnay flow is defined on the domain [-0.5, 1.5]^2, which we
  // create by passing the min and max values to GridGenerator::hyper_cube.
  template <int dim>
  void StokesProblem<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, -0.5, 1.5);
    triangulation.refine_global(3);
  }

  // @sect3{System Setup}
  //
  // The construction of the block matrices and vectors is new compared to
  // step-40 and is different compared to serial codes like step-22, because
  // we need to supply the set of rows that belong to our processor.
  template <int dim>
  void StokesProblem<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup");

    dof_handler.distribute_dofs(fe);

    // Put all dim velocities into block 0 and the pressure into block 1,
    // then reorder the unknowns by block. Finally count how many unknowns
    // we have per block.
    std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
    stokes_sub_blocks[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks);

    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
          << n_u << '+' << n_p << ')' << std::endl;

    // We split up the IndexSet for locally owned and locally relevant DoFs
    // into two IndexSets based on how we want to create the block matrices
    // and vectors.
    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u);
    owned_partitioning[1] =
      dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_p);

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    // Setting up the constraints for boundary conditions and hanging nodes
    // is identical to step-40. Rven though we don't have any hanging nodes
    // because we only perform global refinement, it is still a good idea
    // to put this function call in, in case adaptive refinement gets
    // introduced later.
    {
      constraints.reinit(locally_relevant_dofs);

      FEValuesExtractors::Vector velocities(0);
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               ExactSolution<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
      constraints.close();
    }

    // Now we create the system matrix based on a BlockDynamicSparsityPattern.
    // We know that we won't have coupling between different velocity
    // components (because we use the laplace and not the deformation tensor)
    // and no coupling between pressure with its test functions, so we use
    // a Table to communicate this coupling information to
    // DoFTools::make_sparsity_pattern.
    {
      system_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
      for (unsigned int c = 0; c < dim + 1; ++c)
        for (unsigned int d = 0; d < dim + 1; ++d)
          if (c == dim && d == dim)
            coupling[c][d] = DoFTools::none;
          else if (c == dim || d == dim || c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        dof_handler, coupling, dsp, constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        dof_handler.locally_owned_dofs(),
        mpi_communicator,
        locally_relevant_dofs);

      system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);
    }

    // The preconditioner matrix has a different coupling (we only fill in
    // the 1,1 block with the mass matrix), otherwise this code is identical
    // to the construction of the system_matrix above.
    {
      preconditioner_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
      for (unsigned int c = 0; c < dim + 1; ++c)
        for (unsigned int d = 0; d < dim + 1; ++d)
          if (c == dim && d == dim)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        dof_handler, coupling, dsp, constraints, false);
      SparsityTools::distribute_sparsity_pattern(
        dsp,
        Utilities::MPI::all_gather(mpi_communicator,
                                   dof_handler.locally_owned_dofs()),
        mpi_communicator,
        locally_relevant_dofs);
      preconditioner_matrix.reinit(owned_partitioning,
                                   //      owned_partitioning,
                                   dsp,
                                   mpi_communicator);
    }

    // Finally, we construct the block vectors with the right sizes. The
    // function call with two std::vector<IndexSet> will create a ghosted
    // vector.
    locally_relevant_solution.reinit(owned_partitioning,
                                     relevant_partitioning,
                                     mpi_communicator);
    system_rhs.reinit(owned_partitioning, mpi_communicator);
  }



  // @sect3{Assembly}
  //
  // This function assembles the system matrix, the preconditioner matrix,
  // and the right hand side. The code is pretty standard.
  template <int dim>
  void StokesProblem<dim>::assemble_system()
  {
    TimerOutput::Scope t(computing_timer, "assembly");

    system_matrix         = 0;
    preconditioner_matrix = 0;
    system_rhs            = 0;

    const QGauss<dim> quadrature_formula(velocity_degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    const RightHandSide<dim>    right_hand_side;
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1));

    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double>         div_phi_u(dofs_per_cell);
    std::vector<double>         phi_p(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Vector     velocities(0);
    const FEValuesExtractors::Scalar     pressure(dim);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix  = 0;
          cell_matrix2 = 0;
          cell_rhs     = 0;

          fe_values.reinit(cell);
          right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                            rhs_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                  phi_p[k]      = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) +=
                        (viscosity *
                           scalar_product(grad_phi_u[i], grad_phi_u[j]) -
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                        fe_values.JxW(q);

                      cell_matrix2(i, j) += 1.0 / viscosity * phi_p[i] *
                                            phi_p[j] * fe_values.JxW(q);
                    }

                  const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                  cell_rhs(i) += fe_values.shape_value(i, q) *
                                 rhs_values[q](component_i) * fe_values.JxW(q);
                }
            }


          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);

          constraints.distribute_local_to_global(cell_matrix2,
                                                 local_dof_indices,
                                                 preconditioner_matrix);
        }

    system_matrix.compress(VectorOperation::add);
    preconditioner_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }



  // @sect3{Solving}
  //
  // This function solves the linear system with MINRES with a block diagonal
  // preconditioner and AMG for the two diagonal blocks as described in the
  // introduction. The preconditioner applies a v cycle to the 0,0 block
  // and a CG with the mass matrix for the 1,1 block (the Schur complement).
  template <int dim>
  void StokesProblem<dim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "solve");

    LA::MPI::PreconditionAMG prec_A;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_A.initialize(system_matrix.block(0, 0), data);
    }

    LA::MPI::PreconditionAMG prec_S;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_S.initialize(preconditioner_matrix.block(1, 1), data);
    }

    // The InverseMatrix is used to solve for the mass matrix:
    using mp_inverse_t = LinearSolvers::InverseMatrix<LA::MPI::SparseMatrix,
                                                      LA::MPI::PreconditionAMG>;
    const mp_inverse_t mp_inverse(preconditioner_matrix.block(1, 1), prec_S);

    // This constructs the block preconditioner based on the preconditioners
    // for the individual blocks defined above.
    const LinearSolvers::BlockDiagonalPreconditioner<LA::MPI::PreconditionAMG,
                                                     mp_inverse_t>
      preconditioner(prec_A, mp_inverse);

    // With that, we can finally set up the linear solver and solve the system:
    SolverControl solver_control(system_matrix.m(),
                                 1e-10 * system_rhs.l2_norm());

    SolverMinRes<LA::MPI::BlockVector> solver(solver_control);

    LA::MPI::BlockVector distributed_solution(owned_partitioning,
                                              mpi_communicator);

    constraints.set_zero(distributed_solution);

    solver.solve(system_matrix,
                 distributed_solution,
                 system_rhs,
                 preconditioner);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(distributed_solution);

    // Like in step-56, we subtract the mean pressure to allow error
    // computations against our reference solution, which has a mean value
    // of zero.
    locally_relevant_solution = distributed_solution;
    const double mean_pressure =
      VectorTools::compute_mean_value(dof_handler,
                                      QGauss<dim>(velocity_degree + 2),
                                      locally_relevant_solution,
                                      dim);
    distributed_solution.block(1).add(-mean_pressure);
    locally_relevant_solution.block(1) = distributed_solution.block(1);
  }



  // @sect3{The rest}
  //
  // The remainder of the code that deals with mesh refinement, output, and
  // the main loop is pretty standard.
  template <int dim>
  void StokesProblem<dim>::refine_grid()
  {
    TimerOutput::Scope t(computing_timer, "refine");

    triangulation.refine_global();
  }



  template <int dim>
  void StokesProblem<dim>::output_results(const unsigned int cycle) const
  {
    {
      const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
      const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                       dim + 1);

      Vector<double> cellwise_errors(triangulation.n_active_cells());
      QGauss<dim>    quadrature(velocity_degree + 2);

      VectorTools::integrate_difference(dof_handler,
                                        locally_relevant_solution,
                                        ExactSolution<dim>(),
                                        cellwise_errors,
                                        quadrature,
                                        VectorTools::L2_norm,
                                        &velocity_mask);

      const double error_u_l2 =
        VectorTools::compute_global_error(triangulation,
                                          cellwise_errors,
                                          VectorTools::L2_norm);

      VectorTools::integrate_difference(dof_handler,
                                        locally_relevant_solution,
                                        ExactSolution<dim>(),
                                        cellwise_errors,
                                        quadrature,
                                        VectorTools::L2_norm,
                                        &pressure_mask);

      const double error_p_l2 =
        VectorTools::compute_global_error(triangulation,
                                          cellwise_errors,
                                          VectorTools::L2_norm);

      pcout << "error: u_0: " << error_u_l2 << " p_0: " << error_p_l2
            << std::endl;
    }


    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.emplace_back("pressure");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    LA::MPI::BlockVector interpolated;
    interpolated.reinit(owned_partitioning, MPI_COMM_WORLD);
    VectorTools::interpolate(dof_handler, ExactSolution<dim>(), interpolated);

    LA::MPI::BlockVector interpolated_relevant(owned_partitioning,
                                               relevant_partitioning,
                                               MPI_COMM_WORLD);
    interpolated_relevant = interpolated;
    {
      std::vector<std::string> solution_names(dim, "ref_u");
      solution_names.emplace_back("ref_p");
      data_out.add_data_vector(interpolated_relevant,
                               solution_names,
                               DataOut<dim>::type_dof_data,
                               data_component_interpretation);
    }


    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, mpi_communicator, 2);
  }



  template <int dim>
  void StokesProblem<dim>::run()
  {
#ifdef USE_PETSC_LA
    pcout << "Running using PETSc." << std::endl;
#else
    pcout << "Running using Trilinos." << std::endl;
#endif
    const unsigned int n_cycles = 5;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          make_grid();
        else
          refine_grid();

        setup_system();

        assemble_system();
        solve();

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results(cycle);
          }

        computing_timer.print_summary();
        computing_timer.reset();

        pcout << std::endl;
      }
  }
} // namespace Step55



int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step55;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      StokesProblem<2> problem(2);
      problem.run();
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
