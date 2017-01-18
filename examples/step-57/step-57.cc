/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Liang Zhao and Timo Heister, Clemson University, 2016
 */

// @sect3{Include files}

// As usual, we start by including some well-known files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>

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

// To transfer solutions between meshes, this file is included:
#include <deal.II/numerics/solution_transfer.h>

// This file includes UMFPACK: the direct solver:
#include <deal.II/lac/sparse_direct.h>

// And the one for ILU preconditioner:
#include <deal.II/lac/sparse_ilu.h>


#include <fstream>
#include <iostream>
#include <sstream>

namespace Step57
{
  using namespace dealii;

  // @sect3{The <code>NavierStokesProblem</code> class template}

  // This is the main function and its member functions.
  // As explained in the introduction, what we obtain at each step is the
  // Newton's update, so we define two variables: the present solution
  // and the update. Additionally, the evaluation point is
  // for temporarily holding Newton update in line search. A sparse matrix
  // for the pressure mass matrix is created for the operator of a block Schur
  // complement preconditioner. We use one ConstraintMatrix for Dirichlet boundary
  // conditions at the initial step and a zero ConstraintMatrix for the Newton
  // is defined by 1/Re which has been discussed in the introduction.

  template <int dim>
  class StationaryNavierStokes
  {
  public:
    StationaryNavierStokes(const unsigned int degree);
    void run(const unsigned int refinement);

  private:
    void setup_dofs();
    void initialize_system();
    void assemble_system(const bool initial_step,
                         const bool assemble_matrix,
                         const bool assemble_rhs);
    void assemble_matrix(const bool initial_step);
    void assemble_rhs(const bool initial_step);
    void solve(const bool initial_step);
    void refine_mesh();
    void process_solution(unsigned int refinement);
    void output_results (const unsigned int refinement_cycle) const;
    void newton_iteration(const double tolerance,
                          const unsigned int max_iteration,
                          const unsigned int n_refinements,
                          const bool is_initial_step,
                          const bool output_result);
    void compute_initial_guess(double step_size);

    double viscosity;
    double gamma;
    const unsigned int           degree;
    std::vector<types::global_dof_index> dofs_per_block;

    Triangulation<dim>           triangulation;
    FESystem<dim>                fe;
    DoFHandler<dim>              dof_handler;

    ConstraintMatrix             zero_constraints;
    ConstraintMatrix             nonzero_constraints;

    BlockSparsityPattern         sparsity_pattern;
    BlockSparseMatrix<double>    system_matrix;
    SparseMatrix<double>         pressure_mass_matrix;

    BlockVector<double>          present_solution;
    BlockVector<double>          newton_update;
    BlockVector<double>          system_rhs;
    BlockVector<double>          evaluation_point;

  };

  // @sect3{Boundary values and right hand side}
  // In this problem we set the velocity along the upper surface of the cavity
  // to be one and zero on the other three walls. The right hand side function
  // is zero so we do not need to set the right hand side function in this
  // tutorial. The number of components of the boundary function is dim+1.
  // In practice, the boundary values are
  // applied to our solution through ConstraintMatrix which is obtained by using
  // VectorTools::interpolate_boundary_values. The components of boundary value
  // functions are required to be chosen according to the finite element space.
  // Therefore we have to define the boundary value of pressure even though we
  // actually do not need it.

  // The following function represents the boundary values:
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues() : Function<dim>(dim+1) {}
    virtual double value(const Point<dim> &p,
                         const unsigned int component) const;

    virtual void   vector_value(const Point <dim>    &p,
                                Vector<double> &values) const;
  };

  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> &p,
                                    const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));
    if (component == 0 && std::abs(p[dim-1]-1.0) < 1e-10)
      return 1.0;

    return 0;
  }

  template <int dim>
  void BoundaryValues<dim>::vector_value ( const Point<dim> &p,
                                           Vector<double>   &values ) const
  {
    for (unsigned int c = 0; c < this->n_components; ++c)
      values(c) = BoundaryValues<dim>::value (p, c);
  }


  // @sect3{BlockSchurPreconditioner for Navier Stokes equations}
  //
  // The block
  // Schur complement preconditioner is defined in this part. As discussed in
  // the introduction, the preconditioner in Krylov iterative methods is
  // implemented as a matrix-vector product operator. In practice, the Schur
  // complement preconditioner is decomposed as a product of three matrices (as
  // presented in the first section).  The $\tilde{A}^{-1}$ in the first factor
  // involves a solve for the linear system $\tilde{A}x=b$. Here we solve
  // this system via a direct solver for simplicity. The computation involved
  // in the second factor is a simple matrix-vector multiplication. The Schur
  // complement $\tilde{S}$ can be well approximated by the pressure mass
  // matrix and its inverse can be obtained through an inexact solver. Because
  // the pressure mass matrix is symmetric and positive definite, we can use
  // CG to solve the corresponding linear system.
  //
  template <class PreconditionerMp>
  class BlockSchurPreconditioner : public Subscriptor
  {
  public:
    BlockSchurPreconditioner (double                                     gamma,
                              double                                     viscosity,
                              const BlockSparseMatrix<double>            &S,
                              const SparseMatrix<double>                 &P,
                              const PreconditionerMp                     &Mppreconditioner);

    void vmult (BlockVector<double>       &dst,
                const BlockVector<double> &src) const;

  private:
    const double gamma;
    const double viscosity;
    const BlockSparseMatrix<double> &stokes_matrix;
    const SparseMatrix<double> &pressure_mass_matrix;
    const PreconditionerMp &mp_preconditioner;
    SparseDirectUMFPACK A_inverse;
  };

  // We can notice that the initialization of the inverse of the matrix at (0,0) corner
  // is completed in the constructor. If so, every application of the preconditioner then
  // no longer requires the computation of the matrix factors.

  template <class PreconditionerMp>
  BlockSchurPreconditioner<PreconditionerMp>::
  BlockSchurPreconditioner (double                           gamma,
                            double                           viscosity,
                            const BlockSparseMatrix<double>  &S,
                            const SparseMatrix<double>       &P,
                            const PreconditionerMp           &Mppreconditioner)
    :
    gamma                (gamma),
    viscosity            (viscosity),
    stokes_matrix        (S),
    pressure_mass_matrix (P),
    mp_preconditioner    (Mppreconditioner)
  {
    A_inverse.initialize(stokes_matrix.block(0,0));
  }

  template <class PreconditionerMp>
  void
  BlockSchurPreconditioner<PreconditionerMp>::
  vmult (BlockVector<double>       &dst,
         const BlockVector<double> &src) const
  {
    Vector<double> utmp(src.block(0));

    {
      SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
      SolverCG<>    cg (solver_control);

      dst.block(1) = 0.0;
      cg.solve(pressure_mass_matrix,
               dst.block(1), src.block(1),
               mp_preconditioner);
      dst.block(1) *= -(viscosity+gamma);
    }

    {
      stokes_matrix.block(0,1).vmult(utmp, dst.block(1));
      utmp *= -1.0;
      utmp += src.block(0);
    }

    A_inverse.vmult (dst.block(0), utmp);
  }

  // @sect3{StationaryNavierStokes class implementation}
  // @sect4{StationaryNavierStokes::StationaryNavierStokes}
  // The constructor of this class looks very similar to the one in step-22. The only
  // difference is the viscosity and the Augmented Lagrangian coefficient gamma.
  //

  template <int dim>
  StationaryNavierStokes<dim>::StationaryNavierStokes(const unsigned int degree)
    :
    viscosity(1.0/7500.0),
    gamma(1.0),
    degree(degree),
    triangulation(Triangulation<dim>::maximum_smoothing),
    fe(FE_Q<dim>(degree+1), dim,
       FE_Q<dim>(degree),   1),
    dof_handler(triangulation)
  {}

  // @sect4{StationaryNavierStokes::setup_dofs}
  // This function initializes the DoFHandler enumerating the degrees of freedom
  // and constraints on the current mesh.

  template <int dim>
  void StationaryNavierStokes<dim>::setup_dofs()
  {
    system_matrix.clear();
    pressure_mass_matrix.clear();

    // The first step is to associate DoFs with a given mesh.
    dof_handler.distribute_dofs (fe);

    // We renumber the components to have all velocity DoFs come before
    // the pressure DoFs to be able to split the solution vector in two blocks
    // which are separately accessed in the block preconditioner.
    //
    std::vector<unsigned int> block_component(dim+1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise (dof_handler, block_component);

    dofs_per_block.resize (2);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
    unsigned int dof_u = dofs_per_block[0];
    unsigned int dof_p = dofs_per_block[1];

    // In Newton's scheme, we first apply the boundary condition on the solution
    // obtained from the initial step. To make sure the boundary conditions remain
    // satisfied during Newton's iteration, zero boundary conditions are used for
    // the update $\delta u^k$. Therefore we set up two different constraint objects.

    FEValuesExtractors::Vector velocities(0);
    {
      nonzero_constraints.clear();

      DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               BoundaryValues<dim>(),
                                               nonzero_constraints,
                                               fe.component_mask(velocities));
    }
    nonzero_constraints.close();

    {
      zero_constraints.clear();

      DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               ZeroFunction<dim>(dim+1),
                                               zero_constraints,
                                               fe.component_mask(velocities));
    }
    zero_constraints.close();

    std::cout << "   Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (" << dof_u << '+' << dof_p << ')'
              << std::endl;
  }

  // @sect4{StationaryNavierStokes::initialize_system}
  // On each mesh the sparsity pattern and the size of the linear system
  // are different. This function initializes them after mesh refinement.

  template <int dim>
  void StationaryNavierStokes<dim>::initialize_system()
  {
    {
      BlockDynamicSparsityPattern dsp (dofs_per_block, dofs_per_block);
      DoFTools::make_sparsity_pattern (dof_handler, dsp, nonzero_constraints);
      sparsity_pattern.copy_from (dsp);
    }

    system_matrix.reinit (sparsity_pattern);

    present_solution.reinit (dofs_per_block);
    newton_update.reinit (dofs_per_block);
    system_rhs.reinit (dofs_per_block);
  }

  // @sect4{StationaryNavierStokes::assemble_system}

  // This function builds the system matrix and right hand side that we
  // actually work on. "initial_step" is given for applying different
  // constraints (nonzero for the initial step and zero for the others). The
  // other two flags are to determine whether to assemble the system matrix
  // or the right hand side vector, respectively.

  template <int dim>
  void StationaryNavierStokes<dim>::assemble_system(const bool initial_step,
                                                    const bool assemble_matrix,
                                                    const bool assemble_rhs)
  {
    if (assemble_matrix)
      system_matrix    = 0;

    if (assemble_rhs)
      system_rhs       = 0;

    QGauss<dim>   quadrature_formula(degree+2);

    FEValues<dim> fe_values (fe,
                             quadrature_formula,
                             update_values |
                             update_quadrature_points |
                             update_JxW_values |
                             update_gradients );

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs    (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    // For the linearized system, we create temporary storage for present velocity
    // and gradient, and present pressure. In practice, they are
    // all obtained through their shape functions at quadrature points.

    std::vector<Tensor<1, dim> >  present_velocity_values    (n_q_points);
    std::vector<Tensor<2, dim> >  present_velocity_gradients (n_q_points);
    std::vector<double>           present_pressure_values    (n_q_points);

    std::vector<double>           div_phi_u                 (dofs_per_cell);
    std::vector<Tensor<1, dim> >  phi_u                     (dofs_per_cell);
    std::vector<Tensor<2, dim> >  grad_phi_u                (dofs_per_cell);
    std::vector<double>           phi_p                     (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        fe_values.reinit(cell);

        local_matrix = 0;
        local_rhs    = 0;

        fe_values[velocities].get_function_values(evaluation_point,
                                                  present_velocity_values);

        fe_values[velocities].get_function_gradients(evaluation_point,
                                                     present_velocity_gradients);

        fe_values[pressure].get_function_values(evaluation_point,
                                                present_pressure_values);

        // The assembly is similar to step-22. An additional term with gamma as a coefficient
        // is the Augmented Lagrangian (AL), which is assembled via grad-div stabilization.
        // As we discussed in the introduction, the bottom right block of the system matrix should be
        // zero. Since the pressure mass matrix is used while creating the preconditioner,
        // we assemble it here and then move it into a separate SparseMatrix at the end (same as in step-22).

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                phi_u[k]      =  fe_values[velocities].value(k, q);
                phi_p[k]      =  fe_values[pressure]  .value(k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                if (assemble_matrix)
                  {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      {
                        local_matrix(i, j) += (  viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])
                                                 + present_velocity_gradients[q]*phi_u[j]*phi_u[i]
                                                 + grad_phi_u[j]*present_velocity_values[q]*phi_u[i]
                                                 - div_phi_u[i]*phi_p[j]
                                                 - phi_p[i]*div_phi_u[j]
                                                 + gamma*div_phi_u[j]*div_phi_u[i]
                                                 + phi_p[i] * phi_p[j])
                                              * fe_values.JxW(q);
                      }
                  }

                if (assemble_rhs)
                  {
                    double present_velocity_divergence =  trace(present_velocity_gradients[q]);
                    local_rhs(i) += ( - viscosity*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                                      - present_velocity_gradients[q]*present_velocity_values[q]*phi_u[i]
                                      + present_pressure_values[q]*div_phi_u[i]
                                      + present_velocity_divergence*phi_p[i]
                                      - gamma*present_velocity_divergence*div_phi_u[i])
                                    * fe_values.JxW(q);
                  }
              }
          }

        cell->get_dof_indices (local_dof_indices);

        const ConstraintMatrix &constraints_used = initial_step ? nonzero_constraints : zero_constraints;
        // Finally we move pressure mass matrix into a separate matrix:

        if (assemble_matrix)
          {
            constraints_used.distribute_local_to_global(local_matrix,
                                                        local_dof_indices,
                                                        system_matrix);
          }

        if (assemble_rhs)
          {
            constraints_used.distribute_local_to_global(local_rhs,
                                                        local_dof_indices,
                                                        system_rhs);
          }

      }

    if (assemble_matrix)
      {
        pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
        pressure_mass_matrix.copy_from(system_matrix.block(1,1));

        // Note that settings this pressure block to zero is not identical to
        // not assembling anything in this block, because this operation here
        // will (incorrectly) delete diagonal entries that come in from
        // hanging node constraints for pressure DoFs. This means that our
        // whole system matrix will have rows that are completely
        // zero. Luckily, FGMRES handles these rows without any problem.
        system_matrix.block(1,1) = 0;
      }
  }

  template <int dim>
  void StationaryNavierStokes<dim>::assemble_matrix(const bool initial_step)
  {
    assemble_system(initial_step, true, false);
  }

  template <int dim>
  void StationaryNavierStokes<dim>::assemble_rhs(const bool initial_step)
  {
    assemble_system(initial_step, false, true);
  }
  // @sect4{StationaryNavierStokes::solve}
  // In this function, we use FGMRES together with the block preconditioner,
  // which is defined at the beginning of the program, to solve the linear
  // system. What we obtain at this step is the solution vector. If this is
  // the initial step, the solution vector gives us an initial guess for the
  // Navier Stokes equations. For the initial step, nonzero constraints are
  // applied in order to make sure boundary conditions are satisfied. In the
  // following steps, we will solve for the Newton update so zero
  // constraints are used.
  template <int dim>
  void StationaryNavierStokes<dim>::solve (const bool initial_step)
  {
    const ConstraintMatrix &constraints_used = initial_step ? nonzero_constraints : zero_constraints;

    SolverControl solver_control (system_matrix.m(), 1e-4*system_rhs.l2_norm(), true);
    SolverFGMRES<BlockVector<double> > gmres(solver_control);

    SparseILU<double> pmass_preconditioner;
    pmass_preconditioner.initialize (pressure_mass_matrix,
                                     SparseILU<double>::AdditionalData());

    const BlockSchurPreconditioner<SparseILU<double> >
    preconditioner (gamma,
                    viscosity,
                    system_matrix,
                    pressure_mass_matrix,
                    pmass_preconditioner);

    gmres.solve (system_matrix,
                 newton_update,
                 system_rhs,
                 preconditioner);
    std::cout << " ****FGMRES steps: " << solver_control.last_step() << std::endl;

    constraints_used.distribute(newton_update);
  }

  // @sect4{StationaryNavierStokes::refine_mesh}
  //
  // After finding a good initial guess on the coarse mesh, we hope to
  // decrease the error through refining the mesh. Here we do adaptive
  // refinement similar to step-15 except that we use the Kelly estimator on
  // the velocity only. We also need to transfer the current solution to the
  // next mesh using the SolutionTransfer class.
  template <int dim>
  void StationaryNavierStokes<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    FEValuesExtractors::Vector velocity(0);
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(degree+1),
                                        typename FunctionMap<dim>::type(),
                                        present_solution,
                                        estimated_error_per_cell,
                                        fe.component_mask(velocity));

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.0);

    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim, BlockVector<double> > solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
    triangulation.execute_coarsening_and_refinement ();

    // First the DoFHandler is set up and constraints are generated. Then we
    // create a temporary vector "tmp", whose size is according with the
    // solution on the new mesh.
    setup_dofs();

    BlockVector<double> tmp (dofs_per_block);

    // Transfer solution from coarse to fine mesh and apply boundary value
    // constraints to the new transfered solution. Note that present_solution
    // is still a vector corresponding to the old mesh.
    solution_transfer.interpolate(present_solution, tmp);
    nonzero_constraints.distribute(tmp);

    // Finally set up matrix and vectors and set the present_solution to the
    // interpolated data.
    initialize_system();
    present_solution = tmp;
  }

  // @sect4{StationaryNavierStokes<dim>::newton_iteration}
  //
  // This function implements the Newton iteration with given tolerance, maximum number of iterations,
  // and the number of mesh refinements to do. "is_initial_step" is the flag to tell us whether
  // "setup_system" is necessary, and which part,
  // system matrix or right hand side vector, should be assembled. If we do a line search,
  // the right hand side is already assembled while checking the residual norm in the last iteration.
  // Therefore, we just need to assemble the system matrix at the current iteration. The last
  // argument "output_result" is whether output should be produced.

  template <int dim>
  void StationaryNavierStokes<dim>::newton_iteration(const double tolerance,
                                                     const unsigned int max_iteration,
                                                     const unsigned int max_refinement,
                                                     const bool  is_initial_step,
                                                     const bool  output_result)
  {
    double current_res;
    double last_res;
    bool   first_step = is_initial_step;

    for (unsigned int refinement = 0; refinement < max_refinement+1; ++refinement)
      {
        unsigned int outer_iteration = 0;
        last_res = 1.0;
        current_res = 1.0;
        std::cout << "*****************************************" << std::endl;
        std::cout << "************  refinement = " << refinement << " ************ " << std::endl;
        std::cout << "viscosity= " << viscosity << std::endl;
        std::cout << "*****************************************" << std::endl;

        while ((first_step || (current_res > tolerance)) && outer_iteration < max_iteration)
          {
            if (first_step)
              {
                setup_dofs();
                initialize_system();
                evaluation_point = present_solution;
                assemble_matrix(first_step);
                assemble_rhs(first_step);
                solve(first_step);
                present_solution = newton_update;
                nonzero_constraints.distribute(present_solution);
                first_step = false;
                evaluation_point = present_solution;
                assemble_rhs(first_step);
                current_res = system_rhs.l2_norm();
                std::cout << "******************************" << std::endl;
                std::cout << " The residual of initial guess is " << current_res << std::endl;
                std::cout << " Initialization complete!  " << std::endl;
                last_res = current_res;
              }

            else
              {
                evaluation_point = present_solution;
                assemble_matrix(first_step);
                if (outer_iteration == 0)
                  assemble_rhs(first_step);
                solve(first_step);

                // To make sure our solution is getting close to the exact solution, we
                // let the solution be updated with a weight alpha such
                // that the new residual is smaller than the one of last step, which is done
                // in the following loop. Also the line search method can be located in step-15.

                for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5)
                  {
                    evaluation_point = present_solution;
                    evaluation_point.add(alpha, newton_update);
                    nonzero_constraints.distribute(evaluation_point);
                    assemble_rhs(first_step);
                    current_res = system_rhs.l2_norm();
                    std::cout << " alpha = " << std::setw(6) << alpha << std::setw(0)
                              << " res = " << current_res << std::endl;
                    if (current_res < last_res)
                      break;
                  }
                {
                  present_solution = evaluation_point;
                  std::cout << " ---- Iteration " << outer_iteration << " residual: " << current_res << std::endl;
                  last_res = current_res;
                }

              }
            ++outer_iteration;

            if (output_result)
              {
                output_results (max_iteration*refinement+outer_iteration);

                if (current_res <= tolerance)
                  process_solution(refinement);
              }
          }

        if (refinement < max_refinement)
          {
            refine_mesh();
            std::cout << "*****************************************" << std::endl
                      << "        Do refinement ------   " << std::endl;
          }
      }


  }

  // @sect4{StationaryNavierStokes::compute_initial_guess}
  //
  // This function will provide us with an initial guess by using a
  // continuation method as we discussed in the introduction. The Reynolds
  // number is increased step-by-step until we reach the target value. By
  // experiment, the solution to Stokes is good enough to be the initial guess
  // of NSE with Reynolds number 1000 so we start there.  To make sure the
  // solution from previous problem is close enough to the next one, the step
  // size must be small enough.
  template <int dim>
  void StationaryNavierStokes<dim>::compute_initial_guess(double step_size)
  {
    const double target_Re = 1.0/viscosity;

    bool is_initial_step = true;

    for (double Re=1000.0; Re < target_Re; Re = std::min(Re+step_size, target_Re))
      {
        viscosity = 1.0/Re;
        std::cout << "*****************************************" << std::endl;
        std::cout << " Searching for initial guess with Re = " << Re << std::endl;
        std::cout << "*****************************************" << std::endl;

        newton_iteration(1e-12, 50, 0, is_initial_step, false);
        is_initial_step = false;
      }
  }

  // @sect4{StationaryNavierStokes::output_results}
  // This function is the same as in step-22.
  template <int dim>
  void StationaryNavierStokes<dim>::output_results (const unsigned int output_index)  const
  {
    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (present_solution, solution_names,
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);
    data_out.build_patches ();

    std::ostringstream filename;
    filename << 1.0/viscosity
             << "-solution-"
             << Utilities::int_to_string (output_index, 4)
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
  }

  // @sect4{StationaryNavierStokes::process_solution}
  // In our test case, we do not know the analytical solution. This function
  // outputs the velocity components along x=0.5 and y from 0 to 1 so they
  // can be compared with data from the literature.
  template <int dim>
  void StationaryNavierStokes<dim>::process_solution(unsigned int refinement)
  {
    std::ostringstream filename;
    filename << (1.0/viscosity) << "-line-" << refinement << ".txt";

    std::ofstream f (filename.str().c_str());
    f << "# y u_x u_y" << std::endl;

    Point<dim> p;
    p(0)= 0.5;
    p(1)= 0.5;

    f << std::scientific;

    for (unsigned int i=0; i<=100; ++i)
      {
        p(dim-1) = i/100.0;

        Vector<double> tmp_vector(dim+1);
        VectorTools::point_value(dof_handler, present_solution, p, tmp_vector);
        f << p(dim-1);

        for (int j=0; j<dim; j++)
          f << " " << tmp_vector(j);
        f << std::endl;
      }
  }


  // @sect4{StationaryNavierStokes::run}
  // This is the last step of this program. In this part, we generate the grid and run
  // the other functions respectively. The max refinement can be set by the argument.
  template <int dim>
  void StationaryNavierStokes<dim>::run(const unsigned int refinement)
  {

    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(5);

    const double Re =  1.0/viscosity;

    // If the viscosity is smaller than 1/1000, we have to first search for an
    // initial guess via a continuation method. What we should notice is the
    // search is always on the initial mesh, that is the $8 \times 8$ mesh in
    // this program. After that, we just do the same as we did when viscosity
    // is larger than 1/1000: run Newton's iteration, refine the mesh,
    // transfer solutions, and repeat.
    if (Re > 1000.0)
      {
        std::cout << "       Searching for initial guess ... " << std::endl;
        const double step_size = 2000.0;
        compute_initial_guess(step_size);
        std::cout << "*****************************************" << std::endl
                  << "       Initial guess obtained            " << std::endl
                  << "                  *                      " << std::endl
                  << "                  *                      " << std::endl
                  << "                  *                      " << std::endl
                  << "                  *                      " << std::endl
                  << "*****************************************" << std::endl;

        std::cout << "       Computing solution with target Re = " << Re << std::endl;
        viscosity = 1.0/Re;
        newton_iteration(1e-12, 50, refinement, false, true);
      }
    else
      {
        // When the viscosity is larger than 1/1000, the solution to Stokes
        // equations is good enough as an initial guess. If so, we do not need
        // to search for the initial guess using a continuation
        // method. Newton's iteration can be started directly.

        newton_iteration(1e-12, 50, refinement, true, true);
      }
  }
}

int main()
{
  try
    {
      using namespace dealii;
      using namespace Step57;

      StationaryNavierStokes<2> flow(/* degree = */1);
      flow.run(4);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
