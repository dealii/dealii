/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2013
 */


#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>


namespace Step26
{
  using namespace dealii;



  template<int dim>
  class HeatEquation
  {
  public:
    HeatEquation();
    void run();

  private:
    void setup_system();
    void solve_time_step();
    void output_results() const;
    void refine_mesh (const unsigned int min_grid_level,
                      const unsigned int max_grid_level);

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       old_solution;
    Vector<double>       system_rhs;

    double               time;
    double               time_step;
    unsigned int         timestep_number;

    const double         theta;
  };



  template<int dim>
  class RightHandSide: public Function<dim>
  {
  public:
    RightHandSide()
      :
      Function<dim>(),
      period (0.2)
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;

  private:
    const double period;
  };



  template<int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int component) const
  {
    Assert (component == 0, ExcInternalError());
    Assert (dim == 2, ExcNotImplemented());

    const double time = this->get_time();
    const double point_within_period = (time/period - std::floor(time/period));

    if ((point_within_period >= 0.0) && (point_within_period <= 0.2))
      {
        if ((p[0] > 0.5) && (p[1] > -0.5))
          return 1;
        else
          return 0;
      }
    else if ((point_within_period >= 0.5) && (point_within_period <= 0.7))
      {
        if ((p[0] > -0.5) && (p[1] > 0.5))
          return 1;
        else
          return 0;
      }
    else
      return 0;
  }



  template<int dim>
  class BoundaryValues: public Function<dim>
  {
  public:
    BoundaryValues()
      :
      Function<dim>()
    {
    }

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
  };

  template<int dim>
  double BoundaryValues<dim>::value(const Point<dim> &/*p*/,
                                    const unsigned int component) const
  {
    Assert(component == 0, ExcInternalError());
    return 0;
  }



  template<int dim>
  HeatEquation<dim>::HeatEquation()
    :
    fe(1),
    dof_handler(triangulation),
    time_step(1. / 500),
    theta(0.5)
  {
  }



  template<int dim>
  void HeatEquation<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::cout << std::endl
              << "==========================================="
              << std::endl
              << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    constraints.close();

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    c_sparsity,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from(c_sparsity);

    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);
    system_matrix.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(fe.degree+1),
                                      mass_matrix,
                                      (const Function<dim> *)0,
                                      constraints);
    MatrixCreator::create_laplace_matrix(dof_handler,
                                         QGauss<dim>(fe.degree+1),
                                         laplace_matrix,
                                         (const Function<dim> *)0,
                                         constraints);

    solution.reinit(dof_handler.n_dofs());
    old_solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }



  template<int dim>
  void HeatEquation<dim>::solve_time_step()
  {
    SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<> cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.0);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    std::cout << "     " << solver_control.last_step()
              << " CG iterations." << std::endl;
  }



  template<int dim>
  void HeatEquation<dim>::output_results() const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "U");

    data_out.build_patches();

    const std::string filename = "solution-"
                                 + Utilities::int_to_string(timestep_number, 3) +
                                 ".vtk";
    std::ofstream output(filename.c_str());
    data_out.write_vtk(output);
  }


  // @sect4{BoussinesqFlowProblem::refine_mesh}
  //
  // This function takes care of the adaptive mesh refinement. The three tasks
  // this function performs is to first find out which cells to
  // refine/coarsen, then to actually do the refinement and eventually
  // transfer the solution vectors between the two different grids. The first
  // task is simply achieved by using the well-established Kelly error
  // estimator on the temperature (it is the temperature we're mainly
  // interested in for this program, and we need to be accurate in regions of
  // high temperature gradients, also to not have too much numerical
  // diffusion). The second task is to actually do the remeshing. That
  // involves only basic functions as well, such as the
  // <code>refine_and_coarsen_fixed_fraction</code> that refines those cells
  // with the largest estimated error that together make up 80 per cent of the
  // error, and coarsens those cells with the smallest error that make up for
  // a combined 10 per cent of the error.
  //
  // If implemented like this, we would get a program that will not make much
  // progress: Remember that we expect temperature fields that are nearly
  // discontinuous (the diffusivity $\kappa$ is very small after all) and
  // consequently we can expect that a freely adapted mesh will refine further
  // and further into the areas of large gradients. This decrease in mesh size
  // will then be accompanied by a decrease in time step, requiring an
  // exceedingly large number of time steps to solve to a given final time. It
  // will also lead to meshes that are much better at resolving
  // discontinuities after several mesh refinement cycles than in the
  // beginning.
  //
  // In particular to prevent the decrease in time step size and the
  // correspondingly large number of time steps, we limit the maximal
  // refinement depth of the mesh. To this end, after the refinement indicator
  // has been applied to the cells, we simply loop over all cells on the
  // finest level and unselect them from refinement if they would result in
  // too high a mesh level.
  template <int dim>
  void HeatEquation<dim>::refine_mesh (const unsigned int min_grid_level,
                                       const unsigned int max_grid_level)
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(fe.degree+1),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                       estimated_error_per_cell,
                                                       0.6, 0.4);
    if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
           cell != triangulation.end(); ++cell)
        cell->clear_refine_flag ();
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active(min_grid_level);
         cell != triangulation.end_active(min_grid_level); ++cell)
      cell->clear_coarsen_flag ();


    // As part of mesh refinement we need to transfer the solution vectors
    // from the old mesh to the new one. To this end we use the
    // SolutionTransfer class and we have to prepare the solution vectors that
    // should be transferred to the new grid (we will lose the old grid once
    // we have done the refinement so the transfer has to happen concurrently
    // with refinement). What we definitely need are the current and the old
    // temperature (BDF-2 time stepping requires two old solutions). Since the
    // SolutionTransfer objects only support to transfer one object per dof
    // handler, we need to collect the two temperature solutions in one data
    // structure. Moreover, we choose to transfer the Stokes solution, too,
    // since we need the velocity at two previous time steps, of which only
    // one is calculated on the fly.
    //
    // Consequently, we initialize two SolutionTransfer objects for the Stokes
    // and temperature DoFHandler objects, by attaching them to the old dof
    // handlers. With this at place, we can prepare the triangulation and the
    // data vectors for refinement (in this order).
    std::vector<Vector<double> > x_solution (2);
    x_solution[0] = solution;
    x_solution[1] = old_solution;

    SolutionTransfer<dim> solution_trans(dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    solution_trans.prepare_for_coarsening_and_refinement(x_solution);

    // Now everything is ready, so do the refinement and recreate the dof
    // structure on the new grid, and initialize the matrix structures and the
    // new vectors in the <code>setup_dofs</code> function. Next, we actually
    // perform the interpolation of the solutions between the grids. We create
    // another copy of temporary vectors for temperature (now corresponding to
    // the new grid), and let the interpolate function do the job. Then, the
    // resulting array of vectors is written into the respective vector member
    // variables. For the Stokes vector, everything is just the same &ndash;
    // except that we do not need another temporary vector since we just
    // interpolate a single vector. In the end, we have to tell the program
    // that the matrices and preconditioners need to be regenerated, since the
    // mesh has changed.
    triangulation.execute_coarsening_and_refinement ();
    setup_system ();

    std::vector<Vector<double> > tmp (2);
    tmp[0].reinit (solution);
    tmp[1].reinit (solution);
    solution_trans.interpolate(x_solution, tmp);

    solution = tmp[0];
    old_solution = tmp[1];
  }



  template<int dim>
  void HeatEquation<dim>::run()
  {
    const unsigned int initial_global_refinement = (dim == 2 ? 1 : 2);
    const unsigned int n_adaptive_pre_refinement_steps = 1;

    GridGenerator::hyper_L (triangulation);
    triangulation.refine_global (initial_global_refinement);

    setup_system();

    unsigned int pre_refinement_step = 0;

    Vector<double> tmp;
    Vector<double> forcing_terms;

start_time_iteration:

    VectorTools::interpolate(dof_handler,
                             ZeroFunction<dim>(),
                             old_solution);
    solution = old_solution;

    timestep_number = 0;
    time            = 0;

    output_results();

    while (time <= 0.5)
      {
        time += time_step;
        ++timestep_number;

        std::cout << "Time step " << timestep_number << " at t=" << time
                  << std::endl;

        tmp.reinit (solution.size());
        forcing_terms.reinit (solution.size());

        mass_matrix.vmult(system_rhs, old_solution);

        laplace_matrix.vmult(tmp, old_solution);
        system_rhs.add(-(1 - theta) * time_step, tmp);

        RightHandSide<dim> rhs_function;
        rhs_function.set_time(time);
        VectorTools::create_right_hand_side(dof_handler,
                                            QGauss<dim>(fe.degree+1),
                                            rhs_function,
                                            tmp);
        forcing_terms = tmp;
        forcing_terms *= time_step * theta;

        rhs_function.set_time(time - time_step);
        VectorTools::create_right_hand_side(dof_handler,
                                            QGauss<dim>(fe.degree+1),
                                            rhs_function,
                                            tmp);

        forcing_terms.add(time_step * (1 - theta), tmp);

        system_rhs += forcing_terms;

        {
          BoundaryValues<dim> boundary_values_function;
          boundary_values_function.set_time(time);

          std::map<types::global_dof_index, double> boundary_values;
          VectorTools::interpolate_boundary_values(dof_handler,
                                                   0,
                                                   boundary_values_function,
                                                   boundary_values);

          system_matrix.copy_from(mass_matrix);
          system_matrix.add(theta * time_step, laplace_matrix);
          MatrixTools::apply_boundary_values(boundary_values,
                                             system_matrix,
                                             solution,
                                             system_rhs);
        }

        constraints.condense (system_rhs);

        solve_time_step();

        output_results();

        if ((timestep_number == 1) &&
            (pre_refinement_step < n_adaptive_pre_refinement_steps))
          {
            refine_mesh (initial_global_refinement,
                         initial_global_refinement + n_adaptive_pre_refinement_steps);
            ++pre_refinement_step;

            std::cout << std::endl;

            goto start_time_iteration;
          }
        else if ((timestep_number > 0) && (timestep_number % 5 == 0))
          refine_mesh (initial_global_refinement,
                       initial_global_refinement + n_adaptive_pre_refinement_steps);

        old_solution = solution;
      }
  }
}

int main()
{
  try
    {
      using namespace dealii;
      using namespace Step26;

      deallog.depth_console(0);

      HeatEquation<2> heat_equation_solver;
      heat_equation_solver.run();

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
                << std::endl << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
