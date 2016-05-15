/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2011 - 2015 by the deal.II authors
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
 * Author: Katharina Kormann, Martin Kronbichler, Uppsala University, 2011-2012
 */


// The necessary files from the deal.II library.
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/distributed/tria.h>

// This includes the data structures for the efficient implementation of
// matrix-free methods.
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <fstream>
#include <iostream>
#include <iomanip>


namespace Step48
{
  using namespace dealii;

  // We start by defining two global variables to collect all parameters
  // subject to changes at one place: One for the dimension and one for the
  // finite element degree. The dimension is used in the main function as a
  // template argument for the actual classes (like in all other deal.II
  // programs), whereas the degree of the finite element is more crucial, as
  // it is passed as a template argument to the implementation of the
  // Sine-Gordon operator. Therefore, it needs to be a compile-time constant.
  const unsigned int dimension = 2;
  const unsigned int fe_degree = 4;


  // @sect3{SineGordonOperation}

  // The <code>SineGordonOperation</code> class implements the cell-based
  // operation that is needed in each time step. This nonlinear operation can
  // be implemented straight-forwardly based on the <code>MatrixFree</code>
  // class, in the same way as a linear operation would be treated by this
  // implementation of the finite element operator application. We apply two
  // template arguments to the class, one for the dimension and one for the
  // degree of the finite element. This is a difference to other functions in
  // deal.II where only the dimension is a template argument. This is
  // necessary to provide the inner loops in @p FEEvaluation with information
  // about loop lengths etc., which is essential for efficiency. On the other
  // hand, it makes it more challenging to implement the degree as a run-time
  // parameter.
  template <int dim, int fe_degree>
  class SineGordonOperation
  {
  public:
    SineGordonOperation(const MatrixFree<dim,double> &data_in,
                        const double                  time_step);

    void apply (parallel::distributed::Vector<double>                     &dst,
                const std::vector<parallel::distributed::Vector<double>*> &src) const;

  private:
    const MatrixFree<dim,double>         &data;
    const VectorizedArray<double>         delta_t_sqr;
    parallel::distributed::Vector<double> inv_mass_matrix;

    void local_apply (const MatrixFree<dim,double>               &data,
                      parallel::distributed::Vector<double>      &dst,
                      const std::vector<parallel::distributed::Vector<double>*> &src,
                      const std::pair<unsigned int,unsigned int> &cell_range) const;
  };



  // @sect4{SineGordonOperation::SineGordonOperation}

  // This is the constructor of the SineGordonOperation class. It receives a
  // reference to the MatrixFree holding the problem information and the time
  // step size as input parameters. The initialization routine sets up the
  // mass matrix. Since we use Gauss-Lobatto elements, the mass matrix is a
  // diagonal matrix and can be stored as a vector. The computation of the
  // mass matrix diagonal is simple to achieve with the data structures
  // provided by FEEvaluation: Just loop over all (macro-) cells and integrate
  // over the function that is constant one on all quadrature points by using
  // the <code>integrate</code> function with @p true argument at the slot for
  // values. Finally, we invert the diagonal entries since we have to multiply
  // by the inverse mass matrix in each time step.
  template <int dim, int fe_degree>
  SineGordonOperation<dim,fe_degree>::
  SineGordonOperation(const MatrixFree<dim,double> &data_in,
                      const double                  time_step)
    :
    data(data_in),
    delta_t_sqr(make_vectorized_array(time_step *time_step))
  {
    VectorizedArray<double> one = make_vectorized_array (1.);

    data.initialize_dof_vector (inv_mass_matrix);

    FEEvaluation<dim,fe_degree>   fe_eval(data);
    const unsigned int            n_q_points = fe_eval.n_q_points;

    for (unsigned int cell=0; cell<data.n_macro_cells(); ++cell)
      {
        fe_eval.reinit(cell);
        for (unsigned int q=0; q<n_q_points; ++q)
          fe_eval.submit_value(one,q);
        fe_eval.integrate (true,false);
        fe_eval.distribute_local_to_global (inv_mass_matrix);
      }

    inv_mass_matrix.compress(VectorOperation::add);
    for (unsigned int k=0; k<inv_mass_matrix.local_size(); ++k)
      if (inv_mass_matrix.local_element(k)>1e-15)
        inv_mass_matrix.local_element(k) = 1./inv_mass_matrix.local_element(k);
      else
        inv_mass_matrix.local_element(k) = 0;
  }



  // @sect4{SineGordonOperation::local_apply}

  // This operator implements the core operation of the program, the
  // integration over a range of cells for the nonlinear operator of the
  // Sine-Gordon problem. The implementation is based on the FEEvaluation
  // class as in step-37. Due to the special structure in Gauss-Lobatto
  // elements, certain operations become simpler, in particular the evaluation
  // of shape function values on quadrature points which is simply the
  // injection of the values of cell degrees of freedom. The MatrixFree class
  // detects possible structure of the finite element at quadrature points
  // when initializing, which is then used by FEEvaluation for selecting the
  // most appropriate numerical kernel.

  // The nonlinear function that we have to evaluate for the time stepping
  // routine includes the value of the function at the present time @p current
  // as well as the value at the previous time step @p old. Both values are
  // passed to the operator in the collection of source vectors @p src, which is
  // simply a <tt>std::vector</tt> of pointers to the actual solution
  // vectors. This construct of collecting several source vectors into one is
  // necessary as the cell loop in @p MatrixFree takes exactly one source and
  // one destination vector, even if we happen to use many vectors like the two
  // in this case. Note that the cell loop accepts any valid class for input and
  // output, which does not only include vectors but general data types.
  // However, only in case it encounters a parallel::distributed::Vector<Number>
  // or a <tt>std::vector</tt> collecting these vectors, it calls functions that
  // exchange data at the beginning and the end of the loop. In the loop over
  // the cells, we first have to read in the values in the vectors related to
  // the local values. Then, we evaluate the value and the gradient of the
  // current solution vector and the values of the old vector at the quadrature
  // points. Then, we combine the terms in the scheme in the loop over the
  // quadrature points. Finally, we integrate the result against the test
  // function and accumulate the result to the global solution vector @p dst.
  template <int dim, int fe_degree>
  void SineGordonOperation<dim, fe_degree>::
  local_apply (const MatrixFree<dim>                      &data,
               parallel::distributed::Vector<double>      &dst,
               const std::vector<parallel::distributed::Vector<double>*> &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    AssertDimension (src.size(), 2);
    FEEvaluation<dim,fe_degree> current (data), old (data);
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        current.reinit (cell);
        old.reinit (cell);

        current.read_dof_values (*src[0]);
        old.read_dof_values     (*src[1]);

        current.evaluate (true, true, false);
        old.evaluate (true, false, false);

        for (unsigned int q=0; q<current.n_q_points; ++q)
          {
            const VectorizedArray<double> current_value = current.get_value(q);
            const VectorizedArray<double> old_value     = old.get_value(q);

            current.submit_value (2.*current_value - old_value -
                                  delta_t_sqr * std::sin(current_value),q);
            current.submit_gradient (- delta_t_sqr *
                                     current.get_gradient(q), q);
          }

        current.integrate (true,true);
        current.distribute_local_to_global (dst);
      }
  }



  //@sect4{SineGordonOperation::apply}

  // This function performs the time stepping routine based on the cell-local
  // strategy. First the destination vector is set to zero, then the cell-loop
  // is called, and finally the solution is multiplied by the inverse mass
  // matrix. The structure of the cell loop is implemented in the cell finite
  // element operator class. On each cell it applies the routine defined as
  // the <code>local_apply()</code> method of the class
  // <code>SineGordonOperation</code>, i.e., <code>this</code>. One could also
  // provide a function with the same signature that is not part of a class.
  template <int dim, int fe_degree>
  void SineGordonOperation<dim, fe_degree>::
  apply (parallel::distributed::Vector<double>                     &dst,
         const std::vector<parallel::distributed::Vector<double>*> &src) const
  {
    dst = 0;
    data.cell_loop (&SineGordonOperation<dim,fe_degree>::local_apply,
                    this, dst, src);
    dst.scale(inv_mass_matrix);
  }


  //@sect3{Equation data}

  // We define a time-dependent function that is used as initial
  // value. Different solutions can be obtained by varying the starting
  // time. This function has already been explained in step-25.
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution (const unsigned int n_components = 1,
                   const double time = 0.) : Function<dim>(n_components, time) {}
    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
  };

  template <int dim>
  double ExactSolution<dim>::value (const Point<dim> &p,
                                    const unsigned int /* component */) const
  {
    double t = this->get_time ();

    const double m = 0.5;
    const double c1 = 0.;
    const double c2 = 0.;
    const double factor = (m / std::sqrt(1.-m*m) *
                           std::sin(std::sqrt(1.-m*m)*t+c2));
    double result = 1.;
    for (unsigned int d=0; d<dim; ++d)
      result *= -4. * std::atan (factor / std::cosh(m*p[d]+c1));
    return result;
  }



  // @sect3{SineGordonProblem class}

  // This is the main class that builds on the class in step-25.  However, we
  // replaced the SparseMatrix<double> class by the MatrixFree class to store
  // the geometry data. Also, we use a distributed triangulation in this
  // example.
  template <int dim>
  class SineGordonProblem
  {
  public:
    SineGordonProblem ();
    void run ();

  private:
    ConditionalOStream pcout;

    void make_grid_and_dofs ();
    void output_results (const unsigned int timestep_number);

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim>   triangulation;
#else
    Triangulation<dim>   triangulation;
#endif
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    ConstraintMatrix     constraints;
    IndexSet             locally_relevant_dofs;

    MatrixFree<dim,double> matrix_free_data;

    parallel::distributed::Vector<double> solution, old_solution, old_old_solution;

    const unsigned int n_global_refinements;
    double time, time_step;
    const double final_time;
    const double cfl_number;
    const unsigned int output_timestep_skip;
  };


  //@sect4{SineGordonProblem::SineGordonProblem}

  // This is the constructor of the SineGordonProblem class. The time interval
  // and time step size are defined here. Moreover, we use the degree of the
  // finite element that we defined at the top of the program to initialize a
  // FE_Q finite element based on Gauss-Lobatto support points. These points
  // are convenient because in conjunction with a QGaussLobatto quadrature
  // rule of the same order they give a diagonal mass matrix without
  // compromising accuracy too much (note that the integration is inexact,
  // though), see also the discussion in the introduction. Note that FE_Q
  // selects the Gauss-Lobatto nodal points by default due to their improved
  // conditioning versus equidistant points. To make things more explicit, we
  // choose to state the selection of the nodal points nonetheless.
  template <int dim>
  SineGordonProblem<dim>::SineGordonProblem ()
    :
    pcout (std::cout,
           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
#ifdef DEAL_II_WITH_P4EST
    triangulation (MPI_COMM_WORLD),
#endif
    fe (QGaussLobatto<1>(fe_degree+1)),
    dof_handler (triangulation),
    n_global_refinements (10-2*dim),
    time (-10),
    final_time (10),
    cfl_number (.1/fe_degree),
    output_timestep_skip (200)
  {}

  //@sect4{SineGordonProblem::make_grid_and_dofs}

  // As in step-25 this functions sets up a cube grid in <code>dim</code>
  // dimensions of extent $[-15,15]$. We refine the mesh more in the center of
  // the domain since the solution is concentrated there. We first refine all
  // cells whose center is within a radius of 11, and then refine once more
  // for a radius 6.  This simple ad hoc refinement could be done better by
  // adapting the mesh to the solution using error estimators during the time
  // stepping as done in other example programs, and using
  // parallel::distributed::SolutionTransfer to transfer the solution to the
  // new mesh.
  template <int dim>
  void SineGordonProblem<dim>::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation, -15, 15);
    triangulation.refine_global (n_global_refinements);
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active(),
      end_cell = triangulation.end();
      for ( ; cell != end_cell; ++cell)
        if (cell->is_locally_owned())
          if (cell->center().norm() < 11)
            cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();

      cell = triangulation.begin_active();
      end_cell = triangulation.end();
      for ( ; cell != end_cell; ++cell)
        if (cell->is_locally_owned())
          if (cell->center().norm() < 6)
            cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }

    pcout << "   Number of global active cells: "
#ifdef DEAL_II_WITH_P4EST
          << triangulation.n_global_active_cells()
#else
          << triangulation.n_active_cells()
#endif
          << std::endl;

    dof_handler.distribute_dofs (fe);

    pcout << "   Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;


    // We generate hanging node constraints for ensuring continuity of the
    // solution. As in step-40, we need to equip the constraint matrix with
    // the IndexSet of locally relevant degrees of freedom to avoid it to
    // consume too much memory for big problems. Next, the <code> MatrixFree
    // </code> for the problem is set up. Note that we specify the MPI
    // communicator which we are going to use, and that we also want to use
    // shared-memory parallelization (hence one would use multithreading for
    // intra-node parallelism and not MPI; note that we here choose the
    // standard option &mdash; if we wanted to disable shared memory
    // parallelization, we would choose @p none). Finally, three solution
    // vectors are initialized. MatrixFree stores the layout that is to be
    // used by distributed vectors, so we just ask it to initialize the
    // vectors.
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    constraints.clear();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close();

    QGaussLobatto<1> quadrature (fe_degree+1);
    typename MatrixFree<dim>::AdditionalData additional_data;
    additional_data.mpi_communicator = MPI_COMM_WORLD;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim>::AdditionalData::partition_partition;

    matrix_free_data.reinit (dof_handler, constraints,
                             quadrature, additional_data);

    matrix_free_data.initialize_dof_vector (solution);
    old_solution.reinit (solution);
    old_old_solution.reinit (solution);
  }



  //@sect4{SineGordonProblem::output_results}

  // This function prints the norm of the solution and writes the solution
  // vector to a file. The norm is standard (except for the fact that we need
  // to be sure to only count norms on locally owned cells), and the second is
  // similar to what we did in step-40. Note that we can use the same vector
  // for output as we used for computation: The vectors in the matrix-free
  // framework always provide full information on all locally owned cells
  // (this is what is needed in the local evaluations, too), including ghost
  // vector entries on these cells. This is the only data that is needed in
  // the integrate_difference function as well as in DataOut. We only need to
  // make sure that we tell the vector to update its ghost values before we
  // read them. This is a feature present only in the
  // parallel::distributed::Vector class. Distributed vectors with PETSc and
  // Trilinos, on the other hand, need to be copied to special vectors
  // including ghost values (see the relevant section in step-40). If we
  // wanted to access all degrees of freedom on ghost cells, too (e.g. when
  // computing error estimators that use the jump of solution over cell
  // boundaries), we would need more information and create a vector
  // initialized with locally relevant dofs just as in step-40. Observe also
  // that we need to distribute constraints for output - they are not filled
  // during computations (rather, they are distributed on the fly in the
  // matrix-free method read_dof_values).
  template <int dim>
  void
  SineGordonProblem<dim>::output_results (const unsigned int timestep_number)
  {
    constraints.distribute (solution);

    Vector<float> norm_per_cell (triangulation.n_active_cells());
    solution.update_ghost_values();
    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       ZeroFunction<dim>(),
                                       norm_per_cell,
                                       QGauss<dim>(fe_degree+1),
                                       VectorTools::L2_norm);
    const double solution_norm =
      VectorTools::compute_global_error(triangulation,
                                        norm_per_cell,
                                        VectorTools::L2_norm);

    pcout << "   Time:"
          << std::setw(8) << std::setprecision(3) << time
          << ", solution norm: "
          << std::setprecision(5) << std::setw(7) << solution_norm
          << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();

    const std::string filename =
      "solution-" + Utilities::int_to_string (timestep_number, 3);

    std::ofstream output ((filename +
                           "." + Utilities::int_to_string (Utilities::MPI::
                                                           this_mpi_process(MPI_COMM_WORLD),4) + ".vtu").c_str());
    data_out.write_vtu (output);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
          filenames.push_back ("solution-" +
                               Utilities::int_to_string (timestep_number, 3) +
                               "." +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");

        std::ofstream master_output ((filename + ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }
  }


  // @sect4{SineGordonProblem::run}

  // This function is called by the main function and calls the subroutines of
  // the class.
  //
  // The first step is to set up the grid and the cell operator. Then, the
  // time step is computed from the CFL number given in the constructor and
  // the finest mesh size. The finest mesh size is computed as the diameter of
  // the last cell in the triangulation, which is the last cell on the finest
  // level of the mesh. This is only possible for Cartesian meshes, otherwise,
  // one needs to loop over all cells. Note that we need to query all the
  // processors for their finest cell since the not all processors might hold
  // a region where the mesh is at the finest level. Then, we readjust the
  // time step a little to hit the final time exactly.
  template <int dim>
  void
  SineGordonProblem<dim>::run ()
  {
    make_grid_and_dofs();

    const double local_min_cell_diameter =
      triangulation.last()->diameter()/std::sqrt(dim);
    const double global_min_cell_diameter
      = -Utilities::MPI::max(-local_min_cell_diameter, MPI_COMM_WORLD);
    time_step = cfl_number * global_min_cell_diameter;
    time_step = (final_time-time)/(int((final_time-time)/time_step));
    pcout << "   Time step size: " << time_step << ", finest cell: "
          << global_min_cell_diameter << std::endl << std::endl;

    // Next the initial value is set. Since we have a two-step time stepping
    // method, we also need a value of the solution at time-time_step. For
    // accurate results, one would need to compute this from the time
    // derivative of the solution at initial time, but here we ignore this
    // difficulty and just set it to the initial value function at that
    // artificial time.

    // We create an output of the initial value. Then we also need to collect
    // the two starting solutions in a <tt>std::vector</tt> of pointers field
    // and to set up an instance of the <code> SineGordonOperation class </code>
    // based on the finite element degree specified at the top of this file.
    VectorTools::interpolate (dof_handler,
                              ExactSolution<dim> (1, time),
                              solution);
    VectorTools::interpolate (dof_handler,
                              ExactSolution<dim> (1, time-time_step),
                              old_solution);
    output_results (0);

    std::vector<parallel::distributed::Vector<double>*> previous_solutions;
    previous_solutions.push_back(&old_solution);
    previous_solutions.push_back(&old_old_solution);

    SineGordonOperation<dim,fe_degree> sine_gordon_op (matrix_free_data,
                                                       time_step);

    // Now loop over the time steps. In each iteration, we shift the solution
    // vectors by one and call the <code> apply </code> function of the <code>
    // SineGordonOperator </code>. Then, we write the solution to a file. We
    // clock the wall times for the computational time needed as wall as the
    // time needed to create the output and report the numbers when the time
    // stepping is finished.
    //
    // Note how this shift is implemented: We simply call the swap method on
    // the two vectors which swaps only some pointers without the need to copy
    // data around. Obviously, this is a more efficient way to update the
    // vectors during time stepping. Let us see what happens in more detail:
    // First, we exchange <code>old_solution</code> with
    // <code>old_old_solution</code>, which means that
    // <code>old_old_solution</code> gets <code>old_solution</code>, which is
    // what we expect. Similarly, <code>old_solution</code> gets the content
    // from <code>solution</code> in the next step. Afterward,
    // <code>solution</code> holds <code>old_old_solution</code>, but that
    // will be overwritten during this step.
    unsigned int timestep_number = 1;

    Timer timer;
    double wtime = 0;
    double output_time = 0;
    for (time+=time_step; time<=final_time; time+=time_step, ++timestep_number)
      {
        timer.restart();
        old_old_solution.swap (old_solution);
        old_solution.swap (solution);
        sine_gordon_op.apply (solution, previous_solutions);
        wtime += timer.wall_time();

        timer.restart();
        if (timestep_number % output_timestep_skip == 0)
          output_results(timestep_number / output_timestep_skip);

        output_time += timer.wall_time();
      }
    timer.restart();
    output_results(timestep_number / output_timestep_skip + 1);
    output_time += timer.wall_time();

    pcout << std::endl
          << "   Performed " << timestep_number << " time steps."
          << std::endl;

    pcout << "   Average wallclock time per time step: "
          << wtime / timestep_number << "s" << std::endl;

    pcout << "   Spent " << output_time << "s on output and "
          << wtime << "s on computations." << std::endl;
  }
}



// @sect3{The <code>main</code> function}

// As in step-40, we initialize MPI at the start of the program. Since we will
// in general mix MPI parallelization with threads, we also set the third
// argument in MPI_InitFinalize that controls the number of threads to an
// invalid number, which means that the TBB library chooses the number of
// threads automatically, typically to the number of available cores in the
// system. As an alternative, you can also set this number manually if you
// want to set a specific number of threads (e.g. when MPI-only is required).
int main (int argc, char **argv)
{
  using namespace Step48;
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv,
                                                       numbers::invalid_unsigned_int);

  try
    {
      SineGordonProblem<dimension> sg_problem;
      sg_problem.run ();
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
