/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2014 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2009, 2010
 *         Timo Heister, University of Goettingen, 2009, 2010
 */


// @sect3{Include files}
//
// Most of the include files we need for this program have already been
// discussed in previous programs. In particular, all of the following should
// already be familiar friends:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>

#define USE_PETSC_LA

namespace LA
{
#ifdef USE_PETSC_LA
  using namespace dealii::LinearAlgebraPETSc;
#else
  using namespace dealii::LinearAlgebraTrilinos;
#endif
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// The following, however, will be new or be used in new roles. Let's walk
// through them. The first of these will provide the tools of the
// Utilities::System namespace that we will use to query things like the
// number of processors associated with the current MPI universe, or the
// number within this universe the processor this job runs on is:
#include <deal.II/base/utilities.h>
// The next one provides a class, ConditionOStream that allows us to write
// code that would output things to a stream (such as <code>std::cout</code>
// on every processor but throws the text away on all but one of them. We
// could achieve the same by simply putting an <code>if</code> statement in
// front of each place where we may generate output, but this doesn't make the
// code any prettier. In addition, the condition whether this processor should
// or should not produce output to the screen is the same every time -- and
// consequently it should be simple enough to put it into the statements that
// generate output itself.
#include <deal.II/base/conditional_ostream.h>
// After these preliminaries, here is where it becomes more interesting. As
// mentioned in the @ref distributed module, one of the fundamental truths of
// solving problems on large numbers of processors is that there is no way for
// any processor to store everything (e.g. information about all cells in the
// mesh, all degrees of freedom, or the values of all elements of the solution
// vector). Rather, every processor will <i>own</i> a few of each of these
// and, if necessary, may <i>know</i> about a few more, for example the ones
// that are located on cells adjacent to the ones this processor owns
// itself. We typically call the latter <i>ghost cells</i>, <i>ghost nodes</i>
// or <i>ghost elements of a vector</i>. The point of this discussion here is
// that we need to have a way to indicate which elements a particular
// processor owns or need to know of. This is the realm of the IndexSet class:
// if there are a total of $N$ cells, degrees of freedom, or vector elements,
// associated with (non-negative) integral indices $[0,N)$, then both the set
// of elements the current processor owns as well as the (possibly larger) set
// of indices it needs to know about are subsets of the set $[0,N)$. IndexSet
// is a class that stores subsets of this set in an efficient format:
#include <deal.II/base/index_set.h>
// The next header file is necessary for a single function,
// SparsityTools::distribute_sparsity_pattern. The role of this function will
// be explained below.
#include <deal.II/lac/sparsity_tools.h>
// The final two, new header files provide the class
// parallel::distributed::Triangulation that provides meshes distributed
// across a potentially very large number of processors, while the second
// provides the namespace parallel::distributed::GridRefinement that offers
// functions that can adaptively refine such distributed meshes:
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

namespace Step40
{
  using namespace dealii;

  // @sect3{The <code>LaplaceProblem</code> class template}

  // Next let's declare the main class of this program. Its structure is
  // almost exactly that of the step-6 tutorial program. The only significant
  // differences are:
  // - The <code>mpi_communicator</code> variable that
  //   describes the set of processors we want this code to run on. In practice,
  //   this will be MPI_COMM_WORLD, i.e. all processors the batch scheduling
  //   system has assigned to this particular job.
  // - The presence of the <code>pcout</code> variable of type ConditionOStream.
  // - The obvious use of parallel::distributed::Triangulation instead of Triangulation.
  // - The presence of two IndexSet objects that denote which sets of degrees of
  //   freedom (and associated elements of solution and right hand side vectors)
  //   we own on the current processor and which we need (as ghost elements) for
  //   the algorithms in this program to work.
  // - The fact that all matrices and vectors are now distributed. We use
  //   their PETScWrapper versions for this since deal.II's own classes do not
  //   provide %parallel functionality. Note that as part of this class, we
  //   store a solution vector that does not only contain the degrees of freedom
  //   the current processor owns, but also (as ghost elements) all those vector
  //   elements that correspond to "locally relevant" degrees of freedom
  //   (i.e. all those that live on locally owned cells or the layer of ghost
  //   cells that surround it).
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    MPI_Comm                                  mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    DoFHandler<dim>                           dof_handler;
    FE_Q<dim>                                 fe;

    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;

    ConstraintMatrix                          constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector locally_relevant_solution;
    LA::MPI::Vector system_rhs;

    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
  };


  // @sect3{The <code>LaplaceProblem</code> class implementation}

  // @sect4{Constructors and destructors}

  // Constructors and destructors are rather trivial. In addition to what we
  // do in step-6, we set the set of processors we want to work on to all
  // machines available (MPI_COMM_WORLD); ask the triangulation to ensure that
  // the mesh remains smooth and free to refined islands, for example; and
  // initialize the <code>pcout</code> variable to only allow processor zero
  // to output anything:
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    dof_handler (triangulation),
    fe (2),
    pcout (std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0)),
    computing_timer (pcout,
                     TimerOutput::summary,
                     TimerOutput::wall_times)
  {}



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem ()
  {
    dof_handler.clear ();
  }


  // @sect4{LaplaceProblem::setup_system}

  // The following function is, arguably, the most interesting one in the
  // entire program since it goes to the heart of what distinguishes %parallel
  // step-40 from sequential step-6.
  //
  // At the top we do what we always do: tell the DoFHandler object to
  // distribute degrees of freedom. Since the triangulation we use here is
  // distributed, the DoFHandler object is smart enough to recognize that on
  // each processor it can only distribute degrees of freedom on cells it
  // owns; this is followed by an exchange step in which processors tell each
  // other about degrees of freedom on ghost cell. The result is a DoFHandler
  // that knows about the degrees of freedom on locally owned cells and ghost
  // cells (i.e. cells adjacent to locally owned cells) but nothing about
  // cells that are further away, consistent with the basic philosophy of
  // distributed computing that no processor can know everything.
  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    TimerOutput::Scope t(computing_timer, "setup");

    dof_handler.distribute_dofs (fe);

    // The next two lines extract some information we will need later on,
    // namely two index sets that provide information about which degrees of
    // freedom are owned by the current processor (this information will be
    // used to initialize solution and right hand side vectors, and the system
    // matrix, indicating which elements to store on the current processor and
    // which to expect to be stored somewhere else); and an index set that
    // indicates which degrees of freedom are locally relevant (i.e. live on
    // cells that the current processor owns or on the layer of ghost cells
    // around the locally owned cells; we need all of these degrees of
    // freedom, for example, to estimate the error on the local cells).
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);

    // Next, let us initialize the solution and right hand side vectors. As
    // mentioned above, the solution vector we seek does not only store
    // elements we own, but also ghost entries; on the other hand, the right
    // hand side vector only needs to have the entries the current processor
    // owns since all we will ever do is write into it, never read from it on
    // locally owned cells (of course the linear solvers will read from it,
    // but they do not care about the geometric location of degrees of
    // freedom).
    locally_relevant_solution.reinit (locally_owned_dofs,
                                      locally_relevant_dofs, mpi_communicator);
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);

    // The next step is to compute hanging node and boundary value
    // constraints, which we combine into a single object storing all
    // constraints.
    //
    // As with all other things in %parallel, the mantra must be that no
    // processor can store all information about the entire universe. As a
    // consequence, we need to tell the constraints object for which degrees
    // of freedom it can store constraints and for which it may not expect any
    // information to store. In our case, as explained in the @ref distributed
    // module, the degrees of freedom we need to care about on each processor
    // are the locally relevant ones, so we pass this to the
    // ConstraintMatrix::reinit function. As a side note, if you forget to
    // pass this argument, the ConstraintMatrix class will allocate an array
    // with length equal to the largest DoF index it has seen so far. For
    // processors with high MPI process number, this may be very large --
    // maybe on the order of billions. The program would then allocate more
    // memory than for likely all other operations combined for this single
    // array.
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    constraints.close ();

    // The last part of this function deals with initializing the matrix with
    // accompanying sparsity pattern. As in previous tutorial programs, we use
    // the CompressedSimpleSparsityPattern as an intermediate with which we
    // then initialize the PETSc matrix. To do so we have to tell the sparsity
    // pattern its size but as above there is no way the resulting object will
    // be able to store even a single pointer for each global degree of
    // freedom; the best we can hope for is that it stores information about
    // each locally relevant degree of freedom, i.e. all those that we may
    // ever touch in the process of assembling the matrix (the @ref
    // distributed_paper "distributed computing paper" has a long discussion
    // why one really needs the locally relevant, and not the small set of
    // locally active degrees of freedom in this context).
    //
    // So we tell the sparsity pattern its size and what DoFs to store
    // anything for and then ask DoFTools::make_sparsity_pattern to fill it
    // (this function ignores all cells that are not locally owned, mimicking
    // what we will do below in the assembly process). After this, we call a
    // function that exchanges entries in these sparsity pattern between
    // processors so that in the end each processor really knows about all the
    // entries that will exist in that part of the finite element matrix that
    // it will own. The final step is to initialize the matrix with the
    // sparsity pattern.
    CompressedSimpleSparsityPattern csp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, csp,
                                     constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp,
                                                dof_handler.n_locally_owned_dofs_per_processor(),
                                                mpi_communicator,
                                                locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs,
                          locally_owned_dofs,
                          csp,
                          mpi_communicator);
  }



  // @sect4{LaplaceProblem::assemble_system}

  // The function that then assembles the linear system is comparatively
  // boring, being almost exactly what we've seen before. The points to watch
  // out for are:
  // - Assembly must only loop over locally owned cells. There
  //   are multiple ways to test that; for example, we could compare a cell's
  //   subdomain_id against information from the triangulation as in
  //   <code>cell->subdomain_id() ==
  //   triangulation.locally_owned_subdomain()</code>, or skip all cells for
  //   which the condition <code>cell->is_ghost() ||
  //   cell->is_artificial()</code> is true. The simplest way, however, is to
  //   simply ask the cell whether it is owned by the local processor.
  // - Copying local contributions into the global matrix must include
  //   distributing constraints and boundary values. In other words, we can now
  //   (as we did in step-6) first copy every local contribution into the global
  //   matrix and only in a later step take care of hanging node constraints and
  //   boundary values. The reason is, as discussed in step-17, that PETSc does
  //   not provide access to arbitrary elements of the matrix once they have
  //   been assembled into it -- in parts because they may simple no longer
  //   reside on the current processor but have instead been shipped to a
  //   different machine.
  // - The way we compute the right hand side (given the
  //   formula stated in the introduction) may not be the most elegant but will
  //   do for a program whose focus lies somewhere entirely different.
  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    TimerOutput::Scope t(computing_timer, "assembly");

    const QGauss<dim>  quadrature_formula(3);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit (cell);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            {
              const double
              rhs_value
                = (fe_values.quadrature_point(q_point)[1]
                   >
                   0.5+0.25*std::sin(4.0 * numbers::PI *
                                     fe_values.quadrature_point(q_point)[0])
                   ? 1 : -1);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                         fe_values.shape_grad(j,q_point) *
                                         fe_values.JxW(q_point));

                  cell_rhs(i) += (rhs_value *
                                  fe_values.shape_value(i,q_point) *
                                  fe_values.JxW(q_point));
                }
            }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix,
                                                  cell_rhs,
                                                  local_dof_indices,
                                                  system_matrix,
                                                  system_rhs);
        }

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }



  // @sect4{LaplaceProblem::solve}

  // Even though solving linear systems on potentially tens of thousands of
  // processors is by far not a trivial job, the function that does this is --
  // at least at the outside -- relatively simple. Most of the parts you've
  // seen before. There are really only two things worth mentioning:
  // - Solvers and preconditioners are built on the deal.II wrappers of PETSc
  //   functionality. It is relatively well known that the primary bottleneck of
  //   massively %parallel linear solvers is not actually the communication
  //   between processors, but the fact that it is difficult to produce
  //   preconditioners that scale well to large numbers of processors. Over the
  //   second half of the first decade of the 21st century, it has become clear
  //   that algebraic multigrid (AMG) methods turn out to be extremely efficient
  //   in this context, and we will use one of them -- the BoomerAMG
  //   implementation of the Hypre package that can be interfaced to through
  //   PETSc -- for the current program. The rest of the solver itself is
  //   boilerplate and has been shown before. Since the linear system is
  //   symmetric and positive definite, we can use the CG method as the outer
  //   solver.
  // - Ultimately, we want a vector that stores not only the elements
  //   of the solution for degrees of freedom the current processor owns, but
  //   also all other locally relevant degrees of freedom. On the other hand,
  //   the solver itself needs a vector that is uniquely split between
  //   processors, without any overlap. We therefore create a vector at the
  //   beginning of this function that has these properties, use it to solve the
  //   linear system, and only assign it to the vector we want at the very
  //   end. This last step ensures that all ghost elements are also copied as
  //   necessary.
  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector
    completely_distributed_solution (locally_owned_dofs, mpi_communicator);

    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);

    LA::SolverCG solver(solver_control, mpi_communicator);
    LA::MPI::PreconditionAMG preconditioner;

    LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#else
    /* Trilinos defaults are good */
#endif
    preconditioner.initialize(system_matrix, data);

    solver.solve (system_matrix, completely_distributed_solution, system_rhs,
                  preconditioner);

    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;

    constraints.distribute (completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }



  // @sect4{LaplaceProblem::refine_grid}

  // The function that estimates the error and refines the grid is again
  // almost exactly like the one in step-6. The only difference is that the
  // function that flags cells to be refined is now in namespace
  // parallel::distributed::GridRefinement -- a namespace that has functions
  // that can communicate between all involved processors and determine global
  // thresholds to use in deciding which cells to refine and which to coarsen.
  //
  // Note that we didn't have to do anything special about the
  // KellyErrorEstimator class: we just give it a vector with as many elements
  // as the local triangulation has cells (locally owned cells, ghost cells,
  // and artificial ones), but it only fills those entries that correspond to
  // cells that are locally owned.
  template <int dim>
  void LaplaceProblem<dim>::refine_grid ()
  {
    TimerOutput::Scope t(computing_timer, "refine");

    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        locally_relevant_solution,
                                        estimated_error_per_cell);
    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (triangulation,
                                     estimated_error_per_cell,
                                     0.3, 0.03);
    triangulation.execute_coarsening_and_refinement ();
  }



  // @sect4{LaplaceProblem::output_results}

  // Compared to the corresponding function in step-6, the one here is a tad
  // more complicated. There are two reasons: the first one is that we do not
  // just want to output the solution but also for each cell which processor
  // owns it (i.e. which "subdomain" it is in). Secondly, as discussed at
  // length in step-17 and step-18, generating graphical data can be a
  // bottleneck in parallelizing. In step-18, we have moved this step out of
  // the actual computation but shifted it into a separate program that later
  // combined the output from various processors into a single file. But this
  // doesn't scale: if the number of processors is large, this may mean that
  // the step of combining data on a single processor later becomes the
  // longest running part of the program, or it may produce a file that's so
  // large that it can't be visualized any more. We here follow a more
  // sensible approach, namely creating individual files for each MPI process
  // and leaving it to the visualization program to make sense of that.
  //
  // To start, the top of the function looks like always. In addition to
  // attaching the solution vector (the one that has entries for all locally
  // relevant, not only the locally owned, elements), we attach a data vector
  // that stores, for each cell, the subdomain the cell belongs to. This is
  // slightly tricky, because of course not every processor knows about every
  // cell. The vector we attach therefore has an entry for every cell that the
  // current processor has in its mesh (locally owned ones, ghost cells, and
  // artificial cells), but the DataOut class will ignore all entries that
  // correspond to cells that are not owned by the current processor. As a
  // consequence, it doesn't actually matter what values we write into these
  // vector entries: we simply fill the entire vector with the number of the
  // current MPI process (i.e. the subdomain_id of the current process); this
  // correctly sets the values we care for, i.e. the entries that correspond
  // to locally owned cells, while providing the wrong value for all other
  // elements -- but these are then ignored anyway.
  template <int dim>
  void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (locally_relevant_solution, "u");

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches ();

    // The next step is to write this data to disk. We choose file names of
    // the form <code>solution-XX-PPPP.vtu</code> where <code>XX</code>
    // indicates the refinement cycle, <code>PPPP</code> refers to the
    // processor number (enough for up to 10,000 processors, though we hope
    // that nobody ever tries to generate this much data -- you would likely
    // overflow all file system quotas), and <code>.vtu</code> indicates the
    // XML-based Visualization Toolkit (VTK) file format.
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);

    // The last step is to write a "master record" that lists for the
    // visualization program the names of the various files that combined
    // represents the graphical data for the entire domain. The
    // DataOutBase::write_pvtu_record does this, and it needs a list of
    // filenames that we create first. Note that only one processor needs to
    // generate this file; we arbitrarily choose processor zero to take over
    // this job.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          filenames.push_back ("solution-" +
                               Utilities::int_to_string (cycle, 2) +
                               "." +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");

        std::ofstream master_output ((filename + ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }
  }



  // @sect4{LaplaceProblem::run}

  // The function that controls the overall behavior of the program is again
  // like the one in step-6. The minor difference are the use of
  // <code>pcout</code> instead of <code>std::cout</code> for output to the
  // console (see also step-17) and that we only generate graphical output if
  // at most 32 processors are involved. Without this limit, it would be just
  // too easy for people carelessly running this program without reading it
  // first to bring down the cluster interconnect and fill any file system
  // available :-)
  //
  // A functional difference to step-6 is the use of a square domain and that
  // we start with a slightly finer mesh (5 global refinement cycles) -- there
  // just isn't much of a point showing a massively %parallel program starting
  // on 4 cells (although admittedly the point is only slightly stronger
  // starting on 1024).
  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    const unsigned int n_cycles = 8;
    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation);
            triangulation.refine_global (5);
          }
        else
          refine_grid ();

        setup_system ();

        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells()
              << std::endl
              << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

        assemble_system ();
        solve ();

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results (cycle);
          }

        computing_timer.print_summary ();
        computing_timer.reset ();

        pcout << std::endl;
      }
  }
}



// @sect4{main()}

// The final function, <code>main()</code>, again has the same structure as in
// all other programs, in particular step-6. Like in the other programs that
// use PETSc, we have to initialize and finalize PETSc, which is done using the
// helper object MPI_InitFinalize.
//
// Note how we enclose the use the use of the LaplaceProblem class in a pair
// of braces. This makes sure that all member variables of the object are
// destroyed by the time we destroy the mpi_initialization object. Not doing
// this will lead to strange and hard to debug errors when
// <code>PetscFinalize</code> first deletes all PETSc vectors that are still
// around, and the destructor of the LaplaceProblem class then tries to delete
// them again.
int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step40;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      deallog.depth_console (0);

      {
        LaplaceProblem<2> laplace_problem_2d;
        laplace_problem_2d.run ();
      }
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
