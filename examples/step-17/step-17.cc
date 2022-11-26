/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2022 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004
 *         Wolfgang Bangerth, Texas A&M University, 2016
 */

// @sect3{Include files}

// First the usual assortment of header files we have already used in previous
// example programs:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// And here come the things that we need particularly for this example program
// and that weren't in step-8. First, we replace the standard output
// <code>std::cout</code> by a new stream <code>pcout</code> which is used in
// parallel computations for generating output only on one of the MPI
// processes.
#include <deal.II/base/conditional_ostream.h>
// We are going to query the number of processes and the number of the present
// process by calling the respective functions in the Utilities::MPI
// namespace.
#include <deal.II/base/mpi.h>
// Then, we are going to replace all linear algebra components that involve
// the (global) linear system by classes that wrap interfaces similar to our
// own linear algebra classes around what PETSc offers (PETSc is a library
// written in C, and deal.II comes with wrapper classes that provide the PETSc
// functionality with an interface that is similar to the interface we already
// had for our own linear algebra classes). In particular, we need vectors and
// matrices that are distributed across several
// @ref GlossMPIProcess "processes" in MPI programs (and
// simply map to sequential, local vectors and matrices if there is only a
// single process, i.e., if you are running on only one machine, and without
// MPI support):
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
// Then we also need interfaces for solvers and preconditioners that PETSc
// provides:
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
// And in addition, we need some algorithms for partitioning our meshes so
// that they can be efficiently distributed across an MPI network. The
// partitioning algorithm is implemented in the <code>GridTools</code>
// namespace, and we need an additional include file for a function in
// <code>DoFRenumbering</code> that allows to sort the indices associated with
// degrees of freedom so that they are numbered according to the subdomain
// they are associated with:
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

// And this is simply C++ again:
#include <fstream>
#include <iostream>

// The last step is as in all previous programs:
namespace Step17
{
  using namespace dealii;

  // @sect3{The <code>ElasticProblem</code> class template}

  // The first real part of the program is the declaration of the main
  // class.  As mentioned in the introduction, almost all of this has
  // been copied verbatim from step-8, so we only comment on the few
  // differences between the two tutorials.  There is one (cosmetic)
  // change in that we let <code>solve</code> return a value, namely
  // the number of iterations it took to converge, so that we can
  // output this to the screen at the appropriate place.
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem();
    void run();

  private:
    void         setup_system();
    void         assemble_system();
    unsigned int solve();
    void         refine_grid();
    void         output_results(const unsigned int cycle) const;

    // The first change is that we have to declare a variable that
    // indicates the @ref GlossMPICommunicator "MPI communicator" over
    // which we are supposed to distribute our computations.
    MPI_Comm mpi_communicator;

    // Then we have two variables that tell us where in the parallel
    // world we are. The first of the following variables,
    // <code>n_mpi_processes</code>, tells us how many MPI processes
    // there exist in total, while the second one,
    // <code>this_mpi_process</code>, indicates which is the number of
    // the present process within this space of processes (in MPI
    // language, this corresponds to the @ref GlossMPIRank "rank" of
    // the process). The latter will have a unique value for each
    // process between zero and (less than)
    // <code>n_mpi_processes</code>. If this program is run on a
    // single machine without MPI support, then their values are
    // <code>1</code> and <code>0</code>, respectively.
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;

    // Next up is a stream-like variable <code>pcout</code>. It is, in essence,
    // just something we use for convenience: in a parallel program,
    // if each process outputs status information, then there quickly
    // is a lot of clutter. Rather, we would want to only have one
    // @ref GlossMPIProcess "process" output everything once, for
    // example the one with @ref GlossMPIRank "rank" zero. At the same
    // time, it seems silly to prefix <i>every</i> place where we
    // create output with an <code>if (my_rank==0)</code> condition.
    //
    // To make this simpler, the ConditionalOStream class does exactly
    // this under the hood: it acts as if it were a stream, but only
    // forwards to a real, underlying stream if a flag is set. By
    // setting this condition to <code>this_mpi_process==0</code>
    // (where <code>this_mpi_process</code> corresponds to the rank of
    // an MPI process), we make sure that output is only generated
    // from the first process and that we don't get the same lines of
    // output over and over again, once per process. Thus, we can use
    // <code>pcout</code> everywhere and in every process, but on all
    // but one process nothing will ever happen to the information
    // that is piped into the object via
    // <code>operator&lt;&lt;</code>.
    ConditionalOStream pcout;

    // The remainder of the list of member variables is fundamentally the
    // same as in step-8. However, we change the declarations of matrix
    // and vector types to use parallel PETSc objects instead. Note that
    // we do not use a separate sparsity pattern, since PETSc manages this
    // internally as part of its matrix data structures.
    Triangulation<dim> triangulation;
    FESystem<dim>      fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> hanging_node_constraints;

    PETScWrappers::MPI::SparseMatrix system_matrix;

    PETScWrappers::MPI::Vector solution;
    PETScWrappers::MPI::Vector system_rhs;
  };


  // @sect3{Right hand side values}

  // The following is taken from step-8 without change:
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &  values) const override
    {
      AssertDimension(values.size(), dim);
      Assert(dim >= 2, ExcInternalError());

      Point<dim> point_1, point_2;
      point_1(0) = 0.5;
      point_2(0) = -0.5;

      if (((p - point_1).norm_square() < 0.2 * 0.2) ||
          ((p - point_2).norm_square() < 0.2 * 0.2))
        values(0) = 1;
      else
        values(0) = 0;

      if (p.square() < 0.2 * 0.2)
        values(1) = 1;
      else
        values(1) = 0;
    }

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  value_list) const override
    {
      const unsigned int n_points = points.size();

      AssertDimension(value_list.size(), n_points);

      for (unsigned int p = 0; p < n_points; ++p)
        RightHandSide<dim>::vector_value(points[p], value_list[p]);
    }
  };



  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{ElasticProblem::ElasticProblem}

  // The first step in the actual implementation is the constructor of
  // the main class. Apart from initializing the same member variables
  // that we already had in step-8, we here initialize the MPI
  // communicator variable we shall use with the global MPI
  // communicator linking all processes together (in more complex
  // applications, one could here use a communicator object that only
  // links a subset of all processes), and call the Utilities::MPI
  // helper functions to determine the number of processes and where
  // the present one fits into this picture. In addition, we make sure
  // that output is only generated by the (globally) first process.
  // We do so by passing the stream we want to output to
  // (<code>std::cout</code>) and a true/false flag as arguments where
  // the latter is determined by testing whether the process currently
  // executing the constructor call is the first in the MPI universe.
  template <int dim>
  ElasticProblem<dim>::ElasticProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
    , pcout(std::cout, (this_mpi_process == 0))
    , fe(FE_Q<dim>(1), dim)
    , dof_handler(triangulation)
  {}



  // @sect4{ElasticProblem::setup_system}

  // Next, the function in which we set up the various variables
  // for the global linear system to be solved needs to be implemented.
  //
  // However, before we proceed with this, there is one thing to do for a
  // parallel program: we need to determine which MPI process is
  // responsible for each of the cells. Splitting cells among
  // processes, commonly called "partitioning the mesh", is done by
  // assigning a @ref GlossSubdomainId "subdomain id" to each cell. We
  // do so by calling into the METIS library that does this in a very
  // efficient way, trying to minimize the number of nodes on the
  // interfaces between subdomains. Rather than trying to call METIS
  // directly, we do this by calling the
  // GridTools::partition_triangulation() function that does this at a
  // much higher level of programming.
  //
  // @note As mentioned in the introduction, we could avoid this manual
  //   partitioning step if we used the parallel::shared::Triangulation
  //   class for the triangulation object instead (as we do in step-18).
  //   That class does, in essence, everything a regular triangulation
  //   does, but it then also automatically partitions the mesh after
  //   every mesh creation or refinement operation.
  //
  // Following partitioning, we need to enumerate all degrees of
  // freedom as usual.  However, we would like to enumerate the
  // degrees of freedom in a way so that all degrees of freedom
  // associated with cells in subdomain zero (which resides on process
  // zero) come before all DoFs associated with cells on subdomain
  // one, before those on cells on process two, and so on. We need
  // this since we have to split the global vectors for right hand
  // side and solution, as well as the matrix into contiguous chunks
  // of rows that live on each of the processors, and we will want to
  // do this in a way that requires minimal communication. This
  // particular enumeration can be obtained by re-ordering degrees of
  // freedom indices using DoFRenumbering::subdomain_wise().
  //
  // The final step of this initial setup is that we get ourselves an
  // IndexSet that indicates the subset of the global number of unknowns
  // this
  // process is responsible for. (Note that a degree of freedom is not
  // necessarily owned by the process that owns a cell just because
  // the degree of freedom lives on this cell: some degrees of freedom
  // live on interfaces between subdomains, and are consequently only owned by
  // one of the processes adjacent to this interface.)
  //
  // Before we move on, let us recall a fact already discussed in the
  // introduction: The triangulation we use here is replicated across
  // all processes, and each process has a complete copy of the entire
  // triangulation, with all cells. Partitioning only provides a way
  // to identify which cells out of all each process "owns", but it
  // knows everything about all of them. Likewise, the DoFHandler
  // object knows everything about every cell, in particular the
  // degrees of freedom that live on each cell, whether it is one that
  // the current process owns or not. This can not scale to large
  // problems because eventually just storing the entire mesh, and
  // everything that is associated with it, on every process will
  // become infeasible if the problem is large enough. On the other
  // hand, if we split the triangulation into parts so that every
  // process stores only those cells it "owns" but nothing else (or,
  // at least a sufficiently small fraction of everything else), then
  // we can solve large problems if only we throw a large enough
  // number of MPI processes at them. This is what we are going to in
  // step-40, for example, using the
  // parallel::distributed::Triangulation class.  On the other hand,
  // most of the rest of what we demonstrate in the current program
  // will actually continue to work whether we have the entire
  // triangulation available, or only a piece of it.
  template <int dim>
  void ElasticProblem<dim>::setup_system()
  {
    GridTools::partition_triangulation(n_mpi_processes, triangulation);

    dof_handler.distribute_dofs(fe);
    DoFRenumbering::subdomain_wise(dof_handler);

    // We need to initialize the objects denoting hanging node constraints for
    // the present grid. As with the triangulation and DoFHandler objects, we
    // will simply store <i>all</i> constraints on each process; again, this
    // will not scale, but we show in step-40 how one can work around this by
    // only storing on each MPI process the constraints for degrees of freedom
    // that actually matter on this particular process.
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    // Now we create the sparsity pattern for the system matrix. Note that we
    // again compute and store all entries and not only the ones relevant
    // to this process (see step-18 or step-40 for a more efficient way to
    // handle this).
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    hanging_node_constraints,
                                    false);

    // Now we determine the set of locally owned DoFs and use that to
    // initialize parallel vectors and matrix. Since the matrix and vectors
    // need to work in parallel, we have to pass them an MPI communication
    // object, as well as information about the partitioning contained in the
    // IndexSet @p locally_owned_dofs.  The IndexSet contains information
    // about the global size (the <i>total</i> number of degrees of freedom)
    // and also what subset of rows is to be stored locally.  Note that the
    // system matrix needs that partitioning information for the rows and
    // columns. For square matrices, as it is the case here, the columns
    // should be partitioned in the same way as the rows, but in the case of
    // rectangular matrices one has to partition the columns in the same way
    // as vectors are partitioned with which the matrix is multiplied, while
    // rows have to partitioned in the same way as destination vectors of
    // matrix-vector multiplications:
    const std::vector<IndexSet> locally_owned_dofs_per_proc =
      DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
    const IndexSet locally_owned_dofs =
      locally_owned_dofs_per_proc[this_mpi_process];

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);

    solution.reinit(locally_owned_dofs, mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  }



  // @sect4{ElasticProblem::assemble_system}

  // We now assemble the matrix and right hand side of the
  // problem. There are some things worth mentioning before we go into
  // detail. First, we will be assembling the system in parallel,
  // i.e., each process will be responsible for assembling on cells
  // that belong to this particular process. Note that the degrees of
  // freedom are split in a way such that all DoFs in the interior of
  // cells and between cells belonging to the same subdomain belong to
  // the process that <code>owns</code> the cell. However, even then
  // we sometimes need to assemble on a cell with a neighbor that
  // belongs to a different process, and in these cases when we add up
  // the local contributions into the global matrix or right hand side
  // vector, we have to transfer these entries to the process that
  // owns these elements. Fortunately, we don't have to do this by
  // hand: PETSc does all this for us by caching these elements
  // locally, and sending them to the other processes as necessary
  // when we call the <code>compress()</code> functions on the matrix
  // and vector at the end of this function.
  //
  // The second point is that once we have handed over matrix and
  // vector contributions to PETSc, it is a) hard, and b) very
  // inefficient to get them back for modifications. This is not only
  // the fault of PETSc, it is also a consequence of the distributed
  // nature of this program: if an entry resides on another processor,
  // then it is necessarily expensive to get it. The consequence of
  // this is that we should not try to first assemble the matrix and
  // right hand side as if there were no hanging node constraints and
  // boundary values, and then eliminate these in a second step
  // (using, for example, AffineConstraints::condense()). Rather, we
  // should try to eliminate hanging node constraints before handing
  // these entries over to PETSc. This is easy: instead of copying
  // elements by hand into the global matrix (as we do in step-4), we
  // use the AffineConstraints::distribute_local_to_global() functions
  // to take care of hanging nodes at the same time. We also already
  // did this in step-6. The second step, elimination of boundary
  // nodes, could also be done this way by putting the boundary values
  // into the same AffineConstraints object as hanging nodes (see the
  // way it is done in step-6, for example); however, it is not
  // strictly necessary to do this here because eliminating boundary
  // values can be done with only the data stored on each process
  // itself, and consequently we use the approach used before in
  // step-4, i.e., via MatrixTools::apply_boundary_values().
  //
  // All of this said, here is the actual implementation starting with
  // the general setup of helper variables.  (Note that we still use
  // the deal.II full matrix and vector types for the local systems as
  // these are small and need not be shared across processes.)
  template <int dim>
  void ElasticProblem<dim>::assemble_system()
  {
    QGauss<dim>   quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);

    Functions::ConstantFunction<dim> lambda(1.), mu(1.);

    RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim));


    // The next thing is the loop over all elements. Note that we do
    // not have to do <i>all</i> the work on every process: our job
    // here is only to assemble the system on cells that actually
    // belong to this MPI process, all other cells will be taken care
    // of by other processes. This is what the if-clause immediately
    // after the for-loop takes care of: it queries the subdomain
    // identifier of each cell, which is a number associated with each
    // cell that tells us about the owner process. In more generality,
    // the subdomain id is used to split a domain into several parts
    // (we do this above, at the beginning of
    // <code>setup_system()</code>), and which allows to identify
    // which subdomain a cell is living on. In this application, we
    // have each process handle exactly one subdomain, so we identify
    // the terms <code>subdomain</code> and <code>MPI process</code>.
    //
    // Apart from this, assembling the local system is relatively uneventful
    // if you have understood how this is done in step-8. As mentioned above,
    // distributing local contributions into the global matrix
    // and right hand sides also takes care of hanging node constraints in the
    // same way as is done in step-6.
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->subdomain_id() == this_mpi_process)
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);

          lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
          mu.value_list(fe_values.get_quadrature_points(), mu_values);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;

              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const unsigned int component_j =
                    fe.system_to_component_index(j).first;

                  for (unsigned int q_point = 0; q_point < n_q_points;
                       ++q_point)
                    {
                      cell_matrix(i, j) +=
                        ((fe_values.shape_grad(i, q_point)[component_i] *
                          fe_values.shape_grad(j, q_point)[component_j] *
                          lambda_values[q_point]) +
                         (fe_values.shape_grad(i, q_point)[component_j] *
                          fe_values.shape_grad(j, q_point)[component_i] *
                          mu_values[q_point]) +
                         ((component_i == component_j) ?
                            (fe_values.shape_grad(i, q_point) *
                             fe_values.shape_grad(j, q_point) *
                             mu_values[q_point]) :
                            0)) *
                        fe_values.JxW(q_point);
                    }
                }
            }

          right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                            rhs_values);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                cell_rhs(i) += fe_values.shape_value(i, q_point) *
                               rhs_values[q_point](component_i) *
                               fe_values.JxW(q_point);
            }

          cell->get_dof_indices(local_dof_indices);
          hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                                              cell_rhs,
                                                              local_dof_indices,
                                                              system_matrix,
                                                              system_rhs);
        }

    // The next step is to "compress" the vector and the system matrix. This
    // means that each process sends the additions that were made to those
    // entries of the matrix and vector that the process did not own itself to
    // the process that owns them. After receiving these additions from other
    // processes, each process then adds them to the values it already
    // has. These additions are combining the integral contributions of shape
    // functions living on several cells just as in a serial computation, with
    // the difference that the cells are assigned to different processes.
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    // The global matrix and right hand side vectors have now been
    // formed. We still have to apply boundary values, in the same way as we
    // did, for example, in step-3, step-4, and a number of other programs.
    //
    // The last argument to the call to
    // MatrixTools::apply_boundary_values() below allows for some
    // optimizations. It controls whether we should also delete
    // entries (i.e., set them to zero) in the matrix columns
    // corresponding to boundary nodes, or to keep them (and passing
    // <code>true</code> means: yes, do eliminate the columns). If we
    // do eliminate columns, then the resulting matrix will be
    // symmetric again if it was before; if we don't, then it
    // won't. The solution of the resulting system should be the same,
    // though. The only reason why we may want to make the system
    // symmetric again is that we would like to use the CG method,
    // which only works with symmetric matrices. The reason why we may
    // <i>not</i> want to make the matrix symmetric is because this
    // would require us to write into column entries that actually
    // reside on other processes, i.e., it involves communicating
    // data.
    //
    // Experience tells us that CG also works (and works almost as
    // well) if we don't remove the columns associated with boundary
    // nodes, which can be explained by the special structure of this
    // particular non-symmetry. To avoid the expense of communication,
    // we therefore do not eliminate the entries in the affected
    // columns.
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             boundary_values);
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, false);
  }



  // @sect4{ElasticProblem::solve}

  // Having assembled the linear system, we next need to solve
  // it. PETSc offers a variety of sequential and parallel solvers,
  // for which we have written wrappers that have almost the same
  // interface as is used for the deal.II solvers used in all previous
  // example programs. The following code should therefore look rather
  // familiar.
  //
  // At the top of the function, we set up a convergence monitor, and
  // assign it the accuracy to which we would like to solve the linear
  // system. Next, we create an actual solver object using PETSc's CG
  // solver which also works with parallel (distributed) vectors and
  // matrices. And finally a preconditioner; we choose to use a block
  // Jacobi preconditioner which works by computing an incomplete LU
  // decomposition on each diagonal block of the matrix.  (In other
  // words, each MPI process computes an ILU from the rows it stores
  // by throwing away columns that correspond to row indices not
  // stored locally; this yields a square matrix block from which we
  // can compute an ILU. That means that if you run the program with
  // only one process, then you will use an ILU(0) as a
  // preconditioner, while if it is run on many processes, then we
  // will have a number of blocks on the diagonal and the
  // preconditioner is the ILU(0) of each of these blocks. In the
  // extreme case of one degree of freedom per processor, this
  // preconditioner is then simply a Jacobi preconditioner since the
  // diagonal matrix blocks consist of only a single entry. Such a
  // preconditioner is relatively easy to compute because it does not
  // require any kind of communication between processors, but it is
  // in general not very efficient for large numbers of processors.)
  //
  // Following this kind of setup, we then solve the linear system:
  template <int dim>
  unsigned int ElasticProblem<dim>::solve()
  {
    SolverControl solver_control(solution.size(), 1e-8 * system_rhs.l2_norm());
    PETScWrappers::SolverCG cg(solver_control);

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    // The next step is to distribute hanging node constraints. This is a
    // little tricky, since to fill in the value of a constrained node you
    // need access to the values of the nodes to which it is constrained (for
    // example, for a Q1 element in 2d, we need access to the two nodes on the
    // big side of a hanging node face, to compute the value of the
    // constrained node in the middle).
    //
    // The problem is that we have built our vectors (in
    // <code>setup_system()</code>) in such a way that every process
    // is responsible for storing only those elements of the solution
    // vector that correspond to the degrees of freedom this process
    // "owns". There are, however, cases where in order to compute the
    // value of the vector entry for a constrained degree of freedom
    // on one process, we need to access vector entries that are
    // stored on other processes.  PETSc (and, for that matter, the
    // MPI model on which it is built) does not allow to simply query
    // vector entries stored on other processes, so what we do here is
    // to get a copy of the "distributed" vector where we store all
    // elements locally. This is simple, since the deal.II wrappers
    // have a conversion constructor for the deal.II Vector
    // class. (This conversion of course requires communication, but
    // in essence every process only needs to send its data to every
    // other process once in bulk, rather than having to respond to
    // queries for individual elements):
    Vector<double> localized_solution(solution);

    // Of course, as in previous discussions, it is clear that such a
    // step cannot scale very far if you wanted to solve large
    // problems on large numbers of processes, because every process
    // now stores <i>all elements</i> of the solution vector. (We will
    // show how to do this better in step-40.)  On the other hand,
    // distributing hanging node constraints is simple on this local
    // copy, using the usual function
    // AffineConstraints::distributed(). In particular, we can compute
    // the values of <i>all</i> constrained degrees of freedom,
    // whether the current process owns them or not:
    hanging_node_constraints.distribute(localized_solution);

    // Then transfer everything back into the global vector. The
    // following operation copies those elements of the localized
    // solution that we store locally in the distributed solution, and
    // does not touch the other ones. Since we do the same operation
    // on all processors, we end up with a distributed vector (i.e., a
    // vector that on every process only stores the vector entries
    // corresponding to degrees of freedom that are owned by this
    // process) that has all the constrained nodes fixed.
    //
    // We end the function by returning the number of iterations it
    // took to converge, to allow for some output.
    solution = localized_solution;

    return solver_control.last_step();
  }


  // @sect4{ElasticProblem::refine_grid}

  // Using some kind of refinement indicator, the mesh can be
  // refined. The problem is basically the same as with distributing
  // hanging node constraints: in order to compute the error indicator
  // (even if we were just interested in the indicator on the cells
  // the current process owns), we need access to more elements of the
  // solution vector than just those the current processor stores. To
  // make this happen, we do essentially what we did in
  // <code>solve()</code> already, namely get a <i>complete</i> copy
  // of the solution vector onto every process, and use that to
  // compute. This is in itself expensive as explained above and it
  // is particular unnecessary since we had just created and then
  // destroyed such a vector in <code>solve()</code>, but efficiency
  // is not the point of this program and so let us opt for a design
  // in which every function is as self-contained as possible.
  //
  // Once we have such a "localized" vector that contains <i>all</i>
  // elements of the solution vector, we can compute the indicators
  // for the cells that belong to the present process. In fact, we
  // could of course compute <i>all</i> refinement indicators since
  // our Triangulation and DoFHandler objects store information about
  // all cells, and since we have a complete copy of the solution
  // vector. But in the interest in showing how to operate in
  // %parallel, let us demonstrate how one would operate if one were
  // to only compute <i>some</i> error indicators and then exchange
  // the remaining ones with the other processes. (Ultimately, each
  // process needs a complete set of refinement indicators because
  // every process needs to refine their mesh, and needs to refine it
  // in exactly the same way as all of the other processes.)
  //
  // So, to do all of this, we need to:
  // - First, get a local copy of the distributed solution vector.
  // - Second, create a vector to store the refinement indicators.
  // - Third, let the KellyErrorEstimator compute refinement
  //   indicators for all cells belonging to the present
  //   subdomain/process. The last argument of the call indicates
  //   which subdomain we are interested in. The three arguments
  //   before it are various other default arguments that one usually
  //   does not need (and does not state values for, but rather uses the
  //   defaults), but which we have to state here explicitly since we
  //   want to modify the value of a following argument (i.e., the one
  //   indicating the subdomain).
  template <int dim>
  void ElasticProblem<dim>::refine_grid()
  {
    const Vector<double> localized_solution(solution);

    Vector<float> local_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       localized_solution,
                                       local_error_per_cell,
                                       ComponentMask(),
                                       nullptr,
                                       MultithreadInfo::n_threads(),
                                       this_mpi_process);

    // Now all processes have computed error indicators for their own
    // cells and stored them in the respective elements of the
    // <code>local_error_per_cell</code> vector. The elements of this
    // vector for cells not owned by the present process are
    // zero. However, since all processes have a copy of the entire
    // triangulation and need to keep these copies in sync, they need
    // the values of refinement indicators for all cells of the
    // triangulation. Thus, we need to distribute our results. We do
    // this by creating a distributed vector where each process has
    // its share and sets the elements it has computed. Consequently,
    // when you view this vector as one that lives across all
    // processes, then every element of this vector has been set
    // once. We can then assign this parallel vector to a local,
    // non-parallel vector on each process, making <i>all</i> error
    // indicators available on every process.
    //
    // So in the first step, we need to set up a parallel vector. For
    // simplicity, every process will own a chunk with as many
    // elements as this process owns cells, so that the first chunk of
    // elements is stored with process zero, the next chunk with
    // process one, and so on. It is important to remark, however,
    // that these elements are not necessarily the ones we will write
    // to. This is a consequence of the order in which cells are arranged,
    // i.e., the order in which the elements of the vector correspond
    // to cells is not ordered according to the subdomain these cells
    // belong to. In other words, if on this process we compute
    // indicators for cells of a certain subdomain, we may write the
    // results to more or less random elements of the distributed
    // vector; in particular, they may not necessarily lie within the
    // chunk of vector we own on the present process. They will
    // subsequently have to be copied into another process' memory
    // space, an operation that PETSc does for us when we call the
    // <code>compress()</code> function. This inefficiency could be
    // avoided with some more code, but we refrain from it since it is
    // not a major factor in the program's total runtime.
    //
    // So here is how we do it: count how many cells belong to this
    // process, set up a distributed vector with that many elements to
    // be stored locally, copy over the elements we computed
    // locally, and finally compress the result. In fact, we really only copy
    // the elements that are nonzero, so we may miss a few that we
    // computed to zero, but this won't hurt since the original values
    // of the vector are zero anyway.
    const unsigned int n_local_cells =
      GridTools::count_cells_with_subdomain_association(triangulation,
                                                        this_mpi_process);
    PETScWrappers::MPI::Vector distributed_all_errors(
      mpi_communicator, triangulation.n_active_cells(), n_local_cells);

    for (unsigned int i = 0; i < local_error_per_cell.size(); ++i)
      if (local_error_per_cell(i) != 0)
        distributed_all_errors(i) = local_error_per_cell(i);
    distributed_all_errors.compress(VectorOperation::insert);


    // So now we have this distributed vector that contains the
    // refinement indicators for all cells. To use it, we need to
    // obtain a local copy and then use it to mark cells for
    // refinement or coarsening, and actually do the refinement and
    // coarsening. It is important to recognize that <i>every</i>
    // process does this to its own copy of the triangulation, and
    // does it in exactly the same way.
    const Vector<float> localized_all_errors(distributed_all_errors);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    localized_all_errors,
                                                    0.3,
                                                    0.03);
    triangulation.execute_coarsening_and_refinement();
  }


  // @sect4{ElasticProblem::output_results}

  // The final function of significant interest is the one that
  // creates graphical output. This works the same way as in step-8,
  // with two small differences. Before discussing these, let us state
  // the general philosophy this function will work: we intend for all
  // of the data to be generated on a single process, and subsequently
  // written to a file. This is, as many other parts of this program
  // already discussed, not something that will scale. Previously, we
  // had argued that we will get into trouble with triangulations,
  // DoFHandlers, and copies of the solution vector where every
  // process has to store all of the data, and that there will come to
  // be a point where each process simply doesn't have enough memory
  // to store that much data. Here, the situation is different: it's
  // not only the memory, but also the run time that's a problem. If
  // one process is responsible for processing <i>all</i> of the data
  // while all of the other processes do nothing, then this one
  // function will eventually come to dominate the overall run time of
  // the program.  In particular, the time this function takes is
  // going to be proportional to the overall size of the problem
  // (counted in the number of cells, or the number of degrees of
  // freedom), independent of the number of processes we throw at it.
  //
  // Such situations need to be avoided, and we will show in step-18
  // and step-40 how to address this issue. For the current problem,
  // the solution is to have each process generate output data only
  // for its own local cells, and write them to separate files, one
  // file per process. This is how step-18 operates. Alternatively,
  // one could simply leave everything in a set of independent files
  // and let the visualization software read all of them (possibly
  // also using multiple processors) and create a single visualization
  // out of all of them; this is the path step-40, step-32, and all
  // other parallel programs developed later on take.
  //
  // More specifically for the current function, all processes call
  // this function, but not all of them need to do the work associated
  // with generating output. In fact, they shouldn't, since we would
  // try to write to the same file multiple times at once. So we let
  // only the first process do this, and all the other ones idle
  // around during this time (or start their work for the next
  // iteration, or simply yield their CPUs to other jobs that happen
  // to run at the same time). The second thing is that we not only
  // output the solution vector, but also a vector that indicates
  // which subdomain each cell belongs to. This will make for some
  // nice pictures of partitioned domains.
  //
  // To implement this, process zero needs a complete set of solution
  // components in a local vector. Just as with the previous function,
  // the efficient way to do this would be to re-use the vector
  // already created in the <code>solve()</code> function, but to keep
  // things more self-contained, we simply re-create one here from the
  // distributed solution vector.
  //
  // An important thing to realize is that we do this localization operation
  // on all processes, not only the one that actually needs the data. This
  // can't be avoided, however, with the simplified communication model of MPI
  // we use for vectors in this tutorial program: MPI does not have a way to
  // query data on another process, both sides have to initiate a
  // communication at the same time. So even though most of the processes do
  // not need the localized solution, we have to place the statement
  // converting the distributed into a localized vector so that all processes
  // execute it.
  //
  // (Part of this work could in fact be avoided. What we do is
  // send the local parts of all processes to all other processes. What we
  // would really need to do is to initiate an operation on all processes
  // where each process simply sends its local chunk of data to process
  // zero, since this is the only one that actually needs it, i.e., we need
  // something like a gather operation. PETSc can do this, but for
  // simplicity's sake we don't attempt to make use of this here. We don't,
  // since what we do is not very expensive in the grand scheme of things:
  // it is one vector communication among all processes, which has to be
  // compared to the number of communications we have to do when solving the
  // linear system, setting up the block-ILU for the preconditioner, and
  // other operations.)
  template <int dim>
  void ElasticProblem<dim>::output_results(const unsigned int cycle) const
  {
    const Vector<double> localized_solution(solution);

    // This being done, process zero goes ahead with setting up the
    // output file as in step-8, and attaching the (localized)
    // solution vector to the output object.
    if (this_mpi_process == 0)
      {
        std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);

        std::vector<std::string> solution_names;
        switch (dim)
          {
            case 1:
              solution_names.emplace_back("displacement");
              break;
            case 2:
              solution_names.emplace_back("x_displacement");
              solution_names.emplace_back("y_displacement");
              break;
            case 3:
              solution_names.emplace_back("x_displacement");
              solution_names.emplace_back("y_displacement");
              solution_names.emplace_back("z_displacement");
              break;
            default:
              Assert(false, ExcInternalError());
          }

        data_out.add_data_vector(localized_solution, solution_names);

        // The only other thing we do here is that we also output one
        // value per cell indicating which subdomain (i.e., MPI
        // process) it belongs to. This requires some conversion work,
        // since the data the library provides us with is not the one
        // the output class expects, but this is not difficult. First,
        // set up a vector of integers, one per cell, that is then
        // filled by the subdomain id of each cell.
        //
        // The elements of this vector are then converted to a
        // floating point vector in a second step, and this vector is
        // added to the DataOut object, which then goes off creating
        // output in VTK format:
        std::vector<unsigned int> partition_int(triangulation.n_active_cells());
        GridTools::get_subdomain_association(triangulation, partition_int);

        const Vector<double> partitioning(partition_int.begin(),
                                          partition_int.end());

        data_out.add_data_vector(partitioning, "partitioning");

        data_out.build_patches();
        data_out.write_vtk(output);
      }
  }


  // @sect4{ElasticProblem::run}

  // Lastly, here is the driver function. It is almost completely
  // unchanged from step-8, with the exception that we replace
  // <code>std::cout</code> by the <code>pcout</code> stream. Apart
  // from this, the only other cosmetic change is that we output how
  // many degrees of freedom there are per process, and how many
  // iterations it took for the linear solver to converge:
  template <int dim>
  void ElasticProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 10; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation, -1, 1);
            triangulation.refine_global(3);
          }
        else
          refine_grid();

        pcout << "   Number of active cells:       "
              << triangulation.n_active_cells() << std::endl;

        setup_system();

        pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << " (by partition:";
        for (unsigned int p = 0; p < n_mpi_processes; ++p)
          pcout << (p == 0 ? ' ' : '+')
                << (DoFTools::count_dofs_with_subdomain_association(dof_handler,
                                                                    p));
        pcout << ')' << std::endl;

        assemble_system();
        const unsigned int n_iterations = solve();

        pcout << "   Solver converged in " << n_iterations << " iterations."
              << std::endl;

        output_results(cycle);
      }
  }
} // namespace Step17


// @sect3{The <code>main</code> function}

// The <code>main()</code> works the same way as most of the main
// functions in the other example programs, i.e., it delegates work to
// the <code>run</code> function of a managing object, and only wraps
// everything into some code to catch exceptions:
int main(int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step17;

      // Here is the only real difference: MPI and PETSc both require that we
      // initialize these libraries at the beginning of the program, and
      // un-initialize them at the end. The class MPI_InitFinalize takes care
      // of all of that. The trailing argument `1` means that we do want to
      // run each MPI process with a single thread, a prerequisite with the
      // PETSc parallel linear algebra.
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      ElasticProblem<2> elastic_problem;
      elastic_problem.run();
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
