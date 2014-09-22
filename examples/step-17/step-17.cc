/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004
 */


// First the usual assortment of header files we have already used in previous
// example programs:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
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
// %parallel computations for generating output only on one of the MPI
// processes.
#include <deal.II/base/conditional_ostream.h>
// We are going to query the number of processes and the number of the present
// process by calling the respective functions in the Utilities::MPI
// namespace.
#include <deal.II/base/utilities.h>
// Then, we are going to replace all linear algebra components that involve
// the (global) linear system by classes that wrap interfaces similar to our
// own linear algebra classes around what PETSc offers (PETSc is a library
// written in C, and deal.II comes with wrapper classes that provide the PETSc
// functionality with an interface that is similar to the interface we already
// had for our own linear algebra classes). In particular, we need vectors and
// matrices that are distributed across several processes in MPI programs (and
// simply map to sequential, local vectors and matrices if there is only a
// single process, i.e. if you are running on only one machine, and without
// MPI support):
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
// Then we also need interfaces for solvers and preconditioners that PETSc
// provides:
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
// And in addition, we need some algorithms for partitioning our meshes so
// that they can be efficiently distributed across an MPI network. The
// partitioning algorithm is implemented in the <code>GridTools</code> class,
// and we need an additional include file for a function in
// <code>DoFRenumbering</code> that allows to sort the indices associated with
// degrees of freedom so that they are numbered according to the subdomain
// they are associated with:
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

// And this is simply C++ again:
#include <fstream>
#include <iostream>
#include <sstream>

// The last step is as in all previous programs:
namespace Step17
{
  using namespace dealii;

  // Now, here comes the declaration of the main class and of various other
  // things below it. As mentioned in the introduction, almost all of this has
  // been copied verbatim from step-8, so we only comment on the few things
  // that are different. There is one (cosmetic) change in that we let
  // <code>solve</code> return a value, namely the number of iterations it
  // took to converge, so that we can output this to the screen at the
  // appropriate place. In addition, we introduce a stream-like variable
  // <code>pcout</code>, explained below:
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem ();
    ~ElasticProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    unsigned int solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    // The first variable is basically only for convenience: in %parallel
    // program, if each process outputs status information, then there quickly
    // is a lot of clutter. Rather, we would want to only have one process
    // output everything once, for example the one with process number
    // zero. <code>ConditionalOStream</code> does exactly this: it acts as if
    // it were a stream, but only forwards to a real, underlying stream if a
    // flag is set. By setting this condition to
    // <code>this_mpi_process==0</code>, we make sure that output is only
    // generated from the first process and that we don't get the same lines
    // of output over and over again, once per process.
    //
    // With this simple trick, we make sure that we don't have to guard each
    // and every write to <code>std::cout</code> by a prefixed
    // <code>if(this_mpi_process==0)</code>.
    ConditionalOStream pcout;

    // The next few variables are taken verbatim from step-8:
    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    ConstraintMatrix     hanging_node_constraints;

    // In step-8, this would have been the place where we would have declared
    // the member variables for the sparsity pattern, the system matrix, right
    // hand, and solution vector. We change these declarations to use
    // %parallel PETSc objects instead (note that the fact that we use the
    // %parallel versions is denoted the fact that we use the classes from the
    // <code>PETScWrappers::MPI</code> namespace; sequential versions of these
    // classes are in the <code>PETScWrappers</code> namespace, i.e. without
    // the <code>MPI</code> part). Note also that we do not use a separate
    // sparsity pattern, since PETSc manages that as part of its matrix data
    // structures.
    PETScWrappers::MPI::SparseMatrix system_matrix;

    PETScWrappers::MPI::Vector       solution;
    PETScWrappers::MPI::Vector       system_rhs;

    // The next change is that we have to declare a variable that indicates
    // the MPI communicator over which we are supposed to distribute our
    // computations. Note that if this is a sequential job without support by
    // MPI, then PETSc provides some dummy type for <code>MPI_Comm</code>, so
    // we do not have to care here whether the job is really a %parallel one:
    MPI_Comm mpi_communicator;

    // Then we have two variables that tell us where in the %parallel world we
    // are. The first of the following variables, <code>n_mpi_processes</code>
    // tells us how many MPI processes there exist in total, while the second
    // one, <code>this_mpi_process</code>, indicates which is the number of
    // the present process within this space of processes. The latter variable
    // will have a unique value for each process between zero and (less than)
    // <code>n_mpi_processes</code>. If this program is run on a single
    // machine without MPI support, then their values are <code>1</code> and
    // <code>0</code>, respectively.
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;
  };


  // The following is again taken from step-8 without change:
  template <int dim>
  class RightHandSide :  public Function<dim>
  {
  public:
    RightHandSide ();

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &value_list) const;
  };


  template <int dim>
  RightHandSide<dim>::RightHandSide () :
    Function<dim> (dim)
  {}


  template <int dim>
  inline
  void RightHandSide<dim>::vector_value (const Point<dim> &p,
                                         Vector<double>   &values) const
  {
    Assert (values.size() == dim,
            ExcDimensionMismatch (values.size(), dim));
    Assert (dim >= 2, ExcInternalError());

    Point<dim> point_1, point_2;
    point_1(0) = 0.5;
    point_2(0) = -0.5;

    if (((p-point_1).square() < 0.2*0.2) ||
        ((p-point_2).square() < 0.2*0.2))
      values(0) = 1;
    else
      values(0) = 0;

    if (p.square() < 0.2*0.2)
      values(1) = 1;
    else
      values(1) = 0;
  }



  template <int dim>
  void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                              std::vector<Vector<double> >   &value_list) const
  {
    const unsigned int n_points = points.size();

    Assert (value_list.size() == n_points,
            ExcDimensionMismatch (value_list.size(), n_points));

    for (unsigned int p=0; p<n_points; ++p)
      RightHandSide<dim>::vector_value (points[p],
                                        value_list[p]);
  }


  // The first step in the actual implementation of things is the constructor
  // of the main class. Apart from initializing the same member variables that
  // we already had in step-8, we here initialize the MPI communicator
  // variable we shall use with the global MPI communicator linking all
  // processes together (in more complex applications, one could here use a
  // communicator object that only links a subset of all processes), and call
  // the Utilities helper functions to determine the number of processes and
  // where the present one fits into this picture. In addition, we make sure
  // that output is only generated by the (globally) first process. As
  // <code>this_mpi_process</code> is determined after creation of pcout, we
  // cannot set the condition through the constructor, i.e. by
  // <code>pcout(std::cout, this_mpi_process==0)</code>, but set the
  // condition separately.
  template <int dim>
  ElasticProblem<dim>::ElasticProblem ()
    :
    pcout (std::cout),
    dof_handler (triangulation),
    fe (FE_Q<dim>(1), dim),
    mpi_communicator (MPI_COMM_WORLD),
    n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
    this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
  {
    pcout.set_condition(this_mpi_process == 0);
  }



  template <int dim>
  ElasticProblem<dim>::~ElasticProblem ()
  {
    dof_handler.clear ();
  }


  // The second step is the function in which we set up the various variables
  // for the global linear system to be solved.
  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    // Before we even start out setting up the system, there is one thing to
    // do for a %parallel program: we need to assign cells to each of the
    // processes. We do this by splitting (<code>partitioning</code>) the mesh
    // cells into as many chunks (<code>subdomains</code>) as there are
    // processes in this MPI job (if this is a sequential job, then there is
    // only one job and all cells will get a zero as subdomain
    // indicator). This is done using an interface to the METIS library that
    // does this in a very efficient way, trying to minimize the number of
    // nodes on the interfaces between subdomains. All this is hidden behind
    // the following call to a deal.II library function:
    GridTools::partition_triangulation (n_mpi_processes, triangulation);

    // As for the linear system: First, we need to generate an enumeration for
    // the degrees of freedom in our problem. Further below, we will show how
    // we assign each cell to one of the MPI processes before we even get
    // here. What we then need to do is to enumerate the degrees of freedom in
    // a way so that all degrees of freedom associated with cells in subdomain
    // zero (which resides on process zero) come before all DoFs associated
    // with cells on subdomain one, before those on cells on process two, and
    // so on. We need this since we have to split the global vectors for right
    // hand side and solution, as well as the matrix into contiguous chunks of
    // rows that live on each of the processors, and we will want to do this
    // in a way that requires minimal communication. This is done using the
    // following two functions, which first generates an initial ordering of
    // all degrees of freedom, and then re-sort them according to above
    // criterion:
    dof_handler.distribute_dofs (fe);
    DoFRenumbering::subdomain_wise (dof_handler);

    // While we're at it, let us also count how many degrees of freedom there
    // exist on the present process:
    const types::global_dof_index n_local_dofs
      = DoFTools::count_dofs_with_subdomain_association (dof_handler,
                                                         this_mpi_process);

    // Then we initialize the system matrix, solution, and right hand side
    // vectors. Since they all need to work in %parallel, we have to pass them
    // an MPI communication object, as well as their global sizes (both
    // dimensions are equal to the number of degrees of freedom), and also how
    // many rows out of this global size are to be stored locally
    // (<code>n_local_dofs</code>). In addition, PETSc needs to know how to
    // partition the columns in the chunk of the matrix that is stored
    // locally; for square matrices, the columns should be partitioned in the
    // same way as the rows (indicated by the second <code>n_local_dofs</code>
    // in the call) but in the case of rectangular matrices one has to
    // partition the columns in the same way as vectors are partitioned with
    // which the matrix is multiplied, while rows have to partitioned in the
    // same way as destination vectors of matrix-vector multiplications:
    system_matrix.reinit (mpi_communicator,
                          dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          n_local_dofs,
                          n_local_dofs,
                          dof_handler.max_couplings_between_dofs());

    solution.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
    system_rhs.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

    // Finally, we need to initialize the objects denoting hanging node
    // constraints for the present grid. Note that since PETSc handles the
    // sparsity pattern internally to the matrix, there is no need to set up
    // an independent sparsity pattern here, and to condense it for
    // constraints, as we have done in all other example programs.
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             hanging_node_constraints);
    hanging_node_constraints.close ();
  }


  // The third step is to actually assemble the matrix and right hand side of
  // the problem. There are some things worth mentioning before we go into
  // detail. First, we will be assembling the system in %parallel, i.e. each
  // process will be responsible for assembling on cells that belong to this
  // particular processor. Note that the degrees of freedom are split in a way
  // such that all DoFs in the interior of cells and between cells belonging
  // to the same subdomain belong to the process that <code>owns</code> the
  // cell. However, even then we sometimes need to assemble on a cell with a
  // neighbor that belongs to a different process, and in these cases when we
  // write the local contributions into the global matrix or right hand side
  // vector, we actually have to transfer these entries to the other
  // process. Fortunately, we don't have to do this by hand, PETSc does all
  // this for us by caching these elements locally, and sending them to the
  // other processes as necessary when we call the <code>compress()</code>
  // functions on the matrix and vector at the end of this function.
  //
  // The second point is that once we have handed over matrix and vector
  // contributions to PETSc, it is a) hard, and b) very inefficient to get
  // them back for modifications. This is not only the fault of PETSc, it is
  // also a consequence of the distributed nature of this program: if an entry
  // resides on another processor, then it is necessarily expensive to get
  // it. The consequence of this is that where we previously first assembled
  // the matrix and right hand side as if there were no hanging node
  // constraints and boundary values, and then eliminated these in a second
  // step, we should now try to do that while still assembling the local
  // systems, and before handing these entries over to PETSc. At least as far
  // as eliminating hanging nodes is concerned, this is actually possible,
  // though removing boundary nodes isn't that simple. deal.II provides
  // functions to do this first part: instead of copying elements by hand into
  // the global matrix, we use the <code>distribute_local_to_global</code>
  // functions below to take care of hanging nodes at the same time. The
  // second step, elimination of boundary nodes, is then done in exactly the
  // same way as in all previous example programs.
  //
  // So, here is the actual implementation:
  template <int dim>
  void ElasticProblem<dim>::assemble_system ()
  {
    // The infrastructure to assemble linear systems is the same as in all the
    // other programs, and in particular unchanged from step-8. Note that we
    // still use the deal.II full matrix and vector types for the local
    // systems.
    QGauss<dim>  quadrature_formula(2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     lambda_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    ConstantFunction<dim> lambda(1.), mu(1.);

    RightHandSide<dim>      right_hand_side;
    std::vector<Vector<double> > rhs_values (n_q_points,
                                             Vector<double>(dim));


    // The next thing is the loop over all elements. Note that we do not have
    // to do all the work: our job here is only to assemble the system on
    // cells that actually belong to this MPI process, all other cells will be
    // taken care of by other processes. This is what the if-clause
    // immediately after the for-loop takes care of: it queries the subdomain
    // identifier of each cell, which is a number associated with each cell
    // that tells which process handles it. In more generality, the subdomain
    // id is used to split a domain into several parts (we do this above, at
    // the beginning of <code>setup_system</code>), and which allows to
    // identify which subdomain a cell is living on. In this application, we
    // have each process handle exactly one subdomain, so we identify the
    // terms <code>subdomain</code> and <code>MPI process</code> with each
    // other.
    //
    // Apart from this, assembling the local system is relatively uneventful
    // if you have understood how this is done in step-8, and only becomes
    // interesting again once we start distributing it into the global matrix
    // and right hand sides.
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->subdomain_id() == this_mpi_process)
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit (cell);

          lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
          mu.value_list     (fe_values.get_quadrature_points(), mu_values);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int
              component_i = fe.system_to_component_index(i).first;

              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  const unsigned int
                  component_j = fe.system_to_component_index(j).first;

                  for (unsigned int q_point=0; q_point<n_q_points;
                       ++q_point)
                    {
//TODO investigate really small values here
                      cell_matrix(i,j)
                      +=
                        (
                          (fe_values.shape_grad(i,q_point)[component_i] *
                           fe_values.shape_grad(j,q_point)[component_j] *
                           lambda_values[q_point])
                          +
                          (fe_values.shape_grad(i,q_point)[component_j] *
                           fe_values.shape_grad(j,q_point)[component_i] *
                           mu_values[q_point])
                          +
                          ((component_i == component_j) ?
                           (fe_values.shape_grad(i,q_point) *
                            fe_values.shape_grad(j,q_point) *
                            mu_values[q_point])  :
                           0)
                        )
                        *
                        fe_values.JxW(q_point);
                    }
                }
            }

          right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
                                             rhs_values);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int
              component_i = fe.system_to_component_index(i).first;

              for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                cell_rhs(i) += fe_values.shape_value(i,q_point) *
                               rhs_values[q_point](component_i) *
                               fe_values.JxW(q_point);
            }

          // Now we have the local system, and need to transfer it into the
          // global objects. However, as described in the introduction to this
          // function, we want to avoid any operations to matrix and vector
          // entries after handing them off to PETSc (i.e. after distributing
          // to the global objects). Therefore, we will take care of hanging
          // node constraints already here. This is not quite trivial since
          // the rows and columns of constrained nodes have to be distributed
          // to the rows and columns of those nodes to which they are
          // constrained. This can't be done on a purely local basis (because
          // the degrees of freedom to which hanging nodes are constrained may
          // not be associated with the cell we are presently treating, and
          // are therefore not represented in the local matrix and vector),
          // but it can be done while distributing the local system to the
          // global one. This is what the following call does, i.e. we
          // distribute to the global objects and at the same time make sure
          // that hanging node constraints are taken care of:
          cell->get_dof_indices (local_dof_indices);
          hanging_node_constraints
          .distribute_local_to_global(cell_matrix, cell_rhs,
                                      local_dof_indices,
                                      system_matrix, system_rhs);
        }

    // Now compress the vector and the system matrix:
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    // The global matrix and right hand side vectors have now been
    // formed. Note that since we took care of this already above, we do not
    // have to condense away hanging node constraints any more.
    //
    // However, we still have to apply boundary values, in the same way as we
    // always do:
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(dim),
                                              boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix, solution,
                                        system_rhs, false);
    // The last argument to the call just performed allows for some
    // optimizations. It controls whether we should also delete the column
    // corresponding to a boundary node, or keep it (and passing
    // <code>true</code> means: yes, do eliminate the column). If we
    // do, then the resulting matrix will be symmetric again if it was before;
    // if we don't, then it won't. The solution of the resulting system should
    // be the same, though. The only reason why we may want to make the system
    // symmetric again is that we would like to use the CG method, which only
    // works with symmetric matrices.  Experience tells that CG also works
    // (and works almost as well) if we don't remove the columns associated
    // with boundary nodes, which can be easily explained by the special
    // structure of the non-symmetry. Since eliminating columns from dense
    // matrices is not expensive, though, we let the function do it; not doing
    // so is more important if the linear system is either non-symmetric
    // anyway, or we are using the non-local version of this function (as in
    // all the other example programs before) and want to save a few cycles
    // during this operation.
  }



  // The fourth step is to solve the linear system, with its distributed
  // matrix and vector objects. Fortunately, PETSc offers a variety of
  // sequential and %parallel solvers, for which we have written wrappers that
  // have almost the same interface as is used for the deal.II solvers used in
  // all previous example programs.
  template <int dim>
  unsigned int ElasticProblem<dim>::solve ()
  {
    // First, we have to set up a convergence monitor, and assign it the
    // accuracy to which we would like to solve the linear system. Next, an
    // actual solver object using PETSc's CG solver which also works with
    // %parallel (distributed) vectors and matrices. And finally a
    // preconditioner; we choose to use a block Jacobi preconditioner which
    // works by computing an incomplete LU decomposition on each block
    // (i.e. the chunk of matrix that is stored on each MPI process). That
    // means that if you run the program with only one process, then you will
    // use an ILU(0) as a preconditioner, while if it is run on many
    // processes, then we will have a number of blocks on the diagonal and the
    // preconditioner is the ILU(0) of each of these blocks.
    SolverControl           solver_control (solution.size(),
                                            1e-8*system_rhs.l2_norm());
    PETScWrappers::SolverCG cg (solver_control,
                                mpi_communicator);

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

    // Then solve the system:
    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);

    // The next step is to distribute hanging node constraints. This is a
    // little tricky, since to fill in the value of a constrained node you
    // need access to the values of the nodes to which it is constrained (for
    // example, for a Q1 element in 2d, we need access to the two nodes on the
    // big side of a hanging node face, to compute the value of the
    // constrained node in the middle). Since PETSc (and, for that matter, the
    // MPI model on which it is built) does not allow to query the value of
    // another node in a simple way if we should need it, what we do here is
    // to get a copy of the distributed vector where we keep all elements
    // locally. This is simple, since the deal.II wrappers have a conversion
    // constructor for the non-MPI vector class:
    PETScWrappers::Vector localized_solution (solution);

    // Then we distribute hanging node constraints on this local copy, i.e. we
    // compute the values of all constrained nodes:
    hanging_node_constraints.distribute (localized_solution);

    // Then transfer everything back into the global vector. The following
    // operation copies those elements of the localized solution that we store
    // locally in the distributed solution, and does not touch the other
    // ones. Since we do the same operation on all processors, we end up with
    // a distributed vector that has all the constrained nodes fixed.
    solution = localized_solution;

    // Finally return the number of iterations it took to converge, to allow
    // for some output:
    return solver_control.last_step();
  }



  // Step five is to output the results we computed in this iteration. This is
  // actually the same as done in step-8 before, with two small
  // differences. First, all processes call this function, but not all of them
  // need to do the work associated with generating output. In fact, they
  // shouldn't, since we would try to write to the same file multiple times at
  // once. So we let only the first job do this, and all the other ones idle
  // around during this time (or start their work for the next iteration, or
  // simply yield their CPUs to other jobs that happen to run at the same
  // time). The second thing is that we not only output the solution vector,
  // but also a vector that indicates which subdomain each cell belongs
  // to. This will make for some nice pictures of partitioned domains.
  //
  // In practice, the present implementation of the output function is a major
  // bottleneck of this program, since generating graphical output is
  // expensive and doing so only on one process does, of course, not scale if
  // we significantly increase the number of processes. In effect, this
  // function will consume most of the run-time if you go to very large
  // numbers of unknowns and processes, and real applications should limit the
  // number of times they generate output through this function.
  //
  // The solution to this is to have each process generate output data only
  // for it's own local cells, and write them to separate files, one file per
  // process. This would distribute the work of generating the output to all
  // processes equally. In a second step, separate from running this program,
  // we would then take all the output files for a given cycle and merge these
  // parts into one single output file. This has to be done sequentially, but
  // can be done on a different machine, and should be relatively
  // cheap. However, the necessary functionality for this is not yet
  // implemented in the library, and since we are too close to the next
  // release, we do not want to do such major destabilizing changes any
  // more. This has been fixed in the meantime, though, and a better way to do
  // things is explained in the step-18 example program.
  template <int dim>
  void ElasticProblem<dim>::output_results (const unsigned int cycle) const
  {
    // One point to realize is that when we want to generate output on process
    // zero only, we need to have access to all elements of the solution
    // vector. So we need to get a local copy of the distributed vector, which
    // is in fact simple:
    const PETScWrappers::Vector localized_solution (solution);
    // The thing to notice, however, is that we do this localization operation
    // on all processes, not only the one that actually needs the data. This
    // can't be avoided, however, with the communication model of MPI: MPI
    // does not have a way to query data on another process, both sides have
    // to initiate a communication at the same time. So even though most of
    // the processes do not need the localized solution, we have to place the
    // call here so that all processes execute it.
    //
    // (In reality, part of this work can in fact be avoided. What we do is
    // send the local parts of all processes to all other processes. What we
    // would really need to do is to initiate an operation on all processes
    // where each process simply sends its local chunk of data to process
    // zero, since this is the only one that actually needs it, i.e. we need
    // something like a gather operation. PETSc can do this, but for
    // simplicity's sake we don't attempt to make use of this here. We don't,
    // since what we do is not very expensive in the grand scheme of things:
    // it is one vector communication among all processes , which has to be
    // compared to the number of communications we have to do when solving the
    // linear system, setting up the block-ILU for the preconditioner, and
    // other operations.)

    // This being done, process zero goes ahead with setting up the output
    // file as in step-8, and attaching the (localized) solution vector to the
    // output object:. (The code to generate the output file name is stolen
    // and slightly modified from step-5, since we expect that we can do a
    // number of cycles greater than 10, which is the maximum of what the code
    // in step-8 could handle.)
    if (this_mpi_process == 0)
      {
        std::ostringstream filename;
        filename << "solution-" << cycle << ".gmv";

        std::ofstream output (filename.str().c_str());

        DataOut<dim> data_out;
        data_out.attach_dof_handler (dof_handler);

        std::vector<std::string> solution_names;
        switch (dim)
          {
          case 1:
            solution_names.push_back ("displacement");
            break;
          case 2:
            solution_names.push_back ("x_displacement");
            solution_names.push_back ("y_displacement");
            break;
          case 3:
            solution_names.push_back ("x_displacement");
            solution_names.push_back ("y_displacement");
            solution_names.push_back ("z_displacement");
            break;
          default:
            Assert (false, ExcInternalError());
          }

        data_out.add_data_vector (localized_solution, solution_names);

        // The only thing we do here additionally is that we also output one
        // value per cell indicating which subdomain (i.e. MPI process) it
        // belongs to. This requires some conversion work, since the data the
        // library provides us with is not the one the output class expects,
        // but this is not difficult. First, set up a vector of integers, one
        // per cell, that is then filled by the number of subdomain each cell
        // is in:
        std::vector<unsigned int> partition_int (triangulation.n_active_cells());
        GridTools::get_subdomain_association (triangulation, partition_int);

        // Then convert this integer vector into a floating point vector just
        // as the output functions want to see:
        const Vector<double> partitioning(partition_int.begin(),
                                          partition_int.end());

        // And finally add this vector as well:
        data_out.add_data_vector (partitioning, "partitioning");

        // This all being done, generate the intermediate format and write it
        // out in GMV output format:
        data_out.build_patches ();
        data_out.write_gmv (output);
      }
  }



  // The sixth step is to take the solution just computed, and evaluate some
  // kind of refinement indicator to refine the mesh. The problem is basically
  // the same as with distributing hanging node constraints: in order to
  // compute the error indicator, we need access to all elements of the
  // solution vector. We then compute the indicators for the cells that belong
  // to the present process, but then we need to distribute the refinement
  // indicators into a distributed vector so that all processes have the
  // values of the refinement indicator for all cells. But then, in order for
  // each process to refine its copy of the mesh, they need to have access to
  // all refinement indicators locally, so they have to copy the global vector
  // back into a local one. That's a little convoluted, but thinking about it
  // quite straightforward nevertheless. So here's how we do it:
  template <int dim>
  void ElasticProblem<dim>::refine_grid ()
  {
    // So, first part: get a local copy of the distributed solution
    // vector. This is necessary since the error estimator needs to get at the
    // value of neighboring cells even if they do not belong to the subdomain
    // associated with the present MPI process:
    const PETScWrappers::Vector localized_solution (solution);

    // Second part: set up a vector of error indicators for all cells and let
    // the Kelly class compute refinement indicators for all cells belonging
    // to the present subdomain/process. Note that the last argument of the
    // call indicates which subdomain we are interested in. The three
    // arguments before it are various other default arguments that one
    // usually doesn't need (and doesn't state values for, but rather uses the
    // defaults), but which we have to state here explicitly since we want to
    // modify the value of a following argument (i.e. the one indicating the
    // subdomain):
    Vector<float> local_error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(2),
                                        typename FunctionMap<dim>::type(),
                                        localized_solution,
                                        local_error_per_cell,
                                        ComponentMask(),
                                        0,
                                        multithread_info.n_threads(),
                                        this_mpi_process);

    // Now all processes have computed error indicators for their own cells
    // and stored them in the respective elements of the
    // <code>local_error_per_cell</code> vector. The elements of this vector
    // for cells not on the present process are zero. However, since all
    // processes have a copy of a copy of the entire triangulation and need to
    // keep these copies in sync, they need the values of refinement
    // indicators for all cells of the triangulation. Thus, we need to
    // distribute our results. We do this by creating a distributed vector
    // where each process has its share, and sets the elements it has
    // computed. We will then later generate a local sequential copy of this
    // distributed vector to allow each process to access all elements of this
    // vector.
    //
    // So in the first step, we need to set up a %parallel vector. For
    // simplicity, every process will own a chunk with as many elements as
    // this process owns cells, so that the first chunk of elements is stored
    // with process zero, the next chunk with process one, and so on. It is
    // important to remark, however, that these elements are not necessarily
    // the ones we will write to. This is so, since the order in which cells
    // are arranged, i.e. the order in which the elements of the vector
    // correspond to cells, is not ordered according to the subdomain these
    // cells belong to. In other words, if on this process we compute
    // indicators for cells of a certain subdomain, we may write the results
    // to more or less random elements if the distributed vector, that do not
    // necessarily lie within the chunk of vector we own on the present
    // process. They will subsequently have to be copied into another
    // process's memory space then, an operation that PETSc does for us when
    // we call the <code>compress</code> function. This inefficiency could be
    // avoided with some more code, but we refrain from it since it is not a
    // major factor in the program's total runtime.
    //
    // So here's how we do it: count how many cells belong to this process,
    // set up a distributed vector with that many elements to be stored
    // locally, and copy over the elements we computed locally, then compress
    // the result. In fact, we really only copy the elements that are nonzero,
    // so we may miss a few that we computed to zero, but this won't hurt
    // since the original values of the vector is zero anyway.
    const unsigned int n_local_cells
      = GridTools::count_cells_with_subdomain_association (triangulation,
                                                           this_mpi_process);
    PETScWrappers::MPI::Vector
    distributed_all_errors (mpi_communicator,
                            triangulation.n_active_cells(),
                            n_local_cells);

    for (unsigned int i=0; i<local_error_per_cell.size(); ++i)
      if (local_error_per_cell(i) != 0)
        distributed_all_errors(i) = local_error_per_cell(i);
    distributed_all_errors.compress (VectorOperation::insert);


    // So now we have this distributed vector out there that contains the
    // refinement indicators for all cells. To use it, we need to obtain a
    // local copy...
    const Vector<float> localized_all_errors (distributed_all_errors);

    // ...which we can the subsequently use to finally refine the grid:
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     localized_all_errors,
                                                     0.3, 0.03);
    triangulation.execute_coarsening_and_refinement ();
  }



  // Lastly, here is the driver function. It is almost unchanged from step-8,
  // with the exception that we replace <code>std::cout</code> by the
  // <code>pcout</code> stream. Apart from this, the only other cosmetic
  // change is that we output how many degrees of freedom there are per
  // process, and how many iterations it took for the linear solver to
  // converge:
  template <int dim>
  void ElasticProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<10; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation, -1, 1);
            triangulation.refine_global (3);
          }
        else
          refine_grid ();

        pcout << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;

        setup_system ();

        pcout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (by partition:";
        for (unsigned int p=0; p<n_mpi_processes; ++p)
          pcout << (p==0 ? ' ' : '+')
                << (DoFTools::
                    count_dofs_with_subdomain_association (dof_handler,
                                                           p));
        pcout << ")" << std::endl;

        assemble_system ();
        const unsigned int n_iterations = solve ();

        pcout << "   Solver converged in " << n_iterations
              << " iterations." << std::endl;

        output_results (cycle);
      }
  }
}


// So that's it, almost. <code>main()</code> works the same way as most of the
// main functions in the other example programs, i.e. it delegates work to the
// <code>run</code> function of a master object, and only wraps everything
// into some code to catch exceptions:
int main (int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step17;

      // Here is the only real difference: PETSc requires that we initialize
      // it at the beginning of the program, and un-initialize it at the
      // end. The class MPI_InitFinalize takes care of that. The original code
      // sits in between, enclosed in braces to make sure that the
      // <code>elastic_problem</code> variable goes out of scope (and is
      // destroyed) before PETSc is closed with
      // <code>PetscFinalize</code>. (If we wouldn't use braces, the
      // destructor of <code>elastic_problem</code> would run after
      // <code>PetscFinalize</code>; since the destructor involves calls to
      // PETSc functions, we would get strange error messages from PETSc.)
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        deallog.depth_console (0);

        ElasticProblem<2> elastic_problem;
        elastic_problem.run ();
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
