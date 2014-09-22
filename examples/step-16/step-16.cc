/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2003 - 2014 by the deal.II authors
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
 * Authors: Guido Kanschat, University of Heidelberg, 2003
 *          Baerbel Janssen, University of Heidelberg, 2010
 *          Wolfgang Bangerth, Texas A&M University, 2010
 */


// @sect3{Include files}

// Again, the first few include files are already known, so we won't comment
// on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// These, now, are the include necessary for the multilevel methods. The first
// one declares how to handle Dirichlet boundary conditions on each of the
// levels of the multigrid method. For the actual description of the degrees
// of freedom, we do not need any new include file because DoFHandler already
// has all necessary methods implemented. We will only need to distribute the
// DoFs for the levels further down.
//
// The rest of the include files deals with the mechanics of multigrid as a
// linear operator (solver or preconditioner).
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

// Finally we include the MeshWorker framework. This framework through
// its function loop() and integration_loop(), automates loops over
// cells and assembling of data into vectors, matrices, etc. It obeys
// constraints automatically. Since we have to build
// several matrices and have to be aware of several sets of
// constraints, this will save us a lot of headache.
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/output.h>
#include <deal.II/meshworker/loop.h>

// In order to save effort, we use the pre-implemented Laplacian found in
#include <deal.II/integrators/laplace.h>
#include <deal.II/integrators/l2.h>

// This is C++:
#include <fstream>
#include <sstream>

using namespace dealii;

namespace Step16
{
  // @sect3{The integrator on each cell}

  // The MeshWorker::integration_loop() expects a class that provides
  // functions for integration on cells and boundary and interior
  // faces. This is done by the following class. In the constructor,
  // we tell the loop that cell integrals should be computed (the
  // 'true'), but integrals should not be computed on boundary and
  // interior faces (the two 'false'). Accordingly, we only need a
  // cell function, but none for the faces.

  template <int dim>
  class LaplaceIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    LaplaceIntegrator();
    virtual void cell(MeshWorker::DoFInfo<dim> &dinfo, MeshWorker::IntegrationInfo<dim> &info) const;
  };


  template <int dim>
  LaplaceIntegrator<dim>::LaplaceIntegrator()
    :
    MeshWorker::LocalIntegrator<dim>(true, false, false)
  {}


  // Next the actual integrator on each cell. We solve a Poisson problem with a
  // coefficient one in the right half plane and one tenth in the left
  // half plane.

  // The MeshWorker::LocalResults base class of MeshWorker::DoFInfo
  // contains objects that can be filled in this local integrator. How
  // many objects is determined inside the MeshWorker framework by the
  // assembler class. Here, we test for instance that one matrix is
  // required (MeshWorker::LocalResults::n_matrices()). The matrices are accessed
  // through MeshWorker::LocalResults::matrix(), which takes the number of the
  // matrix as its first argument. The second argument is only used
  // for integrals over faces, where there are two matrices for each
  // test function set. In such a case, a second matrix with indicator
  // 'true' would exist with the same index.

  // MeshWorker::IntegrationInfo provides one or several FEValues
  // objects, which below are used by
  // LocalIntegrators::Laplace::cell_matrix() or
  // LocalIntegrators::L2::L2(). Since we are assembling only a single
  // PDE, there is also only one of these objects with index zero.

  // In addition, we note that this integrator serves to compute the
  // matrices for the multilevel preconditioner as well as the matrix
  // and the right hand side for the global system. Since the
  // assembler for a system requires an additional vector, this is
  // indicated by MeshWorker::LocalResults::n_vectors() returning a nonzero
  // value. Accordingly we fill a right hand side vector at the end of
  // this function. Since LocalResults can deal with several
  // BlockVector objects, but we are again in the simplest case here,
  // we enter the information into block zero of vector zero.
  template <int dim>
  void LaplaceIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, MeshWorker::IntegrationInfo<dim> &info) const
  {
    AssertDimension (dinfo.n_matrices(), 1);
    const double coefficient = (dinfo.cell->center()(0) > 0.)
                               ? .1 : 1.;

    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0), coefficient);

    if (dinfo.n_vectors() > 0)
      {
        std::vector<double> rhs(info.fe_values(0).n_quadrature_points, 1.);
        LocalIntegrators::L2::L2(dinfo.vector(0).block(0), info.fe_values(0), rhs);
      }
  }


  // @sect3{The <code>LaplaceProblem</code> class template}

  // This main class is basically the same class as in step-6. As far as
  // member functions is concerned, the only addition is the
  // <code>assemble_multigrid</code> function that assembles the matrices that
  // correspond to the discrete operators on intermediate levels:
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem (const unsigned int degree);
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    // We need an additional object for the hanging nodes constraints. They
    // are handed to the transfer object in the multigrid. Since we call a
    // compress inside the multigrid these constraints are not allowed to be
    // inhomogeneous so we store them in different ConstraintMatrix objects.
    ConstraintMatrix     hanging_node_constraints;
    ConstraintMatrix     constraints;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    const unsigned int degree;

    // The following members are the essential data structures for the
    // multigrid method. The first two represent the sparsity patterns
    // and the matrices on individual levels of the multilevel
    // hierarchy, very much like the objects for the global mesh above.

    // Then we have two new matrices only needed for multigrid
    // methods with local smoothing on adaptive meshes. They convey
    // data between the interior part of the refined region and the
    // refinement edge, as outline in detail in @ref mg_paper.

    // The last object stores information about the boundary indices
    // on each level and information about indices lying on a
    // refinement edge between two different refinement levels. It
    // thus serves a similar purpose as ConstraintMatrix, but on each
    // level.
    MGLevelObject<SparsityPattern>       mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double> > mg_matrices;
    MGLevelObject<SparseMatrix<double> > mg_interface_in;
    MGLevelObject<SparseMatrix<double> > mg_interface_out;
    MGConstrainedDoFs                    mg_constrained_dofs;
  };


  // @sect3{The <code>LaplaceProblem</code> class implementation}

  // Just one short remark about the constructor of the Triangulation:
  // by convention, all adaptively refined triangulations in deal.II never
  // change by more than one level across a face between cells. For our
  // multigrid algorithms, however, we need a slightly stricter guarantee,
  // namely that the mesh also does not change by more than refinement level
  // across vertices that might connect two cells. In other words, we must
  // prevent the following situation:
  //
  // @image html limit_level_difference_at_vertices.png ""
  //
  // This is achieved by passing the
  // Triangulation::limit_level_difference_at_vertices flag to the constructor
  // of the triangulation class.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem (const unsigned int degree)
    :
    triangulation (Triangulation<dim>::
                   limit_level_difference_at_vertices),
    fe (degree),
    dof_handler (triangulation),
    degree(degree)
  {}



  // @sect4{LaplaceProblem::setup_system}

  // In addition to just distributing the degrees of freedom in
  // the DoFHandler, we do the same on each level. Then, we follow the
  // same procedure as before to set up the system on the leaf mesh.
  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
    dof_handler.distribute_mg_dofs (fe);

    deallog << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (by level: ";
    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
      deallog << dof_handler.n_dofs(level)
              << (level == triangulation.n_levels()-1
                  ? ")" : ", ");
    deallog << std::endl;

    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    constraints.clear ();
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);

    typename FunctionMap<dim>::type      dirichlet_boundary_functions;
    ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
    dirichlet_boundary_functions[0] = &homogeneous_dirichlet_bc;
    VectorTools::interpolate_boundary_values (static_cast<const DoFHandler<dim>&>(dof_handler),
                                              dirichlet_boundary_functions,
                                              constraints);
    constraints.close ();
    hanging_node_constraints.close ();
    constraints.condense (sparsity_pattern);
    sparsity_pattern.compress();
    system_matrix.reinit (sparsity_pattern);

    // The multigrid constraints have to be initialized. They need to know
    // about the boundary values as well, so we pass the
    // <code>dirichlet_boundary</code> here as well.
    mg_constrained_dofs.clear();
    mg_constrained_dofs.initialize(dof_handler, dirichlet_boundary_functions);


    // Now for the things that concern the multigrid data structures. First,
    // we resize the multilevel objects to hold matrices and sparsity
    // patterns for every level. The coarse level is zero (this is mandatory
    // right now but may change in a future revision). Note that these
    // functions take a complete, inclusive range here (not a starting index
    // and size), so the finest level is <code>n_levels-1</code>.  We first
    // have to resize the container holding the SparseMatrix classes, since
    // they have to release their SparsityPattern before the can be destroyed
    // upon resizing.
    const unsigned int n_levels = triangulation.n_levels();

    mg_interface_in.resize(0, n_levels-1);
    mg_interface_in.clear ();
    mg_interface_out.resize(0, n_levels-1);
    mg_interface_out.clear ();
    mg_matrices.resize(0, n_levels-1);
    mg_matrices.clear ();
    mg_sparsity_patterns.resize(0, n_levels-1);

    // Now, we have to provide a matrix on each level. To this end, we first
    // use the MGTools::make_sparsity_pattern function to first generate a
    // preliminary compressed sparsity pattern on each level (see the @ref
    // Sparsity module for more information on this topic) and then copy it
    // over to the one we really want. The next step is to initialize both
    // kinds of level matrices with these sparsity patterns.
    //
    // It may be worth pointing out that the interface matrices only have
    // entries for degrees of freedom that sit at or next to the interface
    // between coarser and finer levels of the mesh. They are therefore even
    // sparser than the matrices on the individual levels of our multigrid
    // hierarchy. If we were more concerned about memory usage (and possibly
    // the speed with which we can multiply with these matrices), we should
    // use separate and different sparsity patterns for these two kinds of
    // matrices.
    for (unsigned int level=0; level<n_levels; ++level)
      {
        CompressedSparsityPattern csp;
        csp.reinit(dof_handler.n_dofs(level),
                   dof_handler.n_dofs(level));
        MGTools::make_sparsity_pattern(dof_handler, csp, level);

        mg_sparsity_patterns[level].copy_from (csp);

        mg_matrices[level].reinit(mg_sparsity_patterns[level]);
        mg_interface_in[level].reinit(mg_sparsity_patterns[level]);
        mg_interface_out[level].reinit(mg_sparsity_patterns[level]);
      }
  }


  // @sect4{LaplaceProblem::assemble_system}

  // The following function assembles the linear system on the finest level of
  // the mesh. Since we want to reuse the code here for the level
  // assembling below, we use the local integrator class
  // LaplaceIntegrator and leave the loops to the MeshWorker
  // framework. Thus, this function first sets up the objects
  // necessary for this framework, namely
  // <ol>
  // <li>an MeshWorker::IntegrationInfoBox, which will provide all the required
  // data in quadrature points on the cell. This object can be seen as
  // an extension of FEValues, providing a lot more useful
  // information,</li>
  // <li>a MeshWorker::DoFInfo object, which on the one hand side extends the
  // functionality of cell iterators, but also provides space for
  // return values in its base class LocalResults,</li>
  // <li>an assembler, in this case for the whole system. The term
  // 'simple' here refers to the fact that the global system does not
  // have a block structure,</li>
  // <li>an the local integrator, which implements the actual forms.
  // </ol>
  //
  // After the loop has combined all of these into a matrix and a
  // right hand side, there is one thing left to do: the assemblers
  // leave matrix rows and columns of constrained degrees of freedom
  // untouched. Therefore, we put a one on the diagonal to make the
  // whole system well posed. The value one, or any fixed value has
  // the advantage, that its effect on the spectrum of the matrix is
  // easily understood. Since the corresponding eigenvectors form an
  // invariant subspace, the value chosen does not affect the
  // convergence of Krylov space solvers.
  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    MappingQ1<dim> mapping;
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients | update_hessians;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> > assembler;
    assembler.initialize(constraints);
    assembler.initialize(system_matrix, system_rhs);

    LaplaceIntegrator<dim> matrix_integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box, matrix_integrator, assembler);

    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        system_matrix.set(i, i, 1.);
  }


  // @sect4{LaplaceProblem::assemble_multigrid}

  // The next function is the one that builds the linear operators (matrices)
  // that define the multigrid method on each level of the mesh. The
  // integration core is the same as above, but the loop below will go over
  // all existing cells instead of just the active ones, and the results must
  // be entered into the correct level matrices. Fortunately,
  // MeshWorker hides most of that from us, and thus the difference
  // between this function and the previous lies only in the setup of
  // the assembler and the different iterators in the loop.
  // Also, fixing up the matrices in the end is a little more comlicated.
  template <int dim>
  void LaplaceProblem<dim>::assemble_multigrid ()
  {
    MappingQ1<dim> mapping;
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients | update_hessians;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(mg_constrained_dofs);
    assembler.initialize(mg_matrices);
    assembler.initialize_interfaces(mg_interface_in, mg_interface_out);

    LaplaceIntegrator<dim> matrix_integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_mg(), dof_handler.end_mg(),
      dof_info, info_box, matrix_integrator, assembler);

    const unsigned int nlevels = triangulation.n_levels();
    for (unsigned int level=0; level<nlevels; ++level)
      {
        for (unsigned int i=0; i<dof_handler.n_dofs(level); ++i)
          if (mg_constrained_dofs.is_boundary_index(level,i) ||
              mg_constrained_dofs.at_refinement_edge(level,i))
            mg_matrices[level].set(i, i, 1.);
      }
  }



  // @sect4{LaplaceProblem::solve}

  // This is the other function that is significantly different in support of
  // the multigrid solver (or, in fact, the preconditioner for which we use
  // the multigrid method).
  //
  // Let us start out by setting up two of the components of multilevel
  // methods: transfer operators between levels, and a solver on the coarsest
  // level. In finite element methods, the transfer operators are derived from
  // the finite element function spaces involved and can often be computed in
  // a generic way independent of the problem under consideration. In that
  // case, we can use the MGTransferPrebuilt class that, given the constraints
  // of the final linear system and the MGConstrainedDoFs object that knows
  // about the boundary conditions on the each level and the degrees of
  // freedom on interfaces between different refinement level can build the
  // matrices for those transfer operations from a DoFHandler object with
  // level degrees of freedom.
  //
  // The second part of the following lines deals with the coarse grid
  // solver. Since our coarse grid is very coarse indeed, we decide for a
  // direct solver (a Householder decomposition of the coarsest level matrix),
  // even if its implementation is not particularly sophisticated. If our
  // coarse mesh had many more cells than the five we have here, something
  // better suited would obviously be necessary here.
  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    MGTransferPrebuilt<Vector<double> > mg_transfer(hanging_node_constraints, mg_constrained_dofs);
    mg_transfer.build_matrices(dof_handler);

    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from (mg_matrices[0]);
    MGCoarseGridHouseholder<> coarse_grid_solver;
    coarse_grid_solver.initialize (coarse_matrix);

    // The next component of a multilevel solver or preconditioner is that we
    // need a smoother on each level. A common choice for this is to use the
    // application of a relaxation method (such as the SOR, Jacobi or
    // Richardson method) or a small number of iterations of a solver method
    // (such as CG or GMRES). The mg::SmootherRelaxation and
    // MGSmootherPrecondition classes provide support for these two kinds of
    // smoothers. Here, we opt for the application of a single SOR
    // iteration. To this end, we define an appropriate <code>typedef</code>
    // and then setup a smoother object.
    //
    // Since this smoother needs temporary vectors to store intermediate
    // results, we need to provide a VectorMemory object. Since these vectors
    // will be reused over and over, the GrowingVectorMemory is more time
    // efficient than the PrimitiveVectorMemory class in the current case.
    //
    // The last step is to initialize the smoother object with our level
    // matrices and to set some smoothing parameters.  The
    // <code>initialize()</code> function can optionally take additional
    // arguments that will be passed to the smoother object on each level. In
    // the current case for the SOR smoother, this could, for example, include
    // a relaxation parameter. However, we here leave these at their default
    // values. The call to <code>set_steps()</code> indicates that we will use
    // two pre- and two post-smoothing steps on each level; to use a variable
    // number of smoother steps on different levels, more options can be set
    // in the constructor call to the <code>mg_smoother</code> object.
    //
    // The last step results from the fact that we use the SOR method as a
    // smoother - which is not symmetric - but we use the conjugate gradient
    // iteration (which requires a symmetric preconditioner) below, we need to
    // let the multilevel preconditioner make sure that we get a symmetric
    // operator even for nonsymmetric smoothers:
    typedef PreconditionSOR<SparseMatrix<double> > Smoother;
    mg::SmootherRelaxation<Smoother, Vector<double> > mg_smoother;
    mg_smoother.initialize(mg_matrices);
    mg_smoother.set_steps(2);
    mg_smoother.set_symmetric(true);

    // The next preparatory step is that we must wrap our level and interface
    // matrices in an object having the required multiplication functions. We
    // will create two objects for the interface objects going from coarse to
    // fine and the other way around; the multigrid algorithm will later use
    // the transpose operator for the latter operation, allowing us to
    // initialize both up and down versions of the operator with the matrices
    // we already built:
    mg::Matrix<Vector<double> > mg_matrix(mg_matrices);
    mg::Matrix<Vector<double> > mg_interface_up(mg_interface_in);
    mg::Matrix<Vector<double> > mg_interface_down(mg_interface_out);

    // Now, we are ready to set up the V-cycle operator and the multilevel
    // preconditioner.
    Multigrid<Vector<double> > mg(dof_handler,
                                  mg_matrix,
                                  coarse_grid_solver,
                                  mg_transfer,
                                  mg_smoother,
                                  mg_smoother);
    mg.set_edge_matrices(mg_interface_down, mg_interface_up);

    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double> > >
    preconditioner(dof_handler, mg, mg_transfer);

    // With all this together, we can finally get about solving the linear
    // system in the usual way:
    SolverControl solver_control (1000, 1e-12);
    SolverCG<>    solver (solver_control);

    solution = 0;

    solver.solve (system_matrix, solution, system_rhs,
                  preconditioner);
    constraints.distribute (solution);
  }



  // @sect4{Postprocessing}

  // The following two functions postprocess a solution once it is
  // computed. In particular, the first one refines the mesh at the beginning
  // of each cycle while the second one outputs results at the end of each
  // such cycle. The functions are almost unchanged from those in step-6, with
  // the exception of one minor difference: we generate output in VTK
  // format, to use the more modern visualization programs available today
  // compared to those that were available when step-6 was written.
  template <int dim>
  void LaplaceProblem<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);
    triangulation.execute_coarsening_and_refinement ();
  }



  template <int dim>
  void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();

    std::ostringstream filename;
    filename << "solution-"
             << cycle
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
  }


  // @sect4{LaplaceProblem::run}

  // Like several of the functions above, this is almost exactly a copy of of
  // the corresponding function in step-6. The only difference is the call to
  // <code>assemble_multigrid</code> that takes care of forming the matrices
  // on every level that we need in the multigrid method.
  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<8; ++cycle)
      {
        deallog << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_ball (triangulation);

            static const HyperBallBoundary<dim> boundary;
            triangulation.set_boundary (0, boundary);

            triangulation.refine_global (1);
          }
        else
          refine_grid ();

        deallog << "   Number of active cells:       "
                << triangulation.n_active_cells()
                << std::endl;

        setup_system ();

        assemble_system ();
        assemble_multigrid ();

        solve ();
        output_results (cycle);
      }
  }
}


// @sect3{The main() function}
//
// This is again the same function as in step-6:
int main ()
{
  try
    {
      using namespace Step16;

      LaplaceProblem<2> laplace_problem(1);
      laplace_problem.run ();
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
