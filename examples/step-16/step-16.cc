/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2003 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Authors: Guido Kanschat, University of Heidelberg, 2003
 *          Baerbel Janssen, University of Heidelberg, 2010
 *          Wolfgang Bangerth, Texas A&M University, 2010
 *          Timo Heister, Clemson University, 2018
 */


// @sect3{Include files}

// Again, the first few include files are already known, so we won't comment
// on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

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

// We will be using MeshWorker::mesh_loop to loop over the cells, so include it
// here:
#include <deal.II/meshworker/mesh_loop.h>


// This is C++:
#include <iostream>
#include <fstream>

namespace Step16
{
  using namespace dealii;

  // @sect3{The Scratch and Copy objects}
  //
  // We use MeshWorker::mesh_loop() to assemble our matrices. For this, we
  // need a ScratchData object to store temporary data on each cell (this is
  // just the FEValues object) and a CopyData object that will contain the
  // output of each cell assembly. For more details about the usage of scratch
  // and copy objects, see the WorkStream namespace.
  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim>       &mapping,
                const FiniteElement<dim> &fe,
                const unsigned int        quadrature_degree,
                const UpdateFlags         update_flags)
      : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
    {}

    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
    {}

    FEValues<dim> fe_values;
  };

  struct CopyData
  {
    unsigned int                         level;
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    template <class Iterator>
    void reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_active_or_mg_dof_indices(local_dof_indices);
      level = cell->level();
    }
  };

  // @sect3{The <code>LaplaceProblem</code> class template}

  // This main class is similar to the same class in step-6. As far as
  // member functions is concerned, the only additions are:
  // - The <code>assemble_multigrid</code> function that assembles the matrices
  // that correspond to the discrete operators on intermediate levels.
  // - The <code>cell_worker</code> function that assembles our PDE on a single
  // cell.
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem(const unsigned int degree);
    void run();

  private:
    template <class Iterator>
    void cell_worker(const Iterator   &cell,
                     ScratchData<dim> &scratch_data,
                     CopyData         &copy_data);

    void setup_system();
    void assemble_system();
    void assemble_multigrid();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

    Triangulation<dim> triangulation;
    const FE_Q<dim>    fe;
    DoFHandler<dim>    dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    AffineConstraints<double> constraints;

    Vector<double> solution;
    Vector<double> system_rhs;

    const unsigned int degree;

    // The following members are the essential data structures for the multigrid
    // method. The first four represent the sparsity patterns and the matrices
    // on individual levels of the multilevel hierarchy, very much like the
    // objects for the global mesh above.
    //
    // Then we have two new matrices only needed for multigrid methods with
    // local smoothing on adaptive meshes. They convey data between the interior
    // part of the refined region and the refinement edge, as outlined in detail
    // in the @ref mg_paper "multigrid paper".
    //
    // The last object stores information about the boundary indices on each
    // level and information about indices lying on a refinement edge between
    // two different refinement levels. It thus serves a similar purpose as
    // AffineConstraints, but on each level.
    MGLevelObject<SparsityPattern> mg_sparsity_patterns;
    MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns;

    MGLevelObject<SparseMatrix<double>> mg_matrices;
    MGLevelObject<SparseMatrix<double>> mg_interface_matrices;
    MGConstrainedDoFs                   mg_constrained_dofs;
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
  LaplaceProblem<dim>::LaplaceProblem(const unsigned int degree)
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    , fe(degree)
    , dof_handler(triangulation)
    , degree(degree)
  {}



  // @sect4{LaplaceProblem::setup_system}

  // In addition to just distributing the degrees of freedom in
  // the DoFHandler, we do the same on each level. Then, we follow the
  // same procedure as before to set up the system on the leaf mesh.
  template <int dim>
  void LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << " (by level: ";
    for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
      std::cout << dof_handler.n_dofs(level)
                << (level == triangulation.n_levels() - 1 ? ")" : ", ");
    std::cout << std::endl;


    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    std::set<types::boundary_id> dirichlet_boundary_ids = {0};
    Functions::ZeroFunction<dim> homogeneous_dirichlet_bc;
    const std::map<types::boundary_id, const Function<dim> *>
      dirichlet_boundary_functions = {
        {types::boundary_id(0), &homogeneous_dirichlet_bc}};
    VectorTools::interpolate_boundary_values(dof_handler,
                                             dirichlet_boundary_functions,
                                             constraints);
    constraints.close();

    {
      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
      sparsity_pattern.copy_from(dsp);
    }
    system_matrix.reinit(sparsity_pattern);

    // The multigrid constraints have to be initialized. They need to know
    // where Dirichlet boundary conditions are prescribed.
    mg_constrained_dofs.clear();
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                       dirichlet_boundary_ids);


    // Now for the things that concern the multigrid data structures. First, we
    // resize the multilevel objects to hold matrices and sparsity patterns for
    // every level. The coarse level is zero (this is mandatory right now but
    // may change in a future revision). Note that these functions take a
    // complete, inclusive range here (not a starting index and size), so the
    // finest level is <code>n_levels-1</code>. We first have to resize the
    // container holding the SparseMatrix classes, since they have to release
    // their SparsityPattern before they can be destroyed upon resizing.
    const unsigned int n_levels = triangulation.n_levels();

    mg_interface_matrices.resize(0, n_levels - 1);
    mg_matrices.resize(0, n_levels - 1);
    mg_sparsity_patterns.resize(0, n_levels - 1);
    mg_interface_sparsity_patterns.resize(0, n_levels - 1);

    // Now, we have to provide a matrix on each level. To this end, we first use
    // the MGTools::make_sparsity_pattern function to generate a preliminary
    // compressed sparsity pattern on each level (see the @ref Sparsity topic
    // for more information on this topic) and then copy it over to the one we
    // really want. The next step is to initialize the interface matrices with
    // the fitting sparsity pattern.
    //
    // It may be worth pointing out that the interface matrices only have
    // entries for degrees of freedom that sit at or next to the interface
    // between coarser and finer levels of the mesh. They are therefore even
    // sparser than the matrices on the individual levels of our multigrid
    // hierarchy. Therefore, we use a function specifically build for this
    // purpose to generate it.
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        {
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
                                     dof_handler.n_dofs(level));
          MGTools::make_sparsity_pattern(dof_handler, dsp, level);

          mg_sparsity_patterns[level].copy_from(dsp);
          mg_matrices[level].reinit(mg_sparsity_patterns[level]);
        }
        {
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
                                     dof_handler.n_dofs(level));
          MGTools::make_interface_sparsity_pattern(dof_handler,
                                                   mg_constrained_dofs,
                                                   dsp,
                                                   level);
          mg_interface_sparsity_patterns[level].copy_from(dsp);
          mg_interface_matrices[level].reinit(
            mg_interface_sparsity_patterns[level]);
        }
      }
  }


  // @sect4{LaplaceProblem::cell_worker}

  // The cell_worker function is used to assemble the matrix and right-hand side
  // on the given cell. This function is used for the active cells to generate
  // the system_matrix and on each level to build the level matrices.
  //
  // Note that we also assemble a right-hand side when called from
  // assemble_multigrid() even though it is not used.
  template <int dim>
  template <class Iterator>
  void LaplaceProblem<dim>::cell_worker(const Iterator   &cell,
                                        ScratchData<dim> &scratch_data,
                                        CopyData         &copy_data)
  {
    FEValues<dim> &fe_values = scratch_data.fe_values;
    fe_values.reinit(cell);

    const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell();
    const unsigned int n_q_points    = fe_values.get_quadrature().size();

    copy_data.reinit(cell, dofs_per_cell);

    const std::vector<double> &JxW = fe_values.get_JxW_values();

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const double coefficient =
          (fe_values.get_quadrature_points()[q][0] < 0.0) ? 1.0 : 0.1;

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) +=
                  coefficient *
                  (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)) *
                  JxW[q];
              }
            copy_data.cell_rhs(i) += 1.0 * fe_values.shape_value(i, q) * JxW[q];
          }
      }
  }



  // @sect4{LaplaceProblem::assemble_system}

  // The following function assembles the linear system on the active cells of
  // the mesh. For this, we pass two lambda functions to the mesh_loop()
  // function. The cell_worker function redirects to the class member function
  // of the same name, while the copier is specific to this function and copies
  // local matrix and vector to the corresponding global ones using the
  // constraints.
  template <int dim>
  void LaplaceProblem<dim>::assemble_system()
  {
    const MappingQ1<dim> mapping;

    auto cell_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          ScratchData<dim>                                     &scratch_data,
          CopyData                                             &copy_data) {
        this->cell_worker(cell, scratch_data, copy_data);
      };

    auto copier = [&](const CopyData &cd) {
      this->constraints.distribute_local_to_global(cd.cell_matrix,
                                                   cd.cell_rhs,
                                                   cd.local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
    };

    const unsigned int n_gauss_points = degree + 1;

    ScratchData<dim> scratch_data(mapping,
                                  fe,
                                  n_gauss_points,
                                  update_values | update_gradients |
                                    update_JxW_values |
                                    update_quadrature_points);

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          CopyData(),
                          MeshWorker::assemble_own_cells);
  }


  // @sect4{LaplaceProblem::assemble_multigrid}

  // The next function is the one that builds the matrices
  // that define the multigrid method on each level of the mesh. The integration
  // core is the same as above, but the loop below will go over all existing
  // cells instead of just the active ones, and the results must be entered into
  // the correct level matrices. Fortunately, MeshWorker hides most of that from
  // us, and thus the difference between this function and the previous lies
  // only in the setup of the assembler and the different iterators in the loop.
  //
  // We generate an AffineConstraints object for each level containing the
  // boundary and interface dofs as constrained entries. The corresponding
  // object is then used to generate the level matrices.
  template <int dim>
  void LaplaceProblem<dim>::assemble_multigrid()
  {
    const MappingQ1<dim> mapping;
    const unsigned int   n_levels = triangulation.n_levels();

    std::vector<AffineConstraints<double>> boundary_constraints(n_levels);
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        boundary_constraints[level].reinit(
          dof_handler.locally_owned_mg_dofs(level),
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level));

        for (const types::global_dof_index dof_index :
             mg_constrained_dofs.get_refinement_edge_indices(level))
          boundary_constraints[level].constrain_dof_to_zero(dof_index);
        for (const types::global_dof_index dof_index :
             mg_constrained_dofs.get_boundary_indices(level))
          boundary_constraints[level].constrain_dof_to_zero(dof_index);
        boundary_constraints[level].close();
      }

    auto cell_worker =
      [&](const typename DoFHandler<dim>::level_cell_iterator &cell,
          ScratchData<dim>                                    &scratch_data,
          CopyData                                            &copy_data) {
        this->cell_worker(cell, scratch_data, copy_data);
      };

    auto copier = [&](const CopyData &cd) {
      boundary_constraints[cd.level].distribute_local_to_global(
        cd.cell_matrix, cd.local_dof_indices, mg_matrices[cd.level]);

      const unsigned int dofs_per_cell = cd.local_dof_indices.size();

      // Interface entries are ignored by the boundary_constraints object
      // above when filling the mg_matrices[cd.level]. Instead, we copy these
      // entries into the interface matrix of the current level manually:
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          if (mg_constrained_dofs.is_interface_matrix_entry(
                cd.level, cd.local_dof_indices[i], cd.local_dof_indices[j]))
            {
              mg_interface_matrices[cd.level].add(cd.local_dof_indices[i],
                                                  cd.local_dof_indices[j],
                                                  cd.cell_matrix(i, j));
            }
    };

    const unsigned int n_gauss_points = degree + 1;

    ScratchData<dim> scratch_data(mapping,
                                  fe,
                                  n_gauss_points,
                                  update_values | update_gradients |
                                    update_JxW_values |
                                    update_quadrature_points);

    MeshWorker::mesh_loop(dof_handler.begin_mg(),
                          dof_handler.end_mg(),
                          cell_worker,
                          copier,
                          scratch_data,
                          CopyData(),
                          MeshWorker::assemble_own_cells);
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
  void LaplaceProblem<dim>::solve()
  {
    MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs);
    mg_transfer.build(dof_handler);

    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from(mg_matrices[0]);
    MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
    coarse_grid_solver.initialize(coarse_matrix);

    // The next component of a multilevel solver or preconditioner is that we
    // need a smoother on each level. A common choice for this is to use the
    // application of a relaxation method (such as the SOR, Jacobi or Richardson
    // method) or a small number of iterations of a solver method (such as CG or
    // GMRES). The mg::SmootherRelaxation and MGSmootherPrecondition classes
    // provide support for these two kinds of smoothers. Here, we opt for the
    // application of a single SOR iteration. To this end, we define an
    // appropriate alias and then setup a smoother object.
    //
    // The last step is to initialize the smoother object with our level
    // matrices and to set some smoothing parameters. The
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
    using Smoother = PreconditionSOR<SparseMatrix<double>>;
    mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother;
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
    mg::Matrix<Vector<double>> mg_matrix(mg_matrices);
    mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices);
    mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices);

    // Now, we are ready to set up the V-cycle operator and the multilevel
    // preconditioner.
    Multigrid<Vector<double>> mg(
      mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface_down, mg_interface_up);

    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
      preconditioner(dof_handler, mg, mg_transfer);

    // With all this together, we can finally get about solving the linear
    // system in the usual way:
    SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);

    solution = 0;

    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    std::cout << "   Number of CG iterations: " << solver_control.last_step()
              << '\n'
              << std::endl;
    constraints.distribute(solution);
  }



  // @sect4{Postprocessing}

  // The following two functions postprocess a solution once it is
  // computed. In particular, the first one refines the mesh at the beginning
  // of each cycle while the second one outputs results at the end of each
  // such cycle. The functions are almost unchanged from those in step-6.
  template <int dim>
  void LaplaceProblem<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(degree + 2),
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);
    triangulation.execute_coarsening_and_refinement();
  }



  template <int dim>
  void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
    data_out.write_vtk(output);
  }


  // @sect4{LaplaceProblem::run}

  // Like several of the functions above, this is almost exactly a copy of
  // the corresponding function in step-6. The only difference is the call to
  // <code>assemble_multigrid</code> that takes care of forming the matrices
  // on every level that we need in the multigrid method.
  template <int dim>
  void LaplaceProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 8; ++cycle)
      {
        std::cout << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_ball(triangulation);
            triangulation.refine_global(2);
          }
        else
          refine_grid();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl;

        setup_system();

        assemble_system();
        assemble_multigrid();

        solve();
        output_results(cycle);
      }
  }
} // namespace Step16


// @sect3{The main() function}
//
// This is again the same function as in step-6:
int main()
{
  try
    {
      using namespace Step16;

      LaplaceProblem<2> laplace_problem(1);
      laplace_problem.run();
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
