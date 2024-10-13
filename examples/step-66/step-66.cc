/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
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
 * Authors: Fabian Castelli, Karlsruhe Institute of Technology (KIT)
 */



// First we include the typical headers of the deal.II library needed for this
// tutorial:
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// In particular, we need to include the headers for the matrix-free framework:
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

// And since we want to use a geometric multigrid preconditioner, we need also
// the multilevel headers:
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>


// Finally some common C++ headers for in and output:
#include <fstream>
#include <iostream>



namespace Step66
{
  using namespace dealii;



  // @sect3{Matrix-free JacobianOperator}

  // In the beginning we define the matrix-free operator for the Jacobian. As a
  // guideline we follow the tutorials step-37 and step-48, where the precise
  // interface of the MatrixFreeOperators::Base class was extensively
  // documented.
  //
  // Since we want to use the Jacobian as system matrix and pass it to the
  // linear solver as well as to the multilevel preconditioner classes, we
  // derive the <code>JacobianOperator</code> class from the
  // MatrixFreeOperators::Base class, such that we have already the right
  // interface. The two functions we need to override from the base class are
  // the MatrixFreeOperators::Base::apply_add() and the
  // MatrixFreeOperators::Base::compute_diagonal() function. To allow
  // preconditioning with float precision we define the number type as template
  // argument.
  //
  // As mentioned already in the introduction, we need to evaluate the Jacobian
  // $F'$ at the last Newton step $u_h^n$ for the computation of the Newton
  // update $s_h^n$. To get the information of the last Newton step $u_h^n$ we
  // do pretty much the same as in step-37, where we stored the values of a
  // coefficient function in a table <code>nonlinear_values</code> once before
  // we use the matrix-free operator. Instead of a function
  // <code>evaluate_coefficient()</code>, we here implement a function
  // <code>evaluate_newton_step()</code>.
  //
  // As additional private member functions of the <code>JacobianOperator</code>
  // we implement the <code>local_apply()</code> and the
  // <code>local_compute_diagonal()</code> function. The first one is the actual
  // worker function for the matrix-vector application, which we pass to the
  // MatrixFree::cell_loop() in the <code>apply_add()</code> function. The later
  // one is the worker function to compute the diagonal, which we pass to the
  // MatrixFreeTools::compute_diagonal() function.
  //
  // For better readability of the source code we further define an alias for
  // the FEEvaluation object.
  template <int dim, int fe_degree, typename number>
  class JacobianOperator
    : public MatrixFreeOperators::
        Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    using value_type = number;

    using FECellIntegrator =
      FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number>;

    JacobianOperator();

    virtual void clear() override;

    void evaluate_newton_step(
      const LinearAlgebra::distributed::Vector<number> &newton_step);

    virtual void compute_diagonal() override;

  private:
    virtual void apply_add(
      LinearAlgebra::distributed::Vector<number>       &dst,
      const LinearAlgebra::distributed::Vector<number> &src) const override;

    void
    local_apply(const MatrixFree<dim, number>                    &data,
                LinearAlgebra::distributed::Vector<number>       &dst,
                const LinearAlgebra::distributed::Vector<number> &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

    void local_compute_diagonal(FECellIntegrator &integrator) const;

    Table<2, VectorizedArray<number>> nonlinear_values;
  };



  // The constructor of the <code>JacobianOperator</code> just calls the
  // constructor of the base class MatrixFreeOperators::Base, which is itself
  // derived from the EnableObserverPointer class.
  template <int dim, int fe_degree, typename number>
  JacobianOperator<dim, fe_degree, number>::JacobianOperator()
    : MatrixFreeOperators::Base<dim,
                                LinearAlgebra::distributed::Vector<number>>()
  {}



  // The <code>clear()</code> function resets the table holding the values for
  // the nonlinearity and call the <code>clear()</code> function of the base
  // class.
  template <int dim, int fe_degree, typename number>
  void JacobianOperator<dim, fe_degree, number>::clear()
  {
    nonlinear_values.reinit(0, 0);
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
      clear();
  }



  // @sect4{Evaluation of the old Newton step}

  // The following <code>evaluate_newton_step()</code> function is based on the
  // <code>evaluate_coefficient()</code> function from step-37. However, it does
  // not evaluate a function object, but evaluates a vector representing a
  // finite element function, namely the last Newton step needed for the
  // Jacobian. Therefore we set up a FEEvaluation object and evaluate the finite
  // element function in the quadrature points with the
  // FEEvaluation::read_dof_values_plain() and FEEvaluation::evaluate()
  // functions. We store the evaluated values of the finite element function
  // directly in the <code>nonlinear_values</code> table.
  //
  // This will work well and in the <code>local_apply()</code> function we can
  // use the values stored in the table to apply the matrix-vector product.
  // However, we can also optimize the implementation of the Jacobian at this
  // stage. We can directly evaluate the nonlinear function
  // <code>std::exp(newton_step[q])</code> and store these values in the table.
  // This skips all evaluations of the nonlinearity in each call of the
  // <code>vmult()</code> function.
  //
  // Note that we need to manually call the functions to exchange the ghost
  // data here, by calling
  // LinearAlgebra::distributed::Vector::update_ghost_values(), to ensure all
  // data from neighboring processes is available for evaluating the
  // finite-element interpolation on cells. In the other functions of this
  // tutorial program, MatrixFree::cell_loop() made sure to call this
  // function. Note that we clear the ghost state again at the end of the
  // function, in order to avoid mixing ghosted and non-ghosted vectors in
  // other parts of the solver.
  template <int dim, int fe_degree, typename number>
  void JacobianOperator<dim, fe_degree, number>::evaluate_newton_step(
    const LinearAlgebra::distributed::Vector<number> &newton_step)
  {
    const unsigned int n_cells = this->data->n_cell_batches();
    FECellIntegrator   phi(*this->data);

    newton_step.update_ghost_values();

    nonlinear_values.reinit(n_cells, phi.n_q_points);

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values_plain(newton_step);
        phi.evaluate(EvaluationFlags::values);

        for (const unsigned int q : phi.quadrature_point_indices())
          {
            nonlinear_values(cell, q) = std::exp(phi.get_value(q));
          }
      }
    newton_step.zero_out_ghost_values();
  }



  // @sect4{Nonlinear matrix-free operator application}

  // Now in the <code>local_apply()</code> function, which actually implements
  // the cell wise action of the system matrix, we can use the information of
  // the last Newton step stored in the table <code>nonlinear_values</code>. The
  // rest of this function is basically the same as in step-37. We set up the
  // FEEvaluation object, gather and evaluate the values and gradients of the
  // input vector <code>src</code>, submit the values and gradients according to
  // the form of the Jacobian and finally call FEEvaluation::integrate_scatter()
  // to perform the cell integration and distribute the local contributions into
  // the global vector <code> dst</code>.
  template <int dim, int fe_degree, typename number>
  void JacobianOperator<dim, fe_degree, number>::local_apply(
    const MatrixFree<dim, number>                    &data,
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        AssertDimension(nonlinear_values.size(0),
                        phi.get_matrix_free().n_cell_batches());
        AssertDimension(nonlinear_values.size(1), phi.n_q_points);


        phi.reinit(cell);

        phi.gather_evaluate(src,
                            EvaluationFlags::values |
                              EvaluationFlags::gradients);

        for (const unsigned int q : phi.quadrature_point_indices())
          {
            phi.submit_value(-nonlinear_values(cell, q) * phi.get_value(q), q);
            phi.submit_gradient(phi.get_gradient(q), q);
          }

        phi.integrate_scatter(EvaluationFlags::values |
                                EvaluationFlags::gradients,
                              dst);
      }
  }



  // Next we use MatrixFree::cell_loop() to perform the actual loop over all
  // cells computing the cell contribution to the matrix-vector product.
  template <int dim, int fe_degree, typename number>
  void JacobianOperator<dim, fe_degree, number>::apply_add(
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src) const
  {
    this->data->cell_loop(&JacobianOperator::local_apply, this, dst, src);
  }



  // @sect4{Diagonal of the JacobianOperator}

  // The internal worker function <code>local_compute_diagonal()</code> for the
  // computation of the diagonal is similar to the above worker function
  // <code>local_apply()</code>. However, as major difference we do not read
  // values from a input vector or distribute any local results to an output
  // vector. Instead the only input argument is the used FEEvaluation object.
  template <int dim, int fe_degree, typename number>
  void JacobianOperator<dim, fe_degree, number>::local_compute_diagonal(
    FECellIntegrator &phi) const
  {
    AssertDimension(nonlinear_values.size(0),
                    phi.get_matrix_free().n_cell_batches());
    AssertDimension(nonlinear_values.size(1), phi.n_q_points);

    const unsigned int cell = phi.get_current_cell_index();

    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    for (const unsigned int q : phi.quadrature_point_indices())
      {
        phi.submit_value(-nonlinear_values(cell, q) * phi.get_value(q), q);
        phi.submit_gradient(phi.get_gradient(q), q);
      }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }



  // Finally we override the MatrixFreeOperators::Base::compute_diagonal()
  // function of the base class of the <code>JacobianOperator</code>. Although
  // the name of the function suggests just the computation of the diagonal,
  // this function does a bit more. Because we only really need the inverse of
  // the matrix diagonal elements for the Chebyshev smoother of the multigrid
  // preconditioner, we compute the diagonal and store the inverse elements.
  // Therefore we first initialize the <code>inverse_diagonal_entries</code>.
  // Then we compute the diagonal by passing the worker function
  // <code>local_compute_diagonal()</code> to the
  // MatrixFreeTools::compute_diagonal() function. In the end we loop over the
  // diagonal and invert the elements by hand. Note, that during this loop we
  // catch the constrained DOFs and set them manually to one.
  template <int dim, int fe_degree, typename number>
  void JacobianOperator<dim, fe_degree, number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);

    MatrixFreeTools::compute_diagonal(*this->data,
                                      inverse_diagonal,
                                      &JacobianOperator::local_compute_diagonal,
                                      this);

    for (auto &diagonal_element : inverse_diagonal)
      {
        diagonal_element = (std::abs(diagonal_element) > 1.0e-10) ?
                             (1.0 / diagonal_element) :
                             1.0;
      }
  }



  // @sect3{GelfandProblem class}

  // After implementing the matrix-free operators we can now define the solver
  // class for the <i>Gelfand problem</i>. This class is based on the common
  // structure of all previous tutorial programs, in particular it is based on
  // step-15, solving also a nonlinear problem. Since we are using the
  // matrix-free framework, we no longer need an assemble_system function any
  // more, instead the information of the matrix is rebuilt in every call of the
  // <code>vmult()</code> function. However, for the application of the Newton
  // scheme we need to assemble the right hand side of the linearized problems
  // and compute the residuals. Therefore, we implement an additional function
  // <code>evaluate_residual()</code>, which we later call in the
  // <code>assemble_rhs()</code> and the <code>compute_residual()</code>
  // function. Finally, the typical <code>solve()</code> function here
  // implements the Newton method, whereas the solution of the linearized system
  // is computed in the function <code>compute_update()</code>. As the
  // MatrixFree framework handles the polynomial degree of the Lagrangian finite
  // element method as a template parameter, we declare it also as a template
  // parameter for the problem solver class.
  template <int dim, int fe_degree>
  class GelfandProblem
  {
  public:
    GelfandProblem();

    void run();

  private:
    void make_grid();

    void setup_system();

    void evaluate_residual(
      LinearAlgebra::distributed::Vector<double>       &dst,
      const LinearAlgebra::distributed::Vector<double> &src) const;

    void local_evaluate_residual(
      const MatrixFree<dim, double>                    &data,
      LinearAlgebra::distributed::Vector<double>       &dst,
      const LinearAlgebra::distributed::Vector<double> &src,
      const std::pair<unsigned int, unsigned int>      &cell_range) const;

    void assemble_rhs();

    double compute_residual(const double alpha);

    void compute_update();

    void solve();

    double compute_solution_norm() const;

    void output_results(const unsigned int cycle) const;


    // For the parallel computation we define a
    // parallel::distributed::Triangulation. As the computational domain is a
    // circle in 2d and a ball in 3d, we assign in addition to the
    // SphericalManifold for boundary cells a TransfiniteInterpolationManifold
    // object for the mapping of the inner cells, which takes care of the inner
    // cells. In this example we use an isoparametric finite element approach
    // and thus use the MappingQ class. For further details you may read the
    // detailed description of this class.
    parallel::distributed::Triangulation<dim> triangulation;
    const MappingQ<dim>                       mapping;


    // As usual we then define the Lagrangian finite elements FE_Q and a
    // DoFHandler.
    const FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;


    // For the linearized discrete system we define an AffineConstraints objects
    // and the <code>system_matrix</code>, which is in this example represented
    // as a matrix-free operator.
    AffineConstraints<double> constraints;
    using SystemMatrixType = JacobianOperator<dim, fe_degree, double>;
    SystemMatrixType system_matrix;


    // The multilevel object is also based on the matrix-free operator for the
    // Jacobian. Since we need to evaluate the Jacobian with the last Newton
    // step, we also need to evaluate the level operator with the last Newton
    // step for the preconditioner. Thus in addition to
    // <code>mg_matrices</code>, we also need a MGLevelObject to store the
    // interpolated solution vector on each level. As in step-37 we use float
    // precision for the preconditioner. Moreover, we define the
    // MGTransferMatrixFree object as a class variable, since we need to set it
    // up only once when the triangulation has changed and can then use it again
    // in each Newton step.
    MGConstrainedDoFs mg_constrained_dofs;
    using LevelMatrixType = JacobianOperator<dim, fe_degree, float>;
    MGLevelObject<LevelMatrixType>                           mg_matrices;
    MGLevelObject<LinearAlgebra::distributed::Vector<float>> mg_solution;
    MGTransferMatrixFree<dim, float>                         mg_transfer;


    // Of course we also need vectors holding the <code>solution</code>, the
    // <code>newton_update</code> and the <code>system_rhs</code>. In that way
    // we can always store the last Newton step in the solution vector and just
    // add the update to get the next Newton step.
    LinearAlgebra::distributed::Vector<double> solution;
    LinearAlgebra::distributed::Vector<double> newton_update;
    LinearAlgebra::distributed::Vector<double> system_rhs;


    // Finally we have a variable for the number of iterations of the linear
    // solver.
    unsigned int linear_iterations;


    // For the output in programs running in parallel with MPI, we use the
    // ConditionalOStream class to avoid multiple output of the same data by
    // different MPI ranks.
    ConditionalOStream pcout;


    // Finally for the time measurement we use a TimerOutput object, which
    // prints the elapsed CPU and wall times for each function in a nicely
    // formatted table after the program has finished.
    TimerOutput computing_timer;
  };



  // The constructor of the <code>GelfandProblem</code> initializes the class
  // variables. In particular, we set up the multilevel support for the
  // parallel::distributed::Triangulation, set the mapping degree equal to the
  // finite element degree, initialize the ConditionalOStream and tell the
  // TimerOutput that we want to see the wall times only on demand.
  template <int dim, int fe_degree>
  GelfandProblem<dim, fe_degree>::GelfandProblem()
    : triangulation(MPI_COMM_WORLD,
                    Triangulation<dim>::limit_level_difference_at_vertices,
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
    , mapping(fe_degree)
    , fe(fe_degree)
    , dof_handler(triangulation)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    , computing_timer(MPI_COMM_WORLD,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
  {}



  // @sect4{GelfandProblem::make_grid}

  // As the computational domain we use the <code>dim</code>-dimensional unit
  // ball. We follow the instructions for the TransfiniteInterpolationManifold
  // class and also assign a SphericalManifold for the boundary. Finally, we
  // refine the initial mesh 3 - <code>dim</code> times globally.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::make_grid()
  {
    TimerOutput::Scope t(computing_timer, "make grid");

    SphericalManifold<dim>                boundary_manifold;
    TransfiniteInterpolationManifold<dim> inner_manifold;

    GridGenerator::hyper_ball(triangulation);

    triangulation.set_all_manifold_ids(1);
    triangulation.set_all_manifold_ids_on_boundary(0);

    triangulation.set_manifold(0, boundary_manifold);

    inner_manifold.initialize(triangulation);
    triangulation.set_manifold(1, inner_manifold);

    triangulation.refine_global(3 - dim);
  }



  // @sect4{GelfandProblem::setup_system}

  // The <code>setup_system()</code> function is quasi identical to the one in
  // step-37. The only differences are obviously the time measurement with only
  // one TimerOutput::Scope instead of measuring each part individually, and
  // more importantly the initialization of the MGLevelObject for the
  // interpolated solution vector of the previous Newton step. Another important
  // change is the setup of the MGTransferMatrixFree object, which we can reuse
  // in each Newton step as the <code>triangulation</code> will not be not
  // changed.
  //
  // Note how we can use the same MatrixFree object twice, for the
  // <code>JacobianOperator</code> and the multigrid preconditioner.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup system");

    system_matrix.clear();
    mg_matrices.clear_elements();

    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    constraints.clear();
    constraints.reinit(dof_handler.locally_owned_dofs(),
                       DoFTools::extract_locally_relevant_dofs(dof_handler));
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    {
      typename MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::partition_color;
      additional_data.mapping_update_flags =
        (update_values | update_gradients | update_JxW_values |
         update_quadrature_points);
      auto system_mf_storage = std::make_shared<MatrixFree<dim, double>>();
      system_mf_storage->reinit(mapping,
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);

      system_matrix.initialize(system_mf_storage);
    }

    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(newton_update);
    system_matrix.initialize_dof_vector(system_rhs);


    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels - 1);
    mg_solution.resize(0, nlevels - 1);

    const std::set<types::boundary_id> dirichlet_boundary_ids = {0};
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                       dirichlet_boundary_ids);

    mg_transfer.initialize_constraints(mg_constrained_dofs);
    mg_transfer.build(dof_handler);

    for (unsigned int level = 0; level < nlevels; ++level)
      {
        AffineConstraints<double> level_constraints(
          dof_handler.locally_owned_mg_dofs(level),
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level));

        for (const types::global_dof_index dof_index :
             mg_constrained_dofs.get_boundary_indices(level))
          level_constraints.constrain_dof_to_zero(dof_index);
        level_constraints.close();

        typename MatrixFree<dim, float>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, float>::AdditionalData::partition_color;
        additional_data.mapping_update_flags =
          (update_values | update_gradients | update_JxW_values |
           update_quadrature_points);
        additional_data.mg_level = level;
        auto mg_mf_storage_level = std::make_shared<MatrixFree<dim, float>>();
        mg_mf_storage_level->reinit(mapping,
                                    dof_handler,
                                    level_constraints,
                                    QGauss<1>(fe.degree + 1),
                                    additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level,
                                      mg_constrained_dofs,
                                      level);
        mg_matrices[level].initialize_dof_vector(mg_solution[level]);
      }
  }



  // @sect4{GelfandProblem::evaluate_residual}

  // Next we implement a function which evaluates the nonlinear discrete
  // residual for a given input vector ($\texttt{dst} = F(\texttt{src})$). This
  // function is then used for the assembly of the right hand side of the
  // linearized system and later for the computation of the residual of the next
  // Newton step to check if we already reached the error tolerance. As this
  // function should not affect any class variable we define it as a constant
  // function. Internally we exploit the fast finite element evaluation through
  // the FEEvaluation class and the MatrixFree::cell_loop(), similar to
  // <code>apply_add()</code> function of the <code>JacobianOperator</code>.
  //
  // First we create a pointer to the MatrixFree object, which is stored in the
  // <code>system_matrix</code>. Then we pass the worker function
  // <code>local_evaluate_residual()</code> for the cell wise evaluation of the
  // residual together with the input and output vector to the
  // MatrixFree::cell_loop(). In addition, we enable the zero out of the output
  // vector in the loop, which is more efficient than calling <code>dst =
  // 0.0</code> separately before.
  //
  // Note that with this approach we do not have to take care about the MPI
  // related data exchange, since all the bookkeeping is done by the
  // MatrixFree::cell_loop().
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::evaluate_residual(
    LinearAlgebra::distributed::Vector<double>       &dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    auto matrix_free = system_matrix.get_matrix_free();

    matrix_free->cell_loop(
      &GelfandProblem::local_evaluate_residual, this, dst, src, true);
  }



  // @sect4{GelfandProblem::local_evaluate_residual}

  // This is the internal worker function for the evaluation of the residual.
  // Essentially it has the same structure as the <code>local_apply()</code>
  // function of the <code>JacobianOperator</code> and evaluates the residual
  // for the input vector <code>src</code> on the given set of cells
  // <code>cell_range</code>. The difference to the above mentioned
  // <code>local_apply()</code> function is, that we split the
  // FEEvaluation::gather_evaluate() function into
  // FEEvaluation::read_dof_values_plain() and FEEvaluation::evaluate(), since
  // the input vector might have constrained DOFs.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::local_evaluate_residual(
    const MatrixFree<dim, double>                    &data,
    LinearAlgebra::distributed::Vector<double>       &dst,
    const LinearAlgebra::distributed::Vector<double> &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        phi.read_dof_values_plain(src);
        phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        for (const unsigned int q : phi.quadrature_point_indices())
          {
            phi.submit_value(-std::exp(phi.get_value(q)), q);
            phi.submit_gradient(phi.get_gradient(q), q);
          }

        phi.integrate_scatter(EvaluationFlags::values |
                                EvaluationFlags::gradients,
                              dst);
      }
  }



  // @sect4{GelfandProblem::assemble_rhs}

  // Using the above function <code>evaluate_residual()</code> to evaluate the
  // nonlinear residual, the assembly of the right hand side of the linearized
  // system becomes now a very easy task. We just call the
  // <code>evaluate_residual()</code> function and multiply the result with
  // minus one.
  //
  // Experiences show that using the FEEvaluation class is much faster than a
  // classical implementation with FEValues and co.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::assemble_rhs()
  {
    TimerOutput::Scope t(computing_timer, "assemble right hand side");

    evaluate_residual(system_rhs, solution);

    system_rhs *= -1.0;
  }



  // @sect4{GelfandProblem::compute_residual}

  // According to step-15 the following function computes the norm of the
  // nonlinear residual for the solution $u_h^n + \alpha s_h^n$ with the help of
  // the <code>evaluate_residual()</code> function. The Newton step length
  // $\alpha$ becomes important if we would use an adaptive version of the
  // Newton method. Then for example we would compute the residual for different
  // step lengths and compare the residuals. However, for our problem the full
  // Newton step with $\alpha=1$ is the best we can do. An adaptive version of
  // Newton's method becomes interesting if we have no good initial value. Note
  // that in theory Newton's method converges with quadratic order, but only if
  // we have an appropriate initial value. For unsuitable initial values the
  // Newton method diverges even with quadratic order. A common way is then to
  // use a damped version $\alpha<1$ until the Newton step is good enough and
  // the full Newton step can be performed. This was also discussed in step-15.
  template <int dim, int fe_degree>
  double GelfandProblem<dim, fe_degree>::compute_residual(const double alpha)
  {
    TimerOutput::Scope t(computing_timer, "compute residual");

    LinearAlgebra::distributed::Vector<double> residual;
    LinearAlgebra::distributed::Vector<double> evaluation_point;

    system_matrix.initialize_dof_vector(residual);
    system_matrix.initialize_dof_vector(evaluation_point);

    evaluation_point = solution;
    if (alpha > 1e-12)
      {
        evaluation_point.add(alpha, newton_update);
      }

    evaluate_residual(residual, evaluation_point);

    return residual.l2_norm();
  }



  // @sect4{GelfandProblem::compute_update}

  // In order to compute the Newton updates in each Newton step we solve the
  // linear system with the CG algorithm together with a geometric multigrid
  // preconditioner. For this we first set up the PreconditionMG object with a
  // Chebyshev smoother like we did in step-37.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::compute_update()
  {
    TimerOutput::Scope t(computing_timer, "compute update");

    // We remember that the Jacobian depends on the last Newton step stored in
    // the solution vector. So we update the ghost values of the Newton step and
    // pass it to the <code>JacobianOperator</code> to store the information.
    solution.update_ghost_values();

    system_matrix.evaluate_newton_step(solution);


    // Next we also have to pass the last Newton step to the multilevel
    // operators. Therefore, we need to interpolate the Newton step to all
    // levels of the triangulation. This is done with the
    // MGTransferMatrixFree::interpolate_to_mg().
    mg_transfer.interpolate_to_mg(dof_handler, mg_solution, solution);


    // Now we can set up the preconditioner. We define the smoother and pass the
    // interpolated vectors of the Newton step to the multilevel operators.
    using SmootherType =
      PreconditionChebyshev<LevelMatrixType,
                            LinearAlgebra::distributed::Vector<float>>;
    mg::SmootherRelaxation<SmootherType,
                           LinearAlgebra::distributed::Vector<float>>
                                                         mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        if (level > 0)
          {
            smoother_data[level].smoothing_range     = 15.;
            smoother_data[level].degree              = 4;
            smoother_data[level].eig_cg_n_iterations = 10;
          }
        else
          {
            smoother_data[0].smoothing_range = 1e-3;
            smoother_data[0].degree          = numbers::invalid_unsigned_int;
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
          }

        mg_matrices[level].evaluate_newton_step(mg_solution[level]);
        mg_matrices[level].compute_diagonal();

        smoother_data[level].preconditioner =
          mg_matrices[level].get_matrix_diagonal_inverse();
      }
    mg_smoother.initialize(mg_matrices, smoother_data);

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
      mg_coarse;
    mg_coarse.initialize(mg_smoother);

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(
      mg_matrices);

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        mg_interface_matrices[level].initialize(mg_matrices[level]);
      }
    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface(
      mg_interface_matrices);

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);

    PreconditionMG<dim,
                   LinearAlgebra::distributed::Vector<float>,
                   MGTransferMatrixFree<dim, float>>
      preconditioner(dof_handler, mg, mg_transfer);


    // Finally we set up the SolverControl and the SolverCG to solve the
    // linearized problem for the current Newton update. An important fact of
    // the implementation of SolverCG or also SolverGMRES is, that the vector
    // holding the solution of the linear system (here
    // <code>newton_update</code>) can be used to pass a starting value. In
    // order to start the iterative solver always with a zero vector we reset
    // the <code>newton_update</code> explicitly before calling
    // SolverCG::solve(). Afterwards we distribute the Dirichlet boundary
    // conditions stored in <code>constraints</code> and store the number of
    // iteration steps for the later output.
    SolverControl solver_control(100, 1.e-12);
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);

    newton_update = 0.0;

    cg.solve(system_matrix, newton_update, system_rhs, preconditioner);

    constraints.distribute(newton_update);

    linear_iterations = solver_control.last_step();


    // Then for bookkeeping we zero out the ghost values.
    solution.zero_out_ghost_values();
  }



  // @sect4{GelfandProblem::solve}

  // Now we implement the actual Newton solver for the nonlinear problem.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::solve()
  {
    TimerOutput::Scope t(computing_timer, "solve");


    // We define a maximal number of Newton steps and tolerances for the
    // convergence criterion. Usually, with good starting values, the Newton
    // method converges in three to six steps, so maximal ten steps should be
    // totally sufficient. As tolerances we use $\|F(u^n_h)\|<\text{TOL}_f =
    // 10^{-12}$ for the norm of the residual and $\|s_h^n\| < \text{TOL}_x =
    // 10^{-10}$ for the norm of the Newton update. This seems a bit over the
    // top, but we will see that, for our example, we will achieve these
    // tolerances after a few steps.
    const unsigned int itmax = 10;
    const double       TOLf  = 1e-12;
    const double       TOLx  = 1e-10;


    Timer solver_timer;
    solver_timer.start();


    // Now we start the actual Newton iteration.
    for (unsigned int newton_step = 1; newton_step <= itmax; ++newton_step)
      {
        // We assemble the right hand side of the linearized problem and compute
        // the Newton update.
        assemble_rhs();
        compute_update();


        // Then we compute the errors, namely the norm of the Newton update and
        // the residual. Note that at this point one could incorporate a step
        // size control for the Newton method by varying the input parameter
        // $\alpha$ for the compute_residual function. However, here we just use
        // $\alpha$ equal to one for a plain Newton iteration.
        const double ERRx = newton_update.l2_norm();
        const double ERRf = compute_residual(1.0);


        // Next we advance the Newton step by adding the Newton update to the
        // current Newton step.
        solution.add(1.0, newton_update);


        // A short output will inform us on the current Newton step.
        pcout << "   Nstep " << newton_step << ", errf = " << ERRf
              << ", errx = " << ERRx << ", it = " << linear_iterations
              << std::endl;


        // After each Newton step we check the convergence criteria. If at least
        // one of those is fulfilled we are done and end the loop. If we haven't
        // found a satisfying solution after the maximal amount of Newton
        // iterations, we inform the user about this shortcoming.
        if (ERRf < TOLf || ERRx < TOLx)
          {
            solver_timer.stop();

            pcout << "Convergence step " << newton_step << " value " << ERRf
                  << " (used wall time: " << solver_timer.wall_time() << " s)"
                  << std::endl;

            break;
          }
        else if (newton_step == itmax)
          {
            solver_timer.stop();
            pcout << "WARNING: No convergence of Newton's method after "
                  << newton_step << " steps." << std::endl;

            break;
          }
      }
  }



  // @sect4{GelfandProblem::compute_solution_norm}

  // The computation of the H1-seminorm of the solution can be done in the same
  // way as in step-59. We update the ghost values and use the function
  // VectorTools::integrate_difference(). In the end we gather all computations
  // from all MPI ranks and return the norm.
  template <int dim, int fe_degree>
  double GelfandProblem<dim, fe_degree>::compute_solution_norm() const
  {
    solution.update_ghost_values();

    Vector<float> norm_per_cell(triangulation.n_active_cells());

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Functions::ZeroFunction<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::H1_seminorm);

    solution.zero_out_ghost_values();

    return VectorTools::compute_global_error(triangulation,
                                             norm_per_cell,
                                             VectorTools::H1_seminorm);
  }



  // @sect4{GelfandProblem::output_results}

  // We generate the graphical output files in vtu format together with a pvtu
  // master file at once by calling the DataOut::write_vtu_with_pvtu_record()
  // function in the same way as in step-37. In addition, as in step-40, we
  // query the types::subdomain_id of each cell and write the distribution of
  // the triangulation among the MPI ranks into the output file. Finally, we
  // generate the patches of the solution by calling DataOut::build_patches().
  // However, since we have a computational domain with a curved boundary, we
  // additionally pass the <code>mapping</code> and the finite element degree as
  // number of subdivision. But this is still not enough for the correct
  // representation of the solution, for example in ParaView, because we
  // attached a TransfiniteInterpolationManifold to the inner cells, which
  // results in curved cells in the interior. Therefore we pass as third
  // argument the DataOut::curved_inner_cells option, such that also the inner
  // cells use the corresponding manifold description to build the patches.
  //
  // Note that we could handle the higher order elements with the flag
  // DataOutBase::VtkFlags::write_higher_order_cells. However, due to the
  // limited compatibility to previous version of ParaView and the missing
  // support by VisIt, we left this option for a future version.
  template <int dim, int fe_degree>
  void
  GelfandProblem<dim, fe_degree>::output_results(const unsigned int cycle) const
  {
    if (triangulation.n_global_active_cells() > 1e6)
      return;

    solution.update_ghost_values();

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      {
        subdomain(i) = triangulation.locally_owned_subdomain();
      }
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(mapping,
                           fe.degree,
                           DataOut<dim>::curved_inner_cells);

    DataOutBase::VtkFlags flags;
    flags.compression_level = DataOutBase::CompressionLevel::best_speed;
    data_out.set_flags(flags);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution_" + std::to_string(dim) + "d", cycle, MPI_COMM_WORLD, 3);

    solution.zero_out_ghost_values();
  }



  // @sect4{GelfandProblem::run}

  // The last missing function of the solver class for the <i>Gelfand
  // problem</i> is the run function. In the beginning we print information
  // about the system specifications and the finite element space we use. The
  // problem is solved several times on a successively refined mesh.
  template <int dim, int fe_degree>
  void GelfandProblem<dim, fe_degree>::run()
  {
    {
      const unsigned int n_ranks =
        Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      const unsigned int n_vect_doubles = VectorizedArray<double>::size();
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

      std::string DAT_header = "START DATE: " + Utilities::System::get_date() +
                               ", TIME: " + Utilities::System::get_time();
      std::string MPI_header = "Running with " + std::to_string(n_ranks) +
                               " MPI process" + (n_ranks > 1 ? "es" : "");
      std::string VEC_header =
        "Vectorization over " + std::to_string(n_vect_doubles) +
        " doubles = " + std::to_string(n_vect_bits) + " bits (" +
        Utilities::System::get_current_vectorization_level() +
        "), VECTORIZATION_LEVEL=" +
        std::to_string(DEAL_II_COMPILER_VECTORIZATION_LEVEL);
      std::string SOL_header = "Finite element space: " + fe.get_name();

      pcout << std::string(80, '=') << std::endl;
      pcout << DAT_header << std::endl;
      pcout << std::string(80, '-') << std::endl;

      pcout << MPI_header << std::endl;
      pcout << VEC_header << std::endl;
      pcout << SOL_header << std::endl;

      pcout << std::string(80, '=') << std::endl;
    }


    for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle)
      {
        pcout << std::string(80, '-') << std::endl;
        pcout << "Cycle " << cycle << std::endl;
        pcout << std::string(80, '-') << std::endl;


        // The first task in actually solving the problem is to generate or
        // refine the triangulation.
        if (cycle == 0)
          {
            make_grid();
          }
        else
          {
            triangulation.refine_global(1);
          }


        // Now we set up the system and solve the problem. These steps are
        // accompanied by time measurement and textual output.
        Timer timer;

        pcout << "Set up system..." << std::endl;
        setup_system();

        pcout << "   Triangulation: " << triangulation.n_global_active_cells()
              << " cells" << std::endl;
        pcout << "   DoFHandler:    " << dof_handler.n_dofs() << " DoFs"
              << std::endl;
        pcout << std::endl;


        pcout << "Solve using Newton's method..." << std::endl;
        solve();
        pcout << std::endl;


        timer.stop();
        pcout << "Time for setup+solve (CPU/Wall) " << timer.cpu_time() << '/'
              << timer.wall_time() << " s" << std::endl;
        pcout << std::endl;


        // After the problem was solved we compute the norm of the solution and
        // generate the graphical output files.
        pcout << "Output results..." << std::endl;
        const double norm = compute_solution_norm();
        output_results(cycle);

        pcout << "  H1 seminorm: " << norm << std::endl;
        pcout << std::endl;


        // Finally after each cycle we print the timing information.
        computing_timer.print_summary();
        computing_timer.reset();
      }
  }
} // namespace Step66



// @sect3{The <code>main</code> function}

// As typical for programs running in parallel with MPI we set up the MPI
// framework and disable shared-memory parallelization by limiting the number of
// threads to one. Finally to run the solver for the <i>Gelfand problem</i> we
// create an object of the <code>GelfandProblem</code> class and call the run
// function. Exemplarily we solve the problem once in 2d and once in 3d each
// with fourth-order Lagrangian finite elements.
int main(int argc, char *argv[])
{
  try
    {
      using namespace Step66;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      {
        GelfandProblem<2, 4> gelfand_problem;
        gelfand_problem.run();
      }

      {
        GelfandProblem<3, 4> gelfand_problem;
        gelfand_problem.run();
      }
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
