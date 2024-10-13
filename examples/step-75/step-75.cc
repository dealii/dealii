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
 * Authors: Marc Fehling, Colorado State University, 2021
 *          Peter Munch, Technical University of Munich and Helmholtz-Zentrum
 *                       hereon, 2021
 *          Wolfgang Bangerth, Colorado State University, 2021
 */


// @sect3{Include files}
//
// The following include files have been used and discussed in previous tutorial
// programs, especially in step-27 and step-40.
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <fstream>
#include <iostream>

// For load balancing we will assign individual weights on cells, and for that
// we will use the class parallel::CellWeights.
#include <deal.II/distributed/cell_weights.h>

// The solution function requires a transformation from Cartesian to polar
// coordinates. The GeometricUtilities::Coordinates namespace provides the
// necessary tools.
#include <deal.II/base/function.h>
#include <deal.II/base/geometric_utilities.h>

// The following include files will enable the MatrixFree functionality.
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/tools.h>

// We will use LinearAlgebra::distributed::Vector for linear algebra operations.
#include <deal.II/lac/la_parallel_vector.h>

// We are left to include the files needed by the multigrid solver.
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

namespace Step75
{
  using namespace dealii;

  // @sect3{The <code>Solution</code> class template}

  // We have an analytic solution to the scenario at our disposal. We will use
  // this solution to impose boundary conditions for the numerical solution of
  // the problem. The formulation of the solution requires a transformation to
  // polar coordinates. To transform from Cartesian to spherical coordinates, we
  // will use a helper function from the GeometricUtilities::Coordinates
  // namespace. The first two coordinates of this transformation correspond to
  // polar coordinates in the x-y-plane.
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution()
      : Function<dim>()
    {
      Assert(dim > 1, ExcNotImplemented());
    }

    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/) const override
    {
      const std::array<double, 2> polar =
        GeometricUtilities::Coordinates::to_spherical(Point<2>(p[0], p[1]));

      constexpr const double alpha = 2. / 3.;
      return std::pow(polar[0], alpha) * std::sin(alpha * polar[1]);
    }
  };



  // @sect3{Parameters}

  // For this tutorial, we will use a simplified set of parameters. It is also
  // possible to use a ParameterHandler class here, but to keep this tutorial
  // short we decided on using simple structs. The actual intention of all these
  // parameters will be described in the upcoming classes at their respective
  // location where they are used.
  //
  // The following parameter set controls the coarse-grid solver, the smoothers,
  // and the inter-grid transfer scheme of the multigrid mechanism.
  // We populate it with default parameters.
  struct MultigridParameters
  {
    struct
    {
      std::string  type            = "cg_with_amg";
      unsigned int maxiter         = 10000;
      double       abstol          = 1e-20;
      double       reltol          = 1e-4;
      unsigned int smoother_sweeps = 1;
      unsigned int n_cycles        = 1;
      std::string  smoother_type   = "ILU";
    } coarse_solver;

    struct
    {
      std::string  type                = "chebyshev";
      double       smoothing_range     = 20;
      unsigned int degree              = 5;
      unsigned int eig_cg_n_iterations = 20;
    } smoother;

    struct
    {
      MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType
        p_sequence = MGTransferGlobalCoarseningTools::
          PolynomialCoarseningSequenceType::decrease_by_one;
      bool perform_h_transfer = true;
    } transfer;
  };



  // This is the general parameter struct for the problem class. You will find
  // this struct divided into several categories, including general runtime
  // parameters, level limits, refine and coarsen fractions, as well as
  // parameters for cell weighting. It also contains an instance of the above
  // struct for multigrid parameters which will be passed to the multigrid
  // algorithm.
  struct Parameters
  {
    unsigned int n_cycles         = 8;
    double       tolerance_factor = 1e-12;

    MultigridParameters mg_data;

    unsigned int min_h_level            = 5;
    unsigned int max_h_level            = 12;
    unsigned int min_p_degree           = 2;
    unsigned int max_p_degree           = 6;
    unsigned int max_p_level_difference = 1;

    double refine_fraction    = 0.3;
    double coarsen_fraction   = 0.03;
    double p_refine_fraction  = 0.9;
    double p_coarsen_fraction = 0.9;

    double weighting_factor   = 1.;
    double weighting_exponent = 1.;
  };



  // @sect3{Matrix-free Laplace operator}

  // This is a matrix-free implementation of the Laplace operator that will
  // basically take over the part of the `assemble_system()` function from other
  // tutorials. The meaning of all member functions will be explained at their
  // definition later.
  //
  // We will use the FEEvaluation class to evaluate the solution vector
  // at the quadrature points and to perform the integration. In contrast to
  // other tutorials, the template arguments `degree` is set to $-1$ and
  // `number of quadrature in 1d` to $0$. In this case, FEEvaluation selects
  // dynamically the correct polynomial degree and number of quadrature
  // points. Here, we introduce an alias to FEEvaluation with the correct
  // template parameters so that we do not have to worry about them later on.
  template <int dim, typename number>
  class LaplaceOperator : public EnableObserverPointer
  {
  public:
    using VectorType = LinearAlgebra::distributed::Vector<number>;

    using FECellIntegrator = FEEvaluation<dim, -1, 0, 1, number>;

    LaplaceOperator() = default;

    LaplaceOperator(const hp::MappingCollection<dim> &mapping,
                    const DoFHandler<dim>            &dof_handler,
                    const hp::QCollection<dim>       &quad,
                    const AffineConstraints<number>  &constraints,
                    VectorType                       &system_rhs);

    void reinit(const hp::MappingCollection<dim> &mapping,
                const DoFHandler<dim>            &dof_handler,
                const hp::QCollection<dim>       &quad,
                const AffineConstraints<number>  &constraints,
                VectorType                       &system_rhs);

    types::global_dof_index m() const;

    number el(unsigned int, unsigned int) const;

    void initialize_dof_vector(VectorType &vec) const;

    void vmult(VectorType &dst, const VectorType &src) const;

    void Tvmult(VectorType &dst, const VectorType &src) const;

    const TrilinosWrappers::SparseMatrix &get_system_matrix() const;

    void compute_inverse_diagonal(VectorType &diagonal) const;

  private:
    void do_cell_integral_local(FECellIntegrator &integrator) const;

    void do_cell_integral_global(FECellIntegrator &integrator,
                                 VectorType       &dst,
                                 const VectorType &src) const;


    void do_cell_integral_range(
      const MatrixFree<dim, number>               &matrix_free,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &range) const;

    MatrixFree<dim, number> matrix_free;

    // To solve the equation system on the coarsest level with an AMG
    // preconditioner, we need an actual system matrix on the coarsest level.
    // For this purpose, we provide a mechanism that optionally computes a
    // matrix from the matrix-free formulation, for which we introduce a
    // dedicated SparseMatrix object. In the default case, this matrix stays
    // empty. Once `get_system_matrix()` is called, this matrix is filled (lazy
    // allocation). Since this is a `const` function, we need the "mutable"
    // keyword here. We also need a the constraints object to build the matrix.
    AffineConstraints<number>              constraints;
    mutable TrilinosWrappers::SparseMatrix system_matrix;
  };



  // The following section contains functions to initialize and reinitialize
  // the class. In particular, these functions initialize the internal
  // MatrixFree instance. For sake of simplicity, we also compute the system
  // right-hand-side vector.
  template <int dim, typename number>
  LaplaceOperator<dim, number>::LaplaceOperator(
    const hp::MappingCollection<dim> &mapping,
    const DoFHandler<dim>            &dof_handler,
    const hp::QCollection<dim>       &quad,
    const AffineConstraints<number>  &constraints,
    VectorType                       &system_rhs)
  {
    this->reinit(mapping, dof_handler, quad, constraints, system_rhs);
  }



  template <int dim, typename number>
  void LaplaceOperator<dim, number>::reinit(
    const hp::MappingCollection<dim> &mapping,
    const DoFHandler<dim>            &dof_handler,
    const hp::QCollection<dim>       &quad,
    const AffineConstraints<number>  &constraints,
    VectorType                       &system_rhs)
  {
    // Clear internal data structures (in the case that the operator is reused).
    this->system_matrix.clear();

    // Copy the constraints, since they might be needed for computation of the
    // system matrix later on.
    this->constraints.copy_from(constraints);

    // Set up MatrixFree. At the quadrature points, we only need to evaluate
    // the gradient of the solution and test with the gradient of the shape
    // functions so that we only need to set the flag `update_gradients`.
    typename MatrixFree<dim, number>::AdditionalData data;
    data.mapping_update_flags = update_gradients;

    matrix_free.reinit(mapping, dof_handler, constraints, quad, data);

    // Compute the right-hand side vector. For this purpose, we set up a second
    // MatrixFree instance that uses a modified AffineConstraints not containing
    // the constraints due to Dirichlet-boundary conditions. This modified
    // operator is applied to a vector with only the Dirichlet values set. The
    // result is the negative right-hand-side vector.
    {
      AffineConstraints<number> constraints_without_dbc(
        dof_handler.locally_owned_dofs(),
        DoFTools::extract_locally_relevant_dofs(dof_handler));

      DoFTools::make_hanging_node_constraints(dof_handler,
                                              constraints_without_dbc);
      constraints_without_dbc.close();

      VectorType b, x;

      this->initialize_dof_vector(system_rhs);

      MatrixFree<dim, number> matrix_free;
      matrix_free.reinit(
        mapping, dof_handler, constraints_without_dbc, quad, data);

      matrix_free.initialize_dof_vector(b);
      matrix_free.initialize_dof_vector(x);

      constraints.distribute(x);

      matrix_free.cell_loop(&LaplaceOperator::do_cell_integral_range,
                            this,
                            b,
                            x);

      constraints.set_zero(b);

      system_rhs -= b;
    }
  }



  // The following functions are implicitly needed by the multigrid algorithm,
  // including the smoothers.

  // Since we do not have a matrix, query the DoFHandler for the number of
  // degrees of freedom.
  template <int dim, typename number>
  types::global_dof_index LaplaceOperator<dim, number>::m() const
  {
    return matrix_free.get_dof_handler().n_dofs();
  }



  // Access a particular element in the matrix. This function is neither
  // needed nor implemented, however, is required to compile the program.
  template <int dim, typename number>
  number LaplaceOperator<dim, number>::el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }



  // Initialize the given vector. We simply delegate the task to the
  // MatrixFree function with the same name.
  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }



  // Perform an operator evaluation by looping with the help of MatrixFree
  // over all cells and evaluating the effect of the cell integrals (see also:
  // `do_cell_integral_local()` and `do_cell_integral_global()`).
  template <int dim, typename number>
  void LaplaceOperator<dim, number>::vmult(VectorType       &dst,
                                           const VectorType &src) const
  {
    this->matrix_free.cell_loop(
      &LaplaceOperator::do_cell_integral_range, this, dst, src, true);
  }



  // Perform the transposed operator evaluation. Since we are considering
  // symmetric "matrices", this function can simply delegate it task to vmult().
  template <int dim, typename number>
  void LaplaceOperator<dim, number>::Tvmult(VectorType       &dst,
                                            const VectorType &src) const
  {
    this->vmult(dst, src);
  }



  // Since we do not have a system matrix, we cannot loop over the the
  // diagonal entries of the matrix. Instead, we compute the diagonal by
  // performing a sequence of operator evaluations to unit basis vectors.
  // For this purpose, an optimized function from the MatrixFreeTools
  // namespace is used. The inversion is performed manually afterwards.
  template <int dim, typename number>
  void LaplaceOperator<dim, number>::compute_inverse_diagonal(
    VectorType &diagonal) const
  {
    this->matrix_free.initialize_dof_vector(diagonal);
    MatrixFreeTools::compute_diagonal(matrix_free,
                                      diagonal,
                                      &LaplaceOperator::do_cell_integral_local,
                                      this);

    for (auto &i : diagonal)
      i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0;
  }



  // In the matrix-free context, no system matrix is set up during
  // initialization of this class. As a consequence, it has to be computed
  // here if it should be requested. Since the matrix is only computed in
  // this tutorial for linear elements (on the coarse grid), this is
  // acceptable.
  // The matrix entries are obtained via sequence of operator evaluations.
  // For this purpose, the optimized function MatrixFreeTools::compute_matrix()
  // is used. The matrix will only be computed if it has not been set up yet
  // (lazy allocation).
  template <int dim, typename number>
  const TrilinosWrappers::SparseMatrix &
  LaplaceOperator<dim, number>::get_system_matrix() const
  {
    if (system_matrix.m() == 0 && system_matrix.n() == 0)
      {
        const auto &dof_handler = this->matrix_free.get_dof_handler();

        TrilinosWrappers::SparsityPattern dsp(
          dof_handler.locally_owned_dofs(),
          dof_handler.get_triangulation().get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints);

        dsp.compress();
        system_matrix.reinit(dsp);

        MatrixFreeTools::compute_matrix(
          matrix_free,
          constraints,
          system_matrix,
          &LaplaceOperator::do_cell_integral_local,
          this);
      }

    return this->system_matrix;
  }



  // Perform cell integral on a cell batch without gathering and scattering
  // the values. This function is needed for the MatrixFreeTools functions
  // since these functions operate directly on the buffers of FEEvaluation.
  template <int dim, typename number>
  void LaplaceOperator<dim, number>::do_cell_integral_local(
    FECellIntegrator &integrator) const
  {
    integrator.evaluate(EvaluationFlags::gradients);

    for (const unsigned int q : integrator.quadrature_point_indices())
      integrator.submit_gradient(integrator.get_gradient(q), q);

    integrator.integrate(EvaluationFlags::gradients);
  }



  // Same as above but with access to the global vectors.
  template <int dim, typename number>
  void LaplaceOperator<dim, number>::do_cell_integral_global(
    FECellIntegrator &integrator,
    VectorType       &dst,
    const VectorType &src) const
  {
    integrator.gather_evaluate(src, EvaluationFlags::gradients);

    for (const unsigned int q : integrator.quadrature_point_indices())
      integrator.submit_gradient(integrator.get_gradient(q), q);

    integrator.integrate_scatter(EvaluationFlags::gradients, dst);
  }



  // This function loops over all cell batches within a cell-batch range and
  // calls the above function.
  template <int dim, typename number>
  void LaplaceOperator<dim, number>::do_cell_integral_range(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegrator integrator(matrix_free, range);

    for (unsigned cell = range.first; cell < range.second; ++cell)
      {
        integrator.reinit(cell);

        do_cell_integral_global(integrator, dst, src);
      }
  }



  // @sect3{Solver and preconditioner}

  // @sect4{Conjugate-gradient solver with multigrid preconditioner}

  // This function solves the equation system with a sequence of provided
  // multigrid objects. It is meant to be treated as general as possible, hence
  // the multitude of template parameters.
  template <typename VectorType,
            int dim,
            typename SystemMatrixType,
            typename LevelMatrixType,
            typename MGTransferType>
  static void
  mg_solve(SolverControl             &solver_control,
           VectorType                &dst,
           const VectorType          &src,
           const MultigridParameters &mg_data,
           const DoFHandler<dim>     &dof,
           const SystemMatrixType    &fine_matrix,
           const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices,
           const MGTransferType                                  &mg_transfer)
  {
    AssertThrow(mg_data.coarse_solver.type == "cg_with_amg",
                ExcNotImplemented());
    AssertThrow(mg_data.smoother.type == "chebyshev", ExcNotImplemented());

    const unsigned int min_level = mg_matrices.min_level();
    const unsigned int max_level = mg_matrices.max_level();

    using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
    using SmootherType               = PreconditionChebyshev<LevelMatrixType,
                                               VectorType,
                                               SmootherPreconditionerType>;
    using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;

    // We initialize level operators and Chebyshev smoothers here.
    mg::Matrix<VectorType> mg_matrix(mg_matrices);

    MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
      min_level, max_level);

    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        smoother_data[level].preconditioner =
          std::make_shared<SmootherPreconditionerType>();
        mg_matrices[level]->compute_inverse_diagonal(
          smoother_data[level].preconditioner->get_vector());
        smoother_data[level].smoothing_range = mg_data.smoother.smoothing_range;
        smoother_data[level].degree          = mg_data.smoother.degree;
        smoother_data[level].eig_cg_n_iterations =
          mg_data.smoother.eig_cg_n_iterations;
      }

    MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType>
      mg_smoother;
    mg_smoother.initialize(mg_matrices, smoother_data);

    // Next, we initialize the coarse-grid solver. We use conjugate-gradient
    // method with AMG as preconditioner.
    ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter,
                                                mg_data.coarse_solver.abstol,
                                                mg_data.coarse_solver.reltol,
                                                false,
                                                false);
    SolverCG<VectorType> coarse_grid_solver(coarse_grid_solver_control);

    std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

    TrilinosWrappers::PreconditionAMG                 precondition_amg;
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
    amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps;
    amg_data.n_cycles        = mg_data.coarse_solver.n_cycles;
    amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str();

    precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(),
                                amg_data);

    mg_coarse =
      std::make_unique<MGCoarseGridIterativeSolver<VectorType,
                                                   SolverCG<VectorType>,
                                                   LevelMatrixType,
                                                   decltype(precondition_amg)>>(
        coarse_grid_solver, *mg_matrices[min_level], precondition_amg);

    // Finally, we create the Multigrid object, convert it to a preconditioner,
    // and use it inside of a conjugate-gradient solver to solve the linear
    // system of equations.
    Multigrid<VectorType> mg(
      mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother);

    PreconditionerType preconditioner(dof, mg, mg_transfer);

    SolverCG<VectorType>(solver_control)
      .solve(fine_matrix, dst, src, preconditioner);
  }



  // @sect4{Hybrid polynomial/geometric-global-coarsening multigrid preconditioner}

  // The above function deals with the actual solution for a given sequence of
  // multigrid objects. This functions creates the actual multigrid levels, in
  // particular the operators, and the transfer operator as a
  // MGTransferGlobalCoarsening object.
  template <typename VectorType, typename OperatorType, int dim>
  void solve_with_gmg(SolverControl                    &solver_control,
                      const OperatorType               &system_matrix,
                      VectorType                       &dst,
                      const VectorType                 &src,
                      const MultigridParameters        &mg_data,
                      const hp::MappingCollection<dim> &mapping_collection,
                      const DoFHandler<dim>            &dof_handler,
                      const hp::QCollection<dim>       &quadrature_collection)
  {
    // Create a DoFHandler and operator for each multigrid level,
    // as well as, create transfer operators. To be able to
    // set up the operators, we need a set of DoFHandler that we create
    // via global coarsening of p or h. For latter, we need also a sequence
    // of Triangulation objects that are obtained by
    // Triangulation::coarsen_global().
    //
    // In case no h-transfer is requested, we provide an empty deleter for the
    // `emplace_back()` function, since the Triangulation of our DoFHandler is
    // an external field and its destructor is called somewhere else.
    MGLevelObject<DoFHandler<dim>>                     dof_handlers;
    MGLevelObject<std::unique_ptr<OperatorType>>       operators;
    MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;

    std::vector<std::shared_ptr<const Triangulation<dim>>>
      coarse_grid_triangulations;
    if (mg_data.transfer.perform_h_transfer)
      coarse_grid_triangulations =
        MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
          dof_handler.get_triangulation());
    else
      coarse_grid_triangulations.emplace_back(
        &(dof_handler.get_triangulation()), [](auto *) {});

    // Determine the total number of levels for the multigrid operation and
    // allocate sufficient memory for all levels.
    const unsigned int n_h_levels = coarse_grid_triangulations.size() - 1;

    const auto get_max_active_fe_degree = [&](const auto &dof_handler) {
      unsigned int max = 0;

      for (auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          max =
            std::max(max, dof_handler.get_fe(cell->active_fe_index()).degree);

      return Utilities::MPI::max(max, MPI_COMM_WORLD);
    };

    const unsigned int n_p_levels =
      MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
        get_max_active_fe_degree(dof_handler), mg_data.transfer.p_sequence)
        .size();

    std::map<unsigned int, unsigned int> fe_index_for_degree;
    for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i)
      {
        const unsigned int degree = dof_handler.get_fe(i).degree;
        Assert(fe_index_for_degree.find(degree) == fe_index_for_degree.end(),
               ExcMessage("FECollection does not contain unique degrees."));
        fe_index_for_degree[degree] = i;
      }

    unsigned int minlevel = 0;
    unsigned int maxlevel = n_h_levels + n_p_levels - 1;

    dof_handlers.resize(minlevel, maxlevel);
    operators.resize(minlevel, maxlevel);
    transfers.resize(minlevel, maxlevel);

    // Loop from the minimum (coarsest) to the maximum (finest) level and set up
    // DoFHandler accordingly. We start with the h-levels, where we distribute
    // on increasingly finer meshes linear elements.
    for (unsigned int l = 0; l < n_h_levels; ++l)
      {
        dof_handlers[l].reinit(*coarse_grid_triangulations[l]);
        dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection());
      }

    // After we reached the finest mesh, we will adjust the polynomial degrees
    // on each level. We reverse iterate over our data structure and start at
    // the finest mesh that contains all information about the active FE
    // indices. We then lower the polynomial degree of each cell level by level.
    for (unsigned int i = 0, l = maxlevel; i < n_p_levels; ++i, --l)
      {
        dof_handlers[l].reinit(dof_handler.get_triangulation());

        if (l == maxlevel) // finest level
          {
            auto &dof_handler_mg = dof_handlers[l];

            auto cell_other = dof_handler.begin_active();
            for (auto &cell : dof_handler_mg.active_cell_iterators())
              {
                if (cell->is_locally_owned())
                  cell->set_active_fe_index(cell_other->active_fe_index());
                ++cell_other;
              }
          }
        else // coarse level
          {
            auto &dof_handler_fine   = dof_handlers[l + 1];
            auto &dof_handler_coarse = dof_handlers[l + 0];

            auto cell_other = dof_handler_fine.begin_active();
            for (auto &cell : dof_handler_coarse.active_cell_iterators())
              {
                if (cell->is_locally_owned())
                  {
                    const unsigned int next_degree =
                      MGTransferGlobalCoarseningTools::
                        create_next_polynomial_coarsening_degree(
                          cell_other->get_fe().degree,
                          mg_data.transfer.p_sequence);
                    Assert(fe_index_for_degree.find(next_degree) !=
                             fe_index_for_degree.end(),
                           ExcMessage("Next polynomial degree in sequence "
                                      "does not exist in FECollection."));

                    cell->set_active_fe_index(fe_index_for_degree[next_degree]);
                  }
                ++cell_other;
              }
          }

        dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection());
      }

    // Next, we will create all data structures additionally needed on each
    // multigrid level. This involves determining constraints with homogeneous
    // Dirichlet boundary conditions, and building the operator just like on the
    // active level.
    MGLevelObject<AffineConstraints<typename VectorType::value_type>>
      constraints(minlevel, maxlevel);

    for (unsigned int level = minlevel; level <= maxlevel; ++level)
      {
        const auto &dof_handler = dof_handlers[level];
        auto       &constraint  = constraints[level];

        constraint.reinit(dof_handler.locally_owned_dofs(),
                          DoFTools::extract_locally_relevant_dofs(dof_handler));

        DoFTools::make_hanging_node_constraints(dof_handler, constraint);
        VectorTools::interpolate_boundary_values(mapping_collection,
                                                 dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraint);
        constraint.close();

        VectorType dummy;

        operators[level] = std::make_unique<OperatorType>(mapping_collection,
                                                          dof_handler,
                                                          quadrature_collection,
                                                          constraint,
                                                          dummy);
      }

    // Set up intergrid operators and collect transfer operators within a single
    // operator as needed by the Multigrid solver class.
    for (unsigned int level = minlevel; level < maxlevel; ++level)
      transfers[level + 1].reinit(dof_handlers[level + 1],
                                  dof_handlers[level],
                                  constraints[level + 1],
                                  constraints[level]);

    MGTransferGlobalCoarsening<dim, VectorType> transfer(
      transfers, [&](const auto l, auto &vec) {
        operators[l]->initialize_dof_vector(vec);
      });

    // Finally, proceed to solve the problem with multigrid.
    mg_solve(solver_control,
             dst,
             src,
             mg_data,
             dof_handler,
             system_matrix,
             operators,
             transfer);
  }



  // @sect3{The <code>LaplaceProblem</code> class template}

  // Now we will finally declare the main class of this program, which solves
  // the Laplace equation on subsequently refined function spaces. Its structure
  // will look familiar as it is similar to the main classes of step-27 and
  // step-40. There are basically just two additions:
  // - The SparseMatrix object that would hold the system matrix has been
  //   replaced by an object of the LaplaceOperator class for the MatrixFree
  //   formulation.
  // - An object of parallel::CellWeights, which will help us with load
  //   balancing, has been added.
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem(const Parameters &parameters);

    void run();

  private:
    void initialize_grid();
    void setup_system();
    void print_diagnostics();
    void solve_system();
    void compute_indicators();
    void adapt_resolution();
    void output_results(const unsigned int cycle);

    MPI_Comm mpi_communicator;

    const Parameters prm;

    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;

    hp::MappingCollection<dim> mapping_collection;
    hp::FECollection<dim>      fe_collection;
    hp::QCollection<dim>       quadrature_collection;
    hp::QCollection<dim - 1>   face_quadrature_collection;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    LaplaceOperator<dim, double>               laplace_operator;
    LinearAlgebra::distributed::Vector<double> locally_relevant_solution;
    LinearAlgebra::distributed::Vector<double> system_rhs;

    std::unique_ptr<FESeries::Legendre<dim>>    legendre;
    std::unique_ptr<parallel::CellWeights<dim>> cell_weights;

    Vector<float> estimated_error_per_cell;
    Vector<float> hp_decision_indicators;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
  };



  // @sect3{The <code>LaplaceProblem</code> class implementation}

  // @sect4{Constructor}

  // The constructor starts with an initializer list that looks similar to the
  // one of step-40. We again prepare the ConditionalOStream object to allow
  // only the first process to output anything over the console, and initialize
  // the computing timer properly.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem(const Parameters &parameters)
    : mpi_communicator(MPI_COMM_WORLD)
    , prm(parameters)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
  {
    Assert(prm.min_h_level <= prm.max_h_level,
           ExcMessage(
             "Triangulation level limits have been incorrectly set up."));
    Assert(prm.min_p_degree <= prm.max_p_degree,
           ExcMessage("FECollection degrees have been incorrectly set up."));

    // We need to prepare the data structures for the hp-functionality in the
    // actual body of the constructor, and create corresponding objects for
    // every degree in the specified range from the parameter struct. As we are
    // only dealing with non-distorted rectangular cells, a linear mapping
    // object is sufficient in this context.
    //
    // In the Parameters struct, we provide ranges for levels on which the
    // function space is operating with a reasonable resolution. The multigrid
    // algorithm requires linear elements on the coarsest possible level. So we
    // start with the lowest polynomial degree and fill the collection with
    // consecutively higher degrees until the user-specified maximum is
    // reached.
    mapping_collection.push_back(MappingQ1<dim>());

    for (unsigned int degree = 1; degree <= prm.max_p_degree; ++degree)
      {
        fe_collection.push_back(FE_Q<dim>(degree));
        quadrature_collection.push_back(QGauss<dim>(degree + 1));
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
      }

    // As our FECollection contains more finite elements than we want to use for
    // the finite element approximation of our solution, we would like to limit
    // the range on which active FE indices can operate on. For this, the
    // FECollection class allows to register a hierarchy that determines the
    // succeeding and preceding finite element in case of of p-refinement and
    // p-coarsening, respectively. All functions in the hp::Refinement namespace
    // consult this hierarchy to determine future FE indices. We will register
    // such a hierarchy that only works on finite elements with polynomial
    // degrees in the proposed range <code>[min_p_degree, max_p_degree]</code>.
    const unsigned int min_fe_index = prm.min_p_degree - 1;
    fe_collection.set_hierarchy(
      /*next_index=*/
      [](const typename hp::FECollection<dim> &fe_collection,
         const unsigned int                    fe_index) -> unsigned int {
        return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 :
                                                         fe_index;
      },
      /*previous_index=*/
      [min_fe_index](const typename hp::FECollection<dim> &,
                     const unsigned int fe_index) -> unsigned int {
        Assert(fe_index >= min_fe_index,
               ExcMessage("Finite element is not part of hierarchy!"));
        return (fe_index > min_fe_index) ? fe_index - 1 : fe_index;
      });

    // We initialize the FESeries::Legendre object in the default configuration
    // for smoothness estimation.
    legendre = std::make_unique<FESeries::Legendre<dim>>(
      SmoothnessEstimator::Legendre::default_fe_series(fe_collection));

    // The next part is going to be tricky. During execution of refinement, a
    // few hp-algorithms need to interfere with the actual refinement process on
    // the Triangulation object. We do this by connecting several functions to
    // Triangulation::Signals: signals will be called at different stages during
    // the actual refinement process and trigger all connected functions. We
    // require this functionality for load balancing and to limit the polynomial
    // degrees of neighboring cells.
    //
    // For the former, we would like to assign a weight to every cell that is
    // proportional to the number of degrees of freedom of its future finite
    // element. The library offers a class parallel::CellWeights that allows to
    // easily attach individual weights at the right place during the refinement
    // process, i.e., after all refine and coarsen flags have been set correctly
    // for hp-adaptation and right before repartitioning for load balancing is
    // about to happen. Functions can be registered that will attach weights in
    // the form that $a (n_\text{dofs})^b$ with a provided pair of parameters
    // $(a,b)$. We register such a function in the following.
    //
    // For load balancing, efficient solvers like the one we use should scale
    // linearly with the number of degrees of freedom owned. We set the
    // parameters for cell weighting correspondingly: A weighting factor of $1$
    // and an exponent of $1$ (see the definitions of the `weighting_factor` and
    // `weighting_exponent` above).
    cell_weights = std::make_unique<parallel::CellWeights<dim>>(
      dof_handler,
      parallel::CellWeights<dim>::ndofs_weighting(
        {prm.weighting_factor, prm.weighting_exponent}));

    // In h-adaptive applications, we ensure a 2:1 mesh balance by limiting the
    // difference of refinement levels of neighboring cells to one. With the
    // second call in the following code snippet, we will ensure the same for
    // p-levels on neighboring cells: levels of future finite elements are not
    // allowed to differ by more than a specified difference. The function
    // hp::Refinement::limit_p_level_difference takes care of this, but needs to
    // be connected to a very specific signal in the parallel context. The issue
    // is that we need to know how the mesh will be actually refined to set
    // future FE indices accordingly. As we ask the p4est oracle to perform
    // refinement, we need to ensure that the Triangulation has been updated
    // with the adaptation flags of the oracle first. An instantiation of
    // parallel::distributed::TemporarilyMatchRefineFlags does exactly
    // that for the duration of its life. Thus, we will create an object of this
    // class right before limiting the p-level difference, and connect the
    // corresponding lambda function to the signal
    // Triangulation::Signals::post_p4est_refinement, which will be triggered
    // after the oracle got refined, but before the Triangulation is refined.
    // Furthermore, we specify that this function will be connected to the front
    // of the signal, to ensure that the modification is performed before any
    // other function connected to the same signal.
    triangulation.signals.post_p4est_refinement.connect(
      [&, min_fe_index]() {
        const parallel::distributed::TemporarilyMatchRefineFlags<dim>
          refine_modifier(triangulation);
        hp::Refinement::limit_p_level_difference(dof_handler,
                                                 prm.max_p_level_difference,
                                                 /*contains=*/min_fe_index);
      },
      boost::signals2::at_front);
  }



  // @sect4{LaplaceProblem::initialize_grid}

  // For a L-shaped domain, we could use the function GridGenerator::hyper_L()
  // as demonstrated in step-50. However in the 2d case, that particular
  // function removes the first quadrant, while we need the fourth quadrant
  // removed in our scenario. Thus, we will use a different function
  // GridGenerator::subdivided_hyper_L() which gives us more options to create
  // the mesh. Furthermore, we formulate that function in a way that it also
  // generates a 3d mesh: the 2d L-shaped domain will basically elongated by 1
  // in the positive z-direction.
  //
  // We first pretend to build a GridGenerator::subdivided_hyper_rectangle().
  // The parameters that we need to provide are Point objects for the lower left
  // and top right corners, as well as the number of repetitions that the base
  // mesh will have in each direction. We provide them for the first two
  // dimensions and treat the higher third dimension separately.
  //
  // To create a L-shaped domain, we need to remove the excess cells. For this,
  // we specify the <code>cells_to_remove</code> accordingly. We would like to
  // remove one cell in every cell from the negative direction, but remove one
  // from the positive x-direction.
  //
  // On the coarse grid, we set the initial active FE indices and distribute the
  // degrees of freedom once. We do that in order to assign the hp::FECollection
  // to the DoFHandler, so that all cells know how many DoFs they are going to
  // have. This step is mandatory for the weighted load balancing algorithm,
  // which will be called implicitly in
  // parallel::distributed::Triangulation::refine_global().
  template <int dim>
  void LaplaceProblem<dim>::initialize_grid()
  {
    TimerOutput::Scope t(computing_timer, "initialize grid");

    std::vector<unsigned int> repetitions(dim);
    Point<dim>                bottom_left, top_right;
    for (unsigned int d = 0; d < dim; ++d)
      if (d < 2)
        {
          repetitions[d] = 2;
          bottom_left[d] = -1.;
          top_right[d]   = 1.;
        }
      else
        {
          repetitions[d] = 1;
          bottom_left[d] = 0.;
          top_right[d]   = 1.;
        }

    std::vector<int> cells_to_remove(dim, 1);
    cells_to_remove[0] = -1;

    GridGenerator::subdivided_hyper_L(
      triangulation, repetitions, bottom_left, top_right, cells_to_remove);

    const unsigned int min_fe_index = prm.min_p_degree - 1;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index(min_fe_index);

    dof_handler.distribute_dofs(fe_collection);

    triangulation.refine_global(prm.min_h_level);
  }



  // @sect4{LaplaceProblem::setup_system}

  // This function looks exactly the same to the one of step-40, but you will
  // notice the absence of the system matrix as well as the scaffold that
  // surrounds it. Instead, we will initialize the MatrixFree formulation of the
  // <code>laplace_operator</code> here. For boundary conditions, we will use
  // the Solution class introduced earlier in this tutorial.
  template <int dim>
  void LaplaceProblem<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup system");

    dof_handler.distribute_dofs(fe_collection);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping_collection, dof_handler, 0, Solution<dim>(), constraints);
    constraints.close();

    laplace_operator.reinit(mapping_collection,
                            dof_handler,
                            quadrature_collection,
                            constraints,
                            system_rhs);
  }



  // @sect4{LaplaceProblem::print_diagnostics}

  // This is a function that prints additional diagnostics about the equation
  // system and its partitioning. In addition to the usual global number of
  // active cells and degrees of freedom, we also output their local
  // equivalents. For a regulated output, we will communicate the local
  // quantities with a Utilities::MPI::gather operation to the first process
  // which will then output all information. Output of local quantities is
  // limited to the first 8 processes to avoid cluttering the terminal.
  //
  // On all other processes, the containers for the collected data remain empty.
  // To ensure that we do not access invalid memory with the insertion operator
  // (`<<`) on these processes, we need to check that the containers are not
  // empty.
  //
  // Furthermore, we would like to print the frequencies of the polynomial
  // degrees in the numerical discretization. Since this information is only
  // stored locally, we will count the finite elements on locally owned cells
  // and later communicate them via Utilities::MPI::sum.
  template <int dim>
  void LaplaceProblem<dim>::print_diagnostics()
  {
    const unsigned int first_n_processes =
      std::min<unsigned int>(8,
                             Utilities::MPI::n_mpi_processes(mpi_communicator));
    const bool output_cropped =
      first_n_processes < Utilities::MPI::n_mpi_processes(mpi_communicator);

    {
      pcout << "   Number of active cells:       "
            << triangulation.n_global_active_cells() << std::endl
            << "     by partition:              ";

      std::vector<unsigned int> n_active_cells_per_subdomain =
        Utilities::MPI::gather(mpi_communicator,
                               triangulation.n_locally_owned_active_cells());
      for (unsigned int i = 0; i < first_n_processes; ++i)
        if (n_active_cells_per_subdomain.size() > 0)
          pcout << ' ' << n_active_cells_per_subdomain[i];
      if (output_cropped)
        pcout << " ...";
      pcout << std::endl;
    }

    {
      pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl
            << "     by partition:              ";

      std::vector<types::global_dof_index> n_dofs_per_subdomain =
        Utilities::MPI::gather(mpi_communicator,
                               dof_handler.n_locally_owned_dofs());
      for (unsigned int i = 0; i < first_n_processes; ++i)
        if (n_dofs_per_subdomain.size() > 0)
          pcout << ' ' << n_dofs_per_subdomain[i];
      if (output_cropped)
        pcout << " ...";
      pcout << std::endl;
    }

    {
      std::vector<types::global_dof_index> n_constraints_per_subdomain =
        Utilities::MPI::gather(mpi_communicator, constraints.n_constraints());

      pcout << "   Number of constraints:        "
            << std::accumulate(n_constraints_per_subdomain.begin(),
                               n_constraints_per_subdomain.end(),
                               0)
            << std::endl
            << "     by partition:              ";
      for (unsigned int i = 0; i < first_n_processes; ++i)
        if (n_constraints_per_subdomain.size() > 0)
          pcout << ' ' << n_constraints_per_subdomain[i];
      if (output_cropped)
        pcout << " ...";
      pcout << std::endl;
    }

    {
      std::vector<unsigned int> n_fe_indices(fe_collection.size(), 0);
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          n_fe_indices[cell->active_fe_index()]++;

      Utilities::MPI::sum(n_fe_indices, mpi_communicator, n_fe_indices);

      pcout << "   Frequencies of poly. degrees:";
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        if (n_fe_indices[i] > 0)
          pcout << ' ' << fe_collection[i].degree << ':' << n_fe_indices[i];
      pcout << std::endl;
    }
  }



  // @sect4{LaplaceProblem::solve_system}

  // The scaffold around the solution is similar to the one of step-40. We
  // prepare a vector that matches the requirements of MatrixFree and collect
  // the locally-relevant degrees of freedoms we solved the equation system. The
  // solution happens with the function introduced earlier.
  template <int dim>
  void LaplaceProblem<dim>::solve_system()
  {
    TimerOutput::Scope t(computing_timer, "solve system");

    LinearAlgebra::distributed::Vector<double> completely_distributed_solution;
    laplace_operator.initialize_dof_vector(completely_distributed_solution);

    SolverControl solver_control(system_rhs.size(),
                                 prm.tolerance_factor * system_rhs.l2_norm());

    solve_with_gmg(solver_control,
                   laplace_operator,
                   completely_distributed_solution,
                   system_rhs,
                   prm.mg_data,
                   mapping_collection,
                   dof_handler,
                   quadrature_collection);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution.copy_locally_owned_data_from(
      completely_distributed_solution);
    locally_relevant_solution.update_ghost_values();
  }



  // @sect4{LaplaceProblem::compute_indicators}

  // This function contains only a part of the typical <code>refine_grid</code>
  // function from other tutorials and is new in that sense. Here, we will only
  // calculate all indicators for adaptation with actually refining the grid. We
  // do this for the purpose of writing all indicators to the file system, so we
  // store them for later.
  //
  // Since we are dealing the an elliptic problem, we will make use of the
  // KellyErrorEstimator again, but with a slight difference. Modifying the
  // scaling factor of the underlying face integrals to be dependent on the
  // actual polynomial degree of the neighboring elements is favorable in
  // hp-adaptive applications @cite davydov2017hp. We can do this by specifying
  // the very last parameter from the additional ones you notices. The others
  // are actually just the defaults.
  //
  // For the purpose of hp-adaptation, we will calculate smoothness estimates
  // with the strategy presented in the tutorial introduction and use the
  // implementation in SmoothnessEstimator::Legendre. In the Parameters struct,
  // we set the minimal polynomial degree to 2 as it seems that the smoothness
  // estimation algorithms have trouble with linear elements.
  template <int dim>
  void LaplaceProblem<dim>::compute_indicators()
  {
    TimerOutput::Scope t(computing_timer, "compute indicators");

    estimated_error_per_cell.grow_or_shrink(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      locally_relevant_solution,
      estimated_error_per_cell,
      /*component_mask=*/ComponentMask(),
      /*coefficients=*/nullptr,
      /*n_threads=*/numbers::invalid_unsigned_int,
      /*subdomain_id=*/numbers::invalid_subdomain_id,
      /*material_id=*/numbers::invalid_material_id,
      /*strategy=*/
      KellyErrorEstimator<dim>::Strategy::face_diameter_over_twice_max_degree);

    hp_decision_indicators.grow_or_shrink(triangulation.n_active_cells());
    SmoothnessEstimator::Legendre::coefficient_decay(*legendre,
                                                     dof_handler,
                                                     locally_relevant_solution,
                                                     hp_decision_indicators);
  }



  // @sect4{LaplaceProblem::adapt_resolution}

  // With the previously calculated indicators, we will finally flag all cells
  // for adaptation and also execute refinement in this function. As in previous
  // tutorials, we will use the "fixed number" strategy, but now for
  // hp-adaptation.
  template <int dim>
  void LaplaceProblem<dim>::adapt_resolution()
  {
    TimerOutput::Scope t(computing_timer, "adapt resolution");

    // First, we will set refine and coarsen flags based on the error estimates
    // on each cell. There is nothing new here.
    //
    // We will use general refine and coarsen fractions that have been
    // elaborated in the other deal.II tutorials: using the fixed number
    // strategy, we will flag 30% of all cells for refinement and 3% for
    // coarsening, as provided in the Parameters struct.
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation,
      estimated_error_per_cell,
      prm.refine_fraction,
      prm.coarsen_fraction);

    // Next, we will make all adjustments for hp-adaptation. We want to refine
    // and coarsen those cells flagged in the previous step, but need to decide
    // if we would like to do it by adjusting the grid resolution or the
    // polynomial degree.
    //
    // The next function call sets future FE indices according to the previously
    // calculated smoothness indicators as p-adaptation indicators. These
    // indices will only be set on those cells that have refine or coarsen flags
    // assigned.
    //
    // For the p-adaptation fractions, we will take an educated guess. Since we
    // only expect a single singularity in our scenario, i.e., in the origin of
    // the domain, and a smooth solution anywhere else, we would like to
    // strongly prefer to use p-adaptation over h-adaptation. This reflects in
    // our choice of a fraction of 90% for both p-refinement and p-coarsening.
    hp::Refinement::p_adaptivity_fixed_number(dof_handler,
                                              hp_decision_indicators,
                                              prm.p_refine_fraction,
                                              prm.p_coarsen_fraction);

    // After setting all indicators, we will remove those that exceed the
    // specified limits of the provided level ranges in the Parameters struct.
    // This limitation naturally arises for p-adaptation as the number of
    // supplied finite elements is limited. In addition, we registered a custom
    // hierarchy for p-adaptation in the constructor. Now, we need to do this
    // manually in the h-adaptive context like in step-31.
    //
    // We will iterate over all cells on the designated min and max levels and
    // remove the corresponding flags. As an alternative, we could also flag
    // these cells for p-adaptation by setting future FE indices accordingly
    // instead of simply clearing the refine and coarsen flags.
    Assert(triangulation.n_levels() >= prm.min_h_level + 1 &&
             triangulation.n_levels() <= prm.max_h_level + 1,
           ExcInternalError());

    if (triangulation.n_levels() > prm.max_h_level)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(prm.max_h_level))
        cell->clear_refine_flag();

    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(prm.min_h_level))
      cell->clear_coarsen_flag();

    // At this stage, we have both the future FE indices and the classic refine
    // and coarsen flags set. The latter will be interpreted by
    // Triangulation::execute_coarsening_and_refinement() for h-adaptation, and
    // our previous modification ensures that the resulting Triangulation stays
    // within the specified level range.
    //
    // Now, we would like to only impose one type of adaptation on cells, which
    // is what the next function will sort out for us. In short, on cells which
    // have both types of indicators assigned, we will favor the p-adaptation
    // one and remove the h-adaptation one.
    hp::Refinement::choose_p_over_h(dof_handler);

    // In the end, we are left to execute coarsening and refinement. Here, not
    // only the grid will be updated, but also all previous future FE indices
    // will become active.
    //
    // Remember that we have attached functions to triangulation signals in the
    // constructor, will be triggered in this function call. So there is even
    // more happening: weighted repartitioning will be performed to ensure load
    // balancing, as well as we will limit the difference of p-levels between
    // neighboring cells.
    triangulation.execute_coarsening_and_refinement();
  }



  // @sect4{LaplaceProblem::output_results}

  // Writing results to the file system in parallel applications works exactly
  // like in step-40. In addition to the data containers that we prepared
  // throughout the tutorial, we would also like to write out the polynomial
  // degree of each finite element on the grid as well as the subdomain each
  // cell belongs to. We prepare necessary containers for this in the scope of
  // this function.
  template <int dim>
  void LaplaceProblem<dim>::output_results(const unsigned int cycle)
  {
    TimerOutput::Scope t(computing_timer, "output results");

    Vector<float> fe_degrees(triangulation.n_active_cells());
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        fe_degrees(cell->active_cell_index()) = cell->get_fe().degree;

    Vector<float> subdomain(triangulation.n_active_cells());
    for (auto &subd : subdomain)
      subd = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, "solution");
    data_out.add_data_vector(fe_degrees, "fe_degree");
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.add_data_vector(estimated_error_per_cell, "error");
    data_out.add_data_vector(hp_decision_indicators, "hp_indicator");
    data_out.build_patches(mapping_collection);

    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, mpi_communicator, 2, 1);
  }



  // @sect4{LaplaceProblem::run}

  // The actual run function again looks very familiar to step-40. The only
  // addition is the bracketed section that precedes the actual cycle loop.
  // Here, we will pre-calculate the Legendre transformation matrices. In
  // general, these will be calculated on the fly via lazy allocation whenever a
  // certain matrix is needed. For timing purposes however, we would like to
  // calculate them all at once before the actual time measurement begins. We
  // will thus designate their calculation to their own scope.
  template <int dim>
  void LaplaceProblem<dim>::run()
  {
    pcout << "Running with Trilinos on "
          << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    {
      pcout << "Calculating transformation matrices..." << std::endl;
      TimerOutput::Scope t(computing_timer, "calculate transformation");
      legendre->precalculate_all_transformation_matrices();
    }

    for (unsigned int cycle = 0; cycle < prm.n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          initialize_grid();
        else
          adapt_resolution();

        setup_system();

        print_diagnostics();

        solve_system();

        compute_indicators();

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          output_results(cycle);

        computing_timer.print_summary();
        computing_timer.reset();

        pcout << std::endl;
      }
  }
} // namespace Step75



// @sect4{main()}

// The final function is the <code>main</code> function that will ultimately
// create and run a LaplaceOperator instantiation. Its structure is similar to
// most other tutorial programs.
int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step75;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      Parameters        prm;
      LaplaceProblem<2> laplace_problem(prm);
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
