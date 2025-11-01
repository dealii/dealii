/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2022 - 2025 by the deal.II authors
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
 * Author: Michał Wichrowski, Heidelberg University, 2025
 */

// @sect3{Include files}

// First include the necessary files from the deal.II library.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/tensor_product_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_patch_evaluation.h>



#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// This includes the utilities for the efficient implementation of
// matrix-free patch-smoothing methods.
#include <deal.II/matrix_free/patch_storage.h>
#include <deal.II/matrix_free/patch_distributors.h>
#include <deal.II/matrix_free/patch_smoother_base.h>

// Assembly of 1D matrices needed for local inverse
#include <deal.II/numerics/tensor_product_matrix_creator.h>

#include <fstream>
#include <iostream>



// @sect3{Operators}

namespace Operators
{
  using namespace dealii;

  /**
   *
   * This class implements the Laplace operator $-\Delta u$ as in step-37 using
   * matrix-free techniques based on deal.II's MatrixFree framework. It
   * derives from MatrixFreeOperators::Base, providing the necessary interface
   * for matrix-vector products (`vmult`, `Tvmult`) and diagonal computation.
   *
   *
   * @tparam dim The spatial dimension.
   * @tparam fe_degree The polynomial degree of the finite element.
   * @tparam number The scalar type for computations (e.g., double, float).
   */
  template <int dim, int fe_degree, typename number>
  class LaplaceOperator
    : public MatrixFreeOperators::
        Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    using value_type = number;
    using VectorType = LinearAlgebra::distributed::Vector<number>;


    /**
     * @brief Constructor.
     */
    LaplaceOperator() = default;

    /**
     * Clears internal data structures, including the diagonal entries and the
     * base class state.
     */
    void clear() override;

    /**
     * @brief Compute the diagonal entries of the operator.
     *
     * This function is not need for patch smoothing.
     */
    virtual void compute_diagonal() override
    {}

  private:
    /**
     * @brief Perform the matrix-vector product `dst += A*src`.
     *
     * This is the core function called by the public `vmult` and `Tvmult`
     * methods (and their variants) provided by the base class. It run
     * the cell loop using MatrixFree::cell_loop.
     *
     * @param dst The destination vector, to which the result is added.
     * @param src The source vector.
     */
    virtual void apply_add(VectorType       &dst,
                           const VectorType &src) const override;

    /**
     * @brief Local worker routine for the matrix-vector product.
     *
     */
    void
    local_apply(const MatrixFree<dim, number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;
  };


  /**
   *
   * Calls the base class clear function. We intentionally keep this
   * trivial override as an obvious placeholder so users extending the
   * operator (for example by adding more advanced local solvers,
   * internal buffers, or other per-operator state) remember to call
   * the base-class clear(). You can remove this override if you prefer,
   * but keeping it here acts as a convenient reminder when building on
   * this code.
   */
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::clear()
  {
    MatrixFreeOperators::Base<dim, VectorType>::clear();
  }

  /**
   * @brief Local apply implementation.
   *
   * Uses FEEvaluation to perform the cell-local computations for the
   * Laplace operator:
   * 1. Reinitialize FEEvaluation for the current cell batch.
   * 2. Read DoF values from the source vector `src`.
   * 3. Evaluate gradients at quadrature points.
   * 4. Submit the gradients (multiplied by coefficient, if any -- here
   * just 1.0) back for integration.
   * 5. Integrate the gradients against test function gradients.
   * 6. Distribute the local results to the global destination vector `dst`.
   */
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::local_apply(
    const MatrixFree<dim, number>               &data,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  /**
   * @brief Apply_add implementation.
   *
   * Uses the MatrixFree::cell_loop mechanism to execute the `local_apply`
   * function in parallel over cell batches, handling MPI communication
   * (ghost value updates and result compression) implicitly.
   */
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::apply_add(
    VectorType       &dst,
    const VectorType &src) const
  {
    this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
  }

  /**
   * @brief Implements a patch smoother for the Laplace operator.
   *
   * This class defines a smoother operation suitable for multigrid methods,
   * specifically designed for matrix-free operators. It operates on patches
   * of cells, applying a local inverse (provided by PatchTensorInverse) to
   * approximate the action of the global operator's inverse on the residual.
   *
   * It derives from  SmootherBase, which provides the framework
   * for patch-based smoothers, including handling patch iteration and data
   * distribution.
   *
   * @tparam dim The spatial dimension.
   * @tparam fe_degree The polynomial degree of the finite element.
   * @tparam number The scalar type for computations (e.g., double).
   */
  template <int dim, int fe_degree, typename number>
  class LaplacePatchSmoother : public PatchSmootherBase<dim, number>
  {
  public:
    using BaseType       = PatchSmootherBase<dim, number>;
    using VectorType     = LinearAlgebra::distributed::Vector<number>;
    using MatrixFreeType = MatrixFree<dim, number>;

    /**
     * @brief Type alias for the storage of managing patch data (DoF indices, etc.).
     */
    using PatchStorageType = PatchStorage<MatrixFreeType>;

    using VectorizedNumber = VectorizedArray<number>;

    /**
     * @brief The size (number of rows/columns) of the 1D matrices that form
     * the tensor product patch matrix. For a patch of two adjacent reference
     * cells in 1D and FE_Q elements of degree `fe_degree`, we consider the
     * unique 1D basis functions after imposing Dirichlet boundary conditions
     * on the ends of the patch. Imposing Dirichlet conditions on both patch
     * ends removes the two boundary degrees of freedom, yielding
     * 2*fe_degree - 1 unknown interior basis functions per axis.
     */
    const constexpr static unsigned int n_rows_1d = 2 * fe_degree - 1;

    const constexpr static unsigned int n_dofs_interior =
      Utilities::pow<unsigned int>(n_rows_1d, dim);

    /**
     * @brief Type aliases for the class providing the local patch inverse
     * and standard FEEvaluation on a single cell batch
     */
    using InverseType = TensorProductMatrixSymmetricSum<dim, number, n_rows_1d>;
    using FEEval      = FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number>;


    /**
     * @brief Type alias for additional data required by the smoother (the patch storage).
     */
    using AdditionalData = std::shared_ptr<PatchStorageType>;

    /**
     * @brief Local worker routine for applying the smoother on a range of patches.
     *
     * This function is called by the base class framework (e.g., within a
     * patch loop). For each patch in the given range, it performs the
     * following steps:
     * 1. Gather the relevant portion of the global right-hand side vector

     * `src` to a local patch vector `local_rhs`.
     * 2. Apply the global operator (implicitly via FEEvaluation) to the
     *    current global solution estimate `dst` and gathers the result to a
     *    local patch vector `local_residual`.
     * 3. Compute the actual residual on the patch: `local_residual = local_rhs
     * - A_patch * current_solution_on_patch`.
     * 4. Apply the inverse patch operator (`inverse.vmult`) to the residual
     *    to get a local correction: `local_correction = A_patch^{-1} *
     * local_residual`.
     * 5. Update the global solution `dst` by distributing the negative of the
     *    local correction.
     *
     * @param patch_storage The storage object containing patch information.
     * @param dst The global solution vector (input/output).
     * @param src The global right-hand side vector (input).
     * @param patch_range The range of patch indices to process.
     * @param do_forward Flag indicating the direction of iteration (forward or backward sweep).
     */
    void local_apply(const PatchStorageType                      &patch_storage,
                     VectorType                                  &dst,
                     const VectorType                            &src,
                     const typename PatchStorageType::PatchRange &patch_range,
                     const bool &do_forward) const override;

    /**
     *
     * Initialize the base class and the local patch inverse operator.
     *
     * @param patch_storage Shared pointer to the patch storage object.
     */
    void initialize(std::shared_ptr<PatchStorageType> &patch_storage) override;

    /**
     * @brief Clear internal data.
     *
     * Clears the base class state and the local patch inverse.
     */
    void clear() override;

  private:
    /**
     * @brief The tensor product matrix representing the patch operator.
     * Provides the `apply_inverse` method.
     */
    InverseType patch_inverse_operator;
  };

  /**
   * @brief Initialize implementation.
   *
   * Calls the base class initialize and initializes the `inverse` member if
   * the patch storage contains any patches.
   */
  template <int dim, int fe_degree, typename number>
  void LaplacePatchSmoother<dim, fe_degree, number>::initialize(
    std::shared_ptr<PatchStorageType> &patch_storage)
  {
    // After calling the initialize function of the base class to set up
    // the connection to the patch storage object, we check whether there
    // are any patches to be worked on the current processor. Then, we
    // construct an FEEvaluation object to provide access to the finite element
    // and the geometry of the cell.
    BaseType::initialize(patch_storage);

    if (patch_storage->n_patches() == 0)
      return;

    FEEval fe_eval(*patch_storage->get_matrix_free());

    // The patch operator is constructed using a tensor product of 1D matrices.
    // The deal.II FE_Q elements use a hierarchical ordering of the Lagrangian
    // basis functions, but the tensor product construction assumes a
    // lexicographic ordering of degrees of freedom. We compute a permutation
    // vector that maps from the lexicographic ordering (used for the matrix) to
    // the hierarchical ordering (used by the FE_Q element). This vector is then
    // passed to the matrix creation functions.
    std::vector<unsigned int> numbering =
      FETools::lexicographic_to_hierarchic_numbering<1>(fe_degree);

    // Get the finite element being used on the full-dimensional mesh.
    const auto &fe =
      patch_storage->get_matrix_free()->get_dof_handler().get_fe();

    // From the full-dimensional FE, we create its 1D equivalent. For example,
    // if the original is FE_Q<dim>, this creates an FE_Q<1>. This is done by
    // getting the name of the FE as a string and replacing the dimension.
    std::string name = fe.get_name();
    name.replace(name.find('<') + 1, 1, "1");
    std::unique_ptr<FiniteElement<1>> fe_1d = FETools::get_fe_by_name<1>(name);

    // Reinitialize the FEEvaluation object for the first cell (index 0) to
    // access its geometric properties. We assume the mesh is uniform, so the
    // properties of one cell are representative of all others.
    fe_eval.reinit(0);

    // Get the inverse of the Jacobian matrix for this cell. The Jacobian maps
    // from the reference cell to the real cell in physical space.
    auto inverse_jacobian = fe_eval.inverse_jacobian(0);
    // This implementation assumes a simple Cartesian grid where the coordinate
    // axes are aligned with the grid lines. We verify this by checking that
    // the off-diagonal entries of the inverse Jacobian are zero.
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        if (d != e)
          AssertThrow(inverse_jacobian[d][e][0] == 0., ExcNotImplemented());

    // Assuming a uniform mesh, we can determine the mesh size 'h' from the
    // first diagonal element of the inverse Jacobian.
    const number h = inverse_jacobian[0][0][0];



    // Now, we create the 1D matrices for a
    // single cell. These are the fundamental building blocks.
    // `create_1d_cell_mass_matrix` assembles the matrix for int(v*u) dx.
    // {true true} indicates that we want to include the right and left boundary
    // DoFs in the matrix.
    FullMatrix<number> cell_mass_matrix =
      TensorProductMatrixCreator::create_1d_cell_mass_matrix<number>(
        *fe_1d, h, {true, true}, numbering);
    // `create_1d_cell_laplace_matrix` assembles the matrix for
    // int(grad(v)*grad(u)) dx.
    FullMatrix<number> cell_laplace_matrix =
      TensorProductMatrixCreator::create_1d_cell_laplace_matrix<number>(
        *fe_1d, h, {true, true}, numbering);


    // A patch consists of 2x1 cells in 1D. We now assemble the 1D *patch*
    // matrices from the 1D *cell* matrices we just created. This helper
    // function combines the cell matrices to form the larger matrix that
    // represents the operator over the entire 1D patch.
    // The patch inverse operator only acts on the interior DoFs, so we
    // create the patch matrices with `false` for the boundary DoFs.
    Table<2, number> patch_mass_matrix =
      TensorProductMatrixCreator::create_1D_discretization_matrix<number>(
        cell_mass_matrix, 2, 1, {false, false});

    Table<2, number> patch_laplace_matrix =
      TensorProductMatrixCreator::create_1D_discretization_matrix<number>(
        cell_laplace_matrix, 2, 1, {false, false});

    // Finally, we initialize the patch inverse operator. This object will be
    // used during the smoothing steps to solve the local problem on each
    // patch. It takes the 1D mass and Laplace matrices and uses them to
    // construct the full `dim`-dimensional inverse operator via tensor
    // products.
    patch_inverse_operator.reinit(patch_mass_matrix, patch_laplace_matrix);
  }

  /**
   * @brief Clear implementation.
   *
   * Calls clear on the base class and the `inverse` member.
   */
  template <int dim, int fe_degree, typename number>
  void LaplacePatchSmoother<dim, fe_degree, number>::clear()
  {
    BaseType::clear();
  }

  /**
   * @brief Local apply implementation.
   *
   * This function performs one smoother step (forward or backward Gauss-Seidel
   * type iteration over patches) for a given range of patches. It calculates
   * the local residual on each patch, applies the inverse of the local patch
   * operator to compute a correction, and updates the global solution vector.
   *
   * Uses  FEPatchEvaluation to manage data transfer between
   * global vectors and local patch vectors. Iterates through the specified
   * patch range, computes the local residual, applies the patch inverse, and
   * updates the global solution vector. At the  moment, only do_forward==true
   * is supported by PatchStorage that manages the loop.
   */
  template <int dim, int fe_degree, typename number>
  void LaplacePatchSmoother<dim, fe_degree, number>::local_apply(
    const PatchStorageType &patch_storage,
    VectorType             &dst, // Global solution vector (input/output)
    const VectorType       &src, // Global right-hand side vector (input)
    const typename PatchStorageType::PatchRange
               &patch_range,      // Range of patches to process
    const bool &do_forward) const // Direction of iteration
  {
    // Type alias for the FEPatchEvaluation, which handles operations across
    // the cells within a patch, including data gathering and scattering.
    // Moving values between patches vector and individual cell vectors, handled
    // by the last template parameter, is crucial for efficiency. For the
    // reference, all operations with distributor take around 10% of smoothing
    // time, that is a comparable amount of time to inverse application via
    // TensorProductMatrixSymmetricSum. If the residual is also  computed via
    // TensorProductMatrixSymmetricSum this fraction is significantly higher.
    //
    // Micro-benchmarking can help determine the optimal method for
    // distributing DoFs, as performance depends on several factors. For
    // instance, with Clang, lookup tables are highly efficient; the overhead
    // of copying data between cell vectors and patch vectors is minimal,
    // comparable to a standard copy operation of the same data size. With GCC,
    // an index-computing implementation (not included here) outperforms lookup
    // tables, though it's still significantly slower than Clang. Using
    // `PatchDistributors::Lookup` is crucial for efficiency, as it provides the
    // mapping between patch DoFs and cell DoFs.
    using PatchEval =
      FEPatchEvaluation<FEEval, PatchDistributors::Lookup<dim, fe_degree>>;



    // Set up the FEPatchEvaluation object. It takes the patch storage
    // and a configured FEEvaluation object (based on the MatrixFree data
    // associated with the patch storage).
    PatchEval patch_eval(patch_storage,
                         FEEval(*patch_storage.get_matrix_free()));

    // Determine the start, end, and step direction for the patch loop based
    // on whether it's a forward or backward sweep.
    const auto begin = do_forward ? patch_range.first : patch_range.second - 1;
    const auto end   = do_forward ? patch_range.second : patch_range.first - 1;
    const auto step  = do_forward ? 1 : -1;


    // Iterate over the assigned range of patches.
    for (auto patch_index = begin; patch_index != end; patch_index += step)
      {
        // Reinitialize the patch evaluation helper for the current patch index.
        // This sets up internal structures for the specific cells in this
        // patch.
        patch_eval.reinit(patch_index);

        // Allocate aligned vectors for local computations on the patch.
        // The size is determined by the number of unique DoFs within the patch.
        AlignedVector<number> local_rhs(patch_eval.n_patch_dofs());
        AlignedVector<number> local_residual(patch_eval.n_patch_dofs());
        AlignedVector<number> local_correction(patch_eval.n_patch_dofs());

        // --- Step 1: Gather the right-hand side (RHS) for the patch ---
        // Read the relevant global RHS values ('src') into the FEEvaluation
        // objects associated with the cells of the current patch.
        patch_eval.read_dof_values(src);
        // Gather the values from the individual cells' FEEvaluation objects
        // into the contiguous 'local_rhs' vector for the patch.
        // 'false' indicates this is a simple gather, we are not accumulating
        // results In other words, values that are share between cell will be
        // imported once,only for the cell with lowest index (cells are ordered
        // lexicographically).
        patch_eval.gather_local_to_patch(ArrayView<double>(local_rhs), false);

        // --- Step 2: Apply the local operator A_patch * u_patch ---
        // Read the relevant global solution values ('dst') into the
        // FEEvaluation objects.
        patch_eval.read_dof_values(dst);
        // Loop through the FEEvaluation objects (one for each cell batch in the
        // patch).
        for (auto &phi : patch_eval.fe_evaluations)
          {
            // Evaluate gradients of the current solution at quadrature points.
            phi.evaluate(EvaluationFlags::gradients);
            // Submit these gradients for integration (standard Laplace weak
            // form).
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            // Perform the integration (test function gradients against
            // submitted gradients).
            phi.integrate(EvaluationFlags::gradients);
            // The results of the local operator application (A_cell * u_cell)
            // are now stored within each 'phi' object.
          }

        // Gather the results of the local operator application from the cells
        // into the 'local_residual' vector for the patch.
        // 'true' indicates accumulation is needed if multiple cell DoFs map
        // to the same patch DoF (though typically not the case here).
        patch_eval.gather_local_to_patch(ArrayView<double>(local_residual),
                                         true);

        // --- Step 3: Compute the actual residual on the patch ---
        // residual = rhs - A*solution (Note: local_residual currently holds
        // A*solution)
        for (unsigned int i = 0; i < local_residual.size(); ++i)
          local_residual[i] = local_rhs[i] - local_residual[i];


        // --- Step 4: Apply the inverse patch operator ---
        // Solve the local system A_patch * correction = residual for the
        // correction. 'inverse' is the PatchTensorInverse object initialized
        // earlier.

        patch_inverse_operator.apply_inverse(
          ArrayView<number>(local_correction.data(), n_dofs_interior),
          ArrayView<const number>(local_residual.data(), n_dofs_interior));

        // inverse.vmult(local_correction, local_residual);

        // --- Step 5: Update the solution ---
        // The standard smoother update is: u_new = u_old + correction
        // However, the distribute functions from underlying FEEvaluations
        // perform an *addition* to the target vector. Since there are DoF that
        // are shared on multiple cells and we want to add the correction only
        // once, we need zero out the "duplicates".
        //
        // Distribute the patch correction vector back to the individual cells'
        // FEEvaluation objects. 'false' indicates that DoFs shared between
        // multiple cell are only written the first owning cell, not
        // accumulated. The values of shared DoFs in subsequent cells will be
        // set to 0.
        patch_eval.distribute_patch_to_local(
          ArrayView<const number>(local_correction), false);

        // Add the corrections stored in the local FEEvaluation objects
        // to the global solution vector 'dst'.
        patch_eval.distribute_local_to_global(dst);
      }
  }


} // namespace Operators

// @sect3{The main class}

namespace LaplaceSolver
{
  using namespace dealii;

  /**
   * This class implements a Laplace solver using matrix-free techniques.
   * There are two main differences between this class and the one from step-37:
   *  - MatrixFree object for level operators  have to enable evaluation on
   * ghosted cells
   *  - The smoother is a patch smoother, which is more efficient for high order
   * elements.  It requires building PatchStorage to store information necessary
   * for loops. The smoother is good enough on its own, there is no need to
   * combine it with PreconditionChebyshev
   *
   */
  template <int dim, int fe_degree>
  class LaplaceSolver
  {
  public:
    using Number              = double;
    using LevelNumber         = double;
    using VectorType          = LinearAlgebra::distributed::Vector<Number>;
    using LevelVectorType     = LinearAlgebra::distributed::Vector<LevelNumber>;
    using LevelMatrixFreeType = MatrixFree<dim, LevelNumber>;
    LaplaceSolver();
    void initialize(const unsigned int n_refinements = 4);

    void solve();

    void output_results(const unsigned int cycle) const;

  private:
    void setup_system();
    void assemble_rhs();

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif

    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;

    MappingQ1<dim> mapping;

    AffineConstraints<Number> constraints;
    using SystemMatrixType = Operators::LaplaceOperator<dim, fe_degree, Number>;
    SystemMatrixType system_matrix;

    MGConstrainedDoFs mg_constrained_dofs;
    using LevelMatrixType =
      Operators::LaplaceOperator<dim, fe_degree, LevelNumber>;
    MGLevelObject<LevelMatrixType> mg_matrices;

    VectorType solution;
    VectorType system_rhs;

    double             setup_time;
    ConditionalOStream pcout;
    ConditionalOStream time_details;
  };



  template <int dim, int fe_degree>
  LaplaceSolver<dim, fe_degree>::LaplaceSolver()
    :
#ifdef DEAL_II_WITH_P4EST
    triangulation(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
    ,
#else
    triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    ,
#endif
    fe(fe_degree)
    , dof_handler(triangulation)
    , setup_time(0.)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    , time_details(std::cout,
                   Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {}


  template <int dim, int fe_degree>
  void LaplaceSolver<dim, fe_degree>::setup_system()
  {
    Timer time;
    setup_time = 0;

    system_matrix.clear();
    mg_matrices.clear_elements();

    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

    IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs, locally_owned_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();
    setup_time += time.wall_time();
    time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
    time.restart();

    {
      typename MatrixFree<dim, Number>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, Number>::AdditionalData::none;
      additional_data.mapping_update_flags =
        (update_gradients | update_JxW_values | update_quadrature_points);
      std::shared_ptr<MatrixFree<dim, Number>> system_mf_storage(
        new MatrixFree<dim, Number>());
      system_mf_storage->reinit(mapping,
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);


      system_mf_storage->reinit(mapping,
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);

      system_matrix.initialize(system_mf_storage);
    }

    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(system_rhs);

    setup_time += time.wall_time();
    time_details << "Setup matrix-free system   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
    time.restart();

    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels - 1);

    std::set<types::boundary_id> dirichlet_boundary;
    dirichlet_boundary.insert(0);
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                       dirichlet_boundary);

    for (unsigned int level = 0; level < nlevels; ++level)
      {
        IndexSet locally_owned_dofs = dof_handler.locally_owned_mg_dofs(level);
        IndexSet relevant_dofs =
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
        AffineConstraints<Number> level_constraints;
        level_constraints.reinit(locally_owned_dofs, relevant_dofs);
        level_constraints.add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        level_constraints.close();

        typename MatrixFree<dim, LevelNumber>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, LevelNumber>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        additional_data.mg_level = level;

        // We need to set the level of the matrix-free object to
        // enable evaluation on ghosted cell.
        additional_data.store_ghost_cells = true;
        std::shared_ptr<MatrixFree<dim, LevelNumber>> mg_mf_storage_level(
          new MatrixFree<dim, LevelNumber>());
        mg_mf_storage_level->reinit(mapping,
                                    dof_handler,
                                    level_constraints,
                                    QGauss<1>(fe.degree + 1),
                                    additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level,
                                      mg_constrained_dofs,
                                      level);
      }
    setup_time += time.wall_time();
    time_details << "Setup matrix-free levels   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
  }



  template <int dim, int fe_degree>
  void LaplaceSolver<dim, fe_degree>::assemble_rhs()
  {
    Timer time;

    system_rhs = 0;
    FEEvaluation<dim, fe_degree> phi(*system_matrix.get_matrix_free());
    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches();
         ++cell)
      {
        phi.reinit(cell);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_value(make_vectorized_array<Number>(1.0), q);
        phi.integrate(EvaluationFlags::values);
        phi.distribute_local_to_global(system_rhs);
      }
    system_rhs.compress(VectorOperation::add);

    setup_time += time.wall_time();
    time_details << "Assemble right hand side   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
  }

  template <int dim, int fe_degree>
  void LaplaceSolver<dim, fe_degree>::solve()
  {
    using StorageType = PatchStorage<MatrixFree<dim, double>>;

    // Get the number of levels in the multigrid hierarchy.
    const unsigned max_lvl = triangulation.n_global_levels() - 1;

    // Initialize the transfer operator between multigrid levels.
    MGTransferMatrixFree<dim, LevelNumber> mg_transfer(mg_constrained_dofs);
    mg_transfer.build(dof_handler);

    // Define the smoother type based on the LaplacePatchSmoother.
    using SmootherType =
      Operators::LaplacePatchSmoother<dim, fe_degree, double>;
    // Initialize the relaxation smoother object for multigrid.
    mg::SmootherRelaxation<SmootherType, LevelVectorType> mg_smoother;

    // Initialize the patch smoother for each level.
    mg_smoother.resize(0, max_lvl);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        // Create patch storage for the current level using the matrix-free
        // data.
        std::shared_ptr<StorageType> storage =
          std::make_shared<StorageType>(mg_matrices[level].get_matrix_free());
        // Initialize the patch storage,  this generates patches
        //  and initializes loop structures
        storage->initialize();

        // Initialize the smoother for the current level with the patch storage.
        mg_smoother[level].initialize(storage);
      }

    // Optionally output patch information for debugging/visualization.
    mg_smoother[max_lvl].get_storage()->output_patches("patches");
    mg_smoother[max_lvl].get_storage()->output_centerpoints("centerpoints");

    // One smoothing step is enough, we are using an effective but also
    // expesinsive smoother
    mg_smoother.set_steps(1);

    // Initialize the coarse grid solver using the smoother defined for level 0.
    MGCoarseGridApplySmoother<LevelVectorType> mg_coarse;
    mg_coarse.initialize(mg_smoother);

    // Create the multigrid matrix object using the level matrices.
    mg::Matrix<LevelVectorType> mg_matrix(mg_matrices);

    // Standard  Multigrid initialization.
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<LevelVectorType> mg_interface(mg_interface_matrices);
    Multigrid<LevelVectorType>  mg(
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);

    // Initialize the multigrid preconditioner.
    PreconditionMG<dim, LevelVectorType, MGTransferMatrixFree<dim, LevelNumber>>
      preconditioner(dof_handler, mg, mg_transfer);

    // Initialize the Conjugate Gradient solver.
    SolverControl        solver_control(100, 1e-8 * system_rhs.l2_norm());
    SolverCG<VectorType> cg(solver_control);

    // Start timer for the solve phase.
    Timer time;
    time.restart();
    // Ensure the initial solution satisfies boundary conditions.
    constraints.set_zero(solution);

    // Solve the linear system Ax = b using CG with the multigrid
    // preconditioner.
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
    // Distribute hanging node constraints to the final solution.
    constraints.distribute(solution);

    // Print timing information and iteration count.
    pcout << "Solver finished after " << solver_control.last_step()
          << " iterations." << std::endl;
    time_details << "Time solve                 (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
  }



  template <int dim, int fe_degree>
  void
  LaplaceSolver<dim, fe_degree>::output_results(const unsigned int cycle) const
  {
    Timer time;
    if (triangulation.n_global_active_cells() > 1000000)
      return;

    DataOut<dim> data_out;

    solution.update_ghost_values();
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches(mapping, fe_degree);

    DataOutBase::VtkFlags flags;
    data_out.set_flags(flags);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, MPI_COMM_WORLD, 3);

    time_details << "Time write output          (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s\n";
  }


  template <int dim, int fe_degree>
  void
  LaplaceSolver<dim, fe_degree>::initialize(const unsigned int n_refinements)
  {
    {
      const unsigned int n_vect_doubles = VectorizedArray<double>::size();
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

      pcout << "Vectorization over " << n_vect_doubles
            << " doubles = " << n_vect_bits << " bits ("
            << Utilities::System::get_current_vectorization_level() << ")"
            << std::endl;
    }


    GridGenerator::subdivided_hyper_cube(triangulation, 2);

    triangulation.refine_global(n_refinements);
    setup_system();
    assemble_rhs();
    pcout << std::endl;
  }


} // namespace LaplaceSolver

int main(int argc, char *argv[])
{
  using namespace dealii;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  const unsigned int dim       = 2;
  const unsigned int fe_degree = 3;

  LaplaceSolver::LaplaceSolver<dim, fe_degree> solver;
  solver.initialize(2);
  solver.solve();
  solver.output_results(0);

  return 0;
}
