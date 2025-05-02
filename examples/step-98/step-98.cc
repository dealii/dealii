/* ------------------------------------------------------------------------
 *
 *
 * Authors: Michal Wichrowski, Heidelberg University,
 * 2025
 */


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



#include <fstream>
#include <iostream>



namespace Operators
{
  using namespace dealii;

  /**
   *
   * This class implements the Laplace operator $-\Delta u$ as in step-37 using
   * matrix-free techniques based on the deal.II::MatrixFree framework. It
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
    LaplaceOperator();

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
    virtual void compute_diagonal() override;

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

    /**
     * @brief Local worker routine for computing the diagonal entries.
     *
     * This function is called by MatrixFree::cell_loop during the diagonal
     * computation. It computes the diagonal entries for the cells within the
     * specified range by applying the local operator to unit vectors.
     *
     * @param data The MatrixFree object holding cached data.
     * @param dst The vector where the diagonal entries are accumulated.
     * @param dummy A placeholder argument required by the cell_loop interface
     *              when computing the diagonal.
     * @param cell_range The pair of start and end indices for the cell batches
     *                   to be processed.
     */
    void local_compute_diagonal(
      const MatrixFree<dim, number>               &data,
      VectorType                                  &dst,
      const unsigned int                          &dummy,
      const std::pair<unsigned int, unsigned int> &cell_range) const;
  };

  /**
   * @brief Constructor implementation.
   *
   * Initializes the base class.
   */
  template <int dim, int fe_degree, typename number>
  LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
    : MatrixFreeOperators::Base<dim, VectorType>()
  {}


  /**
   * @brief Clear implementation.
   *
   * Calls the base class clear function.
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
   * 4. Submit the gradients (multiplied by coefficient, if any - here just 1.0)
   *    back for integration.
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
   * @brief Compute diagonal implementation, as in step-37.
   *
   */
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(new DiagonalMatrix<VectorType>());
    VectorType &inverse_diagonal = this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    unsigned int dummy = 0;
    this->data->cell_loop(&LaplaceOperator::local_compute_diagonal,
                          this,
                          inverse_diagonal,
                          dummy);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1. / inverse_diagonal.local_element(i);
      }
  }

  /**
   * Computes the diagonal entries for a range of cell batches, as in step-37.
   */
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
    const MatrixFree<dim, number> &data,
    VectorType                    &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

    AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.submit_dof_value(VectorizedArray<number>(), j);
            phi.submit_dof_value(make_vectorized_array<number>(1.), i);

            phi.evaluate(EvaluationFlags::gradients);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(EvaluationFlags::gradients);
            diagonal[i] = phi.get_dof_value(i);
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.submit_dof_value(diagonal[i], i);
        phi.distribute_local_to_global(dst);
      }
  }



  /**
   * @brief Implements the inverse of the local patch operator using tensor
   * product techniques.
   *
   * This class computes and applies the inverse of the Laplace operator
   * assembled on a single patch (typically consisting of $2^d$ cells). It
   * exploits the tensor product structure of the finite element basis and
   * quadrature on the reference patch to efficiently compute the inverse.
   *
   * The patch operator is assumed to be a sum of Kronecker products involving
   * 1D mass and stiffness matrices. This class pre-computes these 1D matrices
   * and uses the TensorProductMatrixSymmetricSum class to apply the inverse.
   *
   */
  template <int dim, int fe_degree, typename number, typename VectorizedValues>
  class PatchTensorInverse
  {
  public:
    //  using value_type = number;


    /**
     * @brief The size (number of rows/columns) of the 1D matrices that form
     * the tensor product patch matrix. For a patch of $2^d$ cells with DGQ
     * elements of degree `fe_degree`, this corresponds to the number of unique
     * 1D basis functions on the two adjacent reference cells along one axis.
     */
    const constexpr static unsigned int n_rows_1d = 2 * fe_degree - 1;

    /**
     * @brief Clear internal data structures.
     */
    void clear()
    {}

    /**
     * @brief Initialize the patch inverse operator.
     *
     * This function computes the 1D mass and stiffness matrices on the
     * reference cell, scales them according to the mapping information
     * obtained from the MatrixFree object, assembles the patch matrices
     * (patch_mass, patch_laplace) by combining contributions from adjacent
     * cells, and initializes the TensorProductMatrixSymmetricSum object used
     * for applying the inverse.
     *
     * @param mf A shared pointer to the MatrixFree object containing geometric
     *           and finite element information. Assumes a uniform grid patch.
     */
    void initialize(std::shared_ptr<const MatrixFree<dim, number>> mf);

    /**
     * @brief Apply the inverse of the patch matrix to a vector.
     *
     * Solves the linear system `patch_matrix * dst = src` for `dst`.
     *
     * @tparam VectorType The type of the source and destination vectors
     *                    (typically AlignedVector<VectorizedValues>).
     * @param dst The destination vector (solution).
     * @param src The source vector (right-hand side).
     */
    template <typename VectorType>
    void vmult(VectorType &dst, const VectorType &src) const;

    /**
     * @brief Apply the transpose of the inverse patch matrix (equivalent to
     * vmult for symmetric matrices).
     *
     * @param dst The destination vector.
     * @param src The source vector.
     * @param tmp_memory Temporary storage (unused in this implementation).
     */
    void Tvmult(LinearAlgebra::distributed::Vector<number>       &dst,
                const LinearAlgebra::distributed::Vector<number> &src,
                AlignedVector<VectorizedValues> &tmp_memory) const
    {
      vmult(dst, src, tmp_memory);
    }

  private:
    /**
     * @brief Shared pointer to the MatrixFree data associated with this patch
     * operator. Used to extract mapping information.
     */
    std::shared_ptr<const MatrixFree<dim, number>> data;

    /**
     * @brief The tensor product matrix representing the patch operator.
     * Provides the `apply_inverse` method.
     */
    TensorProductMatrixSymmetricSum<dim, VectorizedValues, n_rows_1d>
      patch_matrix;

    /**
     * @brief Array storing the 1D mass matrices for each dimension, assembled
     * over the two reference cells forming the patch dimension.
     */
    std::array<Table<2, VectorizedValues>, dim> patch_mass;
    /**
     * @brief Array storing the 1D Laplace matrices for each dimension,
     * assembled over the two reference cells forming the patch dimension and
     * scaled by mapping factors.
     */
    std::array<Table<2, VectorizedValues>, dim> patch_laplace;
  };

  /**
   * @brief Initialize implementation.
   *
   * Computes 1D reference mass and stiffness matrices, scales them using
   * Jacobian information from the first cell batch in the MatrixFree object
   * (assuming a uniform grid), assembles the 1D patch matrices by summing
   * contributions from the two cells along each dimension, and initializes
   * the `patch_matrix` member.
   *
   * The implementation should be replaced by calling function from
   * https://github.com/dealii/dealii/pull/18361
   */
  template <int dim, int fe_degree, typename number, typename VectorizedValues>
  void PatchTensorInverse<dim, fe_degree, number, VectorizedValues>::initialize(
    std::shared_ptr<const MatrixFree<dim, number>> mf)
  {
    data = mf;

    auto fe_1d = std::make_unique<FE_DGQ<1>>(fe_degree);

    const unsigned int                          N = fe_degree + 1;
    FullMatrix<number>                          laplace_unscaled(N, N);
    std::array<Table<2, VectorizedValues>, dim> mass_matrices;
    std::array<Table<2, VectorizedValues>, dim> laplace_matrices;
    for (unsigned int d = 0; d < dim; ++d)
      {
        mass_matrices[d].reinit(N, N);
        laplace_matrices[d].reinit(N, N);
      }

    QGauss<1> quadrature(N);
    // assemble 1d matrices on reference cell
    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < N; ++j)
        {
          double sum_mass = 0, sum_laplace = 0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              sum_mass += (fe_1d->shape_value(i, quadrature.point(q)) *
                           fe_1d->shape_value(j, quadrature.point(q))) *
                          quadrature.weight(q);
              sum_laplace += (fe_1d->shape_grad(i, quadrature.point(q))[0] *
                              fe_1d->shape_grad(j, quadrature.point(q))[0]) *
                             quadrature.weight(q);
            }
          for (unsigned int d = 0; d < dim; ++d)
            mass_matrices[d](i, j) = sum_mass;

          laplace_unscaled(i, j) = sum_laplace;
        }

    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data);

    phi.reinit(0);

    auto inverse_jacobian = phi.inverse_jacobian(0);

    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        if (d != e)
          AssertThrow(inverse_jacobian[d][e][0] == 0., ExcNotImplemented());

    // number -> VectorizedArray<Number> for vectorization
    // Fixme ? not sure if number ...
    number jacobian_determinant = inverse_jacobian[0][0][0];
    for (unsigned int e = 1; e < dim; ++e)
      jacobian_determinant *= inverse_jacobian[e][e][0];
    jacobian_determinant = 1. / jacobian_determinant;


    for (unsigned int d = 0; d < dim; ++d)
      {
        // FIXME later:
        const number scaling_factor = inverse_jacobian[d][d][0] *
                                      inverse_jacobian[d][d][0] *
                                      jacobian_determinant;

        for (unsigned int i = 0; i < N; ++i)
          for (unsigned int j = 0; j < N; ++j)
            laplace_matrices[d](i, j) =
              scaling_factor * VectorizedValues(laplace_unscaled(i, j));
      }

    const unsigned int NN = n_rows_1d;

    for (unsigned int d = 0; d < dim; ++d)
      {
        patch_mass[d].reinit(NN, NN);
        patch_laplace[d].reinit(NN, NN);
      }

    // assemble tensor product on one patch
    for (unsigned int p = 0; p < 2; ++p)
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int i = 0; i < fe_degree; ++i)
          for (unsigned int j = 0; j < fe_degree; ++j)
            {
              patch_mass[d](i + p * (fe_degree - 1), j + p * (fe_degree - 1)) +=
                mass_matrices[d](i + 1 - p, j + 1 - p);
              patch_laplace[d](i + p * (fe_degree - 1),
                               j + p * (fe_degree - 1)) +=
                laplace_matrices[d](i + 1 - p, j + 1 - p);
            }

    patch_matrix.reinit(patch_mass, patch_laplace);

#ifdef DEBUG_INVERSE
    for (unsigned int i = 0; i < NN; ++i)
      {
        for (unsigned int j = 0; j < NN; ++j)
          std::cout << patch_mass[0](i, j) << " ";
      }
    std::cout << std::endl;
#endif
  }

  /**
   * @brief vmult implementation.
   *
   * Calls the `apply_inverse` method of the underlying
   * TensorProductMatrixSymmetricSum object.
   */
  template <int dim, int fe_degree, typename number, typename VectorizedValues>
  template <typename VectorType>
  void PatchTensorInverse<dim, fe_degree, number, VectorizedValues>::vmult(
    VectorType       &dst,
    const VectorType &src) const
  {
    const unsigned int n_dofs_interior =
      Utilities::pow<unsigned int>(n_rows_1d, dim);
    patch_matrix.apply_inverse(
      ArrayView<VectorizedValues>(dst.data(), n_dofs_interior),
      ArrayView<const VectorizedValues>(src.data(), n_dofs_interior));
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
  class LaplacePatchSmoother : public SmootherBase<dim, number>
  {
  public:
    using BaseType       = SmootherBase<dim, number>;
    using VectorType     = LinearAlgebra::distributed::Vector<number>;
    using MatrixFreeType = MatrixFree<dim, number>;
    /**
     * @brief Type alias for the storage managing patch data (DoF indices, etc.).
     */
    using PatchStorageType = PatchStorage<MatrixFreeType>;
    /**
     * @brief Type alias for the class providing the local patch inverse.
     */
    using InverseType = PatchTensorInverse<dim, fe_degree, number, number>;

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
     * 1. Gathers the relevant portion of the global right-hand side vector
     * `src` to a local patch vector `local_rhs`.
     * 2. Applies the global operator (implicitly via FEEvaluation) to the
     *    current global solution estimate `dst` and gathers the result to a
     *    local patch vector `local_residual`.
     * 3. Computes the actual residual on the patch: `local_residual = local_rhs
     * - A_patch * current_solution_on_patch`.
     * 4. Applies the inverse patch operator (`inverse.vmult`) to the residual
     *    to get a local correction: `local_correction = A_patch^{-1} *
     * local_residual`.
     * 5. Updates the global solution `dst` by distributing the negative of the
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
     * @brief Initialize the smoother.
     *
     * Initializes the base class and the local patch inverse operator.
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
     * @brief The object responsible for applying the inverse of the operator
     * on a single patch.
     */
    InverseType inverse;
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
    BaseType::initialize(patch_storage);
    if (patch_storage->n_patches() != 0)
      inverse.initialize(patch_storage->get_matrix_free());
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
    inverse.clear();
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
    // Type alias for the standard FEEvaluation on a single cell batch.
    using FEEval = FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number>;


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
    // `DistributorLookup` is crucial for efficiency, as it provides the
    // mapping between patch DoFs and cell DoFs.
    using PatchEval =
      FEPatchEvaluation<FEEval,
                        PatchDistributors::DistributorLookup<dim, fe_degree>>;



    // Instantiate the FEPatchEvaluation object. It takes the patch storage
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
        inverse.vmult(local_correction, local_residual);

        // --- Step 5: Update the solution ---
        // The standard smoother update is: u_new = u_old + correction
        // However, the distribute functions from underlying FEEvaluations
        // perform an *addition* to the target vector. So, we negate the
        // correction and add it. dst_new = dst_old +
        // (-correction) is equivalent to dst_new = dst_old - correction
        //
        // FIXME: Gemini thinks that minus here is a mistake and to be fair it
        // has a point. It is not straighforward to explain why we need to
        // negate the correction. I am sure that his version is correct
        // (veryfied), so explanation is needed.

        // Negate the computed correction before distributing.
        for (unsigned int i = 0; i < local_correction.size(); ++i)
          local_correction[i] *= -1;

        // Distribute the patch correction vector back to the individual cells'
        // FEEvaluation objects. 'false' indicates that DoFs shared beteween
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


namespace LaplaceSolver
{
  using namespace dealii;


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

    void print_timings() const;

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


    mutable std::vector<std::vector<
      std::pair<double, std::chrono::time_point<std::chrono::system_clock>>>>
      all_mg_timers;
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
                   false &&
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

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
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
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                      level,
                                                      relevant_dofs);
        AffineConstraints<Number> level_constraints;
        level_constraints.reinit(relevant_dofs);
        level_constraints.add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        level_constraints.close();

        typename MatrixFree<dim, LevelNumber>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, LevelNumber>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        additional_data.mg_level          = level;
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

        mg_matrices[level].compute_diagonal();
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
    pcout << "Time solve (" << solver_control.last_step() << " iterations)"
          << (solver_control.last_step() < 10 ? "  " : " ") << "(CPU/wall) "
          << time.cpu_time() << "s/" << time.wall_time() << "s\n";
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
    // flags.compression_level = DataOutBase::VtkFlags::best_speed;
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

      // pcout << "Running with  "
      // << n_refinements
      // << "  refinements and task_chunk_size = " << task_chunk_size
      // << std::endl;
    }


    {
      GridGenerator::subdivided_hyper_cube(triangulation, 2);
    }
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
  const unsigned int fe_degree = 1;

  LaplaceSolver::LaplaceSolver<dim, fe_degree> solver;
  solver.initialize(2);
  solver.solve();
  solver.output_results(0);

  return 0;
}
