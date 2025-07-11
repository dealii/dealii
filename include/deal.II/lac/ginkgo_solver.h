// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_ginkgo_solver_h
#define dealii_ginkgo_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/lac/block_sparse_matrix.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/solver_control.h>
#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <ginkgo/ginkgo.hpp>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  /**
   * This class forms the base class for all of Ginkgo's iterative solvers.
   * The various derived classes only take
   * the additional data that is specific to them and solve the given linear
   * system. The entire collection of solvers that Ginkgo implements is
   * available at <a Ginkgo
   * href="https://ginkgo-project.github.io/ginkgo/doc/develop/"> documentation
   * and manual pages</a>.
   */
  template <typename ValueType, typename IndexType>
  class SolverBase
  {
  public:
    /**
     * Constructor.
     *
     * The @p exec_type defines the paradigm where the solution is computed.
     * It is a string and the choices are "omp" , "reference" or "cuda".
     * The respective strings create the respective executors as given below.
     *
     * Ginkgo currently supports three different executor types:
     *
     * +    OmpExecutor specifies that the data should be stored and the
     * associated operations executed on an OpenMP-supporting device (e.g. host
     * CPU);
     * ```
     * auto omp = gko::create<gko::OmpExecutor>();
     * ```
     * +    CudaExecutor specifies that the data should be stored and the
     *      operations executed on the NVIDIA GPU accelerator;
     * ```
     * if(gko::CudaExecutor::get_num_devices() > 0 ) {
     *    auto cuda = gko::create<gko::CudaExecutor>();
     * }
     * ```
     * +    ReferenceExecutor executes a non-optimized reference implementation,
     *      which can be used to debug the library.
     * ```
     * auto ref = gko::create<gko::ReferenceExecutor>();
     * ```
     *
     * The following code snippet demonstrates the using of the OpenMP executor
     * to create a solver which would use the OpenMP paradigm to the solve the
     * system on the CPU.
     *
     * ```
     * auto omp = gko::create<gko::OmpExecutor>();
     * using cg = gko::solver::Cg<>;
     * auto solver_gen =
     *     cg::build()
     *          .with_criteria(
     *              gko::stop::Iteration::build().with_max_iters(20u).on(omp),
     *              gko::stop::ResidualNormReduction<>::build()
     *                  .with_reduction_factor(1e-6)
     *                  .on(omp))
     *          .on(omp);
     * auto solver = solver_gen->generate(system_matrix);
     *
     * solver->apply(lend(rhs), lend(solution));
     * ```
     *
     *
     * The @p solver_control object is the same as for other
     * deal.II iterative solvers.
     */
    SolverBase(SolverControl &solver_control, const std::string &exec_type);

    /**
     * Destructor.
     */
    virtual ~SolverBase() = default;

    /**
     * Initialize the matrix and copy over its data to Ginkgo's data structures.
     */
    void
    initialize(const SparseMatrix<ValueType> &matrix);

    /**
     * Solve the linear system <tt>Ax=b</tt>. Dependent on the information
     * provided by derived classes one of Ginkgo's linear solvers is
     * chosen.
     */
    void
    apply(Vector<ValueType> &solution, const Vector<ValueType> &rhs);

    /**
     * Solve the linear system <tt>Ax=b</tt>. Dependent on the information
     * provided by derived classes one of Ginkgo's linear solvers is
     * chosen.
     */
    void
    solve(const SparseMatrix<ValueType> &matrix,
          Vector<ValueType>             &solution,
          const Vector<ValueType>       &rhs);

    /**
     * Access to the object that controls convergence.
     */
    SolverControl &
    control() const;


  protected:
    /**
     * Reference to the object that controls convergence of the iterative
     * solvers.
     */
    SolverControl &solver_control;

    /**
     * The Ginkgo generated solver factory object.
     */
    std::shared_ptr<gko::LinOpFactory> solver_gen;

    /**
     * The residual criterion object that controls the reduction of the residual
     * based on the tolerance set in the solver_control member.
     */
    std::shared_ptr<gko::stop::ResidualNormReduction<>::Factory>
      residual_criterion;

    /**
     * The Ginkgo convergence logger used to check for convergence and other
     * solver data if needed.
     */
    std::shared_ptr<gko::log::Convergence<>> convergence_logger;

    /**
     * The Ginkgo combined factory object is used to create a combined stopping
     * criterion to be passed to the solver.
     */
    std::shared_ptr<gko::stop::Combined::Factory> combined_factory;

    /**
     * The execution paradigm in Ginkgo. The choices are between
     * `gko::OmpExecutor`, `gko::CudaExecutor` and `gko::ReferenceExecutor`
     * and more details can be found in Ginkgo's documentation.
     */
    std::shared_ptr<gko::Executor> executor;

  private:
    /**
     * Initialize the Ginkgo logger object with event masks. Refer to
     * <a
     * href="https://github.com/ginkgo-project/ginkgo/blob/develop/include/ginkgo/core/log/logger.hpp">Ginkgo's
     * logging event masks.</a>
     */
    void
    initialize_ginkgo_log();

    /**
     * Ginkgo matrix data structure. First template parameter is for storing the
     * array of the non-zeros of the matrix. The second is for the row pointers
     * and the column indices.
     *
     * @todo Templatize based on Matrix type.
     */
    std::shared_ptr<gko::matrix::Csr<ValueType, IndexType>> system_matrix;

    /**
     * The execution paradigm as a string to be set by the user. The choices
     * are between `omp`, `cuda` and `reference` and more details can be found
     * in Ginkgo's documentation.
     */
    const std::string exec_type;
  };


  /**
   * An implementation of the solver interface using the Ginkgo CG solver.
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverCG : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * A standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the CG solver from the CG factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the CG solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverCG(SolverControl        &solver_control,
             const std::string    &exec_type,
             const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the CG solver from the CG factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the CG solver.
     *
     * @param[in] preconditioner The preconditioner for the solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverCG(SolverControl                            &solver_control,
             const std::string                        &exec_type,
             const std::shared_ptr<gko::LinOpFactory> &preconditioner,
             const AdditionalData                     &data = AdditionalData());

  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };


  /**
   * An implementation of the solver interface using the Ginkgo Bicgstab solver.
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverBicgstab : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * A standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the Bicgstab solver from the Bicgstab
     * factory which solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the Bicgstab solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverBicgstab(SolverControl        &solver_control,
                   const std::string    &exec_type,
                   const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the Bicgstab solver from the Bicgstab
     * factory which solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the Bicgstab solver.
     *
     * @param[in] preconditioner The preconditioner for the solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverBicgstab(SolverControl                            &solver_control,
                   const std::string                        &exec_type,
                   const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                   const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the solver interface using the Ginkgo CGS solver.
   *
   * CGS or the conjugate gradient square method is an iterative type Krylov
   * subspace method which is suitable for general systems.
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverCGS : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * A standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the CGS solver from the CGS factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the CGS solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverCGS(SolverControl        &solver_control,
              const std::string    &exec_type,
              const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the CGS solver from the CGS factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the CGS solver.
     *
     * @param[in] preconditioner The preconditioner for the solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverCGS(SolverControl                            &solver_control,
              const std::string                        &exec_type,
              const std::shared_ptr<gko::LinOpFactory> &preconditioner,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the solver interface using the Ginkgo FCG solver.
   *
   * FCG or the flexible conjugate gradient method is an iterative type Krylov
   * subspace method which is suitable for symmetric positive definite methods.
   *
   * Though this method performs very well for symmetric positive definite
   * matrices, it is in general not suitable for general matrices.
   *
   * In contrast to the standard CG based on the Polack-Ribiere formula, the
   * flexible CG uses the Fletcher-Reeves formula for creating the orthonormal
   * vectors spanning the Krylov subspace. This increases the computational cost
   * of every Krylov solver iteration but allows for non-constant
   * preconditioners.
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverFCG : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * A standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the FCG solver from the FCG factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the FCG solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverFCG(SolverControl        &solver_control,
              const std::string    &exec_type,
              const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the FCG solver from the FCG factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the FCG solver.
     *
     * @param[in] preconditioner The preconditioner for the solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverFCG(SolverControl                            &solver_control,
              const std::string                        &exec_type,
              const std::shared_ptr<gko::LinOpFactory> &preconditioner,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the solver interface using the Ginkgo GMRES solver.
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverGMRES : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * A standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the number of temporary vectors to 30,
       * i.e. do a restart every 30 iterations.
       */
      AdditionalData(const unsigned int restart_parameter = 30);

      /**
       * Maximum number of tmp vectors.
       */
      unsigned int restart_parameter;
    };

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the GMRES solver from the GMRES factory
     * which solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the GMRES solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverGMRES(SolverControl        &solver_control,
                const std::string    &exec_type,
                const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the GMRES solver from the GMRES factory
     * which solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the GMRES solver.
     *
     * @param[in] preconditioner The preconditioner for the solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverGMRES(SolverControl                            &solver_control,
                const std::string                        &exec_type,
                const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the solver interface using the Ginkgo IR solver.
   *
   * Iterative refinement (IR) is an iterative method that uses another coarse
   * method to approximate the error of the current solution via the current
   * residual.
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverIR : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * A standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the IR solver from the IR factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the IR solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverIR(SolverControl        &solver_control,
             const std::string    &exec_type,
             const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param[in,out] solver_control The solver control object is then used to
     * set the parameters and set up the IR solver from the IR factory which
     * solves the linear system.
     *
     * @param[in] exec_type The execution paradigm for the IR solver.
     *
     * @param[in] inner_solver The Inner solver for the IR solver.
     *
     * @param[in] data The additional data required by the solver.
     */
    SolverIR(SolverControl                            &solver_control,
             const std::string                        &exec_type,
             const std::shared_ptr<gko::LinOpFactory> &inner_solver,
             const AdditionalData                     &data = AdditionalData());

  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };


} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_GINKGO

#endif
