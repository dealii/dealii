// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_ginkgo_solver_h
#  define dealii_ginkgo_solver_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_GINKGO

#    include <deal.II/lac/block_sparse_matrix.h>
#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/solver_control.h>
#    include <deal.II/lac/sparse_matrix.h>
#    include <deal.II/lac/vector.h>

#    include <ginkgo/ginkgo.hpp>

#    include <memory>

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
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType, typename IndexType>
  class SolverBase
  {
  public:
    /**
     * Constructor.
     *
     * The @p executor defines the paradigm where the solution is computed.
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
    SolverBase(SolverControl &                solver_control,
               std::shared_ptr<gko::Executor> executor);

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
          Vector<ValueType> &            solution,
          const Vector<ValueType> &      rhs);

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
  };


  /**
   * An implementation of the solver interface using the Ginkgo CG solver.
   *
   * @ingroup GinkgoWrappers
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
     * @p solver_control The solver control object is then used to set the
     * parameters and setup the CG solver from the CG factory which solves the
     * linear system.
     *
     * @p executor The execution paradigm for the CG solver.
     */
    SolverCG(SolverControl &                solver_control,
             std::shared_ptr<gko::Executor> executor,
             const AdditionalData &         data = AdditionalData());

  /**
   * Constructor.
   *
   * @p solver_control The solver control object is then used to set the
   * parameters and setup the CG solver from the CG factory which solves the
   * linear system.
   *
   * @p executor The execution paradigm for the CG solver.
   */
  SolverCG(SolverControl &                solver_control,
           std::shared_ptr<gko::Executor> executor,
           std::shared_ptr<gko::LinOpFactory> preconditioner,
           const AdditionalData &         data = AdditionalData());
  protected:
    /**
     * Store a copy of the settings for this particular solver.
     */
    const AdditionalData additional_data;
  };


} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_GINKGO

#endif
/*----------------------------   ginkgo_solver.h ---------------------------*/
