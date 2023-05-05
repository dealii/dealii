// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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
#define dealii_ginkgo_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/lac/block_sparse_matrix.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/ginkgo_preconditioner.h>
#  include <deal.II/lac/ginkgo_sparse_matrix.h>
#  include <deal.II/lac/solver_control.h>
#  include <deal.II/lac/vector.h>

#  include <ginkgo/core/log/convergence.hpp>
#  include <ginkgo/core/solver/bicgstab.hpp>
#  include <ginkgo/core/solver/cg.hpp>
#  include <ginkgo/core/solver/cgs.hpp>
#  include <ginkgo/core/solver/fcg.hpp>
#  include <ginkgo/core/solver/gmres.hpp>
#  include <ginkgo/core/solver/ir.hpp>
#  include <ginkgo/core/stop/iteration.hpp>
#  include <ginkgo/core/stop/residual_norm.hpp>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  namespace detail
  {
    // move to .cc when deprecated constructors are removed
    std::shared_ptr<gko::Executor>
    exec_from_name(const std::string &exec_type)
    {
      if (exec_type == "reference")
        {
          return gko::ReferenceExecutor::create();
        }
      if (exec_type == "omp")
        {
          return gko::OmpExecutor::create();
        }
      if (exec_type == "cuda" && gko::CudaExecutor::get_num_devices() > 0)
        {
          return gko::CudaExecutor::create(0, gko::OmpExecutor::create());
        }
      if (exec_type == "hip" && gko::HipExecutor::get_num_devices() > 0)
        {
          return gko::HipExecutor::create(0, gko::OmpExecutor::create());
        }
      if (exec_type == "dpcpp" &&
          gko::DpcppExecutor::get_num_devices("all") > 0)
        {
          return gko::DpcppExecutor::create(0, gko::OmpExecutor::create());
        }
      Assert(
        false,
        ExcMessage(
          " exec_type needs to be one of the following strings: \"cuda\", \"dpcpp\", \"hip\", \"omp\", or \"reference\" "));
    }
  } // namespace detail
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
    SolverBase(std::shared_ptr<const gko::Executor> exec,
               SolverControl &                      solver_control);

    /**
     * Destructor.
     */
    virtual ~SolverBase() = default;

    /**
     * Initialize the matrix and copy over its data to Ginkgo's data structures.
     */
    [[deprecated(
      "Create GinkgoWrapper objects directly and call solve on them")]] void
    initialize(const ::dealii::SparseMatrix<ValueType> &matrix);

    /**
     * Solve the linear system <tt>Ax=b</tt>. Dependent on the information
     * provided by derived classes one of Ginkgo's linear solvers is
     * chosen.
     */
    [[deprecated(
      "Create GinkgoWrapper objects directly and call solve on them")]] void
    apply(::dealii::Vector<ValueType> &      solution,
          const ::dealii::Vector<ValueType> &rhs);

    /**
     * Solve the linear system <tt>Ax=b</tt>. Dependent on the information
     * provided by derived classes one of Ginkgo's linear solvers is
     * chosen.
     */
    [[deprecated(
      "Create GinkgoWrapper objects directly and call solve on them")]] void
    solve(const ::dealii::SparseMatrix<ValueType> &matrix,
          ::dealii::Vector<ValueType> &            solution,
          const ::dealii::Vector<ValueType> &      rhs);

    void
    solve(const AbstractMatrix<ValueType, IndexType> &  matrix,
          Vector<ValueType> &                           solution,
          const Vector<ValueType> &                     rhs,
          const PreconditionBase<ValueType, IndexType> &preconditioner =
            PreconditionIdentity<ValueType, IndexType>());

    /**
     * Access to the object that controls convergence.
     */
    SolverControl &
    control() const;

  protected:
    virtual void
    solve_impl(
      const AbstractMatrix<ValueType, IndexType> &  matrix,
      Vector<ValueType> &                           solution,
      const Vector<ValueType> &                     rhs,
      const PreconditionBase<ValueType, IndexType> &preconditioner) = 0;

    /**
     * Reference to the object that controls convergence of the iterative
     * solvers.
     */
    SolverControl &solver_control;

    /**
     * The execution paradigm in Ginkgo. The choices are between
     * `gko::OmpExecutor`, `gko::CudaExecutor` and `gko::ReferenceExecutor`
     * and more details can be found in Ginkgo's documentation.
     */
    std::shared_ptr<const gko::Executor> executor;

  private:
    /**
     * Ginkgo matrix data structure. First template parameter is for storing the
     * array of the non-zeros of the matrix. The second is for the row pointers
     * and the column indices.
     *
     * @todo Templatize based on Matrix type.
     */
    Csr<ValueType, IndexType> system_matrix;
  };


  namespace detail
  {
    template <template <class, class> class Solver>
    struct solver_traits;
  }


  template <typename ValueType,
            typename IndexType,
            template <class>
            class GinkgoType,
            typename AdditionalDataType>
  class EnableSolverBase : public SolverBase<ValueType, IndexType>
  {
  public:
    using ginkgo_type     = GinkgoType<ValueType>;
    using parameters_type = decltype(std::declval<ginkgo_type>().build());
    using AdditionalData  = AdditionalDataType;

    EnableSolverBase(std::shared_ptr<const gko::Executor> exec,
                     SolverControl &                      solver_control,
                     const AdditionalData &               data = {})
      : SolverBase<ValueType, IndexType>(exec, solver_control)
      , parameters_(ginkgo_type::build())
    {
      data.enhance_parameters(parameters_);
      combined_factory =
        gko::stop::Combined::build()
          .with_criteria(gko::stop::ResidualNormReduction<ValueType>::build()
                           .with_reduction_factor(
                             static_cast<gko::remove_complex<ValueType>>(
                               solver_control.tolerance()))
                           .on(this->executor),
                         gko::stop::Iteration::build()
                           .with_max_iters(solver_control.max_steps())
                           .on(this->executor))
          .on(this->executor);
    }

  protected:
    EnableSolverBase(std::shared_ptr<const gko::Executor>      exec,
                     SolverControl &                           solver_control,
                     const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                     const AdditionalData &                    data = {});

    void
    solve_impl(
      const AbstractMatrix<ValueType, IndexType> &  matrix,
      Vector<ValueType> &                           solution,
      const Vector<ValueType> &                     rhs,
      const PreconditionBase<ValueType, IndexType> &preconditioner) override;

    /**
     * Initialize the Ginkgo logger object with event masks. Refer to
     * <a
     * href="https://github.com/ginkgo-project/ginkgo/blob/develop/include/ginkgo/core/log/logger.hpp">Ginkgo's
     * logging event masks.</a>
     */
    void
    initialize_ginkgo_log();

    parameters_type parameters_;

    /**
     * The Ginkgo generated solver factory object.
     */
    std::shared_ptr<gko::LinOpFactory> solver_gen;

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
  };


  namespace detail
  {
    struct EmptyAdditionalData
    {
      template <typename ParametersType>
      void
      enhance_parameters(ParametersType &) const
      {}
    };

    struct GmresAdditionalData
    {
      /**
       * Maximum number of tmp vectors.
       */
      unsigned int restart_parameter = 30;

      /**
       * Should flexible GMRES be used?
       */
      bool flexible = false;

      template <typename ParametersType>
      void
      enhance_parameters(ParametersType &parameters) const
      {
        parameters.with_krylov_dim(restart_parameter).with_flexible(flexible);
      }
    };
  } // namespace detail


  /**
   * An implementation of the solver interface using the Ginkgo CG solver.
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverCG : public EnableSolverBase<ValueType,
                                           IndexType,
                                           gko::solver::Cg,
                                           detail::EmptyAdditionalData>
  {
    using Base = EnableSolverBase<ValueType,
                                  IndexType,
                                  gko::solver::Cg,
                                  detail::EmptyAdditionalData>;

  public:
    using Base::Base;
    using typename Base::AdditionalData;

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
    [[deprecated("Use constructor with ginkgo executor")]] SolverCG(
      SolverControl &       solver_control,
      const std::string &   exec_type,
      const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {}

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
    [[deprecated("Provide preconditioner to solve method directly")]] SolverCG(
      SolverControl &                           solver_control,
      const std::string &                       exec_type,
      const std::shared_ptr<gko::LinOpFactory> &preconditioner,
      const AdditionalData &                    data = AdditionalData())
      : Base(detail::exec_from_name(exec_type),
             solver_control,
             preconditioner,
             data)
    {}
  };


  /**
   * An implementation of the solver interface using the Ginkgo Bicgstab solver.
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverBicgstab : public EnableSolverBase<ValueType,
                                                 IndexType,
                                                 gko::solver::Bicgstab,
                                                 detail::EmptyAdditionalData>
  {
    using Base = EnableSolverBase<ValueType,
                                  IndexType,
                                  gko::solver::Bicgstab,
                                  detail::EmptyAdditionalData>;

  public:
    using Base::Base;
    using typename Base::AdditionalData;

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
    SolverBicgstab(SolverControl &       solver_control,
                   const std::string &   exec_type,
                   const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {}

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
    SolverBicgstab(SolverControl &                           solver_control,
                   const std::string &                       exec_type,
                   const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                   const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type),
             solver_control,
             preconditioner,
             data)
    {}
  };

  /**
   * An implementation of the solver interface using the Ginkgo CGS solver.
   *
   * CGS or the conjugate gradient square method is an iterative type Krylov
   * subspace method which is suitable for general systems.
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverCGS : public EnableSolverBase<ValueType,
                                            IndexType,
                                            gko::solver::Cgs,
                                            detail::EmptyAdditionalData>
  {
    using Base = EnableSolverBase<ValueType,
                                  IndexType,
                                  gko::solver::Cgs,
                                  detail::EmptyAdditionalData>;

  public:
    using Base::Base;
    using typename Base::AdditionalData;

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
    SolverCGS(SolverControl &       solver_control,
              const std::string &   exec_type,
              const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {}

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
    SolverCGS(SolverControl &                           solver_control,
              const std::string &                       exec_type,
              const std::shared_ptr<gko::LinOpFactory> &preconditioner,
              const AdditionalData &                    data = AdditionalData())
      : Base(detail::exec_from_name(exec_type),
             solver_control,
             preconditioner,
             data)
    {}
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
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverFCG : public EnableSolverBase<ValueType,
                                            IndexType,
                                            gko::solver::Fcg,
                                            detail::EmptyAdditionalData>
  {
    using Base = EnableSolverBase<ValueType,
                                  IndexType,
                                  gko::solver::Fcg,
                                  detail::EmptyAdditionalData>;

  public:
    using Base::Base;
    using typename Base::AdditionalData;

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
    SolverFCG(SolverControl &       solver_control,
              const std::string &   exec_type,
              const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {}

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
    SolverFCG(SolverControl &                           solver_control,
              const std::string &                       exec_type,
              const std::shared_ptr<gko::LinOpFactory> &preconditioner,
              const AdditionalData &                    data = AdditionalData())
      : Base(detail::exec_from_name(exec_type),
             solver_control,
             preconditioner,
             data)
    {}
  };

  /**
   * An implementation of the solver interface using the Ginkgo GMRES solver.
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverGMRES : public EnableSolverBase<ValueType,
                                              IndexType,
                                              gko::solver::Gmres,
                                              detail::GmresAdditionalData>
  {
    using Base = EnableSolverBase<ValueType,
                                  IndexType,
                                  gko::solver::Gmres,
                                  detail::GmresAdditionalData>;

  public:
    using Base::Base;
    using typename Base::AdditionalData;

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
    SolverGMRES(SolverControl &       solver_control,
                const std::string &   exec_type,
                const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {}

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
    SolverGMRES(SolverControl &                           solver_control,
                const std::string &                       exec_type,
                const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type),
             solver_control,
             preconditioner,
             data)
    {}
  };

  /**
   * An implementation of the solver interface using the Ginkgo IR solver.
   *
   * Iterative refinement (IR) is an iterative method that uses another coarse
   * method to approximate the error of the current solution via the current
   * residual.
   *
   * @ingroup GinkgoWrappers
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverIR : public EnableSolverBase<ValueType,
                                           IndexType,
                                           gko::solver::Ir,
                                           detail::EmptyAdditionalData>
  {
    using Base = EnableSolverBase<ValueType,
                                  IndexType,
                                  gko::solver::Ir,
                                  detail::EmptyAdditionalData>;

  public:
    using Base::Base;
    using typename Base::AdditionalData;

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
    SolverIR(SolverControl &       solver_control,
             const std::string &   exec_type,
             const AdditionalData &data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {}

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
    SolverIR(SolverControl &                           solver_control,
             const std::string &                       exec_type,
             const std::shared_ptr<gko::LinOpFactory> &inner_solver,
             const AdditionalData &                    data = AdditionalData())
      : Base(detail::exec_from_name(exec_type), solver_control, data)
    {
      this->solver_gen = this->parameters_.with_criteria(this->combined_factory)
                           .with_solver(inner_solver)
                           .on(this->executor);
    }
  };


} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_GINKGO

#endif
