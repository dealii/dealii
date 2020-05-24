// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


#include <deal.II/base/logstream.h>

#include <deal.II/lac/ginkgo_solver.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/lac/exceptions.h>

#  include <cmath>
#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  template <typename ValueType, typename IndexType>
  SolverBase<ValueType, IndexType>::SolverBase(SolverControl &solver_control,
                                               const std::string &exec_type)
    : solver_control(solver_control)
    , exec_type(exec_type)
  {
    if (exec_type == "reference")
      {
        executor = gko::ReferenceExecutor::create();
      }
    else if (exec_type == "omp")
      {
        executor = gko::OmpExecutor::create();
      }
    else if (exec_type == "cuda" && gko::CudaExecutor::get_num_devices() > 0)
      {
        executor = gko::CudaExecutor::create(0, gko::OmpExecutor::create());
      }
    else
      {
        Assert(
          false,
          ExcMessage(
            " exec_type needs to be one of the three strings: \"reference\", \"cuda\" or \"omp\" "));
      }
    using ResidualCriterionFactory = gko::stop::ResidualNormReduction<>;
    residual_criterion             = ResidualCriterionFactory::build()
                           .with_reduction_factor(solver_control.tolerance())
                           .on(executor);

    combined_factory =
      gko::stop::Combined::build()
        .with_criteria(residual_criterion,
                       gko::stop::Iteration::build()
                         .with_max_iters(solver_control.max_steps())
                         .on(executor))
        .on(executor);
  }



  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::initialize_ginkgo_log()
  {
    // Add the logger object. See the different masks available in Ginkgo's
    // documentation
    convergence_logger = gko::log::Convergence<>::create(
      executor, gko::log::Logger::criterion_check_completed_mask);
  }



  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::apply(Vector<ValueType> &      solution,
                                          const Vector<ValueType> &rhs)
  {
    // some shortcuts.
    using val_array = gko::Array<ValueType>;
    using vec       = gko::matrix::Dense<ValueType>;

    Assert(system_matrix, ExcNotInitialized());
    Assert(executor, ExcNotInitialized());
    Assert(rhs.size() == solution.size(),
           ExcDimensionMismatch(rhs.size(), solution.size()));

    // Generate the solver from the solver using the system matrix.
    auto solver = solver_gen->generate(system_matrix);

    // Create the rhs vector in Ginkgo's format.
    std::vector<ValueType> f(rhs.size());
    std::copy(rhs.begin(), rhs.begin() + rhs.size(), f.begin());
    auto b =
      vec::create(executor,
                  gko::dim<2>(rhs.size(), 1),
                  val_array::view(executor->get_master(), rhs.size(), f.data()),
                  1);

    // Create the solution vector in Ginkgo's format.
    std::vector<ValueType> u(solution.size());
    std::copy(solution.begin(), solution.begin() + solution.size(), u.begin());
    auto x = vec::create(executor,
                         gko::dim<2>(solution.size(), 1),
                         val_array::view(executor->get_master(),
                                         solution.size(),
                                         u.data()),
                         1);

    // Create the logger object to log some data from the solvers to confirm
    // convergence.
    initialize_ginkgo_log();

    Assert(convergence_logger, ExcNotInitialized());
    // Add the convergence logger object to the combined factory to retrieve the
    // solver and other data
    combined_factory->add_logger(convergence_logger);

    // Finally, apply the solver to b and get the solution in x.
    solver->apply(gko::lend(b), gko::lend(x));

    // The convergence_logger object contains the residual vector after the
    // solver has returned. use this vector to compute the residual norm of the
    // solution. Get the residual norm from the logger. As the convergence
    // logger returns a `linop`, it is necessary to convert it to a Dense
    // matrix. Additionally, if the logger is logging on the gpu, it is
    // necessary to copy the data to the host and hence the
    // `residual_norm_d_master`
    auto residual_norm = convergence_logger->get_residual_norm();
    auto residual_norm_d =
      gko::as<gko::matrix::Dense<ValueType>>(residual_norm);
    auto residual_norm_d_master =
      gko::matrix::Dense<ValueType>::create(executor->get_master(),
                                            gko::dim<2>{1, 1});
    residual_norm_d_master->copy_from(residual_norm_d);

    // Get the number of iterations taken to converge to the solution.
    auto num_iteration = convergence_logger->get_num_iterations();

    // Ginkgo works with a relative residual norm through its
    // ResidualNormReduction criterion. Therefore, to get the normalized
    // residual, we divide by the norm of the rhs.
    auto b_norm = gko::matrix::Dense<ValueType>::create(executor->get_master(),
                                                        gko::dim<2>{1, 1});
    if (executor != executor->get_master())
      {
        auto b_master = vec::create(executor->get_master(),
                                    gko::dim<2>(rhs.size(), 1),
                                    val_array::view(executor->get_master(),
                                                    rhs.size(),
                                                    f.data()),
                                    1);
        b_master->compute_norm2(b_norm.get());
      }
    else
      {
        b->compute_norm2(b_norm.get());
      }

    Assert(b_norm.get()->at(0, 0) != 0.0, ExcDivideByZero());
    // Pass the number of iterations and residual norm to the solver_control
    // object. As both `residual_norm_d_master` and `b_norm` are seen as Dense
    // matrices, we use the `at` function to get the first value here. In case
    // of multiple right hand sides, this will need to be modified.
    const SolverControl::State state =
      solver_control.check(num_iteration,
                           residual_norm_d_master->at(0, 0) / b_norm->at(0, 0));

    // in case of failure: throw exception
    if (state != SolverControl::success)
      AssertThrow(false,
                  SolverControl::NoConvergence(solver_control.last_step(),
                                               solver_control.last_value()));

    // Check if the solution is on a CUDA device, if so, copy it over to the
    // host.
    if (executor != executor->get_master())
      {
        auto x_master = vec::create(executor->get_master(),
                                    gko::dim<2>(solution.size(), 1),
                                    val_array::view(executor,
                                                    solution.size(),
                                                    x->get_values()),
                                    1);
        x.reset(x_master.release());
      }
    // Finally copy over the solution vector to deal.II's solution vector.
    std::copy(x->get_values(),
              x->get_values() + solution.size(),
              solution.begin());
  }



  template <typename ValueType, typename IndexType>
  SolverControl &
  SolverBase<ValueType, IndexType>::control() const
  {
    return solver_control;
  }



  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::initialize(
    const SparseMatrix<ValueType> &matrix)
  {
    // Needs to be a square matrix
    Assert(matrix.m() == matrix.n(), ExcNotQuadratic());

    using size_type   = dealii::types::global_dof_index;
    const size_type N = matrix.m();

    using mtx = gko::matrix::Csr<ValueType, IndexType>;
    std::shared_ptr<mtx> system_matrix_compute;
    system_matrix_compute   = mtx::create(executor->get_master(),
                                        gko::dim<2>(N),
                                        matrix.n_nonzero_elements());
    ValueType *mat_values   = system_matrix_compute->get_values();
    IndexType *mat_row_ptrs = system_matrix_compute->get_row_ptrs();
    IndexType *mat_col_idxs = system_matrix_compute->get_col_idxs();

    // Copy over the data from the matrix to the data structures Ginkgo needs.
    //
    // Final note: if the matrix has entries in the sparsity pattern that are
    // actually occupied by entries that have a zero numerical value, then we
    // keep them anyway. people are supposed to provide accurate sparsity
    // patterns.

    // first fill row lengths array
    mat_row_ptrs[0] = 0;
    for (size_type row = 1; row <= N; ++row)
      mat_row_ptrs[row] =
        mat_row_ptrs[row - 1] + matrix.get_row_length(row - 1);

    // Copy over matrix elements. note that for sparse matrices,
    // iterators are sorted so that they traverse each row from start to end
    // before moving on to the next row. however, this isn't true for block
    // matrices, so we have to do a bit of book keeping
    {
      // Have an array that for each row points to the first entry not yet
      // written to
      std::vector<IndexType> row_pointers(N + 1);
      std::copy(system_matrix_compute->get_row_ptrs(),
                system_matrix_compute->get_row_ptrs() + N + 1,
                row_pointers.begin());

      // Loop over the elements of the matrix row by row, as suggested in the
      // documentation of the sparse matrix iterator class
      for (size_type row = 0; row < N; ++row)
        {
          for (typename SparseMatrix<ValueType>::const_iterator p =
                 matrix.begin(row);
               p != matrix.end(row);
               ++p)
            {
              // Write entry into the first free one for this row
              mat_col_idxs[row_pointers[row]] = p->column();
              mat_values[row_pointers[row]]   = p->value();

              // Then move pointer ahead
              ++row_pointers[row];
            }
        }

      // At the end, we should have written all rows completely
      for (size_type i = 0; i < N - 1; ++i)
        Assert(row_pointers[i] == mat_row_ptrs[i + 1], ExcInternalError());
    }
    system_matrix =
      mtx::create(executor, gko::dim<2>(N), matrix.n_nonzero_elements());
    system_matrix->copy_from(system_matrix_compute.get());
  }



  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::solve(const SparseMatrix<ValueType> &matrix,
                                          Vector<ValueType> &      solution,
                                          const Vector<ValueType> &rhs)
  {
    initialize(matrix);
    apply(solution, rhs);
  }



  /* ---------------------- SolverCG ------------------------ */
  template <typename ValueType, typename IndexType>
  SolverCG<ValueType, IndexType>::SolverCG(SolverControl &       solver_control,
                                           const std::string &   exec_type,
                                           const AdditionalData &data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using cg = gko::solver::Cg<ValueType>;
    this->solver_gen =
      cg::build().with_criteria(this->combined_factory).on(this->executor);
  }



  template <typename ValueType, typename IndexType>
  SolverCG<ValueType, IndexType>::SolverCG(
    SolverControl &                           solver_control,
    const std::string &                       exec_type,
    const std::shared_ptr<gko::LinOpFactory> &preconditioner,
    const AdditionalData &                    data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using cg         = gko::solver::Cg<ValueType>;
    this->solver_gen = cg::build()
                         .with_criteria(this->combined_factory)
                         .with_preconditioner(preconditioner)
                         .on(this->executor);
  }



  /* ---------------------- SolverBicgstab ------------------------ */
  template <typename ValueType, typename IndexType>
  SolverBicgstab<ValueType, IndexType>::SolverBicgstab(
    SolverControl &       solver_control,
    const std::string &   exec_type,
    const AdditionalData &data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using bicgstab   = gko::solver::Bicgstab<ValueType>;
    this->solver_gen = bicgstab::build()
                         .with_criteria(this->combined_factory)
                         .on(this->executor);
  }



  template <typename ValueType, typename IndexType>
  SolverBicgstab<ValueType, IndexType>::SolverBicgstab(
    SolverControl &                           solver_control,
    const std::string &                       exec_type,
    const std::shared_ptr<gko::LinOpFactory> &preconditioner,
    const AdditionalData &                    data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using bicgstab   = gko::solver::Bicgstab<ValueType>;
    this->solver_gen = bicgstab::build()
                         .with_criteria(this->combined_factory)
                         .with_preconditioner(preconditioner)
                         .on(this->executor);
  }



  /* ---------------------- SolverCGS ------------------------ */
  template <typename ValueType, typename IndexType>
  SolverCGS<ValueType, IndexType>::SolverCGS(SolverControl &    solver_control,
                                             const std::string &exec_type,
                                             const AdditionalData &data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using cgs = gko::solver::Cgs<ValueType>;
    this->solver_gen =
      cgs::build().with_criteria(this->combined_factory).on(this->executor);
  }



  template <typename ValueType, typename IndexType>
  SolverCGS<ValueType, IndexType>::SolverCGS(
    SolverControl &                           solver_control,
    const std::string &                       exec_type,
    const std::shared_ptr<gko::LinOpFactory> &preconditioner,
    const AdditionalData &                    data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using cgs        = gko::solver::Cgs<ValueType>;
    this->solver_gen = cgs::build()
                         .with_criteria(this->combined_factory)
                         .with_preconditioner(preconditioner)
                         .on(this->executor);
  }



  /* ---------------------- SolverFCG ------------------------ */
  template <typename ValueType, typename IndexType>
  SolverFCG<ValueType, IndexType>::SolverFCG(SolverControl &    solver_control,
                                             const std::string &exec_type,
                                             const AdditionalData &data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using fcg = gko::solver::Fcg<ValueType>;
    this->solver_gen =
      fcg::build().with_criteria(this->combined_factory).on(this->executor);
  }



  template <typename ValueType, typename IndexType>
  SolverFCG<ValueType, IndexType>::SolverFCG(
    SolverControl &                           solver_control,
    const std::string &                       exec_type,
    const std::shared_ptr<gko::LinOpFactory> &preconditioner,
    const AdditionalData &                    data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using fcg        = gko::solver::Fcg<ValueType>;
    this->solver_gen = fcg::build()
                         .with_criteria(this->combined_factory)
                         .with_preconditioner(preconditioner)
                         .on(this->executor);
  }



  /* ---------------------- SolverGMRES ------------------------ */
  template <typename ValueType, typename IndexType>
  SolverGMRES<ValueType, IndexType>::AdditionalData::AdditionalData(
    const unsigned int restart_parameter)
    : restart_parameter(restart_parameter)
  {}



  template <typename ValueType, typename IndexType>
  SolverGMRES<ValueType, IndexType>::SolverGMRES(SolverControl &solver_control,
                                                 const std::string &exec_type,
                                                 const AdditionalData &data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using gmres      = gko::solver::Gmres<ValueType>;
    this->solver_gen = gmres::build()
                         .with_krylov_dim(additional_data.restart_parameter)
                         .with_criteria(this->combined_factory)
                         .on(this->executor);
  }



  template <typename ValueType, typename IndexType>
  SolverGMRES<ValueType, IndexType>::SolverGMRES(
    SolverControl &                           solver_control,
    const std::string &                       exec_type,
    const std::shared_ptr<gko::LinOpFactory> &preconditioner,
    const AdditionalData &                    data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using gmres      = gko::solver::Gmres<ValueType>;
    this->solver_gen = gmres::build()
                         .with_krylov_dim(additional_data.restart_parameter)
                         .with_criteria(this->combined_factory)
                         .with_preconditioner(preconditioner)
                         .on(this->executor);
  }



  /* ---------------------- SolverIR ------------------------ */
  template <typename ValueType, typename IndexType>
  SolverIR<ValueType, IndexType>::SolverIR(SolverControl &       solver_control,
                                           const std::string &   exec_type,
                                           const AdditionalData &data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using ir = gko::solver::Ir<ValueType>;
    this->solver_gen =
      ir::build().with_criteria(this->combined_factory).on(this->executor);
  }



  template <typename ValueType, typename IndexType>
  SolverIR<ValueType, IndexType>::SolverIR(
    SolverControl &                           solver_control,
    const std::string &                       exec_type,
    const std::shared_ptr<gko::LinOpFactory> &inner_solver,
    const AdditionalData &                    data)
    : SolverBase<ValueType, IndexType>(solver_control, exec_type)
    , additional_data(data)
  {
    using ir         = gko::solver::Ir<ValueType>;
    this->solver_gen = ir::build()
                         .with_criteria(this->combined_factory)
                         .with_solver(inner_solver)
                         .on(this->executor);
  }



  // Explicit instantiations in GinkgoWrappers
#  define DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(_macro) \
    template _macro(float, int32_t);                               \
    template _macro(double, int32_t);                              \
    template _macro(float, int64_t);                               \
    template _macro(double, int64_t);

#  define DECLARE_SOLVER_BASE(ValueType, IndexType) \
    class SolverBase<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_BASE)
#  undef DECLARE_SOLVER_BASE

#  define DECLARE_SOLVER_CG(ValueType, IndexType) \
    class SolverCG<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_CG)
#  undef DECLARE_SOLVER_CG

#  define DECLARE_SOLVER_Bicgstab(ValueType, IndexType) \
    class SolverBicgstab<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_Bicgstab)
#  undef DECLARE_SOLVER_Bicgstab

#  define DECLARE_SOLVER_CGS(ValueType, IndexType) \
    class SolverCGS<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_CGS)
#  undef DECLARE_SOLVER_CGS

#  define DECLARE_SOLVER_FCG(ValueType, IndexType) \
    class SolverFCG<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_FCG)
#  undef DECLARE_SOLVER_FCG

#  define DECLARE_SOLVER_GMRES(ValueType, IndexType) \
    class SolverGMRES<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_GMRES)
#  undef DECLARE_SOLVER_GMRES

#  define DECLARE_SOLVER_IR(ValueType, IndexType) \
    class SolverIR<ValueType, IndexType>
  DEALII_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_SOLVER_IR)
#  undef DECLARE_SOLVER_IR

} // namespace GinkgoWrappers


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_GINKGO
