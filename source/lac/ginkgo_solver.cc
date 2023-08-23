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
    : SolverBase(detail::exec_from_name(exec_type), solver_control)
  {}

  template <typename ValueType, typename IndexType>
  SolverBase<ValueType, IndexType>::SolverBase(
    std::shared_ptr<const gko::Executor> exec,
    SolverControl &                      solver_control)
    : solver_control(solver_control)
    , executor(std::move(exec))
    , system_matrix(executor, 0, 0)
  {}


  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::apply(
    ::dealii::Vector<ValueType> &      solution,
    const ::dealii::Vector<ValueType> &rhs)
  {
    auto view_solution = Vector<ValueType>::create_view(executor, solution);
    auto view_rhs      = Vector<ValueType>::create_view(executor, rhs);

    solve(system_matrix, *view_solution, *view_rhs);
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
    const ::dealii::SparseMatrix<ValueType> &matrix)
  {
    system_matrix = Csr<ValueType, IndexType>{executor, matrix};
  }


  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::solve(
    const ::dealii::SparseMatrix<ValueType> &matrix,
    ::dealii::Vector<ValueType> &            solution,
    const ::dealii::Vector<ValueType> &      rhs)
  {
    initialize(matrix);
    apply(solution, rhs);
  }


  template <typename ValueType, typename IndexType>
  void
  SolverBase<ValueType, IndexType>::solve(
    const AbstractMatrix<ValueType, IndexType> &  matrix,
    Vector<ValueType> &                           solution,
    const Vector<ValueType> &                     rhs,
    const PreconditionBase<ValueType, IndexType> &preconditioner)
  {
    solve_impl(matrix, solution, rhs, preconditioner);
  }


  template <typename ValueType,
            typename IndexType,
            template <class>
            class GinkgoType,
            typename AdditionalDataType>
  void
  EnableSolverBase<ValueType, IndexType, GinkgoType, AdditionalDataType>::
    initialize_ginkgo_log()
  {
    // Add the logger object. See the different masks available in Ginkgo's
    // documentation
    convergence_logger = gko::log::Convergence<>::create(
      this->executor, gko::log::Logger::criterion_check_completed_mask);
  }


  template <typename T, typename = void>
  struct has_with_preconditioner : std::false_type
  {};

  template <typename T>
  struct has_with_preconditioner<
    T,
    boost::void_t<decltype(std::declval<T>().with_generated_preconditioner(
      std::declval<std::shared_ptr<const gko::LinOp>>()))>> : std::true_type
  {};


  template <typename T, typename = void>
  struct has_with_solver : std::false_type
  {};

  template <typename T>
  struct has_with_solver<
    T,
    boost::void_t<decltype(std::declval<T>().with_generated_solver(
      std::declval<std::shared_ptr<const gko::LinOp>>()))>> : std::true_type
  {};


  template <typename ValueType, typename IndexType, typename ParametersType>
  std::enable_if_t<has_with_preconditioner<ParametersType>::value,
                   ParametersType &>
  add_preconditioner(
    ParametersType &                              parameters,
    const PreconditionBase<ValueType, IndexType> &preconditioner)
  {
    return parameters.with_generated_preconditioner(
      preconditioner.get_shared_gko_object());
  }

  template <typename ValueType, typename IndexType, typename ParametersType>
  std::enable_if_t<has_with_solver<ParametersType>::value, ParametersType &>
  add_preconditioner(
    ParametersType &                              parameters,
    const PreconditionBase<ValueType, IndexType> &preconditioner)
  {
    return parameters.with_generated_solver(
      preconditioner.get_shared_gko_object());
  }


  template <typename ValueType,
            typename IndexType,
            template <class>
            class GinkgoType,
            typename AdditionalDataType>
  EnableSolverBase<ValueType, IndexType, GinkgoType, AdditionalDataType>::
    EnableSolverBase(std::shared_ptr<const gko::Executor>      exec,
                     SolverControl &                           solver_control,
                     const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                     const AdditionalData &                    data)
    : EnableSolverBase(exec, solver_control, data)
  {
    solver_gen = ginkgo_type::build()
                   .with_criteria(this->combined_factory)
                   .with_preconditioner(preconditioner)
                   .on(this->executor);
  }


  template <typename ValueType,
            typename IndexType,
            template <class>
            class GinkgoType,
            typename AdditionalDataType>
  void
  EnableSolverBase<ValueType, IndexType, GinkgoType, AdditionalDataType>::
    solve_impl(const AbstractMatrix<ValueType, IndexType> &  matrix,
               Vector<ValueType> &                           solution,
               const Vector<ValueType> &                     rhs,
               const PreconditionBase<ValueType, IndexType> &preconditioner)
  {
    Assert(this->executor, ExcNotInitialized());
    Assert(rhs.size() == solution.size(),
           ExcDimensionMismatch(rhs.size(), solution.size()));

    // Generate the solver from the solver using the system matrix.
    auto solver =
      solver_gen ?
        solver_gen->generate(matrix.get_shared_gko_object()) :
        add_preconditioner(parameters_.with_criteria(combined_factory),
                           preconditioner)
          .on(this->executor)
          ->generate(matrix.get_shared_gko_object());

    // Create the logger object to log some data from the solvers to confirm
    // convergence.
    this->initialize_ginkgo_log();

    Assert(this->convergence_logger, ExcNotInitialized());
    // Add the convergence logger object to the combined factory to retrieve the
    // solver and other data
    this->combined_factory->add_logger(convergence_logger);

    // Finally, apply the solver to b and get the solution in x.
    solver->apply(rhs.get_gko_object(), solution.get_gko_object());

    // The convergence_logger object contains the residual vector after the
    // solver has returned. use this vector to compute the residual norm of the
    // solution. Get the residual norm from the logger. As the convergence
    // logger returns a `linop`, it is necessary to convert it to a Dense
    // matrix. Additionally, if the logger is logging on the gpu, it is
    // necessary to copy the data to the host and hence the
    // `residual_norm_host`
    auto residual_norm = gko::as<typename Vector<ValueType>::ginkgo_type>(
      convergence_logger->get_residual_norm());
    auto residual_norm_host =
      gko::make_temporary_clone(this->executor->get_master(), residual_norm);

    // Get the number of iterations taken to converge to the solution.
    auto num_iteration = convergence_logger->get_num_iterations();

    // Ginkgo works with a relative residual norm through its
    // ResidualNormReduction criterion. Therefore, to get the normalized
    // residual, we divide by the norm of the rhs.
    auto rhs_norm = rhs.l2_norm();

    Assert(rhs_norm != 0.0, ExcDivideByZero());
    // Pass the number of iterations and residual norm to the solver_control
    // object. As both `residual_norm_host` and `rhs_norm` are seen as Dense
    // matrices, we use the `at` function to get the first value here. In case
    // of multiple right hand sides, this will need to be modified.
    const SolverControl::State state =
      this->solver_control.check(num_iteration,
                                 residual_norm_host->at(0, 0) / rhs_norm);

    // in case of failure: throw exception
    if (state != SolverControl::success)
      AssertThrow(
        false,
        SolverControl::NoConvergence(this->solver_control.last_step(),
                                     this->solver_control.last_value()));
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
