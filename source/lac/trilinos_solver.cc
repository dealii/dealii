// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/trilinos_solver.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/conditional_ostream.h>

#  include <deal.II/lac/trilinos_precondition.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#  include <AztecOO_StatusTest.h>
#  include <AztecOO_StatusTestCombo.h>
#  include <AztecOO_StatusTestMaxIters.h>
#  include <AztecOO_StatusTestResNorm.h>
#  include <AztecOO_StatusType.h>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <cmath>
#  include <limits>
#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  SolverBase::AdditionalData::AdditionalData(
    const bool         output_solver_details,
    const unsigned int gmres_restart_parameter)
    : output_solver_details(output_solver_details)
    , gmres_restart_parameter(gmres_restart_parameter)
  {}



  SolverBase::SolverBase(SolverControl &cn, const AdditionalData &data)
    : solver_name(gmres)
    , solver_control(cn)
    , additional_data(data)
  {}



  SolverBase::SolverBase(const SolverBase::SolverName solver_name,
                         SolverControl               &cn,
                         const AdditionalData        &data)
    : solver_name(solver_name)
    , solver_control(cn)
    , additional_data(data)
  {}



  SolverControl &
  SolverBase::control() const
  {
    return solver_control;
  }



  void
  SolverBase::solve(const SparseMatrix     &A,
                    MPI::Vector            &x,
                    const MPI::Vector      &b,
                    const PreconditionBase &preconditioner)
  {
    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()),
      &x.trilinos_vector(),
      const_cast<Epetra_MultiVector *>(&b.trilinos_vector()));

    do_solve(preconditioner);
  }



  // Note: "A" is set as a constant reference so that all patterns for ::solve
  //       can be used by the inverse_operator of LinearOperator
  void
  SolverBase::solve(const Epetra_Operator  &A,
                    MPI::Vector            &x,
                    const MPI::Vector      &b,
                    const PreconditionBase &preconditioner)
  {
    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem =
      std::make_unique<Epetra_LinearProblem>(const_cast<Epetra_Operator *>(&A),
                                             &x.trilinos_vector(),
                                             const_cast<Epetra_MultiVector *>(
                                               &b.trilinos_vector()));

    do_solve(preconditioner);
  }



  // Note: "A" is set as a constant reference so that all patterns for ::solve
  //       can be used by the inverse_operator of LinearOperator
  void
  SolverBase::solve(const Epetra_Operator &A,
                    MPI::Vector           &x,
                    const MPI::Vector     &b,
                    const Epetra_Operator &preconditioner)
  {
    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem =
      std::make_unique<Epetra_LinearProblem>(const_cast<Epetra_Operator *>(&A),
                                             &x.trilinos_vector(),
                                             const_cast<Epetra_MultiVector *>(
                                               &b.trilinos_vector()));

    do_solve(preconditioner);
  }



  // Note: "A" is set as a constant reference so that all patterns for ::solve
  //       can be used by the inverse_operator of LinearOperator
  void
  SolverBase::solve(const Epetra_Operator    &A,
                    Epetra_MultiVector       &x,
                    const Epetra_MultiVector &b,
                    const PreconditionBase   &preconditioner)
  {
    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem =
      std::make_unique<Epetra_LinearProblem>(const_cast<Epetra_Operator *>(&A),
                                             &x,
                                             const_cast<Epetra_MultiVector *>(
                                               &b));

    do_solve(preconditioner);
  }



  // Note: "A" is set as a constant reference so that all patterns for ::solve
  //       can be used by the inverse_operator of LinearOperator
  void
  SolverBase::solve(const Epetra_Operator    &A,
                    Epetra_MultiVector       &x,
                    const Epetra_MultiVector &b,
                    const Epetra_Operator    &preconditioner)
  {
    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem =
      std::make_unique<Epetra_LinearProblem>(const_cast<Epetra_Operator *>(&A),
                                             &x,
                                             const_cast<Epetra_MultiVector *>(
                                               &b));

    do_solve(preconditioner);
  }



  void
  SolverBase::solve(const SparseMatrix           &A,
                    dealii::Vector<double>       &x,
                    const dealii::Vector<double> &b,
                    const PreconditionBase       &preconditioner)
  {
    // In case we call the solver with deal.II vectors, we create views of the
    // vectors in Epetra format.
    Assert(x.size() == A.n(), ExcDimensionMismatch(x.size(), A.n()));
    Assert(b.size() == A.m(), ExcDimensionMismatch(b.size(), A.m()));
    Assert(A.local_range().second == A.m(),
           ExcMessage("Can only work in serial when using deal.II vectors."));
    Assert(A.trilinos_matrix().Filled(),
           ExcMessage("Matrix is not compressed. Call compress() method."));

    Epetra_Vector ep_x(View, A.trilinos_matrix().DomainMap(), x.begin());
    Epetra_Vector ep_b(View,
                       A.trilinos_matrix().RangeMap(),
                       const_cast<double *>(b.begin()));

    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()), &ep_x, &ep_b);

    do_solve(preconditioner);
  }



  void
  SolverBase::solve(Epetra_Operator              &A,
                    dealii::Vector<double>       &x,
                    const dealii::Vector<double> &b,
                    const PreconditionBase       &preconditioner)
  {
    Epetra_Vector ep_x(View, A.OperatorDomainMap(), x.begin());
    Epetra_Vector ep_b(View,
                       A.OperatorRangeMap(),
                       const_cast<double *>(b.begin()));

    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(&A, &ep_x, &ep_b);

    do_solve(preconditioner);
  }



  void
  SolverBase::solve(const SparseMatrix                                       &A,
                    dealii::LinearAlgebra::distributed::Vector<double>       &x,
                    const dealii::LinearAlgebra::distributed::Vector<double> &b,
                    const PreconditionBase &preconditioner)
  {
    // In case we call the solver with deal.II vectors, we create views of the
    // vectors in Epetra format.
    AssertDimension(x.locally_owned_size(),
                    A.trilinos_matrix().DomainMap().NumMyElements());
    AssertDimension(b.locally_owned_size(),
                    A.trilinos_matrix().RangeMap().NumMyElements());

    Epetra_Vector ep_x(View, A.trilinos_matrix().DomainMap(), x.begin());
    Epetra_Vector ep_b(View,
                       A.trilinos_matrix().RangeMap(),
                       const_cast<double *>(b.begin()));

    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()), &ep_x, &ep_b);

    do_solve(preconditioner);
  }



  void
  SolverBase::solve(Epetra_Operator                                          &A,
                    dealii::LinearAlgebra::distributed::Vector<double>       &x,
                    const dealii::LinearAlgebra::distributed::Vector<double> &b,
                    const PreconditionBase &preconditioner)
  {
    AssertDimension(x.locally_owned_size(),
                    A.OperatorDomainMap().NumMyElements());
    AssertDimension(b.locally_owned_size(),
                    A.OperatorRangeMap().NumMyElements());

    Epetra_Vector ep_x(View, A.OperatorDomainMap(), x.begin());
    Epetra_Vector ep_b(View,
                       A.OperatorRangeMap(),
                       const_cast<double *>(b.begin()));

    // We need an Epetra_LinearProblem object to let the AztecOO solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(&A, &ep_x, &ep_b);

    do_solve(preconditioner);
  }


  namespace internal
  {
    namespace
    {
      double
      compute_residual(const Epetra_MultiVector *const residual_vector)
      {
        Assert(residual_vector->NumVectors() == 1,
               ExcMessage("Residual multivector holds more than one vector"));
        TrilinosScalar res_l2_norm = 0.0;
        const int      ierr        = residual_vector->Norm2(&res_l2_norm);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        return res_l2_norm;
      }

      class TrilinosReductionControl : public AztecOO_StatusTest
      {
      public:
        TrilinosReductionControl(const int                   max_steps,
                                 const double                tolerance,
                                 const double                reduction,
                                 const Epetra_LinearProblem &linear_problem);

        virtual ~TrilinosReductionControl() override = default;

        virtual bool
        ResidualVectorRequired() const override
        {
          return status_test_collection->ResidualVectorRequired();
        }

        virtual AztecOO_StatusType
        CheckStatus(int                 CurrentIter,
                    Epetra_MultiVector *CurrentResVector,
                    double              CurrentResNormEst,
                    bool                SolutionUpdated) override
        {
          // Note: CurrentResNormEst is set to -1.0 if no estimate of the
          // residual value is available
          current_residual =
            (CurrentResNormEst < 0.0 ? compute_residual(CurrentResVector) :
                                       CurrentResNormEst);
          if (CurrentIter == 0)
            initial_residual = current_residual;

          return status_test_collection->CheckStatus(CurrentIter,
                                                     CurrentResVector,
                                                     CurrentResNormEst,
                                                     SolutionUpdated);
        }

        virtual AztecOO_StatusType
        GetStatus() const override
        {
          return status_test_collection->GetStatus();
        }

        virtual std::ostream &
        Print(std::ostream &stream, int indent = 0) const override
        {
          return status_test_collection->Print(stream, indent);
        }

        double
        get_initial_residual() const
        {
          return initial_residual;
        }

        double
        get_current_residual() const
        {
          return current_residual;
        }

      private:
        double                                      initial_residual;
        double                                      current_residual;
        std::unique_ptr<AztecOO_StatusTestCombo>    status_test_collection;
        std::unique_ptr<AztecOO_StatusTestMaxIters> status_test_max_steps;
        std::unique_ptr<AztecOO_StatusTestResNorm>  status_test_abs_tol;
        std::unique_ptr<AztecOO_StatusTestResNorm>  status_test_rel_tol;
      };


      TrilinosReductionControl::TrilinosReductionControl(
        const int                   max_steps,
        const double                tolerance,
        const double                reduction,
        const Epetra_LinearProblem &linear_problem)
        : initial_residual(std::numeric_limits<double>::max())
        , current_residual(std::numeric_limits<double>::max())
        // Consider linear problem converged if any of the collection of
        // criterion are met
        , status_test_collection(std::make_unique<AztecOO_StatusTestCombo>(
            AztecOO_StatusTestCombo::OR))
      {
        // Maximum number of iterations
        Assert(max_steps >= 0, ExcInternalError());
        status_test_max_steps =
          std::make_unique<AztecOO_StatusTestMaxIters>(max_steps);
        status_test_collection->AddStatusTest(*status_test_max_steps);

        Assert(linear_problem.GetRHS()->NumVectors() == 1,
               ExcMessage("RHS multivector holds more than one vector"));

        // Residual norm is below some absolute value
        status_test_abs_tol = std::make_unique<AztecOO_StatusTestResNorm>(
          *linear_problem.GetOperator(),
          *(linear_problem.GetLHS()->operator()(0)),
          *(linear_problem.GetRHS()->operator()(0)),
          tolerance);
        status_test_abs_tol->DefineResForm(AztecOO_StatusTestResNorm::Explicit,
                                           AztecOO_StatusTestResNorm::TwoNorm);
        status_test_abs_tol->DefineScaleForm(
          AztecOO_StatusTestResNorm::None, AztecOO_StatusTestResNorm::TwoNorm);
        status_test_collection->AddStatusTest(*status_test_abs_tol);

        // Residual norm, scaled by some initial value, is below some threshold
        status_test_rel_tol = std::make_unique<AztecOO_StatusTestResNorm>(
          *linear_problem.GetOperator(),
          *(linear_problem.GetLHS()->operator()(0)),
          *(linear_problem.GetRHS()->operator()(0)),
          reduction);
        status_test_rel_tol->DefineResForm(AztecOO_StatusTestResNorm::Explicit,
                                           AztecOO_StatusTestResNorm::TwoNorm);
        status_test_rel_tol->DefineScaleForm(
          AztecOO_StatusTestResNorm::NormOfInitRes,
          AztecOO_StatusTestResNorm::TwoNorm);
        status_test_collection->AddStatusTest(*status_test_rel_tol);
      }

    } // namespace
  }   // namespace internal


  template <typename Preconditioner>
  void
  SolverBase::do_solve(const Preconditioner &preconditioner)
  {
    int ierr;

    // Next we can allocate the AztecOO solver...
    solver.SetProblem(*linear_problem);

    // ... and we can specify the solver to be used.
    switch (solver_name)
      {
        case cg:
          solver.SetAztecOption(AZ_solver, AZ_cg);
          break;
        case cgs:
          solver.SetAztecOption(AZ_solver, AZ_cgs);
          break;
        case gmres:
          solver.SetAztecOption(AZ_solver, AZ_gmres);
          solver.SetAztecOption(AZ_kspace,
                                additional_data.gmres_restart_parameter);
          break;
        case bicgstab:
          solver.SetAztecOption(AZ_solver, AZ_bicgstab);
          break;
        case tfqmr:
          solver.SetAztecOption(AZ_solver, AZ_tfqmr);
          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }

    // Set the preconditioner
    set_preconditioner(solver, preconditioner);

    // ... set some options, ...
    solver.SetAztecOption(AZ_output,
                          additional_data.output_solver_details ? AZ_all :
                                                                  AZ_none);
    solver.SetAztecOption(AZ_conv, AZ_noscaled);

    // By default, the Trilinos solver chooses convergence criterion based on
    // the number of iterations made and an absolute tolerance.
    // This implies that the use of the standard Trilinos convergence test
    // actually coincides with dealii::IterationNumberControl because the
    // solver, unless explicitly told otherwise, will Iterate() until a number
    // of max_steps() are taken or an absolute tolerance() is attained.
    // It is therefore suitable for use with both SolverControl or
    // IterationNumberControl. The final check at the end will determine whether
    // failure to converge to the defined residual norm constitutes failure
    // (SolverControl) or is alright (IterationNumberControl).
    // In the case that the SolverControl wants to perform ReductionControl,
    // then we have to do a little extra something by prescribing a custom
    // status test.
    if (!status_test)
      {
        if (const ReductionControl *const reduction_control =
              dynamic_cast<const ReductionControl *>(&solver_control))
          {
            status_test = std::make_unique<internal::TrilinosReductionControl>(
              reduction_control->max_steps(),
              reduction_control->tolerance(),
              reduction_control->reduction(),
              *linear_problem);
            solver.SetStatusTest(status_test.get());
          }
      }

    // ... and then solve!
    ierr =
      solver.Iterate(solver_control.max_steps(), solver_control.tolerance());

    // report errors in more detail than just by checking whether the return
    // status is zero or greater. the error strings are taken from the
    // implementation of the AztecOO::Iterate function
    switch (ierr)
      {
        case -1:
          AssertThrow(false,
                      ExcMessage("AztecOO::Iterate error code -1: "
                                 "option not implemented"));
          break;
        case -2:
          AssertThrow(false,
                      ExcMessage("AztecOO::Iterate error code -2: "
                                 "numerical breakdown"));
          break;
        case -3:
          AssertThrow(false,
                      ExcMessage("AztecOO::Iterate error code -3: "
                                 "loss of precision"));
          break;
        case -4:
          AssertThrow(false,
                      ExcMessage("AztecOO::Iterate error code -4: "
                                 "GMRES Hessenberg ill-conditioned"));
          break;
        default:
          AssertThrow(ierr >= 0, ExcTrilinosError(ierr));
      }

    // Finally, let the deal.II SolverControl object know what has
    // happened. If the solve succeeded, the status of the solver control will
    // turn into SolverControl::success.
    // If the residual is not computed/stored by the solver, as can happen for
    // certain choices of solver or if a custom status test is set, then the
    // result returned by TrueResidual() is equal to -1. In this case we must
    // compute it ourself.
    if (const internal::TrilinosReductionControl
          *const reduction_control_status =
            dynamic_cast<const internal::TrilinosReductionControl *>(
              status_test.get()))
      {
        Assert(dynamic_cast<const ReductionControl *>(&solver_control),
               ExcInternalError());

        // Check to see if solver converged in one step
        // This can happen if the matrix is diagonal and a non-trivial
        // preconditioner is used.
        if (solver.NumIters() > 0)
          {
            // For ReductionControl, we must first register the initial residual
            // value. This is the basis from which it will determine whether the
            // current residual corresponds to a converged state.
            solver_control.check(
              0, reduction_control_status->get_initial_residual());
            solver_control.check(
              solver.NumIters(),
              reduction_control_status->get_current_residual());
          }
        else
          solver_control.check(
            solver.NumIters(),
            reduction_control_status->get_current_residual());
      }
    else
      {
        Assert(solver.TrueResidual() >= 0.0, ExcInternalError());
        solver_control.check(solver.NumIters(), solver.TrueResidual());
      }

    if (solver_control.last_check() != SolverControl::success)
      AssertThrow(false,
                  SolverControl::NoConvergence(solver_control.last_step(),
                                               solver_control.last_value()));
  }



  template <>
  void
  SolverBase::set_preconditioner(AztecOO                &solver,
                                 const PreconditionBase &preconditioner)
  {
    // Introduce the preconditioner, if the identity preconditioner is used,
    // the precondioner is set to none, ...
    if (preconditioner.preconditioner.strong_count() != 0)
      {
        const int ierr = solver.SetPrecOperator(
          const_cast<Epetra_Operator *>(preconditioner.preconditioner.get()));
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }
    else
      solver.SetAztecOption(AZ_precond, AZ_none);
  }


  template <>
  void
  SolverBase::set_preconditioner(AztecOO               &solver,
                                 const Epetra_Operator &preconditioner)
  {
    const int ierr =
      solver.SetPrecOperator(const_cast<Epetra_Operator *>(&preconditioner));
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }


  /* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cg, cn, data)
  {}


  /* ---------------------- SolverGMRES ------------------------ */

  SolverGMRES::SolverGMRES(SolverControl &cn, const AdditionalData &data)
    : SolverBase(gmres, cn, data)
  {}


  /* ---------------------- SolverBicgstab ------------------------ */

  SolverBicgstab::SolverBicgstab(SolverControl &cn, const AdditionalData &data)
    : SolverBase(bicgstab, cn, data)
  {}


  /* ---------------------- SolverCGS ------------------------ */

  SolverCGS::SolverCGS(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cgs, cn, data)
  {}


  /* ---------------------- SolverTFQMR ------------------------ */

  SolverTFQMR::SolverTFQMR(SolverControl &cn, const AdditionalData &data)
    : SolverBase(tfqmr, cn, data)
  {}



  /* ---------------------- SolverDirect ------------------------ */

  SolverDirect::AdditionalData::AdditionalData(const bool output_solver_details,
                                               const std::string &solver_type)
    : output_solver_details(output_solver_details)
    , solver_type(solver_type)
  {}



  SolverDirect::SolverDirect(const AdditionalData &data)
    : solver_control(solver_control_own)
    , additional_data(data.output_solver_details, data.solver_type)
  {}



  SolverDirect::SolverDirect(SolverControl &cn, const AdditionalData &data)
    : solver_control(cn)
    , additional_data(data.output_solver_details, data.solver_type)
  {}



  SolverControl &
  SolverDirect::control() const
  {
    return solver_control;
  }



  void
  SolverDirect::initialize(const SparseMatrix &A)
  {
    // We need an Epetra_LinearProblem object to let the Amesos solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>();

    // Assign the matrix operator to the Epetra_LinearProblem object
    linear_problem->SetOperator(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()));

    // Fetch return value of Amesos Solver functions
    int ierr;

    // First set whether we want to print the solver information to screen or
    // not.
    ConditionalOStream verbose_cout(std::cout,
                                    additional_data.output_solver_details);

    // Next allocate the Amesos solver, this is done in two steps, first we
    // create a solver Factory and generate with that the concrete Amesos
    // solver, if possible.
    Amesos Factory;

    AssertThrow(Factory.Query(additional_data.solver_type.c_str()),
                ExcMessage(
                  "You tried to select the solver type <" +
                  additional_data.solver_type +
                  "> but this solver is not supported by Trilinos either "
                  "because it does not exist, or because Trilinos was not "
                  "configured for its use."));

    solver.reset(
      Factory.Create(additional_data.solver_type.c_str(), *linear_problem));

    verbose_cout << "Starting symbolic factorization" << std::endl;
    ierr = solver->SymbolicFactorization();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    verbose_cout << "Starting numeric factorization" << std::endl;
    ierr = solver->NumericFactorization();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SolverDirect::initialize(const SparseMatrix &A, const AdditionalData &data)
  {
    this->additional_data = data;

    this->initialize(A);
  }


  void
  SolverDirect::solve(MPI::Vector &x, const MPI::Vector &b)
  {
    this->vmult(x, b);
  }



  void
  SolverDirect::solve(
    dealii::LinearAlgebra::distributed::Vector<double>       &x,
    const dealii::LinearAlgebra::distributed::Vector<double> &b)
  {
    this->vmult(x, b);
  }


  void
  SolverDirect::vmult(MPI::Vector &x, const MPI::Vector &b) const
  {
    // Assign the empty LHS vector to the Epetra_LinearProblem object
    linear_problem->SetLHS(&x.trilinos_vector());

    // Assign the RHS vector to the Epetra_LinearProblem object
    linear_problem->SetRHS(
      const_cast<Epetra_MultiVector *>(&b.trilinos_vector()));

    // First set whether we want to print the solver information to screen or
    // not.
    ConditionalOStream verbose_cout(std::cout,
                                    additional_data.output_solver_details);


    verbose_cout << "Starting solve" << std::endl;

    // Fetch return value of Amesos Solver functions
    int ierr = solver->Solve();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    // Finally, force the SolverControl object to report convergence
    solver_control.check(0, 0);
  }



  void
  SolverDirect::vmult(
    dealii::LinearAlgebra::distributed::Vector<double>       &x,
    const dealii::LinearAlgebra::distributed::Vector<double> &b) const
  {
    Epetra_Vector ep_x(View,
                       linear_problem->GetOperator()->OperatorDomainMap(),
                       x.begin());
    Epetra_Vector ep_b(View,
                       linear_problem->GetOperator()->OperatorRangeMap(),
                       const_cast<double *>(b.begin()));

    // Assign the empty LHS vector to the Epetra_LinearProblem object
    linear_problem->SetLHS(&ep_x);

    // Assign the RHS vector to the Epetra_LinearProblem object
    linear_problem->SetRHS(&ep_b);

    // First set whether we want to print the solver information to screen or
    // not.
    ConditionalOStream verbose_cout(std::cout,
                                    additional_data.output_solver_details);

    verbose_cout << "Starting solve" << std::endl;

    // Fetch return value of Amesos Solver functions
    int ierr = solver->Solve();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    // Finally, force the SolverControl object to report convergence
    solver_control.check(0, 0);
  }



  void
  SolverDirect::do_solve()
  {
    // Fetch return value of Amesos Solver functions
    int ierr;

    // First set whether we want to print the solver information to screen or
    // not.
    ConditionalOStream verbose_cout(std::cout,
                                    additional_data.output_solver_details);

    // Next allocate the Amesos solver, this is done in two steps, first we
    // create a solver Factory and generate with that the concrete Amesos
    // solver, if possible.
    Amesos Factory;

    AssertThrow(Factory.Query(additional_data.solver_type.c_str()),
                ExcMessage(
                  "You tried to select the solver type <" +
                  additional_data.solver_type +
                  "> but this solver is not supported by Trilinos either "
                  "because it does not exist, or because Trilinos was not "
                  "configured for its use."));

    solver.reset(
      Factory.Create(additional_data.solver_type.c_str(), *linear_problem));

    verbose_cout << "Starting symbolic factorization" << std::endl;
    ierr = solver->SymbolicFactorization();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    verbose_cout << "Starting numeric factorization" << std::endl;
    ierr = solver->NumericFactorization();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    verbose_cout << "Starting solve" << std::endl;
    ierr = solver->Solve();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    // Finally, let the deal.II SolverControl object know what has
    // happened. If the solve succeeded, the status of the solver control will
    // turn into SolverControl::success.
    solver_control.check(0, 0);

    if (solver_control.last_check() != SolverControl::success)
      AssertThrow(false,
                  SolverControl::NoConvergence(solver_control.last_step(),
                                               solver_control.last_value()));
  }


  void
  SolverDirect::solve(const Epetra_Operator    &A,
                      Epetra_MultiVector       &x,
                      const Epetra_MultiVector &b)
  {
    // We need an Epetra_LinearProblem object to let the Amesos solver know
    // about the matrix and vectors.
    linear_problem =
      std::make_unique<Epetra_LinearProblem>(const_cast<Epetra_Operator *>(&A),
                                             &x,
                                             const_cast<Epetra_MultiVector *>(
                                               &b));

    do_solve();
  }


  void
  SolverDirect::solve(const SparseMatrix       &sparse_matrix,
                      FullMatrix<double>       &solution,
                      const FullMatrix<double> &rhs)
  {
    Assert(sparse_matrix.m() == sparse_matrix.n(), ExcInternalError());
    Assert(rhs.m() == sparse_matrix.m(), ExcInternalError());
    Assert(rhs.m() == solution.m(), ExcInternalError());
    Assert(rhs.n() == solution.n(), ExcInternalError());

    const unsigned int m = rhs.m();
    const unsigned int n = rhs.n();

    FullMatrix<double> rhs_t(n, m);
    FullMatrix<double> solution_t(n, m);

    rhs_t.copy_transposed(rhs);
    solution_t.copy_transposed(solution);

    std::vector<double *> rhs_ptrs(n);
    std::vector<double *> sultion_ptrs(n);

    for (unsigned int i = 0; i < n; ++i)
      {
        rhs_ptrs[i]     = &rhs_t[i][0];
        sultion_ptrs[i] = &solution_t[i][0];
      }

    const Epetra_CrsMatrix &mat = sparse_matrix.trilinos_matrix();

    Epetra_MultiVector trilinos_dst(View,
                                    mat.OperatorRangeMap(),
                                    sultion_ptrs.data(),
                                    sultion_ptrs.size());
    Epetra_MultiVector trilinos_src(View,
                                    mat.OperatorDomainMap(),
                                    rhs_ptrs.data(),
                                    rhs_ptrs.size());

    this->initialize(sparse_matrix);
    this->solve(mat, trilinos_dst, trilinos_src);

    solution.copy_transposed(solution_t);
  }


  void
  SolverDirect::solve(const SparseMatrix &A,
                      MPI::Vector        &x,
                      const MPI::Vector  &b)
  {
    // We need an Epetra_LinearProblem object to let the Amesos solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()),
      &x.trilinos_vector(),
      const_cast<Epetra_MultiVector *>(&b.trilinos_vector()));

    do_solve();
  }



  void
  SolverDirect::solve(const SparseMatrix           &A,
                      dealii::Vector<double>       &x,
                      const dealii::Vector<double> &b)
  {
    // In case we call the solver with deal.II vectors, we create views of the
    // vectors in Epetra format.
    Assert(x.size() == A.n(), ExcDimensionMismatch(x.size(), A.n()));
    Assert(b.size() == A.m(), ExcDimensionMismatch(b.size(), A.m()));
    Assert(A.local_range().second == A.m(),
           ExcMessage("Can only work in serial when using deal.II vectors."));
    Epetra_Vector ep_x(View, A.trilinos_matrix().DomainMap(), x.begin());
    Epetra_Vector ep_b(View,
                       A.trilinos_matrix().RangeMap(),
                       const_cast<double *>(b.begin()));

    // We need an Epetra_LinearProblem object to let the Amesos solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()), &ep_x, &ep_b);

    do_solve();
  }



  void
  SolverDirect::solve(
    const SparseMatrix                                       &A,
    dealii::LinearAlgebra::distributed::Vector<double>       &x,
    const dealii::LinearAlgebra::distributed::Vector<double> &b)
  {
    AssertDimension(x.locally_owned_size(),
                    A.trilinos_matrix().DomainMap().NumMyElements());
    AssertDimension(b.locally_owned_size(),
                    A.trilinos_matrix().RangeMap().NumMyElements());
    Epetra_Vector ep_x(View, A.trilinos_matrix().DomainMap(), x.begin());
    Epetra_Vector ep_b(View,
                       A.trilinos_matrix().RangeMap(),
                       const_cast<double *>(b.begin()));

    // We need an Epetra_LinearProblem object to let the Amesos solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()), &ep_x, &ep_b);

    do_solve();
  }
} // namespace TrilinosWrappers


// explicit instantiations
// TODO: put these instantiations into generic file
namespace TrilinosWrappers
{
  template void
  SolverBase::do_solve(const PreconditionBase &preconditioner);

  template void
  SolverBase::do_solve(const Epetra_Operator &preconditioner);
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
