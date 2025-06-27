// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_nox_templates
#define dealii_trilinos_nox_templates

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_NOX

#  include <deal.II/trilinos/nox.h>

#  include <NOX_Abstract_Group.H>
#  include <NOX_Abstract_Vector.H>
#  include <NOX_Solver_Factory.H>
#  include <NOX_Solver_Generic.H>
#  include <NOX_StatusTest_Combo.H>
#  include <NOX_StatusTest_MaxIters.H>
#  include <NOX_StatusTest_NormF.H>
#  include <NOX_StatusTest_RelativeNormF.H>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN

namespace TrilinosWrappers
{
  namespace internal
  {
    namespace NOXWrappers
    {
      /**
       * A function that calls the function object given by its first argument
       * with the set of arguments following at the end. If the call returns
       * regularly, the current function returns zero to indicate success. If
       * the call fails with an exception, then the current function returns
       * with an error code of -1. In that case, the exception thrown by `f`
       * is captured and `eptr` is set to the exception. In case of success,
       * `eptr` is set to `nullptr`.
       */
      template <typename F, typename... Args>
      int
      call_and_possibly_capture_exception(const F            &f,
                                          std::exception_ptr &eptr,
                                          Args &&...args)
      {
        // See whether there is already something in the exception pointer
        // variable. There is no reason why this should be so, and
        // we should probably bail out:
        AssertThrow(eptr == nullptr, ExcInternalError());

        // Call the function and if that succeeds, return zero:
        try
          {
            f(std::forward<Args>(args)...);
            eptr = nullptr;
            return 0;
          }
        // In case of an exception, capture the exception and
        // return -1:
        catch (...)
          {
            eptr = std::current_exception();
            return -1;
          }
      }


      template <typename VectorType>
      class Group;

      /**
       * Implementation of the abstract interface
       * NOX::Abstract::Vector for deal.II vectors. For details,
       * see
       * https://docs.trilinos.org/dev/packages/nox/doc/html/classNOX_1_1Abstract_1_1Vector.html.
       */
      template <typename VectorType>
      class Vector : public NOX::Abstract::Vector
      {
      public:
        /**
         * Create empty vector.
         */
        Vector() = default;

        /**
         * Wrap an existing vector. The ownership is not transferred.
         */
        Vector(VectorType &vector)
        {
          this->vector.reset(&vector, [](auto *) { /*nothing to do*/ });
        }

        /**
         * Initialize every element of this vector with @p gamma .
         */
        NOX::Abstract::Vector &
        init(double gamma) override
        {
          *vector = static_cast<typename VectorType::value_type>(gamma);
          return *this;
        }

        /**
         * Initialize each element of this vector with a random value.
         */
        NOX::Abstract::Vector &
        random(bool useSeed = false, int seed = 1) override
        {
          AssertThrow(false, ExcNotImplemented());

          (void)useSeed;
          (void)seed;

          return *this;
        }

        /**
         * Put element-wise absolute values of source vector @p y into this vector.
         */
        NOX::Abstract::Vector &
        abs(const NOX::Abstract::Vector &y) override
        {
          AssertThrow(false, ExcNotImplemented());

          (void)y;

          return *this;
        }

        /**
         * Copy source vector @p y into this vector.
         */
        NOX::Abstract::Vector &
        operator=(const NOX::Abstract::Vector &y) override
        {
          if (vector == nullptr)
            vector = std::shared_ptr<VectorType>();

          const auto y_ = dynamic_cast<const Vector<VectorType> *>(&y);

          Assert(y_, ExcInternalError());

          vector->reinit(*y_->vector);

          *vector = *y_->vector;

          return *this;
        }

        /**
         * Put element-wise reciprocal of source vector @p y into this vector.
         */
        NOX::Abstract::Vector &
        reciprocal(const NOX::Abstract::Vector &y) override
        {
          AssertThrow(false, ExcNotImplemented());

          (void)y;

          return *this;
        }

        /**
         * Scale each element of this vector by @p gamma .
         */
        NOX::Abstract::Vector &
        scale(double gamma) override
        {
          *vector *= static_cast<typename VectorType::value_type>(gamma);

          return *this;
        }

        /**
         * Scale this vector element-by-element by the vector @p a.
         */
        NOX::Abstract::Vector &
        scale(const NOX::Abstract::Vector &a) override
        {
          const auto a_ = dynamic_cast<const Vector<VectorType> *>(&a);

          Assert(a_, ExcInternalError());

          vector->scale(*a_->vector);

          return *this;
        }

        /**
         * Compute `x = (alpha * a) + (gamma * x)` where `x` is this vector.
         */
        NOX::Abstract::Vector &
        update(double                       alpha,
               const NOX::Abstract::Vector &a,
               double                       gamma = 0.0) override
        {
          const auto a_ = dynamic_cast<const Vector<VectorType> *>(&a);

          Assert(a_, ExcInternalError());

          vector->sadd(static_cast<typename VectorType::value_type>(gamma),
                       static_cast<typename VectorType::value_type>(alpha),
                       *a_->vector);

          return *this;
        }

        /**
         * Compute `x = (alpha * a) + (beta * b) + (gamma * x)` where `x` is
         * this vector.
         */
        NOX::Abstract::Vector &
        update(double                       alpha,
               const NOX::Abstract::Vector &a,
               double                       beta,
               const NOX::Abstract::Vector &b,
               double                       gamma = 0.0) override
        {
          const auto a_ = dynamic_cast<const Vector<VectorType> *>(&a);
          const auto b_ = dynamic_cast<const Vector<VectorType> *>(&b);

          Assert(a_, ExcInternalError());
          Assert(b_, ExcInternalError());

          vector->operator*=(
            static_cast<typename VectorType::value_type>(gamma));
          vector->add(static_cast<typename VectorType::value_type>(alpha),
                      *a_->vector,
                      static_cast<typename VectorType::value_type>(beta),
                      *b_->vector);

          return *this;
        }

        /**
         * Create a new Vector of the same underlying type by cloning "this",
         * and return a pointer to the new vector.
         */
        Teuchos::RCP<NOX::Abstract::Vector>
        clone(NOX::CopyType copy_type) const override
        {
          auto new_vector    = Teuchos::rcp(new Vector<VectorType>());
          new_vector->vector = std::make_shared<VectorType>();
          new_vector->vector->reinit(*this->vector);

          if (copy_type == NOX::CopyType::DeepCopy)
            *new_vector->vector = *this->vector;
          else
            Assert(copy_type == NOX::CopyType::ShapeCopy, ExcInternalError());

          return new_vector;
        }

        /**
         * Return the vector norm.
         */
        double
        norm(NOX::Abstract::Vector::NormType type =
               NOX::Abstract::Vector::TwoNorm) const override
        {
          if (type == NOX::Abstract::Vector::NormType::TwoNorm)
            return vector->l2_norm();
          if (type == NOX::Abstract::Vector::NormType::OneNorm)
            return vector->l1_norm();
          if (type == NOX::Abstract::Vector::NormType::MaxNorm)
            return vector->linfty_norm();

          DEAL_II_ASSERT_UNREACHABLE();

          return 0.0;
        }

        /**
         * Return the vector's weighted 2-norm.
         */
        double
        norm(const NOX::Abstract::Vector &weights) const override
        {
          AssertThrow(false, ExcNotImplemented());

          (void)weights;

          return 0.0;
        }

        /**
         * Return the inner product of the vector with @p y.
         */
        double
        innerProduct(const NOX::Abstract::Vector &y) const override
        {
          const auto y_ = dynamic_cast<const Vector<VectorType> *>(&y);

          Assert(y_, ExcInternalError());

          return (*vector) * (*y_->vector);
        }

        /**
         * Return the length of vector.
         */
        NOX::size_type
        length() const override
        {
          return vector->size();
        }

        /**
         * Return underlying vector.
         */
        operator VectorType &() const
        {
          AssertThrow(vector, ExcInternalError());

          return *vector;
        }

      private:
        /**
         * The underlying deal.II vector.
         */
        std::shared_ptr<VectorType> vector;

        friend Group<VectorType>;
      };

      /**
       * Implementation of the abstract interface
       * NOX::Abstract::Group for deal.II vectors and deal.II solvers. For
       * details, see
       * https://docs.trilinos.org/dev/packages/nox/doc/html/classNOX_1_1Abstract_1_1Group.html.
       */
      template <typename VectorType>
      class Group : public NOX::Abstract::Group
      {
      public:
        /**
         * Constructor. The class is initialized by the solution vector and
         * functions to compute the residual, to set up the Jacobian, and
         * to solve the Jacobian.
         */
        Group(
          VectorType                                                 &solution,
          const std::function<int(const VectorType &, VectorType &)> &residual,
          const std::function<int(const VectorType &)> &setup_jacobian,
          const std::function<int(const VectorType &, VectorType &)>
                                                 &apply_jacobian,
          const std::function<int(const VectorType &,
                                  VectorType &,
                                  const double)> &solve_with_jacobian)
          : x(solution)
          , residual(residual)
          , setup_jacobian(setup_jacobian)
          , apply_jacobian(apply_jacobian)
          , solve_with_jacobian(solve_with_jacobian)
          , is_valid_f(false)
          , is_valid_j(false)
        {}

        /**
         * Copies the source group into this group.
         */
        NOX::Abstract::Group &
        operator=(const NOX::Abstract::Group &source) override
        {
          if (this != &source)
            {
              const auto other =
                dynamic_cast<const Group<VectorType> *>(&source);

              Assert(other, ExcInternalError());

              if (other->x.vector)
                {
                  if (this->x.vector == nullptr)
                    this->x.vector = std::make_shared<VectorType>();

                  *this->x.vector = *other->x.vector;
                }
              else
                {
                  this->x.vector = {};
                }

              if (other->f.vector)
                {
                  if (this->f.vector == nullptr)
                    this->f.vector = std::make_shared<VectorType>();

                  *this->f.vector = *other->f.vector;
                }
              else
                {
                  this->f.vector = {};
                }

              if (other->gradient.vector)
                {
                  if (this->gradient.vector == nullptr)
                    this->gradient.vector = std::make_shared<VectorType>();

                  *this->gradient.vector = *other->gradient.vector;
                }
              else
                {
                  this->gradient.vector = {};
                }

              if (other->newton.vector)
                {
                  if (this->newton.vector == nullptr)
                    this->newton.vector = std::make_shared<VectorType>();

                  *this->newton.vector = *other->newton.vector;
                }
              else
                {
                  this->newton.vector = {};
                }

              this->residual            = other->residual;
              this->setup_jacobian      = other->setup_jacobian;
              this->apply_jacobian      = other->apply_jacobian;
              this->solve_with_jacobian = other->solve_with_jacobian;

              this->is_valid_f = other->is_valid_f;
              this->is_valid_j = other->is_valid_j;
            }

          return *this;
        }

        /**
         * Set the solution vector `x` to @p y.
         */
        void
        setX(const NOX::Abstract::Vector &y) override
        {
          reset();

          x = y;
        }

        /**
         * Compute the solution update `x = grp.x + step * d`.
         */
        void
        computeX(const NOX::Abstract::Group  &grp,
                 const NOX::Abstract::Vector &d,
                 double                       step) override
        {
          reset();

          const auto grp_ = dynamic_cast<const Group *>(&grp);

          Assert(grp_, ExcInternalError());

          x.update(1.0, grp_->x, step, d);
        }

        /**
         * Compute and store the residual update `F(x)`.
         */
        NOX::Abstract::Group::ReturnType
        computeF() override
        {
          if (isF() == false)
            {
              f.vector = std::make_shared<VectorType>();
              f.vector->reinit(*x.vector);

              if (residual(*x.vector, *f.vector) != 0)
                return NOX::Abstract::Group::Failed;

              is_valid_f = true;
            }

          return NOX::Abstract::Group::Ok;
        }

        /**
         * Return `true` if the residual vector `F` is valid.
         */
        bool
        isF() const override
        {
          return is_valid_f;
        }

        /**
         * Compute and store the Jacobian.
         */
        NOX::Abstract::Group::ReturnType
        computeJacobian() override
        {
          if (isJacobian() == false)
            {
              if (setup_jacobian(*x.vector) != 0)
                return NOX::Abstract::Group::Failed;

              is_valid_j = true;
            }

          return NOX::Abstract::Group::Ok;
        }

        /**
         * Return `true` if the Jacobian is valid.
         */
        bool
        isJacobian() const override
        {
          return is_valid_j;
        }

        /**
         * Return the (total) solution vector.
         */
        const NOX::Abstract::Vector &
        getX() const override
        {
          return x;
        }

        /**
         * Return the residual `F(x)`.
         */
        const NOX::Abstract::Vector &
        getF() const override
        {
          return f;
        }

        /**
         * Return the 2-norm of `F(x)`.
         */
        double
        getNormF() const override
        {
          return f.norm();
        }

        /**
         * Return the gradient.
         */
        const NOX::Abstract::Vector &
        getGradient() const override
        {
          return gradient;
        }

        /**
         * Return the Newton descent direction.
         */
        const NOX::Abstract::Vector &
        getNewton() const override
        {
          return newton;
        }

        /**
         * Return a reference counting pointer to solution vector.
         */
        Teuchos::RCP<const NOX::Abstract::Vector>
        getXPtr() const override
        {
          AssertThrow(false, ExcNotImplemented());
          return {};
        }

        /**
         * Return a reference counting pointer to the residual `F(x)`.
         */
        Teuchos::RCP<const NOX::Abstract::Vector>
        getFPtr() const override
        {
          AssertThrow(false, ExcNotImplemented());
          return {};
        }

        /**
         * Return a reference counting pointer to gradient.
         */
        Teuchos::RCP<const NOX::Abstract::Vector>
        getGradientPtr() const override
        {
          AssertThrow(false, ExcNotImplemented());
          return {};
        }

        /**
         * Return a reference counting pointer to the Newton descent direction.
         */
        Teuchos::RCP<const NOX::Abstract::Vector>
        getNewtonPtr() const override
        {
          AssertThrow(false, ExcNotImplemented());
          return {};
        }

        /**
         * Create a new Group of the same derived type as this one by
         * cloning this one, and return a reference counting pointer to
         * the new group.
         */
        Teuchos::RCP<NOX::Abstract::Group>
        clone(NOX::CopyType copy_type) const override
        {
          auto new_group =
            Teuchos::rcp(new Group<VectorType>(*x.vector,
                                               residual,
                                               setup_jacobian,
                                               apply_jacobian,
                                               solve_with_jacobian));

          if (x.vector)
            {
              new_group->x.vector = std::make_shared<VectorType>();
              new_group->x.vector->reinit(*x.vector);
            }

          if (f.vector)
            {
              new_group->f.vector = std::make_shared<VectorType>();
              new_group->f.vector->reinit(*f.vector);
            }

          if (gradient.vector)
            {
              new_group->gradient.vector = std::make_shared<VectorType>();
              new_group->gradient.vector->reinit(*gradient.vector);
            }

          if (newton.vector)
            {
              new_group->newton.vector = std::make_shared<VectorType>();
              new_group->newton.vector->reinit(*newton.vector);
            }

          if (copy_type == NOX::CopyType::DeepCopy)
            {
              if (x.vector)
                *new_group->x.vector = *x.vector;

              if (f.vector)
                *new_group->f.vector = *f.vector;

              if (gradient.vector)
                *new_group->gradient.vector = *gradient.vector;

              if (newton.vector)
                *new_group->newton.vector = *newton.vector;

              new_group->is_valid_f = is_valid_f;
              new_group->is_valid_j = is_valid_j;
            }
          else
            Assert(copy_type == NOX::CopyType::ShapeCopy, ExcInternalError());

          return new_group;
        }

        /**
         * Compute the Newton direction, using the chosen
         * parameters for the linear solve.
         */
        NOX::Abstract::Group::ReturnType
        computeNewton(Teuchos::ParameterList &p) override
        {
          if (isNewton())
            return NOX::Abstract::Group::Ok;

          if (isF() == false || isJacobian() == false)
            return NOX::Abstract::Group::BadDependency;

          if (newton.vector == nullptr)
            newton.vector = std::make_shared<VectorType>();

          newton.vector->reinit(*f.vector, false);

          const double tolerance = p.get<double>("Tolerance");

          if (solve_with_jacobian(*f.vector, *newton.vector, tolerance) != 0)
            return NOX::Abstract::Group::NotConverged;

          newton.scale(-1.0);

          return NOX::Abstract::Group::Ok;
        }

        /**
         * Applies the Jacobian to the given @p input vector and assigns
         * the output to the @p result.
         */
        NOX::Abstract::Group::ReturnType
        applyJacobian(const NOX::Abstract::Vector &input,
                      NOX::Abstract::Vector       &result) const override
        {
          if (apply_jacobian == nullptr)
            return NOX::Abstract::Group::NotDefined;

          if (!isJacobian())
            return NOX::Abstract::Group::BadDependency;

          const auto *input_ = dynamic_cast<const Vector<VectorType> *>(&input);
          const auto *result_ =
            dynamic_cast<const Vector<VectorType> *>(&result);

          if (apply_jacobian(*input_->vector, *result_->vector) != 0)
            return NOX::Abstract::Group::Failed;

          return NOX::Abstract::Group::Ok;
        }

      private:
        /**
         * Reset the state of this object.
         */
        void
        reset()
        {
          is_valid_f = false;
          is_valid_j = false;
        }

        /**
         * An internal vector for the current solution.
         */
        Vector<VectorType> x;

        /**
         * An internal vector for the residual.
         */
        Vector<VectorType> f;

        /**
         * An internal vector for the solution gradient.
         */
        Vector<VectorType> gradient;

        /**
         * An internal vector for the newton step.
         */
        Vector<VectorType> newton;

        /**
         * A helper function to compute residual.
         */
        std::function<int(const VectorType &x, VectorType &f)> residual;

        /**
         * A helper function to set up Jacobian.
         */
        std::function<int(const VectorType &x)> setup_jacobian;

        /**
         * A helper function to apply Jacobian.
         */
        std::function<int(const VectorType &x, VectorType &v)> apply_jacobian;

        /**
         * A helper function to solve Jacobian.
         */
        std::function<
          int(const VectorType &f, VectorType &x, const double tolerance)>
          solve_with_jacobian;

        /**
         * A flag that indicates if the has residual been computed.
         */
        bool is_valid_f;

        /**
         * A flag that indicates if the Jacobian has been computed.
         */
        bool is_valid_j;
      };


      /**
       * Wrapper class around the user function that allows to check
       * convergence.
       */
      template <typename VectorType>
      class NOXCheck : public NOX::StatusTest::Generic
      {
      public:
        /**
         * Constructor.
         */
        NOXCheck(const std::function<SolverControl::State(const unsigned int,
                                                          const double,
                                                          const VectorType &,
                                                          const VectorType &)>
                   check_iteration_status)
          : check_iteration_status(check_iteration_status)
          , status(NOX::StatusTest::Unevaluated)
        {}

        /**
         * Check the status of the nonlinear solver.
         */
        NOX::StatusTest::StatusType
        checkStatus(const NOX::Solver::Generic &problem,
                    NOX::StatusTest::CheckType  checkType) override
        {
          if (checkType == NOX::StatusTest::None)
            {
              status = NOX::StatusTest::Unevaluated;
            }
          else
            {
              if (check_iteration_status == nullptr)
                {
                  status = NOX::StatusTest::Converged;
                }
              else
                {
                  // unwrap the various vectors
                  const VectorType &x = *dynamic_cast<
                    const internal::NOXWrappers::Vector<VectorType> *>(
                    &problem.getSolutionGroup().getX());
                  const VectorType &f = *dynamic_cast<
                    const internal::NOXWrappers::Vector<VectorType> *>(
                    &problem.getSolutionGroup().getF());

                  // forward to the user-provided function and check
                  // convergence
                  const auto state = this->check_iteration_status(
                    problem.getNumIterations(), f.l2_norm(), x, f);

                  // translate the returned value back to Trilinos data
                  // structure
                  switch (state)
                    {
                      case SolverControl::iterate:
                        status = NOX::StatusTest::Unconverged;
                        break;
                      case SolverControl::failure:
                        status = NOX::StatusTest::Failed;
                        break;
                      case SolverControl::success:
                        status = NOX::StatusTest::Converged;
                        break;
                      default:
                        AssertThrow(false, ExcNotImplemented());
                    }
                }
            }

          return status;
        }

        /**
         * Return the last value that was given by `checkStatus()`.
         */
        NOX::StatusTest::StatusType
        getStatus() const override
        {
          return status;
        }

        /**
         * Print the last value that was given by `checkStatus()`
         */
        virtual std::ostream &
        print(std::ostream &stream, int indent = 0) const override
        {
          for (int j = 0; j < indent; ++j)
            stream << ' ';
          stream << status << std::endl;
          return stream;
        }

      private:
        /**
         * The user function that allows the solver to check for
         * convergence.
         */
        const std::function<SolverControl::State(const unsigned int i,
                                                 const double       f_norm,
                                                 const VectorType  &x,
                                                 const VectorType  &f)>
          check_iteration_status;

        /**
         * The last returned value of `checkStatus()`.
         */
        NOX::StatusTest::StatusType status;
      };
    } // namespace NOXWrappers
  }   // namespace internal



  template <typename VectorType>
  NOXSolver<VectorType>::AdditionalData::AdditionalData(
    const unsigned int max_iter,
    const double       abs_tol,
    const double       rel_tol,
    const unsigned int threshold_nonlinear_iterations,
    const unsigned int threshold_n_linear_iterations,
    const bool         reuse_solver)
    : max_iter(max_iter)
    , abs_tol(abs_tol)
    , rel_tol(rel_tol)
    , threshold_nonlinear_iterations(threshold_nonlinear_iterations)
    , threshold_n_linear_iterations(threshold_n_linear_iterations)
    , reuse_solver(reuse_solver)
  {}



  template <typename VectorType>
  NOXSolver<VectorType>::NOXSolver(
    AdditionalData                             &additional_data,
    const Teuchos::RCP<Teuchos::ParameterList> &parameters)
    : additional_data(additional_data)
    , parameters(parameters)
    , n_residual_evaluations(0)
    , n_jacobian_applications(0)
    , n_nonlinear_iterations(0)
    , n_last_linear_iterations(0)
    , pending_exception(nullptr)
  {}


  template <typename VectorType>
  NOXSolver<VectorType>::~NOXSolver()
  {
    Assert(pending_exception == nullptr, ExcInternalError());
  }



  template <typename VectorType>
  void
  NOXSolver<VectorType>::clear()
  {
    // clear internal counters
    n_residual_evaluations   = 0;
    n_jacobian_applications  = 0;
    n_nonlinear_iterations   = 0;
    n_last_linear_iterations = 0;
  }



  template <typename VectorType>
  unsigned int
  NOXSolver<VectorType>::solve(VectorType &solution)
  {
    if (additional_data.reuse_solver == false)
      clear(); // clear state

    // create group
    const auto group = Teuchos::rcp(new internal::NOXWrappers::Group<
                                    VectorType>(
      /* Starting vector */
      solution,

      /* Residual function */
      [&](const VectorType &x, VectorType &f) -> int {
        Assert(
          residual,
          ExcMessage(
            "No residual function has been attached to the NOXSolver object."));

        ++n_residual_evaluations;

        // evaluate residual
        return internal::NOXWrappers::call_and_possibly_capture_exception(
          residual, pending_exception, x, f);
      },

      /* setup_jacobian function */
      [&](const VectorType &x) -> int {
        Assert(
          setup_jacobian,
          ExcMessage(
            "No setup_jacobian function has been attached to the NOXSolver object."));

        // setup Jacobian
        int flag = internal::NOXWrappers::call_and_possibly_capture_exception(
          setup_jacobian, pending_exception, x);

        if (flag != 0)
          return flag;

        if (setup_preconditioner)
          {
            // check if preconditioner needs to be updated
            bool update_preconditioner =
              ((additional_data.threshold_nonlinear_iterations > 0) &&
               ((n_nonlinear_iterations %
                 additional_data.threshold_nonlinear_iterations) == 0)) ||
              (solve_with_jacobian_and_track_n_linear_iterations &&
               (n_last_linear_iterations >
                additional_data.threshold_n_linear_iterations));

            if (update_preconditioner_predicate)
              update_preconditioner |= update_preconditioner_predicate();

            if (update_preconditioner)
              flag = internal::NOXWrappers::call_and_possibly_capture_exception(
                setup_preconditioner, pending_exception, x);
          }

        return flag;
      },

      /* apply_jacobian function */
      [&](const VectorType &x, VectorType &v) -> int {
        Assert(
          apply_jacobian,
          ExcMessage(
            "No apply_jacobian function has been attached to the NOXSolver object."));

        ++n_jacobian_applications;

        // apply Jacobian
        return internal::NOXWrappers::call_and_possibly_capture_exception(
          apply_jacobian, pending_exception, x, v);
      },

      /* solve_with_jacobian function */
      [&](const VectorType &f, VectorType &x, const double tolerance) -> int {
        ++n_nonlinear_iterations;

        int ret_code = 0;

        // invert Jacobian
        if (solve_with_jacobian)
          {
            Assert(
              !solve_with_jacobian_and_track_n_linear_iterations,
              ExcMessage(
                "It does not make sense to provide both solve_with_jacobian and "
                "solve_with_jacobian_and_track_n_linear_iterations!"));

            // without tracking of linear iterations
            ret_code =
              internal::NOXWrappers::call_and_possibly_capture_exception(
                solve_with_jacobian, pending_exception, f, x, tolerance);
          }
        else if (solve_with_jacobian_and_track_n_linear_iterations)
          {
            // With tracking of linear iterations. The callback function we
            // have here is using an integer return type, which is not
            // what internal::NOXWrappers::call_and_possibly_capture_exception
            // knows how to deal with. As a consequence, package the call
            // and assignment of the return value into a lambda function;
            // if the call to the user-defined callback fails, this will
            // trigger an exception that will propagate through the lambda
            // function and be treated correctly by the logic in
            // internal::NOXWrappers::call_and_possibly_capture_exception.
            ret_code =
              internal::NOXWrappers::call_and_possibly_capture_exception(
                [&]() {
                  this->n_last_linear_iterations =
                    solve_with_jacobian_and_track_n_linear_iterations(
                      f, x, tolerance);
                },
                pending_exception);
          }
        else
          {
            Assert(
              false,
              ExcMessage(
                "Neither a solve_with_jacobian or a "
                "solve_with_jacobian_and_track_n_linear_iterations function "
                "has been attached to the NOXSolver object."));

            DEAL_II_NOT_IMPLEMENTED();
            ret_code = 1;
          }

        // NOX has a recovery feature that is enabled by default. In this case,
        // if a solve_with_jacobian or a
        // solve_with_jacobian_and_track_n_linear_iterations function triggers
        // an exception and therefore call_and_possibly_capture_exception
        // returns code different from 0, then NOX does not interrupt the
        // solution process but rather performs a recovery step. To ensure this
        // feature is available to the user, we need to suppress the exception
        // in this case, since it is exactly that, what NOX expects from our
        // callbacks.
        const bool do_rescue =
          parameters->sublist("Newton").get("Rescue Bad Newton Solve", true);
        if (do_rescue && (pending_exception != nullptr))
          {
            try
              {
                std::rethrow_exception(pending_exception);
              }
            catch (const RecoverableUserCallbackError &exc)
              {
                pending_exception = nullptr;

                // If the callback threw a recoverable exception, and if
                // recovery is enabled, then eat the exception and return the
                // error code.
                return ret_code;
              }
            catch (...)
              {
                // If not a recoverable exception, then just re-throw the
                // exception and hope that NOX knows what to do with propagating
                // exceptions (i.e., does not create a resource leak, for
                // example).
                pending_exception = nullptr;
                throw;
              }
          }
        else
          // Rescue not allowed, or there was no exception -> simply return
          // the value produced by the callback.
          return ret_code;
      }));

    // setup solver control
    auto check =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

    if (this->check_iteration_status)
      {
        const auto info =
          Teuchos::rcp(new internal::NOXWrappers::NOXCheck<VectorType>(
            this->check_iteration_status));
        check->addStatusTest(info);
      }

    if (additional_data.abs_tol > 0.0)
      {
        const auto additional_data_norm_f_abs =
          Teuchos::rcp(new NOX::StatusTest::NormF(additional_data.abs_tol));
        check->addStatusTest(additional_data_norm_f_abs);
      }

    if (additional_data.rel_tol > 0.0)
      {
        const auto additional_data_norm_f_rel = Teuchos::rcp(
          new NOX::StatusTest::RelativeNormF(additional_data.rel_tol));
        check->addStatusTest(additional_data_norm_f_rel);
      }

    if (additional_data.max_iter > 0)
      {
        const auto additional_data_max_iterations =
          Teuchos::rcp(new NOX::StatusTest::MaxIters(additional_data.max_iter));
        check->addStatusTest(additional_data_max_iterations);
      }

    // create non-linear solver
    const auto solver = NOX::Solver::buildSolver(group, check, parameters);

    // Solve, then check whether an exception was thrown by one of the user
    // callback functions. If so, exit by rethrowing the exception that
    // we had previously saved. This also calls the destructors of all of
    // the member variables above, so we do not have to clean things up by hand.
    //
    // NOX has the annoying habit of reporting success through a return code,
    // but if a user callback function throws an exception (which we translate
    // into a return code of -1), then it throws its own exception. So we
    // have to check both the return code and also catch exceptions :-(
    try
      {
        const auto status = solver->solve();

        if (status == NOX::StatusTest::Converged)
          return solver->getNumIterations();
        else
          {
            // See if NOX aborted because we had thrown an exception in a user
            // callback:
            if (pending_exception)
              {
                std::exception_ptr this_exception = pending_exception;
                pending_exception                 = nullptr;

                std::rethrow_exception(this_exception);
              }

            // If that was not the case, NOX just didn't converge:
            AssertThrow(status == NOX::StatusTest::Converged,
                        ExcNOXNoConvergence());
          }
      }
      // See if NOX returned by triggering an exception.
#    if DEAL_II_TRILINOS_VERSION_GTE(14, 2, 0)
    // Starting with Trilinos version 14.2.0 NOX started to throw a
    // std::runtime_error instead of a plain char*.
    catch (const std::runtime_error &exc)
      {
        const char *s = exc.what();
#    else
    // In a sign of generally poor software design, NOX prior to Trilinos
    // version 14.2.0 throws an exception that is not of a class derived
    // from std::exception, but just a char*. That's a nuisance -- you just
    // have to know :-(
    catch (const char *s)
      {
#    endif
        // Like above, see if NOX aborted because there was an exception
        // in a user callback. In that case, collate the errors if we can
        // (namely, if the user exception was derived from std::exception),
        // and otherwise just let the user exception propagate (then swallowing
        // the NOX exception):
        if (pending_exception)
          {
            std::exception_ptr this_exception = pending_exception;
            pending_exception                 = nullptr;

            try
              {
                std::rethrow_exception(this_exception);
              }
            catch (const std::exception &e)
              {
                // Collate the exception texts:
                throw ExcMessage(
                  "NOX aborted with an error text of <" + std::string(s) +
                  "> after a user callback function had thrown an exception " +
                  "with the following message:\n" + e.what());
              }
            catch (...)
              {
                // Let user exception propagate if of a different type:
                throw;
              }
          }

        // NOX just happened to throw an exception, but it wasn't because there
        // was a user callback exception before. Convert the char* to something
        // more readable:
        AssertThrow(false,
                    ExcMessage("NOX aborted with an error text of <" +
                               std::string(s) + ">."));
      }

    return 0; // unreachable
  }

} // namespace TrilinosWrappers

#  endif

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
