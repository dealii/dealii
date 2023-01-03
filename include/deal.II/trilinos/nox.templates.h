// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

          Assert(false, ExcInternalError());

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
         * Constructor. The class is intialized by the solution vector and
         * functions to compute the residual, to setup the jacobian, and
         * to solve the Jacobian.
         */
        Group(
          VectorType &                                                solution,
          const std::function<int(const VectorType &, VectorType &)> &residual,
          const std::function<int(const VectorType &)> &setup_jacobian,
          const std::function<int(const VectorType &, VectorType &)>
            &                                     apply_jacobian,
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
        computeX(const NOX::Abstract::Group & grp,
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
                      NOX::Abstract::Vector &      result) const override
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
         * A helper function to setup Jacobian.
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
                  const VectorType &x__ = *dynamic_cast<
                    const internal::NOXWrappers::Vector<VectorType> *>(
                    &problem.getSolutionGroup().getX());
                  const VectorType &f__ = *dynamic_cast<
                    const internal::NOXWrappers::Vector<VectorType> *>(
                    &problem.getSolutionGroup().getF());

                  // forward to the user-provided function and checks
                  // convergence
                  const auto state = this->check_iteration_status(
                    problem.getNumIterations(), f__.l2_norm(), x__, f__);

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
                                                 const VectorType & x,
                                                 const VectorType & f)>
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
    AdditionalData &                            additional_data,
    const Teuchos::RCP<Teuchos::ParameterList> &parameters)
    : additional_data(additional_data)
    , parameters(parameters)
    , n_residual_evaluations(0)
    , n_jacobian_applications(0)
    , n_nonlinear_iterations(0)
    , n_last_linear_iterations(0)
  {}



  template <typename VectorType>
  void
  NOXSolver<VectorType>::clear()
  {
    // clear interal counters
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
      solution,
      [&](const VectorType &x, VectorType &f) -> int {
        Assert(
          residual,
          ExcMessage(
            "No residual function has been attached to the NOXSolver object."));

        n_residual_evaluations++;

        // evalute residual
        return residual(x, f);
      },
      [&](const VectorType &x) -> int {
        Assert(
          setup_jacobian,
          ExcMessage(
            "No setup_jacobian function has been attached to the NOXSolver object."));

        // setup Jacobian
        int flag = setup_jacobian(x);

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

            if ((update_preconditioner == false) &&
                (update_preconditioner_predicate != nullptr))
              update_preconditioner = update_preconditioner_predicate();

            if (update_preconditioner) // update preconditioner
              flag = setup_preconditioner(x);
          }

        return flag;
      },
      [&](const VectorType &x, VectorType &v) -> int {
        Assert(
          apply_jacobian,
          ExcMessage(
            "No apply_jacobian function has been attached to the NOXSolver object."));

        n_jacobian_applications++;

        // apply Jacobian
        return apply_jacobian(x, v);
      },
      [&](const VectorType &f, VectorType &x, const double tolerance) -> int {
        n_nonlinear_iterations++;

        // invert Jacobian
        if (solve_with_jacobian)
          {
            Assert(
              !solve_with_jacobian_and_track_n_linear_iterations,
              ExcMessage(
                "It does not make sense to provide both solve_with_jacobian and "
                "solve_with_jacobian_and_track_n_linear_iterations!"));

            // without tracking of linear iterations
            return solve_with_jacobian(f, x, tolerance);
          }
        else if (solve_with_jacobian_and_track_n_linear_iterations)
          {
            // with tracking of linear iterations
            const int n_linear_iterations =
              solve_with_jacobian_and_track_n_linear_iterations(f,
                                                                x,
                                                                tolerance);

            if (n_linear_iterations == -1)
              return 1;

            this->n_last_linear_iterations = n_linear_iterations;

            return 0;
          }
        else
          {
            Assert(
              false,
              ExcMessage(
                "Neither a solve_with_jacobian or a "
                "solve_with_jacobian_and_track_n_linear_iterations function "
                "has been attached to the NOXSolver object."));

            Assert(false, ExcNotImplemented());
            return 1;
          }
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

    // solve
    const auto status = solver->solve();

    AssertThrow(status == NOX::StatusTest::Converged, ExcNOXNoConvergence());

    return solver->getNumIterations();
  }

} // namespace TrilinosWrappers

#  endif

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
