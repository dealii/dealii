//-----------------------------------------------------------
//
//    Copyright (C) 2022 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------
//
// Author: Stefano Zampini, King Abdullah University of Science and Technology.

#ifndef dealii_petsc_ts_h
#define dealii_petsc_ts_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/parameter_handler.h>
#  include <deal.II/base/smartpointer.h>

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscdm.h>
#  include <petscts.h>

DEAL_II_NAMESPACE_OPEN

// Shorthand notation for PETSc error codes.
#  define AssertTS(code)                                 \
    do                                                   \
      {                                                  \
        PetscErrorCode __ierr = (code);                  \
        AssertThrow(__ierr == 0, ExcPETScError(__ierr)); \
      }                                                  \
    while (0)

namespace PETScWrappers
{
  /**
   * Additional parameters that can be passed to the TimeStepper class.
   */
  class TimeStepperData
  {
  public:
    /**
     * Initialization parameters for TimeStepper.
     *
     * Running parameters:
     *
     * @param tstype Solver type
     * @param initial_time Initial time
     * @param final_time Final time
     * @param initial_step_size Initial step size
     * @param max_steps maximum number of steps
     * @param match_step match requested final time
     *
     * Error parameters:
     *
     * @param minimum_step_size Minimum step size
     * @param maximum_step_size Maximum step size
     * @param absolute_tolerance Absolute error tolerance
     * @param relative_tolerance Relative error tolerance
     *
     * Note that either one of final_time or max_steps must
     * be specified by the user, otherwise PETSc may complain.
     * Adaptive time stepping is disabled by default. To
     * activate it, customize tolerance with nonzero values.
     * Negative values indicate using PETSc's default.
     */
    TimeStepperData(
      // Running parameters
      const std::string &tstype            = std::string(),
      const double       initial_time      = 0.0,
      const double       final_time        = 0.0,
      const double       initial_step_size = 0.0,
      const int          max_steps         = -1,
      const bool         match_step        = false,
      // Error parameters
      const double minimum_step_size  = 0.0,
      const double maximum_step_size  = 0.0,
      const double absolute_tolerance = 0.0,
      const double relative_tolerance = 0.0)
      : tstype(tstype)
      , initial_time(initial_time)
      , final_time(final_time)
      , initial_step_size(initial_step_size)
      , max_steps(max_steps)
      , match_step(match_step)
      , minimum_step_size(minimum_step_size)
      , maximum_step_size(maximum_step_size)
      , absolute_tolerance(absolute_tolerance)
      , relative_tolerance(relative_tolerance)
    {}

    /**
     * Add all TimeStepperData() parameters to the given ParameterHandler
     * object. When the parameters are parsed from a file, the internal
     * parameters are automatically updated.
     *
     * These are one-to-one with the options you can pass at construction
     * time.
     *
     * The options you pass at construction time are set as default values in
     * the ParameterHandler object `prm`. You can later modify them by parsing
     * a parameter file using `prm`. The values of the parameter will be
     * updated whenever the content of `prm` is updated.
     *
     * Make sure that this class lives longer than `prm`. Undefined behavior
     * will occur if you destroy this class, and then parse a parameter file
     * using `prm`.
     */
    void
    add_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Running parameters");
      prm.add_parameter("Solver type", tstype);
      prm.add_parameter("Initial time", initial_time);
      prm.add_parameter("Final time", final_time);
      prm.add_parameter("Initial step size", initial_step_size);
      prm.add_parameter("Maximum number of steps", max_steps);
      prm.add_parameter("Match final time", match_step);
      prm.leave_subsection();

      prm.enter_subsection("Error control");
      prm.add_parameter("Minimum step size", minimum_step_size);
      prm.add_parameter("Maximum step size", maximum_step_size);
      prm.add_parameter("Absolute error tolerance", absolute_tolerance);
      prm.add_parameter("Relative error tolerance", relative_tolerance);
      prm.leave_subsection();
    }

    /**
     * PETSc solver type
     */
    std::string tstype;

    /**
     * Initial time for the DAE.
     */
    double initial_time;

    /**
     * Final time.
     */
    double final_time;

    /**
     * Initial step size.
     */
    double initial_step_size;

    /**
     * Maximum number of steps to be taken.
     * Negative values are ignored.
     */
    int max_steps;

    /**
     * Match final time requested?
     */
    bool match_step;

    /**
     * Minimum allowed step size.
     * Non-positive values indicate to use PETSc's default.
     */
    double minimum_step_size;

    /**
     * Maximum allowed step size.
     * Non-positive values indicate to use PETSc's default.
     */
    double maximum_step_size;

    /**
     * Absolute error tolerance for adaptive time stepping.
     * Negative values indicate to use PETSc's default.
     */
    double absolute_tolerance;

    /**
     * Relative error tolerance for adaptive time stepping.
     * Negative values indicate to use PETSc's default.
     */
    double relative_tolerance;
  };

  /**
   * Interface to PETSc TS solver for Ordinary Differential Equations and
   * Differential-Algebraic Equations.
   */
  template <typename VectorType  = PETScWrappers::VectorBase,
            typename PMatrixType = PETScWrappers::MatrixBase,
            typename AMatrixType = PMatrixType>
  class TimeStepper
  {
  public:
    /**
     * Constructor.
     */
    TimeStepper(const TimeStepperData &data     = TimeStepperData(),
                const MPI_Comm &       mpi_comm = PETSC_COMM_WORLD)
    {
      AssertTS(TSCreate(mpi_comm, &ts));
      AssertTS(TSSetApplicationContext(ts, this));
      reinit(data);
    }

    /**
     * Destructor.
     */
    virtual ~TimeStepper()
    {
      AssertTS(TSDestroy(&ts));
    }

    /**
     * Return the PETSc TS object.
     */
    TS
    petsc_ts()
    {
      return ts;
    }

    /**
     * Return a reference to the MPI communicator object in use with this
     * object.
     */
    const MPI_Comm &
    get_mpi_communicator() const
    {
      static MPI_Comm comm = PETSC_COMM_SELF;
      MPI_Comm pcomm       = PetscObjectComm(reinterpret_cast<PetscObject>(ts));
      if (pcomm != MPI_COMM_NULL)
        comm = pcomm;
      return comm;
    }

    /**
     * Reset the subobjects associated with TS, do not change customization.
     */
    void
    reinit()
    {
      AssertTS(TSReset(ts));
    }

    /**
     * Reset the subobjects associated with TS, changes customization according
     * to data.
     */
    void
    reinit(const TimeStepperData &data)
    {
      reinit();

      // Solver type
      if (data.tstype.size())
        AssertTS(TSSetType(ts, data.tstype.c_str()));

      // Time and steps limits
      AssertTS(TSSetTime(ts, data.initial_time));
      if (data.final_time > data.initial_time)
        AssertTS(TSSetMaxTime(ts, data.final_time));
      if (data.initial_step_size > 0.0)
        AssertTS(TSSetTimeStep(ts, data.initial_step_size));
      if (data.max_steps >= 0)
        AssertTS(TSSetMaxSteps(ts, data.max_steps));
      // Decide how to end the integration. Either stepover the final time or
      // match it.
      AssertTS(TSSetExactFinalTime(ts,
                                   data.match_step ?
                                     TS_EXACTFINALTIME_MATCHSTEP :
                                     TS_EXACTFINALTIME_STEPOVER));

      // Adaptive tolerances
      PetscReal atol =
        data.absolute_tolerance > 0 ? data.absolute_tolerance : PETSC_DEFAULT;
      PetscReal rtol =
        data.relative_tolerance > 0 ? data.relative_tolerance : PETSC_DEFAULT;
      AssertTS(TSSetTolerances(ts, atol, nullptr, rtol, nullptr));

      // Adaptive time stepping
      TSAdapt tsadapt;
      AssertTS(TSGetAdapt(ts, &tsadapt));

      // Use adaptive time stepping only if explicitly requested
      if (data.absolute_tolerance == 0 && data.relative_tolerance == 0)
        {
          AssertTS(TSAdaptSetType(tsadapt, TSADAPTNONE));
        }
      else
        {
          AssertTS(TSAdaptSetType(tsadapt, TSADAPTBASIC));
        }

      // Step limits
      PetscReal hmin =
        data.minimum_step_size > 0 ? data.minimum_step_size : PETSC_DEFAULT;
      PetscReal hmax =
        data.maximum_step_size > 0 ? data.maximum_step_size : PETSC_DEFAULT;
      AssertTS(TSAdaptSetStepLimits(tsadapt, hmin, hmax));
    }


    /**
     * Pass the preconditioning matrix only.
     * When used with 'setup_jacobian' and 'solve_for_jacobian_system', PETSc
     * will approximate the linear system matrix-vector product using an
     * internal matrix-free representation.
     * When used with 'implicit_jacobian' or 'explicit_jacobian',
     * PETSc will use the same matrix for both preconditioning and
     * matrix-vector products.
     */
    void
    reinit_matrices(PMatrixType &P)
    {
      this->A = nullptr;
      this->P = &P;
    }

    /**
     * Pass both the linear system matrix and the preconditioning matrix.
     */
    void
    reinit_matrices(AMatrixType &A, PMatrixType &P)
    {
      this->A = &A;
      this->P = &P;
    }

    /**
     * Return current time
     */
    double
    get_time()
    {
      PetscReal t;
      AssertTS(TSGetTime(ts, &t));
      return t;
    }

    /**
     * Return current time step
     */
    double
    get_time_step()
    {
      PetscReal dt;
      AssertTS(TSGetTimeStep(ts, &dt));
      return dt;
    }

    /**
     * Integrate the differential-algebraic equations.
     * This function returns the final number of computed steps.
     */
    unsigned int
    solve(VectorType &solution)
    {
      auto _dealii_ts_ifunction_ =
        [](TS ts, PetscReal t, Vec x, Vec xdot, Vec f, void *ctx)
        -> PetscErrorCode {
        TimeStepper *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertTS(TSGetApplicationContext(ts, (void *)&myctx));
        auto user = static_cast<TimeStepper *>(myctx);

        VectorType xdealii(x);
        VectorType xdotdealii(xdot);
        VectorType fdealii(f);
        user->implicit_function(t, xdealii, xdotdealii, fdealii);
        PetscFunctionReturn(0);
      };

      auto _dealii_ts_ijacobian_ = [](TS        ts,
                                      PetscReal t,
                                      Vec       x,
                                      Vec       xdot,
                                      PetscReal s,
                                      Mat       A,
                                      Mat       P,
                                      void *    ctx) -> PetscErrorCode {
        TimeStepper *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertTS(TSGetApplicationContext(ts, (void *)&myctx));
        auto user = static_cast<TimeStepper *>(myctx);

        VectorType  xdealii(x);
        VectorType  xdotdealii(xdot);
        AMatrixType Adealii(A);
        PMatrixType Pdealii(P);
        user->implicit_jacobian(t, xdealii, xdotdealii, s, Adealii, Pdealii);
        PetscBool flg;
        AssertTS(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
        if (flg)
          {
            AssertTS(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            AssertTS(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
          }
        PetscFunctionReturn(0);
      };

      auto _dealii_ts_ijacobian_with_setup_ = [](TS        ts,
                                                 PetscReal t,
                                                 Vec       x,
                                                 Vec       xdot,
                                                 PetscReal s,
                                                 Mat       A,
                                                 Mat       P,
                                                 void *ctx) -> PetscErrorCode {
        TimeStepper *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertTS(TSGetApplicationContext(ts, (void *)&myctx));
        auto user = static_cast<TimeStepper *>(myctx);

        VectorType  xdealii(x);
        VectorType  xdotdealii(xdot);
        AMatrixType Adealii(A);
        PMatrixType Pdealii(P);

        user->A = &Adealii;
        user->P = &Pdealii;
        user->setup_jacobian(t, xdealii, xdotdealii, s);

        PetscBool flg;
        AssertTS(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
        if (flg)
          {
            AssertTS(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            AssertTS(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
          }
        PetscFunctionReturn(0);
      };

      auto _dealii_ts_rhsfunction_ =
        [](TS ts, PetscReal t, Vec x, Vec f, void *ctx) -> PetscErrorCode {
        TimeStepper *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertTS(TSGetApplicationContext(ts, (void *)&myctx));
        auto user = static_cast<TimeStepper *>(myctx);

        VectorType xdealii(x);
        VectorType fdealii(f);

        user->explicit_function(t, xdealii, fdealii);
        PetscFunctionReturn(0);
      };

      auto _dealii_ts_rhsjacobian_ =
        [](TS ts, PetscReal t, Vec x, Mat A, Mat P, void *ctx)
        -> PetscErrorCode {
        TimeStepper *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertTS(TSGetApplicationContext(ts, (void *)&myctx));
        auto user = static_cast<TimeStepper *>(myctx);

        VectorType  xdealii(x);
        AMatrixType Adealii(A);
        PMatrixType Pdealii(P);

        user->explicit_jacobian(t, xdealii, Adealii, Pdealii);
        PetscBool flg;
        AssertTS(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
        if (flg)
          {
            AssertTS(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            AssertTS(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
          }
        PetscFunctionReturn(0);
      };

      auto _dealii_ts_monitor_ = [](TS        ts,
                                    PetscInt  it,
                                    PetscReal t,
                                    Vec       x,
                                    void *    ctx) -> PetscErrorCode {
        TimeStepper *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertTS(TSGetApplicationContext(ts, (void *)&myctx));
        auto user = static_cast<TimeStepper *>(myctx);

        VectorType xdealii(x);
        user->monitor(t, xdealii, it);
        PetscFunctionReturn(0);
      };

      Assert(explicit_function || implicit_function,
             ExcFunctionNotProvided("explicit_function || implicit_function"));

      AssertTS(TSSetSolution(ts, solution.petsc_vector()));

      if (explicit_function)
        AssertTS(
          TSSetRHSFunction(ts, nullptr, _dealii_ts_rhsfunction_, nullptr));

      if (implicit_function)
        AssertTS(TSSetIFunction(ts, nullptr, _dealii_ts_ifunction_, nullptr));

      if (setup_jacobian)
        {
          AssertTS(TSSetIJacobian(ts,
                                  A ? A->petsc_matrix() : nullptr,
                                  P ? P->petsc_matrix() : nullptr,
                                  _dealii_ts_ijacobian_with_setup_,
                                  nullptr));
          SNES snes;
          AssertTS(TSGetSNES(ts, &snes));
          if (!A)
            AssertTS(SNESSetUseMatrixFree(snes, PETSC_TRUE, PETSC_FALSE));

          if (!P)
            {
              DM dm;
              AssertTS(SNESGetDM(snes, &dm));
              AssertTS(DMSetMatType(dm, MATSHELL));
            }
        }
      else
        {
          if (explicit_jacobian)
            {
              AssertTS(TSSetRHSJacobian(ts,
                                        A ? A->petsc_matrix() :
                                            (P ? P->petsc_matrix() : nullptr),
                                        P ? P->petsc_matrix() : nullptr,
                                        _dealii_ts_rhsjacobian_,
                                        nullptr));
            }

          if (implicit_jacobian)
            {
              AssertTS(TSSetIJacobian(ts,
                                      A ? A->petsc_matrix() :
                                          (P ? P->petsc_matrix() : nullptr),
                                      P ? P->petsc_matrix() : nullptr,
                                      _dealii_ts_ijacobian_,
                                      nullptr));
            }
        }


      // In case solve_for_jacobian_system is provided, create a shell
      // preconditioner wrapping the user call. The internal Krylov
      // solver will apply the preconditioner only once. This choice
      // can be overriden by command line and users can use any other
      // Krylov method if their solve is not accurate enough.
      PreconditionShell precond(
        PetscObjectComm(reinterpret_cast<PetscObject>(ts)));
      if (solve_for_jacobian_system)
        {
          precond.apply = [&](VectorBase &      indst,
                              const VectorBase &insrc) -> int {
            VectorType       dst(static_cast<const Vec &>(indst));
            const VectorType src(static_cast<const Vec &>(insrc));
            return solve_for_jacobian_system(src, dst);
          };

          SNES snes;
          KSP  ksp;
          AssertTS(TSGetSNES(ts, &snes));
          AssertTS(SNESGetKSP(snes, &ksp));
          AssertTS(KSPSetType(ksp, KSPPREONLY));
          AssertTS(KSPSetPC(ksp, precond.get_pc()));
        }

      // Attach user monitoring routine
      if (monitor)
        AssertTS(TSMonitorSet(ts, _dealii_ts_monitor_, nullptr, nullptr));

      AssertTS(TSSetFromOptions(ts));

      AssertTS(TSSolve(ts, nullptr));

      PetscInt nt;
      AssertTS(TSGetStepNumber(ts, &nt));
      return nt;
    }

    /**
     * Integrate the differential-algebraic equations.
     * Here we also pass the matrices to handle Jacobians.
     */
    unsigned int
    solve(VectorType &solution, PMatrixType &P)
    {
      reinit_matrices(P);
      return solve(solution);
    }

    /**
     * Integrate the differential-algebraic equations.
     * Here we also pass the matrices to handle Jacobians.
     */
    unsigned int
    solve(VectorType &solution, AMatrixType &A, PMatrixType &P)
    {
      reinit_matrices(A, P);
      return solve(solution);
    }

    /**
     * Compute implicit residual $F(t, y, \dot y)$.
     */
    std::function<int(const double      t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      VectorType &      res)>
      implicit_function;

    /**
     * Compute implicit Jacobian
     * $\dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot
     * y}$
     */
    std::function<int(const double      t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      const double      alpha,
                      AMatrixType &     A,
                      PMatrixType &     P)>
      implicit_jacobian;

    /**
     * Compute explicit residual $G(t, y)$.
     */
    std::function<int(const double t, const VectorType &y, VectorType &res)>
      explicit_function;

    /**
     * Compute explicit Jacobian $\dfrac{\partial G}{\partial y}$.
     */
    std::function<
      int(const double t, const VectorType &y, AMatrixType &A, PMatrixType &P)>
      explicit_jacobian;

    /**
     * Process solution. This function is called by TimeStepper at the beginning
     * of each time step.
     */
    std::function<
      int(const double t, const VectorType &y, const unsigned int step_number)>
      monitor;

    /**
     * Setup Jacobian without matrices. This callback gives full control to
     * users to setup the linearized equations
     * $\dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot
     * y}$
     */
    std::function<int(const double      t,
                      const VectorType &y,
                      const VectorType &ydot,
                      const double      shift)>
      setup_jacobian;

    /**
     * Solution of the Jacobian system setup with setup_jacobian
     */
    std::function<int(const VectorType &src, VectorType &dst)>
      solve_for_jacobian_system;

  protected:
    TS                                     ts;
    SmartPointer<AMatrixType, TimeStepper> A;
    SmartPointer<PMatrixType, TimeStepper> P;

  private:
    DeclException1(ExcFunctionNotProvided,
                   std::string,
                   << "Please provide an implementation for the function \""
                   << arg1 << "\"");
  };
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
