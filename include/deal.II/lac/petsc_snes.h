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

#ifndef dealii_petsc_snes_h
#define dealii_petsc_snes_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/parameter_handler.h>
#  include <deal.II/base/smartpointer.h>

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscdm.h>
#  include <petscsnes.h>

DEAL_II_NAMESPACE_OPEN

// Shorthand notation for PETSc error codes.
#  define AssertSNES(code)                               \
    do                                                   \
      {                                                  \
        PetscErrorCode __ierr = (code);                  \
        AssertThrow(__ierr == 0, ExcPETScError(__ierr)); \
      }                                                  \
    while (0)

namespace PETScWrappers
{
  /**
   * Additional parameters that can be passed to the NonlinearSolver class.
   */
  class NonlinearSolverData
  {
  public:
    /**
     * Initialization parameters for NonlinearSolverData.
     *
     * Running parameters:
     *
     * @param snestype Solver type
     * @param absolute_tolerance Absolute error tolerance
     * @param relative_tolerance Relative error tolerance
     * @param step_tolerance Step tolerance
     * @param max_it Maximum number of iterations
     * @param max_fe Maximum number of function evaluations
     *
     */
    NonlinearSolverData(
      // Running parameters
      const std::string &snestype           = std::string(),
      const double       absolute_tolerance = 0,
      const double       relative_tolerance = 0,
      const double       step_tolerance     = 0,
      const unsigned int max_it             = 0,
      const unsigned int max_fe             = 0)
      : snestype(snestype)
      , absolute_tolerance(absolute_tolerance)
      , relative_tolerance(relative_tolerance)
      , step_tolerance(step_tolerance)
      , max_it(max_it)
      , max_fe(max_fe)
    {}

    /**
     * Add all NonlinearSolverData() parameters to the given ParameterHandler
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
      prm.add_parameter("Solver type", snestype);
      prm.add_parameter("Absolute error tolerance", absolute_tolerance);
      prm.add_parameter("Relative error tolerance", relative_tolerance);
      prm.add_parameter("Step tolerance", step_tolerance);
      prm.add_parameter("Maximum iterations", max_it);
      prm.add_parameter("Maximum function evaluations", max_fe);
      prm.leave_subsection();
    }

    /**
     * PETSc solver type
     */
    std::string snestype;

    /**
     * Absolute error tolerance for function evaluation.
     */
    double absolute_tolerance;

    /**
     * Relative error tolerance for function evaluation.
     */
    double relative_tolerance;

    /**
     * Step tolerance for solution update.
     */
    double step_tolerance;

    /**
     * Maximum number of nonlinear iterations.
     */
    unsigned int max_it;

    /**
     * Maximum number of function evaluations.
     */
    unsigned int max_fe;
  };

  /**
   * Interface to PETSc SNES solver for nonlinear equations.
   */
  template <typename VectorType  = PETScWrappers::VectorBase,
            typename PMatrixType = PETScWrappers::MatrixBase,
            typename AMatrixType = PMatrixType>
  class NonlinearSolver
  {
  public:
    /**
     * Constructor.
     */
    NonlinearSolver(const NonlinearSolverData &data     = NonlinearSolverData(),
                    const MPI_Comm &           mpi_comm = PETSC_COMM_WORLD)
    {
      AssertSNES(SNESCreate(mpi_comm, &snes));
      AssertSNES(SNESSetApplicationContext(snes, this));
      reinit(data);
    }

    /**
     * Destructor.
     */
    virtual ~NonlinearSolver()
    {
      AssertSNES(SNESDestroy(&snes));
    }

    /**
     * Return the PETSc SNES object.
     */
    SNES
    petsc_snes()
    {
      return snes;
    }

    /**
     * Return a reference to the MPI communicator object in use with this
     * object.
     */
    const MPI_Comm &
    get_mpi_communicator() const
    {
      static MPI_Comm comm = PETSC_COMM_SELF;
      MPI_Comm pcomm = PetscObjectComm(reinterpret_cast<PetscObject>(snes));
      if (pcomm != MPI_COMM_NULL)
        comm = pcomm;
      return comm;
    }

    /**
     * Reset the subobjects associated with SNES, do not change customization.
     */
    void
    reinit()
    {
      AssertSNES(SNESReset(snes));
    }

    /**
     * Reset the subobjects associated with SNES, changes customization
     * according to data.
     */
    void
    reinit(const NonlinearSolverData &data)
    {
      reinit();
      if (data.snestype.size())
        AssertSNES(SNESSetType(snes, data.snestype.c_str()));

      PetscReal atol  = (PetscReal)data.absolute_tolerance > 0 ?
                          data.absolute_tolerance :
                          PETSC_DEFAULT;
      PetscReal rtol  = (PetscReal)data.relative_tolerance > 0 ?
                          data.relative_tolerance :
                          PETSC_DEFAULT;
      PetscReal stol  = (PetscReal)data.step_tolerance > 0 ?
                          data.step_tolerance :
                          PETSC_DEFAULT;
      PetscInt  maxit = (PetscInt)data.max_it > 0 ? data.max_it : PETSC_DEFAULT;
      PetscInt  maxfe = (PetscInt)data.max_fe > 0 ? data.max_fe : PETSC_DEFAULT;
      AssertSNES(SNESSetTolerances(snes, atol, rtol, stol, maxit, maxfe));
    }


    /**
     * Pass the preconditioning matrix only.
     * When used with 'setup_jacobian' and 'solve_for_jacobian_system', PETSc
     * will approximate the linear system matrix-vector product using an
     * internal matrix-free representation.
     * When used with 'jacobian', PETSc will use the same matrix for both
     * preconditioning and matrix-vector products.
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
     * Solve the nonlinear system of equations F(u) = 0
     */
    unsigned int
    solve(VectorType &solution)
    {
      auto _dealii_snes_function_ =
        [](SNES snes, Vec x, Vec f, void *ctx) -> PetscErrorCode {
        NonlinearSolver *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertSNES(SNESGetApplicationContext(snes, (void *)&myctx));
        auto user = static_cast<NonlinearSolver *>(myctx);

        VectorType xdealii(x);
        VectorType fdealii(f);
        user->residual(xdealii, fdealii);
        PetscFunctionReturn(0);
      };

      auto _dealii_snes_jacobian_ =
        [](SNES snes, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
        NonlinearSolver *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertSNES(SNESGetApplicationContext(snes, (void *)&myctx));
        auto user = static_cast<NonlinearSolver *>(myctx);

        VectorType  xdealii(x);
        AMatrixType Adealii(A);
        PMatrixType Pdealii(P);
        user->jacobian(xdealii, Adealii, Pdealii);

        PetscBool flg;
        AssertSNES(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
        if (flg)
          {
            AssertSNES(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            AssertSNES(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
          }
        PetscFunctionReturn(0);
      };

      auto _dealii_snes_jacobian_with_setup_ =
        [](SNES snes, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
        NonlinearSolver *myctx;

        PetscFunctionBeginUser;
        (void)ctx;
        AssertSNES(SNESGetApplicationContext(snes, (void *)&myctx));
        auto user = static_cast<NonlinearSolver *>(myctx);

        VectorType  xdealii(x);
        AMatrixType Adealii(A);
        PMatrixType Pdealii(P);

        user->A = &Adealii;
        user->P = &Pdealii;
        user->setup_jacobian(xdealii);

        PetscBool flg;
        AssertSNES(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
        if (flg)
          {
            AssertSNES(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            AssertSNES(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
          }
        PetscFunctionReturn(0);
      };

      auto _dealii_snes_monitor_ =
        [](SNES snes, PetscInt it, PetscReal f, void *ctx) -> PetscErrorCode {
        NonlinearSolver *myctx;
        PetscFunctionBeginUser;
        (void)ctx;
        AssertSNES(SNESGetApplicationContext(snes, (void *)&myctx));
        auto user = static_cast<NonlinearSolver *>(myctx);

        Vec x;
        AssertSNES(SNESGetSolution(snes, &x));
        VectorType xdealii(x);
        user->monitor(xdealii, it, f);
        PetscFunctionReturn(0);
      };

      Assert(residual, ExcFunctionNotProvided("residual"));

      AssertSNES(SNESSetSolution(snes, solution.petsc_vector()));

      AssertSNES(
        SNESSetFunction(snes, nullptr, _dealii_snes_function_, nullptr));

      if (setup_jacobian)
        {
          AssertSNES(SNESSetJacobian(snes,
                                     A ? A->petsc_matrix() : nullptr,
                                     P ? P->petsc_matrix() : nullptr,
                                     _dealii_snes_jacobian_with_setup_,
                                     nullptr));
          if (!A)
            AssertSNES(SNESSetUseMatrixFree(snes, PETSC_TRUE, PETSC_FALSE));

          if (!P)
            {
              DM dm;
              AssertSNES(SNESGetDM(snes, &dm));
              AssertSNES(DMSetMatType(dm, MATSHELL));
            }
        }
      else if (jacobian)
        {
          AssertSNES(SNESSetJacobian(snes,
                                     A ? A->petsc_matrix() :
                                         (P ? P->petsc_matrix() : nullptr),
                                     P ? P->petsc_matrix() : nullptr,
                                     _dealii_snes_jacobian_,
                                     nullptr));
        }


      // In case solve_for_jacobian_system is provided, create a shell
      // preconditioner wrapping the user call. The internal Krylov
      // solver will apply the preconditioner only once. This choice
      // can be overriden by command line and users can use any other
      // Krylov method if their solve is not accurate enough.
      PreconditionShell precond(
        PetscObjectComm(reinterpret_cast<PetscObject>(snes)));
      if (solve_for_jacobian_system)
        {
          precond.apply = [&](VectorBase &      indst,
                              const VectorBase &insrc) -> int {
            VectorType       dst(static_cast<const Vec &>(indst));
            const VectorType src(static_cast<const Vec &>(insrc));
            return solve_for_jacobian_system(src, dst);
          };

          KSP ksp;
          AssertSNES(SNESGetKSP(snes, &ksp));
          AssertSNES(KSPSetType(ksp, KSPPREONLY));
          AssertSNES(KSPSetPC(ksp, precond.get_pc()));
        }

      // Attach user monitoring routine
      if (monitor)
        AssertSNES(
          SNESMonitorSet(snes, _dealii_snes_monitor_, nullptr, nullptr));

      AssertSNES(SNESSetFromOptions(snes));

      AssertSNES(
        SNESSolve(snes, B ? static_cast<const Vec &>(*B) : nullptr, nullptr));

      PetscInt nt;
      AssertSNES(SNESGetIterationNumber(snes, &nt));
      return nt;
    }

    /**
     * Solve the nonlinear system of equations F(u) = 0.
     * Here we also pass the matrices to handle Jacobians.
     */
    unsigned int
    solve(VectorType &solution, PMatrixType &P)
    {
      reinit_matrices(P);
      return solve(solution);
    }

    /**
     * Solve the nonlinear system of equations F(u) = 0.
     * Here we also pass the matrices to handle Jacobians.
     */
    unsigned int
    solve(VectorType &solution, AMatrixType &A, PMatrixType &P)
    {
      reinit_matrices(A, P);
      return solve(solution);
    }

    /**
     * Solve the nonlinear system of equations F(u) = rhs
     */
    unsigned int
    solve(VectorType &solution, const VectorType &rhs)
    {
      this->B = &rhs;
      auto nt = solve(solution);
      this->B = nullptr;
      return nt;
    }

    /**
     * Solve the nonlinear system of equations F(u) = rhs.
     * Here we also pass the matrices to handle Jacobians.
     */
    unsigned int
    solve(VectorType &solution, const VectorType &rhs, PMatrixType &P)
    {
      this->B = &rhs;
      auto nt = solve(solution, P);
      this->B = nullptr;
      return nt;
    }

    /**
     * Solve the nonlinear system of equations F(u) = rhs.
     * Here we also pass the matrices to handle Jacobians.
     */
    unsigned int
    solve(VectorType &      solution,
          const VectorType &rhs,
          AMatrixType &     A,
          PMatrixType &     P)
    {
      this->B = &rhs;
      auto nt = solve(solution, A, P);
      this->B = nullptr;
      return nt;
    }

    /**
     * Compute nonlinear residual $F(x)$.
     */
    std::function<int(const VectorType &x, VectorType &res)> residual;

    /**
     * Compute Jacobian
     * $\dfrac{\partial F}{\partial x}$
     */
    std::function<int(const VectorType &x, AMatrixType &A, PMatrixType &P)>
      jacobian;

    /**
     * Process solution. This function is called by NonlinearSolver at the
     * beginning of each nonlinear step.
     */
    std::function<
      int(const VectorType &x, const unsigned int step_number, double fval)>
      monitor;

    /**
     * Setup Jacobian without matrices. This callback gives full control to
     * users to setup the linearized equations
     * $\dfrac{\partial F}{\partial x}$
     */
    std::function<int(const VectorType &x)> setup_jacobian;

    /**
     * Solution of the Jacobian system setup with setup_jacobian
     */
    std::function<int(const VectorType &src, VectorType &dst)>
      solve_for_jacobian_system;

  protected:
    SNES                                            snes;
    SmartPointer<AMatrixType, NonlinearSolver>      A;
    SmartPointer<PMatrixType, NonlinearSolver>      P;
    SmartPointer<const VectorType, NonlinearSolver> B;

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
