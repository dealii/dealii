//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//    Author: Toby D. Young, Polish Academy of Sciences, 2008, 2009
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__slepc_solver_h
#define __deal2__slepc_solver_h


#include <base/config.h>

#ifdef DEAL_II_USE_SLEPC

#  include <base/std_cxx1x/shared_ptr.hpp>
#  include <lac/exceptions.h>
#  include <lac/solver_control.h>
#  include <lac/petsc_matrix_base.h>
#  include <lac/slepc_spectral_transformation.h>

#  include <petscksp.h>
#  include <slepceps.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Base class for solver classes using the SLEPc solvers which are
 * selected based on flags passed to the eigenvalue problem solver
 * context. Derived classes set the right flags to set the right
 * solver. On the other hand, note that: the
 * <code>AdditionalData</code> structure is a dummy structure that
 * currently exists for backward/forward compatibility only.
 *
 * The SLEPc solvers are intended to be used for solving the
 * generalized eigenspectrum problem $(A-\lambda M)x=0$, for $x\neq0$;
 * where $A$ is a system matrix, $M$ is a mass matrix, and $\lambda,
 * x$ are a set of eigenvalues and eigenvectors respectively. The
 * emphasis is on methods and techniques appropriate for problems in
 * which the associated matrices are sparse and therefore, most of the
 * methods offered by the SLEPc library are projection methods or
 * other methods with similar properties. SLEPc implements these basic
 * methods as well as more sophisticated algorithms. On the other
 * hand, SLEPc is a general library in the sense that it covers
 * standard and generalized eigenvalue problems, and wrappers are
 * provided to interface to SLEPc solvers that handle both of these
 * problem sets.
 *
 * SLEPcWrappers can be implemented in application codes in the
 * following way:
 @verbatim
  SolverControl solver_control (1000, 1e-9);
  SolverArnoldi system (solver_control,
                        mpi_communicator);
  system.solve (A, B, eigenvalues, eigenvectors,
                size_of_spectrum);
 @endverbatim

 * for the generalized eigenvalue problem $Ax=B\lambda x$, where the
 * variable <code>const unsigned int size_of_spectrum</code> tells
 * SLEPc the number of eigenvector/eigenvalue pairs to solve for: See
 * also step-36 for a hands-on example.
 *
 * An alternative implementation to the one above is to use the API
 * internals directly within the application code. In this way the
 * calling sequence requires calling several of SolverBase functions
 * rather than just one. This freedom is intended for use of the
 * SLEPcWrappers that require a greater handle on the eigenvalue
 * problem solver context. See also the API of:
 @verbatim
  template <typename OutputVector>
  void
  SolverBase::solve (const PETScWrappers::MatrixBase &A,
                     const PETScWrappers::MatrixBase &B,
                     std::vector<double>             &kr,
                     std::vector<OutputVector>       &vr,
                     const unsigned int               n_eigenvectors)
  { ... }
 @endverbatim
 * as an example on how to do this.
 *
 * For further information and explanations on handling the @ref
 * SLEPcWrappers "SLEPcWrappers", see also the @ref PETScWrappers
 * "PETScWrappers", on which they depend.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2008, 2009, 2010; and Rickard Armiento 2008
 *
 * Various tweaks to the SLEPcWrappers have been contributed by Eloy
 * Romeoro and Jose Roman.
 */
namespace SLEPcWrappers
{

                                   /**
                                    * Base class for solver classes
                                    * using the SLEPc solvers. Since
                                    * solvers in SLEPc are selected
                                    * based on flags passed to a
                                    * generic solver object, basically
                                    * all the actual solver calls
                                    * happen in this class, and
                                    * derived classes simply set the
                                    * right flags to select one solver
                                    * or another, or to set certain
                                    * parameters for individual
                                    * solvers.
                                    */
  class SolverBase
    {
    public:

                                   /**
                                    * Constructor. Takes the MPI
                                    * communicator over which parallel
                                    * computations are to happen.
                                    */
      SolverBase (SolverControl &cn,
		  const MPI_Comm &mpi_communicator);

                                   /**
                                    * Destructor.
                                    */
      virtual ~SolverBase ();

                                   /**
                                    * Composite method that solves the
                                    * eigensystem $Ax=\lambda x$. The
                                    * eigenvector sent in has to have
                                    * at least one element that we can
                                    * use as a template when resizing,
                                    * since we do not know the
                                    * parameters of the specific
                                    * vector class used
                                    * (i.e. local_dofs for MPI
                                    * vectors). However, while copying
                                    * eigenvectors, at least twice the
                                    * memory size of <tt>vr</tt> is
                                    * being used (and can be more). To
                                    * avoid doing this, the fairly
                                    * standard calling sequence
                                    * excecuted here is used:
                                    * Initialise; Set up matrices for
                                    * solving; Actually solve the
                                    * system; Gather the solution(s);
                                    * and reset.
                                    *
                                    * Note that the number of
                                    * converged eigenvectors can be
                                    * larger than the number of
                                    * eigenvectors requested; this is
                                    * due to a round off error
                                    * (success) of the eigenproblem
                                    * solver context. If this is found
                                    * to be the case we simply do not
                                    * bother with more eigenpairs than
                                    * requested, but handle that it
                                    * may be more than specified by
                                    * ignoring any extras. By default
                                    * one eigenvector/eigenvalue pair
                                    * is computed.
				    */
      template <typename OutputVector>
	void
	solve (const PETScWrappers::MatrixBase &A,
	       std::vector<double>             &kr,
	       std::vector<OutputVector>       &vr,
	       const unsigned int              n_eigenvectors);

                                   /**
				    * Same as above, but here a
                                    * composite method for solving the
                                    * system $A x=\lambda B x$.
                                    */
      template <typename OutputVector>
	void
	solve (const PETScWrappers::MatrixBase &A,
	       const PETScWrappers::MatrixBase &B,
	       std::vector<double>             &kr,
	       std::vector<OutputVector>       &vr,
	       const unsigned int               n_eigenvectors);

                                   /**
                                    * Initialize solver for the linear
                                    * system $Ax=\lambda x$. (Note:
                                    * this is required before calling
                                    * solve ())
                                    */
      void
	set_matrices (const PETScWrappers::MatrixBase &A);

                                   /**
				    * Same as above, but here
                                    * initialize solver for the linear
                                    * system $A x=\lambda B x$.
                                    */
      void
	set_matrices (const PETScWrappers::MatrixBase &A,
		      const PETScWrappers::MatrixBase &B);

                                   /**
                                    * Set the initial vector for the
                                    * solver.
                                    */
      void
	set_initial_vector (const PETScWrappers::VectorBase &initial_vec);

                                   /**
				    * Set the spectral transformation
				    * to be used.
                                    */
      void
	set_transformation (SLEPcWrappers::TransformationBase &trans);

                                   /**
				    * Indicate which part of the
				    * spectrum is to be computed. By
				    * default largest magnitude
				    * eigenvalues are computed.  For
				    * other allowed values see the
				    * SLEPc documentation.
                                    */
      void
	set_which_eigenpairs (EPSWhich set_which);

                                   /**
                                    * Solve the linear system for
				    * n_eigenvectors
				    * eigenstates. Parameter
				    * <code>n_converged</code>
				    * contains the actual number of
				    * eigenstates that have .
				    * converged; this can be both
				    * fewer or more than
				    * n_eigenvectors, depending on the
				    * SLEPc eigensolver used.
                                    */
      void
	solve (const unsigned int n_eigenvectors, unsigned int *n_converged);


                                   /**
                                    * Access the solutions for a
                                    * solved eigenvector problem, pair
                                    * index solutions,
                                    * $\text{index}\,\in\,0\hdots
                                    * \text{n\_converged}-1$.
                                    */
      void
	get_eigenpair (const unsigned int index,
#ifndef DEAL_II_USE_PETSC_COMPLEX
		       double                    &kr,
#else
		       std::complex<double>      &kr,
#endif
		       PETScWrappers::VectorBase &vr);

                                   /**
                                    * Reset the solver, and return
                                    * memory for eigenvectors
                                    */
      void
	reset();

                                   /**
				    * Retrieve the SLEPc solver object
				    * that is internally used.
                                    */
      EPS *
	get_EPS ();


                                   /**
                                    * Access to the object that
                                    * controls convergence.
                                    */
      SolverControl &control () const;

                                   /**
                                    * Exceptions.
                                    */
      DeclException0 (ExcSLEPcWrappersUsageError);

      DeclException1 (ExcSLEPcError,
		      int,
		      << "    An error with error number " << arg1
		      << " occured while calling a SLEPc function");

      DeclException2 (ExcSLEPcEigenvectorConvergenceMismatchError,
		      int, int,
		      << "    The number of converged eigenvectors is " << arg1
		      << " but " << arg2 << " were requested. ");

    protected:

                                   /**
                                    * Reference to the object that
                                    * controls convergence of the
                                    * iterative solver.
                                    */
      SolverControl &solver_control;

                                   /**
                                    * Copy of the MPI communicator
                                    * object to be used for the
                                    * solver.
                                    */
      const MPI_Comm mpi_communicator;

                                   /**
                                    * Function that takes an
                                    * Eigenvalue Problem Solver
                                    * context object, and sets the
                                    * type of solver that is requested
                                    * by the derived class.
                                    */
      virtual void set_solver_type (EPS &eps) const = 0;

                                   /**
                                    * Attributes that store the
				    * relevant information for the
				    * eigenproblem solver context.
                                    */
      EPSWhich                           set_which;

      const PETScWrappers::MatrixBase   *opA;
      const PETScWrappers::MatrixBase   *opB;
      const PETScWrappers::VectorBase   *ini_vec;

      SLEPcWrappers::TransformationBase *transform;

    private:

                                   /**
                                    * A function that is used in SLEPc
                                    * as a callback to check on
                                    * convergence. It takes the
                                    * information provided from SLEPc
                                    * and checks it against deal.II's
                                    * own SolverControl objects to see
                                    * if convergence has been reached.
                                    */
      static
	int
	convergence_test (EPS                 eps,
#ifdef PETSC_USE_64BIT_INDICES
			  const PetscInt      iteration,
#else
			  const int           iteration,
#endif
			  const PetscReal     residual_norm,
			  EPSConvergedReason *reason,
			  void               *solver_control);

                                   /**
                                    * Objects of this type are
                                    * explicitly created, but are
                                    * destroyed when the surrounding
                                    * solver object goes out of scope,
                                    * or when we assign a new value to
                                    * the pointer to this object. The
                                    * respective Destroy functions are
                                    * therefore written into the
                                    * destructor of this object, even
                                    * though the object does not have
                                    * a constructor.
                                    */
      struct SolverData
      {

                                   /**
                           	    * Destructor.
                              	    */
	~SolverData ();

                                   /**
                                    * Objects for Eigenvalue Problem
                                    * Solver.
                                    */
	EPS eps;
      };

      std_cxx1x::shared_ptr<SolverData> solver_data;
    };

/**
 * An implementation of the solver interface using the SLEPc
 * Krylov-Schur solver. Usage: All spectrum, all problem types,
 * complex.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2008
 */
  class SolverKrylovSchur : public SolverBase
    {
    public:

                                   /**
                                    * Standardized data struct to pipe
                                    * additional data to the solver,
                                    * should it be needed.
                                    */
      struct AdditionalData
      {};

                                   /**
                                    * SLEPc solvers will want to have
                                    * an MPI communicator context over
                                    * which computations are
                                    * parallelized. By default, this
                                    * carries the same behaviour has
                                    * the PETScWrappers, but you can
                                    * change that.
                                    */
       SolverKrylovSchur (SolverControl        &cn,
			  const MPI_Comm       &mpi_communicator = PETSC_COMM_WORLD,
			  const AdditionalData &data = AdditionalData());

    protected:

                                   /**
                                    * Store a copy of the flags for
                                    * this particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Eigenvalue
                                    * Problem Solver context object,
                                    * and sets the type of solver that
                                    * is appropriate for this class.
                                    */
      virtual void set_solver_type (EPS &eps) const;
    };

/**
 * An implementation of the solver interface using the SLEPc Arnoldi
 * solver. Usage: All spectrum, all problem types, complex.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2008
 */
  class SolverArnoldi : public SolverBase
    {
    public:
                                   /**
                                    * Standardized data struct to pipe
                                    * additional data to the solver,
                                    * should it be needed.
                                    */
      struct AdditionalData
      {};

                                   /**
                                    * SLEPc solvers will want to have
                                    * an MPI communicator context over
                                    * which computations are
                                    * parallelized. By default, this
                                    * carries the same behaviour has
                                    * the PETScWrappers, but you can
                                    * change that.
                                    */
      SolverArnoldi (SolverControl        &cn,
		     const MPI_Comm       &mpi_communicator = PETSC_COMM_WORLD,
		     const AdditionalData &data = AdditionalData());

    protected:

                                   /**
                                    * Store a copy of the flags for
                                    * this particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Eigenvalue
                                    * Problem Solver context object,
                                    * and sets the type of solver that
                                    * is appropriate for this class.
                                    */
      virtual void set_solver_type (EPS &eps) const;

    };

/**
 * An implementation of the solver interface using the SLEPc Lanczos
 * solver. Usage: All spectrum, all problem types, complex.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2009
 */
  class SolverLanczos : public SolverBase
    {
    public:
                                   /**
                                    * Standardized data struct to pipe
                                    * additional data to the solver,
                                    * should it be needed.
                                    */
      struct AdditionalData
      {};

                                   /**
                                    * SLEPc solvers will want to have
                                    * an MPI communicator context over
                                    * which computations are
                                    * parallelized. By default, this
                                    * carries the same behaviour has
                                    * the PETScWrappers, but you can
                                    * change that.
                                    */
      SolverLanczos (SolverControl        &cn,
		     const MPI_Comm       &mpi_communicator = PETSC_COMM_WORLD,
		     const AdditionalData &data = AdditionalData());

    protected:

                                   /**
                                    * Store a copy of the flags for
                                    * this particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Eigenvalue
                                    * Problem Solver context object,
                                    * and sets the type of solver that
                                    * is appropriate for this class.
                                    */
      virtual void set_solver_type (EPS &eps) const;

    };


  // --------------------------- inline and template functions -----------
  /**
   * This is declared here to make it possible to take a std::vector
   * of different PETScWrappers vector types
   */

  template <typename OutputVector>
    void
    SolverBase::solve (const PETScWrappers::MatrixBase &A,
                       std::vector<double>             &kr,
                       std::vector<OutputVector>       &vr,
                       const unsigned int               n_eigenvectors = 1)
    {
      unsigned int n_converged;

      set_matrices (A);
      solve (n_eigenvectors,&n_converged);

      if (n_converged > n_eigenvectors)
	n_converged = n_eigenvectors;
      AssertThrow (n_converged == n_eigenvectors,
		   ExcSLEPcEigenvectorConvergenceMismatchError(n_converged, n_eigenvectors));

      AssertThrow (vr.size() >= 1, ExcSLEPcWrappersUsageError());
      vr.resize (n_converged, vr.front());
      kr.resize (n_converged);

      for (unsigned int index=0; index < n_converged;
	   +index)
	get_eigenpair (index, kr[index], vr[index]);
    }

  template <typename OutputVector>
    void
    SolverBase::solve (const PETScWrappers::MatrixBase &A,
                       const PETScWrappers::MatrixBase &B,
                       std::vector<double>             &kr,
                       std::vector<OutputVector>       &vr,
                       const unsigned int               n_eigenvectors = 1)
    {
      unsigned int n_converged;

      set_matrices (A,B);
      solve (n_eigenvectors, &n_converged);

      if (n_converged >= n_eigenvectors)
	n_converged = n_eigenvectors;

      AssertThrow (n_converged == n_eigenvectors,
		   ExcSLEPcEigenvectorConvergenceMismatchError(n_converged, n_eigenvectors));
      AssertThrow (vr.size() >= 1, ExcSLEPcWrappersUsageError());

      vr.resize (n_converged, vr.front());
      kr.resize (n_converged);

      for (unsigned int index=0; index < n_converged;
	   ++index)
	get_eigenpair (index, kr[index], vr[index]);
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_SLEPC

/*----------------------------   slepc_solver.h  ---------------------------*/

#endif

/*----------------------------   slepc_solver.h  ---------------------------*/

