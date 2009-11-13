//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_precondition_block_h
#define __deal2__trilinos_precondition_block_h


#include <base/config.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <lac/exceptions.h>
#  include <lac/trilinos_precondition.h>
#  include <lac/trilinos_sparse_matrix.h>

#  include <Teuchos_RCP.hpp>
#  include <Thyra_LinearOperatorDecl.hpp>
#  include <Epetra_Operator.h>

// some forward declarations
class Ifpack_Preconditioner;


DEAL_II_NAMESPACE_OPEN

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{
				   // forward declarations
  class BlockSparseMatrix;

/**
 * This implements an interface class of preconditioners for block
 * matrices based on Trilinos, which are intended to be used with the
 * solver classes SolverBlockCG and SolverBlockGMRES. This is based on
 * the Trilinos Thyra abstract interface of linear operators. This
 * class does not implement any method, and derived classes will need
 * to do that job.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionBlockBase : public Subscriptor
  {
    public:
				       /**
					* Empty constructor.
					*/
      PreconditionBlockBase();

				       /**
					* Destructor.
					*/
      ~PreconditionBlockBase();

    protected:
				       /**
					* Preconditioner object.
					*/
      Teuchos::RCP<Thyra::ConstLinearOperator<double> > preconditioner;

    friend class SolverBlockBase;
};

				       /**
					* Internal function that sets
					* up an inverse matrix.
					*/
  inline
  Thyra::ConstLinearOperator<double>
  inverse_matrix (const SparseMatrix    &M,
		  const Epetra_Operator *P,
		  const bool             is_symmetric,
		  const unsigned int     n_iterations,
		  const double           solve_tolerance,
		  const bool             output_details);


					// Forward declarations.

/**
 * This class implements a black box preconditioner for saddle points
 * systems arising from the Stokes or Navier&ndash;Stokes equations as
 * specified by the papers <em>D. Silvester, A. Wathen, Fast iterative
 * solution of stabilised Stokes systems part II. Using general block
 * preconditioners, SIAM J. Numer. Anal. 31:1352&ndash;1367 (1994)</em>
 * and <em> D. Kay, D. Loghin, A. Wathen, A preconditioner for the
 * steady-state Navier&ndash;Stokes equations, SIAM
 * J. Sci. Comput. 24(1):237&ndash;256 (2002)</em>, respectively.
 *
 * The preconditioner is based an approximation to the Schur
 * complement of the block matrix. The Schur complement $S=B
 * A_{\mathbf u}^{-1} B^T$ is approximated by a mass matrix $M_p$ on
 * the pressure space in the case of the Stokes equations, and as a
 * product $S^{-1} = L_p^{-1} F_p M_p^{-1}$ with pressure Laplace
 * matrix $L_p$, pressure convection-diffusion operator $F_p$
 * (corresponding to the sum of time derivative, convection and
 * diffusion), and pressure mass matrix $M_p$ in the case of the
 * Navier&ndash;Stokes equations.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionStokes : public PreconditionBlockBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {
                                       /**
					* Constructor. This sets the
					* additional data to its
					* default values. For lazy
					* initializations, the user
					* will just set the parameters
					* in here and let the
					* <tt>initialize</tt> function
					* set up good preconditioners,
					* where some good values have
					* to be provided here. The
					* default ones mean that we
					* assume the velocity-velocity
					* matrix to be symmetric (only
					* true for Stokes problems!),
					* we use a relative inner
					* tolerance of 1e-9 for the
					* inversion of matrices, and
					* take in IC preconditioner
					* for the inversion of the
					* pressure mass
					* matrix. Moreover, there is
					* no inner CG/GMRES solve
					* performed on the
					* velocity-velocity block by
					* default. Otherwise, we do a
					* solve with the specified
					* solver and the number of
					* iterations and
					* tolerance. Note that too
					* many iterations on Av can
					* slow down the overall
					* procedure considerably.
					*/
	AdditionalData (const bool          right_preconditioning = false,
			const bool          Av_is_symmetric = true,
			const std::vector<std::vector<bool> > &Av_constant_modes = std::vector<std::vector<bool> > (1),
			const double        inner_solve_tolerance = 1e-9,
			const bool          use_ssor_on_mass_matrix = false,
			const bool          velocity_uses_higher_order_elements = false,
			const bool          pressure_uses_higher_order_elements = false,
			const bool          Av_do_solve   = false,
			const unsigned int  Av_n_iters    = 1,
			const double        Av_tolerance  = 1e-3,
			const bool          ouput_details = false);

				       /**
					* This flag specifies whether
					* the preconditioner should be
					* build as a left
					* preconditioner or right
					* preconditioner. Note that
					* this setting must be the
					* same as the one used in the
					* solver call.
					*/
	bool right_preconditioning;

				       /**
					* This flag determines whether
					* the velocity-velocity matrix
					* is symmetric (and positive
					* definite) or not. This
					* choice will influence the
					* AMG preconditioner that is
					* built for that system as
					* well as the inner solve on
					* Av (if that flag is
					* enabled).
					*/
	bool Av_is_symmetric;

				       /**
					* This determines the near
					* null space for the operator
					* Av, the vector-valued
					* velocity-velocity space. In
					* order to get good
					* performance, this vector has
					* to be specified using the
					* DoFTools::extract_constant_modes()
					* function.
					*/
	std::vector<std::vector<bool> > Av_constant_modes;

				       /**
					* Tolerance level that the
					* inner solvers on the
					* pressure mass matrix (and
					* the pressure Laplace matrix
					* in case of a Navier-Stokes
					* problem) should be solve
					* to. It is essential that
					* this tolerance is in the
					* same order of magnitude than
					* the one chosen for the outer
					* GMRES solver (or little
					* less). Otherwise unexpected
					* variations in the number of
					* iterations can occur from
					* solve to solve. 
					*
					* Default value: <tt>1e-9</tt>
					* residual relative to initial
					* residual.
					*/
	double inner_solve_tolerance;

				       /**
					* Determines whether the inner
					* solve for the pressure mass
					* matrix should use an IC
					* preconditioner (incomplete
					* Cholesky factorization) or
					* an SSOR preconditioner. The
					* former tends to be faster in
					* most situations, whereas the
					* latter uses almost no
					* additional memory besides
					* the matrix storage.
					*/
	bool use_ssor_on_mass_matrix;

				       /**
					* Determines whether the
					* underlying velocity
					* discretization uses linear
					* elements or higher order
					* elements, which has
					* consequences for the set up
					* of the internal setup of the
					* ML AMG preconditioner.
					*/
	bool velocity_uses_higher_order_elements;

				       /**
					* Determines whether the
					* underlying pressure
					* discretization uses linear
					* elements or higher order
					* elements, which has
					* consequences for the set up
					* of the internal setup of the
					* ML AMG preconditioner.
					*/
	bool pressure_uses_higher_order_elements;

                                       /**
					* Flag to determine whether we
					* should do an inner solve on
					* the velocity-velocity
					* matrix. By default, this
					* option is not on, since it
					* is not necessary for a good
					* performance of the outer
					* solver. However, it can be
					* beneficial to perform some
					* iterations on this block
					* only and not on the whole
					* block, so that the number of
					* outer iterations is
					* reduced. See the discussion
					* in the @ref step_22
					* "step-22" tutorial program.
					*/
	unsigned int Av_do_solve;
	  
				       /**
					* The number of iterations in
					* the inner solve.
					*/
	unsigned int Av_n_iters;
	  
				       /**
					* The residual tolerance that
					* is going to be used for the
					* inner solve on the
					* velocity-velocity matrix.
					*/
	double Av_tolerance;

				       /**
					* This defines whether
					* internal details about the
					* solve process should be
					* written to screen. This can
					* be a lot of information.
					*/
	bool output_details;
      };

				       /**
					* Lazy setup of the block
					* preconditioner for the
					* (Navier-) Stokes
					* system. This function takes
					* the matrices given here and
					* firs calculates good
					* preconditioners, i.e.,
					* algebraic multigrid
					* preconditioners for the
					* velocity-velocity matrix
					* <tt>Av</tt>, that can be
					* specified to be different
					* from the (0,0) block in the
					* system matrix, IC/SSOR on
					* the pressure mass matrix
					* <tt>Mp_matrix</tt>, and AMG
					* on the pressure Laplace
					* matrix
					* <tt>Lp_matrix</tt>. Next,
					* these preconditioners are
					* used to generate a block
					* operator.
					*/
      void initialize (const BlockSparseMatrix  &system_matrix,
		       const SparseMatrix       &Mp_matrix,
		       const AdditionalData     &additional_data,
		       const SparseMatrix       &Av_matrix = SparseMatrix(),
		       const SparseMatrix       &Lp_matrix = SparseMatrix(),
		       const SparseMatrix       &Fp_matrix = SparseMatrix());

				       /**
					* Set up the block
					* preconditioner for the
					* (Navier-) Stokes system. In
					* contrast to the other
					* <tt>initialize</tt>
					* function, this one expects
					* the user to specify
					* preconditioner objects for
					* the individual matrices,
					* that will then be used for
					* the inner inversion of the
					* matrices when building the
					* Schur complement
					* preconditioner. If no
					* preconditioner is specified,
					* the inner solvers will use
					* an identity preconditioner
					* object.
					*/
      void initialize (const BlockSparseMatrix  &system_matrix,
		       const SparseMatrix       &Mp_matrix,
		       const AdditionalData     &additional_data,
		       const PreconditionBase   &Mp_preconditioner = PreconditionBase(),
		       const PreconditionBase   &Av_preconditioner = PreconditionBase(),
		       const SparseMatrix       &Lp_matrix = SparseMatrix(),
		       const PreconditionBase   &Lp_preconditioner = PreconditionBase(),
		       const SparseMatrix       &Fp_matrix = SparseMatrix());

				       /**
					* This function can be used
				        * for a faster recalculation
				        * of the preconditioner
				        * construction when the system
				        * matrix entries underlying
				        * the preconditioner have
				        * changed but not the sparsity
				        * pattern (this means that the
				        * grid has remained the same
				        * and the matrix structure has
				        * not been
				        * reinitialized). Moreover, it
				        * is assumed that the PDE
				        * parameters have not changed
				        * drastically. This function
				        * is only needed for the lazy
				        * variant of the
				        * preconditioner. In the other
				        * case, the preconditioner can
				        * be modified outside this
				        * function, which will be
				        * recongnized in this
				        * preconditioner as well when
				        * doing a solve, without the
				        * need to restructure anything
				        * here.
					*/
      void reinit_lazy ();

                                       /**
					* Exception.
					*/
      DeclException1 (ExcNonMatchingMaps,
		      std::string,
		      << "The sparse matrix the preconditioner is based on "
		      << "uses a map that is not compatible to the one in vector"
		      << arg1
		      << ". Check preconditioner and matrix setup.");

    private:

				       /**
					* Pointer to preconditioner
					* for the mass matrix.
					*/
      Teuchos::RCP<PreconditionSSOR> Mp_precondition_ssor;

				       /**
					* Pointer to preconditioner
					* for the mass matrix.
					*/
      Teuchos::RCP<PreconditionIC>   Mp_precondition_ic;

				       /**
					* Pointer to preconditioner
					* for the velocity-velocity
					* matrix.
					*/
      Teuchos::RCP<PreconditionAMG> Av_precondition;

				       /**
					* Pointer to preconditioner
					* for the pressure Laplace
					* matrix.
					*/
      Teuchos::RCP<PreconditionAMG> Lp_precondition;
  };
}

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*--------------------------   trilinos_precondition_block.h   ---------------------------*/

#endif
/*--------------------------   trilinos_precondition_block.h   ---------------------------*/
