//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/trilinos_precondition_block.h>

#ifdef DEAL_II_USE_TRILINOS

#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>

#include <Thyra_VectorSpaceImpl.hpp>
#include <Thyra_LinearOperatorDecl.hpp>
#include <Thyra_DefaultBlockedLinearOpDecl.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_InverseLinearOperator.hpp>
#include <Thyra_DefaultInverseLinearOp.hpp>
#include <Thyra_AztecOOLinearOpWithSolveFactory.hpp>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionBlockBase::PreconditionBlockBase()
  {}



  PreconditionBlockBase::~PreconditionBlockBase()
  {}



				       /**
					* Internal function that sets
					* up an inverse matrix.
					*/
  inline
  Thyra::ConstLinearOperator<double>
  inverse_matrix (const SparseMatrix    &input_M,
		  const Epetra_Operator *input_P,
		  const bool             is_symmetric,
		  const unsigned int     n_iterations,
		  const double           solve_tolerance,
		  const bool             output_details)
  {

					// Create discrete
					// Epetra operator
					// (Thyra wrapper around the
					// matrix).
    Teuchos::RCP<const Thyra::LinearOpBase<double> > tmpM = 
      Thyra::epetraLinearOp(Teuchos::rcp(&input_M.trilinos_matrix(), false));

    Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > M 
      = Teuchos::rcp(new Thyra::DefaultLinearOpSource<double>(tmpM), true);

    //Thyra::ConstLinearOperator<double> opM = tmpM;

					// Create a Thyra version of
					// the preconditioner.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > tmpP = 
      Thyra::epetraLinearOp(Teuchos::rcp(const_cast<Epetra_Operator*>(input_P),
					 false), 
			    Thyra::NOTRANS,
			    Thyra::EPETRA_OP_APPLY_APPLY_INVERSE,
			    Thyra::EPETRA_OP_ADJOINT_SUPPORTED);
	
    //Thyra::ConstLinearOperator<double> opPM = tmpP;

    Teuchos::RCP<Thyra::PreconditionerBase<double> > PM = 
      Teuchos::rcp (new Thyra::DefaultPreconditioner<double>(), true);

    {
      Thyra::DefaultPreconditioner<double>
	*defaultPrec = &Teuchos::dyn_cast<Thyra::DefaultPreconditioner<double> >
	               (*PM);

      (*defaultPrec).initializeUnspecified(tmpP);
    }

					// Set up the solver.
    Teuchos::RCP<Teuchos::ParameterList> aztecParams 
      = Teuchos::rcp(new Teuchos::ParameterList("aztecOOFSolverFactory"), true);

    aztecParams->sublist("Forward Solve").set("Max Iterations", 
					      (int)n_iterations);
    aztecParams->sublist("Forward Solve").set("Tolerance", 
					      solve_tolerance);

					// aztecOO solver settings
    if (is_symmetric == true)
      {
	aztecParams->sublist("Forward Solve")
	  .sublist("AztecOO Settings").set("Aztec Solver", "CG");
      }
    else
      {
	aztecParams->sublist("Forward Solve")
	  .sublist("AztecOO Settings").set("Aztec Solver", "GMRES");
	aztecParams->sublist("Forward Solve")
	  .sublist("AztecOO Settings").set("Size of Krylov Subspace", 
					   (int)n_iterations);
      }

    aztecParams->sublist("Forward Solve")
      .sublist("AztecOO Settings").set("Aztec Preconditioner", "none");
    aztecParams->sublist("Forward Solve")
      .sublist("AztecOO Settings").set("Convergence Test", "r0");

					// This sets if solver details
					// (the residual in each
					// iterations) should be
					// written to screen.
    if (output_details == true)
      {
	aztecParams->sublist("Forward Solve")
	  .sublist("AztecOO Settings").set("Output Frequency", 1);
      }
    else
      {
	aztecParams->sublist("VerboseObject").set("Verbosity Level", "none");
      }

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      aztecLowsFactory = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory(), true);

    aztecLowsFactory->setParameterList(aztecParams);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > rcpAztec 
      = aztecLowsFactory->createOp();

					// Does not work like
					// this. The preconditioner
					// will only appear in the
					// operator, but won't be set
					// correctly to the right hand
					// side vector, which makes
					// the results wrong. So try
					// something else...
    //Meros::LinearSolveStrategy<double> azM 
    //  = new Meros::AztecSolveStrategy(*(aztecParams.get()));
    //Thyra::LinearOperator<double> Minverse 
    //  = Thyra::inverse(*aztecLowsFactory, opM);
    //return new Meros::InverseOperator<double>(opM * opPM, azM);

					// Attach the preconditioner
					// to the solver.
    aztecLowsFactory->initializePreconditionedOp(M, PM, &*rcpAztec, 
						 Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
 
					// Create an inverse matrix
					// object and return.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > Minverse = 
      Teuchos::rcp(new Thyra::DefaultInverseLinearOp<double>(rcpAztec));

    return Minverse;
  }



  PreconditionStokes::AdditionalData::
  AdditionalData (const bool          right_preconditioning,
		  const bool          Av_is_symmetric,
		  const std::vector<std::vector<bool> > &Av_constant_modes,
		  const double        inner_solve_tolerance,
		  const bool          use_ssor_on_mass_matrix,
		  const bool          velocity_uses_higher_order_elements,
		  const bool          pressure_uses_higher_order_elements,
		  const bool          Av_do_solve,
		  const unsigned int  Av_n_iters,
		  const double        Av_tolerance,
		  const bool          output_details)
                  :
                  right_preconditioning (right_preconditioning),
                  Av_is_symmetric (Av_is_symmetric),
		  Av_constant_modes (Av_constant_modes),
		  inner_solve_tolerance (inner_solve_tolerance),
		  use_ssor_on_mass_matrix (use_ssor_on_mass_matrix),
		  velocity_uses_higher_order_elements (velocity_uses_higher_order_elements),
		  pressure_uses_higher_order_elements (pressure_uses_higher_order_elements),
		  Av_do_solve (Av_do_solve),
		  Av_n_iters (Av_n_iters),
		  Av_tolerance (Av_tolerance),
		  output_details (output_details)
  {}



  void
  PreconditionStokes::initialize(const BlockSparseMatrix  &system_matrix,
				 const SparseMatrix       &Mp_matrix,
				 const AdditionalData     &additional_data,
				 const SparseMatrix       &Av_matrix,
				 const SparseMatrix       &Lp_matrix,
				 const SparseMatrix       &Fp_matrix)
  {
    AssertThrow (system_matrix.n_block_rows() == 2,
		 ExcDimensionMismatch(system_matrix.n_block_rows(),
				      2));
    AssertThrow (system_matrix.n_block_cols() == 2,
		 ExcDimensionMismatch(system_matrix.n_block_cols(),
				      2));

					// Create Mp preconditioner.
    Teuchos::RCP<PreconditionBase> Mp_precondition;
    {
      if (additional_data.use_ssor_on_mass_matrix == true)
	{
	  Mp_precondition_ssor = Teuchos::rcp (new PreconditionSSOR(), true);
	  Mp_precondition_ssor->initialize (Mp_matrix,
					 PreconditionSSOR::AdditionalData (1.2,0,0));
	  Mp_precondition = Teuchos::rcp(&*Mp_precondition_ssor, false);
	}
      else
	{
	  Mp_precondition_ic = Teuchos::rcp (new PreconditionIC(), true);
	  Mp_precondition_ic->initialize (Mp_matrix);
	  Mp_precondition = Teuchos::rcp(&*Mp_precondition_ic, false);
	}
    }

					// Create Av AMG
					// preconditioner.
    {
      Av_precondition = Teuchos::rcp (new PreconditionAMG(), true);
      PreconditionAMG::AdditionalData av_amg_data;
      av_amg_data.elliptic = additional_data.Av_is_symmetric;
      av_amg_data.higher_order_elements = additional_data.velocity_uses_higher_order_elements;
      av_amg_data.constant_modes = additional_data.Av_constant_modes;

      if (Av_matrix.m() == system_matrix.block(0,0).m())
	Av_precondition->initialize (Av_matrix, av_amg_data);
      else
	Av_precondition->initialize (system_matrix.block(0,0), av_amg_data);
    }

					// Create Lp AMG
					// preconditioner.
    if (Lp_matrix.m() == system_matrix.block(1,1).m())
      {
	Lp_precondition = Teuchos::rcp (new PreconditionAMG(), true);
	PreconditionAMG::AdditionalData lp_amg_data;
	lp_amg_data.elliptic = true;
	lp_amg_data.higher_order_elements = additional_data.pressure_uses_higher_order_elements;

	Lp_precondition->initialize (Lp_matrix, lp_amg_data);
      }

    initialize (system_matrix, Mp_matrix, additional_data, *Mp_precondition,
		*Av_precondition, Lp_matrix, *Lp_precondition, Fp_matrix);

  }

  void 
  PreconditionStokes::initialize (const BlockSparseMatrix  &system_matrix,
				  const SparseMatrix       &Mp_matrix,
				  const AdditionalData     &additional_data,
				  const PreconditionBase   &Mp_preconditioner,
				  const PreconditionBase   &Av_preconditioner,
				  const SparseMatrix       &Lp_matrix,
				  const PreconditionBase   &Lp_preconditioner,
				  const SparseMatrix       &Fp_matrix)
  {
    AssertThrow (system_matrix.n_block_rows() == 2,
		 ExcDimensionMismatch(system_matrix.n_block_rows(),
				      2));
    AssertThrow (system_matrix.n_block_cols() == 2,
		 ExcDimensionMismatch(system_matrix.n_block_cols(),
				      2));

					// Create Av inverse matrix
    Thyra::ConstLinearOperator<double> Av_inverse;
    if (additional_data.Av_do_solve)
      Av_inverse = inverse_matrix(system_matrix.block(0,0),
				  &*Av_preconditioner.preconditioner,
				  additional_data.Av_is_symmetric,
				  additional_data.Av_n_iters,
				  additional_data.Av_tolerance,
				  additional_data.output_details);
    else
      {
	Teuchos::RCP<const Thyra::LinearOpBase<double> >
	  tmpAvp = Thyra::epetraLinearOp(Teuchos::rcp(&*Av_preconditioner.preconditioner,
						      false), 
					 Thyra::NOTRANS,
					 Thyra::EPETRA_OP_APPLY_APPLY_INVERSE,
					 Thyra::EPETRA_OP_ADJOINT_SUPPORTED);
	Av_inverse = tmpAvp;
      }


					// Create Mp inverse matrix.
    Thyra::ConstLinearOperator<double> Mp_inverse 
      = inverse_matrix(Mp_matrix, &*Mp_preconditioner.preconditioner, true, 
		       Mp_matrix.m(),
		       additional_data.inner_solve_tolerance,
		       additional_data.output_details);

					// Create Lp preconditioner
					// and inverse matrix in case
					// there is a matrix given.
    Thyra::ConstLinearOperator<double> Lp_inverse;

    if (Lp_matrix.m() == system_matrix.block(1,1).m())
      {
					// Create Lp inverse matrix.
	Lp_inverse = inverse_matrix(Lp_matrix, &*Lp_preconditioner.preconditioner,
				    true, Lp_matrix.m(), 
				    additional_data.inner_solve_tolerance,
				    additional_data.output_details);
      }
    else
      Lp_inverse = Thyra::identity<double> (Mp_inverse.range());

					// Create Fp Thyra operator.
    Thyra::ConstLinearOperator<double> Fp_op;
    if (Fp_matrix.m() != system_matrix.block(1,1).m())
      {
	Fp_op = Thyra::identity<double> (Mp_inverse.range());
      }
    else
      {
	Teuchos::RCP<const Thyra::LinearOpBase<double> >
	  tmpFp = Thyra::epetraLinearOp(Teuchos::rcp(&Fp_matrix.trilinos_matrix(),
						     false));
	Fp_op = tmpFp;
      }

					// Build the PCD block
					// preconditioner factory.
					// Build identity matrices on
					// the velocity and pressure
					// spaces
    Thyra::ConstLinearOperator<double> Ivel 
      = Thyra::identity<double>(Av_inverse.range());
    Thyra::ConstLinearOperator<double> Ipress 
      = Thyra::identity<double>(Mp_inverse.domain());

					// Build zero operators. Need
					// one that is pressure x
					// velocity and one that is
					// velocity x pressure
    Thyra::ConstLinearOperator<double> zero;
    
					// Build the composed Schur
					// complement approximation
					// inverse
    Thyra::ConstLinearOperator<double> Xinv;
    if (additional_data.right_preconditioning == true)
      Xinv = Mp_inverse * Fp_op * Lp_inverse;
    else
      Xinv = Lp_inverse * Fp_op * Mp_inverse;

					// Create discrete pressure
					// gradient and divergence
					// operators in Thyra format.
    Teuchos::RCP<const Thyra::LinearOpBase<double> > tmpS01 = 
      Thyra::epetraLinearOp(Teuchos::rcp(&system_matrix.block(0,1).trilinos_matrix(), 
					 false));
    Thyra::ConstLinearOperator<double> S01 = tmpS01;

    Teuchos::RCP<const Thyra::LinearOpBase<double> > tmpS10 = 
      Thyra::epetraLinearOp(Teuchos::rcp(&system_matrix.block(1,0).trilinos_matrix(), 
					 false));
    Thyra::ConstLinearOperator<double> S10 = tmpS10;

					// Build the 3 block operators
					// for the
					// preconditioner. Here we
					// have to set different
					// matrices depending on
					// whether we do right or left
					// preconditioning.
    Thyra::ConstLinearOperator<double> P1 
      = block2x2( Av_inverse, zero, zero, Ipress );

    Thyra::ConstLinearOperator<double> P2;
    if (additional_data.right_preconditioning == true)
      P2 = block2x2( Ivel, (-1.0)*S01, zero, Ipress  );
    else
      P2 = block2x2( Ivel, zero, (-1.0)*S10, Ipress );

    Thyra::ConstLinearOperator<double> P3 
      = block2x2( Ivel,   zero, zero, (-1.0)*Xinv );

    if (additional_data.right_preconditioning == true)    
      preconditioner = Teuchos::rcp (new 
		       Thyra::ConstLinearOperator<double>(P1 * P2 * P3));
    else
      preconditioner = Teuchos::rcp (new 
		       Thyra::ConstLinearOperator<double>(P3 * P2 * P1));

  }



  void
  PreconditionStokes::reinit_lazy ()
  {
    Assert (&*Av_precondition != 0,
	    ExcMessage ("Can only be done when the outer preconditioner owns"
			" the AMG object!"));

    Av_precondition->reinit();
  }

} /* end of namespace TrilinosWrappers */

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
