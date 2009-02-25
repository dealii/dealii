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

#include <lac/trilinos_solver_block.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_precondition_block.h>

#include <base/utilities.h>

#include <cmath>

#ifdef DEAL_II_USE_TRILINOS

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>

#include <Thyra_VectorDecl.hpp>
#include <Thyra_VectorImpl.hpp>
#include <Thyra_VectorSpaceImpl.hpp>
#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_MultiVectorDefaultBase.hpp>
#include <Thyra_LinearOperatorDecl.hpp>
#include <Thyra_LinearOperatorImpl.hpp>
#include <Thyra_DefaultBlockedLinearOpDecl.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_AztecOOLinearOpWithSolveFactory.hpp>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  SolverBlockBase::AdditionalData::
  AdditionalData (const unsigned int gmres_restart_parameter,
		  const bool         right_preconditioning,
		  const bool         output_details)
                  :
                  gmres_restart_parameter (gmres_restart_parameter),
		  right_preconditioning   (right_preconditioning),
		  output_details          (output_details)
  {}

  
  
  SolverBlockBase::SolverBlockBase (SolverControl  &cn)
                  :
                  solver_name    (gmres),
                  solver_control (cn)
  {}

  

  SolverBlockBase::SolverBlockBase (SolverBlockBase::SolverBlockName solver_name,
				    SolverControl                   &cn)
                  :
                  solver_name    (solver_name),
                  solver_control (cn)
  {}

  

  SolverBlockBase::~SolverBlockBase ()
  {}

  

  SolverControl &
  SolverBlockBase::control() const
  {
    return solver_control;
  }



  void
  SolverBlockBase::solve (const BlockSparseMatrix     &A,
			  BlockVector                 &x,
			  const BlockVector           &b,
			  const PreconditionBlockBase &preconditioner)
  {

    Assert (x.n_blocks() == b.n_blocks(),
	    ExcDimensionMismatch(x.n_blocks(), b.n_blocks()));

					// Just copy everything over
					// to a distributed vector...
    MPI::BlockVector tmpx (x.n_blocks());
    MPI::BlockVector tmpb (b.n_blocks());
    for (unsigned int i=0; i<x.n_blocks(); ++i)
      {
	tmpx.block(i).reinit (A.block(i,i).domain_partitioner(), false);
	tmpx.block(i) = x.block(i);
	tmpb.block(i).reinit (A.block(i,i).range_partitioner(), false);
	tmpb.block(i) = b.block(i);
      }
    tmpx.collect_sizes();
    tmpb.collect_sizes();

    solve (A, tmpx, tmpb, preconditioner);
    x = tmpx;
  }
  
    

  void
  SolverBlockBase::solve (const BlockSparseMatrix     &input_A,
			  MPI::BlockVector            &input_x,
			  const MPI::BlockVector      &input_b,
			  const PreconditionBlockBase &preconditioner)
  {

					// Get block information about
					// the preconditioner.
    Thyra::ConstLinearOperator<double> P = *preconditioner.preconditioner;
    const unsigned int n_rows = P.numBlockRows();
    const unsigned int n_cols = P.numBlockCols();

					// Check the dimensions of the
					// block vectors and matrix
					// and check if the parallel
					// distribution of the
					// matrices and vectors is the
					// same.
    AssertThrow (n_rows == n_cols,
		 ExcDimensionMismatch (n_rows, n_cols));

    AssertThrow (input_x.n_blocks() == input_A.n_block_cols(),
		 ExcDimensionMismatch(input_x.n_blocks(),
				      input_A.n_block_cols()));
    AssertThrow (input_b.n_blocks() == input_A.n_block_rows(),
		 ExcDimensionMismatch(input_b.n_blocks(),
				      input_A.n_block_rows()));
    AssertThrow (input_A.n_block_rows() == n_rows,
		 ExcDimensionMismatch(input_A.n_block_rows(), n_rows));
    AssertThrow (input_A.n_block_cols() == n_cols,
		 ExcDimensionMismatch(input_A.n_block_cols(), n_cols));

    for (unsigned int i=0; i<n_rows; ++i)
      {
	AssertThrow (input_x.block(i).vector_partitioner().UniqueGIDs() == true,
		     ExcOverlappingMaps("vector", "x"));
	AssertThrow (input_b.block(i).vector_partitioner().UniqueGIDs() == true,
		     ExcOverlappingMaps("vector", "b"));
	for (unsigned int j=0; j<n_cols; ++j)
	  {
	    {
	      std::ostringstream i_str;
	      i_str << i;
	      std::ostringstream error_component;
	      error_component << "x.block(" << i_str.str() << ")";
	      AssertThrow (input_x.block(j).vector_partitioner().SameAs(
			     input_A.block(i,j).domain_partitioner()) == true,
			   ExcNonMatchingMaps (error_component.str()));
	    }
	    {
	      std::ostringstream i_str;
	      i_str << j;
	      std::ostringstream error_component;
	      error_component << "b.block(" << i_str.str() << ")";
	      AssertThrow (input_b.block(i).vector_partitioner().SameAs(
			     input_A.block(i,j).range_partitioner()) == true,
			   ExcNonMatchingMaps (error_component.str()));
	    }
	  }
      }


					// Wrap Epetra vectors into
					// Thyra vectors in order to
					// be able to work with the
					// abstract Thyra/Stratimikos
					// interfaces to the AztecOO
					// solver.
					// 
					// As a first step, create a
					// copy of the vector spaces.
    std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > 
      epetra_vector_spaces;
    for (unsigned int i=0; i<n_rows; ++i)
      {
	Teuchos::RCP<const Thyra::VectorSpaceBase<double> > tmp_space
	  = Thyra::create_VectorSpace(
				      Teuchos::rcp(&input_A.block(i,i).domain_partitioner(), 
						   false));

	epetra_vector_spaces.push_back(tmp_space);
      }

					// Convert the block matrix to
					// a Thyra linear operator
					// that acts on blocks.
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> >
      tmpA = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<double>());
    tmpA->beginBlockFill(n_rows,n_cols);
    for (unsigned int i=0; i<n_rows; ++i)
      for (unsigned int j=0; j<n_cols; ++j)
	{
	  Teuchos::RCP<const Thyra::LinearOpBase<double> > 
	    A_ij = Thyra::epetraLinearOp(Teuchos::rcp(&input_A.block(i,j).trilinos_matrix(), 
						      false));
	  tmpA->setBlock(i, j, A_ij);
	}
    tmpA->endBlockFill();
    tmpA->setObjectLabel("Blocked system matrix");
    Thyra::RCP<const Thyra::LinearOpBase<double> > baseA = tmpA;

    Thyra::ConstLinearOperator<double> A = baseA;

					// Convert block vectors for
					// right hand side and
					// solution to Thyra objects.
    Thyra::Vector<double> rhs = A.range().createMember();
    Thyra::Vector<double> sol = A.domain().createMember();
    for (unsigned int i=0; i<n_rows; ++i)
      {
	Epetra_Vector *block_ptr = const_cast<Epetra_Vector*>
	  ((input_b.block(i).trilinos_vector())(0));
	Teuchos::RCP<Thyra::VectorBase<double> > tmp_rhs_i
	  = Thyra::create_Vector(Teuchos::rcp(block_ptr, false), 
				 epetra_vector_spaces[i]);
	const Thyra::Vector<double> rhs_i = tmp_rhs_i;
	rhs.setBlock(i, rhs_i);

	block_ptr = (input_x.block(i).trilinos_vector())(0);
	Teuchos::RCP<Thyra::VectorBase<double> > tmp_sol_i
	  = Thyra::create_Vector(Teuchos::rcp(block_ptr, false), 
				 epetra_vector_spaces[i]);
	Thyra::Vector<double> sol_i = tmp_sol_i;
	sol.setBlock(i, sol_i);
      }

    // ------- Now build a solver factory for the block problem ------- 

					// Set up a parameter list for
					// the AztecOO solver.
    Teuchos::RCP<Teuchos::ParameterList> aztecBlockParams 
      = rcp(new Teuchos::ParameterList("AztecOO block solver parameters"));
      
					// Set the number of
					// iterations and the solver
					// tolerance from the solver
					// control object.
    aztecBlockParams->sublist("Forward Solve").set("Max Iterations", 
						   (int)solver_control.max_steps());
    aztecBlockParams->sublist("Forward Solve").set("Tolerance",
						    solver_control.tolerance());
    aztecBlockParams->sublist("Forward Solve")
      .sublist("AztecOO Settings").set("Convergence Test", "no scaling");

					// The settings about the
					// actual solver that is going
					// to be used - let this be
					// determined by the
					// <tt>solver_name</tt> set in
					// a derived class.
    switch (solver_name)
      {
        case cg:
	  aztecBlockParams->sublist("Forward Solve")
	    .sublist("AztecOO Settings").set("Aztec Solver", "CG");
	  break;
        case gmres:
	  aztecBlockParams->sublist("Forward Solve")
		    .sublist("AztecOO Settings").set("Aztec Solver", "GMRES");
	  aztecBlockParams->sublist("Forward Solve")
	    .sublist("AztecOO Settings").set("Size of Krylov Subspace", 
					     (int)additional_data.gmres_restart_parameter);
	  break;
        default:
	  ExcNotImplemented();
      }

					// We use an externally
					// provided preconditioner, so
					// do not set anything here.
    aztecBlockParams->sublist("Forward Solve")
      .sublist("AztecOO Settings").set("Aztec Preconditioner", "none");

    if (additional_data.output_details == true)
      aztecBlockParams->sublist("Forward Solve")
	.sublist("AztecOO Settings").set("Output Frequency", 1);
    else
      aztecBlockParams->sublist("VerboseObject").set("Verbosity Level", 
						     "none");

					// Create solve factory from
					// parameter list.
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      aztecBlockLowsFactory = Teuchos::rcp
        (new Thyra::AztecOOLinearOpWithSolveFactory());

    aztecBlockLowsFactory->setParameterList(aztecBlockParams);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > rcpBlockAztec 
      = aztecBlockLowsFactory->createOp();


    // ---- Set up the preconditioned inverse object and do the solve ----

					// Extract base operators from
					// block matrix A and block
					// preconditioner P.
    Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > A_ptr 
      = Teuchos::rcp(new Thyra::DefaultLinearOpSource<double>(A.constPtr()), 
		     true);
    Teuchos::RCP<Thyra::PreconditionerBase<double> > P_ptr = 
      Teuchos::rcp (new Thyra::DefaultPreconditioner<double>(), true);
    {
      Thyra::DefaultPreconditioner<double>
	*defaultPrec = &Teuchos::dyn_cast<Thyra::DefaultPreconditioner<double> >
	              (*P_ptr);
      if (additional_data.right_preconditioning == true)
	(*defaultPrec).initializeRight(P.constPtr());
      else
	(*defaultPrec).initializeUnspecified(P.constPtr());
    }

					// Attach the matrix and the
					// preconditioner to the
					// solver object.
    aztecBlockLowsFactory->initializePreconditionedOp(A_ptr, P_ptr, 
						      &*rcpBlockAztec, 
						      Thyra::SUPPORT_SOLVE_FORWARD_ONLY);

					// Now do the solve.
    Thyra::SolveStatus<double>
      solveStatus = Thyra::solve( *rcpBlockAztec, Thyra::NOTRANS, 
				  *rhs.constPtr().get(), &*sol.ptr().get() );

					// Extract the number of
					// needed iterations from the
					// solver message.
    unsigned int n_iterations;
    {
      std::string output_text = solveStatus.message;
      unsigned int position = output_text.find ("using ");
      const std::pair<int,unsigned int> tmp 
	= Utilities::get_integer_at_position (output_text, position+6);
      n_iterations = tmp.first;
    }

					// Do a dummy check of
					// convergence to tell the
					// SolverControl object the
					// number of iterations. We do
					// not impose the achieved
					// tolerance here since that
					// would require additional
					// work (computing a
					// residual).
    solver_control.check (n_iterations, 0);

  }
  
    


  
/* ------------------------ SolverCG -------------------------- */

  SolverBlockCG::AdditionalData::
  AdditionalData (const bool         right_preconditioning,
		  const bool         output_details)
                  :
                  right_preconditioning   (right_preconditioning),
		  output_details          (output_details)
  {}

  
  
  SolverBlockCG::SolverBlockCG (SolverControl        &cn,
				const AdditionalData &data)
                  :
                  SolverBlockBase (cn),
                  data (data)
  {
    solver_name = cg;
    additional_data.right_preconditioning = data.right_preconditioning;
    additional_data.output_details = data.output_details;
  }
  

/* ---------------------- SolverGMRES ------------------------- */

  SolverBlockGMRES::AdditionalData::
  AdditionalData (const unsigned int gmres_restart_parameter,
		  const bool         right_preconditioning,
		  const bool         output_details)
                  :
                  gmres_restart_parameter (gmres_restart_parameter),
		  right_preconditioning   (right_preconditioning),
		  output_details          (output_details)
  {}

  

  SolverBlockGMRES::SolverBlockGMRES (SolverControl        &cn,
				      const AdditionalData &data)
                  :
                  SolverBlockBase (cn),
                  data (data)
  {
    solver_name = gmres;
    additional_data.gmres_restart_parameter = data.gmres_restart_parameter;
    additional_data.right_preconditioning = data.right_preconditioning;
    additional_data.output_details = data.output_details;
  }


}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
