//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
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
 * @author Martin Kronbichler, 2008, 2009
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
}

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*--------------------------   trilinos_precondition_block.h   ---------------------------*/

#endif
/*--------------------------   trilinos_precondition_block.h   ---------------------------*/
