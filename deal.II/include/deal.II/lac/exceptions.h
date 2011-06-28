//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__lac_exceptions_h
#define __deal2__lac_exceptions_h

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace LACExceptions 
{
				   /**
				    * @addtogroup Exceptions
				    */
				   //@{

				   /**
				    * This function only works for
				    * quadratic matrices.
				    */
  DeclException0 (ExcNotQuadratic);

				   /**
				    * The operation cannot be finished
				    * since the matrix is singular.
				    */
  DeclException0 (ExcSingular);
  
				   /**
				    * Block indices of two block
				    * objects are different.
				    */
  DeclException0 (ExcDifferentBlockIndices);
  
				   /**
				    * An error of a PETSc function was
				    * encountered. Check the PETSc
				    * documentation for details.
				    */
  DeclException1 (ExcPETScError,
		  int,
		  << "An error with error number " << arg1
		  << " occured while calling a PETSc function");
  
				   /**
				    * An error of a Trilinos function was
				    * encountered. Check the Trilinos
				    * documentation for details.
				    */
  DeclException1 (ExcTrilinosError,
		  int,
		  << "An error with error number " << arg1
		  << " occured while calling a Trilinos function");
  
				   //@}
}


using namespace LACExceptions;


DEAL_II_NAMESPACE_CLOSE

#endif
