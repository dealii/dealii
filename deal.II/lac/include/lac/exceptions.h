//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__lac_exceptions_h
#define __deal2__lac_exceptions_h

#include <base/exceptions.h>

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
				    * Block indices of two block
				    * objects are different.
				    */
  DeclException0(ExcDifferentBlockIndices);
  
				   /**
				    * An error of a PETSc function was
				    * encountered. Check the PETSc
				    * documentation for details.
				    */
  DeclException1 (ExcPETScError,
		  int,
		  << "An error with error number " << arg1
		  << " occured while calling a PETSc function");
  
				   //@}
}


using namespace LACExceptions;


#endif
