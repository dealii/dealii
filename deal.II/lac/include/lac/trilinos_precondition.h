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
#ifndef __deal2__trilinos_precondition_h
#define __deal2__trilinos_precondition_h


#include <base/config.h>
#include <boost/shared_ptr.hpp>


#ifdef DEAL_II_USE_TRILINOS

// forward declarations
class Ifpack_Preconditioner;


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  class Vector;
  class SparseMatrix;
  

  class PreconditionSSOR
  {
    public:
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {
                                           /**
                                            * Constructor. By default,
                                            * set the damping parameter
                                            * to one.
                                            */
          AdditionalData (const double omega = 1);
	
                                           /**
                                            * Relaxation parameter.
                                            */
          double omega;
      };

                                       /**
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        */
      PreconditionSSOR (const SparseMatrix   &matrix,
                        const AdditionalData &additional_data = AdditionalData());
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (Vector       &dst,
		  const Vector &src) const;

    private:
      boost::shared_ptr<Ifpack_Preconditioner> preconditioner;
  };
  

}



DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_precondition_base.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition_base.h     ---------------------------*/
