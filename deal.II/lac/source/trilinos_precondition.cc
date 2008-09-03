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


#include <lac/trilinos_precondition.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_sparse_matrix.h>

#include <Ifpack.h>



#ifdef DEAL_II_USE_TRILINOS


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionSSOR::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}

  
  PreconditionSSOR::PreconditionSSOR (const SparseMatrix   &matrix,
				      const AdditionalData &additional_data)
		  :
		  preconditioner (Ifpack().Create ("point relaxation",
						   &*matrix.matrix, 0))
  {
    Assert (preconditioner != 0, ExcMessage ("Trilinos could not create this preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    
    ierr = preconditioner->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }


  void
  PreconditionSSOR::vmult (Vector       &dst,
			   const Vector &src) const
  {
    preconditioner->ApplyInverse (*src.vector, *dst.vector);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
