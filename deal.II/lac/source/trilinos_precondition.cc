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
#include <lac/trilinos_vector_base.h>
#include <lac/trilinos_sparse_matrix.h>

#ifdef DEAL_II_USE_TRILINOS

#include <Ifpack.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionJacobi::AdditionalData::
  AdditionalData (const double omega,
		  const double min_diagonal)
                  :
                  omega (omega),
		  min_diagonal (min_diagonal)
  {}

  

  PreconditionJacobi::PreconditionJacobi ()
  {}



  void
  PreconditionJacobi::initialize (const SparseMatrix   &matrix,
				  const AdditionalData &additional_data)
  {

    preconditioner = Teuchos::rcp (Ifpack().Create ("point relaxation", &*matrix.matrix, 0));
    Assert (&*preconditioner != 0, ExcMessage ("Trilinos could not create this "
					       "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: type", "Jacobi");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value", 
			additional_data.min_diagonal);
    
    ierr = preconditioner->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  PreconditionJacobi::vmult (VectorBase       &dst,
			     const VectorBase &src) const
  {
    const int ierr = preconditioner->ApplyInverse (*src.vector, *dst.vector);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  PreconditionSSOR::AdditionalData::
  AdditionalData (const double       omega,
		  const double       min_diagonal,
		  const unsigned int overlap)
                  :
                  omega        (omega),
		  min_diagonal (min_diagonal),
		  overlap      (overlap)
  {}

  

  PreconditionSSOR::PreconditionSSOR ()
  {}



  void
  PreconditionSSOR::initialize (const SparseMatrix   &matrix,
				const AdditionalData &additional_data)
  {

    preconditioner = Teuchos::rcp (Ifpack().Create ("point relaxation",
						    &*matrix.matrix, 
						    additional_data.overlap));
    Assert (&*preconditioner != 0, ExcMessage ("Trilinos could not create this "
					       "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value", 
			additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");
    
    ierr = preconditioner->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  PreconditionSSOR::vmult (VectorBase       &dst,
			   const VectorBase &src) const
  {
    const int ierr = preconditioner->ApplyInverse (*src.vector, *dst.vector);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  PreconditionSOR::AdditionalData::
  AdditionalData (const double       omega,
		  const double       min_diagonal,
		  const unsigned int overlap)
                  :
                  omega        (omega),
		  min_diagonal (min_diagonal),
		  overlap      (overlap)
  {}

  

  PreconditionSOR::PreconditionSOR ()
  {}



  void
  PreconditionSOR::initialize (const SparseMatrix   &matrix,
			       const AdditionalData &additional_data)
  {

    preconditioner = Teuchos::rcp (Ifpack().Create ("point relaxation",
						    &*matrix.matrix, 
						    additional_data.overlap));
    Assert (&*preconditioner != 0, ExcMessage ("Trilinos could not create this "
					       "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: type", "Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value", 
			additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");
    
    ierr = preconditioner->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  PreconditionSOR::vmult (VectorBase       &dst,
			  const VectorBase &src) const
  {
    const int ierr = preconditioner->ApplyInverse (*src.vector, *dst.vector);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  PreconditionIC::AdditionalData::
  AdditionalData (const double       ic_fill,
		  const double       ic_atol,
		  const double       ic_rtol,
		  const unsigned int overlap)
                  :
                  ic_fill (ic_fill),
		  ic_atol (ic_atol),
		  ic_rtol (ic_rtol),
		  overlap (overlap)
  {}



  PreconditionIC::PreconditionIC ()
  {}



  void
  PreconditionIC::initialize (const SparseMatrix   &matrix,
			      const AdditionalData &additional_data)
  {
    preconditioner = Teuchos::rcp (Ifpack().Create ("IC", &*matrix.matrix, 
						    additional_data.overlap));
    Assert (&*preconditioner != 0, ExcMessage ("Trilinos could not create this "
					       "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: level-of-fill",additional_data.ic_fill); 
    parameter_list.set ("fact: absolute threshold",additional_data.ic_atol); 
    parameter_list.set ("fact: relative threshold",additional_data.ic_rtol); 
    parameter_list.set ("schwarz: combine mode", "Add");
    
    ierr = preconditioner->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  PreconditionIC::vmult (VectorBase       &dst,
			 const VectorBase &src) const
  {
    const int ierr = preconditioner->ApplyInverse (*src.vector, *dst.vector);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  PreconditionILU::AdditionalData::
  AdditionalData (const double       ilu_fill,
		  const double       ilu_atol,
		  const double       ilu_rtol,
		  const unsigned int overlap)
                  :
                  ilu_fill (ilu_fill),
		  ilu_atol (ilu_atol),
		  ilu_rtol (ilu_rtol),
		  overlap  (overlap)
  {}



  PreconditionILU::PreconditionILU ()
  {}



  void
  PreconditionILU::initialize (const SparseMatrix   &matrix,
			       const AdditionalData &additional_data)
  {
    preconditioner = Teuchos::rcp (Ifpack().Create ("ILU", &*matrix.matrix, 
						    additional_data.overlap));
    Assert (&*preconditioner != 0, ExcMessage ("Trilinos could not create this "
					       "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: level-of-fill",additional_data.ilu_fill); 
    parameter_list.set ("fact: absolute threshold",additional_data.ilu_atol); 
    parameter_list.set ("fact: relative threshold",additional_data.ilu_rtol); 
    parameter_list.set ("schwarz: combine mode", "Add");
    
    ierr = preconditioner->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = preconditioner->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  PreconditionILU::vmult (VectorBase       &dst,
			  const VectorBase &src) const
  {
    const int ierr = preconditioner->ApplyInverse (*src.vector, *dst.vector);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
