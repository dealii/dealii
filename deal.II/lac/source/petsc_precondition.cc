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


#include <lac/petsc_matrix_base.h>
#include <lac/petsc_vector_base.h>
#include <lac/petsc_precondition.h>

#include <cmath>

#ifdef DEAL_II_USE_PETSC


namespace PETScWrappers
{
  PreconditionerBase::PreconditionerBase (const MatrixBase &matrix)
                  :
                  matrix (matrix)
  {}


  
  PreconditionerBase::~PreconditionerBase ()
  {}

  

  PreconditionerBase::operator const Mat () const
  {
    return matrix;
  }


/* ----------------- PreconditionJacobi -------------------- */


  PreconditionJacobi::PreconditionJacobi (const MatrixBase     &matrix,
                                          const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionJacobi::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCJACOBI));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionJacobi -------------------- */


  PreconditionBlockJacobi::
  PreconditionBlockJacobi (const MatrixBase     &matrix,
                           const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionBlockJacobi::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCBJACOBI));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionSOR -------------------- */

  PreconditionSOR::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}

  
  
  PreconditionSOR::PreconditionSOR (const MatrixBase     &matrix,
                                    const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionSOR::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCSOR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags as given
    ierr = PCSORSetOmega (pc, additional_data.omega);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionSSOR -------------------- */

  PreconditionSSOR::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}

  
  
  PreconditionSSOR::PreconditionSSOR (const MatrixBase     &matrix,
                                      const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionSSOR::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCSOR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags as given
    ierr = PCSORSetOmega (pc, additional_data.omega);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // convert SOR to SSOR
    ierr = PCSORSetSymmetric (pc, SOR_SYMMETRIC_SWEEP);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionEisenstat -------------------- */

  PreconditionEisenstat::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}

  
  
  PreconditionEisenstat::PreconditionEisenstat (const MatrixBase     &matrix,
                                                const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionEisenstat::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCEISENSTAT));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags as given
    ierr = PCEisenstatSetOmega (pc, additional_data.omega);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionICC -------------------- */


  PreconditionICC::AdditionalData::
  AdditionalData (const unsigned int levels)
                  :
                  levels (levels)
  {}

  
  
  PreconditionICC::PreconditionICC (const MatrixBase     &matrix,
                                    const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionICC::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCICC));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags
    PCICCSetLevels (pc, additional_data.levels);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionILU -------------------- */

  PreconditionILU::AdditionalData::
  AdditionalData (const unsigned int levels)
                  :
                  levels (levels)
  {}

  
  
  PreconditionILU::PreconditionILU (const MatrixBase     &matrix,
                                    const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionILU::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCILU));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags
    PCILUSetLevels (pc, additional_data.levels);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionLU -------------------- */

  PreconditionLU::AdditionalData::
  AdditionalData (const double pivoting,
		  const double zero_pivot,
		  const double damping)
                  :
                  pivoting (pivoting),
                  zero_pivot (zero_pivot),
		  damping (damping)
  {}

  
  
  PreconditionLU::PreconditionLU (const MatrixBase     &matrix,
				  const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionLU::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCLU));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // set flags as given
    ierr = PCLUSetPivoting (pc, additional_data.pivoting);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCLUSetZeroPivot (pc, additional_data.zero_pivot);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCLUSetDamping (pc, additional_data.damping);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
