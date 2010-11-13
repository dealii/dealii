//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2006, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_precondition.h>

#ifdef DEAL_II_USE_PETSC

#  include <base/utilities.h>
#  include <lac/petsc_matrix_base.h>
#  include <lac/petsc_vector_base.h>
#  include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  PreconditionerBase::PreconditionerBase (const MatrixBase &matrix)
                  :
                  matrix (matrix)
  {}


  
  PreconditionerBase::~PreconditionerBase ()
  {}

  

  PreconditionerBase::operator Mat () const
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

    ierr = PCSetFromOptions (pc);
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

    ierr = PCSetFromOptions (pc);
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

    ierr = PCSetFromOptions (pc);
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

    ierr = PCSetFromOptions (pc);
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

    ierr = PCSetFromOptions (pc);
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
#if DEAL_II_PETSC_VERSION_LT(2,3,1)
    PCICCSetLevels (pc, additional_data.levels);
#else
    PCFactorSetLevels (pc, additional_data.levels);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
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
#if DEAL_II_PETSC_VERSION_LT(2,3,1)
    PCILUSetLevels (pc, additional_data.levels);  
#else
    PCFactorSetLevels (pc, additional_data.levels);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ----------------- PreconditionBoomerAMG -------------------- */

  PreconditionBoomerAMG::AdditionalData::
  AdditionalData(const bool symmetric_operator,
		 const double strong_threshold,
		 const double max_row_sum,
		 const unsigned int aggressive_coarsening_num_levels,
		 const bool output_details
  )
		  :
		  symmetric_operator(symmetric_operator),
		  strong_threshold(strong_threshold),
		  max_row_sum(max_row_sum),
		  aggressive_coarsening_num_levels(aggressive_coarsening_num_levels),
		  output_details(output_details)
  {}

  
  PreconditionBoomerAMG::PreconditionBoomerAMG (const MatrixBase     &matrix,
						const AdditionalData &additional_data)
                  :
                  PreconditionerBase (matrix),
                  additional_data (additional_data)
  {}

  
  void
  PreconditionBoomerAMG::set_preconditioner_type (PC &pc) const
  {
                                     // set the right type for the
                                     // preconditioner
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCHYPRE));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    ierr = PCHYPRESetType(pc, "boomeramg");
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    if (additional_data.output_details)
      PetscOptionsSetValue("-pc_hypre_boomeramg_print_statistics","1");
    
    PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl",
			 Utilities::int_to_string(
			   additional_data.aggressive_coarsening_num_levels
			 ).c_str());

    std::stringstream ssStream;
    ssStream << additional_data.max_row_sum;
    PetscOptionsSetValue("-pc_hypre_boomeramg_max_row_sum", ssStream.str().c_str());
    
    ssStream.str(""); // empty the stringstream
    ssStream << additional_data.strong_threshold;    
    PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", ssStream.str().c_str());

    if (additional_data.symmetric_operator)
      {
	PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_up", "symmetric-SOR/Jacobi");
	PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_down", "symmetric-SOR/Jacobi");
	PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_coarse", "Gaussian-elimination");
      }
    else
      {
	PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_up", "SOR/Jacobi");
	PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_down", "SOR/Jacobi");
	PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_coarse", "Gaussian-elimination");
      }
    
    ierr = PCSetFromOptions (pc);
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
#if DEAL_II_PETSC_VERSION_LT(2,3,1)
    ierr = PCLUSetPivoting (pc, additional_data.pivoting);
#elif DEAL_II_PETSC_VERSION_LT(3,0,1)
    ierr = PCFactorSetPivoting (pc, additional_data.pivoting);
#else
    ierr = PCFactorSetColumnPivot (pc, additional_data.pivoting);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(2,3,0)
    ierr = PCLUSetZeroPivot (pc, additional_data.zero_pivot);
#else
    ierr = PCFactorSetZeroPivot (pc, additional_data.zero_pivot);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(2,3,0)
    ierr = PCLUSetDamping (pc, additional_data.damping);
#elif DEAL_II_PETSC_VERSION_LT(3,0,1)
    ierr = PCFactorSetShiftNonzero (pc, additional_data.damping);
#else
    ierr = PCFactorSetShiftAmount (pc, additional_data.damping);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
