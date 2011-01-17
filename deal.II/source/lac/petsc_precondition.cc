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
#  include <lac/petsc_solver.h>
#  include <petscconf.h>
#  include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  PreconditionerBase::PreconditionerBase ()
		  :
		  pc(NULL), matrix(NULL)
  {}


  PreconditionerBase::~PreconditionerBase ()
  {
    if (pc!=NULL)
      {
	int ierr = PCDestroy(pc);
	AssertThrow (ierr == 0, ExcPETScError(ierr));
      }
  }
  
  
  void
  PreconditionerBase::vmult (VectorBase       &dst,
			     const VectorBase &src) const
  {
    AssertThrow (pc != NULL, StandardExceptions::ExcInvalidState ());

    int ierr;
    ierr = PCApply(pc, src, dst);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }  

  
  void
  PreconditionerBase::create_pc ()
  {
				     // only allow the creation of the
				     // preconditioner once
    AssertThrow (pc == NULL, StandardExceptions::ExcInvalidState ());
    
    MPI_Comm comm;
    int ierr;
				     // this ugly cast is necessay because the
				     // type Mat and PETScObject are
				     // unrelated.
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    ierr = PCCreate(comm, &pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetOperators(pc , matrix, matrix, SAME_PRECONDITIONER);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


  const PC &
  PreconditionerBase::get_pc () const
  {
    return pc;
  }

  
  PreconditionerBase::operator Mat () const
  {
    return matrix;
  }


/* ----------------- PreconditionJacobi -------------------- */

  PreconditionJacobi::PreconditionJacobi ()
  {}
  
    
  PreconditionJacobi::PreconditionJacobi (const MatrixBase     &matrix,
                                          const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionJacobi::initialize (const MatrixBase     &matrix_,
				  const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();
    
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCJACOBI));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


/* ----------------- PreconditionBlockJacobi -------------------- */

  PreconditionBlockJacobi::PreconditionBlockJacobi ()
  {}
  
    
  PreconditionBlockJacobi::
  PreconditionBlockJacobi (const MatrixBase     &matrix,
                           const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }

  void
  PreconditionBlockJacobi::initialize (const MatrixBase     &matrix_,
				       const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCBJACOBI));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


/* ----------------- PreconditionSOR -------------------- */

  PreconditionSOR::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}


  PreconditionSOR::PreconditionSOR ()
  {}
  

  PreconditionSOR::PreconditionSOR (const MatrixBase     &matrix,
                                    const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionSOR::initialize (const MatrixBase     &matrix_,
			       const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCSOR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags as given
    ierr = PCSORSetOmega (pc, additional_data.omega);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


/* ----------------- PreconditionSSOR -------------------- */

  PreconditionSSOR::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}

  
  PreconditionSSOR::PreconditionSSOR ()
  {}


  PreconditionSSOR::PreconditionSSOR (const MatrixBase     &matrix,
                                      const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionSSOR::initialize (const MatrixBase     &matrix_,
				const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

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

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


/* ----------------- PreconditionEisenstat -------------------- */

  PreconditionEisenstat::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}


  PreconditionEisenstat::PreconditionEisenstat ()
  {}
  

  PreconditionEisenstat::PreconditionEisenstat (const MatrixBase     &matrix,
                                                const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionEisenstat::initialize (const MatrixBase     &matrix_,
				     const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCEISENSTAT));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // then set flags as given
    ierr = PCEisenstatSetOmega (pc, additional_data.omega);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


/* ----------------- PreconditionICC -------------------- */


  PreconditionICC::AdditionalData::
  AdditionalData (const unsigned int levels)
                  :
                  levels (levels)
  {}


  PreconditionICC::PreconditionICC ()
  {}

  
  PreconditionICC::PreconditionICC (const MatrixBase     &matrix,
                                    const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionICC::initialize (const MatrixBase     &matrix_,
			       const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

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

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


/* ----------------- PreconditionILU -------------------- */

  PreconditionILU::AdditionalData::
  AdditionalData (const unsigned int levels)
                  :
                  levels (levels)
  {}


  PreconditionILU::PreconditionILU ()
  {}
  

  PreconditionILU::PreconditionILU (const MatrixBase     &matrix,
                                    const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionILU::initialize (const MatrixBase     &matrix_,
			       const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

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

    ierr = PCSetUp (pc);
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


  PreconditionBoomerAMG::PreconditionBoomerAMG ()
  {}

  
  PreconditionBoomerAMG::PreconditionBoomerAMG (const MatrixBase     &matrix,
						const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionBoomerAMG::initialize (const MatrixBase     &matrix_,
				     const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;

#ifdef PETSC_HAVE_HYPRE
    create_pc();

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

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#else // PETSC_HAVE_HYPRE
    (void)pc;
    Assert (false,
	    ExcMessage ("Your PETSc installation does not include a copy of "
			"the hypre package necessary for this preconditioner."));
#endif
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


  PreconditionLU::PreconditionLU ()
  {}

  
  PreconditionLU::PreconditionLU (const MatrixBase     &matrix,
				  const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);    
  }


  void
  PreconditionLU::initialize (const MatrixBase     &matrix_,
			      const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;
    
    create_pc();

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

    ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
