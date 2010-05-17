//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//    Author: Toby D. Young, Polish Academy of Sciences, 2009
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/slepc_spectral_transformation.h>

#ifdef DEAL_II_USE_SLEPC

#  include <lac/slepc_solver.h>
#  include <lac/petsc_matrix_base.h>
#  include <lac/petsc_vector_base.h>
#  include <lac/petsc_vector.h>

#  include <cmath>
#  include <vector>

#  include <petscversion.h>

DEAL_II_NAMESPACE_OPEN

namespace SLEPcWrappers
{
  TransformationBase::TransformationData::~TransformationData ()
  {
  }

  TransformationBase::TransformationBase ()
  {
  }

  TransformationBase::~TransformationBase ()
  {
  }

  void TransformationBase::set_context (EPS &eps)
  {
    AssertThrow (transformation_data.get() == 0,
		 SolverBase::ExcSLEPcWrappersUsageError());
    transformation_data.reset (new TransformationData());

    int ierr = EPSGetST(eps, &transformation_data->st);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    set_transformation_type(transformation_data->st);
  }

  /* ------------------- TransformationShift --------------------- */
  TransformationShift::AdditionalData::
  AdditionalData (const double shift_parameter)
                  :
                  shift_parameter (shift_parameter)
  {}

  TransformationShift::TransformationShift (const AdditionalData &data)
    :
    additional_data (data)
  {}

  void
  TransformationShift::set_transformation_type (ST &st) const
  {
    int ierr;
    ierr = STSetType (st, const_cast<char *>(STSHIFT));
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift (st, additional_data.shift_parameter);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* ---------------- TransformationShiftInvert ------------------ */
  TransformationShiftInvert::AdditionalData::
  AdditionalData (const double shift_parameter)
                  :
                  shift_parameter (shift_parameter)
  {}

  TransformationShiftInvert::TransformationShiftInvert (const AdditionalData &data)
    :
    additional_data (data)
  {}

  void
  TransformationShiftInvert::set_transformation_type (ST &st) const
  {
    int ierr;
    ierr = STSetType (st, const_cast<char *>(STSINV));
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift (st, additional_data.shift_parameter);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* --------------- TransformationSpectrumFolding ----------------- */
  TransformationSpectrumFolding::AdditionalData::
  AdditionalData (const double shift_parameter)
                  :
                  shift_parameter (shift_parameter)
  {}

  TransformationSpectrumFolding::TransformationSpectrumFolding (const AdditionalData &data)
    :
    additional_data (data)
  {}


  void
  TransformationSpectrumFolding::set_transformation_type (ST &st) const
  {
    int ierr;
    ierr = STSetType (st, const_cast<char *>(STFOLD));
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift (st, additional_data.shift_parameter);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* ------------------- TransformationCayley --------------------- */
  TransformationCayley::AdditionalData::
  AdditionalData (const double shift_parameter,
		  const double antishift_parameter)
                  :
                  shift_parameter (shift_parameter),
		  antishift_parameter (antishift_parameter)
  {
  }

  TransformationCayley::TransformationCayley (const double shift,
					      const double antishift)
    :
    additional_data (shift, antishift)
  {}

  void
  TransformationCayley::set_transformation_type (ST &st) const
  {
    int ierr = STSetType (st, const_cast<char *>(STCAYLEY));
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift (st, additional_data.shift_parameter);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STCayleySetAntishift (st, additional_data.antishift_parameter);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

}

DEAL_II_NAMESPACE_CLOSE

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_SLEPC

