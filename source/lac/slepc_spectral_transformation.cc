// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/slepc_spectral_transformation.h>

#ifdef DEAL_II_WITH_SLEPC

#  include <deal.II/lac/slepc_solver.h>
#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/lac/petsc_vector.h>

#  include <cmath>
#  include <vector>

#  include <petscversion.h>

DEAL_II_NAMESPACE_OPEN

namespace SLEPcWrappers
{
  TransformationBase::TransformationData::~TransformationData ()
  {}

  TransformationBase::TransformationBase ()
  {}

  TransformationBase::~TransformationBase ()
  {}

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
#if DEAL_II_PETSC_VERSION_LT(3,1,0)
    ierr = STSetType (st, const_cast<char *>(STSINV));
#else
    ierr = STSetType (st, const_cast<char *>(STSINVERT));
#endif
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
#if DEAL_II_PETSC_VERSION_LT(3,5,0)
    int ierr;
    ierr = STSetType (st, const_cast<char *>(STFOLD));
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift (st, additional_data.shift_parameter);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));
#else
    // PETSc/SLEPc version must be < 3.5.0.
    Assert ((false),
            ExcMessage ("Folding transformation has been removed in SLEPc 3.5.0 and newer."
                        "You cannot use this transformation anymore."));
#endif
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

#endif // DEAL_II_WITH_SLEPC

