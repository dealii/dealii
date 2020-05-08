// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/slepc_spectral_transformation.h>

#ifdef DEAL_II_WITH_SLEPC

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/slepc_solver.h>

#  include <petscversion.h>

#  include <cmath>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace SLEPcWrappers
{
  TransformationBase::TransformationBase(const MPI_Comm &mpi_communicator)
  {
    const PetscErrorCode ierr = STCreate(mpi_communicator, &st);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  TransformationBase::~TransformationBase()
  {
    if (st != nullptr)
      {
        const PetscErrorCode ierr = STDestroy(&st);
        (void)ierr;
        AssertNothrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
      }
  }

  void
  TransformationBase::set_matrix_mode(const STMatMode mode)
  {
    const PetscErrorCode ierr = STSetMatMode(st, mode);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  void
  TransformationBase::set_solver(const PETScWrappers::SolverBase &solver)
  {
    PetscErrorCode ierr = STSetKSP(st, solver.solver_data->ksp);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* ------------------- TransformationShift --------------------- */

  TransformationShift::AdditionalData::AdditionalData(
    const double shift_parameter)
    : shift_parameter(shift_parameter)
  {}

  TransformationShift::TransformationShift(const MPI_Comm &mpi_communicator,
                                           const AdditionalData &data)
    : TransformationBase(mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = STSetType(st, const_cast<char *>(STSHIFT));
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift(st, additional_data.shift_parameter);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* ---------------- TransformationShiftInvert ------------------ */

  TransformationShiftInvert::AdditionalData::AdditionalData(
    const double shift_parameter)
    : shift_parameter(shift_parameter)
  {}

  TransformationShiftInvert::TransformationShiftInvert(
    const MPI_Comm &      mpi_communicator,
    const AdditionalData &data)
    : TransformationBase(mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = STSetType(st, const_cast<char *>(STSINVERT));
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift(st, additional_data.shift_parameter);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* --------------- TransformationSpectrumFolding ----------------- */

  TransformationSpectrumFolding::AdditionalData::AdditionalData(
    const double shift_parameter)
    : shift_parameter(shift_parameter)
  {}

  TransformationSpectrumFolding::TransformationSpectrumFolding(
    const MPI_Comm &      mpi_communicator,
    const AdditionalData &data)
    : TransformationBase(mpi_communicator)
    , additional_data(data)
  {
#  if DEAL_II_PETSC_VERSION_LT(3, 5, 0)
    PetscErrorCode ierr = STSetType(st, const_cast<char *>(STFOLD));
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift(st, additional_data.shift_parameter);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
#  else
    // PETSc/SLEPc version must be < 3.5.0.
    (void)st;
    Assert((false),
           ExcMessage(
             "Folding transformation has been removed in SLEPc 3.5.0 and newer."
             " You cannot use this transformation anymore."));
#  endif
  }

  /* ------------------- TransformationCayley --------------------- */

  TransformationCayley::AdditionalData::AdditionalData(
    const double shift_parameter,
    const double antishift_parameter)
    : shift_parameter(shift_parameter)
    , antishift_parameter(antishift_parameter)
  {}

  TransformationCayley::TransformationCayley(const MPI_Comm &mpi_communicator,
                                             const AdditionalData &data)
    : TransformationBase(mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = STSetType(st, const_cast<char *>(STCAYLEY));
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift(st, additional_data.shift_parameter);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STCayleySetAntishift(st, additional_data.antishift_parameter);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SLEPC
