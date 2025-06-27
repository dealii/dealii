// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
  TransformationBase::TransformationBase(const MPI_Comm mpi_communicator)
  {
    const PetscErrorCode ierr = STCreate(mpi_communicator, &st);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  TransformationBase::~TransformationBase()
  {
    if (st != nullptr)
      {
        const PetscErrorCode ierr = STDestroy(&st);
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
    PetscErrorCode ierr = STSetKSP(st, solver);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* ------------------- TransformationShift --------------------- */

  TransformationShift::AdditionalData::AdditionalData(
    const double shift_parameter)
    : shift_parameter(shift_parameter)
  {}

  TransformationShift::TransformationShift(const MPI_Comm mpi_communicator,
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
    const MPI_Comm        mpi_communicator,
    const AdditionalData &data)
    : TransformationBase(mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = STSetType(st, const_cast<char *>(STSINVERT));
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

    ierr = STSetShift(st, additional_data.shift_parameter);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  /* ------------------- TransformationCayley --------------------- */

  TransformationCayley::AdditionalData::AdditionalData(
    const double shift_parameter,
    const double antishift_parameter)
    : shift_parameter(shift_parameter)
    , antishift_parameter(antishift_parameter)
  {}

  TransformationCayley::TransformationCayley(const MPI_Comm mpi_communicator,
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
