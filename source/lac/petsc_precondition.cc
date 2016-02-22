// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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

#include <deal.II/lac/petsc_precondition.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/utilities.h>
#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/lac/petsc_solver.h>
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
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
        int ierr = PCDestroy(pc);
#else
        int ierr = PCDestroy(&pc);
#endif
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
    // this ugly cast is necessary because the
    // type Mat and PETScObject are
    // unrelated.
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCCreate(comm, &pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(3, 5, 0)
    ierr = PCSetOperators(pc , matrix, matrix, SAME_PRECONDITIONER);
#else
    ierr = PCSetOperators(pc , matrix, matrix);
#endif
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
  PreconditionJacobi::PreconditionJacobi (const MPI_Comm comm,
                                          const AdditionalData &additional_data_)
  {
    additional_data = additional_data_;

    int ierr = PCCreate(comm, &pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    initialize();
  }


  PreconditionJacobi::PreconditionJacobi ()
  {}


  PreconditionJacobi::PreconditionJacobi (const MatrixBase     &matrix,
                                          const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);
  }

  void
  PreconditionJacobi::initialize()
  {
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCJACOBI));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  void
  PreconditionJacobi::initialize (const MatrixBase     &matrix_,
                                  const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;

    create_pc();
    initialize();

    int ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


  /* ----------------- PreconditionBlockJacobi -------------------- */
  PreconditionBlockJacobi::PreconditionBlockJacobi (const MPI_Comm comm,
                                                    const AdditionalData &additional_data_)
  {
    additional_data = additional_data_;

    int ierr = PCCreate(comm, &pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    initialize();
  }


  PreconditionBlockJacobi::PreconditionBlockJacobi ()
  {}


  PreconditionBlockJacobi::
  PreconditionBlockJacobi (const MatrixBase     &matrix,
                           const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);
  }

  void
  PreconditionBlockJacobi::initialize()
  {
    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCBJACOBI));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


  void
  PreconditionBlockJacobi::initialize (const MatrixBase     &matrix_,
                                       const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;

    create_pc();
    initialize();

    int ierr = PCSetUp (pc);
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
    PCFactorSetLevels (pc, additional_data.levels);
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
    PCFactorSetLevels (pc, additional_data.levels);
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

  PreconditionBoomerAMG::PreconditionBoomerAMG (const MPI_Comm comm,
                                                const AdditionalData &additional_data_)
  {
    additional_data = additional_data_;

    int ierr = PCCreate(comm, &pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#ifdef PETSC_HAVE_HYPRE
    initialize();
#else // PETSC_HAVE_HYPRE
    (void)pc;
    Assert (false,
            ExcMessage ("Your PETSc installation does not include a copy of "
                        "the hypre package necessary for this preconditioner."));
#endif
  }


  PreconditionBoomerAMG::PreconditionBoomerAMG (const MatrixBase     &matrix,
                                                const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);
  }

  void
  PreconditionBoomerAMG::initialize ()
  {
#ifndef PETSC_USE_COMPLEX
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
#else
    Assert(false, ExcMessage("Complex-valued PETSc does not support BoomerAMG preconditioner."));
#endif
  }

  void
  PreconditionBoomerAMG::initialize (const MatrixBase     &matrix_,
                                     const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;

#ifdef PETSC_HAVE_HYPRE
    create_pc();
    initialize ();

    int ierr = PCSetUp (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#else // PETSC_HAVE_HYPRE
    (void)pc;
    Assert (false,
            ExcMessage ("Your PETSc installation does not include a copy of "
                        "the hypre package necessary for this preconditioner."));
#endif
  }


  /* ----------------- PreconditionParaSails -------------------- */

  PreconditionParaSails::AdditionalData::
  AdditionalData(const unsigned int symmetric,
                 const unsigned int n_levels,
                 const double threshold,
                 const double filter,
                 const bool output_details)
    :
    symmetric(symmetric),
    n_levels(n_levels),
    threshold(threshold),
    filter(filter),
    output_details(output_details)
  {}


  PreconditionParaSails::PreconditionParaSails ()
  {}


  PreconditionParaSails::PreconditionParaSails (const MatrixBase     &matrix,
                                                const AdditionalData &additional_data)
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionParaSails::initialize (const MatrixBase     &matrix_,
                                     const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;

#ifdef PETSC_HAVE_HYPRE
    create_pc();

    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCHYPRE));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCHYPRESetType(pc, "parasails");
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    if (additional_data.output_details)
      PetscOptionsSetValue("-pc_hypre_parasails_logging","1");

    Assert ((additional_data.symmetric == 0 ||
             additional_data.symmetric == 1 ||
             additional_data.symmetric == 2),
            ExcMessage("ParaSails parameter symmetric can only be equal to 0, 1, 2!"));

    std::stringstream ssStream;

    switch (additional_data.symmetric)
      {
      case 0:
      {
        ssStream << "nonsymmetric";
        break;
      }

      case 1:
      {
        ssStream << "SPD";
        break;
      }

      case 2:
      {
        ssStream << "nonsymmetric,SPD";
        break;
      }

      default:
        Assert (false,
                ExcMessage("ParaSails parameter symmetric can only be equal to 0, 1, 2!"));
      };

    PetscOptionsSetValue("-pc_hypre_parasails_sym",ssStream.str().c_str());

    PetscOptionsSetValue("-pc_hypre_parasails_nlevels",
                         Utilities::int_to_string(
                           additional_data.n_levels
                         ).c_str());

    ssStream.str(""); // empty the stringstream
    ssStream << additional_data.threshold;
    PetscOptionsSetValue("-pc_hypre_parasails_thresh", ssStream.str().c_str());

    ssStream.str(""); // empty the stringstream
    ssStream << additional_data.filter;
    PetscOptionsSetValue("-pc_hypre_parasails_filter", ssStream.str().c_str());

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


  /* ----------------- PreconditionNone ------------------------- */

  PreconditionNone::PreconditionNone ()
  {}


  PreconditionNone::PreconditionNone (const MatrixBase     &matrix,
                                      const AdditionalData &additional_data)
  {
    initialize (matrix, additional_data);
  }


  void
  PreconditionNone::initialize (const MatrixBase     &matrix_,
                                const AdditionalData &additional_data_)
  {
    matrix = static_cast<Mat>(matrix_);
    additional_data = additional_data_;

    create_pc();

    int ierr;
    ierr = PCSetType (pc, const_cast<char *>(PCNONE));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions (pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCSetUp (pc);
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
#if DEAL_II_PETSC_VERSION_LT(3,0,1)
    ierr = PCFactorSetPivoting (pc, additional_data.pivoting);
#else
    ierr = PCFactorSetColumnPivot (pc, additional_data.pivoting);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = PCFactorSetZeroPivot (pc, additional_data.zero_pivot);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(3,0,1)
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

#endif // DEAL_II_WITH_PETSC
