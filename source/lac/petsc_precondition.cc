// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/petsc_precondition.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_solver.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscconf.h>

#  include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  PreconditionBase::PreconditionBase(const MPI_Comm comm)
    : pc(nullptr)
  {
    create_pc_with_comm(comm);
  }

  PreconditionBase::PreconditionBase()
    : pc(nullptr)
  {}

  PreconditionBase::~PreconditionBase()
  {
    try
      {
        clear();
      }
    catch (...)
      {}
  }

  void
  PreconditionBase::clear()
  {
    if (pc)
      {
        PetscErrorCode ierr = PCDestroy(&pc);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
  }

  void
  PreconditionBase::vmult(VectorBase &dst, const VectorBase &src) const
  {
    AssertThrow(pc != nullptr, StandardExceptions::ExcInvalidState());

    PetscErrorCode ierr = PCApply(pc, src, dst);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  PreconditionBase::Tvmult(VectorBase &dst, const VectorBase &src) const
  {
    AssertThrow(pc != nullptr, StandardExceptions::ExcInvalidState());

    PetscErrorCode ierr = PCApplyTranspose(pc, src, dst);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  PreconditionBase::setup()
  {
    AssertThrow(pc != nullptr, StandardExceptions::ExcInvalidState());

    PetscErrorCode ierr = PCSetUp(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  MPI_Comm
  PreconditionBase::get_mpi_communicator() const
  {
    return PetscObjectComm(reinterpret_cast<PetscObject>(pc));
  }

  void
  PreconditionBase::create_pc_with_mat(const MatrixBase &matrix)
  {
    // only allow the creation of the
    // preconditioner once
    AssertThrow(pc == nullptr, StandardExceptions::ExcInvalidState());

    MPI_Comm       comm;
    PetscErrorCode ierr = PetscObjectGetComm(
      reinterpret_cast<PetscObject>(static_cast<const Mat &>(matrix)), &comm);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    create_pc_with_comm(comm);

    ierr = PCSetOperators(pc, matrix, matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  PreconditionBase::create_pc_with_comm(const MPI_Comm comm)
  {
    clear();
    PetscErrorCode ierr = PCCreate(comm, &pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  const PC &
  PreconditionBase::get_pc() const
  {
    return pc;
  }


  /* ----------------- PreconditionJacobi -------------------- */

  PreconditionJacobi::PreconditionJacobi()
    : PreconditionBase()
  {}



  PreconditionJacobi::PreconditionJacobi(const MPI_Comm        comm,
                                         const AdditionalData &additional_data_)
    : PreconditionBase(comm)
  {
    additional_data = additional_data_;

    PetscErrorCode ierr = PCCreate(comm, &pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    initialize();
  }



  PreconditionJacobi::PreconditionJacobi(const MatrixBase     &matrix,
                                         const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }



  void
  PreconditionJacobi::initialize()
  {
    AssertThrow(pc != nullptr, StandardExceptions::ExcInvalidState());

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCJACOBI));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  PreconditionJacobi::initialize(const MatrixBase     &matrix_,
                                 const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);
    initialize();
  }


  /* ----------------- PreconditionBlockJacobi -------------------- */

  PreconditionBlockJacobi::PreconditionBlockJacobi()
    : PreconditionBase()
  {}

  PreconditionBlockJacobi::PreconditionBlockJacobi(
    const MPI_Comm        comm,
    const AdditionalData &additional_data_)
    : PreconditionBase(comm)
  {
    additional_data = additional_data_;

    PetscErrorCode ierr = PCCreate(comm, &pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    initialize();
  }



  PreconditionBlockJacobi::PreconditionBlockJacobi(
    const MatrixBase     &matrix,
    const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }



  void
  PreconditionBlockJacobi::initialize()
  {
    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCBJACOBI));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  PreconditionBlockJacobi::initialize(const MatrixBase     &matrix_,
                                      const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);
    initialize();
  }


  /* ----------------- PreconditionSOR -------------------- */

  PreconditionSOR::PreconditionSOR()
    : PreconditionBase()
  {}



  PreconditionSOR::AdditionalData::AdditionalData(const double omega)
    : omega(omega)
  {}



  PreconditionSOR::PreconditionSOR(const MatrixBase     &matrix,
                                   const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionSOR::initialize(const MatrixBase     &matrix_,
                              const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCSOR));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then set flags as given
    ierr = PCSORSetOmega(pc, additional_data.omega);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ----------------- PreconditionSSOR -------------------- */

  PreconditionSSOR::PreconditionSSOR()
    : PreconditionBase()
  {}



  PreconditionSSOR::AdditionalData::AdditionalData(const double omega)
    : omega(omega)
  {}



  PreconditionSSOR::PreconditionSSOR(const MatrixBase     &matrix,
                                     const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionSSOR::initialize(const MatrixBase     &matrix_,
                               const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCSOR));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then set flags as given
    ierr = PCSORSetOmega(pc, additional_data.omega);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // convert SOR to SSOR
    ierr = PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ----------------- PreconditionICC -------------------- */

  PreconditionICC::PreconditionICC()
    : PreconditionBase()
  {}



  PreconditionICC::AdditionalData::AdditionalData(const unsigned int levels)
    : levels(levels)
  {}



  PreconditionICC::PreconditionICC(const MatrixBase     &matrix,
                                   const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionICC::initialize(const MatrixBase     &matrix_,
                              const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCICC));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then set flags
    ierr = PCFactorSetLevels(pc, additional_data.levels);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ----------------- PreconditionILU -------------------- */

  PreconditionILU::PreconditionILU()
    : PreconditionBase()
  {}



  PreconditionILU::AdditionalData::AdditionalData(const unsigned int levels)
    : levels(levels)
  {}



  PreconditionILU::PreconditionILU(const MatrixBase     &matrix,
                                   const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionILU::initialize(const MatrixBase     &matrix_,
                              const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCILU));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then set flags
    ierr = PCFactorSetLevels(pc, additional_data.levels);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ----------------- PreconditionBoomerAMG -------------------- */

  PreconditionBoomerAMG::AdditionalData::AdditionalData(
    const bool           symmetric_operator,
    const double         strong_threshold,
    const double         max_row_sum,
    const unsigned int   aggressive_coarsening_num_levels,
    const bool           output_details,
    const RelaxationType relaxation_type_up,
    const RelaxationType relaxation_type_down,
    const RelaxationType relaxation_type_coarse,
    const unsigned int   n_sweeps_coarse,
    const double         tol,
    const unsigned int   max_iter,
    const bool           w_cycle)
    : symmetric_operator(symmetric_operator)
    , strong_threshold(strong_threshold)
    , max_row_sum(max_row_sum)
    , aggressive_coarsening_num_levels(aggressive_coarsening_num_levels)
    , output_details(output_details)
    , relaxation_type_up(relaxation_type_up)
    , relaxation_type_down(relaxation_type_down)
    , relaxation_type_coarse(relaxation_type_coarse)
    , n_sweeps_coarse(n_sweeps_coarse)
    , tol(tol)
    , max_iter(max_iter)
    , w_cycle(w_cycle)
  {}



#  ifdef DEAL_II_PETSC_WITH_HYPRE
  namespace
  {
    /**
     * Converts the enums for the different relaxation types to the respective
     * strings for PETSc.
     */
    std::string
    to_string(
      PreconditionBoomerAMG::AdditionalData::RelaxationType relaxation_type)
    {
      std::string string_type;

      switch (relaxation_type)
        {
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::Jacobi:
            string_type = "Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            sequentialGaussSeidel:
            string_type = "sequential-Gauss-Seidel";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            seqboundaryGaussSeidel:
            string_type = "seqboundary-Gauss-Seidel";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::SORJacobi:
            string_type = "SOR/Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            backwardSORJacobi:
            string_type = "backward-SOR/Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            symmetricSORJacobi:
            string_type = "symmetric-SOR/Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            l1scaledSORJacobi:
            string_type = " l1scaled-SOR/Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            GaussianElimination:
            string_type = "Gaussian-elimination";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            l1GaussSeidel:
            string_type = "l1-Gauss-Seidel";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            backwardl1GaussSeidel:
            string_type = "backward-l1-Gauss-Seidel";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::CG:
            string_type = "CG";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::Chebyshev:
            string_type = "Chebyshev";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::FCFJacobi:
            string_type = "FCF-Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::
            l1scaledJacobi:
            string_type = "l1scaled-Jacobi";
            break;
          case PreconditionBoomerAMG::AdditionalData::RelaxationType::None:
            string_type = "None";
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }
      return string_type;
    }
  } // namespace
#  endif



  PreconditionBoomerAMG::PreconditionBoomerAMG()
    : PreconditionBase()
  {}



  PreconditionBoomerAMG::PreconditionBoomerAMG(
    const MPI_Comm        comm,
    const AdditionalData &additional_data_)
    : PreconditionBase(comm)
  {
    additional_data = additional_data_;

    PetscErrorCode ierr = PCCreate(comm, &pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

#  ifdef DEAL_II_PETSC_WITH_HYPRE
    initialize();
#  else // DEAL_II_PETSC_WITH_HYPRE
    Assert(false,
           ExcMessage("Your PETSc installation does not include a copy of "
                      "the hypre package necessary for this preconditioner."));
#  endif
  }



  PreconditionBoomerAMG::PreconditionBoomerAMG(
    const MatrixBase     &matrix,
    const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }



  void
  PreconditionBoomerAMG::initialize()
  {
#  ifdef DEAL_II_PETSC_WITH_HYPRE
    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCHYPRE));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCHYPRESetType(pc, "boomeramg");
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    if (additional_data.output_details)
      {
        set_option_value("-pc_hypre_boomeramg_print_statistics", "1");
      }

    set_option_value("-pc_hypre_boomeramg_agg_nl",
                     std::to_string(
                       additional_data.aggressive_coarsening_num_levels));

    set_option_value("-pc_hypre_boomeramg_max_row_sum",
                     std::to_string(additional_data.max_row_sum));

    set_option_value("-pc_hypre_boomeramg_strong_threshold",
                     std::to_string(additional_data.strong_threshold));

    // change to symmetric SOR/Jacobi when using a symmetric operator for
    // backward compatibility
    if (additional_data.relaxation_type_up ==
          AdditionalData::RelaxationType::SORJacobi &&
        additional_data.symmetric_operator)
      {
        additional_data.relaxation_type_up =
          AdditionalData::RelaxationType::symmetricSORJacobi;
      }

    if (additional_data.relaxation_type_down ==
          AdditionalData::RelaxationType::SORJacobi &&
        additional_data.symmetric_operator)
      {
        additional_data.relaxation_type_down =
          AdditionalData::RelaxationType::symmetricSORJacobi;
      }

    if (additional_data.relaxation_type_coarse ==
          AdditionalData::RelaxationType::SORJacobi &&
        additional_data.symmetric_operator)
      {
        additional_data.relaxation_type_coarse =
          AdditionalData::RelaxationType::symmetricSORJacobi;
      }

    auto relaxation_type_is_symmetric =
      [](AdditionalData::RelaxationType relaxation_type) {
        return relaxation_type == AdditionalData::RelaxationType::Jacobi ||
               relaxation_type ==
                 AdditionalData::RelaxationType::symmetricSORJacobi ||
               relaxation_type ==
                 AdditionalData::RelaxationType::GaussianElimination ||
               relaxation_type == AdditionalData::RelaxationType::None ||
               relaxation_type ==
                 AdditionalData::RelaxationType::l1scaledJacobi ||
               relaxation_type == AdditionalData::RelaxationType::CG ||
               relaxation_type == AdditionalData::RelaxationType::Chebyshev;
      };

    if (additional_data.symmetric_operator &&
        !relaxation_type_is_symmetric(additional_data.relaxation_type_up))
      Assert(false,
             ExcMessage("Use a symmetric smoother for relaxation_type_up"));

    if (additional_data.symmetric_operator &&
        !relaxation_type_is_symmetric(additional_data.relaxation_type_down))
      Assert(false,
             ExcMessage("Use a symmetric smoother for relaxation_type_down"));

    if (additional_data.symmetric_operator &&
        !relaxation_type_is_symmetric(additional_data.relaxation_type_coarse))
      Assert(false,
             ExcMessage("Use a symmetric smoother for relaxation_type_coarse"));

    set_option_value("-pc_hypre_boomeramg_relax_type_up",
                     to_string(additional_data.relaxation_type_up));
    set_option_value("-pc_hypre_boomeramg_relax_type_down",
                     to_string(additional_data.relaxation_type_down));
    set_option_value("-pc_hypre_boomeramg_relax_type_coarse",
                     to_string(additional_data.relaxation_type_coarse));
    set_option_value("-pc_hypre_boomeramg_grid_sweeps_coarse",
                     std::to_string(additional_data.n_sweeps_coarse));

    set_option_value("-pc_hypre_boomeramg_tol",
                     std::to_string(additional_data.tol));
    set_option_value("-pc_hypre_boomeramg_max_iter",
                     std::to_string(additional_data.max_iter));

    if (additional_data.w_cycle)
      {
        set_option_value("-pc_hypre_boomeramg_cycle_type", "W");
      }

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
#  else
    Assert(false,
           ExcMessage("Your PETSc installation does not include a copy of "
                      "the hypre package necessary for this preconditioner."));
#  endif
  }



  void
  PreconditionBoomerAMG::initialize(const MatrixBase     &matrix_,
                                    const AdditionalData &additional_data_)
  {
#  ifdef DEAL_II_PETSC_WITH_HYPRE
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);
    initialize();

#  else // DEAL_II_PETSC_WITH_HYPRE
    (void)matrix_;
    (void)additional_data_;
    Assert(false,
           ExcMessage("Your PETSc installation does not include a copy of "
                      "the hypre package necessary for this preconditioner."));
#  endif
  }


  /* ----------------- PreconditionParaSails -------------------- */

  PreconditionParaSails::AdditionalData::AdditionalData(
    const unsigned int symmetric,
    const unsigned int n_levels,
    const double       threshold,
    const double       filter,
    const bool         output_details)
    : symmetric(symmetric)
    , n_levels(n_levels)
    , threshold(threshold)
    , filter(filter)
    , output_details(output_details)
  {}



  PreconditionParaSails::PreconditionParaSails()
    : PreconditionBase()
  {}



  PreconditionParaSails::PreconditionParaSails(
    const MatrixBase     &matrix,
    const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionParaSails::initialize(const MatrixBase     &matrix_,
                                    const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

#  ifdef DEAL_II_PETSC_WITH_HYPRE
    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCHYPRE));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCHYPRESetType(pc, "parasails");
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    if (additional_data.output_details)
      {
        set_option_value("-pc_hypre_parasails_logging", "1");
      }

    Assert((additional_data.symmetric == 0 || additional_data.symmetric == 1 ||
            additional_data.symmetric == 2),
           ExcMessage(
             "ParaSails parameter symmetric can only be equal to 0, 1, 2!"));

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
          Assert(
            false,
            ExcMessage(
              "ParaSails parameter symmetric can only be equal to 0, 1, 2!"));
      }

    set_option_value("-pc_hypre_parasails_sym", ssStream.str());

    set_option_value("-pc_hypre_parasails_nlevels",
                     std::to_string(additional_data.n_levels));

    ssStream.str(""); // empty the stringstream
    ssStream << additional_data.threshold;
    set_option_value("-pc_hypre_parasails_thresh", ssStream.str());

    ssStream.str(""); // empty the stringstream
    ssStream << additional_data.filter;
    set_option_value("-pc_hypre_parasails_filter", ssStream.str());

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

#  else // DEAL_II_PETSC_WITH_HYPRE
    (void)matrix_;
    Assert(false,
           ExcMessage("Your PETSc installation does not include a copy of "
                      "the hypre package necessary for this preconditioner."));
#  endif
  }


  /* ----------------- PreconditionNone ------------------------- */

  PreconditionNone::PreconditionNone()
    : PreconditionBase()
  {}



  PreconditionNone::PreconditionNone(const MatrixBase     &matrix,
                                     const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionNone::initialize(const MatrixBase     &matrix_,
                               const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCNONE));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ----------------- PreconditionLU -------------------- */

  PreconditionLU::AdditionalData::AdditionalData(const double pivoting,
                                                 const double zero_pivot,
                                                 const double damping)
    : pivoting(pivoting)
    , zero_pivot(zero_pivot)
    , damping(damping)
  {}



  PreconditionLU::PreconditionLU()
    : PreconditionBase()
  {}



  PreconditionLU::PreconditionLU(const MatrixBase     &matrix,
                                 const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }


  void
  PreconditionLU::initialize(const MatrixBase     &matrix_,
                             const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);

    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCLU));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set flags as given
    ierr = PCFactorSetColumnPivot(pc, additional_data.pivoting);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCFactorSetZeroPivot(pc, additional_data.zero_pivot);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCFactorSetShiftAmount(pc, additional_data.damping);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  /* ----------------- PreconditionBDDC -------------------- */

  template <int dim>
  PreconditionBDDC<dim>::AdditionalData::AdditionalData(
    const bool                    use_vertices,
    const bool                    use_edges,
    const bool                    use_faces,
    const bool                    symmetric,
    const std::vector<Point<dim>> coords)
    : use_vertices(use_vertices)
    , use_edges(use_edges)
    , use_faces(use_faces)
    , symmetric(symmetric)
    , coords(coords)
  {}



  template <int dim>
  PreconditionBDDC<dim>::PreconditionBDDC()
    : PreconditionBase()
  {}



  template <int dim>
  PreconditionBDDC<dim>::PreconditionBDDC(
    const MPI_Comm        comm,
    const AdditionalData &additional_data_)
    : PreconditionBase(comm)
  {
    additional_data = additional_data_;

    PetscErrorCode ierr = PCCreate(comm, &pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    initialize();
  }



  template <int dim>
  PreconditionBDDC<dim>::PreconditionBDDC(const MatrixBase     &matrix,
                                          const AdditionalData &additional_data)
    : PreconditionBase(matrix.get_mpi_communicator())
  {
    initialize(matrix, additional_data);
  }



  template <int dim>
  void
  PreconditionBDDC<dim>::initialize()
  {
#  if DEAL_II_PETSC_VERSION_GTE(3, 10, 0)
    PetscErrorCode ierr = PCSetType(pc, const_cast<char *>(PCBDDC));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // The matrix must be of IS type. We check for this to avoid the PETSc error
    // in order to suggest the correct matrix reinit method.
    {
      MatType   current_type;
      Mat       A, P;
      PetscBool flg;

      ierr = PCGetOperators(pc, &A, &P);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      ierr = PCGetUseAmat(pc, &flg);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      ierr = MatGetType(flg ? A : P, &current_type);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      AssertThrow(
        strcmp(current_type, MATIS) == 0,
        ExcMessage(
          "Matrix must be of IS type. For this, the variant of reinit that includes the active dofs must be used."));
    }


    std::stringstream ssStream;

    if (additional_data.use_vertices)
      set_option_value("-pc_bddc_use_vertices", "true");
    else
      set_option_value("-pc_bddc_use_vertices", "false");
    if (additional_data.use_edges)
      set_option_value("-pc_bddc_use_edges", "true");
    else
      set_option_value("-pc_bddc_use_edges", "false");
    if (additional_data.use_faces)
      set_option_value("-pc_bddc_use_faces", "true");
    else
      set_option_value("-pc_bddc_use_faces", "false");
    if (additional_data.symmetric)
      set_option_value("-pc_bddc_symmetric", "true");
    else
      set_option_value("-pc_bddc_symmetric", "false");
    if (additional_data.coords.size() > 0)
      {
        set_option_value("-pc_bddc_corner_selection", "true");
        // Convert coords vector to PETSc data array
        std::vector<PetscReal> coords_petsc(additional_data.coords.size() *
                                            dim);
        for (unsigned int i = 0, j = 0; i < additional_data.coords.size(); ++i)
          {
            for (j = 0; j < dim; ++j)
              coords_petsc[dim * i + j] = additional_data.coords[i][j];
          }

        ierr = PCSetCoordinates(pc,
                                dim,
                                additional_data.coords.size(),
                                coords_petsc.data());
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    else
      {
        set_option_value("-pc_bddc_corner_selection", "false");
        ierr = PCSetCoordinates(pc, 0, 0, nullptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }


    ierr = PCSetFromOptions(pc);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
#  else
    AssertThrow(
      false, ExcMessage("BDDC preconditioner requires PETSc 3.10.0 or newer"));
#  endif
  }



  template <int dim>
  void
  PreconditionBDDC<dim>::initialize(const MatrixBase     &matrix_,
                                    const AdditionalData &additional_data_)
  {
    clear();

    additional_data = additional_data_;

    create_pc_with_mat(matrix_);
    initialize();
  }

  /* ----------------- PreconditionShell -------------------- */

  PreconditionShell::PreconditionShell(const MatrixBase &matrix)
  {
    initialize(matrix);
  }

  PreconditionShell::PreconditionShell(const MPI_Comm comm)
  {
    initialize(comm);
  }

  void
  PreconditionShell::initialize(const MPI_Comm comm)
  {
    PetscErrorCode ierr;
    if (pc)
      {
        ierr = PCDestroy(&pc);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    create_pc_with_comm(comm);

    ierr = PCSetType(pc, PCSHELL);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = PCShellSetContext(pc, static_cast<void *>(this));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = PCShellSetSetUp(pc, PreconditionShell::pcsetup);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = PCShellSetApply(pc, PreconditionShell::pcapply);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = PCShellSetApplyTranspose(pc, PreconditionShell::pcapply_transpose);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = PCShellSetName(pc, "deal.II user solve");
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  PreconditionShell::initialize(const MatrixBase &matrix)
  {
    initialize(matrix.get_mpi_communicator());
    PetscErrorCode ierr;
    ierr = PCSetOperators(pc, matrix, matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

#  ifndef PetscCall
#    define PetscCall(code)             \
      do                                \
        {                               \
          PetscErrorCode ierr = (code); \
          CHKERRQ(ierr);                \
        }                               \
      while (false)
#  endif

  PetscErrorCode
  PreconditionShell::pcsetup(PC ppc)
  {
    PetscFunctionBeginUser;
    // Failed reason is not reset uniformly within the
    // interface code of PCSetUp in PETSc.
    // We handle it here.
    PetscCall(pc_set_failed_reason(ppc, PC_NOERROR));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscErrorCode
  PreconditionShell::pcapply(PC ppc, Vec x, Vec y)
  {
    void *ctx;

    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(ppc, &ctx));

    auto *user = static_cast<PreconditionShell *>(ctx);
    if (!user->vmult)
      SETERRQ(
        PetscObjectComm((PetscObject)ppc),
        PETSC_ERR_LIB,
        "Failure in dealii::PETScWrappers::PreconditionShell::pcapply. Missing std::function vmult");

    VectorBase src(x);
    VectorBase dst(y);
    const int  lineno = __LINE__;
    try
      {
        user->vmult(dst, src);
      }
    catch (const RecoverableUserCallbackError &)
      {
        PetscCall(pc_set_failed_reason(ppc, PC_SUBPC_ERROR));
      }
    catch (...)
      {
        return PetscError(
          PetscObjectComm((PetscObject)ppc),
          lineno + 3,
          "vmult",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in pcapply from dealii::PETScWrappers::NonlinearSolver");
      }
    petsc_increment_state_counter(y);
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscErrorCode
  PreconditionShell::pcapply_transpose(PC ppc, Vec x, Vec y)
  {
    void *ctx;

    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(ppc, &ctx));

    auto *user = static_cast<PreconditionShell *>(ctx);
    if (!user->vmultT)
      SETERRQ(
        PetscObjectComm((PetscObject)ppc),
        PETSC_ERR_LIB,
        "Failure in dealii::PETScWrappers::PreconditionShell::pcapply_transpose. Missing std::function vmultT");

    VectorBase src(x);
    VectorBase dst(y);
    const int  lineno = __LINE__;
    try
      {
        user->vmultT(dst, src);
      }
    catch (const RecoverableUserCallbackError &)
      {
        PetscCall(pc_set_failed_reason(ppc, PC_SUBPC_ERROR));
      }
    catch (...)
      {
        return PetscError(
          PetscObjectComm((PetscObject)ppc),
          lineno + 3,
          "vmultT",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in pcapply_transpose from dealii::PETScWrappers::NonlinearSolver");
      }
    petsc_increment_state_counter(y);
    PetscFunctionReturn(PETSC_SUCCESS);
  }


} // namespace PETScWrappers

template class PETScWrappers::PreconditionBDDC<2>;
template class PETScWrappers::PreconditionBDDC<3>;

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
