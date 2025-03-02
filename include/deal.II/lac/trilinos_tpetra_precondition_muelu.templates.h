// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_precondition_muelu_templates_h
#define dealii_trilinos_tpetra_precondition_muelu_templates_h


#include <deal.II/base/config.h>

#include "deal.II/base/exceptions.h"
#include <deal.II/base/index_set.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#include <MueLu_TpetraOperator_decl.hpp>
#include <Teuchos_ParameterList.hpp>

#include <string>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  include <deal.II/lac/trilinos_tpetra_precondition.h>
#  include <deal.II/lac/trilinos_tpetra_precondition.templates.h>

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU
#    include <MueLu_CreateTpetraPreconditioner.hpp>
#  endif

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU

    /* ---------------------- PreconditionAMGMueLu ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionAMGMueLu<Number, MemorySpace>::AdditionalData::AdditionalData(
      const bool         elliptic,
      const bool         symmetric,
      const bool         w_cycle,
      const double       aggregation_threshold,
      const int          smoother_sweeps,
      const int          smoother_overlap,
      const bool         output_details,
      const std::string &smoother_type,
      const std::string &coarse_type)
      : elliptic(elliptic)
      , symmetric(symmetric)
      , w_cycle(w_cycle)
      , aggregation_threshold(aggregation_threshold)
      , smoother_sweeps(smoother_sweeps)
      , smoother_overlap(smoother_overlap)
      , output_details(output_details)
      , smoother_type(smoother_type)
      , coarse_type(coarse_type)
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionAMGMueLu<Number, MemorySpace>::initialize(
      SparseMatrix<Number, MemorySpace> &A,
      Teuchos::ParameterList            &parameters)
    {
      this->parameter_list.setParameters(parameters);
      Teuchos::RCP<Tpetra::Operator<Number,
                                    TpetraTypes::LO,
                                    TpetraTypes::GO,
                                    TpetraTypes::NodeType<MemorySpace>>>
        op = A.trilinos_rcp();
      Teuchos::RCP<MueLu::TpetraOperator<Number,
                                         TpetraTypes::LO,
                                         TpetraTypes::GO,
                                         TpetraTypes::NodeType<MemorySpace>>>
        precon             = MueLu::CreateTpetraPreconditioner(op, parameters);
      this->preconditioner = precon;
    }



    template <typename Number, typename MemorySpace>
    void
    PreconditionAMGMueLu<Number, MemorySpace>::initialize(
      SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData              &ad)
    {
      Teuchos::ParameterList parameter_list;

      if (ad.elliptic)
        {
          parameter_list.set("multigrid algorithm", "sa");
        }
      else
        {
          parameter_list.set("multigrid algorithm", "unsmoothed");
        }
      parameter_list.set("problem: symmetric", ad.symmetric);

      Teuchos::ParameterList &smootherParams =
        parameter_list.sublist("smoother: params");

      if (ad.smoother_type == "Jacobi")
        {
          parameter_list.set("smoother: type", "RELAXATION");
          smootherParams.set("relaxation: sweeps", ad.smoother_sweeps);
          smootherParams.set("relaxation: type", "Jacobi");
        }
      else if (ad.smoother_type == "l1 Jacobi")
        {
          parameter_list.set("smoother: type", "RELAXATION");
          smootherParams.set("relaxation: sweeps", ad.smoother_sweeps);
          smootherParams.set("relaxation: type", "Jacobi");
          smootherParams.set("relaxation: use l1", true);
        }
      else if (ad.smoother_type == "Gauss Seidel")
        {
          parameter_list.set("smoother: type", "RELAXATION");
          smootherParams.set("relaxation: sweeps", ad.smoother_sweeps);
          smootherParams.set("relaxation: type", "Gauss-Seidel");
        }
      else if (ad.smoother_type == "l1 Gauss Seidel")
        {
          parameter_list.set("smoother: type", "RELAXATION");
          smootherParams.set("relaxation: sweeps", ad.smoother_sweeps);
          smootherParams.set("relaxation: type", "Gauss-Seidel");
          smootherParams.set("relaxation: use l1", true);
        }
      else if (ad.smoother_type == "Symmetric Gauss Seidel")
        {
          parameter_list.set("smoother: type", "RELAXATION");
          smootherParams.set("relaxation: sweeps", ad.smoother_sweeps);
          smootherParams.set("relaxation: type", "Symmetric Gauss-Seidel");
        }
      else if (ad.smoother_type == "Chebyshev")
        {
          parameter_list.set("smoother: type", "CHEBYSHEV");
          smootherParams.set("smoother: sweeps", ad.smoother_sweeps);
        }
      else if (ad.smoother_type == "ILU")
        {
          parameter_list.set("smoother: type", "RILUK");
        }
      else if (ad.smoother_type == "ILUT")
        {
          parameter_list.set("smoother: type", "ILUT");
        }
      else if (ad.smoother_type == "KLU2")
        {
          parameter_list.set("smoother: type", "RELAXATION");
        }
      else if (ad.smoother_type == "SuperLU")
        {
          parameter_list.set("smoother: type", "RELAXATION");
        }
      else if (ad.smoother_type == "SuperLU_dist")
        {
          parameter_list.set("smoother: type", "RELAXATION");
        }
      else
        {
          AssertThrow(false,
                      ExcTrilinosMueLuSmootherUnsupported(ad.smoother_type));
        }
      parameter_list.set("smoother: overlap", ad.smoother_overlap);


      Teuchos::ParameterList &coarseParams =
        parameter_list.sublist("coarse: params");
      if (ad.coarse_type == "Jacobi")
        {
          parameter_list.set("coarse: type", "RELAXATION");
          coarseParams.set("relaxation: type", "Jacobi");
        }
      else if (ad.coarse_type == "l1 Jacobi")
        {
          parameter_list.set("coarse: type", "RELAXATION");
          coarseParams.set("relaxation: type", "Jacobi");
          coarseParams.set("relaxation: use l1", true);
        }
      else if (ad.coarse_type == "Gauss Seidel")
        {
          parameter_list.set("coarse: type", "RELAXATION");
          coarseParams.set("relaxation: type", "Gauss-Seidel");
        }
      else if (ad.coarse_type == "l1 Gauss Seidel")
        {
          parameter_list.set("coarse: type", "RELAXATION");
          coarseParams.set("relaxation: type", "Gauss-Seidel");
          coarseParams.set("relaxation: use l1", true);
        }
      else if (ad.coarse_type == "Symmetric Gauss Seidel")
        {
          parameter_list.set("coarse: type", "RELAXATION");
          coarseParams.set("relaxation: type", "Symmetric Gauss-Seidel");
        }
      else if (ad.coarse_type == "Chebyshev")
        {
          parameter_list.set("coarse: type", "CHEBYSHEV");
        }
      else if (ad.coarse_type == "ILU")
        {
          parameter_list.set("coarse: type", "RILUK");
        }
      else if (ad.coarse_type == "ILUT")
        {
          parameter_list.set("coarse: type", "ILUT");
        }
      else if (ad.coarse_type == "KLU2")
        {
          parameter_list.set("coarse: type", "RELAXATION");
        }
      else if (ad.coarse_type == "SuperLU")
        {
          parameter_list.set("coarse: type", "RELAXATION");
        }
      else if (ad.coarse_type == "SuperLU_dist")
        {
          parameter_list.set("coarse: type", "RELAXATION");
        }
      else
        {
          AssertThrow(false,
                      ExcTrilinosMueLuCoarseSolverUnsupported(ad.coarse_type));
        }

      if (ad.w_cycle)
        parameter_list.set("cycle type", "W");
      else
        parameter_list.set("cycle type", "V");

      parameter_list.set("aggregation: drop tol", ad.aggregation_threshold);

      parameter_list.set("coarse: max size", 2000);

      if (ad.output_details == true)
        parameter_list.set("verbosity", "high");
      else
        parameter_list.set("verbosity", "none");

      initialize(A, parameter_list);
    }

#    endif // DEAL_II_TRILINOS_WITH_TPETRA_MUELU
#  endif   // DEAL_II_TRILINOS_WITH_IFPACK2
  }        // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif
