// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_precondition_muelu_templates_h
#define dealii_trilinos_tpetra_precondition_muelu_templates_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/index_set.h>

#  include <deal.II/lac/trilinos_tpetra_precondition.h>
#  include <deal.II/lac/trilinos_tpetra_precondition.templates.h>
#  include <deal.II/lac/trilinos_tpetra_types.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <MueLu_CreateTpetraPreconditioner.hpp>
#  include <MueLu_TpetraOperator_decl.hpp>
#  include <Teuchos_ParameterList.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <string>

#endif // DEAL_II_TRILINOS_WITH_TPETRA_MUELU

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2

    /* ---------------------- PreconditionAMGMueLu ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionAMGMueLu<Number, MemorySpace>::AdditionalData::AdditionalData(
      const bool                            elliptic,
      const bool                            symmetric,
      const bool                            w_cycle,
      const double                          aggregation_threshold,
      const std::vector<std::vector<bool>> &constant_modes,
      const int                             smoother_sweeps,
      const int                             smoother_overlap,
      const bool                            output_details,
      const std::string                    &smoother_type,
      const std::string                    &coarse_type)
      : elliptic(elliptic)
      , symmetric(symmetric)
      , w_cycle(w_cycle)
      , aggregation_threshold(aggregation_threshold)
      , constant_modes(constant_modes)
      , smoother_sweeps(smoother_sweeps)
      , smoother_overlap(smoother_overlap)
      , output_details(output_details)
      , smoother_type(smoother_type)
      , coarse_type(coarse_type)
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionAMGMueLu<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const Teuchos::ParameterList            &parameters)
    {
      this->parameter_list = parameters;
      Teuchos::RCP<const Tpetra::Operator<Number,
                                          TpetraTypes::LO,
                                          TpetraTypes::GO,
                                          TpetraTypes::NodeType<MemorySpace>>>

        op = A.trilinos_rcp();

      // FIXME We shouldn't need to use const_cast here but up to at least
      // Trilinos 17.0.0 MueLu::CreateTpetraPreconditioner only works with
      // non-const arguments
      this->preconditioner = MueLu::CreateTpetraPreconditioner(
        Teuchos::rcp_const_cast<
          Tpetra::Operator<Number,
                           TpetraTypes::LO,
                           TpetraTypes::GO,
                           TpetraTypes::NodeType<MemorySpace>>>(op),
        this->parameter_list);
    }



    template <typename Number, typename MemorySpace>
    void
    PreconditionAMGMueLu<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
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
          parameter_list.set("coarse: type", "KLU2");
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

      using size_type = typename SparseMatrix<Number, MemorySpace>::size_type;
      const auto domain_map               = A.trilinos_rcp()->getDomainMap();
      const auto constant_modes_dimension = ad.constant_modes.size();

      // Create the nullspace vector
      Teuchos::RCP<TpetraTypes::MultiVectorType<Number, MemorySpace>>
        distributed_constant_modes = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MultiVectorType<Number, MemorySpace>>(
          domain_map,
          constant_modes_dimension > 0 ? constant_modes_dimension : 1);

      if (constant_modes_dimension > 0)
        {
          // MueLu benefits from knowing how many solution components
          // a vector-valued problem has. We do not have this information
          // readily available, except if the constant modes are given by
          // the user. In that case one constant mode corresponds to each
          // solution component. Provide a good guess if possible, otherwise use
          // the default of 1.
          parameter_list.set("number of equations",
                             static_cast<int>(constant_modes_dimension));

          const size_type n_rows = A.trilinos_rcp()->getGlobalNumRows();
          const bool      constant_modes_are_global =
            ad.constant_modes[0].size() == n_rows;
          const size_type n_relevant_rows =
            constant_modes_are_global ? n_rows : ad.constant_modes[0].size();
          const size_type my_size = domain_map->getLocalNumElements();
          if (constant_modes_are_global == false)
            Assert(n_relevant_rows == my_size,
                   ExcDimensionMismatch(n_relevant_rows, my_size));
          Assert(n_rows == static_cast<size_type>(
                             distributed_constant_modes->getGlobalLength()),
                 ExcDimensionMismatch(
                   n_rows, distributed_constant_modes->getGlobalLength()));

          auto vector_2d = distributed_constant_modes
                             ->template getLocalView<Kokkos::HostSpace>(
                               Tpetra::Access::ReadWriteStruct{});

          // Convert the given constant modes into a Tpetra multi-vector
          for (size_type d = 0; d < constant_modes_dimension; ++d)
            {
              for (size_type row = 0; row < my_size; ++row)
                {
                  const auto global_row_id =
                    constant_modes_are_global ?
                      domain_map->getGlobalElement(row) :
                      row;
                  vector_2d(row, d) =
                    static_cast<Number>(ad.constant_modes[d][global_row_id]);
                }
            }

          Teuchos::ParameterList &userParams =
            parameter_list.sublist("user data");
          userParams.set("Nullspace", distributed_constant_modes);
        }

      if (ad.output_details == true)
        parameter_list.set("verbosity", "high");
      else
        parameter_list.set("verbosity", "none");

      initialize(A, parameter_list);
    }

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
  }      // namespace TpetraWrappers
} // namespace LinearAlgebra

#endif // DEAL_II_TRILINOS_WITH_TPETRA_MUELU

DEAL_II_NAMESPACE_CLOSE

#endif
