// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/trilinos_index_access.h>
#include <deal.II/lac/trilinos_precondition.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#  include <Epetra_MultiVector.h>
#  include <Ifpack.h>
#  include <Ifpack_Chebyshev.h>
#  include <Teuchos_ParameterList.hpp>
#  include <Teuchos_RCP.hpp>
#  include <ml_MultiLevelPreconditioner.h>
#  include <ml_include.h>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  /* -------------------------- PreconditionAMG -------------------------- */

  PreconditionAMG::AdditionalData::AdditionalData(
    const bool                            elliptic,
    const bool                            higher_order_elements,
    const unsigned int                    n_cycles,
    const bool                            w_cycle,
    const double                          aggregation_threshold,
    const std::vector<std::vector<bool>> &constant_modes,
    const unsigned int                    smoother_sweeps,
    const unsigned int                    smoother_overlap,
    const bool                            output_details,
    const char                           *smoother_type,
    const char                           *coarse_type)
    : elliptic(elliptic)
    , higher_order_elements(higher_order_elements)
    , n_cycles(n_cycles)
    , w_cycle(w_cycle)
    , aggregation_threshold(aggregation_threshold)
    , constant_modes(constant_modes)
    , smoother_sweeps(smoother_sweeps)
    , smoother_overlap(smoother_overlap)
    , output_details(output_details)
    , smoother_type(smoother_type)
    , coarse_type(coarse_type)
  {}



  void
  PreconditionAMG::AdditionalData::set_parameters(
    Teuchos::ParameterList              &parameter_list,
    std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
    const Epetra_RowMatrix              &matrix) const
  {
    if (elliptic == true)
      {
        ML_Epetra::SetDefaults("SA", parameter_list);

        // uncoupled mode can give a lot of warnings or even fail when there
        // are too many entries per row and aggregation gets complicated, but
        // MIS does not work if too few elements are located on one
        // processor. work around these warnings by choosing the different
        // strategies in different situations: for low order, always use the
        // standard choice uncoupled. if higher order, right now we also just
        // use Uncoupled, but we should be aware that maybe MIS might be
        // needed
        if (higher_order_elements)
          parameter_list.set("aggregation: type", "Uncoupled");
      }
    else
      {
        ML_Epetra::SetDefaults("NSSA", parameter_list);
        parameter_list.set("aggregation: type", "Uncoupled");
        parameter_list.set("aggregation: block scaling", true);
      }

    parameter_list.set("smoother: type", smoother_type);
    parameter_list.set("coarse: type", coarse_type);

    // Force re-initialization of the random seed to make ML deterministic
    parameter_list.set("initialize random seed", true);

    parameter_list.set("smoother: sweeps", static_cast<int>(smoother_sweeps));
    parameter_list.set("cycle applications", static_cast<int>(n_cycles));
    if (w_cycle == true)
      parameter_list.set("prec type", "MGW");
    else
      parameter_list.set("prec type", "MGV");

    parameter_list.set("smoother: Chebyshev alpha", 10.);
    parameter_list.set("smoother: ifpack overlap",
                       static_cast<int>(smoother_overlap));
    parameter_list.set("aggregation: threshold", aggregation_threshold);
    parameter_list.set("coarse: max size", 2000);

    if (output_details)
      parameter_list.set("ML output", 10);
    else
      parameter_list.set("ML output", 0);

    set_operator_null_space(parameter_list, distributed_constant_modes, matrix);
  }



  void
  PreconditionAMG::AdditionalData::set_operator_null_space(
    Teuchos::ParameterList              &parameter_list,
    std::unique_ptr<Epetra_MultiVector> &ptr_distributed_constant_modes,
    const Epetra_RowMatrix              &matrix) const
  {
    const auto run = [&](const auto &constant_modes) {
      const Epetra_Map &domain_map = matrix.OperatorDomainMap();

      const size_type constant_modes_dimension = constant_modes.size();
      ptr_distributed_constant_modes =
        std::make_unique<Epetra_MultiVector>(domain_map,
                                             constant_modes_dimension > 0 ?
                                               constant_modes_dimension :
                                               1);
      Assert(ptr_distributed_constant_modes, ExcNotInitialized());
      Epetra_MultiVector &distributed_constant_modes =
        *ptr_distributed_constant_modes;

      if (constant_modes_dimension > 0)
        {
          const size_type global_size = TrilinosWrappers::n_global_rows(matrix);
          Assert(global_size ==
                   static_cast<size_type>(TrilinosWrappers::global_length(
                     distributed_constant_modes)),
                 ExcDimensionMismatch(global_size,
                                      TrilinosWrappers::global_length(
                                        distributed_constant_modes)));
          const bool constant_modes_are_global =
            constant_modes[0].size() == global_size;
          const size_type my_size = domain_map.NumMyElements();

          // Reshape null space as a contiguous vector of doubles so that
          // Trilinos can read from it.
          const size_type expected_mode_size =
            constant_modes_are_global ? global_size : my_size;
          for (size_type d = 0; d < constant_modes_dimension; ++d)
            {
              Assert(constant_modes[d].size() == expected_mode_size,
                     ExcDimensionMismatch(constant_modes[d].size(),
                                          expected_mode_size));
              for (size_type row = 0; row < my_size; ++row)
                {
                  const TrilinosWrappers::types::int_type mode_index =
                    constant_modes_are_global ?
                      TrilinosWrappers::global_index(domain_map, row) :
                      row;
                  distributed_constant_modes[d][row] =
                    static_cast<double>(constant_modes[d][mode_index]);
                }
            }

          parameter_list.set("null space: type", "pre-computed");
          parameter_list.set("null space: dimension",
                             distributed_constant_modes.NumVectors());
          parameter_list.set("null space: vectors",
                             distributed_constant_modes.Values());
        }
    };

    if (!constant_modes_values.empty())
      {
        AssertDimension(constant_modes.size(), 0);
        run(constant_modes_values);
      }
    else
      run(constant_modes);
  }



  void
  PreconditionAMG::AdditionalData::set_parameters(
    Teuchos::ParameterList              &parameter_list,
    std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
    const SparseMatrix                  &matrix) const
  {
    return set_parameters(parameter_list,
                          distributed_constant_modes,
                          matrix.trilinos_matrix());
  }



  void
  PreconditionAMG::AdditionalData::set_operator_null_space(
    Teuchos::ParameterList              &parameter_list,
    std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
    const SparseMatrix                  &matrix) const
  {
    return set_operator_null_space(parameter_list,
                                   distributed_constant_modes,
                                   matrix.trilinos_matrix());
  }



  PreconditionAMG::~PreconditionAMG()
  {
    preconditioner.reset();
    trilinos_matrix.reset();
  }



  void
  PreconditionAMG::initialize(const SparseMatrix   &matrix,
                              const AdditionalData &additional_data)
  {
    initialize(matrix.trilinos_matrix(), additional_data);
  }



  void
  PreconditionAMG::initialize(const Epetra_RowMatrix &matrix,
                              const AdditionalData   &additional_data)
  {
    // Build the AMG preconditioner.
    Teuchos::ParameterList              ml_parameters;
    std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
    additional_data.set_parameters(ml_parameters,
                                   distributed_constant_modes,
                                   matrix);

    initialize(matrix, ml_parameters);

    if (additional_data.output_details)
      {
        ML_Epetra::MultiLevelPreconditioner *multilevel_operator =
          dynamic_cast<ML_Epetra::MultiLevelPreconditioner *>(
            preconditioner.get());
        Assert(multilevel_operator != nullptr,
               ExcMessage("Preconditioner setup failed."));
        multilevel_operator->PrintUnused(0);
      }
  }



  void
  PreconditionAMG::initialize(const SparseMatrix           &matrix,
                              const Teuchos::ParameterList &ml_parameters)
  {
    initialize(matrix.trilinos_matrix(), ml_parameters);
  }



  void
  PreconditionAMG::initialize(const Epetra_RowMatrix       &matrix,
                              const Teuchos::ParameterList &ml_parameters)
  {
    preconditioner.reset(
      new ML_Epetra::MultiLevelPreconditioner(matrix, ml_parameters));
  }



  template <typename number>
  void
  PreconditionAMG::initialize(
    const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
    const AdditionalData                 &additional_data,
    const double                          drop_tolerance,
    const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    preconditioner.reset();
    const size_type n_rows = deal_ii_sparse_matrix.m();

    // Init Epetra Matrix using an equidistributed map; avoid storing the
    // nonzero elements.
    IndexSet           distributor(n_rows);
    const unsigned int n_mpi_processes = communicator.NumProc();
    const unsigned int my_id           = communicator.MyPID();
    distributor.add_range(my_id * n_rows / n_mpi_processes,
                          (my_id + 1) * n_rows / n_mpi_processes);

    if (trilinos_matrix.get() == nullptr)
      trilinos_matrix = std::make_shared<SparseMatrix>();

    trilinos_matrix->reinit(distributor,
                            distributor,
                            deal_ii_sparse_matrix,
                            communicator.Comm(),
                            drop_tolerance,
                            true,
                            use_this_sparsity);

    initialize(*trilinos_matrix, additional_data);
  }



  void
  PreconditionAMG::reinit()
  {
    ML_Epetra::MultiLevelPreconditioner *multilevel_operator =
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner *>(preconditioner.get());
    multilevel_operator->ReComputePreconditioner();
  }



  void
  PreconditionAMG::clear()
  {
    PreconditionBase::clear();
    trilinos_matrix.reset();
  }



  PreconditionAMG::size_type
  PreconditionAMG::memory_consumption() const
  {
    unsigned int memory = sizeof(*this);

    // todo: find a way to read out ML's data
    // sizes
    if (trilinos_matrix.get() != nullptr)
      memory += trilinos_matrix->memory_consumption();
    return memory;
  }



#  ifndef DOXYGEN
  // explicit instantiations
  template void
  PreconditionAMG::initialize(const ::dealii::SparseMatrix<double> &,
                              const AdditionalData &,
                              const double,
                              const ::dealii::SparsityPattern *);
  template void
  PreconditionAMG::initialize(const ::dealii::SparseMatrix<float> &,
                              const AdditionalData &,
                              const double,
                              const ::dealii::SparsityPattern *);
#  endif



} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
