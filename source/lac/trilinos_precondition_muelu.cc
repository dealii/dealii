// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2019 by the deal.II authors
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

#include <deal.II/lac/trilinos_index_access.h>
#include <deal.II/lac/trilinos_precondition.h>

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_TRILINOS_WITH_MUELU
#    include <deal.II/lac/sparse_matrix.h>
#    include <deal.II/lac/trilinos_sparse_matrix.h>

#    include <MueLu_CreateEpetraPreconditioner.hpp>
#    include <ml_MultiLevelPreconditioner.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  PreconditionAMGMueLu::AdditionalData::AdditionalData(
    const bool                            elliptic,
    const unsigned int                    n_cycles,
    const bool                            w_cycle,
    const double                          aggregation_threshold,
    const std::vector<std::vector<bool>> &constant_modes,
    const unsigned int                    smoother_sweeps,
    const unsigned int                    smoother_overlap,
    const bool                            output_details,
    const char *                          smoother_type,
    const char *                          coarse_type)
    : elliptic(elliptic)
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



  PreconditionAMGMueLu::PreconditionAMGMueLu()
  {
    // clang-tidy wants to default the constructor if we disable the check
    // in case we compile without 64-bit indices
#    ifdef DEAL_II_WITH_64BIT_INDICES
    constexpr bool enabled = false;
#    else
    constexpr bool enabled = true;
#    endif
    AssertThrow(enabled,
                ExcMessage(
                  "PreconditionAMGMueLu does not support 64bit-indices!"));
  }



  void
  PreconditionAMGMueLu::initialize(const SparseMatrix &  matrix,
                                   const AdditionalData &additional_data)
  {
    initialize(matrix.trilinos_matrix(), additional_data);
  }



  void
  PreconditionAMGMueLu::initialize(const Epetra_CrsMatrix &matrix,
                                   const AdditionalData &  additional_data)
  {
    // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;

    parameter_list.set("parameterlist: syntax", "ml");

    if (additional_data.elliptic == true)
      ML_Epetra::SetDefaults("SA", parameter_list);
    else
      {
        ML_Epetra::SetDefaults("NSSA", parameter_list);
        parameter_list.set("aggregation: block scaling", true);
      }
    // MIS does not exist anymore, only choice are uncoupled and coupled. When
    // using uncoupled, aggregates cannot span multiple processes. When using
    // coupled aggregates can span multiple processes.
    parameter_list.set("aggregation: type", "Uncoupled");

    parameter_list.set("smoother: type", additional_data.smoother_type);
    parameter_list.set("coarse: type", additional_data.coarse_type);

    parameter_list.set("smoother: sweeps",
                       static_cast<int>(additional_data.smoother_sweeps));
    parameter_list.set("cycle applications",
                       static_cast<int>(additional_data.n_cycles));
    if (additional_data.w_cycle == true)
      parameter_list.set("prec type", "MGW");
    else
      parameter_list.set("prec type", "MGV");

    parameter_list.set("smoother: Chebyshev alpha", 10.);
    parameter_list.set("smoother: ifpack overlap",
                       static_cast<int>(additional_data.smoother_overlap));
    parameter_list.set("aggregation: threshold",
                       additional_data.aggregation_threshold);
    parameter_list.set("coarse: max size", 2000);

    if (additional_data.output_details)
      parameter_list.set("ML output", 10);
    else
      parameter_list.set("ML output", 0);

    const Epetra_Map &domain_map = matrix.OperatorDomainMap();

    const size_type constant_modes_dimension =
      additional_data.constant_modes.size();
    Epetra_MultiVector distributed_constant_modes(
      domain_map, constant_modes_dimension > 0 ? constant_modes_dimension : 1);
    std::vector<double> dummy(constant_modes_dimension);

    if (constant_modes_dimension > 0)
      {
        const size_type n_rows = TrilinosWrappers::n_global_rows(matrix);
        const bool      constant_modes_are_global =
          additional_data.constant_modes[0].size() == n_rows;
        const size_type n_relevant_rows =
          constant_modes_are_global ? n_rows :
                                      additional_data.constant_modes[0].size();
        const size_type my_size = domain_map.NumMyElements();
        if (constant_modes_are_global == false)
          Assert(n_relevant_rows == my_size,
                 ExcDimensionMismatch(n_relevant_rows, my_size));
        Assert(n_rows == static_cast<size_type>(TrilinosWrappers::global_length(
                           distributed_constant_modes)),
               ExcDimensionMismatch(n_rows,
                                    TrilinosWrappers::global_length(
                                      distributed_constant_modes)));

        (void)n_relevant_rows;
        (void)global_length;

        // Reshape null space as a contiguous vector of doubles so that
        // Trilinos can read from it.
        for (size_type d = 0; d < constant_modes_dimension; ++d)
          for (size_type row = 0; row < my_size; ++row)
            {
              TrilinosWrappers::types::int_type global_row_id =
                constant_modes_are_global ?
                  TrilinosWrappers::global_index(domain_map, row) :
                  row;
              distributed_constant_modes[d][row] =
                additional_data.constant_modes[d][global_row_id];
            }

        parameter_list.set("null space: type", "pre-computed");
        parameter_list.set("null space: dimension",
                           distributed_constant_modes.NumVectors());
        if (my_size > 0)
          parameter_list.set("null space: vectors",
                             distributed_constant_modes.Values());
        // We need to set a valid pointer to data even if there is no data on
        // the current processor. Therefore, pass a dummy in that case
        else
          parameter_list.set("null space: vectors", dummy.data());
      }

    initialize(matrix, parameter_list);
  }



  void
  PreconditionAMGMueLu::initialize(const SparseMatrix &    matrix,
                                   Teuchos::ParameterList &muelu_parameters)
  {
    initialize(matrix.trilinos_matrix(), muelu_parameters);
  }



  void
  PreconditionAMGMueLu::initialize(const Epetra_CrsMatrix &matrix,
                                   Teuchos::ParameterList &muelu_parameters)
  {
    const auto teuchos_wrapped_matrix =
      Teuchos::rcp(const_cast<Epetra_CrsMatrix *>(&matrix), false);
    preconditioner = MueLu::CreateEpetraPreconditioner(teuchos_wrapped_matrix,
                                                       muelu_parameters);
  }



  template <typename number>
  void
  PreconditionAMGMueLu::initialize(
    const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
    const AdditionalData &                additional_data,
    const double                          drop_tolerance,
    const ::dealii::SparsityPattern *     use_this_sparsity)
  {
    preconditioner.reset();
    const size_type n_rows = deal_ii_sparse_matrix.m();

    // Init Epetra Matrix using an
    // equidistributed map; avoid
    // storing the nonzero
    // elements.
    vector_distributor = std::make_shared<Epetra_Map>(
      static_cast<TrilinosWrappers::types::int_type>(n_rows), 0, communicator);

    if (trilinos_matrix.get() == nullptr)
      trilinos_matrix = std::make_shared<SparseMatrix>();

    trilinos_matrix->reinit(*vector_distributor,
                            *vector_distributor,
                            deal_ii_sparse_matrix,
                            drop_tolerance,
                            true,
                            use_this_sparsity);

    initialize(*trilinos_matrix, additional_data);
  }



  void
  PreconditionAMGMueLu::clear()
  {
    PreconditionBase::clear();
    trilinos_matrix.reset();
  }



  PreconditionAMGMueLu::size_type
  PreconditionAMGMueLu::memory_consumption() const
  {
    unsigned int memory = sizeof(*this);

    // todo: find a way to read out ML's data
    // sizes
    if (trilinos_matrix.get() != nullptr)
      memory += trilinos_matrix->memory_consumption();
    return memory;
  }



#    ifndef DOXYGEN
  // explicit instantiations
  template void
  PreconditionAMGMueLu::initialize(const ::dealii::SparseMatrix<double> &,
                                   const AdditionalData &,
                                   const double,
                                   const ::dealii::SparsityPattern *);
  template void
  PreconditionAMGMueLu::initialize(const ::dealii::SparseMatrix<float> &,
                                   const AdditionalData &,
                                   const double,
                                   const ::dealii::SparsityPattern *);
#    endif

} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_MUELU
#endif   // DEAL_II_WITH_TRILINOS
