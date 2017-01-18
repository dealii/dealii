// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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

#include <deal.II/lac/trilinos_precondition.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Ifpack.h>
#  include <Ifpack_Chebyshev.h>
#  include <Teuchos_ParameterList.hpp>
#  include <Teuchos_RCP.hpp>
#  include <Epetra_MultiVector.h>
#  include <ml_include.h>
#  include <ml_MultiLevelPreconditioner.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    int n_global_rows (const Epetra_RowMatrix &matrix)
    {
      return matrix.NumGlobalRows();
    }

    int global_length (const Epetra_MultiVector &vector)
    {
      return vector.GlobalLength();
    }

    int gid(const Epetra_Map &map, unsigned int i)
    {
      return map.GID(i);
    }
#else
    long long int n_global_rows (const Epetra_RowMatrix &matrix)
    {
      return matrix.NumGlobalRows64();
    }

    long long int global_length (const Epetra_MultiVector &vector)
    {
      return vector.GlobalLength64();
    }

    long long int gid(const Epetra_Map &map, dealii::types::global_dof_index i)
    {
      return map.GID64(i);
    }
#endif
  }



  /* -------------------------- PreconditionAMG -------------------------- */

  PreconditionAMG::AdditionalData::
  AdditionalData (const bool                             elliptic,
                  const bool                             higher_order_elements,
                  const unsigned int                     n_cycles,
                  const bool                             w_cycle,
                  const double                           aggregation_threshold,
                  const std::vector<std::vector<bool> > &constant_modes,
                  const unsigned int                     smoother_sweeps,
                  const unsigned int                     smoother_overlap,
                  const bool                             output_details,
                  const char                            *smoother_type,
                  const char                            *coarse_type)
    :
    elliptic (elliptic),
    higher_order_elements (higher_order_elements),
    n_cycles (n_cycles),
    w_cycle (w_cycle),
    aggregation_threshold (aggregation_threshold),
    constant_modes (constant_modes),
    smoother_sweeps (smoother_sweeps),
    smoother_overlap (smoother_overlap),
    output_details (output_details),
    smoother_type (smoother_type),
    coarse_type (coarse_type)
  {}


  PreconditionAMG::~PreconditionAMG()
  {
    preconditioner.reset();
    trilinos_matrix.reset();
  }



  void
  PreconditionAMG::initialize (const SparseMatrix   &matrix,
                               const AdditionalData &additional_data)
  {
    initialize(matrix.trilinos_matrix(), additional_data);
  }



  void
  PreconditionAMG::initialize (const Epetra_RowMatrix &matrix,
                               const AdditionalData   &additional_data)
  {
    // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;

    if (additional_data.elliptic == true)
      {
        ML_Epetra::SetDefaults("SA",parameter_list);

        // uncoupled mode can give a lot of warnings or even fail when there
        // are too many entries per row and aggreggation gets complicated, but
        // MIS does not work if too few elements are located on one
        // processor. work around these warnings by choosing the different
        // strategies in different situations: for low order, always use the
        // standard choice uncoupled. if higher order, right now we also just
        // use Uncoupled, but we should be aware that maybe MIS might be
        // needed
        if (additional_data.higher_order_elements)
          parameter_list.set("aggregation: type", "Uncoupled");
      }
    else
      {
        ML_Epetra::SetDefaults("NSSA",parameter_list);
        parameter_list.set("aggregation: type", "Uncoupled");
        parameter_list.set("aggregation: block scaling", true);
      }

    parameter_list.set("smoother: type", additional_data.smoother_type);
    parameter_list.set("coarse: type", additional_data.coarse_type);

    // Force re-initialization of the random seed to make ML deterministic
    // (only supported in trilinos >12.2):
#if DEAL_II_TRILINOS_VERSION_GTE(12,4,0)
    parameter_list.set("initialize random seed", true);
#endif

    parameter_list.set("smoother: sweeps",
                       static_cast<int>(additional_data.smoother_sweeps));
    parameter_list.set("cycle applications",
                       static_cast<int>(additional_data.n_cycles));
    if (additional_data.w_cycle == true)
      parameter_list.set("prec type", "MGW");
    else
      parameter_list.set("prec type", "MGV");

    parameter_list.set("smoother: Chebyshev alpha",10.);
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
    Epetra_MultiVector distributed_constant_modes (domain_map,
                                                   constant_modes_dimension > 0 ?
                                                   constant_modes_dimension : 1);
    std::vector<double> dummy (constant_modes_dimension);

    if (constant_modes_dimension > 0)
      {
        const size_type global_size = n_global_rows(matrix);
        (void)global_length; // work around compiler warning about unused function in release mode
        Assert (global_size ==
                static_cast<size_type>(global_length(distributed_constant_modes)),
                ExcDimensionMismatch(global_size,
                                     global_length(distributed_constant_modes)));
        const bool constant_modes_are_global
          = additional_data.constant_modes[0].size() == global_size;
        const size_type my_size = domain_map.NumMyElements();

        // Reshape null space as a contiguous vector of doubles so that
        // Trilinos can read from it.
        const size_type expected_mode_size =
          constant_modes_are_global ? global_size : my_size;
        for (size_type d=0; d<constant_modes_dimension; ++d)
          {
            Assert (additional_data.constant_modes[d].size() == expected_mode_size,
                    ExcDimensionMismatch(additional_data.constant_modes[d].size(), expected_mode_size));
            for (size_type row=0; row<my_size; ++row)
              {
                const TrilinosWrappers::types::int_type mode_index =
                  constant_modes_are_global ? gid(domain_map,row) : row;
                distributed_constant_modes[d][row] =
                  additional_data.constant_modes[d][mode_index];
              }
          }
        (void)expected_mode_size;

        parameter_list.set("null space: type", "pre-computed");
        parameter_list.set("null space: dimension",
                           distributed_constant_modes.NumVectors());
        if (my_size > 0)
          parameter_list.set("null space: vectors",
                             distributed_constant_modes.Values());
        // We need to set a valid pointer to data even if there is no data on
        // the current processor. Therefore, pass a dummy in that case
        else
          parameter_list.set("null space: vectors",
                             &dummy[0]);
      }

    initialize (matrix, parameter_list);

    if (additional_data.output_details)
      {
        ML_Epetra::MultiLevelPreconditioner *multilevel_operator =
          dynamic_cast<ML_Epetra::MultiLevelPreconditioner *> (preconditioner.get());
        Assert (multilevel_operator != 0,
                ExcMessage ("Preconditioner setup failed."));
        multilevel_operator->PrintUnused(0);
      }
  }



  void
  PreconditionAMG::initialize (const SparseMatrix           &matrix,
                               const Teuchos::ParameterList &ml_parameters)
  {
    initialize(matrix.trilinos_matrix(), ml_parameters);
  }



  void
  PreconditionAMG::initialize (const Epetra_RowMatrix       &matrix,
                               const Teuchos::ParameterList &ml_parameters)
  {
    preconditioner.reset ();
    preconditioner.reset (new ML_Epetra::MultiLevelPreconditioner
                          (matrix, ml_parameters));
  }



  template <typename number>
  void
  PreconditionAMG::
  initialize (const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
              const AdditionalData                 &additional_data,
              const double                          drop_tolerance,
              const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    preconditioner.reset();
    const size_type n_rows = deal_ii_sparse_matrix.m();

    // Init Epetra Matrix using an
    // equidistributed map; avoid
    // storing the nonzero
    // elements.
    vector_distributor.reset (new Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(n_rows),
                                             0, communicator));

    if (trilinos_matrix.get() == 0)
      trilinos_matrix.reset (new SparseMatrix());

    trilinos_matrix->reinit (*vector_distributor, *vector_distributor,
                             deal_ii_sparse_matrix, drop_tolerance, true,
                             use_this_sparsity);

    initialize (*trilinos_matrix, additional_data);
  }



  void PreconditionAMG::reinit ()
  {
    ML_Epetra::MultiLevelPreconditioner *multilevel_operator =
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner *> (preconditioner.get());
    multilevel_operator->ReComputePreconditioner();
  }



  void PreconditionAMG::clear ()
  {
    PreconditionBase::clear();
    trilinos_matrix.reset();
  }



  PreconditionAMG::size_type
  PreconditionAMG::memory_consumption() const
  {
    unsigned int memory = sizeof(this);

    // todo: find a way to read out ML's data
    // sizes
    if (trilinos_matrix.get() != 0)
      memory += trilinos_matrix->memory_consumption();
    return memory;
  }




  // explicit instantiations
  template void PreconditionAMG::initialize (const ::dealii::SparseMatrix<double> &,
                                             const AdditionalData &, const double,
                                             const ::dealii::SparsityPattern *);
  template void PreconditionAMG::initialize (const ::dealii::SparseMatrix<float> &,
                                             const AdditionalData &, const double,
                                             const ::dealii::SparsityPattern *);






}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
