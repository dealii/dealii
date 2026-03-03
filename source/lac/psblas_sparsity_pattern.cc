// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/index_set.h>
#include <deal.II/base/init_finalize.h>

#include "deal.II/lac/psblas_common.h"

#include <psb_c_base.h>


#ifdef DEAL_II_WITH_PSBLAS
#  include <deal.II/lac/psblas_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{
  // SparsityPattern
  SparsityPattern::SparsityPattern(const IndexSet &index_set,
                                   const MPI_Comm  communicator)
  {
    SparsityPatternBase::resize(index_set.size(), index_set.size());

    Assert(communicator != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to SparseMatrix::reinit()."));

    psblas_descriptor.reset(psb_c_new_descriptor(),
                            PSCToolkitWrappers::internal::DescriptorDeleter());

    // Use get_index_vector() from IndexSet to get the indexes
    const std::vector<types::global_dof_index> &indexes =
      index_set.get_index_vector();

    psb_i_t number_of_local_indexes = indexes.size(); // Number of local indexes
    // Copy the indexes into a psb_l_t vector
    std::vector<psb_l_t> vl(number_of_local_indexes);
    for (psb_i_t i = 0; i < number_of_local_indexes; ++i)
      {
        const auto psblas_index = static_cast<psb_l_t>(indexes[i]);
        AssertIntegerConversion(psblas_index, indexes[i]);
        vl[i] = psblas_index;
      }

    // Insert the indexes into the descriptor
    psb_c_cdall_vl(number_of_local_indexes,
                   vl.data(),
                   *InitFinalize::get_psblas_context(),
                   psblas_descriptor.get());
  }



  void
  SparsityPattern::add(const PSCToolkitWrappers::SparsityPattern::size_type i,
                       const PSCToolkitWrappers::SparsityPattern::size_type j)
  {
    add_entries(i, &j, &j + 1);
  }



  template <typename ForwardIterator>
  inline void
  SparsityPattern::add_entries(
    const PSCToolkitWrappers::SparsityPattern::size_type row,
    ForwardIterator                                      begin,
    ForwardIterator                                      end,
    const bool                                           indices_are_sorted)
  {
    Assert(psblas_descriptor.get() != nullptr,
           ExcMessage("PSBLAS descriptor is null."));
    if (begin == end)
      return;

    (void)indices_are_sorted;
    psb_i_t              nz = static_cast<int>(end - begin);
    std::vector<psb_l_t> ia(nz);
    std::vector<psb_l_t> ja(nz);

    for (int k = 0; k < nz; ++k)
      {
        ia[k] = row;          // row index
        ja[k] = *(begin + k); // column index
      }
    int err = psb_c_cdins(nz, ia.data(), ja.data(), psblas_descriptor.get());
    Assert(err == 0, ExcInsertionInPSBLASMatrix(err));
  }



  void
  SparsityPattern::add_row_entries(const size_type                  &row,
                                   const ArrayView<const size_type> &columns,
                                   const bool indices_are_sorted)
  {
    add_entries(row, columns.begin(), columns.end(), indices_are_sorted);
  }



  void
  SparsityPattern::compress()
  {
    Assert(psblas_descriptor.get() != nullptr,
           ExcMessage("PSBLAS descriptor is null."));
    int err = -1;
    if (!psb_c_cd_is_asb(psblas_descriptor.get()))
      {
        err = psb_c_cdasb(psblas_descriptor.get());
        Assert(
          err == 0,
          ExcCallingPSBLASFunction(
            err,
            "Error while finalizing SparsityPattern through psb_c_cdasb."));
      }
  }
} // namespace PSCToolkitWrappers

DEAL_II_NAMESPACE_CLOSE
#endif
