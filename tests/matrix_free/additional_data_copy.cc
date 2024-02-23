// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test copy constructor and copy assignment of MatrixFree::AdditionalData

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <typename T>
bool
operator==(const T &lhs, const T &rhs)
{
  return lhs.tasks_parallel_scheme == rhs.tasks_parallel_scheme &&
         lhs.tasks_block_size == rhs.tasks_block_size &&
         lhs.mapping_update_flags == rhs.mapping_update_flags &&
         lhs.mapping_update_flags_boundary_faces ==
           rhs.mapping_update_flags_boundary_faces &&
         lhs.mapping_update_flags_inner_faces ==
           rhs.mapping_update_flags_inner_faces &&
         lhs.mapping_update_flags_faces_by_cells ==
           rhs.mapping_update_flags_faces_by_cells &&
         lhs.mg_level == rhs.mg_level &&
         lhs.store_plain_indices == rhs.store_plain_indices &&
         lhs.initialize_indices == rhs.initialize_indices &&
         lhs.initialize_mapping == rhs.initialize_mapping &&
         lhs.overlap_communication_computation ==
           rhs.overlap_communication_computation &&
         lhs.hold_all_faces_to_owned_cells ==
           rhs.hold_all_faces_to_owned_cells &&
         lhs.cell_vectorization_categories_strict ==
           rhs.cell_vectorization_categories_strict &&
         lhs.cell_vectorization_category == rhs.cell_vectorization_category;
}

int
main()
{
  initlog();

  using AD =
    typename MatrixFree<2, double, VectorizedArray<double>>::AdditionalData;

  AD ad;
  ad.tasks_parallel_scheme = AD::TasksParallelScheme::partition_color;
  ad.tasks_block_size      = 100;
  ad.mapping_update_flags  = update_values;
  ad.mapping_update_flags_boundary_faces  = update_values;
  ad.mapping_update_flags_inner_faces     = update_values;
  ad.mapping_update_flags_faces_by_cells  = update_values;
  ad.store_plain_indices                  = false;
  ad.initialize_indices                   = false;
  ad.initialize_mapping                   = false;
  ad.overlap_communication_computation    = false;
  ad.hold_all_faces_to_owned_cells        = true;
  ad.cell_vectorization_categories_strict = true;
  ad.cell_vectorization_category          = {1, 2, 3};

  {
    // copy constructor
    AD copy(ad);
    Assert(copy == ad, ExcInternalError());
  }

  {
    // copy assignment
    AD copy;
    copy = ad;
    Assert(copy == ad, ExcInternalError());
  }
}
