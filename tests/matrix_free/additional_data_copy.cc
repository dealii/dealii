// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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
         lhs.level_mg_handler == rhs.level_mg_handler &&
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
  ad.mg_level                             = 10;
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
