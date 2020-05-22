// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.templates.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    template struct ConstraintValues<double>;
    template struct ConstraintValues<float>;

    template void
    DoFInfo::read_dof_indices<double>(
      const std::vector<types::global_dof_index> &,
      const std::vector<unsigned int> &,
      const dealii::AffineConstraints<double> &,
      const unsigned int,
      ConstraintValues<double> &,
      bool &);

    template void
    DoFInfo::read_dof_indices<float>(
      const std::vector<types::global_dof_index> &,
      const std::vector<unsigned int> &,
      const dealii::AffineConstraints<float> &,
      const unsigned int,
      ConstraintValues<double> &,
      bool &);

    template void
    DoFInfo::compute_face_index_compression<1>(
      const std::vector<FaceToCellTopology<1>> &);
    template void
    DoFInfo::compute_face_index_compression<2>(
      const std::vector<FaceToCellTopology<2>> &);
    template void
    DoFInfo::compute_face_index_compression<4>(
      const std::vector<FaceToCellTopology<4>> &);
    template void
    DoFInfo::compute_face_index_compression<8>(
      const std::vector<FaceToCellTopology<8>> &);
    template void
    DoFInfo::compute_face_index_compression<16>(
      const std::vector<FaceToCellTopology<16>> &);

    template void
    DoFInfo::compute_vector_zero_access_pattern<1>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<1>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<2>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<2>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<4>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<4>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<8>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<8>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<16>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<16>> &);

    template void
    DoFInfo::print_memory_consumption<std::ostream>(std::ostream &,
                                                    const TaskInfo &) const;
    template void
    DoFInfo::print_memory_consumption<ConditionalOStream>(
      ConditionalOStream &,
      const TaskInfo &) const;

    template void
    DoFInfo::print<double>(const std::vector<double> &,
                           const std::vector<unsigned int> &,
                           std::ostream &) const;

    template void
    DoFInfo::print<float>(const std::vector<float> &,
                          const std::vector<unsigned int> &,
                          std::ostream &) const;
  } // namespace MatrixFreeFunctions
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
