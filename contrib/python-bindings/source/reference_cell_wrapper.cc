// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/reference_cell.h>

#include <reference_cell_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  ReferenceCellWrapper::ReferenceCellWrapper()
  {}



  ReferenceCellWrapper::ReferenceCellWrapper(const ReferenceCellWrapper &other)
  {
    cell_type = other.cell_type;
  }



  ReferenceCellWrapper::ReferenceCellWrapper(const std::uint8_t &kind)
  {
    cell_type = internal::make_reference_cell_from_int(kind);
  }



  ReferenceCellWrapper::ReferenceCellWrapper(const ReferenceCell &cell_type_in)
  {
    cell_type = cell_type_in;
  }



  ReferenceCellWrapper::~ReferenceCellWrapper()
  {}



  std::uint8_t
  ReferenceCellWrapper::cell_kind() const
  {
    return cell_type;
  }

} // namespace python

DEAL_II_NAMESPACE_CLOSE
