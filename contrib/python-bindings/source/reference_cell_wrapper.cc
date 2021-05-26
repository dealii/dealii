// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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
    cell_type = internal::ReferenceCell::make_reference_cell_from_int(kind);
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
