// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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
    // We need to convert the given integer value to a ReferenceCell object.
    // The issue is that we need to know the dimension of the cell in order to
    // do so, but this is not given as an argument. However, we can infer it
    // from the integer value, since the integer values for the different cell
    // types are defined in such a way that they are unique across dimensions.
    // We can thus just call a function that converts the integer value to a
    // ReferenceCell object, and this function will internally determine the
    // dimension of the cell type based on the integer value.
    switch (kind)
      {
        case ReferenceCells::Vertex:
          cell_type = internal::make_reference_cell_from_int<0>(kind);
          break;
        case ReferenceCells::Line:
          cell_type = internal::make_reference_cell_from_int<1>(kind);
          break;
        case ReferenceCells::Triangle:
        case ReferenceCells::Quadrilateral:
          cell_type = internal::make_reference_cell_from_int<2>(kind);
          break;
        case ReferenceCells::Tetrahedron:
        case ReferenceCells::Hexahedron:
        case ReferenceCells::Wedge:
        case ReferenceCells::Pyramid:
          cell_type = internal::make_reference_cell_from_int<3>(kind);
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
      }
  }



  template <int dim>
  ReferenceCellWrapper::ReferenceCellWrapper(
    const ReferenceCell<dim> &cell_type_in)
  {
    cell_type = cell_type_in;
  }



  ReferenceCellWrapper::~ReferenceCellWrapper()
  {}



  std::uint8_t
  ReferenceCellWrapper::cell_kind() const
  {
    // One of the four variants in 'cell_type' is always active, so we can just
    // return the one that is active.
    switch (cell_type.index())
      {
        case 0:
          return std::get<0>(cell_type);
        case 1:
          return std::get<1>(cell_type);
        case 2:
          return std::get<2>(cell_type);
        case 3:
          return std::get<3>(cell_type);
        default:
          DEAL_II_ASSERT_UNREACHABLE();
          return 0; // silence compiler warning
      }
  }


  template ReferenceCellWrapper::ReferenceCellWrapper<0>(
    const ReferenceCell<0> &);
  template ReferenceCellWrapper::ReferenceCellWrapper<1>(
    const ReferenceCell<1> &);
  template ReferenceCellWrapper::ReferenceCellWrapper<2>(
    const ReferenceCell<2> &);
  template ReferenceCellWrapper::ReferenceCellWrapper<3>(
    const ReferenceCell<3> &);


} // namespace python

DEAL_II_NAMESPACE_CLOSE
