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

#ifndef dealii_reference_cell_wrapper_h
#define dealii_reference_cell_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/grid/reference_cell.h>

#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class ReferenceCellWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    ReferenceCellWrapper(const ReferenceCellWrapper &other);

    /**
     * Constructor. Takes a cell kind id field and creates a Type class.
     */
    ReferenceCellWrapper(const std::uint8_t &kind);

    /**
     * Constructor. Takes a ReferenceCell object and creates a Type class.
     */
    ReferenceCellWrapper(const ReferenceCell &cell_type_in);

    /**
     * Constructor for an empty object.
     */
    ReferenceCellWrapper();

    /**
     * Destructor.
     */
    ~ReferenceCellWrapper();

    std::uint8_t
    cell_kind() const;

  private:
    ReferenceCell cell_type;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
