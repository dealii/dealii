// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
