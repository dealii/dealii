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
  class CellTypeWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    CellTypeWrapper(const CellTypeWrapper &other);

    /**
     * Constructor. Takes a CellKinds enum field and creates a Type class.
     */
    CellTypeWrapper(const ReferenceCell::Type::CellKinds &kind);

    /**
     * Constructor. Takes a CellKinds enum field and creates a Type class.
     */
    CellTypeWrapper(const ReferenceCell::Type &cell_type_in);

    /**
     * Constructor for an empty object.
     */
    CellTypeWrapper();

    /**
     * Destructor.
     */
    ~CellTypeWrapper();

    ReferenceCell::Type::CellKinds
    cell_kind() const;

  private:
    ReferenceCell::Type cell_type;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
