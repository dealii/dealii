// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_reference_cell_wrapper_h
#define dealii_reference_cell_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/grid/reference_cell.h>

#include <boost/python.hpp>

#include <variant>

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
    template <int dim>
    ReferenceCellWrapper(const ReferenceCell<dim> &cell_type_in);

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
    /**
     * A variable that stores the reference cell this object was initialized
     * with -- regardless of what dimension the object has.
     */
    std::variant<ReferenceCell<0>,
                 ReferenceCell<1>,
                 ReferenceCell<2>,
                 ReferenceCell<3>>
      cell_type;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
