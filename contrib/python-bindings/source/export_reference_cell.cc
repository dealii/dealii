// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#include <boost/python.hpp>

#include <reference_cell_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  const char cell_kind_docstring[] =
    "Return geometric type of the corresponding reference cell.               \n";


  void
  export_reference_cell()
  {
    boost::python::enum_<ReferenceCell::Type::CellKinds>("CellKinds")
      .value("Vertex", ReferenceCell::Type::Vertex)
      .value("Line", ReferenceCell::Type::Line)
      .value("Tri", ReferenceCell::Type::Tri)
      .value("Quad", ReferenceCell::Type::Quad)
      .value("Tet", ReferenceCell::Type::Tet)
      .value("Pyramid", ReferenceCell::Type::Pyramid)
      .value("Wedge", ReferenceCell::Type::Wedge)
      .value("Hex", ReferenceCell::Type::Hex)
      .value("Invalid", ReferenceCell::Type::Invalid);

    boost::python::class_<CellTypeWrapper>(
      "CellType",
      boost::python::init<ReferenceCell::Type::CellKinds>(
        boost::python::args("kind")))
      .add_property("kind", &CellTypeWrapper::cell_kind, cell_kind_docstring);
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
