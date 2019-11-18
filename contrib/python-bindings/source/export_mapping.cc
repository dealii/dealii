// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#include <mapping_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  const char transform_unit_to_real_cell_docstring[] =
    " Map the point p on the unit cell to the corresponding point       \n"
    " on the real cell.                                                 \n";

  const char transform_real_to_unit_cell_docstring[] =
    " Map the point p on the real cell to the corresponding point       \n"
    " on the unit cell.                                                 \n";



  void
  export_mapping()
  {
    boost::python::class_<MappingQGenericWrapper>(
      "MappingQGeneric",
      boost::python::init<const int, const int, const int>(
        boost::python::args("dim", "spacedim", "degree")))
      .def("transform_real_to_unit_cell",
           &MappingQGenericWrapper::transform_real_to_unit_cell,
           transform_real_to_unit_cell_docstring,
           boost::python::args("self", "cell_wrapper", "point_wrapper"))
      .def("transform_unit_to_real_cell",
           &MappingQGenericWrapper::transform_unit_to_real_cell,
           transform_unit_to_real_cell_docstring,
           boost::python::args("self", "cell_wrapper", "point_wrapper"));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
