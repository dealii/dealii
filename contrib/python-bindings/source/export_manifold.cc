// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#include <manifold_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(create_cylindrical_overloads,
                                         create_cylindrical,
                                         0,
                                         2)


  const char create_spherical_docstring[] =
    " Create spherical manifold with a given center point.      \n";


  const char create_polar_docstring[] =
    " Create polar manifold with a given center point.          \n";


  const char create_function_docstring[] =
    " Create manifold with the given python push forward and    \n"
    " pull back functions.                                      \n";


  const char create_function_string_docstring[] =
    " Create manifold with given string expression for the push \n"
    " forward and pull back functions.                          \n";


  const char create_cylindrical_fixed_docstring[] =
    " Create cylindrical manifold oriented along a given axis   \n"
    " (0 - x, 1 - y, 2 - z).                                    \n";


  const char create_cylindrical_direction_docstring[] =
    " Create cylindrical manifold with an axis that points in   \n"
    " direction direction and goes through the given point on   \n"
    " axis.                                                     \n";


  void
  export_manifold()
  {
    boost::python::class_<ManifoldWrapper>(
      "Manifold",
      boost::python::init<const int, const int>(
        boost::python::args("dim", "spacedim")))
      .def("create_spherical",
           &ManifoldWrapper::create_spherical,
           create_spherical_docstring,
           boost::python::args("self", "center"))
      .def("create_polar",
           &ManifoldWrapper::create_polar,
           create_polar_docstring,
           boost::python::args("self", "center"))
      .def("create_cylindrical",
           static_cast<void (ManifoldWrapper::*)(const int, const double)>(
             &ManifoldWrapper::create_cylindrical),
           create_cylindrical_overloads(
             boost::python::args("self", "axis", "tolerance"),
             create_cylindrical_fixed_docstring))
      .def("create_cylindrical",
           static_cast<void (ManifoldWrapper::*)(const boost::python::list &,
                                                 const boost::python::list &)>(
             &ManifoldWrapper::create_cylindrical),
           create_cylindrical_direction_docstring,
           boost::python::args("self", "direction", "axial_point"))
      .def("create_function",
           &ManifoldWrapper::create_function,
           create_function_docstring,
           boost::python::args("self", "push_forward", "pull_back"))
      .def("create_function_string",
           &ManifoldWrapper::create_function_string,
           create_function_string_docstring,
           boost::python::args("self", "push_forward", "pull_back"));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
