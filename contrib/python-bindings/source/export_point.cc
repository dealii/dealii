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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <boost/python.hpp>

#include <point_wrapper.h>

// clang complains about explicitly assigning boost::python::self to itself
// below. However, this is the correct way to define the python bindings.
#if defined(__clang__) && __clang_major__ > 6
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif

DEAL_II_NAMESPACE_OPEN

namespace python
{
  const char to_list_docstring[] =
    "Convert point's coordinates to a python list with [x,y] or [x,y,z]     \n"
    "for 2-D or 3-D, respectively.                                          \n";



  const char distance_docstring[] =
    "Return the Euclidean distance of this point to the point p             \n";



  const char norm_docstring[] =
    "Return the L2 norm of the vector connecting the origin to the point    \n";



  const char norm_square_docstring[] =
    "Return the sum of the absolute squares of all entries                  \n";



  const char get_x_docstring[] =
    "Get the x component of the point                                       \n";



  const char get_y_docstring[] =
    "Get the y component of the point                                       \n";



  const char get_z_docstring[] =
    "Get the z component of the point                                       \n";



  void
  export_point()
  {
    boost::python::class_<PointWrapper>(
      "Point", boost::python::init<boost::python::list>())
      .def("to_list",
           &PointWrapper::to_list,
           to_list_docstring,
           boost::python::args("self"))
      .def("distance",
           &PointWrapper::distance,
           distance_docstring,
           boost::python::args("self", "p"))
      .def("norm",
           &PointWrapper::norm,
           norm_docstring,
           boost::python::args("self"))
      .def("norm_square",
           &PointWrapper::norm_square,
           norm_square_docstring,
           boost::python::args("self"))
      .def(boost::python::self != boost::python::self)
      .def(boost::python::self == boost::python::self)
      .def(boost::python::self * boost::python::self)
      .def(boost::python::self + boost::python::self)
      .def(boost::python::self - boost::python::self)
      .def(-boost::python::self)
      .def(boost::python::self / float())
      .def(boost::python::self * float())
      .def(boost::python::self += boost::python::self)
      .def(boost::python::self -= boost::python::self)
      .def(boost::python::self *= float())
      .def(boost::python::self /= float())
      .add_property("x",
                    &PointWrapper::get_x,
                    &PointWrapper::set_x,
                    get_x_docstring)
      .add_property("y",
                    &PointWrapper::get_y,
                    &PointWrapper::set_y,
                    get_y_docstring)
      .add_property("z",
                    &PointWrapper::get_z,
                    &PointWrapper::set_z,
                    get_z_docstring);
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE

#if defined(__clang__) && __clang_major__ > 6
#  pragma GCC diagnostic pop
#endif
