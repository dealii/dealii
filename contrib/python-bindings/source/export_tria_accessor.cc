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

#include <tria_accessor_wrapper.h>
#include <triangulation_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_center_overloads, get_center, 0, 2)

  const char manifold_id_docstring[] =
    "Get/Set the manifold_id of the face                                \n";

  const char boundary_id_docstring[] =
    "Get/Set the boundary_id of the face                                \n";

  const char set_all_boundary_ids_docstring[] =
    "Do as set_boundary_id() but also set the boundary indicators       \n"
    "of the objects that bound the current object.                      \n";

  const char barycenter_docstring[] =
    "Return the barycenter of the current face                          \n";

  const char center_docstring[] =
    "Return the center of the current face taking into account manifold.\n";

  const char set_vertex_docstring[] =
    " Set the ith vertex of the face to point_wrapper                   \n";

  const char get_vertex_docstring[] =
    " Get the ith vertex of the face                                    \n";

  const char at_boundary_docstring[] =
    " Return whether the face is at the boundary                        \n";

  const char measure_docstring[] =
    " Compute the dim-dimensional measure of the object.                 \n";

  void
  export_tria_accessor()
  {
    boost::python::class_<TriaAccessorWrapper>(
      "TriaAccessor",
      boost::python::init<void *, const int, const int, const int>())
      .add_property("boundary_id",
                    &TriaAccessorWrapper::get_boundary_id,
                    &TriaAccessorWrapper::set_boundary_id,
                    boundary_id_docstring)
      .add_property("manifold_id",
                    &TriaAccessorWrapper::get_manifold_id,
                    &TriaAccessorWrapper::set_manifold_id,
                    manifold_id_docstring)
      .def("barycenter",
           &TriaAccessorWrapper::get_barycenter,
           barycenter_docstring,
           boost::python::args("self"))
      .def("center",
           &TriaAccessorWrapper::get_center,
           get_center_overloads(
             boost::python::args("self",
                                 "respect_manifold",
                                 "interpolate_from_surrounding"),
             center_docstring))
      .def("set_vertex",
           &TriaAccessorWrapper::set_vertex,
           set_vertex_docstring,
           boost::python::args("self", "i", "point_wrapper"))
      .def("get_vertex",
           &TriaAccessorWrapper::get_vertex,
           get_vertex_docstring,
           boost::python::args("self", "i"))
      .def("at_boundary",
           &TriaAccessorWrapper::at_boundary,
           at_boundary_docstring,
           boost::python::args("self"))
      .def("set_all_boundary_ids",
           &TriaAccessorWrapper::set_all_boundary_ids,
           set_all_boundary_ids_docstring,
           boost::python::args("self", "boundary_id"))
      .def("measure",
           &TriaAccessorWrapper::measure,
           measure_docstring,
           boost::python::args("self"));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
