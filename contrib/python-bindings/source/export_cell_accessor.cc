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

#include <cell_accessor_wrapper.h>
#include <triangulation_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_center_overloads, get_center, 0, 2)

  const char refine_flag_docstring[] =
    "Get/Set the refine_flag of the cell. In 2D, the possibilities are: \n"
    "  - isotropic                                                      \n"
    "  - no_refinement                                                  \n"
    "  - cut_x                                                          \n"
    "  - cut_y                                                          \n"
    "  - cut_xy                                                         \n"
    "In 3D, the possibilities are:                                      \n"
    "  - isotropic                                                      \n"
    "  - no_refinement                                                  \n"
    "  - cut_x                                                          \n"
    "  - cut_y                                                          \n"
    "  - cut_z                                                          \n"
    "  - cut_xy                                                         \n"
    "  - cut_xz                                                         \n"
    "  - cut_yz                                                         \n"
    "  - cut_xyz                                                        \n";



  const char coarsen_flag_docstring[] =
    "Get/Set the coarsen_flag of the cell.                              \n";



  const char material_id_docstring[] =
    "Get/Set the material_id of the cell.                               \n";



  const char manifold_id_docstring[] =
    "Get/Set the manifold_id of the cell.                               \n";



  const char set_all_manifold_ids_docstring[] =
    "Do as set_manifold_id() but also set the manifold indicators of    \n"
    "the objects that bound the current object.                         \n";



  const char barycenter_docstring[] =
    "Return the barycenter of the current cell.                         \n";



  const char center_docstring[] =
    "Center of the object. This method can take into account manifold   \n"
    "attached to a cell.                                                \n";



  const char set_vertex_docstring[] =
    " Set the ith vertex of the cell to point_wrapper.                  \n";



  const char get_vertex_docstring[] =
    " Get the ith vertex of the cell.                                   \n";



  const char at_boundary_docstring[] =
    " Return whether the cell is at the boundary.                       \n";



  const char has_boundary_lines_docstring[] =
    " This is a slight variation to the at_boundary function:           \n"
    " for two dimensions, it is equivalent, for three                   \n"
    " dimensions it returns whether at least one of the 12              \n"
    " lines of the hexahedron is at a boundary.                         \n";



  const char neighbor_docstring[] =
    " Return the ith neighbor of a cell. If the neighbor does not exist,\n"
    " i.e., if the ith face of the current object is at the boundary,   \n"
    " then an exception is thrown.                                      \n";



  const char faces_docstring[] =
    " Return list of cell's faces.                                      \n";



  const char measure_docstring[] =
    " Compute the dim-dimensional measure of the object.                 \n";



  const char vertex_index_docstring[] =
    " Return the global index of i-th vertex of a cell.                  \n";



  const char neighbor_of_neighbor_docstring[] =
    " Return the how-many'th neighbor this cell is of                    \n"
    " cell.neighbor(neighbor), i.e. return the face_no such that         \n"
    " cell.neighbor(neighbor).neighbor(face_no)==cell.                   \n";



  const char neighbor_is_coarser_docstring[] =
    " Return whether the neighbor is coarser then the present cell.      \n";



  const char index_docstring[] =
    " Return the index of the element presently pointed to               \n"
    " on the present.                                                    \n";



  const char level_docstring[] =
    " Return the level within the mesh hierarchy at which this cell      \n"
    " is located.                                                        \n";



  const char active_docstring[] =
    " Test whether the cell has children (this is the criterion          \n"
    " for activity of a cell).                                           \n";


  void
  export_cell_accessor()
  {
    boost::python::class_<CellAccessorWrapper>(
      "CellAccessor",
      boost::python::init<TriangulationWrapper &, const int, const int>())
      .add_property("refine_flag",
                    &CellAccessorWrapper::get_refine_flag,
                    &CellAccessorWrapper::set_refine_flag,
                    refine_flag_docstring)
      .add_property("coarsen_flag",
                    &CellAccessorWrapper::get_coarsen_flag,
                    &CellAccessorWrapper::set_coarsen_flag,
                    coarsen_flag_docstring)
      .add_property("material_id",
                    &CellAccessorWrapper::get_material_id,
                    &CellAccessorWrapper::set_material_id,
                    material_id_docstring)
      .add_property("manifold_id",
                    &CellAccessorWrapper::get_manifold_id,
                    &CellAccessorWrapper::set_manifold_id,
                    manifold_id_docstring)
      .def("set_all_manifold_ids",
           &CellAccessorWrapper::set_all_manifold_ids,
           set_all_manifold_ids_docstring,
           boost::python::args("self", "number"))
      .def("barycenter",
           &CellAccessorWrapper::get_barycenter,
           barycenter_docstring,
           boost::python::args("self"))
      .def("center",
           &CellAccessorWrapper::get_center,
           get_center_overloads(
             boost::python::args("self",
                                 "respect_manifold",
                                 "interpolate_from_surrounding"),
             center_docstring))
      .def("set_vertex",
           &CellAccessorWrapper::set_vertex,
           set_vertex_docstring,
           boost::python::args("self", "i", "point_wrapper"))
      .def("get_vertex",
           &CellAccessorWrapper::get_vertex,
           get_vertex_docstring,
           boost::python::args("self", "i"))
      .def("at_boundary",
           &CellAccessorWrapper::at_boundary,
           at_boundary_docstring,
           boost::python::args("self"))
      .def("has_boundary_lines",
           &CellAccessorWrapper::has_boundary_lines,
           has_boundary_lines_docstring,
           boost::python::args("self"))
      .def("neighbor",
           &CellAccessorWrapper::neighbor,
           neighbor_docstring,
           boost::python::args("self", "i"))
      .def("faces",
           &CellAccessorWrapper::faces,
           faces_docstring,
           boost::python::args("self"))
      .def("measure",
           &CellAccessorWrapper::measure,
           measure_docstring,
           boost::python::args("self"))
      .def("active",
           &CellAccessorWrapper::active,
           active_docstring,
           boost::python::args("self"))
      .def("level",
           &CellAccessorWrapper::level,
           level_docstring,
           boost::python::args("self"))
      .def("index",
           &CellAccessorWrapper::index,
           index_docstring,
           boost::python::args("self"))
      .def("neighbor_is_coarser",
           &CellAccessorWrapper::neighbor_is_coarser,
           neighbor_is_coarser_docstring,
           boost::python::args("self", "neighbor"))
      .def("neighbor_of_neighbor",
           &CellAccessorWrapper::neighbor_of_neighbor,
           neighbor_of_neighbor_docstring,
           boost::python::args("self", "neighbor"))
      .def("vertex_index",
           &CellAccessorWrapper::vertex_index,
           vertex_index_docstring,
           boost::python::args("self", "vertex"));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
