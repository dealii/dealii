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

#include <cell_accessor_wrapper.h>
#include <triangulation_wrapper.h>
#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  const char refine_flag_docstring [] =
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
    "  - cut_xyz                                                        \n"
    ;



  const char coarsen_flag_docstring [] =
    "Get/Set the coarsen_flag of the cell                               \n"
    ;



  const char material_id_docstring [] =
    "Get/Set the material_id of the cell                                \n"
    ;



  const char manifold_id_docstring [] =
    "Get/Set the manifold_id of the cell                                \n"
    ;



  const char barycenter_docstring [] =
    "Return the barycenter of the current cell                          \n"
    ;



  const char set_vertex_docstring [] =
    " Set the ith vertex of the cell to point_wrapper                   \n"
    ;



  const char get_vertex_docstring [] =
    " Get the ith vertex of the cell                                    \n"
    ;



  void export_cell_accessor()
  {
    boost::python::class_<CellAccessorWrapper>("CellAccessor",
                                               boost::python::init<TriangulationWrapper &, const int, const int> ())
    .add_property("refine_flag", &CellAccessorWrapper::get_refine_flag,
                  &CellAccessorWrapper::set_refine_flag,
                  refine_flag_docstring)
    .add_property("coarsen_flag", &CellAccessorWrapper::get_coarsen_flag,
                  &CellAccessorWrapper::set_coarsen_flag,
                  coarsen_flag_docstring)
    .add_property("material_id", &CellAccessorWrapper::get_material_id,
                  &CellAccessorWrapper::set_material_id,
                  material_id_docstring)
    .add_property("manifold_id", &CellAccessorWrapper::get_manifold_id,
                  &CellAccessorWrapper::set_manifold_id,
                  manifold_id_docstring)
    .def("barycenter", &CellAccessorWrapper::get_barycenter,
         barycenter_docstring,
         boost::python::args("self"))
    .def("set_vertex", &CellAccessorWrapper::set_vertex,
         set_vertex_docstring,
         boost::python::args("self", "i", "point_wrapper"))
    .def("get_vertex", &CellAccessorWrapper::get_vertex,
         get_vertex_docstring,
         boost::python::args("self", "i"));
  }
}

DEAL_II_NAMESPACE_CLOSE
