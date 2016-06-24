// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <boost/python.hpp>
#include <deal.II/python/cell_accessor_wrapper.h>
#include <deal.II/python/triangulation_wrapper.h>

namespace PyDealII
{
  void export_cell_accessor()
  {
    boost::python::class_<CellAccessorWrapper>("CellAccessor",
                                               boost::python::init<TriangulationWrapper &, const int, const int> ())
    .add_property("refine_flag", &CellAccessorWrapper::get_refine_flag,
                  &CellAccessorWrapper::set_refine_flag,
                  "Get the refine_flag of the cell. The possibilities are in 2D: isotropic, no_refinement, cut_x, cut_y, and cut_xy and in 3D: isotropic, no_refinement, cut_x, cut_y, cut_z, cut_xy, cut_xz, cut_yz, and cut_xyz.")
    .add_property("coarsen_flag", &CellAccessorWrapper::get_coarsen_flag,
                  &CellAccessorWrapper::set_coarsen_flag,
                  "Get the coarsen_flag of the cell.")
    .add_property("material_id", &CellAccessorWrapper::get_material_id,
                  &CellAccessorWrapper::set_material_id,
                  "Get the material_id of the cell.")
    .add_property("manifold_id", &CellAccessorWrapper::get_manifold_id,
                  &CellAccessorWrapper::set_manifold_id,
                  "Get the manifold_id of the cell.")
    .def("barycenter", &CellAccessorWrapper::get_barycenter,
         "Return the barycenter of the current cell.",
         boost::python::args("self"))
    .def("set_vertex", &CellAccessorWrapper::set_vertex,
         "Set the ith vertex of the cell to point_wrapper",
         boost::python::args("self", "i", "point_wrapper"))
    .def("get_vertex", &CellAccessorWrapper::get_vertex,
         "Get the ith vertex of the cell",
         boost::python::args("self", "i"));
  }
}

