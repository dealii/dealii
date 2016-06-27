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
#include <deal.II/python/triangulation_wrapper.h>
#include <deal.II/python/cell_accessor_wrapper.h>

namespace PyDealII
{

// Macro to enable default arguments
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_cube_overloads, generate_hyper_cube, 0, 3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_subdivided_hyper_cube_overloads,
                                         generate_subdivided_hyper_cube, 1, 3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_rectangle_overloads,
                                         generate_hyper_rectangle, 2, 3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_subdivided_hyper_rectangle_overloads,
                                         generate_subdivided_hyper_rectangle, 3, 4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_ball_overloads, generate_hyper_ball, 1, 2)


  void export_triangulation()
  {
    boost::python::class_<TriangulationWrapper>("Triangulation", boost::python::init<const std::string &>())
    .def(boost::python::init<const std::string &, const std::string &>())
    .def("__iter__", boost::python::iterator<TriangulationWrapper>())
    .def("n_active_cells",
         &TriangulationWrapper::n_active_cells,
         "Return the number of active cells",
         boost::python::args("self"))
    .def("generate_hyper_cube",
         &TriangulationWrapper::generate_hyper_cube,
         generate_hyper_cube_overloads(
           boost::python::args("self", "left", "right", "colorize"),
           "Generate a hyper_cube (square in 2D and cube in 3D) with exactly one cell."))
    .def("generate_simplex",
         &TriangulationWrapper::generate_simplex,
         "Generate a simplex with (dim+1) vertices and mesh cells.",
         boost::python::args("self", "vertices"))
    .def("generate_subdivided_hyper_cube",
         &TriangulationWrapper::generate_subdivided_hyper_cube,
         generate_subdivided_hyper_cube_overloads(
           boost::python::args("self", "repetitions", "left", "right"),
           "Same as hyper_cube but not only one cells is created but each coordinate direction is subdivdided in repetitions cells."))
    .def("generate_hyper_rectangle",
         &TriangulationWrapper::generate_hyper_rectangle,
         generate_hyper_rectangle_overloads(
           boost::python::args("self", "p1", "p2", "colorize"),
           "Generate a coordinate-parallel brick from the two diagonally opposite corners points p1 and p2."))
    .def("generate_subdivided_hyper_rectangle",
         &TriangulationWrapper::generate_subdivided_hyper_rectangle,
         generate_subdivided_hyper_rectangle_overloads(
           boost::python::args("self", "repetitions",
                               "p1", "p2", "colorize"),
           "Generate a coordinate-parallel brick from the two diagonally opposite corners point @p1 and @p2. In direction i, repetitions[i] cells are created."))
    .def("generate_hyper_ball",
         &TriangulationWrapper::generate_hyper_ball,
         generate_hyper_ball_overloads(
           boost::python::args("self", "center", "radius"),
           "Generate a hyperball, i.e. a circle or a ball around @p center with given @p radius."))
    .def("shift",
         &TriangulationWrapper::shift,
         "Shift every vertex of the Triangulation by the given shift vector.",
         boost::python::args("self", "shift"))
    .def("merge_triangulations",
         &TriangulationWrapper::merge_triangulations,
         "Given two triangulations, create the triangulation that contains the cells of both triangulations.",
         boost::python::args("self", "triangulation_1", "triangulation_2"))
    .def("refine_global",
         &TriangulationWrapper::refine_global,
         "Refine all the cells times times.",
         boost::python::args("self", "times"))
    .def("execute_coarsening_and_refinement",
         &TriangulationWrapper::execute_coarsening_and_refinement,
         "Execute both refinement and coarsening of the Triangulation.",
         boost::python::args("self"))
    .def("write",
         &TriangulationWrapper::write,
         "Write grid to the output file according to the given data format. The possible formats are: none, dx, gnuplot, eps, ucd, xfig, msh, svg, mathgl, vtk, and vtu.",
         boost::python::args("self", "filename", "format"))
    .def("save",
         &TriangulationWrapper::save,
         "Write the Triangulation to a file.",
         boost::python::args("self", "filename"))
    .def("load",
         &TriangulationWrapper::load,
         "Load the Triangulation from a file.",
         boost::python::args("self", "filename"));
  }
}
