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

DEAL_II_NAMESPACE_OPEN

namespace python
{

// Macro to enable default arguments
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_cube_overloads,
                                         generate_hyper_cube, 0, 3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_subdivided_hyper_cube_overloads,
                                         generate_subdivided_hyper_cube, 1, 3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_rectangle_overloads,
                                         generate_hyper_rectangle, 2, 3)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_subdivided_hyper_rectangle_overloads,
                                         generate_subdivided_hyper_rectangle, 3, 4)
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generate_hyper_ball_overloads,
                                         generate_hyper_ball, 1, 2)



  const char n_active_cells_docstring [] =
    "Return the number of active cells                                      \n"
    ;



  const char generate_hyper_cube_docstring [] =
    "Generate a hyper_cube (square in 2D and cube in 3D)                    \n"
    "with exactly one cell                                                  \n"
    ;



  const char generate_simplex_docstring [] =
    "Generate a simplex with (dim+1) vertices and mesh cells                \n"
    ;



  const char generate_subdivided_hyper_cube_docstring [] =
    "Same as hyper_cube but not only one cell is created but each           \n"
    "coordinate direction is subdivided in repetitions cells                \n"
    ;



  const char generate_hyper_rectangle_docstring [] =
    "Generate a coordinate-parallel brick from the two diagonally opposite  \n"
    "corners points p1 and p2                                               \n"
    ;



  const char generate_subdivided_hyper_rectangle_docstring [] =
    "Generate a coordinate-parallel brick from the two diagonally opposite  \n"
    "corners point p1 and p2. In direction i, repetitions[i] cells are      \n"
    "created                                                                \n"
    ;



  const char generate_hyper_ball_docstring [] =
    "Generate a hyperball, i.e., a circle or a ball around center with     \n"
    "a given radius                                                        \n"
    ;



  const char shift_docstring [] =
    "Shift every vertex of the Triangulation by the gien shift vector      \n"
    ;



  const char merge_docstring [] =
    "Given two triangulations, create the triangulation that contains      \n"
    "the cells of both triangulations                                      \n"
    ;



  const char refine_global_docstring [] =
    "Refine all the cells times time                                       \n"
    ;



  const char execute_coarsening_and_refinement_docstring [] =
    "Execute both refinement and coarsening of the Triangulation           \n"
    ;


  const char active_cells_docstring [] =
    "Return the list of active cell accessors of the Triangulation         \n"
    ;



  const char write_docstring [] =
    "Write the mesh to the output file according to the given data format. \n"
    "The possible formats are:                                             \n"
    "  - none                                                              \n"
    "  - dx                                                                \n"
    "  - gnuplot                                                           \n"
    "  - eps                                                               \n"
    "  - ucd                                                               \n"
    "  - xfig                                                              \n"
    "  - msh                                                               \n"
    "  - svg                                                               \n"
    "  - mathgl                                                            \n"
    "  - vtk                                                               \n"
    "  - vtu                                                               \n"
    ;



  const char save_docstring [] =
    "Write the Triangulation to a file                                     \n"
    ;




  const char load_docstring [] =
    "Load the Triangulation from a file                                    \n"
    ;



  void export_triangulation()
  {
    boost::python::class_<TriangulationWrapper>("Triangulation",
                                                boost::python::init<const std::string &>())
    .def(boost::python::init<const std::string &, const std::string &>())
    .def("n_active_cells",
         &TriangulationWrapper::n_active_cells,
         n_active_cells_docstring,
         boost::python::args("self"))
    .def("generate_hyper_cube",
         &TriangulationWrapper::generate_hyper_cube,
         generate_hyper_cube_overloads(
           boost::python::args("self", "left", "right", "colorize"),
           generate_hyper_cube_docstring))
    .def("generate_simplex",
         &TriangulationWrapper::generate_simplex,
         generate_simplex_docstring,
         boost::python::args("self", "vertices"))
    .def("generate_subdivided_hyper_cube",
         &TriangulationWrapper::generate_subdivided_hyper_cube,
         generate_subdivided_hyper_cube_overloads(
           boost::python::args("self", "repetitions", "left", "right"),
           generate_subdivided_hyper_cube_docstring))
    .def("generate_hyper_rectangle",
         &TriangulationWrapper::generate_hyper_rectangle,
         generate_hyper_rectangle_overloads(
           boost::python::args("self", "p1", "p2", "colorize"),
           generate_hyper_rectangle_docstring))
    .def("generate_subdivided_hyper_rectangle",
         &TriangulationWrapper::generate_subdivided_hyper_rectangle,
         generate_subdivided_hyper_rectangle_overloads(
           boost::python::args("self", "repetitions",
                               "p1", "p2", "colorize"),
           generate_subdivided_hyper_rectangle_docstring))
    .def("generate_hyper_ball",
         &TriangulationWrapper::generate_hyper_ball,
         generate_hyper_ball_overloads(
           boost::python::args("self", "center", "radius"),
           generate_hyper_ball_docstring))
    .def("shift",
         &TriangulationWrapper::shift,
         shift_docstring,
         boost::python::args("self", "shift"))
    .def("merge_triangulations",
         &TriangulationWrapper::merge_triangulations,
         merge_docstring,
         boost::python::args("self", "triangulation_1", "triangulation_2"))
    .def("refine_global",
         &TriangulationWrapper::refine_global,
         refine_global_docstring,
         boost::python::args("self", "times"))
    .def("execute_coarsening_and_refinement",
         &TriangulationWrapper::execute_coarsening_and_refinement,
         execute_coarsening_and_refinement_docstring,
         boost::python::args("self"))
    .def("active_cells",
         &TriangulationWrapper::active_cells,
         active_cells_docstring,
         boost::python::args("self"))
    .def("write",
         &TriangulationWrapper::write,
         write_docstring,
         boost::python::args("self", "filename", "format"))
    .def("save",
         &TriangulationWrapper::save,
         save_docstring,
         boost::python::args("self", "filename"))
    .def("load",
         &TriangulationWrapper::load,
         load_docstring,
         boost::python::args("self", "filename"));
  }
}

DEAL_II_NAMESPACE_CLOSE
