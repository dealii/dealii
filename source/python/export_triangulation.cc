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

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <fstream>

namespace PyDealII
{

  template <int dim>
  unsigned int n_active_cells(const dealii::Triangulation<dim> &triangulation)
  {
    return triangulation.n_active_cells();
  }



  template <int dim>
  void generate_hyper_cube(dealii::Triangulation<dim> &triangulation,
                           const double                left = 0.,
                           const double                right = 0.,
                           const bool                  colorize = false)
  {
    dealii::GridGenerator::hyper_cube(triangulation, left, right, colorize);
  }



  void generate_simplex(dealii::Triangulation<2> &triangulation,
                        dealii::Point<2>          vertex_1,
                        dealii::Point<2>          vertex_2,
                        dealii::Point<2>          vertex_3)
  {
    std::vector<dealii::Point<2>> vertices(3);
    vertices[0] = vertex_1;
    vertices[1] = vertex_2;
    vertices[2] = vertex_3;

    dealii::GridGenerator::simplex(triangulation, vertices);
  }



  template <int dim>
  void generate_subdivided_hyper_cube(dealii::Triangulation<dim> &triangulation,
                                      const unsigned int          repetitions,
                                      const double                left = 0.,
                                      const double                right = 1.)
  {
    dealii::GridGenerator::subdivided_hyper_cube(triangulation, repetitions,
                                                 left, right);
  }



  template <int dim>
  void generate_hyper_rectangle(dealii::Triangulation<dim> &triangulation,
                                const dealii::Point<dim>   &p1,
                                const dealii::Point<dim>   &p2,
                                const bool                  colorize = false)
  {
    dealii::GridGenerator::hyper_rectangle(triangulation, p1, p2, colorize);
  }



  void generate_subdivided_hyper_rectangle(dealii::Triangulation<2> &triangulation,
                                           const unsigned int        repetition_x,
                                           const unsigned int        repetition_y,
                                           const dealii::Point<2>   &p1,
                                           const dealii::Point<2>   &p2,
                                           const bool                colorize = false)
  {
    std::vector<unsigned int> repetitions(2);
    repetitions[0] = repetition_x;
    repetitions[1] = repetition_y;
    dealii::GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions,
                                                      p1, p2, colorize);
  }



  template <int dim>
  void generate_hyper_ball(dealii::Triangulation<dim> &triangulation,
                           const dealii::Point<dim>   &center = dealii::Point<dim>(),
                           const double                radius = 1.)
  {
    dealii::GridGenerator::hyper_ball(triangulation, center, radius);
  }



  template <int dim>
  void save(const dealii::Triangulation<dim> &triangulation,
            const std::string                 filename)
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << triangulation;
  }



  template <int dim>
  void load(dealii::Triangulation<dim> &triangulation,
            const std::string           filename)
  {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> triangulation;
  }



// Macro to enable default arguments
  BOOST_PYTHON_FUNCTION_OVERLOADS(generate_hyper_cube_overloads, generate_hyper_cube, 1, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(generate_subdivided_hyper_cube_overloads,
                                  generate_subdivided_hyper_cube, 2, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(generate_hyper_rectangle_overloads,
                                  generate_hyper_rectangle, 3, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(generate_subdivided_hyper_rectangle_overloads,
                                  generate_subdivided_hyper_rectangle, 5, 6)
  BOOST_PYTHON_FUNCTION_OVERLOADS(generate_hyper_ball_overloads, generate_hyper_ball, 1, 3)


  void export_triangulation()
  {
    boost::python::class_<dealii::Triangulation<2>> ("Triangulation")
                                                 .def("n_active_cells", n_active_cells<2>,
                                                      "Return the number of active cells",
                                                      boost::python::args("self"))
                                                 .def("generate_hyper_cube", generate_hyper_cube<2>,
                                                      generate_hyper_cube_overloads(
                                                        boost::python::args("self", "left", "right", "colorize"),
                                                        "Generate a hyper_cube."))
                                                 .def("generate_simplex", generate_simplex,
                                                      "Generate a simplex.",
                                                      boost::python::args("self", "vertex_1", "vertex_2", "vertex_3"))
                                                 .def("generate_subdivided_hyper_cube", generate_subdivided_hyper_cube<2>,
                                                      generate_subdivided_hyper_cube_overloads(
                                                        boost::python::args("self", "repetitions", "left", "right"),
                                                        "Generate a subdivided hyper_cube."))
                                                 .def("generate_hyper_rectangle", generate_hyper_rectangle<2>,
                                                      generate_hyper_rectangle_overloads(
                                                        boost::python::args("self", "p1", "p2", "colorize"),
                                                        "Generate a hyper_rectangle."))
                                                 .def("generate_subdivided_hyper_rectangle",
                                                      generate_subdivided_hyper_rectangle,
                                                      generate_subdivided_hyper_rectangle_overloads(
                                                        boost::python::args("self", "repetition_x",
                                                            "repetition_y","p1", "p2",
                                                            "colorize"),
                                                        "Generate a subdivided hyper_rectangle."))
                                                 .def("generate_hyper_ball",
                                                      generate_hyper_ball<2>,
                                                      generate_hyper_ball_overloads(
                                                        boost::python::args("self", "center", "radius"),
                                                        "Generate a hyper_ball."))
                                                 .def("refine_global", &dealii::Triangulation<2>::refine_global,
                                                      "Refine the mesh uniformly.",
                                                      boost::python::args("self", "times"))
                                                 .def("save", save<2>, "Serialize and save the triangulation.",
                                                      boost::python::args("self", "filename"))
                                                 .def("load", load<2>, "Load and deserialize a triangulation.",
                                                      boost::python::args("self", "filename"));
  }

}
