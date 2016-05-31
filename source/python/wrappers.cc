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

char const *pydealii_docstring =
  "                                                             \n"
  "PyDealII                                                     \n"
  "========                                                     \n"
  "Some interesting doc.                                        \n"
  ;

unsigned int n_active_cells(const dealii::Triangulation<2> &triangulation)
{
  return triangulation.n_active_cells();
}

void generate_hyper_cube(dealii::Triangulation<2> &triangulation, const double left=0.,
                         const double right=0., const bool colorize=false)
{
  dealii::GridGenerator::hyper_cube(triangulation);
}

void save(const dealii::Triangulation<2> &triangulation,
          const std::string filename)
{
  std::ofstream ofs(filename);
  boost::archive::text_oarchive oa(ofs);
  oa << triangulation;
}

void load(dealii::Triangulation<2> &triangulation,
          const std::string filename)
{
  std::ifstream ifs(filename);
  boost::archive::text_iarchive ia(ifs);
  ia >> triangulation;
}

// Macro to enable default arguments
BOOST_PYTHON_FUNCTION_OVERLOADS(generate_hyper_cube_overloads, generate_hyper_cube, 1, 4)

BOOST_PYTHON_MODULE(PyDealII)
{
  boost::python::scope().attr("__doc__") = pydealii_docstring;

  boost::python::docstring_options doc_options;
  doc_options.enable_user_defined();
  doc_options.enable_py_signatures();
  doc_options.disable_cpp_signatures();


  boost::python::class_<dealii::Triangulation<2>> ("Triangulation")
                                               .def("n_active_cells", n_active_cells,
                                                    "Return the number of active cells",
                                                    boost::python::args("self"))
                                               .def("generate_hyper_cube", generate_hyper_cube,
                                                    generate_hyper_cube_overloads(
                                                      boost::python::args("self", "left", "right", "colorize"),
                                                      "Generate a hyper_cube."))
                                               .def("refine_global", &dealii::Triangulation<2>::refine_global,
                                                    "Refine the mesh uniformly",
                                                    boost::python::args("self", "times"))
                                               .def("save", save, "Serialize and save the triangulation",
                                                    boost::python::args("self", "filename"))
                                               .def("load", load, "Load and deserialize a triangulation",
                                                    boost::python::args("self", "filename"));
}
