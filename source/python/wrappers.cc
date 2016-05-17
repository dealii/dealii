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
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

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

void generate_cube(dealii::Triangulation<2> &triangulation)
{
  dealii::GridGenerator::hyper_cube(triangulation);
}

BOOST_PYTHON_MODULE(PyDealII)
{
  boost::python::scope().attr("__doc__") = pydealii_docstring;

  boost::python::docstring_options doc_options;
  doc_options.enable_user_defined();
  doc_options.enable_py_signatures();
  doc_options.disable_cpp_signatures();


  boost::python::class_<dealii::Triangulation<2>> ("Triangulation")
                                               .def("n_active_cells", &n_active_cells)
                                               .def("generate_cube", &generate_cube)
                                               .def("refine_global", &dealii::Triangulation<2>::refine_global);
}
