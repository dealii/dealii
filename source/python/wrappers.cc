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

#include <deal.II/base/config.h>
#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  void export_cell_accessor();
  void export_point();
  void export_triangulation();
}

DEAL_II_NAMESPACE_CLOSE

char const *pydealii_docstring =
  "                                                             \n"
  "PyDealII                                                     \n"
  "========                                                     \n"
  "Some interesting doc.                                        \n"
  ;

BOOST_PYTHON_MODULE(PyDealII)
{
  boost::python::scope().attr("__doc__") = pydealii_docstring;

  boost::python::docstring_options doc_options;
  doc_options.enable_user_defined();
  doc_options.enable_py_signatures();
  doc_options.disable_cpp_signatures();

  dealii::python::export_cell_accessor();
  dealii::python::export_point();
  dealii::python::export_triangulation();
}
