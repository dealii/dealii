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

namespace PyDealII
{
  void export_cell_accessor();
  void export_point();
  void export_triangulation();
}

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

  PyDealII::export_cell_accessor();
  PyDealII::export_point();
  PyDealII::export_triangulation();
}
