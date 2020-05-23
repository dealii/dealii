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

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  void
  export_tria_accessor();
  void
  export_cell_accessor();
  void
  export_point();
  void
  export_triangulation();
  void
  export_mapping();
  void
  export_manifold();
  void
  export_quadrature();
} // namespace python

DEAL_II_NAMESPACE_CLOSE

char const *pydealii_docstring =
  "                                                             \n"
  "PyDealII                                                     \n"
  "========                                                     \n"
  "This module contains the python bindings to deal.II.         \n"
  "The Debug module uses deal.II compiled in Debug mode while   \n"
  "the Release module uses deal.II compiled in Release mode.    \n";

#ifdef DEBUG

BOOST_PYTHON_MODULE(Debug)
{
  boost::python::scope().attr("__doc__") = pydealii_docstring;

  boost::python::docstring_options doc_options;
  doc_options.enable_user_defined();
  doc_options.enable_py_signatures();
  doc_options.disable_cpp_signatures();

  // Switch off call to std::abort when an exception is created using Assert.
  // If the code aborts, the kernel of a Jupyter Notebook is killed and no
  // message is printed.
  dealii::deal_II_exceptions::disable_abort_on_exception();

  dealii::python::export_tria_accessor();
  dealii::python::export_cell_accessor();
  dealii::python::export_point();
  dealii::python::export_triangulation();
  dealii::python::export_mapping();
  dealii::python::export_manifold();
  dealii::python::export_quadrature();
}

#else

BOOST_PYTHON_MODULE(Release)
{
  boost::python::scope().attr("__doc__") = pydealii_docstring;

  boost::python::docstring_options doc_options;
  doc_options.enable_user_defined();
  doc_options.enable_py_signatures();
  doc_options.disable_cpp_signatures();

  dealii::python::export_tria_accessor();
  dealii::python::export_cell_accessor();
  dealii::python::export_point();
  dealii::python::export_triangulation();
  dealii::python::export_mapping();
  dealii::python::export_manifold();
  dealii::python::export_quadrature();
}

#endif
