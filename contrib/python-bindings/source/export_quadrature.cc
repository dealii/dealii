// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <boost/python.hpp>

#include <quadrature_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  const char create_gauss_docstring[] =
    "Create Gauss quadrature with n points                     \n"
    "in each space direction.                                  \n";


  const char create_gauss_lobatto_docstring[] =
    "Create Gauss-Lobatto quadrature with n points             \n"
    "in each space direction.                                  \n";


  const char get_points_docstring[] =
    "Return the list of quadrature points.                     \n";


  const char get_weights_docstring[] =
    "Return the list of quadrature weights.                    \n";


  void
  export_quadrature()
  {
    boost::python::class_<QuadratureWrapper>(
      "Quadrature", boost::python::init<const int>(boost::python::args("dim")))
      .def("create_gauss",
           &QuadratureWrapper::create_gauss,
           create_gauss_docstring,
           boost::python::args("self", "n"))
      .def("create_gauss_lobatto",
           &QuadratureWrapper::create_gauss_lobatto,
           create_gauss_lobatto_docstring,
           boost::python::args("self", "n"))
      .def("points",
           &QuadratureWrapper::get_points,
           get_points_docstring,
           boost::python::args("self"))
      .def("weights",
           &QuadratureWrapper::get_weights,
           get_weights_docstring,
           boost::python::args("self"));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
