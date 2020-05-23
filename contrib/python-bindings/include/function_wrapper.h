// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_function_wrapper_h
#define dealii_function_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>

#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  template <int dim>
  class FunctionWrapper : public Function<dim, double>
  {
  public:
    FunctionWrapper(boost::python::object &python_function,
                    unsigned               n_components)
      : Function<dim, double>(n_components)
      , python_function(python_function)
    {}

    FunctionWrapper(const FunctionWrapper &other)
    {
      python_function = other.python_function;
    }

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      boost::python::list p_list_in;

      for (int d = 0; d < dim; ++d)
        p_list_in.append(p[d]);

      boost::python::list p_list_out =
        boost::python::extract<boost::python::list>(python_function(p_list_in));

      for (size_t i = 0; i < this->n_components; ++i)
        values[i] = boost::python::extract<double>(p_list_out[i]);
    }

  private:
    /**
     * A callback to a python function
     */
    boost::python::object python_function;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
