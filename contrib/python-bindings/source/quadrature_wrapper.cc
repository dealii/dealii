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

#include <deal.II/base/quadrature_lib.h>

#include <quadrature_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  namespace internal
  {
    template <int dim>
    boost::python::list
    get_points(const void *quadrature_ptr)
    {
      const auto quadrature =
        static_cast<const Quadrature<dim> *>(quadrature_ptr);

      const auto points = quadrature->get_points();

      boost::python::list points_list;
      for (const auto &p : points)
        {
          boost::python::list python_p;
          for (int d = 0; d < dim; ++d)
            python_p.append(p[d]);

          points_list.append(python_p);
        }

      return points_list;
    }

    template <int dim>
    boost::python::list
    get_weights(const void *quadrature_ptr)
    {
      const auto quadrature =
        static_cast<const Quadrature<dim> *>(quadrature_ptr);

      const auto weights = quadrature->get_weights();

      boost::python::list weights_list;
      for (const auto &w : weights)
        weights_list.append(w);

      return weights_list;
    }
  } // namespace internal



  QuadratureWrapper::QuadratureWrapper(const int dim)
    : dim(dim)
    , quadrature_ptr(nullptr)
  {
    AssertThrow(dim == 2 || dim == 3, ExcMessage("Unsupported dimension."));
  }



  QuadratureWrapper::QuadratureWrapper(const QuadratureWrapper &other)
  {
    dim = other.dim;

    AssertThrow(other.quadrature_ptr != nullptr,
                ExcMessage("Underlying quadrature does not exist."));

    if (dim == 2)
      {
        const auto quadrature =
          static_cast<const Quadrature<2> *>(other.quadrature_ptr);
        quadrature_ptr = new Quadrature<2>(*quadrature);
      }
    else if (dim == 3)
      {
        const auto quadrature =
          static_cast<const Quadrature<3> *>(other.quadrature_ptr);
        quadrature_ptr = new Quadrature<3>(*quadrature);
      }
    else
      AssertThrow(false, ExcMessage("Given dimension is not implemented."));
  }



  QuadratureWrapper::~QuadratureWrapper()
  {
    if (dim != -1)
      {
        if (dim == 2)
          {
            // We cannot call delete on a void pointer so cast the void pointer
            // back first.
            Quadrature<2> *tmp = static_cast<Quadrature<2> *>(quadrature_ptr);
            delete tmp;
          }
        else if (dim == 3)
          {
            Quadrature<3> *tmp = static_cast<Quadrature<3> *>(quadrature_ptr);
            delete tmp;
          }

        dim            = -1;
        quadrature_ptr = nullptr;
      }
  }



  void
  QuadratureWrapper::create_gauss(const unsigned int n)
  {
    if (dim == 2)
      {
        quadrature_ptr = new QGauss<2>(n);
      }
    else if (dim == 3)
      {
        quadrature_ptr = new QGauss<3>(n);
      }
    else
      AssertThrow(false, ExcMessage("Given dimension is not implemented."));
  }



  void
  QuadratureWrapper::create_gauss_lobatto(const unsigned int n)
  {
    if (dim == 2)
      {
        quadrature_ptr = new QGaussLobatto<2>(n);
      }
    else if (dim == 3)
      {
        quadrature_ptr = new QGaussLobatto<3>(n);
      }
    else
      AssertThrow(false, ExcMessage("Given dimension is not implemented."));
  }



  boost::python::list
  QuadratureWrapper::get_points() const
  {
    if (dim == 2)
      return internal::get_points<2>(quadrature_ptr);
    else if (dim == 3)
      return internal::get_points<3>(quadrature_ptr);
    else
      AssertThrow(false, ExcMessage("Given dimension is not implemented."));
  }



  boost::python::list
  QuadratureWrapper::get_weights() const
  {
    if (dim == 2)
      return internal::get_weights<2>(quadrature_ptr);
    else if (dim == 3)
      return internal::get_weights<3>(quadrature_ptr);
    else
      AssertThrow(false, ExcMessage("Given dimension is not implemented."));
  }



  void *
  QuadratureWrapper::get_quadrature() const
  {
    return quadrature_ptr;
  }



  int
  QuadratureWrapper::get_dim() const
  {
    return dim;
  }

} // namespace python

DEAL_II_NAMESPACE_CLOSE
