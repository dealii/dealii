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


#ifndef dealii_matrix_free_util_h
#define dealii_matrix_free_util_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/simplex/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int dim>
    inline Quadrature<dim - 1>
    get_face_quadrature(const Quadrature<dim> &quad)
    {
      if (dim == 2 || dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == Simplex::QGauss<dim>(i))
            return Simplex::QGauss<dim - 1>(i);

      AssertThrow(false, ExcNotImplemented());

      return Quadrature<dim - 1>();
    }

    template <int dim>
    inline std::vector<Quadrature<dim - 1>>
    get_unique_face_quadratures(const Quadrature<dim> &quad)
    {
      if (dim == 2 || dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == Simplex::QGauss<dim>(i))
            return {{Simplex::QGauss<dim - 1>(i)}};

      if (dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == Simplex::QGaussWedge<dim>(i))
            return {{QGauss<dim - 1>(i), Simplex::QGauss<dim - 1>(i)}};

      if (dim == 3)
        for (unsigned int i = 1; i <= 2; ++i)
          if (quad == Simplex::QGaussPyramid<dim>(i))
            return {{QGauss<dim - 1>(i), Simplex::QGauss<dim - 1>(i)}};

      AssertThrow(false, ExcNotImplemented());

      return {{QGauss<dim - 1>(1)}};
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
