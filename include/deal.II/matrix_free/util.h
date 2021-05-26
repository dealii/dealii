// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/q_collection.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int dim>
    inline std::pair<dealii::ReferenceCell, dealii::hp::QCollection<dim - 1>>
    get_face_quadrature_collection(const Quadrature<dim> &quad,
                                   const bool             do_assert = true)
    {
      if (dim == 2 || dim == 3)
        {
          for (unsigned int i = 1; i <= 4; ++i)
            if (quad == QGaussSimplex<dim>(i))
              {
                QGaussSimplex<dim - 1> tri(i);

                if (dim == 2)
                  return {ReferenceCells::Triangle,
                          dealii::hp::QCollection<dim - 1>(tri, tri, tri)};
                else
                  return {ReferenceCells::Tetrahedron,
                          dealii::hp::QCollection<dim - 1>(tri, tri, tri, tri)};
              }

          for (unsigned int i = 1; i <= 5; ++i)
            if (quad == QWitherdenVincentSimplex<dim>(i))
              {
                QWitherdenVincentSimplex<dim - 1> tri(i);

                if (dim == 2)
                  return {ReferenceCells::Triangle,
                          dealii::hp::QCollection<dim - 1>(tri, tri, tri)};
                else
                  return {ReferenceCells::Tetrahedron,
                          dealii::hp::QCollection<dim - 1>(tri, tri, tri, tri)};
              }
        }

      if (dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == QGaussWedge<dim>(i))
            {
              QGauss<dim - 1>        quad(i);
              QGaussSimplex<dim - 1> tri(i);

              return {
                ReferenceCells::Wedge,
                dealii::hp::QCollection<dim - 1>(tri, tri, quad, quad, quad)};
            }

      if (dim == 3)
        for (unsigned int i = 1; i <= 2; ++i)
          if (quad == QGaussPyramid<dim>(i))
            {
              QGauss<dim - 1>        quad(i);
              QGaussSimplex<dim - 1> tri(i);

              return {
                ReferenceCells::Pyramid,
                dealii::hp::QCollection<dim - 1>(quad, tri, tri, tri, tri)};
            }

      if (do_assert)
        AssertThrow(false, ExcNotImplemented());

      return {ReferenceCells::Invalid, dealii::hp::QCollection<dim - 1>()};
    }



    template <int dim>
    inline std::pair<Quadrature<dim - 1>, Quadrature<dim - 1>>
    get_unique_face_quadratures(const Quadrature<dim> &quad)
    {
      if (dim == 2 || dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == QGaussSimplex<dim>(i))
            {
              if (dim == 2)
                return {QGaussSimplex<dim - 1>(i), Quadrature<dim - 1>()};
              else
                return {Quadrature<dim - 1>(), QGaussSimplex<dim - 1>(i)};
            }

      if (dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == QGaussWedge<dim>(i))
            return {QGauss<dim - 1>(i), QGaussSimplex<dim - 1>(i)};

      if (dim == 3)
        for (unsigned int i = 1; i <= 2; ++i)
          if (quad == QGaussPyramid<dim>(i))
            return {QGauss<dim - 1>(i), QGaussSimplex<dim - 1>(i)};

      AssertThrow(false, ExcNotImplemented());

      return {QGauss<dim - 1>(1), QGauss<dim - 1>(1)};
    }
  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
