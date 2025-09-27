// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_util_h
#define dealii_matrix_free_util_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/q_collection.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * Given a quadrature rule @p quad defined on a cell, return the type of
     * the cell and a collection of lower-dimensional quadrature rules that
     * are defined on each face.
     */
    template <int dim>
    inline std::pair<ReferenceCell, dealii::hp::QCollection<dim - 1>>
    get_face_quadrature_collection(const Quadrature<dim> &quad,
                                   const bool             do_assert = true)
    {
      if (dim == 2 || dim == 3)
        {
          for (unsigned int i = 1; i <= 4; ++i)
            if (quad == QGaussSimplex<dim>(i))
              return {ReferenceCells::get_simplex<dim>(),
                      dealii::hp::QCollection<dim - 1>(
                        QGaussSimplex<dim - 1>(i))};

          for (unsigned int i = 1; i <= 5; ++i)
            if (quad == QWitherdenVincentSimplex<dim>(i))
              return {ReferenceCells::get_simplex<dim>(),
                      dealii::hp::QCollection<dim - 1>(
                        QWitherdenVincentSimplex<dim - 1>(i))};

          for (unsigned int i = 1; i <= 3; ++i)
            {
              const FE_SimplexP<dim> fe(i);
              if (quad == Quadrature<dim>(fe.get_unit_support_points()))
                return {ReferenceCells::get_simplex<dim>(),
                        dealii::hp::QCollection<dim - 1>(Quadrature<dim - 1>(
                          fe.get_unit_face_support_points()))};
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

      // note: handle hypercubes last since normally this function is not
      // called for hypercubes
      for (unsigned int i = 1; i <= 5; ++i)
        if (quad == QGauss<dim>(i))
          return {ReferenceCells::get_hypercube<dim>(),
                  dealii::hp::QCollection<dim - 1>(QGauss<dim - 1>(i))};

      if (do_assert)
        AssertThrow(false, ExcNotImplemented());

      return {ReferenceCells::Invalid, dealii::hp::QCollection<dim - 1>()};
    }



    /**
     * Return face quadrature rules. In contrast to
     * get_face_quadrature_collection(), it does not return the quadrature
     * face for each face but returns one of each type. The first entry
     * of the returned pair might contain a quadrature rule defined on lines and
     * quadrilaterals, while the second entry might contain a quadrature rule
     * defined on a triangle.
     */
    template <int dim>
    inline std::pair<Quadrature<dim - 1>, Quadrature<dim - 1>>
    get_unique_face_quadratures(const Quadrature<dim> &quad)
    {
      AssertThrow(
        quad.size() > 0,
        ExcMessage(
          "There is nothing useful you can do with a MatrixFree/FEEvaluation "
          "object when using a quadrature formula with zero "
          "quadrature points!"));

      if (dim == 2 || dim == 3)
        {
          for (unsigned int i = 1; i <= 4; ++i)
            if (quad == QGaussSimplex<dim>(i))
              {
                if (dim == 2)
                  return {QGaussSimplex<dim - 1>(i), // line!
                          Quadrature<dim - 1>()};
                else
                  return {Quadrature<dim - 1>(), QGaussSimplex<dim - 1>(i)};
              }

          for (unsigned int i = 1; i <= 5; ++i)
            if (quad == QWitherdenVincentSimplex<dim>(i))
              {
                if (dim == 2)
                  return {QWitherdenVincentSimplex<dim - 1>(i), // line!
                          Quadrature<dim - 1>()};
                else
                  return {Quadrature<dim - 1>(),
                          QWitherdenVincentSimplex<dim - 1>(i)};
              }

          for (unsigned int i = 1; i <= 3; ++i)
            {
              const FE_SimplexP<dim> fe(i);
              if (quad == Quadrature<dim>(fe.get_unit_support_points()))
                {
                  if (dim == 2)
                    return {Quadrature<dim - 1>(
                              fe.get_unit_face_support_points()), // line!
                            Quadrature<dim - 1>()};
                  else
                    return {Quadrature<dim - 1>(),
                            Quadrature<dim - 1>(
                              fe.get_unit_face_support_points())};
                }
            }
        }

      if (dim == 3)
        for (unsigned int i = 1; i <= 3; ++i)
          if (quad == QGaussWedge<dim>(i))
            return {QGauss<dim - 1>(i), QGaussSimplex<dim - 1>(i)};

      if (dim == 3)
        for (unsigned int i = 1; i <= 2; ++i)
          if (quad == QGaussPyramid<dim>(i))
            return {QGauss<dim - 1>(i), QGaussSimplex<dim - 1>(i)};

      // note: handle hypercubes last since normally this function is not
      // called for hypercubes
      for (unsigned int i = 1; i <= 5; ++i)
        if (quad == QGauss<dim>(i))
          return {QGauss<dim - 1>(i), Quadrature<dim - 1>()};

      AssertThrow(false, ExcNotImplemented());

      return {Quadrature<dim - 1>(), Quadrature<dim - 1>()};
    }

    inline DEAL_II_ALWAYS_INLINE unsigned int
    indicate_power_of_two(const unsigned int vectorization_length)
    {
      unsigned int vectorization_length_bits = 0;
      unsigned int my_length                 = vectorization_length;
      while (my_length >>= 1)
        ++vectorization_length_bits;
      return 1 << vectorization_length_bits;
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
