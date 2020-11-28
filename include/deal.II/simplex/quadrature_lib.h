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

#ifndef dealii_tet_quadrature_lib_h
#define dealii_tet_quadrature_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN

namespace Simplex
{
  /**
   * Integration rule for simplex entities.
   *
   * Users specify a number `n_points_1D` as an indication of what polynomial
   * degree to be integrated exactly, similarly to the number of points in a
   * QGauss quadrature object, even though the present quadrature formula is not
   * a tensor product. The given value is translated for n_points_1D=1,2,3 to
   * following number of quadrature points for 2D and 3D:
   *   - 2D: 1, 3, 7
   *   - 3D: 1, 4, 10
   *
   * For 1D, the quadrature rule degenerates to a `QGauss<1>(n_points_1D)`.
   *
   * @ingroup simplex
   */
  template <int dim>
  class QGauss : public QSimplex<dim>
  {
  public:
    /**
     * Constructor taking the number of quadrature points in 1D direction
     * @p n_points_1D.
     */
    explicit QGauss(const unsigned int n_points_1D);
  };

  /**
   * Integration rule for wedge entities.
   */
  template <int dim>
  class QGaussWedge : public Quadrature<dim>
  {
  public:
    /**
     * Users specify a number `n_points_1D` as an indication of what polynomial
     * degree to be integrated exactly. For details, see the comments of
     * Simplex::QGauss.
     */
    explicit QGaussWedge(const unsigned int n_points_1D);
  };

} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE

#endif
