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

#ifndef dealii__manifold_tools_h
#define dealii__manifold_tools_h


/*----------------------------   manifold_tools.h     -----------------*/

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/point.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold.h>

DEAL_II_NAMESPACE_OPEN

/**
 * CompositionManifold.  Take two ChartManifold objects, and make
 * their composition. The CompositionManifold object is a
 * ChartManifold going from the chart of the first ChartManifold to
 * the embedding space of the second ChartManifold. If the first
 * ChartManifold is periodic, so is the resulting ChartManifold, with
 * the same periodicity. Notice however that periodicity on the second
 * ChartManifold is lost during the composition.
 *
 * This class only works for dim <= chartdim <= intermediate_spacedim
 * <= spacedim. If you try to instantiate anything different, an
 * Exception will be thrown in one of the ChartManifold classes that
 * violates this condition.
 *
 * Given the ChartManifold F and the ChartManifold G, this class
 * represents the composition of G after F.
 *
 * The template parameters have the following meaning:
 *
 * @tparam dim The dimension of the resulting ChartManifold
 * @tparam spacedim The space dimension of the resulting ChartManifold
 * @tparam chartdim The chart dimension of the resulting ChartManifold
 * @tparam intermediate_dim The space dimension of the first ChartManifold
 * @tparam dim1 The dimension of the first ChartManifold, which coincides also
 * with the chart dimension of the second ChartManifold
 * @tparam dim2 The dimension of the second ChartManifold
 *
 * @ingroup manifold
 * @author Luca Heltai, Timo Heister, 2016
 */
template <int dim, int spacedim=dim, int chartdim=dim,
          int intermediate_dim=dim, int dim1=dim, int dim2=dim>
class CompositionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:

  /**
   * Construct the composition of the two given manifolds.
   */
  CompositionManifold(const ChartManifold<dim1, intermediate_dim, chartdim> &F,
                      const ChartManifold<dim2, spacedim, intermediate_dim> &G);

  /**
   * Pull back the given point in spacedim to the Euclidean chartdim
   * dimensional space. This function calls the pull_back() function
   * of G, and then the pull_back() function of F.
   */
  virtual
  Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const;


  /**
   * Push forward the chartdim dimensional point to a spacedim
   * Euclidean point. The function calls first the push_forward() of
   * F, and then the push_foward() of G.
   */
  virtual
  Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const;

  /**
   * Return the derivative of the composition of G after F.
   */
  virtual
  DerivativeForm<1,chartdim,spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

private:
  /**
   * The first ChartManifold.
   */
  SmartPointer<const ChartManifold<dim1, intermediate_dim, chartdim>,
               CompositionManifold<dim,spacedim,chartdim,dim1,dim2,intermediate_dim> > F;


  /**
   * The second ChartManifold.
   */
  SmartPointer<const ChartManifold<dim2, spacedim, intermediate_dim>,
               CompositionManifold<dim,spacedim,chartdim,dim1,dim2,intermediate_dim> > G;
};

DEAL_II_NAMESPACE_CLOSE

#endif
