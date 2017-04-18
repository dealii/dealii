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
 * CompositionManifold.
 *
 * @ingroup manifold
 * @author Luca Heltai, Timo Heister, 2016
 */
template <int dim, int spacedim=dim, int chartdim=dim, int dim1=dim, int dim2=dim, int spacedim1=dim>
class CompositionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:

  /**
   * Take two ChartManifold() objects, and make their composition. The
   * CompositionManifold() object is a ChartManifold() going from the chart of
   * the first ChartManifold to the embedding space of the second ChartManifold().
   *
   * The dim parameters of both manifolds are simply ignored. This class only works
   * for chartdim <= spacedim1 <= spacedim. If you try to instantiate anything
   * different, an Exception will be thrown in one of the ChartManifold classes
   * that violates this condition.
   *
   * Given the ChartManifold F and the ChartManifold G, this class represents the
   * composition of G after F.
   */
  CompositionManifold(const ChartManifold<dim1, spacedim1, chartdim> &F,
                      const ChartManifold<dim2, spacedim, spacedim1> &G);

  /**
   * Destructor. Does nothing here, but needs to be declared virtual to make
   * class hierarchies derived from this class possible.
   */
  virtual ~CompositionManifold ();


  /**
   * Pull back the given point in spacedim to the Euclidean chartdim
   * dimensional space. This function calls the pull_back() function of
   * G, and then the pull_back() function of F.
   */
  virtual
  Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const;


  /**
   * Push forward the chartdim dimensional point to a spacedim Euclidean
   * point. The function calls first the push_forward() of F, and then the
   * push_foward() of G.
   */
  virtual
  Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const;

  /**
   * Return the derivative of the composition of G after F. This function
   * calls first the push_forward_gradient() and the push_forward() of
   * F. With the resulting Point, it calls the push_forward_gradient() of
   * G, and then it multiplies DF and DG.
   */
  virtual
  DerivativeForm<1,chartdim,spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

private:
  /**
   * The first ChartManifold.
   */
  const ChartManifold<dim1, spacedim1, chartdim> &F;


  /**
   * The second ChartManifold.
   */
  const ChartManifold<dim2, spacedim, spacedim1> &G;
};

DEAL_II_NAMESPACE_CLOSE

#endif
