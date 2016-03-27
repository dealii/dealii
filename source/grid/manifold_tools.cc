// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/manifold_tools.h>
#include <deal.II/fe/fe_q.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN


/* -------------------------- CompositionManifold --------------------- */

template <int dim, int spacedim, int chartdim, int intermediate_dim,
          int dim1, int dim2>
CompositionManifold<dim,spacedim,chartdim,intermediate_dim,dim1,dim2>::CompositionManifold
(const ChartManifold<dim1, intermediate_dim, chartdim> &F,
 const ChartManifold<dim2, spacedim, intermediate_dim> &G) :
  ChartManifold<dim,spacedim,chartdim>(F.get_periodicity()),
  F(&F),
  G(&G)
{
  // We don't know what to do with a periodicity in the second manifold, so
  // throw an assertion if the second manifold is periodic
  Assert(G.get_periodicity().norm() == 0.0,
         ExcMessage("The second manifold cannot be periodic."));
}



template <int dim, int spacedim, int chartdim, int intermediate_dim,
          int dim1, int dim2>
Point<chartdim>
CompositionManifold<dim,spacedim,chartdim,intermediate_dim,dim1,dim2>::pull_back(const Point<spacedim> &space_point) const
{
  return F->pull_back(G->pull_back(space_point));
}



template <int dim, int spacedim, int chartdim, int intermediate_dim,
          int dim1, int dim2>
Point<spacedim>
CompositionManifold<dim,spacedim,chartdim,intermediate_dim,dim1,dim2>::push_forward(const Point<chartdim> &chart_point) const
{
  return G->push_forward(F->push_forward(chart_point));
}



template <int dim, int spacedim, int chartdim, int intermediate_dim,
          int dim1, int dim2>
DerivativeForm<1,chartdim,spacedim>
CompositionManifold<dim,spacedim,chartdim,intermediate_dim,dim1,dim2>::push_forward_gradient(const Point<chartdim> &chart_point) const
{
  DerivativeForm<1,chartdim,intermediate_dim> DF =
    F->push_forward_gradient(chart_point);

  DerivativeForm<1,intermediate_dim,spacedim> DG =
    G->push_forward_gradient(F->push_forward(chart_point));

  DerivativeForm<1,chartdim,spacedim> DF_DG;

  for (unsigned int d=0; d<spacedim; ++d)
    for (unsigned int c=0; c<chartdim; ++c)
      for (unsigned int s=0; s<intermediate_dim; ++s)
        DF_DG[d][c] += DG[d][s]*DF[s][c];

  return DF_DG;
}



// explicit instantiations
#include "manifold_tools.inst"

DEAL_II_NAMESPACE_CLOSE

