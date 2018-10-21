// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#include <deal.II/numerics/kdtree.h>

#ifdef DEAL_II_WITH_NANOFLANN

#  include <deal.II/base/std_cxx14/memory.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
KDTree<dim>::KDTree(const unsigned int             max_leaf_size,
                    const std::vector<Point<dim>> &pts)
  : max_leaf_size(max_leaf_size)
{
  if (pts.size() > 0)
    set_points(pts);
}



template <int dim>
std::vector<std::pair<unsigned int, double>>
KDTree<dim>::get_points_within_ball(const Point<dim> &center,
                                    const double      radius,
                                    bool              sorted) const
{
  DEAL_II_Assert(adaptor, ExcNotInitialized());
  DEAL_II_Assert(kdtree, ExcInternalError());

  DEAL_II_Assert(radius > 0, ExcMessage("Radius is expected to be positive."));

  nanoflann::SearchParams params;
  params.sorted = sorted;

  std::vector<std::pair<unsigned int, double>> matches;
  kdtree->radiusSearch(&center[0], radius, matches, params);

  return matches;
}



template <int dim>
std::vector<std::pair<unsigned int, double>>
KDTree<dim>::get_closest_points(const Point<dim> & target,
                                const unsigned int n_points) const
{
  DEAL_II_Assert(adaptor, ExcNotInitialized());
  DEAL_II_Assert(kdtree, ExcInternalError());

  // get the information out of nanoflann
  std::vector<unsigned int> indices(n_points);
  std::vector<double>       distances(n_points);

  kdtree->knnSearch(&target[0], n_points, &indices[0], &distances[0]);

  // convert it to the format we want to return
  std::vector<std::pair<unsigned int, double>> matches(n_points);
  for (unsigned int i = 0; i < n_points; ++i)
    matches[i] = std::make_pair(indices[i], distances[i]);

  return matches;
}



template <int dim>
void
KDTree<dim>::set_points(const std::vector<Point<dim>> &pts)
{
  DEAL_II_Assert(pts.size() > 0,
                 ExcMessage("Expecting a non zero set of points."));
  adaptor = std_cxx14::make_unique<PointCloudAdaptor>(pts);
  kdtree  = std_cxx14::make_unique<NanoFlannKDTree>(
    dim, *adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));
  kdtree->buildIndex();
}


template class KDTree<1>;
template class KDTree<2>;
template class KDTree<3>;

DEAL_II_NAMESPACE_CLOSE

#endif
