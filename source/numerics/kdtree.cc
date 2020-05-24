// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

  Assert(radius > 0, ExcMessage("Radius is expected to be positive."));

  nanoflann::SearchParams params;
  params.sorted = sorted;

  std::vector<std::pair<unsigned int, double>> matches;
#  if NANOFLANN_VERSION < 0x130
  kdtree->radiusSearch(center.begin_raw(), radius, matches, params);
#  else
  // nanoflann 1.3 performs distance comparisons with squared distances, so
  // square the radius before we query and square root after:
  kdtree->radiusSearch(center.begin_raw(), radius * radius, matches, params);
  for (std::pair<unsigned int, double> &match : matches)
    match.second = std::sqrt(match.second);
#  endif

  return matches;
}



template <int dim>
std::vector<std::pair<unsigned int, double>>
KDTree<dim>::get_closest_points(const Point<dim> & target,
                                const unsigned int n_points) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

  // get the information out of nanoflann
  std::vector<unsigned int> indices(n_points);
  std::vector<double>       distances(n_points);

  kdtree->knnSearch(target.begin_raw(),
                    n_points,
                    indices.data(),
                    distances.data());

  // convert it to the format we want to return
  std::vector<std::pair<unsigned int, double>> matches(n_points);
  for (unsigned int i = 0; i < n_points; ++i)
#  if NANOFLANN_VERSION < 0x130
    matches[i] = std::make_pair(indices[i], distances[i]);
#  else
    // nanoflann 1.3 performs distance comparisons with squared distances, so
    // take a square root:
    matches[i] = std::make_pair(indices[i], std::sqrt(distances[i]));
#  endif

  return matches;
}



template <int dim>
void
KDTree<dim>::set_points(const std::vector<Point<dim>> &pts)
{
  Assert(pts.size() > 0, ExcMessage("Expecting a non zero set of points."));
  adaptor = std::make_unique<PointCloudAdaptor>(pts);
  kdtree  = std::make_unique<NanoFlannKDTree>(
    dim, *adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));
  kdtree->buildIndex();
}


template class KDTree<1>;
template class KDTree<2>;
template class KDTree<3>;

DEAL_II_NAMESPACE_CLOSE

#endif
