#include <deal.II/numerics/kdtree_distance.h>

#ifdef DEAL_II_WITH_NANOFLANN

#include<deal.II/base/std_cxx14/memory.h>

DEAL_II_NAMESPACE_OPEN


template<int dim>
KDTreeDistance<dim>::KDTreeDistance(const unsigned int &max_leaf_size,
                                    const std::vector<Point<dim> > &pts)
  : max_leaf_size(max_leaf_size)
{
  if (pts.size() > 0)
    set_points(pts);
}


template<int dim>
std::vector<std::pair<unsigned int, double> > KDTreeDistance<dim>::get_points_within_ball(const Point<dim> &center,
                                                                                          const double &radius,
                                                                                          bool sorted) const
{
  std::vector<std::pair<unsigned int, double> > matches;
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

  Assert(radius > 0,
         ExcMessage("Radius is expected to be positive."));

  nanoflann::SearchParams params;
  params.sorted = sorted;
  kdtree->radiusSearch(&center[0], radius, matches, params);
  return matches;
}

template<int dim>
std::vector<std::pair<unsigned int, double> > KDTreeDistance<dim>::get_closest_points(const Point<dim> &target,
                                                                                      const unsigned int n_points) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());
  std::vector<unsigned int> indices(n_points);
  std::vector<double> distances(n_points);
  std::vector<std::pair<unsigned int, double> > matches(n_points);

  kdtree->knnSearch(&target[0], n_points, &indices[0], &distances[0]);
  for(unsigned int i=0; i<n_points; ++i)
    matches[i] = std::make_pair(indices[i], distances[i]);
  return matches;
}

template<int dim>
void KDTreeDistance<dim>::set_points(const std::vector<Point<dim> > &pts)
{
  Assert(pts.size() > 0, ExcMessage("Expecting a non zero set of points."));
  adaptor = std_cxx14::make_unique<PointCloudAdaptor>(pts);
  kdtree = std_cxx14::make_unique<KDTree>(dim, *adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));
  kdtree->buildIndex();
}


template class KDTreeDistance<1>;
template class KDTreeDistance<2>;
template class KDTreeDistance<3>;

DEAL_II_NAMESPACE_CLOSE

#endif
