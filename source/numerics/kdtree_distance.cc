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
unsigned int KDTreeDistance<dim>::get_points_within_ball(const Point<dim> &center, const double &radius,
                                                         std::vector<std::pair<unsigned int, double> > &matches,
                                                         bool sorted) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

  Assert(radius > 0,
         ExcMessage("Radius is expected to be positive."));

  nanoflann::SearchParams params;
  params.sorted = sorted;
  return kdtree->radiusSearch(&center[0], radius, matches, params);
}

template<int dim>
void KDTreeDistance<dim>::get_closest_points(const Point<dim> &target,
                                             std::vector<unsigned int> &indices,
                                             std::vector<double> &distances) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());
  AssertDimension(indices.size(), distances.size());

  kdtree->knnSearch(&target[0], indices.size(), &indices[0], &distances[0]);
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
