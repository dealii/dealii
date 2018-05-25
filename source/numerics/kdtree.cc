#include <deal.II/numerics/kdtree.h>

#ifdef DEAL_II_WITH_NANOFLANN

#  include <deal.II/base/std_cxx14/memory.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
KDTree<dim>::KDTree(const unsigned int &           max_leaf_size,
                    const std::vector<Point<dim>> &pts) :
  max_leaf_size(max_leaf_size)
{
  if (pts.size() > 0)
    set_points(pts);
}



template <int dim>
std::vector<std::pair<unsigned int, double>>
KDTree<dim>::get_points_within_ball(const Point<dim> &center,
                                    const double &    radius,
                                    bool              sorted) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

  Assert(radius > 0, ExcMessage("Radius is expected to be positive."));

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
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

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
  Assert(pts.size() > 0, ExcMessage("Expecting a non zero set of points."));
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
