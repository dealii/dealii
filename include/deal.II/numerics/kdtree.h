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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_numerics_kdtree_h
#define dealii_numerics_kdtree_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_NANOFLANN

#  include <deal.II/base/point.h>

#  include <nanoflann.hpp>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

/**
 * A wrapper for the nanoflann library, used to compute the distance from a
 * collection of points, and to efficiently return nearest neighbors to a
 * target point. This class uses nanoflann to efficiently partition the
 * space in a $k$-dimensional tree. The cost of each query is then roughly of
 * order $\log(n)$, where $n$ is the number of points stored in this class.
 *
 * The wrapper provides methods that give access to some of the functionalities
 * of the nanoflann library, like searching the $p$ nearest neighbors of
 * a given point, or
 * searching the points that fall within a radius of a target point.
 *
 * > From wikipedia (https://en.wikipedia.org/wiki/K-d_tree):
 * >
 * > A k-d tree is a binary tree in which every node is a $k$-dimensional point.
 * > Every non-leaf node can be thought of as implicitly generating a splitting
 * > hyperplane that divides the space into two parts, known as half-spaces.
 * > Points to the left of this hyperplane are represented by the left subtree
 * of > that node and points right of the hyperplane are represented by the
 * right > subtree. The hyperplane direction is chosen in the following way:
 * every node > in the tree is associated with one of the $k$-dimensions, with
 * the hyperplane > perpendicular to that dimension's axis. So, for example, if
 * for a particular > split the "x" axis is chosen, all points in the subtree
 * with a smaller "x" > value than the node will appear in the left subtree and
 * all points with > larger "x" value will be in the right subtree. In such a
 * case, the > hyperplane would be set by the $x$-value of the point, and its
 * normal would be > the unit $x$-axis.
 *
 * @author Luca Heltai, 2017.
 */
template <int dim>
class KDTree
{
public:
  /**
   * Constructor.
   *
   * @param[in] max_leaf_size A number denoting how many points per leaf
   * are used in the kdtree algorithm.
   *
   * @param[in] pts A vector of points that are to be represented by
   * the current object. If no points are passed to this constructor
   * (or if the default value of the argument is used), then you have
   * to pass them later to this object by calling the set_points()
   * method.
   *
   * Access to any of the methods without first passing a reference to
   * a vector of points will result in an exception. Only a reference
   * to the points is stored, so you should make sure that the life of
   * the vector you pass is longer than the life of this class, or
   * you will get undefined behaviour.
   *
   * @warning If you change the contents of the vector of points that you
   * passed either to the constructor or to set_points(), remember to call
   * the set_points() method again. The tree and the index are
   * constructed only once when you pass the points (either at
   * construction time, or when you call set_points()). If you update
   * your points, and do not call set_points() again, then all following results
   * will likely be wrong.
   */
  KDTree(const unsigned int &           max_leaf_size = 10,
         const std::vector<Point<dim>> &pts = std::vector<Point<dim>>());


  /**
   * Adaptor class used internally by nanoflann. This class stores a reference
   * to the vector of points, and generates some helper functions for
   * nanoflann.
   */
  struct PointCloudAdaptor
  {
    /**
     * A typedef used by nanoflann.
     */
    typedef double coord_t;


    /**
     * Reference to the vector of points from which we want to compute
     * the distance.
     */
    const std::vector<Point<dim>> &points;


    /**
     * The constructor needs the vector of points from which we want to build
     * the tree.
     */
    PointCloudAdaptor(const std::vector<Point<dim>> &_points);


    /**
     * Return number of points in the data set (required by nanoflann).
     */
    size_t
    kdtree_get_point_count() const;


    /**
     * Return the L2 distance between points
     */
    coord_t
    kdtree_distance(const coord_t *p1,
                    const size_t   idx_p2,
                    const size_t   size) const;


    /**
     * Return the d-th component of the idx-th point in the class.
     */
    coord_t
    kdtree_get_pt(const size_t idx, const int d) const;


    /**
     * Optional bounding-box computation: return false to default to a
     * standard bbox computation loop.  Return true if the BBOX was
     * already computed by the class and returned in "bb" so it can be
     * avoided to redo it again.  Look at bb.size() to find out the
     * expected dimensionality (e.g. 2 or 3 for point clouds).
     */
    template <class BBOX>
    bool
    kdtree_get_bbox(BBOX &) const;
  };


  /**
   * A typedef for the actual KDTree object.
   */
  typedef typename nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloudAdaptor>,
    PointCloudAdaptor,
    dim,
    unsigned int>
    NanoFlannKDTree;


  /**
   * Store a reference to the passed points. After you called this
   * method, you can call the value() method to compute the minimum
   * distance between an evaluation point and the collection of points
   * you passed to this method, or the get_points_within_ball() and
   * the get_closest_points() methods.
   *
   * Notice that the constructor calls this method internally if you
   * pass it a non-empty vector of points.
   *
   * Whenever your points change, you should call this method again,
   * since this is the method responsible for building the index and
   * storing the actual tree internally. If you change your points and
   * don't call again this method, any function you call later will
   * happily return wrong values without you noticing.
   *
   * @param[in] pts A collection of points
   */
  void
  set_points(const std::vector<Point<dim>> &pts);


  /**
   * A const accessor to the @p i'th one among the underlying points.
   */
  const Point<dim> &operator[](const unsigned int i) const;


  /**
   * The number of points currently stored by this class.
   */
  unsigned int
  size() const;


  /**
   * Fill and return a vector with the indices and the distance of the points
   * that are at distance less than or equal to the given radius from
   * the target point.
   *
   * @param[in] target The target point
   * @param[in] radius The radius of the ball
   * @param[in] sorted If @p true, sort the output results in ascending order with respect to distance
   *
   * @return A vector of indices and distances to @p target of the matching points
   */
  std::vector<std::pair<unsigned int, double>>
  get_points_within_ball(const Point<dim> &target,
                         const double &    radius,
                         const bool        sorted = false) const;

  /**
   * Fill and return a vector with the indices and distances of the closest @p n_points
   * points to the given target point.
   *
   * @param[in] target The target point
   * @param[in] n_points The number of requested points
   *
   * @return A vector of pairs of indices and distances of the matching points
   */
  std::vector<std::pair<unsigned int, double>>
  get_closest_points(const Point<dim> & target,
                     const unsigned int n_points) const;

private:
  /**
   * Max number of points per leaf as set in the constructor.
   */
  const unsigned int max_leaf_size;


  /**
   * A point cloud adaptor, to be filled when set points is called.
   */
  std::unique_ptr<PointCloudAdaptor> adaptor;


  /**
   * The actual kdtree.
   */
  std::unique_ptr<NanoFlannKDTree> kdtree;
};


//------------ inline functions -------------
#  ifndef DOXYGEN

template <int dim>
inline unsigned int
KDTree<dim>::size() const
{
  if (adaptor)
    return adaptor->points.size();
  else
    return 0;
}



template <int dim>
inline const Point<dim> &KDTree<dim>::operator[](const unsigned int i) const
{
  AssertIndexRange(i, size());
  return adaptor->points[i];
}



template <int dim>
KDTree<dim>::PointCloudAdaptor::PointCloudAdaptor(
  const std::vector<Point<dim>> &_points) :
  points(_points)
{}



template <int dim>
inline size_t
KDTree<dim>::PointCloudAdaptor::kdtree_get_point_count() const
{
  return points.size();
}



template <int dim>
inline double
KDTree<dim>::PointCloudAdaptor::kdtree_get_pt(const size_t idx, int d) const
{
  AssertIndexRange(d, dim);
  return points[idx][d];
}



template <int dim>
template <class BBOX>
inline bool
KDTree<dim>::PointCloudAdaptor::kdtree_get_bbox(BBOX &) const
{
  return false;
}



template <int dim>
inline double
KDTree<dim>::PointCloudAdaptor::kdtree_distance(const double *p1,
                                                const size_t  idx_p2,
                                                const size_t  size) const
{
  AssertDimension(size, dim);
  double res = 0.0;
  for (size_t d = 0; d < size; ++d)
    res += (p1[d] - points[idx_p2][d]) * (p1[d] - points[idx_p2][d]);
  return std::sqrt(res);
}
#  endif

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_NANO_FLANN
#endif
