// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Check that we can construct a tree with dealii::Point

#include <deal.II/base/patterns.h>

#include <deal.II/boost_adaptors/point.h>

#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>

#include <algorithm>

#include "../tests.h"

using Patterns::Tools::to_string;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

template <int dim, int leaf_size = 16>
using LinearPointRtree = bgi::rtree<Point<dim>, bgi::linear<leaf_size>>;

template <int dim, int leaf_size = 16>
using QuadraticPointRtree = bgi::rtree<Point<dim>, bgi::quadratic<leaf_size>>;

template <int dim, int leaf_size = 16>
using RstarPointRtree = bgi::rtree<Point<dim>, bgi::rstar<leaf_size>>;

int
main(int argc, char **argv)
{
  initlog();

  // N random points in the box.
  const unsigned int    N = 100;
  std::vector<Point<2>> points(N);
  std::generate(points.begin(), points.end(), []() {
    return random_point<2>();
  });

  LinearPointRtree<2> linear_tree;
  linear_tree.insert(points.begin(), points.end());

  QuadraticPointRtree<2> quadratic_tree;
  quadratic_tree.insert(points.begin(), points.end());

  RstarPointRtree<2> rstar_tree;
  rstar_tree.insert(points.begin(), points.end());

  AssertDimension(linear_tree.size(), N);
  AssertDimension(quadratic_tree.size(), N);
  AssertDimension(rstar_tree.size(), N);

  // Now take the middle point, and ask for the k nearest neighbors
  Point<2> p(.5, .5);
  // Return only the first two points, so that the order is unique
  const unsigned int    k = 2;
  std::vector<Point<2>> nearest_linear;
  std::vector<Point<2>> nearest_quadratic;
  std::vector<Point<2>> nearest_rstar;

  linear_tree.query(bgi::nearest(p, k), std::back_inserter(nearest_linear));
  quadratic_tree.query(bgi::nearest(p, k),
                       std::back_inserter(nearest_quadratic));
  rstar_tree.query(bgi::nearest(p, k), std::back_inserter(nearest_rstar));

  AssertDimension(nearest_linear.size(), k);
  AssertDimension(nearest_quadratic.size(), k);
  AssertDimension(nearest_rstar.size(), k);

  if (std::equal(nearest_linear.begin(),
                 nearest_linear.end(),
                 nearest_quadratic.begin()))
    deallog << "OK";
  else
    {
      deallog << "Not OK: check " << std::endl
              << to_string(nearest_linear) << std::endl
              << to_string(nearest_quadratic) << std::endl;
    }
  if (std::equal(nearest_linear.begin(),
                 nearest_linear.end(),
                 nearest_rstar.begin()))
    deallog << "OK";
  else
    {
      deallog << "Not OK: check " << std::endl
              << to_string(nearest_linear) << std::endl
              << to_string(nearest_rstar) << std::endl;
    }
}
