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

#ifndef dealii_boost_adaptor_segment_h
#define dealii_boost_adaptor_segment_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/boost_adaptors/point.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/geometries/segment.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

/**
 * An alias for boost::geometry::model::segment that uses the deal.II
 * Point class.
 *
 * @author Luca Heltai, 2018
 */
template <int dim>
using Segment = boost::geometry::model::segment<Point<dim>>;

DEAL_II_NAMESPACE_CLOSE

#endif
