// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_boost_adaptor_segment_h
#define dealii_boost_adaptor_segment_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/boost_adaptors/point.h>

#include <boost/geometry/geometries/segment.hpp>

DEAL_II_NAMESPACE_OPEN

/**
 * An alias for boost::geometry::model::segment that uses the deal.II
 * Point class.
 */
template <int dim>
using Segment = boost::geometry::model::segment<Point<dim>>;

DEAL_II_NAMESPACE_CLOSE

#endif
