// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__function_map_h
#define __deal2__function_map_h

#include <deal.II/base/config.h>
#include <map>

DEAL_II_NAMESPACE_OPEN

template <int spacedim, typename Number> class Function;



/**
 * This class declares a local typedef that denotes a mapping between a boundary indicator
 * (see @ref GlossBoundaryIndicator) that is used to describe what kind of boundary
 * condition holds on a particular piece of the boundary,
 * and the function describing the actual function that provides the boundary
 * values on this part of the boundary. This type is required in many functions
 * in the library where, for example, we need to know about the functions $h_i(\mathbf x)$
 * used in boundary conditions
 * @f{align*}
 *   \mathbf n \cdot \nabla u = h_i \qquad \qquad \text{on}\ \Gamma_i\subset\partial\Omega.
 * @f}
 * An example is the function KellyErrorEstimator::estimate() that allows us
 * to provide a set of functions $h_i$ for all those boundary indicators $i$ for
 * which the boundary condition is supposed to be of Neumann type. Of course,
 * the same kind of principle can be applied to cases where we care about
 * Dirichlet values, where one needs to provide a map from boundary indicator $i$
 * to Dirichlet function $h_i$ if the boundary conditions are given as
 * @f{align*}
 *   u = h_i \qquad \qquad \text{on}\ \Gamma_i\subset\partial\Omega.
 * @f}
 * This is, for example, the case for the VectorTools::interpolate() functions.
 *
 * Tutorial programs step-6, step-7 and step-8 show examples of how to use
 * function arguments of this type in situations where we actually have an empty
 * map (i.e., we want to describe that <i>no</i> part of the boundary is a
 * Neumann boundary). step-16 actually uses it in a case where one of the
 * parts of the boundary uses a boundary indicator for which we want to use
 * a function object.
 *
 * It seems odd at first to declare this typedef inside a class, rather than
 * declaring a typedef at global scope. The reason is that C++ does not allow
 * to define templated typedefs, where here in fact we want a typdef that
 * depends on the space dimension. (Defining templated typedefs is something that
 * is possible starting with the C++11 standard, but that wasn't possible within
 * the C++98 standard in place when this programming pattern was conceived.)
 *
 * @ingroup functions
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim,typename Number=double>
struct FunctionMap
{
  /**
   * Declare the type as discussed
   * above. Since we can't name it
   * FunctionMap (as that would
   * ambiguate a possible
   * constructor of this class),
   * name it in the fashion of the
   * STL local typedefs.
   */
  typedef std::map<types::boundary_id, const Function<dim,Number>*> type;
};

DEAL_II_NAMESPACE_CLOSE

#endif
