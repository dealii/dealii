// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


/**
 * @defgroup functions Functions
 *
 * Functions are used in various places in deal.II, for example to
 * describe boundary conditions, coefficients in equations, forcing
 * terms, or exact solutions. Since closed form expressions for
 * equations are often hard to pass along as function arguments,
 * deal.II uses the Function base class to describe these
 * objects. Essentially, the interface of this base class requires
 * derived classes to implement the ability to return the value of a
 * function at one or a list of particular locations, and possibly (if
 * needed) of gradients or second derivatives of the function. With
 * this, function objects can then be used by algorithms like
 * VectorTools::interpolate, VectorTools::project_boundary_values, and
 * other functions.
 *
 * Some functions are needed again and again, and are therefore
 * already provided in deal.II. This includes a function with a
 * constant value; a function that is zero everywhere, or a
 * vector-valued function for which only one vector component has a
 * particular value and all other components are zero. Some more
 * specialized functions are also defined in the Functions namespace.
 *
 *
 * <h3>Time dependent functions</h3>
 * 
 * For time dependent computations, boundary conditions and/or right
 * hand side functions may also change with time. Since at a given
 * time step one is usually only interested in the spatial dependence
 * of a function, it would be awkward if one had to pass a value for
 * the time variable to all methods that use function objects. For
 * example, the VectorTools::interpolate_boundary_values function
 * would have to take a time argument which it can use when it wants
 * to query the value of the boundary function at a given time
 * step. However, it would also have to do so if we are considering a
 * stationary problem, for which there is nothing like a time
 * variable.
 *
 * To circumvent this problem, function objects are always considered
 * spatial functions only. However, the Function class is derived from
 * the FunctionTime base class that stores a value for a time
 * variable, if so necessary. This way, one can define a function
 * object that acts as a spatial function but can do so internally by
 * referencing a particular time. In above example, one would set the
 * time of the function object to the present time step before handing
 * it off to the VectorTools::interpolate_boundary_values method.
 *
 * 
 * <h3>Tensor-valued functions</h3>
 *
 * The Function class is the most frequently used, but sometimes one needs a
 * function the values of which are tensors, rather than scalars. The
 * TensorFunction template can do this for you. Apart from the return type,
 * the interface is most the same as that of the Function class.
 */
