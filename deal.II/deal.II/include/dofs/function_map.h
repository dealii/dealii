//----------------------------  function_map.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_map.h  ---------------------------
#ifndef __deal2__function_map_h
#define __deal2__function_map_h

#include <map>

template <int dim> class Function;



/**
 * Declare a data type which denotes a mapping between a boundary
 * indicator and the function denoting the boundary values on this
 * part of the boundary. This type is required in many functions where
 * depending on boundary indicators different functions are used,
 * e.g. for boundary value interpolation, etc. As a variable of this
 * type is often declared in certain contexts where these functions
 * are used, this variable is also used by some functions that are not
 * actually interested in the function pointer itself, but only in the
 * list of selected boundary indicators.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
struct FunctionMap
{
				     /**
				      * Declare the type as discussed
				      * above. Since we can't name it
				      * @p{FunctionMap} (as that would
				      * ambiguate a possible
				      * constructor of this class),
				      * name it in the fashion of the
				      * STL local typedefs.
				      */
    typedef typename std::map<unsigned char, const Function<dim>*> type;
};


#endif
