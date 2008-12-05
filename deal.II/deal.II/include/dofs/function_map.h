//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__function_map_h
#define __deal2__function_map_h

#include <base/config.h>
#include <map>

DEAL_II_NAMESPACE_OPEN

template <int spacedim> class Function;



/**
 * Declare a data type which denotes a mapping between a boundary indicator
 * and the function denoting the boundary values on this part of the
 * boundary. This type is required in many functions where depending on the
 * boundary indicator, different functions are used. An example is boundary
 * value interpolation.
 *
 * It seems odd at first to declare this typedef inside a class, rather than
 * declaring a typedef at global scope. The reason is that C++ does not allow
 * to define templated typedefs, where here in fact we want a typdef that
 * depends on the space dimension.
 *
 * @ingroup functions
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
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
    typedef std::map<unsigned char, const Function<dim>*> type;
};

DEAL_II_NAMESPACE_CLOSE

#endif
