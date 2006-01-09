//----------------------------  select_fe_values.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  select_fe_values.h  ---------------------------
#ifndef __deal2__select_fe_values_h
#define __deal2__select_fe_values_h

#include <base/config.h>

template <int> class DoFHandler;
namespace hp
{
  template <int> class DoFHandler;
}

template <int> class FEValues;
template <int> class FEFaceValues;
template <int> class FESubfaceValues;

namespace hp
{
  template <int> class FEValues;
  template <int> class FEFaceValues;
  template <int> class FESubfaceValues;
}



/**
 * 
 * @ingroup hp
 */  
template <typename> struct SelectFEValues;

/**
 * 
 * @ingroup hp
 */  
template <int dim> struct SelectFEValues<DoFHandler<dim> >
{
    typedef FEValues<dim>        FEValues;
    typedef FEFaceValues<dim>    FEFaceValues;
    typedef FESubfaceValues<dim> FESubfaceValues;
};


/**
 * 
 * @ingroup hp
 */  
template <int dim> struct SelectFEValues<hp::DoFHandler<dim> >
{
    typedef hp::FEValues<dim>        FEValues;
    typedef hp::FEFaceValues<dim>    FEFaceValues;
    typedef hp::FESubfaceValues<dim> FESubfaceValues;
};




#endif
