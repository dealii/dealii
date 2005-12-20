//----------------------------  select_fe_values.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
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
template <int> class hpDoFHandler;

template <int> class FEValues;
template <int> class FEFaceValues;
template <int> class FESubfaceValues;

template <int> class hpFEValues;
template <int> class hpFEFaceValues;
template <int> class hpFESubfaceValues;



template <typename> struct SelectFEValues;

template <int dim> struct SelectFEValues<DoFHandler<dim> >
{
    typedef FEValues<dim>        FEValues;
    typedef FEFaceValues<dim>    FEFaceValues;
    typedef FESubfaceValues<dim> FESubfaceValues;
};


template <int dim> struct SelectFEValues<hpDoFHandler<dim> >
{
    typedef hpFEValues<dim>        FEValues;
    typedef hpFEFaceValues<dim>    FEFaceValues;
    typedef hpFESubfaceValues<dim> FESubfaceValues;
};




#endif
