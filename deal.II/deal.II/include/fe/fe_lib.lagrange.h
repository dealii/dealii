//---------------------------------------------------------------
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
//---------------------------------------------------------------
#ifndef __deal2__fe_lib_lagrange_h
#define __deal2__fe_lib_lagrange_h


#include<fe/fe_q.h>


template<int dim> class FEQ1: public FE_Q<dim>
{
  public:
    FEQ1(): FE_Q<dim>(1) 
      {};
};


template<int dim> class FEQ2: public FE_Q<dim>
{
  public:
    FEQ2(): FE_Q<dim>(2) 
      {};
};


template<int dim> class FEQ3: public FE_Q<dim>
{
  public:
    FEQ3(): FE_Q<dim>(3) 
      {};
};


template<int dim> class FEQ4: public FE_Q<dim>
{
  public:
    FEQ4(): FE_Q<dim>(4) 
      {};
};


#endif
