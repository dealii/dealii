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
#ifndef __deal2__fe_lib_dg_h
#define __deal2__fe_lib_dg_h


#include<fe/fe_dgq.h>


template<int dim> class FEDG_Q0: public FE_DGQ<dim>
{
  public:
    FEDG_Q0(): FE_DGQ<dim>(0) 
      {};
};


template<int dim> class FEDG_Q1: public FE_DGQ<dim>
{
  public:
    FEDG_Q1(): FE_DGQ<dim>(1) 
      {};
};


template<int dim> class FEDG_Q2: public FE_DGQ<dim>
{
  public:
    FEDG_Q2(): FE_DGQ<dim>(2) 
      {};
};


template<int dim> class FEDG_Q3: public FE_DGQ<dim>
{
  public:
    FEDG_Q3(): FE_DGQ<dim>(3) 
      {};
};


template<int dim> class FEDG_Q4: public FE_DGQ<dim>
{
  public:
    FEDG_Q4(): FE_DGQ<dim>(4) 
      {};
};


#endif
