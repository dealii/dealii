//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_lib_dg_h
#define __deal2__fe_lib_dg_h


#include <base/config.h>
#include<fe/fe_dgq.h>





/**
 * This class is an abbreviation for the @ref{FE_DGQ} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_DGQ}) is fixed to zero, i.e. piecewise constant elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_DGQ} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEDG_Q0 : public FE_DGQ<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEDG_Q0(): FE_DGQ<dim>(0)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_DGQ} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_DGQ}) is fixed to one, i.e. piecewise @p{d}-linear elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_DGQ} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEDG_Q1 : public FE_DGQ<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEDG_Q1(): FE_DGQ<dim>(1)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_DGQ} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_DGQ}) is fixed to two, i.e. piecewise @p{d}-quadratic elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_DGQ} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEDG_Q2 : public FE_DGQ<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEDG_Q2(): FE_DGQ<dim>(2)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_DGQ} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_DGQ}) is fixed to three, i.e. piecewise @p{d}-cubic elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_DGQ} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEDG_Q3 : public FE_DGQ<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEDG_Q3(): FE_DGQ<dim>(3)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_DGQ} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_DGQ}) is fixed to four, i.e. piecewise @p{d}-quartic elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_DGQ} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEDG_Q4 : public FE_DGQ<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEDG_Q4(): FE_DGQ<dim>(4)       {};
};


#endif
