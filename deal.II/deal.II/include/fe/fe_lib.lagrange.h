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


#include <base/config.h>
#include<fe/fe_q.h>




/**
 * This class is an abbreviation for the @ref{FE_Q} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_Q}) is fixed to one, i.e. @p{d}-linear elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_Q} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEQ1 : public FE_Q<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEQ1(): FE_Q<dim>(1)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_Q} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_Q}) is fixed to two, i.e. @p{d}-quadratic elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_Q} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEQ2 : public FE_Q<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEQ2(): FE_Q<dim>(2)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_Q} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_Q}) is fixed to three, i.e. @p{d}-cubic elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_Q} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEQ3 : public FE_Q<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEQ3(): FE_Q<dim>(3)       {};
};


/**
 * This class is an abbreviation for the @ref{FE_Q} class, where the
 * polynomial degree (which is passed through the constructor of
 * @ref{FE_Q}) is fixed to four, i.e. @p{d}-quartic elements.
 *
 * This class is only here for backward compatibility and will go away
 * someday. Use the @ref{FE_Q} class instead.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, 2001
 */
template<int dim>
class FEQ4 : public FE_Q<dim>
{
  public:
				     /**
				      * Constructor. Simply pass
				      * correct polynomial degree to
				      * constructor of base class.
				      */
    FEQ4(): FE_Q<dim>(4)       {};
};


#endif
