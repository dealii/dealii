//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_handler_policy_h
#define __deal2__dof_handler_policy_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>

#include <vector>
#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <int, int> class FiniteElement;
template <int, int> class DoFHandler;


namespace internal
{
  namespace DoFHandler
  {
    struct NumberCache;

				     /**
				      * A namespace in which we define
				      * classes that describe how to
				      * distribute and renumber
				      * degrees of freedom.
				      */
    namespace Policy
    {
      struct Implementation;

				       /**
					* A class that implements policies for
					* how the DoFHandler::distribute_dofs
					* and DoFHandler::renumber_dofs
					* functions should work.
					*/
      template <int dim, int spacedim>
      class PolicyBase
      {
	public:
					   /**
					    * Destructor.
					    */
	  virtual ~PolicyBase ();

					   /**
					    * Distribute degrees of freedom on
					    * the object given as last argument.
					    */
	  virtual
	  NumberCache
	  distribute_dofs (const unsigned int        offset,
			   dealii::DoFHandler<dim,spacedim> &dof_handler) const = 0;

					   /**
					    * Renumber degrees of freedom as
					    * specified by the first argument.
					    */
	  virtual
	  NumberCache
	  renumber_dofs (const std::vector<unsigned int>  &new_numbers,
			 dealii::DoFHandler<dim,spacedim> &dof_handler) const = 0;
      };


				       /**
					* This class implements the
					* default policy for sequential
					* operations, i.e. for the case where
					* all cells get degrees of freedom.
					*/
      template <int dim, int spacedim>
      class Sequential : public PolicyBase<dim,spacedim>
      {
	public:
					   /**
					    * Distribute degrees of freedom on
					    * the object given as last argument.
					    */
	  virtual
	  NumberCache
	  distribute_dofs (const unsigned int        offset,
			   dealii::DoFHandler<dim,spacedim> &dof_handler) const;

					   /**
					    * Renumber degrees of freedom as
					    * specified by the first argument.
					    */
	  virtual
	  NumberCache
	  renumber_dofs (const std::vector<unsigned int>  &new_numbers,
			 dealii::DoFHandler<dim,spacedim> &dof_handler) const;
      };


				       /**
					* This class implements the
					* policy for operations when
					* we use a
					* parallel::distributed::Triangulation
					* object.
					*/
      template <int dim, int spacedim>
      class ParallelDistributed : public PolicyBase<dim,spacedim>
      {
	public:
					   /**
					    * Distribute degrees of freedom on
					    * the object given as last argument.
					    */
	  virtual
	  NumberCache
	  distribute_dofs (const unsigned int        offset,
			   dealii::DoFHandler<dim,spacedim> &dof_handler) const;

					   /**
					    * Renumber degrees of freedom as
					    * specified by the first argument.
					    */
	  virtual
	  NumberCache
	  renumber_dofs (const std::vector<unsigned int>  &new_numbers,
			 dealii::DoFHandler<dim,spacedim> &dof_handler) const;
      };
    }
  }
}



DEAL_II_NAMESPACE_CLOSE

/*----------------------------   dof_handler_policy.h     ---------------------------*/
#endif
/*----------------------------   dof_handler_policy.h     ---------------------------*/
