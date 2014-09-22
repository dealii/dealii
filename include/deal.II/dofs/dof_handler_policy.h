// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
        distribute_dofs (dealii::DoFHandler<dim,spacedim> &dof_handler) const = 0;

        /**
         * Distribute the multigrid dofs on each level
         */
        virtual
        void
        distribute_mg_dofs (dealii::DoFHandler<dim,spacedim> &dof_handler,
                            std::vector<NumberCache> &number_caches) const = 0;

        /**
         * Renumber degrees of freedom as
         * specified by the first argument.
         */
        virtual
        NumberCache
        renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
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
        distribute_dofs (dealii::DoFHandler<dim,spacedim> &dof_handler) const;

        /**
         * Distribute multigrid DoFs.
         */
        virtual
        void
        distribute_mg_dofs (dealii::DoFHandler<dim,spacedim> &dof_handler,
                            std::vector<NumberCache> &number_caches) const;

        /**
         * Renumber degrees of freedom as
         * specified by the first argument.
         */
        virtual
        NumberCache
        renumber_dofs (const std::vector<types::global_dof_index>  &new_numbers,
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
        distribute_dofs (dealii::DoFHandler<dim,spacedim> &dof_handler) const;

        /**
         * Distribute multigrid DoFs.
         */
        virtual
        void
        distribute_mg_dofs (dealii::DoFHandler<dim,spacedim> &dof_handler,
                            std::vector<NumberCache> &number_caches) const;

        /**
         * Renumber degrees of freedom as
         * specified by the first argument.
         */
        virtual
        NumberCache
        renumber_dofs (const std::vector<types::global_dof_index>  &new_numbers,
                       dealii::DoFHandler<dim,spacedim> &dof_handler) const;
      };
    }
  }
}



DEAL_II_NAMESPACE_CLOSE

/*----------------------------   dof_handler_policy.h     ---------------------------*/
#endif
/*----------------------------   dof_handler_policy.h     ---------------------------*/
