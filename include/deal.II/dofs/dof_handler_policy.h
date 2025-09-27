// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dof_handler_policy_h
#define dealii_dof_handler_policy_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
#endif

namespace internal
{
  namespace DoFHandlerImplementation
  {
    struct NumberCache;

    /**
     * A namespace in which we define classes that describe how to distribute
     * and renumber degrees of freedom.
     */
    namespace Policy
    {
      struct Implementation;

      /**
       * A class that implements policies for how the
       * DoFHandler::distribute_dofs and DoFHandler::renumber_dofs functions
       * should work.
       */
      template <int dim, int spacedim>
      class PolicyBase
      {
      public:
        /**
         * Destructor.
         */
        virtual ~PolicyBase() = default;

        /**
         * Distribute degrees of freedom on the DoFHandler object associated
         * with this policy object. The argument is a reference to the
         * NumberCache of the DoFHandler object. The function may modify it to
         * make DoFHandler related functions work properly when called within
         * the policies classes. The updated NumberCache is written to that
         * argument.
         */
        virtual NumberCache
        distribute_dofs() const = 0;

        /**
         * Distribute the multigrid dofs on each level of the DoFHandler
         * associated with this policy object. Return a vector of number
         * caches for all of the levels.
         */
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const = 0;

        /**
         * Renumber degrees of freedom as specified by the first argument.
         *
         * Return an updated NumberCache for the DoFHandler after renumbering.
         */
        virtual NumberCache
        renumber_dofs(
          const std::vector<types::global_dof_index> &new_numbers) const = 0;

        /**
         * Renumber multilevel degrees of freedom on one level of a multigrid
         * hierarchy. The second argument specifies the set of new DoF
         * indices.
         *
         * Return an updated NumberCache for the specified level of the
         * DoFHandler after renumbering.
         */
        virtual NumberCache
        renumber_mg_dofs(
          const unsigned int                          level,
          const std::vector<types::global_dof_index> &new_numbers) const = 0;
      };


      /**
       * This class implements the default policy for sequential operations,
       * i.e. for the case where all cells get degrees of freedom.
       */
      template <int dim, int spacedim>
      class Sequential : public PolicyBase<dim, spacedim>
      {
      public:
        /**
         * Constructor.
         * @param dof_handler The DoFHandler object upon which this
         *   policy class is supposed to work.
         */
        Sequential(DoFHandler<dim, spacedim> &dof_handler);

        // documentation is inherited
        virtual NumberCache
        distribute_dofs() const override;

        // documentation is inherited
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const override;

        // documentation is inherited
        virtual NumberCache
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers)
          const override;

        // documentation is inherited
        virtual NumberCache
        renumber_mg_dofs(const unsigned int level,
                         const std::vector<types::global_dof_index>
                           &new_numbers) const override;

      protected:
        /**
         * The DoFHandler object on which this policy object works.
         */
        ObserverPointer<DoFHandler<dim, spacedim>> dof_handler;
      };



      /**
       * This class implements the policy for operations when we use a
       * parallel::shared::Triangulation object.
       */
      template <int dim, int spacedim>
      class ParallelShared : public PolicyBase<dim, spacedim>
      {
      public:
        /**
         * Constructor.
         * @param dof_handler The DoFHandler object upon which this
         *   policy class is supposed to work.
         */
        ParallelShared(DoFHandler<dim, spacedim> &dof_handler);

        /**
         * Distribute degrees of freedom on the object given as first
         * argument.
         *
         * On distribution, DoFs are renumbered subdomain-wise and
         * number_cache.n_locally_owned_dofs_per_processor[i] and
         * number_cache.locally_owned_dofs are updated consistently.
         */
        virtual NumberCache
        distribute_dofs() const override;

        /**
         * This function is not yet implemented.
         */
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const override;

        /**
         * Renumber degrees of freedom as specified by the first argument.
         *
         * The input argument @p new_numbers may either have as many entries
         * as there are global degrees of freedom (i.e. dof_handler.n_dofs() )
         * or dof_handler.locally_owned_dofs().n_elements(). Therefore it can
         * be utilized with renumbering functions implemented for the
         * parallel::distributed case.
         */
        virtual NumberCache
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers)
          const override;

        // documentation is inherited
        virtual NumberCache
        renumber_mg_dofs(const unsigned int level,
                         const std::vector<types::global_dof_index>
                           &new_numbers) const override;

      private:
        /**
         * The DoFHandler object on which this policy object works.
         */
        ObserverPointer<DoFHandler<dim, spacedim>> dof_handler;
      };


      /**
       * This class implements the policy for operations when we use a
       * parallel::DistributedTriangulationBase object.
       */
      template <int dim, int spacedim>
      class ParallelDistributed : public PolicyBase<dim, spacedim>
      {
      public:
        /**
         * Constructor.
         * @param dof_handler The DoFHandler object upon which this
         *   policy class is supposed to work.
         */
        ParallelDistributed(DoFHandler<dim, spacedim> &dof_handler);

        // documentation is inherited
        virtual NumberCache
        distribute_dofs() const override;

        // documentation is inherited
        virtual std::vector<NumberCache>
        distribute_mg_dofs() const override;

        // documentation is inherited
        virtual NumberCache
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers)
          const override;

        // documentation is inherited
        virtual NumberCache
        renumber_mg_dofs(const unsigned int level,
                         const std::vector<types::global_dof_index>
                           &new_numbers) const override;

      private:
        /**
         * The DoFHandler object on which this policy object works.
         */
        ObserverPointer<DoFHandler<dim, spacedim>> dof_handler;
      };
    } // namespace Policy
  }   // namespace DoFHandlerImplementation
} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif
