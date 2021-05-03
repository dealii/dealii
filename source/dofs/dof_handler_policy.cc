// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/geometry_info.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler_policy.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>
#include <memory>
#include <numeric>
#include <set>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DoFHandlerImplementation
  {
    namespace Policy
    {
      // use class dealii::DoFHandler instead
      // of namespace internal::DoFHandler in
      // the following
      using dealii::DoFHandler;

      namespace
      {
        /**
         * We'll use this constant to mark unnumbered degrees of freedom
         * during their enumeration in multiple parts of the code.
         * This constant is necessary to distinguish between valid and
         * invalid degrees of freedom.
         */
        const types::global_dof_index enumeration_dof_index =
          numbers::invalid_dof_index - 1;

        /**
         * Update the cache used for cell dof indices on all (non-artificial)
         * active cells of the given DoFHandler.
         */
        template <int dim, int spacedim>
        void
        update_all_active_cell_dof_indices_caches(
          const DoFHandler<dim, spacedim> &dof_handler)
        {
          const auto worker = [](const auto &cell, void *, void *) {
            if (!cell->is_artificial())
              cell->update_cell_dof_indices_cache();
          };

          // parallelize filling all of the cell caches. by using
          // WorkStream, we make sure that we only run through the
          // range of iterators once, whereas a parallel_for loop
          // for example has to split the range multiple times,
          // which is expensive because cell iterators are not
          // random access iterators with a cheap operator-
          WorkStream::run(dof_handler.begin_active(),
                          dof_handler.end(),
                          worker,
                          /* copier */ std::function<void(void *)>(),
                          /* scratch_data */ nullptr,
                          /* copy_data */ nullptr,
                          2 * MultithreadInfo::n_threads(),
                          /* chunk_size = */ 32);
        }


        /**
         * Update the cache used for cell dof indices on all (non-artificial)
         * level (multigrid) cells of the given DoFHandler.
         */
        template <int dim, int spacedim>
        void
        update_all_level_cell_dof_indices_caches(
          const DoFHandler<dim, spacedim> &dof_handler)
        {
          const auto worker = [](const auto &cell, void *, void *) {
            if (cell->has_children() || !cell->is_artificial())
              cell->update_cell_dof_indices_cache();
          };

          // parallelize filling all of the cell caches. by using
          // WorkStream, we make sure that we only run through the
          // range of iterators once, whereas a parallel_for loop
          // for example has to split the range multiple times,
          // which is expensive because cell iterators are not
          // random access iterators with a cheap operator-
          WorkStream::run(dof_handler.begin(),
                          dof_handler.end(),
                          worker,
                          /* copier */ std::function<void(void *)>(),
                          /* scratch_data */ nullptr,
                          /* copy_data */ nullptr,
                          2 * MultithreadInfo::n_threads(),
                          /* chunk_size = */ 32);
        }


        using DoFIdentities =
          std::vector<std::pair<unsigned int, unsigned int>>;


        /**
         * Make sure that the given @p identities pointer points to a
         * valid array. If the pointer is zero beforehand, create an
         * entry with the correct data. If it is nonzero, don't touch
         * it.
         *
         * @p structdim denotes the dimension of the objects on which
         * identities are to be represented, i.e. zero for vertices,
         * one for lines, etc.
         */
        template <int structdim, int dim, int spacedim>
        const std::unique_ptr<DoFIdentities> &
        ensure_existence_and_return_dof_identities(
          const FiniteElement<dim, spacedim> &fe1,
          const FiniteElement<dim, spacedim> &fe2,
          std::unique_ptr<DoFIdentities> &    identities,
          const unsigned int face_no = numbers::invalid_unsigned_int)
        {
          Assert(structdim == 2 || face_no == numbers::invalid_unsigned_int,
                 ExcInternalError());

          // see if we need to fill this entry, or whether it already
          // exists
          if (identities.get() == nullptr)
            {
              switch (structdim)
                {
                  case 0:
                    {
                      identities = std::make_unique<DoFIdentities>(
                        fe1.hp_vertex_dof_identities(fe2));
                      break;
                    }

                  case 1:
                    {
                      identities = std::make_unique<DoFIdentities>(
                        fe1.hp_line_dof_identities(fe2));
                      break;
                    }

                  case 2:
                    {
                      identities = std::make_unique<DoFIdentities>(
                        fe1.hp_quad_dof_identities(fe2, face_no));
                      break;
                    }

                  default:
                    Assert(false, ExcNotImplemented());
                }

              // double check whether the newly created entries make
              // any sense at all
              for (unsigned int i = 0; i < identities->size(); ++i)
                {
                  Assert((*identities)[i].first <
                           fe1.template n_dofs_per_object<structdim>(face_no),
                         ExcInternalError());
                  Assert((*identities)[i].second <
                           fe2.template n_dofs_per_object<structdim>(face_no),
                         ExcInternalError());
                }
            }

          return identities;
        }
      } // namespace



      struct Implementation
      {
        /* -------------- distribute_dofs functionality ------------- */

        /**
         * Compute identities between DoFs located on vertices. Called from
         * distribute_dofs().
         */
        template <int dim, int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_vertex_dof_identities(
          const DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          std::map<types::global_dof_index, types::global_dof_index>
            dof_identities;

          // Note: we may wish to have something here similar to what
          // we do for lines and quads, namely that we only identify
          // dofs for any FE towards the most dominating one. however,
          // it is not clear whether this is actually necessary for
          // vertices at all, I can't think of a finite element that
          // would make that necessary...
          dealii::Table<2, std::unique_ptr<DoFIdentities>>
            vertex_dof_identities(dof_handler.get_fe_collection().size(),
                                  dof_handler.get_fe_collection().size());

          // loop over all vertices and see which one we need to work on
          for (unsigned int vertex_index = 0;
               vertex_index < dof_handler.get_triangulation().n_vertices();
               ++vertex_index)
            if (dof_handler.get_triangulation()
                  .get_used_vertices()[vertex_index] == true)
              {
                const unsigned int n_active_fe_indices =
                  dealii::internal::DoFAccessorImplementation::Implementation::
                    n_active_fe_indices(dof_handler,
                                        0,
                                        vertex_index,
                                        std::integral_constant<int, 0>());

                if (n_active_fe_indices > 1)
                  {
                    const std::set<unsigned int> fe_indices =
                      dealii::internal::DoFAccessorImplementation::
                        Implementation::get_active_fe_indices(
                          dof_handler,
                          0,
                          vertex_index,
                          std::integral_constant<int, 0>());

                    // find out which is the most dominating finite
                    // element of the ones that are used on this vertex
                    unsigned int most_dominating_fe_index =
                      dof_handler.get_fe_collection().find_dominating_fe(
                        fe_indices,
                        /*codim*/ dim);

                    // if we haven't found a dominating finite element,
                    // choose the very first one to be dominant
                    if (most_dominating_fe_index ==
                        numbers::invalid_unsigned_int)
                      most_dominating_fe_index =
                        dealii::internal::DoFAccessorImplementation::
                          Implementation::nth_active_fe_index(
                            dof_handler,
                            0,
                            vertex_index,
                            0,
                            std::integral_constant<int, 0>());

                    // loop over the indices of all the finite
                    // elements that are not dominating, and
                    // identify their dofs to the most dominating
                    // one
                    for (const auto &other_fe_index : fe_indices)
                      if (other_fe_index != most_dominating_fe_index)
                        {
                          // make sure the entry in the equivalence
                          // table exists
                          const auto &identities =
                            *ensure_existence_and_return_dof_identities<0>(
                              dof_handler.get_fe(most_dominating_fe_index),
                              dof_handler.get_fe(other_fe_index),
                              vertex_dof_identities[most_dominating_fe_index]
                                                   [other_fe_index]);

                          // then loop through the identities we
                          // have. first get the global numbers of the
                          // dofs we want to identify and make sure they
                          // are not yet constrained to anything else,
                          // except for to each other. use the rule that
                          // we will always constrain the dof with the
                          // higher FE index to the one with the lower,
                          // to avoid circular reasoning.
                          for (const auto &identity : identities)
                            {
                              const types::global_dof_index primary_dof_index =
                                dealii::internal::DoFAccessorImplementation::
                                  Implementation::get_dof_index(
                                    dof_handler,
                                    0,
                                    vertex_index,
                                    most_dominating_fe_index,
                                    identity.first,
                                    std::integral_constant<int, 0>());
                              const types::global_dof_index
                                dependent_dof_index =
                                  dealii::internal::DoFAccessorImplementation::
                                    Implementation::get_dof_index(
                                      dof_handler,
                                      0,
                                      vertex_index,
                                      other_fe_index,
                                      identity.second,
                                      std::integral_constant<int, 0>());

                              // on subdomain boundaries, we will
                              // encounter invalid DoFs on ghost cells,
                              // for which we have not yet distributed
                              // valid indices. depending on which finte
                              // element is dominating the other on this
                              // interface, we either have to constrain
                              // the valid to the invalid indices, or vice
                              // versa.
                              //
                              // we only store an identity if we are about
                              // to overwrite a valid DoF. we will skip
                              // constraining invalid DoFs for now, and
                              // consider them later in Phase 5.
                              if (dependent_dof_index !=
                                  numbers::invalid_dof_index)
                                {
                                  // if the DoF indices of both elements
                                  // are already distributed, i.e., both
                                  // of these 'fe_indices' are associated
                                  // with a locally owned cell, then we
                                  // should either not have a dof_identity
                                  // yet, or it must come out here to be
                                  // exactly as we had computed before
                                  if (primary_dof_index !=
                                      numbers::invalid_dof_index)
                                    Assert(
                                      (dof_identities.find(primary_dof_index) ==
                                       dof_identities.end()) ||
                                        (dof_identities[dependent_dof_index] ==
                                         primary_dof_index),
                                      ExcInternalError());

                                  dof_identities[dependent_dof_index] =
                                    primary_dof_index;
                                }
                            }
                        }
                  }
              }

          return dof_identities;
        }


        /**
         * Compute identities between DoFs located on lines. Called from
         * distribute_dofs().
         */
        template <int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_line_dof_identities(const DoFHandler<1, spacedim> &dof_handler)
        {
          (void)dof_handler;
          Assert(dof_handler.hp_capability_enabled == true,
                 (typename DoFHandler<1, spacedim>::ExcOnlyAvailableWithHP()));

          return std::map<types::global_dof_index, types::global_dof_index>();
        }


        template <int dim, int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_line_dof_identities(
          const DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          std::map<types::global_dof_index, types::global_dof_index>
            dof_identities;

          // we will mark lines that we have already treated, so first save and
          // clear the user flags on lines and later restore them
          std::vector<bool> user_flags;
          dof_handler.get_triangulation().save_user_flags_line(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .clear_user_flags_line();

          // An implementation of the algorithm described in the hp-paper,
          // including the modification mentioned later in the "complications in
          // 3-d" subsections
          //
          // as explained there, we do something only if there are exactly 2
          // finite elements associated with an object. if there is only one,
          // then there is nothing to do anyway, and if there are 3 or more,
          // then we can get into trouble. note that this only happens for lines
          // in 3d and higher, and for quads only in 4d and higher, so this
          // isn't a particularly frequent case
          //
          // there is one case, however, that we would like to handle (see, for
          // example, the hp/crash_15 testcase): if we have
          // FESystem(FE_Q(2),FE_DGQ(i)) elements for a bunch of values 'i',
          // then we should be able to handle this because we can simply unify
          // *all* dofs, not only a some. so what we do is to first treat all
          // pairs of finite elements that have *identical* dofs, and then only
          // deal with those that are not identical of which we can handle at
          // most 2
          dealii::Table<2, std::unique_ptr<DoFIdentities>> line_dof_identities(
            dof_handler.fe_collection.size(), dof_handler.fe_collection.size());

          for (const auto &cell : dof_handler.active_cell_iterators())
            for (const auto l : cell->line_indices())
              if (cell->line(l)->user_flag_set() == false)
                {
                  const auto line = cell->line(l);
                  line->set_user_flag();

                  unsigned int unique_sets_of_dofs =
                    line->n_active_fe_indices();

                  // do a first loop over all sets of dofs and do identity
                  // uniquification
                  const unsigned int n_active_fe_indices =
                    line->n_active_fe_indices();
                  for (unsigned int f = 0; f < n_active_fe_indices; ++f)
                    for (unsigned int g = f + 1; g < n_active_fe_indices; ++g)
                      {
                        const unsigned int fe_index_1 =
                                             line->nth_active_fe_index(f),
                                           fe_index_2 =
                                             line->nth_active_fe_index(g);

                        // as described in the hp-paper, we only unify on lines
                        // when there are at most two different FE objects
                        // assigned on it.
                        // however, more than two 'active_fe_indices' can be
                        // attached that still fulfill the above criterion,
                        // i.e. when two different FiniteElement objects are
                        // assigned to neighboring cells that map their degrees
                        // of freedom one-to-one.
                        // we cannot verify with certainty if two dofs each of
                        // separate FiniteElement objects actually map
                        // one-to-one. however, checking for the number of
                        // 'dofs_per_line' turned out to be a reasonable
                        // approach, that also works for e.g. two different
                        // FE_Q objects of the same order, from which one is
                        // enhanced by a bubble function that is zero on the
                        // boundary.
                        if ((dof_handler.get_fe(fe_index_1).n_dofs_per_line() ==
                             dof_handler.get_fe(fe_index_2)
                               .n_dofs_per_line()) &&
                            (dof_handler.get_fe(fe_index_1).n_dofs_per_line() >
                             0))
                          {
                            // the number of dofs per line is identical
                            const unsigned int dofs_per_line =
                              dof_handler.get_fe(fe_index_1).n_dofs_per_line();

                            const auto &identities =
                              *ensure_existence_and_return_dof_identities<1>(
                                dof_handler.get_fe(fe_index_1),
                                dof_handler.get_fe(fe_index_2),
                                line_dof_identities[fe_index_1][fe_index_2]);
                            // see if these sets of dofs are identical. the
                            // first condition for this is that indeed there are
                            // n identities
                            if (identities.size() == dofs_per_line)
                              {
                                unsigned int i = 0;
                                for (; i < dofs_per_line; ++i)
                                  if ((identities[i].first != i) &&
                                      (identities[i].second != i))
                                    // not an identity
                                    break;

                                if (i == dofs_per_line)
                                  {
                                    // The line dofs (i.e., the ones interior to
                                    // a line) of these two finite elements are
                                    // identical. Note that there could be
                                    // situations when one element still
                                    // dominates another, e.g.: FE_Q(2) x
                                    // FE_Nothing(dominate) vs FE_Q(2) x FE_Q(1)

                                    --unique_sets_of_dofs;

                                    // determine which one of both finite
                                    // elements is the dominating one.
                                    const std::set<unsigned int> fe_indices{
                                      fe_index_1, fe_index_2};

                                    unsigned int dominating_fe_index =
                                      dof_handler.get_fe_collection()
                                        .find_dominating_fe(fe_indices,
                                                            /*codim=*/dim - 1);
                                    unsigned int other_fe_index =
                                      numbers::invalid_unsigned_int;

                                    if (dominating_fe_index !=
                                        numbers::invalid_unsigned_int)
                                      other_fe_index =
                                        (dominating_fe_index == fe_index_1) ?
                                          fe_index_2 :
                                          fe_index_1;
                                    else
                                      {
                                        // if we haven't found a dominating
                                        // finite element, choose the one with
                                        // the lower index to be dominating
                                        dominating_fe_index = fe_index_1;
                                        other_fe_index      = fe_index_2;
                                      }

                                    for (unsigned int j = 0; j < dofs_per_line;
                                         ++j)
                                      {
                                        const types::global_dof_index
                                          primary_dof_index = line->dof_index(
                                            j, dominating_fe_index);
                                        const types::global_dof_index
                                          dependent_dof_index =
                                            line->dof_index(j, other_fe_index);

                                        // on subdomain boundaries, we will
                                        // encounter invalid DoFs on ghost
                                        // cells, for which we have not yet
                                        // distributed valid indices. depending
                                        // on which finte element is dominating
                                        // the other on this interface, we
                                        // either have to constrain the valid to
                                        // the invalid indices, or vice versa.
                                        //
                                        // we only store an identity if we are
                                        // about to overwrite a valid DoF. we
                                        // will skip constraining invalid DoFs
                                        // for now, and consider them later in
                                        // Phase 5.
                                        if (dependent_dof_index !=
                                            numbers::invalid_dof_index)
                                          {
                                            if (primary_dof_index !=
                                                numbers::invalid_dof_index)
                                              {
                                                // if primary dof was already
                                                // constrained, constrain to
                                                // that one, otherwise constrain
                                                // dependent to primary
                                                if (dof_identities.find(
                                                      primary_dof_index) !=
                                                    dof_identities.end())
                                                  {
                                                    // if the DoF indices of
                                                    // both elements are already
                                                    // distributed, i.e., both
                                                    // of these 'fe_indices' are
                                                    // associated with a locally
                                                    // owned cell, then we
                                                    // should either not have a
                                                    // dof_identity yet, or it
                                                    // must come out here to be
                                                    // exactly as we had
                                                    // computed before
                                                    Assert(
                                                      dof_identities.find(
                                                        dof_identities
                                                          [primary_dof_index]) ==
                                                        dof_identities.end(),
                                                      ExcInternalError());

                                                    dof_identities
                                                      [dependent_dof_index] =
                                                        dof_identities
                                                          [primary_dof_index];
                                                  }
                                                else
                                                  {
                                                    // see comment above for an
                                                    // explanation of this
                                                    // assertion
                                                    Assert(
                                                      (dof_identities.find(
                                                         primary_dof_index) ==
                                                       dof_identities.end()) ||
                                                        (dof_identities
                                                           [dependent_dof_index] ==
                                                         primary_dof_index),
                                                      ExcInternalError());

                                                    dof_identities
                                                      [dependent_dof_index] =
                                                        primary_dof_index;
                                                  }
                                              }
                                            else
                                              {
                                                // set dependent_dof to
                                                // primary_dof_index, which is
                                                // invalid
                                                dof_identities
                                                  [dependent_dof_index] =
                                                    numbers::invalid_dof_index;
                                              }
                                          }
                                      }
                                  }
                              }
                          }
                      }

                  // if at this point, there is only one unique set of dofs
                  // left, then we have taken care of everything above. if there
                  // are two, then we need to deal with them here. if there are
                  // more, then we punt, as described in the paper (and
                  // mentioned above)
                  // TODO: The check for 'dim==2' was inserted by intuition. It
                  // fixes
                  // the previous problems with step-27 in 3D. But an
                  // explanation for this is still required, and what we do here
                  // is not what we describe in the paper!.
                  if ((unique_sets_of_dofs == 2) && (dim == 2))
                    {
                      const std::set<unsigned int> fe_indices =
                        line->get_active_fe_indices();

                      // find out which is the most dominating finite element of
                      // the ones that are used on this line
                      const unsigned int most_dominating_fe_index =
                        dof_handler.get_fe_collection().find_dominating_fe(
                          fe_indices,
                          /*codim=*/dim - 1);

                      // if we found the most dominating element, then use this
                      // to eliminate some of the degrees of freedom by
                      // identification. otherwise, the code that computes
                      // hanging node constraints will have to deal with it by
                      // computing appropriate constraints along this face/edge
                      if (most_dominating_fe_index !=
                          numbers::invalid_unsigned_int)
                        {
                          // loop over the indices of all the finite elements
                          // that are not dominating, and identify their dofs to
                          // the most dominating one
                          for (const auto &other_fe_index : fe_indices)
                            if (other_fe_index != most_dominating_fe_index)
                              {
                                const auto &identities =
                                  *ensure_existence_and_return_dof_identities<
                                    1>(dof_handler.get_fe(
                                         most_dominating_fe_index),
                                       dof_handler.get_fe(other_fe_index),
                                       line_dof_identities
                                         [most_dominating_fe_index]
                                         [other_fe_index]);

                                for (const auto &identity : identities)
                                  {
                                    const types::global_dof_index
                                      primary_dof_index = line->dof_index(
                                        identity.first,
                                        most_dominating_fe_index);
                                    const types::global_dof_index
                                      dependent_dof_index =
                                        line->dof_index(identity.second,
                                                        other_fe_index);

                                    // on subdomain boundaries, we will
                                    // encounter invalid DoFs on ghost cells,
                                    // for which we have not yet distributed
                                    // valid indices. depending on which finte
                                    // element is dominating the other on this
                                    // interface, we either have to constrain
                                    // the valid to the invalid indices, or vice
                                    // versa.
                                    //
                                    // we only store an identity if we are about
                                    // to overwrite a valid DoF. we will skip
                                    // constraining invalid DoFs for now, and
                                    // consider them later in Phase 5.
                                    if (dependent_dof_index !=
                                        numbers::invalid_dof_index)
                                      {
                                        // if the DoF indices of both elements
                                        // are already distributed, i.e., both
                                        // of these 'fe_indices' are associated
                                        // with a locally owned cell, then we
                                        // should either not have a dof_identity
                                        // yet, or it must come out here to be
                                        // exactly as we had computed before
                                        if (primary_dof_index !=
                                            numbers::invalid_dof_index)
                                          Assert((dof_identities.find(
                                                    primary_dof_index) ==
                                                  dof_identities.end()) ||
                                                   (dof_identities
                                                      [dependent_dof_index] ==
                                                    primary_dof_index),
                                                 ExcInternalError());

                                        dof_identities[dependent_dof_index] =
                                          primary_dof_index;
                                      }
                                  }
                              }
                        }
                    }
                }

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .load_user_flags_line(user_flags);

          return dof_identities;
        }



        /**
         * Compute identities between DoFs located on quads. Called from
         * distribute_dofs().
         */
        template <int dim, int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_quad_dof_identities(
          const DoFHandler<dim, spacedim> &dof_handler)
        {
          (void)dof_handler;
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          // this function should only be called for dim<3 where there are
          // no quad dof identies. for dim==3, the specialization below should
          // take care of it
          Assert(dim < 3, ExcInternalError());

          return std::map<types::global_dof_index, types::global_dof_index>();
        }


        template <int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_quad_dof_identities(const DoFHandler<3, spacedim> &dof_handler)
        {
          Assert(dof_handler.hp_capability_enabled == true,
                 (typename DoFHandler<3, spacedim>::ExcOnlyAvailableWithHP()));

          const int dim = 3;

          std::map<types::global_dof_index, types::global_dof_index>
            dof_identities;


          // we will mark quads that we have already treated, so first
          // save and clear the user flags on quads and later restore
          // them
          std::vector<bool> user_flags;
          dof_handler.get_triangulation().save_user_flags_quad(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .clear_user_flags_quad();

          // An implementation of the algorithm described in the hp-
          // paper, including the modification mentioned later in the
          // "complications in 3-d" subsections
          //
          // as explained there, we do something only if there are
          // exactly 2 finite elements associated with an object. if
          // there is only one, then there is nothing to do anyway,
          // and if there are 3 or more, then we can get into
          // trouble. note that this only happens for lines in 3d and
          // higher, and for quads only in 4d and higher, so this
          // isn't a particularly frequent case
          dealii::Table<3, std::unique_ptr<DoFIdentities>> quad_dof_identities(
            dof_handler.fe_collection.size(),
            dof_handler.fe_collection.size(),
            2 /*triangle (0) or quadrilateral (1)*/);

          for (const auto &cell : dof_handler.active_cell_iterators())
            for (const auto q : cell->face_indices())
              if ((cell->quad(q)->user_flag_set() == false) &&
                  (cell->quad(q)->n_active_fe_indices() == 2))
                {
                  const auto quad = cell->quad(q);
                  quad->set_user_flag();

                  const std::set<unsigned int> fe_indices =
                    quad->get_active_fe_indices();

                  // find out which is the most dominating finite
                  // element of the ones that are used on this quad
                  const unsigned int most_dominating_fe_index =
                    dof_handler.get_fe_collection().find_dominating_fe(
                      fe_indices,
                      /*codim=*/dim - 2);

                  const unsigned int most_dominating_fe_index_face_no =
                    cell->active_fe_index() == most_dominating_fe_index ?
                      q :
                      cell->neighbor_face_no(q);

                  // if we found the most dominating element, then use
                  // this to eliminate some of the degrees of freedom
                  // by identification. otherwise, the code that
                  // computes hanging node constraints will have to
                  // deal with it by computing appropriate constraints
                  // along this face/edge
                  if (most_dominating_fe_index != numbers::invalid_unsigned_int)
                    {
                      // loop over the indices of all the finite
                      // elements that are not dominating, and
                      // identify their dofs to the most dominating
                      // one
                      for (const auto &other_fe_index : fe_indices)
                        if (other_fe_index != most_dominating_fe_index)
                          {
                            const auto &identities =
                              *ensure_existence_and_return_dof_identities<2>(
                                dof_handler.get_fe(most_dominating_fe_index),
                                dof_handler.get_fe(other_fe_index),
                                quad_dof_identities
                                  [most_dominating_fe_index][other_fe_index]
                                  [cell->quad(q)->reference_cell() ==
                                   dealii::ReferenceCells::Quadrilateral],
                                most_dominating_fe_index_face_no);

                            for (const auto &identity : identities)
                              {
                                const types::global_dof_index
                                  primary_dof_index =
                                    quad->dof_index(identity.first,
                                                    most_dominating_fe_index);
                                const types::global_dof_index
                                  dependent_dof_index =
                                    quad->dof_index(identity.second,
                                                    other_fe_index);

                                // we only store an identity if we are about to
                                // overwrite a valid degree of freedom. we will
                                // skip invalid degrees of freedom (that are
                                // associated with ghost cells) for now, and
                                // consider them later in phase 5.
                                if (dependent_dof_index !=
                                    numbers::invalid_dof_index)
                                  {
                                    // if the DoF indices of both elements are
                                    // already distributed, i.e., both of these
                                    // 'fe_indices' are associated with a
                                    // locally owned cell, then we should either
                                    // not have a dof_identity yet, or it must
                                    // come out here to be exactly as we had
                                    // computed before
                                    if (primary_dof_index !=
                                        numbers::invalid_dof_index)
                                      Assert((dof_identities.find(
                                                primary_dof_index) ==
                                              dof_identities.end()) ||
                                               (dof_identities
                                                  [dependent_dof_index] ==
                                                primary_dof_index),
                                             ExcInternalError());

                                    dof_identities[dependent_dof_index] =
                                      primary_dof_index;
                                  }
                              }
                          }
                    }
                }

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .load_user_flags_quad(user_flags);

          return dof_identities;
        }



        /**
         * Compute the constraints that correspond to unifying DoF indices
         * on vertices, lines, and quads. Do so in parallel.
         */
        template <int dim, int spacedim>
        static void
        compute_dof_identities(std::vector<std::map<types::global_dof_index,
                                                    types::global_dof_index>>
                                 &all_constrained_indices,
                               const DoFHandler<dim, spacedim> &dof_handler)
        {
          if (dof_handler.hp_capability_enabled == false)
            return;

          Assert(all_constrained_indices.size() == dim, ExcInternalError());

          Threads::TaskGroup<> tasks;

          unsigned int i = 0;
          tasks += Threads::new_task([&, i]() {
            all_constrained_indices[i] =
              compute_vertex_dof_identities(dof_handler);
          });

          if (dim > 1)
            {
              ++i;
              tasks += Threads::new_task([&, i]() {
                all_constrained_indices[i] =
                  compute_line_dof_identities(dof_handler);
              });
            }

          if (dim > 2)
            {
              ++i;
              tasks += Threads::new_task([&, i]() {
                all_constrained_indices[i] =
                  compute_quad_dof_identities(dof_handler);
              });
            }

          tasks.join_all();
        }



        /**
         * Once degrees of freedom have been distributed on all cells, we may
         * want to eliminate duplicates, and enumerate the remaining ones
         * consecutively. This particular function is responsible for the
         * latter part.
         *
         * This function stores the new indices of all DoFs in ascending order
         * in @p new_dof_indices, which can be used to renumber all DoFs with
         * the renumber_dofs() function later.
         *
         * This vector will contain enumerated indices, skipping invalid indices
         * previously stored in it. Additionally, if a
         * @p all_constrained_indices parameter is provided, DoF identity
         * relations  will be considered as well during the enumeration process
         * by identifying similar DoFs on vertices, lines and quads.
         *
         * Returns the final number of degrees of freedom, which is the number
         * of all valid DoF indices in @p new_dof_indices.
         */
        template <int dim, int spacedim>
        static types::global_dof_index
        enumerate_dof_indices_for_renumbering(
          std::vector<types::global_dof_index> &new_dof_indices,
          const std::vector<
            std::map<types::global_dof_index, types::global_dof_index>>
            &all_constrained_indices,
          const DoFHandler<dim, spacedim> &)
        {
          Assert(all_constrained_indices.size() == dim, ExcInternalError());

          // first preset the new DoF indices that are identities
          for (const auto &constrained_dof_indices : all_constrained_indices)
            for (const auto &p : constrained_dof_indices)
              if (new_dof_indices[p.first] != numbers::invalid_dof_index)
                {
                  Assert(new_dof_indices[p.first] == enumeration_dof_index,
                         ExcInternalError());

                  new_dof_indices[p.first] = p.second;
                }

          // then enumerate the rest
          types::global_dof_index next_free_dof = 0;
          for (auto &new_dof_index : new_dof_indices)
            if (new_dof_index == enumeration_dof_index)
              new_dof_index = next_free_dof++;

          // then loop over all those that are constrained and record the
          // new dof number for those
          for (const auto &constrained_dof_indices : all_constrained_indices)
            for (const auto &p : constrained_dof_indices)
              if (new_dof_indices[p.first] != numbers::invalid_dof_index)
                {
                  Assert(new_dof_indices[p.first] != enumeration_dof_index,
                         ExcInternalError());

                  if (p.second != numbers::invalid_dof_index)
                    new_dof_indices[p.first] = new_dof_indices[p.second];
                }

          for (const types::global_dof_index new_dof_index : new_dof_indices)
            {
              (void)new_dof_index;
              Assert(new_dof_index != enumeration_dof_index,
                     ExcInternalError());
              Assert(new_dof_index < next_free_dof ||
                       new_dof_index == numbers::invalid_dof_index,
                     ExcInternalError());
            }

          return next_free_dof;
        }



        /**
         * Once degrees of freedom have been distributed on all cells, see if
         * we can identify DoFs on neighboring cells. This function does
         * nothing unless the DoFHandler has hp-capabilities.
         *
         * Return the final number of degrees of freedom, which is the old one
         * minus however many were identified.
         */
        template <int dim, int spacedim>
        static types::global_dof_index
        unify_dof_indices(const DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int n_dofs_before_identification,
                          const bool         check_validity)
        {
          if (dof_handler.hp_capability_enabled == false)
            return n_dofs_before_identification;

          std::vector<
            std::map<types::global_dof_index, types::global_dof_index>>
            all_constrained_indices(dim);
          compute_dof_identities(all_constrained_indices, dof_handler);

          std::vector<dealii::types::global_dof_index> renumbering(
            n_dofs_before_identification, enumeration_dof_index);
          const types::global_dof_index n_dofs =
            enumerate_dof_indices_for_renumbering(renumbering,
                                                  all_constrained_indices,
                                                  dof_handler);

          renumber_dofs(renumbering, IndexSet(0), dof_handler, check_validity);

          update_all_active_cell_dof_indices_caches(dof_handler);

          return n_dofs;
        }



        /**
         * Merge invalid DoF indices on vertices located on ghost interfaces
         * by a dominating valid one.
         */
        template <int dim, int spacedim>
        static void
        merge_invalid_vertex_dofs_on_ghost_interfaces(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          // Note: we may wish to have something here similar to what
          // we do for lines and quads, namely that we only identify
          // dofs for any FE towards the most dominating one. however,
          // it is not clear whether this is actually necessary for
          // vertices at all, I can't think of a finite element that
          // would make that necessary...
          dealii::Table<2, std::unique_ptr<DoFIdentities>>
            vertex_dof_identities(dof_handler.get_fe_collection().size(),
                                  dof_handler.get_fe_collection().size());

          // mark all vertices on ghost cells
          std::vector<bool> include_vertex(
            dof_handler.get_triangulation().n_vertices(), false);
          if (dynamic_cast<const dealii::parallel::
                             DistributedTriangulationBase<dim, spacedim> *>(
                &dof_handler.get_triangulation()) != nullptr)
            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->is_ghost())
                for (const unsigned int v : cell->vertex_indices())
                  include_vertex[cell->vertex_index(v)] = true;

          // loop over all vertices and see which one we need to work on
          for (unsigned int vertex_index = 0;
               vertex_index < dof_handler.get_triangulation().n_vertices();
               ++vertex_index)
            if ((dof_handler.get_triangulation()
                   .get_used_vertices()[vertex_index] == true) &&
                (include_vertex[vertex_index] == true))
              {
                const unsigned int n_active_fe_indices =
                  dealii::internal::DoFAccessorImplementation::Implementation::
                    n_active_fe_indices(dof_handler,
                                        0,
                                        vertex_index,
                                        std::integral_constant<int, 0>());

                if (n_active_fe_indices > 1)
                  {
                    const std::set<unsigned int> fe_indices =
                      dealii::internal::DoFAccessorImplementation::
                        Implementation::get_active_fe_indices(
                          dof_handler,
                          0,
                          vertex_index,
                          std::integral_constant<int, 0>());

                    // find out which is the most dominating finite
                    // element of the ones that are used on this vertex
                    unsigned int most_dominating_fe_index =
                      dof_handler.get_fe_collection().find_dominating_fe(
                        fe_indices,
                        /*codim=*/dim);

                    // if we haven't found a dominating finite element,
                    // choose the very first one to be dominant similar
                    // to compute_vertex_dof_identities()
                    if (most_dominating_fe_index ==
                        numbers::invalid_unsigned_int)
                      most_dominating_fe_index =
                        dealii::internal::DoFAccessorImplementation::
                          Implementation::nth_active_fe_index(
                            dof_handler,
                            0,
                            vertex_index,
                            0,
                            std::integral_constant<int, 0>());

                    // loop over the indices of all the finite
                    // elements that are not dominating, and
                    // identify their dofs to the most dominating
                    // one
                    for (const auto &other_fe_index : fe_indices)
                      if (other_fe_index != most_dominating_fe_index)
                        {
                          // make sure the entry in the equivalence
                          // table exists
                          const auto &identities =
                            *ensure_existence_and_return_dof_identities<0>(
                              dof_handler.get_fe(most_dominating_fe_index),
                              dof_handler.get_fe(other_fe_index),
                              vertex_dof_identities[most_dominating_fe_index]
                                                   [other_fe_index]);

                          // then loop through the identities we
                          // have. first get the global numbers of the
                          // dofs we want to identify and make sure they
                          // are not yet constrained to anything else,
                          // except for to each other. use the rule that
                          // we will always constrain the dof with the
                          // higher FE index to the one with the lower,
                          // to avoid circular reasoning.
                          for (const auto &identity : identities)
                            {
                              const types::global_dof_index primary_dof_index =
                                dealii::internal::DoFAccessorImplementation::
                                  Implementation::get_dof_index(
                                    dof_handler,
                                    0,
                                    vertex_index,
                                    most_dominating_fe_index,
                                    identity.first,
                                    std::integral_constant<int, 0>());
                              const types::global_dof_index
                                dependent_dof_index =
                                  dealii::internal::DoFAccessorImplementation::
                                    Implementation::get_dof_index(
                                      dof_handler,
                                      0,
                                      vertex_index,
                                      other_fe_index,
                                      identity.second,
                                      std::integral_constant<int, 0>());

                              // check if we are on an interface between
                              // a locally owned and a ghost cell on which
                              // we need to work on.
                              //
                              // all degrees of freedom belonging to
                              // dominating FE indices or to a processor
                              // with a higher rank have been set at this
                              // point (either in Phase 2, or after the
                              // first ghost exchange in Phase 5). thus,
                              // we only have to set the indices of
                              // degrees of freedom that have been
                              // previously flagged invalid.
                              if ((dependent_dof_index ==
                                   numbers::invalid_dof_index) &&
                                  (primary_dof_index !=
                                   numbers::invalid_dof_index))
                                dealii::internal::DoFAccessorImplementation::
                                  Implementation::set_dof_index(
                                    dof_handler,
                                    0,
                                    vertex_index,
                                    other_fe_index,
                                    identity.second,
                                    std::integral_constant<int, 0>(),
                                    primary_dof_index);
                            }
                        }
                  }
              }
        }



        /**
         * Merge invalid DoF indices on lines located on ghost interfaces
         * by a dominating valid one.
         */
        template <int spacedim>
        static void merge_invalid_line_dofs_on_ghost_interfaces(
          DoFHandler<1, spacedim> &dof_handler)
        {
          (void)dof_handler;
          Assert(dof_handler.hp_capability_enabled == true,
                 (typename DoFHandler<1, spacedim>::ExcOnlyAvailableWithHP()));
        }


        template <int dim, int spacedim>
        static void
        merge_invalid_line_dofs_on_ghost_interfaces(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          // we will mark lines that we have already treated, so first save and
          // clear the user flags on lines and later restore them
          std::vector<bool> user_flags;
          dof_handler.get_triangulation().save_user_flags_line(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .clear_user_flags_line();

          // mark all lines on ghost cells
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_ghost())
              for (const auto l : cell->line_indices())
                cell->line(l)->set_user_flag();

          // An implementation of the algorithm described in the hp-paper,
          // including the modification mentioned later in the "complications in
          // 3-d" subsections
          //
          // as explained there, we do something only if there are exactly 2
          // finite elements associated with an object. if there is only one,
          // then there is nothing to do anyway, and if there are 3 or more,
          // then we can get into trouble. note that this only happens for lines
          // in 3d and higher, and for quads only in 4d and higher, so this
          // isn't a particularly frequent case
          //
          // there is one case, however, that we would like to handle (see, for
          // example, the hp/crash_15 testcase): if we have
          // FESystem(FE_Q(2),FE_DGQ(i)) elements for a bunch of values 'i',
          // then we should be able to handle this because we can simply unify
          // *all* dofs, not only a some. so what we do is to first treat all
          // pairs of finite elements that have *identical* dofs, and then only
          // deal with those that are not identical of which we can handle at
          // most 2
          dealii::Table<2, std::unique_ptr<DoFIdentities>> line_dof_identities(
            dof_handler.fe_collection.size(), dof_handler.fe_collection.size());

          for (const auto &cell : dof_handler.active_cell_iterators())
            for (const auto l : cell->line_indices())
              if ((cell->is_locally_owned()) &&
                  (cell->line(l)->user_flag_set() == true))
                {
                  const auto line = cell->line(l);
                  line->clear_user_flag();

                  unsigned int unique_sets_of_dofs =
                    line->n_active_fe_indices();

                  // do a first loop over all sets of dofs and do identity
                  // uniquification
                  const unsigned int n_active_fe_indices =
                    line->n_active_fe_indices();
                  for (unsigned int f = 0; f < n_active_fe_indices; ++f)
                    for (unsigned int g = f + 1; g < n_active_fe_indices; ++g)
                      {
                        const unsigned int fe_index_1 =
                                             line->nth_active_fe_index(f),
                                           fe_index_2 =
                                             line->nth_active_fe_index(g);

                        if ((dof_handler.get_fe(fe_index_1).n_dofs_per_line() ==
                             dof_handler.get_fe(fe_index_2)
                               .n_dofs_per_line()) &&
                            (dof_handler.get_fe(fe_index_1).n_dofs_per_line() >
                             0))
                          {
                            // the number of dofs per line is identical
                            const unsigned int dofs_per_line =
                              dof_handler.get_fe(fe_index_1).n_dofs_per_line();

                            const auto &identities =
                              *ensure_existence_and_return_dof_identities<1>(
                                dof_handler.get_fe(fe_index_1),
                                dof_handler.get_fe(fe_index_2),
                                line_dof_identities[fe_index_1][fe_index_2]);
                            // see if these sets of dofs are identical. the
                            // first condition for this is that indeed there are
                            // n identities
                            if (identities.size() == dofs_per_line)
                              {
                                unsigned int i = 0;
                                for (; i < dofs_per_line; ++i)
                                  if ((identities[i].first != i) &&
                                      (identities[i].second != i))
                                    // not an identity
                                    break;

                                if (i == dofs_per_line)
                                  {
                                    // The line dofs (i.e., the ones interior to
                                    // a line) of these two finite elements are
                                    // identical. Note that there could be
                                    // situations when one element still
                                    // dominates another, e.g.: FE_Q(2) x
                                    // FE_Nothing(dominate) vs FE_Q(2) x FE_Q(1)

                                    --unique_sets_of_dofs;

                                    // determine which one of both finite
                                    // elements is the dominating one.
                                    const std::set<unsigned int> fe_indices{
                                      fe_index_1, fe_index_2};

                                    unsigned int dominating_fe_index =
                                      dof_handler.get_fe_collection()
                                        .find_dominating_fe(fe_indices,
                                                            /*codim*/ dim - 1);
                                    unsigned int other_fe_index =
                                      numbers::invalid_unsigned_int;

                                    if (dominating_fe_index !=
                                        numbers::invalid_unsigned_int)
                                      other_fe_index =
                                        (dominating_fe_index == fe_index_1) ?
                                          fe_index_2 :
                                          fe_index_1;
                                    else
                                      {
                                        // if we haven't found a dominating
                                        // finite element, choose the one with
                                        // the lower index to be dominating
                                        dominating_fe_index = fe_index_1;
                                        other_fe_index      = fe_index_2;
                                      }

                                    for (unsigned int j = 0; j < dofs_per_line;
                                         ++j)
                                      {
                                        const types::global_dof_index
                                          primary_dof_index = line->dof_index(
                                            j, dominating_fe_index);
                                        const types::global_dof_index
                                          dependent_dof_index =
                                            line->dof_index(j, other_fe_index);

                                        // check if we are on an interface
                                        // between a locally owned and a ghost
                                        // cell on which we need to work on.
                                        //
                                        // all degrees of freedom belonging to
                                        // dominating fe_indices or to a
                                        // processor with a higher rank have
                                        // been set at this point (either in
                                        // Phase 2, or after the first ghost
                                        // exchange in Phase 5). thus, we only
                                        // have to set the indices of degrees
                                        // of freedom that have been previously
                                        // flagged invalid.
                                        if ((dependent_dof_index ==
                                             numbers::invalid_dof_index) &&
                                            (primary_dof_index !=
                                             numbers::invalid_dof_index))
                                          line->set_dof_index(j,
                                                              primary_dof_index,
                                                              fe_index_2);
                                      }
                                  }
                              }
                          }
                      }

                  // if at this point, there is only one unique set of dofs
                  // left, then we have taken care of everything above. if there
                  // are two, then we need to deal with them here. if there are
                  // more, then we punt, as described in the paper (and
                  // mentioned above)
                  // TODO: The check for 'dim==2' was inserted by intuition. It
                  // fixes
                  // the previous problems with step-27 in 3D. But an
                  // explanation for this is still required, and what we do here
                  // is not what we describe in the paper!.
                  if ((unique_sets_of_dofs == 2) && (dim == 2))
                    {
                      const std::set<unsigned int> fe_indices =
                        line->get_active_fe_indices();

                      // find out which is the most dominating finite element of
                      // the ones that are used on this line
                      const unsigned int most_dominating_fe_index =
                        dof_handler.get_fe_collection().find_dominating_fe(
                          fe_indices,
                          /*codim=*/dim - 1);

                      // if we found the most dominating element, then use this
                      // to eliminate some of the degrees of freedom by
                      // identification. otherwise, the code that computes
                      // hanging node constraints will have to deal with it by
                      // computing appropriate constraints along this face/edge
                      if (most_dominating_fe_index !=
                          numbers::invalid_unsigned_int)
                        {
                          // loop over the indices of all the finite elements
                          // that are not dominating, and identify their dofs to
                          // the most dominating one
                          for (const auto &other_fe_index : fe_indices)
                            if (other_fe_index != most_dominating_fe_index)
                              {
                                const auto &identities =
                                  *ensure_existence_and_return_dof_identities<
                                    1>(dof_handler.get_fe(
                                         most_dominating_fe_index),
                                       dof_handler.get_fe(other_fe_index),
                                       line_dof_identities
                                         [most_dominating_fe_index]
                                         [other_fe_index]);

                                for (const auto &identity : identities)
                                  {
                                    const types::global_dof_index
                                      primary_dof_index = line->dof_index(
                                        identity.first,
                                        most_dominating_fe_index);
                                    const types::global_dof_index
                                      dependent_dof_index =
                                        line->dof_index(identity.second,
                                                        other_fe_index);

                                    // check if we are on an interface between
                                    // a locally owned and a ghost cell on which
                                    // we need to work on.
                                    //
                                    // all degrees of freedom belonging to
                                    // dominating FE indices or to a processor
                                    // with a higher rank have been set at this
                                    // point (either in Phase 2, or after the
                                    // first ghost exchange in Phase 5). thus,
                                    // we only have to set the indices of
                                    // degrees of freedom that have been
                                    // previously flagged invalid.
                                    if ((dependent_dof_index ==
                                         numbers::invalid_dof_index) &&
                                        (primary_dof_index !=
                                         numbers::invalid_dof_index))
                                      line->set_dof_index(identity.second,
                                                          primary_dof_index,
                                                          other_fe_index);
                                  }
                              }
                        }
                    }
                }

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .load_user_flags_line(user_flags);
        }



        /**
         * Merge invalid DoF indices on quads located on ghost interfaces
         * by a dominating valid one.
         */
        template <int dim, int spacedim>
        static void
        merge_invalid_quad_dofs_on_ghost_interfaces(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          (void)dof_handler;
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          // this function should only be called for dim<3 where there are
          // no quad dof identies. for dim>=3, the specialization below should
          // take care of it
          Assert(dim < 3, ExcInternalError());
        }


        template <int spacedim>
        static void merge_invalid_quad_dofs_on_ghost_interfaces(
          DoFHandler<3, spacedim> &dof_handler)
        {
          Assert(dof_handler.hp_capability_enabled == true,
                 (typename DoFHandler<3, spacedim>::ExcOnlyAvailableWithHP()));

          const int dim = 3;

          // we will mark quads that we have already treated, so first
          // save and clear the user flags on quads and later restore
          // them
          std::vector<bool> user_flags;
          dof_handler.get_triangulation().save_user_flags_quad(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .clear_user_flags_quad();

          // mark all quads on ghost cells
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_ghost())
              for (const auto q : cell->face_indices())
                cell->quad(q)->set_user_flag();

          // An implementation of the algorithm described in the hp-
          // paper, including the modification mentioned later in the
          // "complications in 3-d" subsections
          //
          // as explained there, we do something only if there are
          // exactly 2 finite elements associated with an object. if
          // there is only one, then there is nothing to do anyway,
          // and if there are 3 or more, then we can get into
          // trouble. note that this only happens for lines in 3d and
          // higher, and for quads only in 4d and higher, so this
          // isn't a particularly frequent case
          dealii::Table<3, std::unique_ptr<DoFIdentities>> quad_dof_identities(
            dof_handler.fe_collection.size(),
            dof_handler.fe_collection.size(),
            2 /*triangle (0) or quadrilateral (1)*/);

          for (const auto &cell : dof_handler.active_cell_iterators())
            for (const auto q : cell->face_indices())
              if ((cell->is_locally_owned()) &&
                  (cell->quad(q)->user_flag_set() == true) &&
                  (cell->quad(q)->n_active_fe_indices() == 2))
                {
                  const auto quad = cell->quad(q);
                  quad->clear_user_flag();

                  const std::set<unsigned int> fe_indices =
                    quad->get_active_fe_indices();

                  // find out which is the most dominating finite
                  // element of the ones that are used on this quad
                  const unsigned int most_dominating_fe_index =
                    dof_handler.get_fe_collection().find_dominating_fe(
                      fe_indices,
                      /*codim=*/dim - 2);

                  const unsigned int most_dominating_fe_index_face_no =
                    cell->active_fe_index() == most_dominating_fe_index ?
                      q :
                      cell->neighbor_face_no(q);

                  // if we found the most dominating element, then use
                  // this to eliminate some of the degrees of freedom
                  // by identification. otherwise, the code that
                  // computes hanging node constraints will have to
                  // deal with it by computing appropriate constraints
                  // along this face/edge
                  if (most_dominating_fe_index != numbers::invalid_unsigned_int)
                    {
                      // loop over the indices of all the finite
                      // elements that are not dominating, and
                      // identify their dofs to the most dominating
                      // one
                      for (const auto &other_fe_index : fe_indices)
                        if (other_fe_index != most_dominating_fe_index)
                          {
                            const auto &identities =
                              *ensure_existence_and_return_dof_identities<2>(
                                dof_handler.get_fe(most_dominating_fe_index),
                                dof_handler.get_fe(other_fe_index),
                                quad_dof_identities
                                  [most_dominating_fe_index][other_fe_index]
                                  [cell->quad(q)->reference_cell() ==
                                   dealii::ReferenceCells::Quadrilateral],
                                most_dominating_fe_index_face_no);

                            for (const auto &identity : identities)
                              {
                                const types::global_dof_index
                                  primary_dof_index =
                                    quad->dof_index(identity.first,
                                                    most_dominating_fe_index);
                                const types::global_dof_index
                                  dependent_dof_index =
                                    quad->dof_index(identity.second,
                                                    other_fe_index);

                                // check if we are on an interface between
                                // a locally owned and a ghost cell on which
                                // we need to work on.
                                //
                                // all degrees of freedom belonging to
                                // dominating FE indices or to a processor with
                                // a higher rank have been set at this point
                                // (either in Phase 2, or after the first ghost
                                // exchange in Phase 5). thus, we only have to
                                // set the indices of degrees of freedom that
                                // have been previously flagged invalid.
                                if ((dependent_dof_index ==
                                     numbers::invalid_dof_index) &&
                                    (primary_dof_index !=
                                     numbers::invalid_dof_index))
                                  quad->set_dof_index(identity.second,
                                                      primary_dof_index,
                                                      other_fe_index);
                              }
                          }
                    }
                }

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation())
            .load_user_flags_quad(user_flags);
        }



        /**
         * After DoF unification in Phase 2, we may still have invalid DoF
         * indices left on ghost interfaces that are dominated by a
         * FiniteElement object of an adjacent cell, which could be either a
         * ghost cell or locally owned, whose DoF indices we did not know at the
         * moment of unification. After the first ghost exchange in Phase 5, we
         * know the indices of all dominating DoFs, and we have to assign those
         * invalid entries to their corresponding global value.
         *
         * This function does nothing unless the DoFHandler has hp-
         * capabilities.
         */
        template <int dim, int spacedim>
        static void
        merge_invalid_dof_indices_on_ghost_interfaces(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          if (dof_handler.hp_capability_enabled == false)
            return;

          {
            Threads::TaskGroup<> tasks;

            tasks += Threads::new_task([&]() {
              merge_invalid_vertex_dofs_on_ghost_interfaces(dof_handler);
            });

            if (dim > 1)
              {
                tasks += Threads::new_task([&]() {
                  merge_invalid_line_dofs_on_ghost_interfaces(dof_handler);
                });
              }

            if (dim > 2)
              {
                tasks += Threads::new_task([&]() {
                  merge_invalid_quad_dofs_on_ghost_interfaces(dof_handler);
                });
              }

            tasks.join_all();
          }

          update_all_active_cell_dof_indices_caches(dof_handler);
        }



        /**
         * Distribute degrees of freedom on all cells, or on cells with the
         * correct subdomain_id if the corresponding argument is not equal to
         * numbers::invalid_subdomain_id. Return the total number of dofs
         * distributed.
         */
        template <int dim, int spacedim>
        static types::global_dof_index
        distribute_dofs(const types::subdomain_id  subdomain_id,
                        DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(dof_handler.get_triangulation().n_levels() > 0,
                 ExcMessage("Empty triangulation"));

          // Step 1: distribute dofs on all cells, but definitely
          // exclude artificial cells
          types::global_dof_index next_free_dof = 0;

          std::vector<types::global_dof_index> dof_indices;

          for (auto cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              if ((subdomain_id == numbers::invalid_subdomain_id) ||
                  (cell->subdomain_id() == subdomain_id))
                {
                  dof_indices.resize(cell->get_fe().n_dofs_per_cell());

                  // circumvent cache
                  internal::DoFAccessorImplementation::Implementation::
                    get_dof_indices(*cell,
                                    dof_indices,
                                    cell->active_fe_index());

                  for (auto &dof_index : dof_indices)
                    if (dof_index == numbers::invalid_dof_index)
                      dof_index = next_free_dof++;

                  cell->set_dof_indices(dof_indices);
                }

          update_all_active_cell_dof_indices_caches(dof_handler);

          return next_free_dof;
        }



        /**
         * During DoF distribution, DoFs on ghost interfaces get different
         * indices assigned by each adjacent subdomain. We need to clarify
         * ownership of those DoFs by imposing a criterion, which is that
         * the adjacent subdomain belonging to the processor of lowest rank
         * prescribes its index.
         *
         * Thus, we have to invalidate all those DoFs on ghost cells that
         * belong to processors of lower rank than the current one, which is
         * the @p subdomain_id parameter. Later during ghost exchange in
         * Phase 5, these values will be overwritten by the correct one. The
         * invalidated indices will be stored in the @p renumbering parameter.
         */
        template <int dim, int spacedim>
        static void
        invalidate_dof_indices_on_weaker_ghost_cells_for_renumbering(
          std::vector<types::global_dof_index> &renumbering,
          const types::subdomain_id             subdomain_id,
          const DoFHandler<dim, spacedim> &     dof_handler)
        {
          std::vector<types::global_dof_index> local_dof_indices;

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_ghost() && (cell->subdomain_id() < subdomain_id))
              {
                // we found a neighboring ghost cell whose subdomain
                // is "stronger" than our own subdomain

                // delete all dofs that live there and that we have
                // previously assigned a number to (i.e. the ones on
                // the interface)
                local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
                cell->get_dof_indices(local_dof_indices);
                for (const auto &local_dof_index : local_dof_indices)
                  if (local_dof_index != numbers::invalid_dof_index)
                    renumbering[local_dof_index] = numbers::invalid_dof_index;
              }
        }



        /* -------------- distribute_mg_dofs functionality ------------- */



        template <int dim, int spacedim>
        static types::global_dof_index
        distribute_dofs_on_level(const types::subdomain_id  level_subdomain_id,
                                 DoFHandler<dim, spacedim> &dof_handler,
                                 const unsigned int         level)
        {
          Assert(dof_handler.hp_capability_enabled == false,
                 ExcInternalError());

          const dealii::Triangulation<dim, spacedim> &tria =
            dof_handler.get_triangulation();
          Assert(tria.n_levels() > 0, ExcMessage("Empty triangulation"));
          if (level >= tria.n_levels())
            return 0; // this is allowed for multigrid

          types::global_dof_index next_free_dof = 0;

          std::vector<types::global_dof_index> dof_indices;

          for (auto cell : dof_handler.cell_iterators_on_level(level))
            if ((level_subdomain_id == numbers::invalid_subdomain_id) ||
                (cell->level_subdomain_id() == level_subdomain_id))
              {
                dof_indices.resize(cell->get_fe().n_dofs_per_cell());

                cell->get_mg_dof_indices(dof_indices);

                for (auto &dof_index : dof_indices)
                  if (dof_index == numbers::invalid_dof_index)
                    dof_index = next_free_dof++;

                cell->set_mg_dof_indices(dof_indices);
              }

          return next_free_dof;
        }



        /* --------------------- renumber_dofs functionality ---------------- */


        /**
         * The part of the renumber_dofs() functionality that operates on faces.
         * This part is dimension dependent and so needs to be implemented in
         * three separate specializations of the function.
         *
         * See renumber_dofs() for the meaning of the arguments.
         */
        template <int dim, int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index> &new_numbers,
          const IndexSet &                            indices_we_care_about,
          DoFHandler<dim, spacedim> &                 dof_handler)
        {
          for (unsigned int d = 1; d < dim; d++)
            for (auto &i : dof_handler.object_dof_indices[0][d])
              if (i != numbers::invalid_dof_index)
                i = ((indices_we_care_about.size() == 0) ?
                       new_numbers[i] :
                       new_numbers[indices_we_care_about.index_within_set(i)]);
        }



        template <int dim, int spacedim>
        static void
        renumber_vertex_dofs(
          const std::vector<types::global_dof_index> &new_numbers,
          const IndexSet &                            indices_we_care_about,
          DoFHandler<dim, spacedim> &                 dof_handler,
          const bool                                  check_validity)
        {
          if (dof_handler.hp_capability_enabled == false)
            {
              // we can not use cell iterators in this function since then
              // we would renumber the dofs on the interface of two cells
              // more than once. Anyway, this way it's not only more
              // correct but also faster; note, however, that dof numbers
              // may be invalid_dof_index, namely when the appropriate
              // vertex/line/etc is unused
              for (std::vector<types::global_dof_index>::iterator i =
                     dof_handler.object_dof_indices[0][0].begin();
                   i != dof_handler.object_dof_indices[0][0].end();
                   ++i)
                if (*i != numbers::invalid_dof_index)
                  *i =
                    (indices_we_care_about.size() == 0) ?
                      (new_numbers[*i]) :
                      (new_numbers[indices_we_care_about.index_within_set(*i)]);
                else if (check_validity)
                  // if index is invalid_dof_index: check if this one
                  // really is unused
                  Assert(dof_handler.get_triangulation().vertex_used(
                           (i - dof_handler.object_dof_indices[0][0].begin()) /
                           dof_handler.get_fe().n_dofs_per_vertex()) == false,
                         ExcInternalError());
              return;
            }


          for (unsigned int vertex_index = 0;
               vertex_index < dof_handler.get_triangulation().n_vertices();
               ++vertex_index)
            {
              const unsigned int n_active_fe_indices =
                dealii::internal::DoFAccessorImplementation::Implementation::
                  n_active_fe_indices(dof_handler,
                                      0,
                                      vertex_index,
                                      std::integral_constant<int, 0>());

              // if this vertex is unused, then we really ought not to have
              // allocated any space for it, i.e., n_active_fe_indices should be
              // zero, and there is no space to actually store dof indices for
              // this vertex
              if (dof_handler.get_triangulation().vertex_used(vertex_index) ==
                  false)
                Assert(n_active_fe_indices == 0, ExcInternalError());

              // otherwise the vertex is used; it may still not hold any dof
              // indices if it is located on an artificial cell and not adjacent
              // to a ghost cell, but in that case there is simply nothing for
              // us to do
              for (unsigned int f = 0; f < n_active_fe_indices; ++f)
                {
                  const unsigned int fe_index =
                    dealii::internal::DoFAccessorImplementation::
                      Implementation::nth_active_fe_index(
                        dof_handler,
                        0,
                        vertex_index,
                        f,
                        std::integral_constant<int, 0>());

                  for (unsigned int d = 0;
                       d < dof_handler.get_fe(fe_index).n_dofs_per_vertex();
                       ++d)
                    {
                      const types::global_dof_index old_dof_index =
                        dealii::internal::DoFAccessorImplementation::
                          Implementation::get_dof_index(
                            dof_handler,
                            0,
                            vertex_index,
                            fe_index,
                            d,
                            std::integral_constant<int, 0>());

                      // if check_validity was set, then we are to verify that
                      // the previous indices were all valid. this really should
                      // be the case: we allocated space for these vertex dofs,
                      // i.e., at least one adjacent cell has a valid
                      // active FE index, so there are DoFs that really live
                      // on this vertex. if check_validity is set, then we
                      // must make sure that they have been set to something
                      // useful
                      if (check_validity)
                        Assert(old_dof_index != numbers::invalid_dof_index,
                               ExcInternalError());

                      if (old_dof_index != numbers::invalid_dof_index)
                        {
                          // In the following blocks, we first check whether
                          // we were given an IndexSet of DoFs to touch. If not
                          // (the first 'if' case here), then we are in the
                          // sequential case and are allowed to touch all DoFs.
                          //
                          // If yes (the 'else' case), then we need to
                          // distinguish whether the DoF whose number we want to
                          // touch is in fact locally owned (i.e., is in the
                          // index set) and then we can actually assign it a new
                          // number; otherwise, we have encountered a
                          // non-locally owned DoF for which we don't know the
                          // new number yet and so set it to an invalid index.
                          // This will later be fixed up after the first ghost
                          // exchange phase when we unify hp-DoFs on neighboring
                          // cells.
                          if (indices_we_care_about.size() == 0)
                            dealii::internal::DoFAccessorImplementation::
                              Implementation::set_dof_index(
                                dof_handler,
                                0,
                                vertex_index,
                                fe_index,
                                d,
                                std::integral_constant<int, 0>(),
                                new_numbers[old_dof_index]);
                          else
                            {
                              if (indices_we_care_about.is_element(
                                    old_dof_index))
                                dealii::internal::DoFAccessorImplementation::
                                  Implementation::set_dof_index(
                                    dof_handler,
                                    0,
                                    vertex_index,
                                    fe_index,
                                    d,
                                    std::integral_constant<int, 0>(),
                                    new_numbers[indices_we_care_about
                                                  .index_within_set(
                                                    old_dof_index)]);
                              else
                                dealii::internal::DoFAccessorImplementation::
                                  Implementation::set_dof_index(
                                    dof_handler,
                                    0,
                                    vertex_index,
                                    fe_index,
                                    d,
                                    std::integral_constant<int, 0>(),
                                    numbers::invalid_dof_index);
                            }
                        }
                    }
                }
            }
        }



        template <int dim, int spacedim>
        static void
        renumber_cell_dofs(
          const std::vector<types::global_dof_index> &new_numbers,
          const IndexSet &                            indices_we_care_about,
          DoFHandler<dim, spacedim> &                 dof_handler)
        {
          if (dof_handler.hp_capability_enabled == false)
            {
              for (unsigned int level = 0;
                   level < dof_handler.object_dof_indices.size();
                   ++level)
                for (auto &i : dof_handler.object_dof_indices[level][dim])
                  if (i != numbers::invalid_dof_index)
                    i = ((indices_we_care_about.size() == 0) ?
                           new_numbers[i] :
                           new_numbers[indices_we_care_about.index_within_set(
                             i)]);
              return;
            }

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              {
                const unsigned int fe_index = cell->active_fe_index();

                for (unsigned int d = 0;
                     d < dof_handler.get_fe(fe_index)
                           .template n_dofs_per_object<dim>();
                     ++d)
                  {
                    const types::global_dof_index old_dof_index =
                      cell->dof_index(d, fe_index);
                    if (old_dof_index != numbers::invalid_dof_index)
                      {
                        // In the following blocks, we first check whether
                        // we were given an IndexSet of DoFs to touch. If not
                        // (the first 'if' case here), then we are in the
                        // sequential case and are allowed to touch all DoFs.
                        //
                        // If yes (the 'else' case), then we need to distinguish
                        // whether the DoF whose number we want to touch is in
                        // fact locally owned (i.e., is in the index set) and
                        // then we can actually assign it a new number;
                        // otherwise, we have encountered a non-locally owned
                        // DoF for which we don't know the new number yet and so
                        // set it to an invalid index. This will later be fixed
                        // up after the first ghost exchange phase when we unify
                        // hp-DoFs on neighboring cells.
                        if (indices_we_care_about.size() == 0)
                          cell->set_dof_index(d,
                                              new_numbers[old_dof_index],
                                              fe_index);
                        else
                          {
                            if (indices_we_care_about.is_element(old_dof_index))
                              cell->set_dof_index(
                                d,
                                new_numbers[indices_we_care_about
                                              .index_within_set(old_dof_index)],
                                fe_index);
                            else
                              cell->set_dof_index(d,
                                                  numbers::invalid_dof_index,
                                                  fe_index);
                          }
                      }
                  }
              }
        }



        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index> & /*new_numbers*/,
          const IndexSet & /*indices_we_care_about*/,
          DoFHandler<1, spacedim> & /*dof_handler*/)
        {
          // nothing to do in 1d since there are no separate faces -- we've
          // already taken care of this when dealing with the vertices
        }



        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index> &new_numbers,
          const IndexSet &                            indices_we_care_about,
          DoFHandler<2, spacedim> &                   dof_handler)
        {
          const unsigned int dim = 2;

          if (dof_handler.hp_capability_enabled == false)
            {
              for (unsigned int d = 1; d < dim; d++)
                for (auto &i : dof_handler.object_dof_indices[0][d])
                  if (i != numbers::invalid_dof_index)
                    i = ((indices_we_care_about.size() == 0) ?
                           new_numbers[i] :
                           new_numbers[indices_we_care_about.index_within_set(
                             i)]);
              return;
            }

          // deal with DoFs on lines
          {
            // save user flags on lines so we can use them to mark lines
            // we've already treated
            std::vector<bool> saved_line_user_flags;
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .save_user_flags_line(saved_line_user_flags);
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .clear_user_flags_line();

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const auto l : cell->line_indices())
                  if (cell->line(l)->user_flag_set() == false)
                    {
                      const auto line = cell->line(l);
                      line->set_user_flag();

                      const unsigned int n_active_fe_indices =
                        line->n_active_fe_indices();

                      for (unsigned int f = 0; f < n_active_fe_indices; ++f)
                        {
                          const unsigned int fe_index =
                            line->nth_active_fe_index(f);

                          for (unsigned int d = 0;
                               d <
                               dof_handler.get_fe(fe_index).n_dofs_per_line();
                               ++d)
                            {
                              const types::global_dof_index old_dof_index =
                                line->dof_index(d, fe_index);
                              if (old_dof_index != numbers::invalid_dof_index)
                                {
                                  // In the following blocks, we first check
                                  // whether we were given an IndexSet of DoFs
                                  // to touch. If not (the first 'if' case
                                  // here), then we are in the sequential case
                                  // and are allowed to touch all DoFs.
                                  //
                                  // If yes (the 'else' case), then we need to
                                  // distinguish whether the DoF whose number we
                                  // want to touch is in fact locally owned
                                  // (i.e., is in the index set) and then we can
                                  // actually assign it a new number; otherwise,
                                  // we have encountered a non-locally owned DoF
                                  // for which we don't know the new number yet
                                  // and so set it to an invalid index. This
                                  // will later be fixed up after the first
                                  // ghost exchange phase when we unify hp-DoFs
                                  // on neighboring cells.
                                  if (indices_we_care_about.size() == 0)
                                    line->set_dof_index(
                                      d, new_numbers[old_dof_index], fe_index);
                                  else
                                    {
                                      if (indices_we_care_about.is_element(
                                            old_dof_index))
                                        line->set_dof_index(
                                          d,
                                          new_numbers[indices_we_care_about
                                                        .index_within_set(
                                                          old_dof_index)],
                                          fe_index);
                                      else
                                        line->set_dof_index(
                                          d,
                                          numbers::invalid_dof_index,
                                          fe_index);
                                    }
                                }
                            }
                        }
                    }

            // at the end, restore the user
            // flags for the lines
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .load_user_flags_line(saved_line_user_flags);
          }
        }



        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index> &new_numbers,
          const IndexSet &                            indices_we_care_about,
          DoFHandler<3, spacedim> &                   dof_handler)
        {
          const unsigned int dim = 3;

          if (dof_handler.hp_capability_enabled == false)
            {
              for (unsigned int d = 1; d < dim; d++)
                for (auto &i : dof_handler.object_dof_indices[0][d])
                  if (i != numbers::invalid_dof_index)
                    i = ((indices_we_care_about.size() == 0) ?
                           new_numbers[i] :
                           new_numbers[indices_we_care_about.index_within_set(
                             i)]);
              return;
            }

          // deal with DoFs on lines
          {
            // save user flags on lines so we can use them to mark lines
            // we've already treated
            std::vector<bool> saved_line_user_flags;
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .save_user_flags_line(saved_line_user_flags);
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .clear_user_flags_line();

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const auto l : cell->line_indices())
                  if (cell->line(l)->user_flag_set() == false)
                    {
                      const auto line = cell->line(l);
                      line->set_user_flag();

                      const unsigned int n_active_fe_indices =
                        line->n_active_fe_indices();

                      for (unsigned int f = 0; f < n_active_fe_indices; ++f)
                        {
                          const unsigned int fe_index =
                            line->nth_active_fe_index(f);

                          for (unsigned int d = 0;
                               d <
                               dof_handler.get_fe(fe_index).n_dofs_per_line();
                               ++d)
                            {
                              const types::global_dof_index old_dof_index =
                                line->dof_index(d, fe_index);
                              if (old_dof_index != numbers::invalid_dof_index)
                                {
                                  // In the following blocks, we first check
                                  // whether we were given an IndexSet of DoFs
                                  // to touch. If not (the first 'if' case
                                  // here), then we are in the sequential case
                                  // and are allowed to touch all DoFs.
                                  //
                                  // If yes (the 'else' case), then we need to
                                  // distinguish whether the DoF whose number we
                                  // want to touch is in fact locally owned
                                  // (i.e., is in the index set) and then we can
                                  // actually assign it a new number; otherwise,
                                  // we have encountered a non-locally owned DoF
                                  // for which we don't know the new number yet
                                  // and so set it to an invalid index. This
                                  // will later be fixed up after the first
                                  // ghost exchange phase when we unify hp-DoFs
                                  // on neighboring cells.
                                  if (indices_we_care_about.size() == 0)
                                    line->set_dof_index(
                                      d, new_numbers[old_dof_index], fe_index);
                                  else if (indices_we_care_about.is_element(
                                             old_dof_index))
                                    line->set_dof_index(
                                      d,
                                      new_numbers[indices_we_care_about
                                                    .index_within_set(
                                                      old_dof_index)],
                                      fe_index);
                                  else
                                    line->set_dof_index(
                                      d, numbers::invalid_dof_index, fe_index);
                                }
                            }
                        }
                    }

            // at the end, restore the user
            // flags for the lines
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .load_user_flags_line(saved_line_user_flags);
          }

          // then deal with dofs on quads
          {
            std::vector<bool> saved_quad_user_flags;
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .save_user_flags_quad(saved_quad_user_flags);
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .clear_user_flags_quad();

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const auto q : cell->face_indices())
                  if (cell->quad(q)->user_flag_set() == false)
                    {
                      const auto quad = cell->quad(q);
                      quad->set_user_flag();

                      const unsigned int n_active_fe_indices =
                        quad->n_active_fe_indices();

                      for (unsigned int f = 0; f < n_active_fe_indices; ++f)
                        {
                          const unsigned int fe_index =
                            quad->nth_active_fe_index(f);

                          for (unsigned int d = 0;
                               d <
                               dof_handler.get_fe(fe_index).n_dofs_per_quad(q);
                               ++d)
                            {
                              const types::global_dof_index old_dof_index =
                                quad->dof_index(d, fe_index);
                              if (old_dof_index != numbers::invalid_dof_index)
                                {
                                  // In the following blocks, we first check
                                  // whether we were given an IndexSet of DoFs
                                  // to touch. If not (the first 'if' case
                                  // here), then we are in the sequential case
                                  // and are allowed to touch all DoFs.
                                  //
                                  // If yes (the 'else' case), then we need to
                                  // distinguish whether the DoF whose number we
                                  // want to touch is in fact locally owned
                                  // (i.e., is in the index set) and then we can
                                  // actually assign it a new number; otherwise,
                                  // we have encountered a non-locally owned DoF
                                  // for which we don't know the new number yet
                                  // and so set it to an invalid index. This
                                  // will later be fixed up after the first
                                  // ghost exchange phase when we unify hp-DoFs
                                  // on neighboring cells.
                                  if (indices_we_care_about.size() == 0)
                                    quad->set_dof_index(
                                      d, new_numbers[old_dof_index], fe_index);
                                  else
                                    {
                                      if (indices_we_care_about.is_element(
                                            old_dof_index))
                                        quad->set_dof_index(
                                          d,
                                          new_numbers[indices_we_care_about
                                                        .index_within_set(
                                                          old_dof_index)],
                                          fe_index);
                                      else
                                        quad->set_dof_index(
                                          d,
                                          numbers::invalid_dof_index,
                                          fe_index);
                                    }
                                }
                            }
                        }
                    }

            // at the end, restore the user flags for the quads
            const_cast<dealii::Triangulation<dim, spacedim> &>(
              dof_handler.get_triangulation())
              .load_user_flags_quad(saved_quad_user_flags);
          }
        }



        /**
         * Implementation of DoFHandler::renumber_dofs()
         *
         * If the second argument has any elements set, elements of
         * the then the vector of new numbers do not relate to the old
         * DoF number but instead to the index of the old DoF number
         * within the set of locally owned DoFs.
         *
         * (The IndexSet argument is not used in 1d because we only need
         * it for parallel meshes and 1d doesn't support that right now.)
         */
        template <int dim, int space_dim>
        static void
        renumber_dofs(const std::vector<types::global_dof_index> &new_numbers,
                      const IndexSet &                  indices_we_care_about,
                      const DoFHandler<dim, space_dim> &dof_handler,
                      const bool                        check_validity)
        {
          if (dim == 1)
            Assert(indices_we_care_about == IndexSet(0), ExcNotImplemented());

          // renumber DoF indices on vertices, cells, and faces. this
          // can be done in parallel because the respective functions
          // work on separate data structures
          Threads::TaskGroup<> tasks;
          tasks += Threads::new_task([&]() {
            renumber_vertex_dofs(new_numbers,
                                 indices_we_care_about,
                                 const_cast<DoFHandler<dim, space_dim> &>(
                                   dof_handler),
                                 check_validity);
          });
          tasks += Threads::new_task([&]() {
            renumber_face_dofs(new_numbers,
                               indices_we_care_about,
                               const_cast<DoFHandler<dim, space_dim> &>(
                                 dof_handler));
          });
          tasks += Threads::new_task([&]() {
            renumber_cell_dofs(new_numbers,
                               indices_we_care_about,
                               const_cast<DoFHandler<dim, space_dim> &>(
                                 dof_handler));
          });
          tasks.join_all();

          // update the cache used for cell dof indices
          update_all_active_cell_dof_indices_caches(
            const_cast<DoFHandler<dim, space_dim> &>(dof_handler));
        }



        /* --------------------- renumber_mg_dofs functionality ----------------
         */

        /**
         * The part of the renumber_mg_dofs() functionality that is dimension
         * independent because it renumbers the DoF indices on vertex interiors
         * (which exist for all dimensions).
         *
         * See renumber_mg_dofs() for the meaning of the arguments.
         */
        template <int dim, int spacedim>
        static void
        renumber_vertex_mg_dofs(
          const std::vector<dealii::types::global_dof_index> &new_numbers,
          const IndexSet &           indices_we_care_about,
          DoFHandler<dim, spacedim> &dof_handler,
          const unsigned int         level,
          const bool                 check_validity)
        {
          (void)check_validity;
          Assert(level < dof_handler.get_triangulation().n_levels(),
                 ExcInternalError());

          for (typename std::vector<
                 typename DoFHandler<dim, spacedim>::MGVertexDoFs>::iterator i =
                 dof_handler.mg_vertex_dofs.begin();
               i != dof_handler.mg_vertex_dofs.end();
               ++i)
            // if the present vertex lives on the current level
            if ((i->get_coarsest_level() <= level) &&
                (i->get_finest_level() >= level))
              for (unsigned int d = 0;
                   d < dof_handler.get_fe().n_dofs_per_vertex();
                   ++d)
                {
                  const dealii::types::global_dof_index idx =
                    i->get_index(level,
                                 d,
                                 dof_handler.get_fe().n_dofs_per_vertex());

                  if (idx != numbers::invalid_dof_index)
                    {
                      Assert(check_validity == false ||
                               (indices_we_care_about.size() > 0 ?
                                  indices_we_care_about.is_element(idx) :
                                  (idx < new_numbers.size())),
                             ExcInternalError());
                      i->set_index(level,
                                   d,
                                   dof_handler.get_fe().n_dofs_per_vertex(),
                                   (indices_we_care_about.size() == 0) ?
                                     (new_numbers[idx]) :
                                     (new_numbers[indices_we_care_about
                                                    .index_within_set(idx)]));
                    }
                }
        }



        /**
         * The part of the renumber_dofs() functionality that is dimension
         * independent because it renumbers the DoF indices on cell interiors
         * (which exist for all dimensions).
         *
         * See renumber_mg_dofs() for the meaning of the arguments.
         */
        template <int dim, int spacedim>
        static void
        renumber_cell_mg_dofs(
          const std::vector<dealii::types::global_dof_index> &new_numbers,
          const IndexSet &           indices_we_care_about,
          DoFHandler<dim, spacedim> &dof_handler,
          const unsigned int         level)
        {
          for (std::vector<types::global_dof_index>::iterator i =
                 dof_handler.mg_levels[level]->dof_object.dofs.begin();
               i != dof_handler.mg_levels[level]->dof_object.dofs.end();
               ++i)
            {
              if (*i != numbers::invalid_dof_index)
                {
                  Assert((indices_we_care_about.size() > 0 ?
                            indices_we_care_about.is_element(*i) :
                            (*i < new_numbers.size())),
                         ExcInternalError());
                  *i =
                    (indices_we_care_about.size() == 0) ?
                      (new_numbers[*i]) :
                      (new_numbers[indices_we_care_about.index_within_set(*i)]);
                }
            }
        }



        /**
         * The part of the renumber_mg_dofs() functionality that operates on
         * faces. This part is dimension dependent and so needs to be
         * implemented in three separate specializations of the function.
         *
         * See renumber_mg_dofs() for the meaning of the arguments.
         */
        template <int spacedim>
        static void
        renumber_face_mg_dofs(
          const std::vector<types::global_dof_index> & /*new_numbers*/,
          const IndexSet & /*indices_we_care_about*/,
          DoFHandler<1, spacedim> & /*dof_handler*/,
          const unsigned int /*level*/,
          const bool /*check_validity*/)
        {
          // nothing to do in 1d because there are no separate faces
        }



        template <int spacedim>
        static void
        renumber_face_mg_dofs(
          const std::vector<dealii::types::global_dof_index> &new_numbers,
          const IndexSet &         indices_we_care_about,
          DoFHandler<2, spacedim> &dof_handler,
          const unsigned int       level,
          const bool               check_validity)
        {
          if (dof_handler.get_fe().n_dofs_per_line() > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<2, spacedim> &>(
                dof_handler.get_triangulation())
                .clear_user_flags();

              // flag all lines adjacent to cells of the current
              // level, as those lines logically belong to the same
              // level as the cell, at least for for isotropic
              // refinement
              for (const auto &cell :
                   dof_handler.cell_iterators_on_level(level))
                if (cell->level_subdomain_id() !=
                    numbers::artificial_subdomain_id)
                  for (const unsigned int line : cell->face_indices())
                    cell->face(line)->set_user_flag();

              for (const auto &cell :
                   dof_handler.cell_iterators_on_level(level))
                for (const auto l : cell->line_indices())
                  if (cell->line(l)->user_flag_set())
                    {
                      for (unsigned int d = 0;
                           d < dof_handler.get_fe().n_dofs_per_line();
                           ++d)
                        {
                          const dealii::types::global_dof_index idx =
                            cell->line(l)->mg_dof_index(level, d);
                          if (check_validity)
                            Assert(idx != numbers::invalid_dof_index,
                                   ExcInternalError());

                          if (idx != numbers::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index(
                              level,
                              d,
                              ((indices_we_care_about.size() == 0) ?
                                 new_numbers[idx] :
                                 new_numbers[indices_we_care_about
                                               .index_within_set(idx)]));
                        }
                      cell->line(l)->clear_user_flag();
                    }
              // finally, restore user flags
              const_cast<dealii::Triangulation<2, spacedim> &>(
                dof_handler.get_triangulation())
                .load_user_flags(user_flags);
            }
        }



        template <int spacedim>
        static void
        renumber_face_mg_dofs(
          const std::vector<dealii::types::global_dof_index> &new_numbers,
          const IndexSet &         indices_we_care_about,
          DoFHandler<3, spacedim> &dof_handler,
          const unsigned int       level,
          const bool               check_validity)
        {
          if (dof_handler.get_fe().n_dofs_per_line() > 0 ||
              dof_handler.get_fe().max_dofs_per_quad() > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<3, spacedim> &>(
                dof_handler.get_triangulation())
                .clear_user_flags();

              // flag all lines adjacent to cells of the current
              // level, as those lines logically belong to the same
              // level as the cell, at least for isotropic refinement
              for (const auto &cell :
                   dof_handler.cell_iterators_on_level(level))
                if (cell->level_subdomain_id() !=
                    numbers::artificial_subdomain_id)
                  for (const auto line : cell->line_indices())
                    cell->line(line)->set_user_flag();

              for (const auto &cell :
                   dof_handler.cell_iterators_on_level(level))
                for (const auto l : cell->line_indices())
                  if (cell->line(l)->user_flag_set())
                    {
                      for (unsigned int d = 0;
                           d < dof_handler.get_fe().n_dofs_per_line();
                           ++d)
                        {
                          const dealii::types::global_dof_index idx =
                            cell->line(l)->mg_dof_index(level, d);
                          if (check_validity)
                            Assert(idx != numbers::invalid_dof_index,
                                   ExcInternalError());

                          if (idx != numbers::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index(
                              level,
                              d,
                              ((indices_we_care_about.size() == 0) ?
                                 new_numbers[idx] :
                                 new_numbers[indices_we_care_about
                                               .index_within_set(idx)]));
                        }
                      cell->line(l)->clear_user_flag();
                    }

              // flag all quads adjacent to cells of the current level, as
              // those quads logically belong to the same level as the cell,
              // at least for isotropic refinement
              for (const auto &cell :
                   dof_handler.cell_iterators_on_level(level))
                if (cell->level_subdomain_id() !=
                    numbers::artificial_subdomain_id)
                  for (const auto quad : cell->face_indices())
                    cell->quad(quad)->set_user_flag();

              for (const auto &cell : dof_handler.cell_iterators())
                for (const auto l : cell->face_indices())
                  if (cell->quad(l)->user_flag_set())
                    {
                      for (unsigned int d = 0;
                           d < dof_handler.get_fe().n_dofs_per_quad(l);
                           ++d)
                        {
                          const dealii::types::global_dof_index idx =
                            cell->quad(l)->mg_dof_index(level, d);
                          if (check_validity)
                            Assert(idx != numbers::invalid_dof_index,
                                   ExcInternalError());

                          if (idx != numbers::invalid_dof_index)
                            cell->quad(l)->set_mg_dof_index(
                              level,
                              d,
                              ((indices_we_care_about.size() == 0) ?
                                 new_numbers[idx] :
                                 new_numbers[indices_we_care_about
                                               .index_within_set(idx)]));
                        }
                      cell->quad(l)->clear_user_flag();
                    }

              // finally, restore user flags
              const_cast<dealii::Triangulation<3, spacedim> &>(
                dof_handler.get_triangulation())
                .load_user_flags(user_flags);
            }
        }



        template <int dim, int spacedim>
        static void
        renumber_mg_dofs(
          const std::vector<dealii::types::global_dof_index> &new_numbers,
          const IndexSet &           indices_we_care_about,
          DoFHandler<dim, spacedim> &dof_handler,
          const unsigned int         level,
          const bool                 check_validity)
        {
          Assert(
            dof_handler.hp_capability_enabled == false,
            (typename DoFHandler<dim, spacedim>::ExcNotImplementedWithHP()));

          Assert(level < dof_handler.get_triangulation().n_global_levels(),
                 ExcInternalError());

          // renumber DoF indices on vertices, cells, and faces. this
          // can be done in parallel because the respective functions
          // work on separate data structures
          Threads::TaskGroup<> tasks;
          tasks += Threads::new_task([&]() {
            renumber_vertex_mg_dofs(new_numbers,
                                    indices_we_care_about,
                                    dof_handler,
                                    level,
                                    check_validity);
          });
          tasks += Threads::new_task([&]() {
            renumber_face_mg_dofs(new_numbers,
                                  indices_we_care_about,
                                  dof_handler,
                                  level,
                                  check_validity);
          });
          tasks += Threads::new_task([&]() {
            renumber_cell_mg_dofs(new_numbers,
                                  indices_we_care_about,
                                  dof_handler,
                                  level);
          });
          tasks.join_all();
        }
      };



      /* --------------------- class Sequential ---------------- */



      template <int dim, int spacedim>
      Sequential<dim, spacedim>::Sequential(
        DoFHandler<dim, spacedim> &dof_handler)
        : dof_handler(&dof_handler)
      {}



      template <int dim, int spacedim>
      NumberCache
      Sequential<dim, spacedim>::distribute_dofs() const
      {
        const types::global_dof_index n_initial_dofs =
          Implementation::distribute_dofs(numbers::invalid_subdomain_id,
                                          *dof_handler);

        const types::global_dof_index n_dofs =
          Implementation::unify_dof_indices(*dof_handler,
                                            n_initial_dofs,
                                            /*check_validity=*/true);

        // return a sequential, complete index set
        return NumberCache(n_dofs);
      }



      template <int dim, int spacedim>
      std::vector<NumberCache>
      Sequential<dim, spacedim>::distribute_mg_dofs() const
      {
        std::vector<bool> user_flags;
        dof_handler->get_triangulation().save_user_flags(user_flags);

        const_cast<dealii::Triangulation<dim, spacedim> &>(
          dof_handler->get_triangulation())
          .clear_user_flags();

        std::vector<NumberCache> number_caches;
        number_caches.reserve(dof_handler->get_triangulation().n_levels());
        for (unsigned int level = 0;
             level < dof_handler->get_triangulation().n_levels();
             ++level)
          {
            // first distribute dofs on this level
            const types::global_dof_index n_level_dofs =
              Implementation::distribute_dofs_on_level(
                numbers::invalid_subdomain_id, *dof_handler, level);

            // then add a complete, sequential index set
            number_caches.emplace_back(n_level_dofs);
          }

        const_cast<dealii::Triangulation<dim, spacedim> &>(
          dof_handler->get_triangulation())
          .load_user_flags(user_flags);

        return number_caches;
      }



      template <int dim, int spacedim>
      NumberCache
      Sequential<dim, spacedim>::renumber_dofs(
        const std::vector<types::global_dof_index> &new_numbers) const
      {
        Implementation::renumber_dofs(new_numbers,
                                      IndexSet(0),
                                      *dof_handler,
                                      /*check_validity=*/true);

        // return a sequential, complete index set. take into account that the
        // number of DoF indices may in fact be smaller than there were before
        // if some previously separately numbered dofs have been identified.
        // this is, for example, what we do when the DoFHandler has hp-
        // capabilities enabled: it first enumerates all DoFs on cells
        // independently, and then unifies some located at vertices or faces;
        // this leaves us with fewer DoFs than there were before, so use the
        // largest index as the one to determine the size of the index space
        return NumberCache(
          *std::max_element(new_numbers.begin(), new_numbers.end()) + 1);
      }



      template <int dim, int spacedim>
      NumberCache
      Sequential<dim, spacedim>::renumber_mg_dofs(
        const unsigned int                          level,
        const std::vector<types::global_dof_index> &new_numbers) const
      {
        Implementation::renumber_mg_dofs(
          new_numbers, IndexSet(0), *dof_handler, level, true);

        // return a sequential, complete index set
        return NumberCache(new_numbers.size());
      }


      /* --------------------- class ParallelShared ---------------- */


      template <int dim, int spacedim>
      ParallelShared<dim, spacedim>::ParallelShared(
        DoFHandler<dim, spacedim> &dof_handler)
        : dof_handler(&dof_handler)
      {}



      namespace
      {
        /**
         * This function is a variation of
         * DoFTools::get_dof_subdomain_association() with the exception that it
         * is (i) specialized for parallel::shared::Triangulation objects, and
         * (ii) does not assume that the internal number cache of the DoFHandler
         * has already been set up. In can, consequently, be called at a point
         * when we are still distributing degrees of freedom.
         */
        template <int dim, int spacedim>
        std::vector<types::subdomain_id>
        get_dof_subdomain_association(
          const DoFHandler<dim, spacedim> &dof_handler,
          const types::global_dof_index    n_dofs,
          const unsigned int               n_procs)
        {
          (void)n_procs;
          std::vector<types::subdomain_id> subdomain_association(
            n_dofs, numbers::invalid_subdomain_id);
          std::vector<types::global_dof_index> local_dof_indices;
          local_dof_indices.reserve(
            dof_handler.get_fe_collection().max_dofs_per_cell());

          // loop over all cells and record which subdomain a DoF belongs to.
          // give to the smaller subdomain_id in case it is on an interface
          for (const auto &cell : dof_handler.active_cell_iterators())
            {
              // get the owner of the cell; note that we have made sure above
              // that all cells are either locally owned or ghosts (not
              // artificial), so this call will always yield the true owner
              const types::subdomain_id subdomain_id = cell->subdomain_id();
              const unsigned int        dofs_per_cell =
                cell->get_fe().n_dofs_per_cell();
              local_dof_indices.resize(dofs_per_cell);
              cell->get_dof_indices(local_dof_indices);

              // set subdomain ids. if dofs already have their values set then
              // they must be on partition interfaces. In that case assign them
              // to the processor with the smaller subdomain id.
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                if (subdomain_association[local_dof_indices[i]] ==
                    numbers::invalid_subdomain_id)
                  subdomain_association[local_dof_indices[i]] = subdomain_id;
                else if (subdomain_association[local_dof_indices[i]] >
                         subdomain_id)
                  {
                    subdomain_association[local_dof_indices[i]] = subdomain_id;
                  }
            }

          Assert(std::find(subdomain_association.begin(),
                           subdomain_association.end(),
                           numbers::invalid_subdomain_id) ==
                   subdomain_association.end(),
                 ExcInternalError());

          Assert(*std::max_element(subdomain_association.begin(),
                                   subdomain_association.end()) < n_procs,
                 ExcInternalError());

          return subdomain_association;
        }


        /**
         * level subdomain association. Similar to the above function only
         * for level meshes. This function assigns boundary dofs in
         * the same way as parallel::DistributedTriangulationBase (proc with
         * smallest index) instead of the coin flip method above.
         */
        template <int dim, int spacedim>
        std::vector<types::subdomain_id>
        get_dof_level_subdomain_association(
          const DoFHandler<dim, spacedim> &dof_handler,
          const types::global_dof_index    n_dofs_on_level,
          const unsigned int               n_procs,
          const unsigned int               level)
        {
          (void)n_procs;
          std::vector<types::subdomain_id> level_subdomain_association(
            n_dofs_on_level, numbers::invalid_subdomain_id);
          std::vector<types::global_dof_index> local_dof_indices;
          local_dof_indices.reserve(
            dof_handler.get_fe_collection().max_dofs_per_cell());

          // loop over all cells and record which subdomain a DoF belongs to.
          // interface goes to proccessor with smaller subdomain id
          for (const auto &cell : dof_handler.cell_iterators_on_level(level))
            {
              // get the owner of the cell; note that we have made sure above
              // that all cells are either locally owned or ghosts (not
              // artificial), so this call will always yield the true owner
              const types::subdomain_id level_subdomain_id =
                cell->level_subdomain_id();
              const unsigned int dofs_per_cell =
                cell->get_fe().n_dofs_per_cell();
              local_dof_indices.resize(dofs_per_cell);
              cell->get_mg_dof_indices(local_dof_indices);

              // set level subdomain ids. if dofs already have their values set
              // then they must be on partition interfaces. In that case assign
              // them to the processor with the smaller subdomain id.
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                if (level_subdomain_association[local_dof_indices[i]] ==
                    numbers::invalid_subdomain_id)
                  level_subdomain_association[local_dof_indices[i]] =
                    level_subdomain_id;
                else if (level_subdomain_association[local_dof_indices[i]] >
                         level_subdomain_id)
                  {
                    level_subdomain_association[local_dof_indices[i]] =
                      level_subdomain_id;
                  }
            }

          Assert(std::find(level_subdomain_association.begin(),
                           level_subdomain_association.end(),
                           numbers::invalid_subdomain_id) ==
                   level_subdomain_association.end(),
                 ExcInternalError());

          Assert(*std::max_element(level_subdomain_association.begin(),
                                   level_subdomain_association.end()) < n_procs,
                 ExcInternalError());

          return level_subdomain_association;
        }
      } // namespace



      template <int dim, int spacedim>
      NumberCache
      ParallelShared<dim, spacedim>::distribute_dofs() const
      {
        const dealii::parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<
            const dealii::parallel::shared::Triangulation<dim, spacedim> *>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        const unsigned int n_procs =
          Utilities::MPI::n_mpi_processes(tr->get_communicator());

        // If an underlying shared::Tria allows artificial cells, we need to
        // restore the true cell owners temporarily.
        // We use the TemporarilyRestoreSubdomainIds class for this purpose: we
        // save the current set of subdomain ids, set subdomain ids to the
        // "true" owner of each cell upon construction of the
        // TemporarilyRestoreSubdomainIds object, and later restore these flags
        // when it is destroyed.
        const internal::parallel::shared::
          TemporarilyRestoreSubdomainIds<dim, spacedim>
            subdomain_modifier(*tr);

        // first let the sequential algorithm do its magic. it is going to
        // enumerate DoFs on all cells, regardless of owner
        const types::global_dof_index n_initial_dofs =
          Implementation::distribute_dofs(numbers::invalid_subdomain_id,
                                          *this->dof_handler);

        const types::global_dof_index n_dofs =
          Implementation::unify_dof_indices(*this->dof_handler,
                                            n_initial_dofs,
                                            /*check_validity=*/true);

        // then re-enumerate them based on their subdomain association.
        // for this, we first have to identify for each current DoF
        // index which subdomain they belong to. ideally, we would
        // like to call DoFRenumbering::subdomain_wise(), but
        // because the NumberCache of the current DoFHandler is not
        // fully set up yet, we can't quite do that. also, that
        // function has to deal with other kinds of triangulations as
        // well, whereas we here know what kind of triangulation
        // we have and can simplify the code accordingly
        std::vector<types::global_dof_index> new_dof_indices(
          n_dofs, enumeration_dof_index);
        {
          // first get the association of each dof with a subdomain and
          // determine the total number of subdomain ids used
          const std::vector<types::subdomain_id> subdomain_association =
            get_dof_subdomain_association(*this->dof_handler, n_dofs, n_procs);

          // then renumber the subdomains by first looking at those belonging
          // to subdomain 0, then those of subdomain 1, etc. note that the
          // algorithm is stable, i.e. if two dofs i,j have i<j and belong to
          // the same subdomain, then they will be in this order also after
          // reordering
          types::global_dof_index next_free_index = 0;
          for (types::subdomain_id subdomain = 0; subdomain < n_procs;
               ++subdomain)
            for (types::global_dof_index i = 0; i < n_dofs; ++i)
              if (subdomain_association[i] == subdomain)
                {
                  Assert(new_dof_indices[i] == enumeration_dof_index,
                         ExcInternalError());
                  new_dof_indices[i] = next_free_index;
                  ++next_free_index;
                }

          // we should have numbered all dofs
          Assert(next_free_index == n_dofs, ExcInternalError());
          Assert(std::find(new_dof_indices.begin(),
                           new_dof_indices.end(),
                           enumeration_dof_index) == new_dof_indices.end(),
                 ExcInternalError());
        }
        // finally do the renumbering. we can use the sequential
        // version of the function because we do things on all
        // cells and all cells have their subdomain ids and DoFs
        // correctly set
        Implementation::renumber_dofs(new_dof_indices,
                                      IndexSet(0),
                                      *this->dof_handler,
                                      /*check_validity=*/true);

        // update the number cache. for this, we first have to find the
        // subdomain association for each DoF again following renumbering, from
        // which we can then compute the IndexSets of locally owned DoFs for all
        // processors. all other fields then follow from this
        //
        // given the way we enumerate degrees of freedom, the locally owned
        // ranges must all be contiguous and consecutive. this makes filling
        // the IndexSets cheap. an assertion at the top verifies that this
        // assumption is true
        const std::vector<types::subdomain_id> subdomain_association =
          get_dof_subdomain_association(*this->dof_handler, n_dofs, n_procs);

        for (unsigned int i = 1; i < n_dofs; ++i)
          Assert(subdomain_association[i] >= subdomain_association[i - 1],
                 ExcInternalError());

        std::vector<IndexSet> locally_owned_dofs_per_processor(
          n_procs, IndexSet(n_dofs));
        {
          // we know that the set of subdomain indices is contiguous from
          // the assertion above; find the start and end index for each
          // processor, taking into account that sometimes a processor
          // may not in fact have any DoFs at all. we do the latter
          // by just identifying contiguous ranges of subdomain_ids
          // and filling IndexSets for those subdomains; subdomains
          // that don't appear will lead to IndexSets that are simply
          // never touched and remain empty as initialized above.
          unsigned int start_index = 0;
          unsigned int end_index   = 0;
          while (start_index < n_dofs)
            {
              while ((end_index) < n_dofs &&
                     (subdomain_association[end_index] ==
                      subdomain_association[start_index]))
                ++end_index;

              // we've now identified a range of same indices. set that
              // range in the corresponding IndexSet
              if (end_index > start_index)
                {
                  const unsigned int subdomain_owner =
                    subdomain_association[start_index];
                  locally_owned_dofs_per_processor[subdomain_owner].add_range(
                    start_index, end_index);
                }

              // then move on to thinking about the next range
              start_index = end_index;
            }
        }

        // return a NumberCache object made up from the sets of locally
        // owned DoFs
        return NumberCache(
          locally_owned_dofs_per_processor,
          this->dof_handler->get_triangulation().locally_owned_subdomain());
      }



      template <int dim, int spacedim>
      std::vector<NumberCache>
      ParallelShared<dim, spacedim>::distribute_mg_dofs() const
      {
        const dealii::parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<
            const dealii::parallel::shared::Triangulation<dim, spacedim> *>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        const unsigned int n_procs =
          Utilities::MPI::n_mpi_processes(tr->get_communicator());
        const unsigned int n_levels = tr->n_global_levels();

        std::vector<NumberCache> number_caches;
        number_caches.reserve(n_levels);

        // We create an index set for each level
        for (unsigned int lvl = 0; lvl < n_levels; ++lvl)
          {
            // If the underlying shared::Tria allows artificial cells,
            // then save the current set of level subdomain ids, and set
            // subdomain ids to the "true" owner of each cell. we later
            // restore these flags
            // Note: "allows_artificial_cells" is currently enforced for
            //       MG computations.
            std::vector<types::subdomain_id> saved_level_subdomain_ids;
            saved_level_subdomain_ids.resize(tr->n_cells(lvl));
            {
              typename dealii::parallel::shared::Triangulation<dim, spacedim>::
                cell_iterator cell =
                                this->dof_handler->get_triangulation().begin(
                                  lvl),
                              endc =
                                this->dof_handler->get_triangulation().end(lvl);

              const std::vector<types::subdomain_id> &true_level_subdomain_ids =
                tr->get_true_level_subdomain_ids_of_cells(lvl);

              for (unsigned int index = 0; cell != endc; ++cell, ++index)
                {
                  saved_level_subdomain_ids[index] = cell->level_subdomain_id();
                  cell->set_level_subdomain_id(true_level_subdomain_ids[index]);
                }
            }

            // Next let the sequential algorithm do its magic. it is going to
            // enumerate DoFs on all cells on the given level, regardless of
            // owner
            const types::global_dof_index n_dofs_on_level =
              Implementation::distribute_dofs_on_level(
                numbers::invalid_subdomain_id, *this->dof_handler, lvl);

            // then re-enumerate them based on their level subdomain
            // association. for this, we first have to identify for each current
            // DoF index which subdomain they belong to. ideally, we would like
            // to call DoFRenumbering::subdomain_wise(), but because the
            // NumberCache of the current DoFHandler is not fully set up yet, we
            // can't quite do that. also, that function has to deal with other
            // kinds of triangulations as well, whereas we here know what kind
            // of triangulation we have and can simplify the code accordingly
            std::vector<types::global_dof_index> new_dof_indices(
              n_dofs_on_level, numbers::invalid_dof_index);
            {
              // first get the association of each dof with a subdomain and
              // determine the total number of subdomain ids used
              const std::vector<types::subdomain_id>
                level_subdomain_association =
                  get_dof_level_subdomain_association(*this->dof_handler,
                                                      n_dofs_on_level,
                                                      n_procs,
                                                      lvl);

              // then renumber the subdomains by first looking at those
              // belonging to subdomain 0, then those of subdomain 1, etc. note
              // that the algorithm is stable, i.e. if two dofs i,j have i<j and
              // belong to the same subdomain, then they will be in this order
              // also after reordering
              types::global_dof_index next_free_index = 0;
              for (types::subdomain_id level_subdomain = 0;
                   level_subdomain < n_procs;
                   ++level_subdomain)
                for (types::global_dof_index i = 0; i < n_dofs_on_level; ++i)
                  if (level_subdomain_association[i] == level_subdomain)
                    {
                      Assert(new_dof_indices[i] == numbers::invalid_dof_index,
                             ExcInternalError());
                      new_dof_indices[i] = next_free_index;
                      ++next_free_index;
                    }

              // we should have numbered all dofs
              Assert(next_free_index == n_dofs_on_level, ExcInternalError());
              Assert(std::find(new_dof_indices.begin(),
                               new_dof_indices.end(),
                               numbers::invalid_dof_index) ==
                       new_dof_indices.end(),
                     ExcInternalError());
            }

            // finally do the renumbering. we can use the sequential
            // version of the function because we do things on all
            // cells and all cells have their subdomain ids and DoFs
            // correctly set
            Implementation::renumber_mg_dofs(
              new_dof_indices, IndexSet(0), *this->dof_handler, lvl, true);

            // update the number cache. for this, we first have to find the
            // level subdomain association for each DoF again following
            // renumbering, from which we can then compute the IndexSets of
            // locally owned DoFs for all processors. all other fields then
            // follow from this
            //
            // given the way we enumerate degrees of freedom, the locally owned
            // ranges must all be contiguous and consecutive. this makes filling
            // the IndexSets cheap. an assertion at the top verifies that this
            // assumption is true
            const std::vector<types::subdomain_id> level_subdomain_association =
              get_dof_level_subdomain_association(*this->dof_handler,
                                                  n_dofs_on_level,
                                                  n_procs,
                                                  lvl);

            for (unsigned int i = 1; i < n_dofs_on_level; ++i)
              Assert(level_subdomain_association[i] >=
                       level_subdomain_association[i - 1],
                     ExcInternalError());

            std::vector<IndexSet> locally_owned_dofs_per_processor(
              n_procs, IndexSet(n_dofs_on_level));
            {
              // we know that the set of subdomain indices is contiguous from
              // the assertion above; find the start and end index for each
              // processor, taking into account that sometimes a processor
              // may not in fact have any DoFs at all. we do the latter
              // by just identifying contiguous ranges of level_subdomain_ids
              // and filling IndexSets for those subdomains; subdomains
              // that don't appear will lead to IndexSets that are simply
              // never touched and remain empty as initialized above.
              unsigned int start_index = 0;
              unsigned int end_index   = 0;
              while (start_index < n_dofs_on_level)
                {
                  while ((end_index) < n_dofs_on_level &&
                         (level_subdomain_association[end_index] ==
                          level_subdomain_association[start_index]))
                    ++end_index;

                  // we've now identified a range of same indices. set that
                  // range in the corresponding IndexSet
                  if (end_index > start_index)
                    {
                      const unsigned int level_subdomain_owner =
                        level_subdomain_association[start_index];
                      locally_owned_dofs_per_processor[level_subdomain_owner]
                        .add_range(start_index, end_index);
                    }

                  // then move on to thinking about the next range
                  start_index = end_index;
                }
            }

            // finally, restore current level subdomain ids
            {
              typename dealii::parallel::shared::Triangulation<dim, spacedim>::
                cell_iterator cell =
                                this->dof_handler->get_triangulation().begin(
                                  lvl),
                              endc =
                                this->dof_handler->get_triangulation().end(lvl);

              for (unsigned int index = 0; cell != endc; ++cell, ++index)
                cell->set_level_subdomain_id(saved_level_subdomain_ids[index]);

              // add NumberCache for current level
              number_caches.emplace_back(
                NumberCache(locally_owned_dofs_per_processor,
                            this->dof_handler->get_triangulation()
                              .locally_owned_subdomain()));
            }
          }

        return number_caches;
      }



      template <int dim, int spacedim>
      NumberCache
      ParallelShared<dim, spacedim>::renumber_dofs(
        const std::vector<types::global_dof_index> &new_numbers) const
      {
#ifndef DEAL_II_WITH_MPI
        (void)new_numbers;
        Assert(false, ExcNotImplemented());
        return NumberCache();
#else
        // Similar to distribute_dofs() we need to have a special treatment in
        // case artificial cells are present.
        const dealii::parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<
            const dealii::parallel::shared::Triangulation<dim, spacedim> *>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        // Set subdomain IDs to the "true" owner of each cell.
        const internal::parallel::shared::
          TemporarilyRestoreSubdomainIds<dim, spacedim>
            subdomain_modifier(*tr);

        std::vector<types::global_dof_index> global_gathered_numbers(
          this->dof_handler->n_dofs(), 0);
        // as we call DoFRenumbering::subdomain_wise (*dof_handler) from
        // distribute_dofs(), we need to support sequential-like input.
        // Distributed-like input from, for example, component_wise renumbering
        // is also supported.
        if (new_numbers.size() == this->dof_handler->n_dofs())
          {
            global_gathered_numbers = new_numbers;
          }
        else
          {
            Assert(new_numbers.size() ==
                     this->dof_handler->locally_owned_dofs().n_elements(),
                   ExcInternalError());
            const unsigned int n_cpu =
              Utilities::MPI::n_mpi_processes(tr->get_communicator());
            std::vector<types::global_dof_index> gathered_new_numbers(
              this->dof_handler->n_dofs(), 0);
            Assert(Utilities::MPI::this_mpi_process(tr->get_communicator()) ==
                     this->dof_handler->get_triangulation()
                       .locally_owned_subdomain(),
                   ExcInternalError())

            // gather new numbers among processors into one vector
            {
              std::vector<types::global_dof_index> new_numbers_copy(
                new_numbers);

              // store the number of elements that are to be received from each
              // process
              std::vector<int> rcounts(n_cpu);

              types::global_dof_index shift = 0;
              // set rcounts based on new_numbers:
              int cur_count = new_numbers_copy.size();
              int ierr      = MPI_Allgather(&cur_count,
                                       1,
                                       MPI_INT,
                                       rcounts.data(),
                                       1,
                                       MPI_INT,
                                       tr->get_communicator());
              AssertThrowMPI(ierr);

              // compute the displacements (relative to recvbuf)
              // at which to place the incoming data from process i
              std::vector<int> displacements(n_cpu);
              for (unsigned int i = 0; i < n_cpu; i++)
                {
                  displacements[i] = shift;
                  shift += rcounts[i];
                }
              Assert(new_numbers_copy.size() ==
                       static_cast<unsigned int>(
                         rcounts[Utilities::MPI::this_mpi_process(
                           tr->get_communicator())]),
                     ExcInternalError());
              ierr = MPI_Allgatherv(new_numbers_copy.data(),
                                    new_numbers_copy.size(),
                                    DEAL_II_DOF_INDEX_MPI_TYPE,
                                    gathered_new_numbers.data(),
                                    rcounts.data(),
                                    displacements.data(),
                                    DEAL_II_DOF_INDEX_MPI_TYPE,
                                    tr->get_communicator());
              AssertThrowMPI(ierr);
            }

            // put new numbers according to the current
            // locally_owned_dofs_per_processor IndexSets
            types::global_dof_index shift = 0;
            // flag_1 and flag_2 are
            // used to control that there is a
            // one-to-one relation between old and new DoFs.
            std::vector<unsigned int> flag_1(this->dof_handler->n_dofs(), 0);
            std::vector<unsigned int> flag_2(this->dof_handler->n_dofs(), 0);
            for (unsigned int i = 0; i < n_cpu; i++)
              {
                const IndexSet iset =
                  this->dof_handler->locally_owned_dofs_per_processor()[i];
                for (types::global_dof_index ind = 0; ind < iset.n_elements();
                     ind++)
                  {
                    const types::global_dof_index target =
                      iset.nth_index_in_set(ind);
                    const types::global_dof_index value =
                      gathered_new_numbers[shift + ind];
                    Assert(target < this->dof_handler->n_dofs(),
                           ExcInternalError());
                    Assert(value < this->dof_handler->n_dofs(),
                           ExcInternalError());
                    global_gathered_numbers[target] = value;
                    flag_1[target]++;
                    flag_2[value]++;
                  }
                shift += iset.n_elements();
              }

            Assert(*std::max_element(flag_1.begin(), flag_1.end()) == 1,
                   ExcInternalError());
            Assert(*std::min_element(flag_1.begin(), flag_1.end()) == 1,
                   ExcInternalError());
            Assert((*std::max_element(flag_2.begin(), flag_2.end())) == 1,
                   ExcInternalError());
            Assert((*std::min_element(flag_2.begin(), flag_2.end())) == 1,
                   ExcInternalError());
          }

        // let the sequential algorithm do its magic; ignore the
        // return type, but reconstruct the number cache based on
        // which DoFs each process owns
        Implementation::renumber_dofs(global_gathered_numbers,
                                      IndexSet(0),
                                      *this->dof_handler,
                                      /*check_validity=*/true);

        const NumberCache number_cache(
          DoFTools::locally_owned_dofs_per_subdomain(*this->dof_handler),
          this->dof_handler->get_triangulation().locally_owned_subdomain());

        return number_cache;
#endif
      }



      template <int dim, int spacedim>
      NumberCache
      ParallelShared<dim, spacedim>::renumber_mg_dofs(
        const unsigned int /*level*/,
        const std::vector<types::global_dof_index> & /*new_numbers*/) const
      {
        // multigrid is not currently implemented for shared triangulations
        Assert(false, ExcNotImplemented());

        return {};
      }



      /* --------------------- class ParallelDistributed ---------------- */

#ifdef DEAL_II_WITH_MPI

      namespace
      {
        template <int dim, int spacedim>
        void
        communicate_mg_ghost_cells(
          const typename dealii::parallel::
            DistributedTriangulationBase<dim, spacedim> &tria,
          DoFHandler<dim, spacedim> &                    dof_handler)
        {
          (void)tria;
          const auto pack = [&](const auto &cell) {
            // why would somebody request a cell that is not ours?
            Assert(cell->level_subdomain_id() == tria.locally_owned_subdomain(),
                   ExcInternalError());

            std::vector<dealii::types::global_dof_index> data(
              cell->get_fe().n_dofs_per_cell());
            cell->get_mg_dof_indices(data);

            return data;
          };

          const auto unpack = [](const auto &cell, const auto &dofs) {
            Assert(cell->get_fe().n_dofs_per_cell() == dofs.size(),
                   ExcInternalError());

            Assert(cell->level_subdomain_id() !=
                     dealii::numbers::artificial_subdomain_id,
                   ExcInternalError());

            std::vector<dealii::types::global_dof_index> dof_indices(
              cell->get_fe().n_dofs_per_cell());
            cell->get_mg_dof_indices(dof_indices);

            bool complete = true;
            for (unsigned int i = 0; i < dof_indices.size(); ++i)
              if (dofs[i] != numbers::invalid_dof_index)
                {
                  Assert((dof_indices[i] == (numbers::invalid_dof_index)) ||
                           (dof_indices[i] == dofs[i]),
                         ExcInternalError());
                  dof_indices[i] = dofs[i];
                }
              else
                complete = false;

            if (!complete)
              const_cast<
                typename DoFHandler<dim, spacedim>::level_cell_iterator &>(cell)
                ->set_user_flag();
            else
              const_cast<
                typename DoFHandler<dim, spacedim>::level_cell_iterator &>(cell)
                ->clear_user_flag();

            const_cast<
              typename DoFHandler<dim, spacedim>::level_cell_iterator &>(cell)
              ->set_mg_dof_indices(dof_indices);
          };

          const auto filter = [](const auto &cell) {
            return cell->user_flag_set();
          };

          GridTools::exchange_cell_data_to_level_ghosts<
            std::vector<types::global_dof_index>,
            DoFHandler<dim, spacedim>>(dof_handler, pack, unpack, filter);
        }



        template <int spacedim>
        void
        communicate_mg_ghost_cells(const typename dealii::parallel::
                                     distributed::Triangulation<1, spacedim> &,
                                   DoFHandler<1, spacedim> &)
        {
          Assert(false, ExcNotImplemented());
        }



        /**
         * A function that communicates the DoF indices from that subset of
         * locally owned cells that have their user indices set to the
         * corresponding ghost cells on other processors.
         *
         * This function makes use of the user flags in the following
         * way:
         * - On locally owned cells, the flag indicates whether we still
         *   need to send the DoF indices to other processors on which
         *   the current cell is a ghost. In phase 1, this is true for
         *   all locally owned cells that are adjacent to ghost cells
         *   in some way. In phase 2, this is only true if before phase
         *   1 we did not know all dof indices yet
         * - On ghost cells, the flag indicates whether we still expect
         *   information to be sent to us. In phase 1, this is true for
         *   all ghost cells. In phase 2, this is only true if we
         *   did not receive a complete set of DoF indices in phase 1.
         */
        template <int dim, int spacedim>
        void
        communicate_dof_indices_on_marked_cells(
          const DoFHandler<dim, spacedim> &dof_handler)
        {
#  ifndef DEAL_II_WITH_MPI
          (void)dof_handler;
          Assert(false, ExcNotImplemented());
#  else

          // define functions that pack data on cells that are ghost cells
          // somewhere else, and unpack data on cells where we get information
          // from elsewhere
          const auto pack = [](const auto &cell) {
            Assert(cell->is_locally_owned(), ExcInternalError());

            std::vector<dealii::types::global_dof_index> data(
              cell->get_fe().n_dofs_per_cell());
            cell->get_dof_indices(data);

            return data;
          };

          const auto unpack = [](const auto &cell, const auto &dofs) {
            Assert(cell->get_fe().n_dofs_per_cell() == dofs.size(),
                   ExcInternalError());

            Assert(cell->is_ghost(), ExcInternalError());

            std::vector<dealii::types::global_dof_index> dof_indices(
              cell->get_fe().n_dofs_per_cell());
            cell->update_cell_dof_indices_cache();
            cell->get_dof_indices(dof_indices);

            bool complete = true;
            for (unsigned int i = 0; i < dof_indices.size(); ++i)
              if (dofs[i] != numbers::invalid_dof_index)
                {
                  Assert((dof_indices[i] == (numbers::invalid_dof_index)) ||
                           (dof_indices[i] == dofs[i]),
                         ExcInternalError());
                  dof_indices[i] = dofs[i];
                }
              else
                complete = false;

            if (!complete)
              const_cast<
                typename DoFHandler<dim, spacedim>::active_cell_iterator &>(
                cell)
                ->set_user_flag();
            else
              const_cast<
                typename DoFHandler<dim, spacedim>::active_cell_iterator &>(
                cell)
                ->clear_user_flag();

            const_cast<
              typename DoFHandler<dim, spacedim>::active_cell_iterator &>(cell)
              ->set_dof_indices(dof_indices);
          };

          const auto filter = [](const auto &cell) {
            return cell->user_flag_set();
          };

          GridTools::exchange_cell_data_to_ghosts<
            std::vector<types::global_dof_index>,
            DoFHandler<dim, spacedim>>(dof_handler, pack, unpack, filter);

          // finally update the cell DoF indices caches to make sure
          // our internal data structures are consistent
          update_all_active_cell_dof_indices_caches(dof_handler);


          // have a barrier so that sends between two calls to this
          // function are not mixed up.
          //
          // this is necessary because above we just see if there are
          // messages and then receive them, without discriminating
          // where they come from and whether they were sent in phase
          // 1 or 2 (the function is called twice in a row). the need
          // for a global communication step like this barrier could
          // be avoided by receiving messages specifically from those
          // processors from which we expect messages, and by using
          // different tags for phase 1 and 2, but the cost of a
          // barrier is negligible compared to everything else we do
          // here
          if (const auto *triangulation =
                dynamic_cast<const dealii::parallel::
                               DistributedTriangulationBase<dim, spacedim> *>(
                  &dof_handler.get_triangulation()))
            {
              const int ierr = MPI_Barrier(triangulation->get_communicator());
              AssertThrowMPI(ierr);
            }
          else
            {
              Assert(false,
                     ExcMessage(
                       "The function communicate_dof_indices_on_marked_cells() "
                       "only works with parallel distributed triangulations."));
            }
#  endif
        }



      } // namespace

#endif // DEAL_II_WITH_MPI



      template <int dim, int spacedim>
      ParallelDistributed<dim, spacedim>::ParallelDistributed(
        DoFHandler<dim, spacedim> &dof_handler)
        : dof_handler(&dof_handler)
      {}



      template <int dim, int spacedim>
      NumberCache
      ParallelDistributed<dim, spacedim>::distribute_dofs() const
      {
#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNotImplemented());
        return NumberCache();
#else

        dealii::parallel::DistributedTriangulationBase<dim, spacedim>
          *triangulation =
            (dynamic_cast<
              dealii::parallel::DistributedTriangulationBase<dim, spacedim> *>(
              const_cast<dealii::Triangulation<dim, spacedim> *>(
                &dof_handler->get_triangulation())));
        Assert(triangulation != nullptr, ExcInternalError());

        const types::subdomain_id subdomain_id =
          triangulation->locally_owned_subdomain();


        /*
           The following algorithm has a number of stages that are all
           documented in the paper that describes the parallel::distributed
           functionality:

           1/ locally enumerate dofs on locally owned cells
           2/ eliminate dof duplicates on all cells.
              un-numerate those that are on interfaces with ghost
              cells and that we don't own based on the tie-breaking
              criterion. unify dofs afterwards.
           3/ unify dofs and re-enumerate the remaining valid ones.
              the end result is that we only enumerate locally owned
              DoFs
           4/ shift indices so that each processor has a unique
              range of indices
           5/ for all locally owned cells that are ghost
              cells somewhere else, send our own DoF indices
              to the appropriate set of other processors.
              overwrite invalid DoF indices on ghost interfaces
              with the corresponding valid ones that we now know.
           6/ send DoF indices again to get the correct indices
              on ghost cells that we may not have known earlier
         */

        // --------- Phase 1: enumerate dofs on locally owned cells
        const types::global_dof_index n_initial_local_dofs =
          Implementation::distribute_dofs(subdomain_id, *dof_handler);

        // --------- Phase 2: eliminate dof duplicates on all cells:
        //                    - un-numerate dofs on interfaces to ghost cells
        //                      that we don't own
        //                    - in case of hp-support, unify dofs
        std::vector<dealii::types::global_dof_index> renumbering(
          n_initial_local_dofs, enumeration_dof_index);

        // first, we invalidate degrees of freedom that belong to processors
        // of a lower rank, from which we will receive the final (and lower)
        // degrees of freedom later.
        Implementation::
          invalidate_dof_indices_on_weaker_ghost_cells_for_renumbering(
            renumbering, subdomain_id, *dof_handler);

        // then, we identify DoF duplicates if the DoFHandler has hp-
        // capabilities
        std::vector<std::map<types::global_dof_index, types::global_dof_index>>
          all_constrained_indices(dim);
        Implementation::compute_dof_identities(all_constrained_indices,
                                               *dof_handler);

        // --------- Phase 3: re-enumerate the valid degrees of freedom
        //                    consecutively. thus, we finally receive the
        //                    correct number of locally owned DoFs after
        //                    this step.
        //
        // the order in which we handle Phases 2 and 3 is important,
        // since we want to clarify ownership of degrees of freedom before
        // we actually unify and enumerate their indices. otherwise, we could
        // end up having a degee of freedom to which only invalid indices will
        // be assigned.
        const types::global_dof_index n_locally_owned_dofs =
          Implementation::enumerate_dof_indices_for_renumbering(
            renumbering, all_constrained_indices, *dof_handler);

        // --------- Phase 4: shift indices so that each processor has a unique
        //                    range of indices
        dealii::types::global_dof_index my_shift = 0;
        const int                       ierr =
          MPI_Exscan(DEAL_II_MPI_CONST_CAST(&n_locally_owned_dofs),
                     &my_shift,
                     1,
                     DEAL_II_DOF_INDEX_MPI_TYPE,
                     MPI_SUM,
                     triangulation->get_communicator());
        AssertThrowMPI(ierr);

        // make dof indices globally consecutive
        for (auto &new_index : renumbering)
          if (new_index != numbers::invalid_dof_index)
            new_index += my_shift;

        // now re-enumerate all dofs to this shifted and condensed
        // numbering form.  we renumber some dofs as invalid, so
        // choose the nocheck-version.
        Implementation::renumber_dofs(renumbering,
                                      IndexSet(0),
                                      *dof_handler,
                                      /*check_validity=*/false);

        // now a little bit of housekeeping
        const dealii::types::global_dof_index n_global_dofs =
          Utilities::MPI::sum(n_locally_owned_dofs,
                              triangulation->get_communicator());

        NumberCache number_cache;
        number_cache.n_global_dofs        = n_global_dofs;
        number_cache.n_locally_owned_dofs = n_locally_owned_dofs;
        number_cache.locally_owned_dofs   = IndexSet(n_global_dofs);
        number_cache.locally_owned_dofs.add_range(my_shift,
                                                  my_shift +
                                                    n_locally_owned_dofs);
        number_cache.locally_owned_dofs.compress();

        // this ends the phase where we enumerate degrees of freedom on
        // each processor. what is missing is communicating DoF indices
        // on ghost cells

        // --------- Phase 5: for all locally owned cells that are ghost
        //                    cells somewhere else, send our own DoF indices
        //                    to the appropriate set of other processors
        {
          std::vector<bool> user_flags;
          triangulation->save_user_flags(user_flags);
          triangulation->clear_user_flags();

          // mark all cells that either have to send data (locally
          // owned cells that are adjacent to ghost neighbors in some
          // way) or receive data (all ghost cells) via the user flags
          for (const auto &cell : dof_handler->active_cell_iterators())
            if (cell->is_ghost())
              cell->set_user_flag();

          // Send and receive cells. After this, only the local cells
          // are marked, that received new data. This has to be
          // communicated in a second communication step.
          //
          // as explained in the 'distributed' paper, this has to be
          // done twice
          communicate_dof_indices_on_marked_cells(*dof_handler);

          // If the DoFHandler has hp-capabilities enabled, then we may have
          // received valid indices of degrees of freedom that are dominated
          // by a FE object adjacent to a ghost interface.  thus, we overwrite
          // the remaining invalid indices with the valid ones in this step.
          Implementation::merge_invalid_dof_indices_on_ghost_interfaces(
            *dof_handler);

          // --------- Phase 6: all locally owned cells have their correct
          //                    DoF indices set. however, some ghost cells
          //                    may still have invalid ones. thus, exchange
          //                    one more time.
          communicate_dof_indices_on_marked_cells(*dof_handler);

          // at this point, we must have taken care of the data transfer
          // on all cells we had previously marked. verify this
#  ifdef DEBUG
          for (const auto &cell : dof_handler->active_cell_iterators())
            Assert(cell->user_flag_set() == false, ExcInternalError());
#  endif

          triangulation->load_user_flags(user_flags);
        }

#  ifdef DEBUG
        // check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          for (const auto &cell : dof_handler->active_cell_iterators())
            if (!cell->is_artificial())
              {
                local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
                cell->get_dof_indices(local_dof_indices);
                if (local_dof_indices.end() !=
                    std::find(local_dof_indices.begin(),
                              local_dof_indices.end(),
                              numbers::invalid_dof_index))
                  {
                    if (cell->is_ghost())
                      {
                        Assert(false,
                               ExcMessage(
                                 "A ghost cell ended up with incomplete "
                                 "DoF index information. This should not "
                                 "have happened!"));
                      }
                    else
                      {
                        Assert(
                          false,
                          ExcMessage(
                            "A locally owned cell ended up with incomplete "
                            "DoF index information. This should not "
                            "have happened!"));
                      }
                  }
              }
        }
#  endif // DEBUG
        return number_cache;
#endif   // DEAL_II_WITH_MPI
      }



      template <int dim, int spacedim>
      std::vector<NumberCache>
      ParallelDistributed<dim, spacedim>::distribute_mg_dofs() const
      {
#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNotImplemented());
        return std::vector<NumberCache>();
#else

        dealii::parallel::DistributedTriangulationBase<dim, spacedim>
          *triangulation =
            (dynamic_cast<
              dealii::parallel::DistributedTriangulationBase<dim, spacedim> *>(
              const_cast<dealii::Triangulation<dim, spacedim> *>(
                &dof_handler->get_triangulation())));
        Assert(triangulation != nullptr, ExcInternalError());

        AssertThrow((triangulation->is_multilevel_hierarchy_constructed()),
                    ExcMessage(
                      "Multigrid DoFs can only be distributed on a parallel "
                      "Triangulation if the flag construct_multigrid_hierarchy "
                      "is set in the constructor."));

        // loop over all levels that exist globally (across all
        // processors), even if the current processor does not in fact
        // have any cells on that level or if the local part of the
        // Triangulation has fewer levels. we need to do this because
        // we need to communicate across all processors on all levels
        const unsigned int       n_levels = triangulation->n_global_levels();
        std::vector<NumberCache> number_caches;
        number_caches.reserve(n_levels);
        for (unsigned int level = 0; level < n_levels; ++level)
          {
            NumberCache level_number_cache;

            //* 1. distribute on own subdomain
            const unsigned int n_initial_local_dofs =
              Implementation::distribute_dofs_on_level(
                triangulation->locally_owned_subdomain(), *dof_handler, level);

            //* 2. iterate over ghostcells and kill dofs that are not
            // owned by us
            std::vector<dealii::types::global_dof_index> renumbering(
              n_initial_local_dofs);
            for (dealii::types::global_dof_index i = 0; i < renumbering.size();
                 ++i)
              renumbering[i] = i;

            if (level < triangulation->n_levels())
              {
                std::vector<dealii::types::global_dof_index> local_dof_indices;

                for (const auto &cell :
                     dof_handler->cell_iterators_on_level(level))
                  if (cell->level_subdomain_id() !=
                        numbers::artificial_subdomain_id &&
                      (cell->level_subdomain_id() <
                       triangulation->locally_owned_subdomain()))
                    {
                      // we found a neighboring ghost cell whose
                      // subdomain is "stronger" than our own
                      // subdomain

                      // delete all dofs that live there and that we
                      // have previously assigned a number to
                      // (i.e. the ones on the interface)
                      local_dof_indices.resize(
                        cell->get_fe().n_dofs_per_cell());
                      cell->get_mg_dof_indices(local_dof_indices);
                      for (unsigned int i = 0;
                           i < cell->get_fe().n_dofs_per_cell();
                           ++i)
                        if (local_dof_indices[i] != numbers::invalid_dof_index)
                          renumbering[local_dof_indices[i]] =
                            numbers::invalid_dof_index;
                    }
              }

            level_number_cache.n_locally_owned_dofs = 0;
            for (types::global_dof_index &index : renumbering)
              if (index != numbers::invalid_dof_index)
                index = level_number_cache.n_locally_owned_dofs++;

            //* 3. communicate local dofcount and shift ids to make
            // them unique
            dealii::types::global_dof_index my_shift = 0;
            int ierr = MPI_Exscan(DEAL_II_MPI_CONST_CAST(
                                    &level_number_cache.n_locally_owned_dofs),
                                  &my_shift,
                                  1,
                                  DEAL_II_DOF_INDEX_MPI_TYPE,
                                  MPI_SUM,
                                  triangulation->get_communicator());
            AssertThrowMPI(ierr);

            // The last processor knows about the total number of dofs, so we
            // can use a cheaper broadcast rather than an MPI_Allreduce via
            // MPI::sum().
            level_number_cache.n_global_dofs =
              my_shift + level_number_cache.n_locally_owned_dofs;
            ierr = MPI_Bcast(&level_number_cache.n_global_dofs,
                             1,
                             DEAL_II_DOF_INDEX_MPI_TYPE,
                             Utilities::MPI::n_mpi_processes(
                               triangulation->get_communicator()) -
                               1,
                             triangulation->get_communicator());
            AssertThrowMPI(ierr);

            // shift indices
            for (types::global_dof_index &index : renumbering)
              if (index != numbers::invalid_dof_index)
                index += my_shift;

            // now re-enumerate all dofs to this shifted and condensed
            // numbering form.  we renumber some dofs as invalid, so
            // choose the nocheck-version of the function
            //
            // of course there is nothing for us to renumber if the
            // level we are currently dealing with doesn't even exist
            // within the current triangulation, so skip renumbering
            // in that case
            if (level < triangulation->n_levels())
              Implementation::renumber_mg_dofs(
                renumbering, IndexSet(0), *dof_handler, level, false);

            // now a little bit of housekeeping
            level_number_cache.locally_owned_dofs =
              IndexSet(level_number_cache.n_global_dofs);
            level_number_cache.locally_owned_dofs.add_range(
              my_shift, my_shift + level_number_cache.n_locally_owned_dofs);
            level_number_cache.locally_owned_dofs.compress();

            number_caches.emplace_back(level_number_cache);
          }


        //* communicate ghost DoFs
        // We mark all ghost cells by setting the user_flag and then request
        // these cells from the corresponding owners. As this information
        // can be incomplete,
        {
          std::vector<bool> user_flags;
          triangulation->save_user_flags(user_flags);
          triangulation->clear_user_flags();

          // mark all ghost cells for transfer
          {
            for (const auto &cell : dof_handler->cell_iterators())
              if (cell->level_subdomain_id() !=
                    dealii::numbers::artificial_subdomain_id &&
                  !cell->is_locally_owned_on_level())
                cell->set_user_flag();
          }

          // Phase 1. Request all marked cells from corresponding owners. If we
          // managed to get every DoF, remove the user_flag, otherwise we
          // will request them again in the step below.
          communicate_mg_ghost_cells(*triangulation, *dof_handler);

          // have a barrier so that sends from above and below this
          // place are not mixed up.
          //
          // this is necessary because above we just see if there are
          // messages and then receive them, without discriminating
          // where they come from and whether they were sent in phase
          // 1 or 2 in communicate_mg_ghost_cells() on another
          // processor. the need for a global communication step like
          // this barrier could be avoided by receiving messages
          // specifically from those processors from which we expect
          // messages, and by using different tags for phase 1 and 2,
          // but the cost of a barrier is negligible compared to
          // everything else we do here
          const int ierr = MPI_Barrier(triangulation->get_communicator());
          AssertThrowMPI(ierr);

          // Phase 2, only request the cells that were not completed
          // in Phase 1.
          communicate_mg_ghost_cells(*triangulation, *dof_handler);

#  ifdef DEBUG
          // make sure we have removed all flags:
          {
            for (const auto &cell : dof_handler->cell_iterators())
              if (cell->level_subdomain_id() !=
                    dealii::numbers::artificial_subdomain_id &&
                  !cell->is_locally_owned_on_level())
                Assert(cell->user_flag_set() == false, ExcInternalError());
          }
#  endif

          triangulation->load_user_flags(user_flags);
        }



#  ifdef DEBUG
        // check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;
          for (const auto &cell : dof_handler->cell_iterators())
            if (cell->level_subdomain_id() !=
                dealii::numbers::artificial_subdomain_id)
              {
                local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
                cell->get_mg_dof_indices(local_dof_indices);
                if (local_dof_indices.end() !=
                    std::find(local_dof_indices.begin(),
                              local_dof_indices.end(),
                              numbers::invalid_dof_index))
                  {
                    Assert(false, ExcMessage("not all DoFs got distributed!"));
                  }
              }
        }
#  endif // DEBUG

        return number_caches;

#endif // DEAL_II_WITH_MPI
      }


      template <int dim, int spacedim>
      NumberCache
      ParallelDistributed<dim, spacedim>::renumber_dofs(
        const std::vector<dealii::types::global_dof_index> &new_numbers) const
      {
        (void)new_numbers;

        Assert(new_numbers.size() == dof_handler->n_locally_owned_dofs(),
               ExcInternalError());

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNotImplemented());
        return NumberCache();
#else

        dealii::parallel::DistributedTriangulationBase<dim, spacedim>
          *triangulation =
            (dynamic_cast<
              dealii::parallel::DistributedTriangulationBase<dim, spacedim> *>(
              const_cast<dealii::Triangulation<dim, spacedim> *>(
                &dof_handler->get_triangulation())));
        Assert(triangulation != nullptr, ExcInternalError());


        // We start by checking whether only the numbering within the MPI
        // ranks changed. In that case, we can apply the renumbering with some
        // local renumbering only (this is similar to the renumber_mg_dofs()
        // function below)
        const bool locally_owned_set_changes =
          std::any_of(new_numbers.cbegin(),
                      new_numbers.cend(),
                      [this](const types::global_dof_index i) {
                        return dof_handler->locally_owned_dofs().is_element(
                                 i) == false;
                      });

        if (Utilities::MPI::sum(static_cast<unsigned int>(
                                  locally_owned_set_changes),
                                triangulation->get_communicator()) == 0)
          {
            // Since only the order within the local subdomains has changed,
            // all we need to do is to propagate the knowledge about the
            // numbers from the locally owned dofs (given by the new_numbers
            // array) to all ghosted dofs on neighboring processors. We can do
            // this by ghost layer exchange routines as in parallel vectors:
            // We create an IndexSet for the relevant dofs and then export
            // into an array of those values via Utilities::MPI::Partitioner.
            IndexSet relevant_dofs;
            DoFTools::extract_locally_relevant_dofs(*dof_handler,
                                                    relevant_dofs);
            std::vector<types::global_dof_index> ghosted_new_numbers(
              relevant_dofs.n_elements());
            {
              Utilities::MPI::Partitioner partitioner(
                dof_handler->locally_owned_dofs(),
                relevant_dofs,
                triangulation->get_communicator());

              // choose some number that makes it unlikely to get conflicts
              // with other ongoing non-blocking communication (there
              // shouldn't be any at this place in most programs).
              const unsigned int                   communication_channel = 19;
              std::vector<types::global_dof_index> temp_array(
                partitioner.n_import_indices());
              std::vector<MPI_Request> requests;
              partitioner.export_to_ghosted_array_start(
                communication_channel,
                make_array_view(new_numbers),
                make_array_view(temp_array),
                ArrayView<types::global_dof_index>(
                  ghosted_new_numbers.data() + new_numbers.size(),
                  partitioner.n_ghost_indices()),
                requests);
              partitioner.export_to_ghosted_array_finish(
                ArrayView<types::global_dof_index>(
                  ghosted_new_numbers.data() + new_numbers.size(),
                  partitioner.n_ghost_indices()),
                requests);

              // we need to fill the indices of the locally owned part into
              // the new numbers array, which is not provided by the parallel
              // partitioner. their right position is somewhere in the middle
              // of the array, so we first copy the ghosted part from smaller
              // ranks to the front, then insert the data in the middle.
              unsigned int n_ghosts_on_smaller_ranks = 0;
              for (std::pair<unsigned int, unsigned int> t :
                   partitioner.ghost_targets())
                {
                  if (t.first > partitioner.this_mpi_process())
                    break;
                  n_ghosts_on_smaller_ranks += t.second;
                }
              if (n_ghosts_on_smaller_ranks > 0)
                {
                  Assert(ghosted_new_numbers.data() != nullptr,
                         ExcInternalError());
                  std::memmove(ghosted_new_numbers.data(),
                               ghosted_new_numbers.data() + new_numbers.size(),
                               sizeof(types::global_dof_index) *
                                 n_ghosts_on_smaller_ranks);
                }
              if (new_numbers.size() > 0)
                {
                  Assert(new_numbers.data() != nullptr, ExcInternalError());
                  std::memcpy(ghosted_new_numbers.data() +
                                n_ghosts_on_smaller_ranks,
                              new_numbers.data(),
                              sizeof(types::global_dof_index) *
                                new_numbers.size());
                }
            }

            // In case we do not carry any relevant dof (but only some remote
            // processor), we do not need to call the renumbering. We call the
            // version without validity check because vertex dofs will be
            // set already in the artificial region.
            if (relevant_dofs.n_elements() > 0)
              Implementation::renumber_dofs(ghosted_new_numbers,
                                            relevant_dofs,
                                            *dof_handler,
                                            /*check_validity=*/false);

            NumberCache number_cache;
            number_cache.locally_owned_dofs = dof_handler->locally_owned_dofs();
            number_cache.n_global_dofs      = dof_handler->n_dofs();
            number_cache.n_locally_owned_dofs =
              number_cache.locally_owned_dofs.n_elements();
            return number_cache;
          }
        else
          {
            // Now back to the more complicated case
            //
            // First figure out the new set of locally owned DoF indices.
            // If we own no DoFs, we still need to go through this function,
            // but we can skip this calculation.
            //
            // The IndexSet::add_indices() function is substantially more
            // efficient if the set of indices is already sorted because
            // it can then insert ranges instead of individual elements.
            // consequently, pre-sort the array of new indices
            IndexSet my_locally_owned_new_dof_indices(dof_handler->n_dofs());
            if (dof_handler->n_locally_owned_dofs() > 0)
              {
                std::vector<dealii::types::global_dof_index>
                  new_numbers_sorted = new_numbers;
                std::sort(new_numbers_sorted.begin(), new_numbers_sorted.end());

                my_locally_owned_new_dof_indices.add_indices(
                  new_numbers_sorted.begin(), new_numbers_sorted.end());
                my_locally_owned_new_dof_indices.compress();

                Assert(my_locally_owned_new_dof_indices.n_elements() ==
                         new_numbers.size(),
                       ExcInternalError());
              }

            // delete all knowledge of DoF indices that are not locally
            // owned. we do so by getting DoF indices on cells, checking
            // whether they are locally owned, if not, setting them to
            // an invalid value, and then setting them again on the current
            // cell
            //
            // DoFs we (i) know about, and (ii) don't own locally must be
            // located either on ghost cells, or on the interface between a
            // locally owned cell and a ghost cell. In any case, it is
            // sufficient to kill them only from the ghost side cell, so loop
            // only over ghost cells
            {
              std::vector<dealii::types::global_dof_index> local_dof_indices;

              for (auto cell : dof_handler->active_cell_iterators())
                if (cell->is_ghost())
                  {
                    local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
                    cell->get_dof_indices(local_dof_indices);

                    for (unsigned int i = 0;
                         i < cell->get_fe().n_dofs_per_cell();
                         ++i)
                      // delete a DoF index if it has not already been deleted
                      // (e.g., by visiting a neighboring cell, if it is on the
                      // boundary), and if we don't own it
                      if ((local_dof_indices[i] !=
                           numbers::invalid_dof_index) &&
                          (!dof_handler->locally_owned_dofs().is_element(
                            local_dof_indices[i])))
                        local_dof_indices[i] = numbers::invalid_dof_index;

                    cell->set_dof_indices(local_dof_indices);
                  }
            }


            // renumber. Skip when there is nothing to do because we own no DoF.
            if (dof_handler->locally_owned_dofs().n_elements() > 0)
              Implementation::renumber_dofs(new_numbers,
                                            dof_handler->locally_owned_dofs(),
                                            *dof_handler,
                                            /*check_validity=*/false);

            // Communicate newly assigned DoF indices to other processors
            // and get the same information for our own ghost cells.
            //
            // This is the same as phase 5+6 in the distribute_dofs() algorithm,
            // taking into account that we have to unify a few DoFs in between
            // then communication phases if we do hp-numbering
            {
              std::vector<bool> user_flags;
              triangulation->save_user_flags(user_flags);
              triangulation->clear_user_flags();

              // mark all own cells for transfer
              for (const auto &cell : dof_handler->active_cell_iterators())
                if (cell->is_ghost())
                  cell->set_user_flag();


              // Send and receive cells. After this, only the local cells
              // are marked, that received new data. This has to be
              // communicated in a second communication step.
              //
              // as explained in the 'distributed' paper, this has to be
              // done twice
              communicate_dof_indices_on_marked_cells(*dof_handler);

              // if the DoFHandler has hp-capabilities then we may have
              // received valid indices of degrees of freedom that are
              // dominated by a FE object adjacent to a ghost interface.
              // thus, we overwrite the remaining invalid indices with the
              // valid ones in this step.
              Implementation::merge_invalid_dof_indices_on_ghost_interfaces(
                *dof_handler);

              communicate_dof_indices_on_marked_cells(*dof_handler);

              triangulation->load_user_flags(user_flags);
            }

            NumberCache number_cache;
            number_cache.locally_owned_dofs = my_locally_owned_new_dof_indices;
            number_cache.n_global_dofs      = dof_handler->n_dofs();
            number_cache.n_locally_owned_dofs =
              number_cache.locally_owned_dofs.n_elements();
            return number_cache;
          }
#endif
      }



      template <int dim, int spacedim>
      NumberCache
      ParallelDistributed<dim, spacedim>::renumber_mg_dofs(
        const unsigned int                          level,
        const std::vector<types::global_dof_index> &new_numbers) const
      {
        // we only implement the case where the multigrid numbers are
        // renumbered within the processor's partition, rather than the most
        // general case
        const IndexSet index_set = dof_handler->locally_owned_mg_dofs(level);

#ifdef DEAL_II_WITH_MPI

        const dealii::parallel::TriangulationBase<dim, spacedim> *tr =
          (dynamic_cast<const dealii::parallel::TriangulationBase<dim, spacedim>
                          *>(&this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        const unsigned int my_rank =
          Utilities::MPI::this_mpi_process(tr->get_communicator());

#  ifdef DEBUG
        for (types::global_dof_index i : new_numbers)
          {
            Assert(index_set.is_element(i),
                   ExcNotImplemented(
                     "Renumberings that change the locally owned mg dofs "
                     "partitioning are currently not implemented for "
                     "the multigrid levels"));
          }
#  endif

        // we need to access all locally relevant degrees of freedom. we
        // use Utilities::MPI::Partitioner for handling the data exchange
        // of the new numbers, which is simply the extraction of ghost data
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(*dof_handler,
                                                      level,
                                                      relevant_dofs);
        std::vector<types::global_dof_index> ghosted_new_numbers(
          relevant_dofs.n_elements());
        {
          Utilities::MPI::Partitioner          partitioner(index_set,
                                                  relevant_dofs,
                                                  tr->get_communicator());
          std::vector<types::global_dof_index> temp_array(
            partitioner.n_import_indices());
          const unsigned int       communication_channel = 17;
          std::vector<MPI_Request> requests;
          partitioner.export_to_ghosted_array_start(
            communication_channel,
            make_array_view(new_numbers),
            make_array_view(temp_array),
            ArrayView<types::global_dof_index>(ghosted_new_numbers.data() +
                                                 new_numbers.size(),
                                               partitioner.n_ghost_indices()),
            requests);
          partitioner.export_to_ghosted_array_finish(
            ArrayView<types::global_dof_index>(ghosted_new_numbers.data() +
                                                 new_numbers.size(),
                                               partitioner.n_ghost_indices()),
            requests);

          // we need to fill the indices of the locally owned part into the
          // new numbers array. their right position is somewhere in the
          // middle of the array, so we first copy the ghosted part from
          // smaller ranks to the front, then insert the data in the middle.
          unsigned int n_ghosts_on_smaller_ranks = 0;
          for (std::pair<unsigned int, unsigned int> t :
               partitioner.ghost_targets())
            {
              if (t.first > my_rank)
                break;
              n_ghosts_on_smaller_ranks += t.second;
            }
          if (n_ghosts_on_smaller_ranks > 0)
            {
              Assert(ghosted_new_numbers.data() != nullptr, ExcInternalError());
              std::memmove(ghosted_new_numbers.data(),
                           ghosted_new_numbers.data() + new_numbers.size(),
                           sizeof(types::global_dof_index) *
                             n_ghosts_on_smaller_ranks);
            }
          if (new_numbers.size() > 0)
            {
              Assert(new_numbers.data() != nullptr, ExcInternalError());
              std::memcpy(ghosted_new_numbers.data() +
                            n_ghosts_on_smaller_ranks,
                          new_numbers.data(),
                          sizeof(types::global_dof_index) * new_numbers.size());
            }
        }

        // in case we do not own any of the given level (but only some remote
        // processor), we do not need to call the renumbering
        if (level < this->dof_handler->get_triangulation().n_levels() &&
            relevant_dofs.n_elements() > 0)
          Implementation::renumber_mg_dofs(
            ghosted_new_numbers, relevant_dofs, *dof_handler, level, true);
#else
        (void)new_numbers;
        Assert(false, ExcNotImplemented());
#endif

        NumberCache number_cache;
        number_cache.locally_owned_dofs = index_set;
        number_cache.n_global_dofs      = dof_handler->n_dofs(level);
        number_cache.n_locally_owned_dofs =
          number_cache.locally_owned_dofs.n_elements();
        return number_cache;
      }
    } // namespace Policy
  }   // namespace DoFHandlerImplementation
} // namespace internal



/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler_policy.inst"


DEAL_II_NAMESPACE_CLOSE
