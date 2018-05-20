// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/device/back_inserter.hpp>
#  include <boost/iostreams/filter/gzip.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#  include <boost/iostreams/stream.hpp>
#  include <boost/serialization/array.hpp>
#endif

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

      namespace hp
      {
        using dealii::hp::DoFHandler;
      }

      namespace
      {
        /**
         * Update the cache used for cell dof indices on all (non-artificial)
         * active cells of the given DoFHandler.
         */
        template <class DoFHandlerType>
        void
        update_all_active_cell_dof_indices_caches(
          const DoFHandlerType& dof_handler)
        {
          typename DoFHandlerType::active_cell_iterator beginc
            = dof_handler.begin_active(),
            endc = dof_handler.end();

          auto worker
            = [](const typename DoFHandlerType::active_cell_iterator& cell,
                 void*,
                 void*) {
                if(!cell->is_artificial())
                  cell->update_cell_dof_indices_cache();
              };

          // parallelize filling all of the cell caches. by using
          // WorkStream, we make sure that we only run through the
          // range of iterators once, whereas a parallel_for loop
          // for example has to split the range multiple times,
          // which is expensive because cell iterators are not
          // random access iterators with a cheap operator-
          WorkStream::run(beginc,
                          endc,
                          worker,
                          /* copier */ std::function<void(void*)>(),
                          /* scratch_data */ nullptr,
                          /* copy_data */ nullptr,
                          2 * MultithreadInfo::n_threads(),
                          /* chunk_size = */ 32);
        }

        /**
         * Update the cache used for cell dof indices on all (non-artificial)
         * level (multigrid) cells of the given DoFHandler.
         */
        template <class DoFHandlerType>
        void
        update_all_level_cell_dof_indices_caches(
          const DoFHandlerType& dof_handler)
        {
          typename DoFHandlerType::level_cell_iterator beginc
            = dof_handler.begin(),
            endc = dof_handler.end();

          auto worker
            = [](const typename DoFHandlerType::level_cell_iterator& cell,
                 void*,
                 void*) {
                if(cell->has_children() || !cell->is_artificial())
                  cell->update_cell_dof_indices_cache();
              };

          // parallelize filling all of the cell caches. by using
          // WorkStream, we make sure that we only run through the
          // range of iterators once, whereas a parallel_for loop
          // for example has to split the range multiple times,
          // which is expensive because cell iterators are not
          // random access iterators with a cheap operator-
          WorkStream::run(beginc,
                          endc,
                          worker,
                          /* copier */ std::function<void(void*)>(),
                          /* scratch_data */ nullptr,
                          /* copy_data */ nullptr,
                          2 * MultithreadInfo::n_threads(),
                          /* chunk_size = */ 32);
        }

        typedef std::vector<std::pair<unsigned int, unsigned int>>
          DoFIdentities;

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
        void
        ensure_existence_of_dof_identities(
          const FiniteElement<dim, spacedim>& fe1,
          const FiniteElement<dim, spacedim>& fe2,
          std::unique_ptr<DoFIdentities>&     identities)
        {
          // see if we need to fill this entry, or whether it already
          // exists
          if(identities.get() == nullptr)
            {
              switch(structdim)
                {
                  case 0:
                    {
                      identities = std_cxx14::make_unique<DoFIdentities>(
                        fe1.hp_vertex_dof_identities(fe2));
                      break;
                    }

                  case 1:
                    {
                      identities = std_cxx14::make_unique<DoFIdentities>(
                        fe1.hp_line_dof_identities(fe2));
                      break;
                    }

                  case 2:
                    {
                      identities = std_cxx14::make_unique<DoFIdentities>(
                        fe1.hp_quad_dof_identities(fe2));
                      break;
                    }

                  default:
                    Assert(false, ExcNotImplemented());
                }

              // double check whether the newly created entries make
              // any sense at all
              for(unsigned int i = 0; i < identities->size(); ++i)
                {
                  Assert((*identities)[i].first
                           < fe1.template n_dofs_per_object<structdim>(),
                         ExcInternalError());
                  Assert((*identities)[i].second
                           < fe2.template n_dofs_per_object<structdim>(),
                         ExcInternalError());
                }
            }
        }

        /**
         * For an object, such as a line or a quad iterator, determine
         * the fe_index of the most dominating finite element that
         * lives on this object.
         *
         * Return numbers::invalid_unsigned_int if we couldn't find one.
         */
        template <int dim, int spacedim, typename iterator>
        unsigned int
        get_most_dominating_fe_index(const iterator& object)
        {
          unsigned int dominating_fe_index = 0;
          for(; dominating_fe_index < object->n_active_fe_indices();
              ++dominating_fe_index)
            {
              const FiniteElement<dim, spacedim>& this_fe = object->get_fe(
                object->nth_active_fe_index(dominating_fe_index));

              FiniteElementDomination::Domination domination
                = FiniteElementDomination::either_element_can_dominate;
              for(unsigned int other_fe_index = 0;
                  other_fe_index < object->n_active_fe_indices();
                  ++other_fe_index)
                if(other_fe_index != dominating_fe_index)
                  {
                    const FiniteElement<dim, spacedim>& that_fe
                      = object->get_fe(
                        object->nth_active_fe_index(other_fe_index));

                    domination = domination
                                 & this_fe.compare_for_face_domination(that_fe);
                  }

              // see if this element is able to dominate all the other
              // ones, and if so take it
              if((domination == FiniteElementDomination::this_element_dominates)
                 || (domination
                     == FiniteElementDomination::either_element_can_dominate)
                 || (domination == FiniteElementDomination::no_requirements))
                break;
            }

          // check that we have found one such fe
          if(dominating_fe_index != object->n_active_fe_indices())
            {
              // return the finite element index used on it. note that
              // only a single fe can be active on such subfaces
              return object->nth_active_fe_index(dominating_fe_index);
            }
          else
            {
              // if we couldn't find the most dominating object
              return numbers::invalid_unsigned_int;
            }
        }
      } // namespace

      struct Implementation
      {
        /* -------------- distribute_dofs functionality ------------- */

        /**
         * Distribute dofs on the given cell, with new dofs starting with index
         * @p next_free_dof. Return the next unused index number.
         *
         * This function is refactored from the main @p distribute_dofs function since
         * it can not be implemented dimension independent.
         */
        template <int spacedim>
        static types::global_dof_index
        distribute_dofs_on_cell(
          const DoFHandler<1, spacedim>& dof_handler,
          const typename DoFHandler<1, spacedim>::active_cell_iterator& cell,
          types::global_dof_index next_free_dof)
        {
          // distribute dofs of vertices
          if(dof_handler.get_fe().dofs_per_vertex > 0)
            for(unsigned int v = 0; v < GeometryInfo<1>::vertices_per_cell; ++v)
              {
                if(cell->vertex_dof_index(v, 0) == numbers::invalid_dof_index)
                  for(unsigned int d = 0;
                      d < dof_handler.get_fe().dofs_per_vertex;
                      ++d)
                    {
                      Assert((cell->vertex_dof_index(v, d)
                              == numbers::invalid_dof_index),
                             ExcInternalError());
                      cell->set_vertex_dof_index(v, d, next_free_dof++);
                    }
                else
                  for(unsigned int d = 0;
                      d < dof_handler.get_fe().dofs_per_vertex;
                      ++d)
                    Assert((cell->vertex_dof_index(v, d)
                            != numbers::invalid_dof_index),
                           ExcInternalError());
              }

          // dofs of line
          for(unsigned int d = 0; d < dof_handler.get_fe().dofs_per_line; ++d)
            cell->set_dof_index(d, next_free_dof++);

          return next_free_dof;
        }

        template <int spacedim>
        static types::global_dof_index
        distribute_dofs_on_cell(
          const DoFHandler<2, spacedim>& dof_handler,
          const typename DoFHandler<2, spacedim>::active_cell_iterator& cell,
          types::global_dof_index next_free_dof)
        {
          if(dof_handler.get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<2>::vertices_per_cell;
                ++vertex)
              // check whether dofs for this vertex have been distributed
              // (checking the first dof should be good enough)
              if(cell->vertex_dof_index(vertex, 0)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0;
                    d < dof_handler.get_fe().dofs_per_vertex;
                    ++d)
                  cell->set_vertex_dof_index(vertex, d, next_free_dof++);

          // for the four sides
          if(dof_handler.get_fe().dofs_per_line > 0)
            for(unsigned int side = 0; side < GeometryInfo<2>::faces_per_cell;
                ++side)
              {
                const typename DoFHandler<2, spacedim>::line_iterator line
                  = cell->line(side);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (checking the first dof should be good enough)
                if(line->dof_index(0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for(unsigned int d = 0;
                      d < dof_handler.get_fe().dofs_per_line;
                      ++d)
                    line->set_dof_index(d, next_free_dof++);
              }

          // dofs of quad
          if(dof_handler.get_fe().dofs_per_quad > 0)
            for(unsigned int d = 0; d < dof_handler.get_fe().dofs_per_quad; ++d)
              cell->set_dof_index(d, next_free_dof++);

          return next_free_dof;
        }

        template <int spacedim>
        static types::global_dof_index
        distribute_dofs_on_cell(
          const DoFHandler<3, spacedim>& dof_handler,
          const typename DoFHandler<3, spacedim>::active_cell_iterator& cell,
          types::global_dof_index next_free_dof)
        {
          if(dof_handler.get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<3>::vertices_per_cell;
                ++vertex)
              // check whether dofs for this vertex have been distributed
              // (checking the first dof should be good enough)
              if(cell->vertex_dof_index(vertex, 0)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0;
                    d < dof_handler.get_fe().dofs_per_vertex;
                    ++d)
                  cell->set_vertex_dof_index(vertex, d, next_free_dof++);

          // for the lines
          if(dof_handler.get_fe().dofs_per_line > 0)
            for(unsigned int l = 0; l < GeometryInfo<3>::lines_per_cell; ++l)
              {
                const typename DoFHandler<3, spacedim>::line_iterator line
                  = cell->line(l);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (checking the first dof should be good enough)
                if(line->dof_index(0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for(unsigned int d = 0;
                      d < dof_handler.get_fe().dofs_per_line;
                      ++d)
                    line->set_dof_index(d, next_free_dof++);
              }

          // for the quads
          if(dof_handler.get_fe().dofs_per_quad > 0)
            for(unsigned int q = 0; q < GeometryInfo<3>::quads_per_cell; ++q)
              {
                const typename DoFHandler<3, spacedim>::quad_iterator quad
                  = cell->quad(q);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (checking the first dof should be good enough)
                if(quad->dof_index(0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for(unsigned int d = 0;
                      d < dof_handler.get_fe().dofs_per_quad;
                      ++d)
                    quad->set_dof_index(d, next_free_dof++);
              }

          // dofs of hex
          if(dof_handler.get_fe().dofs_per_hex > 0)
            for(unsigned int d = 0; d < dof_handler.get_fe().dofs_per_hex; ++d)
              cell->set_dof_index(d, next_free_dof++);

          return next_free_dof;
        }

        // same for the hp::DoFHandler
        template <int spacedim>
        static types::global_dof_index
        distribute_dofs_on_cell(
          const hp::DoFHandler<1, spacedim>&,
          const typename hp::DoFHandler<1, spacedim>::active_cell_iterator&
                                  cell,
          types::global_dof_index next_free_dof)
        {
          const unsigned int dim = 1;

          const FiniteElement<dim, spacedim>& fe = cell->get_fe();
          const unsigned int fe_index            = cell->active_fe_index();

          // number dofs on vertices. to do so, check whether dofs for
          // this vertex have been distributed and for the present fe
          // (only check the first dof), and if this isn't the case
          // distribute new ones there
          if(fe.dofs_per_vertex > 0)
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<1>::vertices_per_cell;
                ++vertex)
              if(cell->vertex_dof_index(vertex, 0, fe_index)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0; d < fe.dofs_per_vertex;
                    ++d, ++next_free_dof)
                  cell->set_vertex_dof_index(
                    vertex, d, next_free_dof, fe_index);

          // finally for the line. this one shouldn't be numbered yet
          if(fe.dofs_per_line > 0)
            {
              Assert(
                (cell->dof_index(0, fe_index) == numbers::invalid_dof_index),
                ExcInternalError());

              for(unsigned int d = 0; d < fe.dofs_per_line;
                  ++d, ++next_free_dof)
                cell->set_dof_index(d, next_free_dof, fe_index);
            }

          // note that this cell has been processed
          cell->set_user_flag();

          return next_free_dof;
        }

        template <int spacedim>
        static types::global_dof_index
        distribute_dofs_on_cell(
          const hp::DoFHandler<2, spacedim>&,
          const typename hp::DoFHandler<2, spacedim>::active_cell_iterator&
                                  cell,
          types::global_dof_index next_free_dof)
        {
          const unsigned int dim = 2;

          const FiniteElement<dim, spacedim>& fe = cell->get_fe();
          const unsigned int fe_index            = cell->active_fe_index();

          // number dofs on vertices. to do so, check whether dofs for
          // this vertex have been distributed and for the present fe
          // (only check the first dof), and if this isn't the case
          // distribute new ones there
          if(fe.dofs_per_vertex > 0)
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<2>::vertices_per_cell;
                ++vertex)
              if(cell->vertex_dof_index(vertex, 0, fe_index)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0; d < fe.dofs_per_vertex;
                    ++d, ++next_free_dof)
                  cell->set_vertex_dof_index(
                    vertex, d, next_free_dof, fe_index);

          // next the sides. do the same as above: check whether the
          // line is already numbered for the present fe_index, and if
          // not do it
          if(fe.dofs_per_line > 0)
            for(unsigned int l = 0; l < GeometryInfo<2>::lines_per_cell; ++l)
              {
                typename hp::DoFHandler<dim, spacedim>::line_iterator line
                  = cell->line(l);

                if(line->dof_index(0, fe_index) == numbers::invalid_dof_index)
                  for(unsigned int d = 0; d < fe.dofs_per_line;
                      ++d, ++next_free_dof)
                    line->set_dof_index(d, next_free_dof, fe_index);
              }

          // finally for the quad. this one shouldn't be numbered yet
          if(fe.dofs_per_quad > 0)
            {
              Assert(
                (cell->dof_index(0, fe_index) == numbers::invalid_dof_index),
                ExcInternalError());

              for(unsigned int d = 0; d < fe.dofs_per_quad;
                  ++d, ++next_free_dof)
                cell->set_dof_index(d, next_free_dof, fe_index);
            }

          // note that this cell has been processed
          cell->set_user_flag();

          return next_free_dof;
        }

        template <int spacedim>
        static types::global_dof_index
        distribute_dofs_on_cell(
          const hp::DoFHandler<3, spacedim>&,
          const typename hp::DoFHandler<3, spacedim>::active_cell_iterator&
                                  cell,
          types::global_dof_index next_free_dof)
        {
          const unsigned int dim = 3;

          const FiniteElement<dim, spacedim>& fe = cell->get_fe();
          const unsigned int fe_index            = cell->active_fe_index();

          // number dofs on vertices. to do so, check whether dofs for
          // this vertex have been distributed and for the present fe
          // (only check the first dof), and if this isn't the case
          // distribute new ones there
          if(fe.dofs_per_vertex > 0)
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<3>::vertices_per_cell;
                ++vertex)
              if(cell->vertex_dof_index(vertex, 0, fe_index)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0; d < fe.dofs_per_vertex;
                    ++d, ++next_free_dof)
                  cell->set_vertex_dof_index(
                    vertex, d, next_free_dof, fe_index);

          // next the four lines. do the same as above: check whether
          // the line is already numbered for the present fe_index,
          // and if not do it
          if(fe.dofs_per_line > 0)
            for(unsigned int l = 0; l < GeometryInfo<3>::lines_per_cell; ++l)
              {
                typename hp::DoFHandler<dim, spacedim>::line_iterator line
                  = cell->line(l);

                if(line->dof_index(0, fe_index) == numbers::invalid_dof_index)
                  for(unsigned int d = 0; d < fe.dofs_per_line;
                      ++d, ++next_free_dof)
                    line->set_dof_index(d, next_free_dof, fe_index);
              }

          // same for quads
          if(fe.dofs_per_quad > 0)
            for(unsigned int q = 0; q < GeometryInfo<3>::quads_per_cell; ++q)
              {
                typename hp::DoFHandler<dim, spacedim>::quad_iterator quad
                  = cell->quad(q);

                if(quad->dof_index(0, fe_index) == numbers::invalid_dof_index)
                  for(unsigned int d = 0; d < fe.dofs_per_quad;
                      ++d, ++next_free_dof)
                    quad->set_dof_index(d, next_free_dof, fe_index);
              }

          // finally for the hex. this one shouldn't be numbered yet
          if(fe.dofs_per_hex > 0)
            {
              Assert(
                (cell->dof_index(0, fe_index) == numbers::invalid_dof_index),
                ExcInternalError());

              for(unsigned int d = 0; d < fe.dofs_per_hex; ++d, ++next_free_dof)
                cell->set_dof_index(d, next_free_dof, fe_index);
            }

          // note that this cell has been processed
          cell->set_user_flag();

          return next_free_dof;
        }

        /**
         * Compute identities between DoFs located on vertices. Called from
         * distribute_dofs().
         */
        template <int dim, int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_vertex_dof_identities(
          hp::DoFHandler<dim, spacedim>& dof_handler)
        {
          std::map<types::global_dof_index, types::global_dof_index>
            dof_identities;

          // Note: we may wish to have something here similar to what
          // we do for lines and quads, namely that we only identify
          // dofs for any fe towards the most dominating one. however,
          // it is not clear whether this is actually necessary for
          // vertices at all, I can't think of a finite element that
          // would make that necessary...
          dealii::Table<2, std::unique_ptr<DoFIdentities>>
            vertex_dof_identities(dof_handler.get_fe_collection().size(),
                                  dof_handler.get_fe_collection().size());

          // first identify vertices we want to exclude from working on.
          // specifically, these are the vertices of artificial and ghost
          // cells because at the time when we get here, we do not yet
          // know DoF indices on ghost cells (and we will never know
          // them for artificial cells). this is, at least the case for
          // parallel::distributed::Triangulations.
          //
          // this means that we will not unify DoF indices between locally
          // owned cells and ghost cells, and this is different from what
          // we would do if the triangulation were not split into subdomains.
          // on the other hand, DoF unification is only an optimization: we
          // will still record these identities when we compute hanging
          // node constraints; we just end up with more DoFs than we would
          // if we unified DoF indices also between locally owned and ghost
          // cells, but we end up with a simpler algorithm in return.
          std::vector<bool> include_vertex
            = dof_handler.get_triangulation().get_used_vertices();
          if(dynamic_cast<
               const parallel::distributed::Triangulation<dim, spacedim>*>(
               &dof_handler.get_triangulation())
             != nullptr)
            for(const auto& cell : dof_handler.active_cell_iterators())
              if(!cell->is_locally_owned())
                for(unsigned int v = 0;
                    v < GeometryInfo<dim>::vertices_per_cell;
                    ++v)
                  include_vertex[cell->vertex_index(v)] = false;

          // loop over all vertices and see which one we need to work
          // on
          for(unsigned int vertex_index = 0;
              vertex_index < dof_handler.get_triangulation().n_vertices();
              ++vertex_index)
            if((dof_handler.get_triangulation()
                  .get_used_vertices()[vertex_index]
                == true)
               && (include_vertex[vertex_index] == true))
              {
                const unsigned int n_active_fe_indices
                  = dealii::internal::DoFAccessorImplementation::
                    Implementation::n_active_vertex_fe_indices(dof_handler,
                                                               vertex_index);
                if(n_active_fe_indices > 1)
                  {
                    const unsigned int first_fe_index
                      = dealii::internal::DoFAccessorImplementation::
                        Implementation::nth_active_vertex_fe_index(
                          dof_handler, vertex_index, 0);

                    // loop over all the other FEs with which we want
                    // to identify the DoF indices of the first FE of
                    for(unsigned int f = 1; f < n_active_fe_indices; ++f)
                      {
                        const unsigned int other_fe_index
                          = dealii::internal::DoFAccessorImplementation::
                            Implementation::nth_active_vertex_fe_index(
                              dof_handler, vertex_index, f);

                        // make sure the entry in the equivalence
                        // table exists
                        ensure_existence_of_dof_identities<0>(
                          dof_handler.get_fe(first_fe_index),
                          dof_handler.get_fe(other_fe_index),
                          vertex_dof_identities[first_fe_index]
                                               [other_fe_index]);

                        // then loop through the identities we
                        // have. first get the global numbers of the
                        // dofs we want to identify and make sure they
                        // are not yet constrained to anything else,
                        // except for to each other. use the rule that
                        // we will always constrain the dof with the
                        // higher fe index to the one with the lower,
                        // to avoid circular reasoning.
                        DoFIdentities& identities
                          = *vertex_dof_identities[first_fe_index]
                                                  [other_fe_index];
                        for(unsigned int i = 0; i < identities.size(); ++i)
                          {
                            const types::global_dof_index lower_dof_index
                              = dealii::internal::DoFAccessorImplementation::
                                Implementation::get_vertex_dof_index(
                                  dof_handler,
                                  vertex_index,
                                  first_fe_index,
                                  identities[i].first);
                            const types::global_dof_index higher_dof_index
                              = dealii::internal::DoFAccessorImplementation::
                                Implementation::get_vertex_dof_index(
                                  dof_handler,
                                  vertex_index,
                                  other_fe_index,
                                  identities[i].second);

                            Assert((dof_identities.find(higher_dof_index)
                                    == dof_identities.end())
                                     || (dof_identities[higher_dof_index]
                                         == lower_dof_index),
                                   ExcInternalError());

                            dof_identities[higher_dof_index] = lower_dof_index;
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
          compute_line_dof_identities(hp::DoFHandler<1, spacedim>&)
        {
          return std::map<types::global_dof_index, types::global_dof_index>();
        }

        template <int dim, int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
        compute_line_dof_identities(hp::DoFHandler<dim, spacedim>& dof_handler)
        {
          std::map<types::global_dof_index, types::global_dof_index>
            dof_identities;

          // we will mark lines that we have already treated, so first save and clear
          // the user flags on lines and later restore them
          std::vector<bool> user_flags;
          dof_handler.get_triangulation().save_user_flags_line(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim>&>(
            dof_handler.get_triangulation())
            .clear_user_flags_line();

          // exclude lines that bound cells we don't locally own, because
          // we do not have information about their dofs at this point.
          // this is, at least the case for parallel::distributed::Triangulations.
          //
          // this means that we will not unify DoF indices between locally
          // owned cells and ghost cells, and this is different from what
          // we would do if the triangulation were not split into subdomains.
          // on the other hand, DoF unification is only an optimization: we
          // will still record these identities when we compute hanging
          // node constraints; we just end up with more DoFs than we would
          // if we unified DoF indices also between locally owned and ghost
          // cells, but we end up with a simpler algorithm in return.
          if(dynamic_cast<
               const parallel::distributed::Triangulation<dim, spacedim>*>(
               &dof_handler.get_triangulation())
             != nullptr)
            for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator
                  cell
                = dof_handler.begin_active();
                cell != dof_handler.end();
                ++cell)
              if(cell->is_locally_owned() == false)
                for(unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell;
                    ++l)
                  cell->line(l)->set_user_flag();

          // An implementation of the algorithm described in the hp paper, including
          // the modification mentioned later in the "complications in 3-d" subsections
          //
          // as explained there, we do something only if there are exactly 2 finite
          // elements associated with an object. if there is only one, then there is
          // nothing to do anyway, and if there are 3 or more, then we can get into
          // trouble. note that this only happens for lines in 3d and higher, and for
          // quads only in 4d and higher, so this isn't a particularly frequent case
          //
          // there is one case, however, that we would like to handle (see, for
          // example, the hp/crash_15 testcase): if we have FESystem(FE_Q(2),FE_DGQ(i))
          // elements for a bunch of values 'i', then we should be able to handle this
          // because we can simply unify *all* dofs, not only a some. so what we do
          // is to first treat all pairs of finite elements that have *identical* dofs,
          // and then only deal with those that are not identical of which we can
          // handle at most 2
          dealii::Table<2, std::unique_ptr<DoFIdentities>> line_dof_identities(
            dof_handler.fe_collection.size(), dof_handler.fe_collection.size());

          for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator cell
              = dof_handler.begin_active();
              cell != dof_handler.end();
              ++cell)
            for(unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
              if(cell->line(l)->user_flag_set() == false)
                {
                  const typename hp::DoFHandler<dim, spacedim>::line_iterator
                    line
                    = cell->line(l);
                  line->set_user_flag();

                  unsigned int unique_sets_of_dofs
                    = line->n_active_fe_indices();

                  // do a first loop over all sets of dofs and do identity
                  // uniquification
                  const unsigned int n_active_fe_indices
                    = line->n_active_fe_indices();
                  for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                    for(unsigned int g = f + 1; g < n_active_fe_indices; ++g)
                      {
                        const unsigned int fe_index_1
                          = line->nth_active_fe_index(f),
                          fe_index_2 = line->nth_active_fe_index(g);

                        if((dof_handler.get_fe(fe_index_1).dofs_per_line
                            == dof_handler.get_fe(fe_index_2).dofs_per_line)
                           && (dof_handler.get_fe(fe_index_1).dofs_per_line
                               > 0))
                          {
                            ensure_existence_of_dof_identities<1>(
                              dof_handler.get_fe(fe_index_1),
                              dof_handler.get_fe(fe_index_2),
                              line_dof_identities[fe_index_1][fe_index_2]);
                            // see if these sets of dofs are identical. the first
                            // condition for this is that indeed there are n identities
                            if(line_dof_identities[fe_index_1][fe_index_2]
                                 ->size()
                               == dof_handler.get_fe(fe_index_1).dofs_per_line)
                              {
                                unsigned int i = 0;
                                for(; i < dof_handler.get_fe(fe_index_1)
                                            .dofs_per_line;
                                    ++i)
                                  if(((*(line_dof_identities[fe_index_1]
                                                            [fe_index_2]))[i]
                                        .first
                                      != i)
                                     && ((*(line_dof_identities[fe_index_1]
                                                               [fe_index_2]))[i]
                                           .second
                                         != i))
                                    // not an identity
                                    break;

                                if(i
                                   == dof_handler.get_fe(fe_index_1)
                                        .dofs_per_line)
                                  {
                                    // The line dofs (i.e., the ones interior to a line) of these two finite elements are identical.
                                    // Note that there could be situations when one element still dominates another, e.g.:
                                    // FE_Q(2) x FE_Nothing(dominate) vs
                                    // FE_Q(2) x FE_Q(1)

                                    --unique_sets_of_dofs;

                                    for(unsigned int j = 0;
                                        j < dof_handler.get_fe(fe_index_1)
                                              .dofs_per_line;
                                        ++j)
                                      {
                                        const types::global_dof_index
                                          master_dof_index
                                          = line->dof_index(j, fe_index_1);
                                        const types::global_dof_index
                                          slave_dof_index
                                          = line->dof_index(j, fe_index_2);

                                        // if master dof was already constrained,
                                        // constrain to that one, otherwise constrain
                                        // slave to master
                                        if(dof_identities.find(master_dof_index)
                                           != dof_identities.end())
                                          {
                                            Assert(dof_identities.find(
                                                     dof_identities
                                                       [master_dof_index])
                                                     == dof_identities.end(),
                                                   ExcInternalError());

                                            dof_identities[slave_dof_index]
                                              = dof_identities
                                                [master_dof_index];
                                          }
                                        else
                                          {
                                            Assert((dof_identities.find(
                                                      master_dof_index)
                                                    == dof_identities.end())
                                                     || (dof_identities
                                                           [slave_dof_index]
                                                         == master_dof_index),
                                                   ExcInternalError());

                                            dof_identities[slave_dof_index]
                                              = master_dof_index;
                                          }
                                      }
                                  }
                              }
                          }
                      }

                  // if at this point, there is only one unique set of dofs left, then
                  // we have taken care of everything above. if there are two, then we
                  // need to deal with them here. if there are more, then we punt, as
                  // described in the paper (and mentioned above)
                  //TODO: The check for 'dim==2' was inserted by intuition. It fixes
                  // the previous problems with step-27 in 3D. But an explanation
                  // for this is still required, and what we do here is not what we
                  // describe in the paper!.
                  if((unique_sets_of_dofs == 2) && (dim == 2))
                    {
                      // find out which is the most dominating finite element of the
                      // ones that are used on this line
                      const unsigned int most_dominating_fe_index
                        = get_most_dominating_fe_index<dim, spacedim>(line);

                      // if we found the most dominating element, then use this to eliminate some of
                      // the degrees of freedom by identification. otherwise, the code that computes
                      // hanging node constraints will have to deal with it by computing
                      // appropriate constraints along this face/edge
                      if(most_dominating_fe_index
                         != numbers::invalid_unsigned_int)
                        {
                          const unsigned int n_active_fe_indices
                            = line->n_active_fe_indices();

                          // loop over the indices of all the finite elements that are not
                          // dominating, and identify their dofs to the most dominating
                          // one
                          for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                            if(line->nth_active_fe_index(f)
                               != most_dominating_fe_index)
                              {
                                const unsigned int other_fe_index
                                  = line->nth_active_fe_index(f);

                                ensure_existence_of_dof_identities<1>(
                                  dof_handler.get_fe(most_dominating_fe_index),
                                  dof_handler.get_fe(other_fe_index),
                                  line_dof_identities[most_dominating_fe_index]
                                                     [other_fe_index]);

                                DoFIdentities& identities
                                  = *line_dof_identities
                                      [most_dominating_fe_index]
                                      [other_fe_index];
                                for(unsigned int i = 0; i < identities.size();
                                    ++i)
                                  {
                                    const types::global_dof_index
                                      master_dof_index
                                      = line->dof_index(
                                        identities[i].first,
                                        most_dominating_fe_index);
                                    const types::global_dof_index
                                      slave_dof_index
                                      = line->dof_index(identities[i].second,
                                                        other_fe_index);

                                    Assert(
                                      (dof_identities.find(master_dof_index)
                                       == dof_identities.end())
                                        || (dof_identities[slave_dof_index]
                                            == master_dof_index),
                                      ExcInternalError());

                                    dof_identities[slave_dof_index]
                                      = master_dof_index;
                                  }
                              }
                        }
                    }
                }

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim>&>(
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
        compute_quad_dof_identities(hp::DoFHandler<dim, spacedim>&)
        {
          // this function should only be called for dim<3 where there are
          // no quad dof identies. for dim>=3, the specialization below should
          // take care of it
          Assert(dim < 3, ExcInternalError());

          return std::map<types::global_dof_index, types::global_dof_index>();
        }

        template <int spacedim>
        static std::map<types::global_dof_index, types::global_dof_index>
          compute_quad_dof_identities(hp::DoFHandler<3, spacedim>& dof_handler)
        {
          const int dim = 3;

          std::map<types::global_dof_index, types::global_dof_index>
            dof_identities;

          // we will mark quads that we have already treated, so first
          // save and clear the user flags on quads and later restore
          // them
          std::vector<bool> user_flags;
          dof_handler.get_triangulation().save_user_flags_quad(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim>&>(
            dof_handler.get_triangulation())
            .clear_user_flags_quad();

          // exclude quads that bound cells we don't locally own, because
          // we do not have information about their dofs at this point.
          // this is, at least the case for parallel::distributed::Triangulations.
          //
          // this means that we will not unify DoF indices between locally
          // owned cells and ghost cells, and this is different from what
          // we would do if the triangulation were not split into subdomains.
          // on the other hand, DoF unification is only an optimization: we
          // will still record these identities when we compute hanging
          // node constraints; we just end up with more DoFs than we would
          // if we unified DoF indices also between locally owned and ghost
          // cells, but we end up with a simpler algorithm in return.
          if(dynamic_cast<
               const parallel::distributed::Triangulation<dim, spacedim>*>(
               &dof_handler.get_triangulation())
             != nullptr)
            for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator
                  cell
                = dof_handler.begin_active();
                cell != dof_handler.end();
                ++cell)
              if(cell->is_locally_owned() == false)
                for(unsigned int q = 0; q < GeometryInfo<dim>::quads_per_cell;
                    ++q)
                  cell->quad(q)->set_user_flag();

          // An implementation of the algorithm described in the hp
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
          dealii::Table<2, std::unique_ptr<DoFIdentities>> quad_dof_identities(
            dof_handler.fe_collection.size(), dof_handler.fe_collection.size());

          for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator cell
              = dof_handler.begin_active();
              cell != dof_handler.end();
              ++cell)
            for(unsigned int q = 0; q < GeometryInfo<dim>::quads_per_cell; ++q)
              if((cell->quad(q)->user_flag_set() == false)
                 && (cell->quad(q)->n_active_fe_indices() == 2))
                {
                  const typename hp::DoFHandler<dim, spacedim>::quad_iterator
                    quad
                    = cell->quad(q);
                  quad->set_user_flag();

                  // find out which is the most dominating finite
                  // element of the ones that are used on this quad
                  const unsigned int most_dominating_fe_index
                    = get_most_dominating_fe_index<dim, spacedim>(quad);

                  // if we found the most dominating element, then use
                  // this to eliminate some of the degrees of freedom
                  // by identification. otherwise, the code that
                  // computes hanging node constraints will have to
                  // deal with it by computing appropriate constraints
                  // along this face/edge
                  if(most_dominating_fe_index != numbers::invalid_unsigned_int)
                    {
                      const unsigned int n_active_fe_indices
                        = quad->n_active_fe_indices();

                      // loop over the indices of all the finite
                      // elements that are not dominating, and
                      // identify their dofs to the most dominating
                      // one
                      for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                        if(quad->nth_active_fe_index(f)
                           != most_dominating_fe_index)
                          {
                            const unsigned int other_fe_index
                              = quad->nth_active_fe_index(f);

                            ensure_existence_of_dof_identities<2>(
                              dof_handler.get_fe(most_dominating_fe_index),
                              dof_handler.get_fe(other_fe_index),
                              quad_dof_identities[most_dominating_fe_index]
                                                 [other_fe_index]);

                            DoFIdentities& identities
                              = *quad_dof_identities[most_dominating_fe_index]
                                                    [other_fe_index];
                            for(unsigned int i = 0; i < identities.size(); ++i)
                              {
                                const types::global_dof_index master_dof_index
                                  = quad->dof_index(identities[i].first,
                                                    most_dominating_fe_index);
                                const types::global_dof_index slave_dof_index
                                  = quad->dof_index(identities[i].second,
                                                    other_fe_index);

                                Assert((dof_identities.find(master_dof_index)
                                        == dof_identities.end())
                                         || (dof_identities[slave_dof_index]
                                             == master_dof_index),
                                       ExcInternalError());

                                dof_identities[slave_dof_index]
                                  = master_dof_index;
                              }
                          }
                    }
                }

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim>&>(
            dof_handler.get_triangulation())
            .load_user_flags_quad(user_flags);

          return dof_identities;
        }

        /**
         * Once degrees of freedom have been distributed on all cells, see if
         * we can identify DoFs on neighboring cells. This function does nothing
         * on regular DoFHandlers, but goes through vertices, lines, and quads
         * for hp::DoFHandler objects.
         *
         * Return the final number of degrees of freedom, which is the old one
         * minus however many were identified
         */
        template <int dim, int spacedim>
        static unsigned int
        unify_dof_indices(const DoFHandler<dim, spacedim>&,
                          const unsigned int n_dofs_before_identification,
                          const bool)
        {
          return n_dofs_before_identification;
        }

        template <int dim, int spacedim>
        static unsigned int
        unify_dof_indices(hp::DoFHandler<dim, spacedim>& dof_handler,
                          const unsigned int n_dofs_before_identification,
                          const bool         check_validity)
        {
          // compute the constraints that correspond to unifying
          // dof indices on vertices, lines, and quads. do so
          // in parallel
          std::map<types::global_dof_index, types::global_dof_index>
            all_constrained_indices[dim];

          {
            Threads::TaskGroup<> tasks;

            unsigned int i = 0;
            tasks += Threads::new_task([&, i]() {
              all_constrained_indices[i]
                = compute_vertex_dof_identities(dof_handler);
            });

            if(dim > 1)
              {
                ++i;
                tasks += Threads::new_task([&, i]() {
                  all_constrained_indices[i]
                    = compute_line_dof_identities(dof_handler);
                });
              }

            if(dim > 2)
              {
                ++i;
                tasks += Threads::new_task([&, i]() {
                  all_constrained_indices[i]
                    = compute_quad_dof_identities(dof_handler);
                });
              }

            tasks.join_all();
          }

          // create a vector that contains the new DoF indices; first preset the
          // ones that are identities as determined above, then enumerate the rest
          std::vector<types::global_dof_index> new_dof_indices(
            n_dofs_before_identification, numbers::invalid_dof_index);

          for(const auto& constrained_dof_indices : all_constrained_indices)
            for(const auto& p : constrained_dof_indices)
              {
                Assert(new_dof_indices[p.first] == numbers::invalid_dof_index,
                       ExcInternalError());
                new_dof_indices[p.first] = p.second;
              }

          types::global_dof_index next_free_dof = 0;
          for(types::global_dof_index i = 0; i < n_dofs_before_identification;
              ++i)
            if(new_dof_indices[i] == numbers::invalid_dof_index)
              {
                new_dof_indices[i] = next_free_dof;
                ++next_free_dof;
              }

          // then loop over all those that are constrained and record the
          // new dof number for those:
          for(const auto& constrained_dof_indices : all_constrained_indices)
            for(const auto& p : constrained_dof_indices)
              {
                Assert(new_dof_indices[p.first] != numbers::invalid_dof_index,
                       ExcInternalError());

                new_dof_indices[p.first] = new_dof_indices[p.second];
              }

          for(types::global_dof_index i = 0; i < n_dofs_before_identification;
              ++i)
            {
              Assert(new_dof_indices[i] != numbers::invalid_dof_index,
                     ExcInternalError());
              Assert(new_dof_indices[i] < next_free_dof, ExcInternalError());
            }

          // finally, do the renumbering. verify that previous dof indices
          // were indeed all valid on all cells that we touch if we were
          // told to do so
          renumber_dofs(
            new_dof_indices, IndexSet(0), dof_handler, check_validity);

          return next_free_dof;
        }

        /**
         * Distribute degrees of freedom on all cells, or on cells with the
         * correct subdomain_id if the corresponding argument is not equal to
         * numbers::invalid_subdomain_id. Return the total number of dofs
         * distributed.
         */
        template <class DoFHandlerType>
        static types::global_dof_index
        distribute_dofs(const types::subdomain_id subdomain_id,
                        DoFHandlerType&           dof_handler)
        {
          Assert(dof_handler.get_triangulation().n_levels() > 0,
                 ExcMessage("Empty triangulation"));

          // Step 1: distribute dofs on all cells, but definitely
          // exclude artificial cells
          types::global_dof_index                       next_free_dof = 0;
          typename DoFHandlerType::active_cell_iterator cell
            = dof_handler.begin_active(),
            endc = dof_handler.end();

          for(; cell != endc; ++cell)
            if(!cell->is_artificial())
              if((subdomain_id == numbers::invalid_subdomain_id)
                 || (cell->subdomain_id() == subdomain_id))
                next_free_dof = Implementation::distribute_dofs_on_cell(
                  dof_handler, cell, next_free_dof);

          // Step 2: unify dof indices in case this is an hp DoFHandler
          //
          // during unification, we need to renumber DoF indices. there,
          // we can check that all previous DoF indices were valid, but
          // this only makes sense if we really distributed DoFs on
          // all (non-artificial) cells above
          next_free_dof = unify_dof_indices(
            dof_handler,
            next_free_dof,
            /* check_validity = */
            (subdomain_id == numbers::invalid_subdomain_id));

          update_all_active_cell_dof_indices_caches(dof_handler);

          return next_free_dof;
        }

        /* -------------- distribute_mg_dofs functionality ------------- */

        /**
         * Distribute multilevel dofs on the given cell, with new dofs starting
         * with index @p next_free_dof. Return the next unused index number.
         *
         * Note that unlike for the usual dofs, here all cells and not
         * only active ones are allowed.
         */
        template <int dim, int spacedim>
        static types::global_dof_index
        distribute_mg_dofs_on_cell(
          const typename DoFHandler<dim, spacedim>::level_cell_iterator& cell,
          types::global_dof_index next_free_dof,
          const std::integral_constant<int, 1>&)
        {
          // distribute dofs of vertices
          if(cell->get_fe().dofs_per_vertex > 0)
            for(unsigned int v = 0; v < GeometryInfo<1>::vertices_per_cell; ++v)
              {
                typename DoFHandler<dim, spacedim>::level_cell_iterator neighbor
                  = cell->neighbor(v);

                if(neighbor.state() == IteratorState::valid)
                  {
                    // has neighbor already been processed?
                    if(neighbor->user_flag_set()
                       && (neighbor->level() == cell->level()))
                      // copy dofs if the neighbor is on the same level (only then are
                      // mg dofs the same)
                      {
                        if(v == 0)
                          for(unsigned int d = 0;
                              d < cell->get_fe().dofs_per_vertex;
                              ++d)
                            cell->set_mg_vertex_dof_index(
                              cell->level(),
                              0,
                              d,
                              neighbor->mg_vertex_dof_index(
                                cell->level(), 1, d));
                        else
                          for(unsigned int d = 0;
                              d < cell->get_fe().dofs_per_vertex;
                              ++d)
                            cell->set_mg_vertex_dof_index(
                              cell->level(),
                              1,
                              d,
                              neighbor->mg_vertex_dof_index(
                                cell->level(), 0, d));

                        // next neighbor
                        continue;
                      }
                  }

                // otherwise: create dofs newly
                for(unsigned int d = 0; d < cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index(
                    cell->level(), v, d, next_free_dof++);
              }

          // dofs of line
          if(cell->get_fe().dofs_per_line > 0)
            for(unsigned int d = 0; d < cell->get_fe().dofs_per_line; ++d)
              cell->set_mg_dof_index(cell->level(), d, next_free_dof++);

          // note that this cell has been processed
          cell->set_user_flag();

          return next_free_dof;
        }

        template <int dim, int spacedim>
        static types::global_dof_index
        distribute_mg_dofs_on_cell(
          const typename DoFHandler<dim, spacedim>::level_cell_iterator& cell,
          types::global_dof_index next_free_dof,
          const std::integral_constant<int, 2>&)
        {
          if(cell->get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<2>::vertices_per_cell;
                ++vertex)
              // check whether dofs for this
              // vertex have been distributed
              // (only check the first dof)
              if(cell->mg_vertex_dof_index(cell->level(), vertex, 0)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0; d < cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index(
                    cell->level(), vertex, d, next_free_dof++);

          // for the four sides
          if(cell->get_fe().dofs_per_line > 0)
            for(unsigned int side = 0; side < GeometryInfo<2>::faces_per_cell;
                ++side)
              {
                typename DoFHandler<dim, spacedim>::line_iterator line
                  = cell->line(side);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (check only first dof)
                if(line->mg_dof_index(cell->level(), 0)
                   == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for(unsigned int d = 0; d < cell->get_fe().dofs_per_line; ++d)
                    line->set_mg_dof_index(cell->level(), d, next_free_dof++);
              }

          // dofs of quad
          if(cell->get_fe().dofs_per_quad > 0)
            for(unsigned int d = 0; d < cell->get_fe().dofs_per_quad; ++d)
              cell->set_mg_dof_index(cell->level(), d, next_free_dof++);

          // note that this cell has been processed
          cell->set_user_flag();

          return next_free_dof;
        }

        template <int dim, int spacedim>
        static types::global_dof_index
        distribute_mg_dofs_on_cell(
          const typename DoFHandler<dim, spacedim>::level_cell_iterator& cell,
          types::global_dof_index next_free_dof,
          const std::integral_constant<int, 3>&)
        {
          if(cell->get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for(unsigned int vertex = 0;
                vertex < GeometryInfo<3>::vertices_per_cell;
                ++vertex)
              // check whether dofs for this vertex have been distributed
              // (only check the first dof)
              if(cell->mg_vertex_dof_index(cell->level(), vertex, 0)
                 == numbers::invalid_dof_index)
                for(unsigned int d = 0; d < cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index(
                    cell->level(), vertex, d, next_free_dof++);

          // for the lines
          if(cell->get_fe().dofs_per_line > 0)
            for(unsigned int l = 0; l < GeometryInfo<3>::lines_per_cell; ++l)
              {
                typename DoFHandler<dim, spacedim>::line_iterator line
                  = cell->line(l);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if(line->mg_dof_index(cell->level(), 0)
                   == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for(unsigned int d = 0; d < cell->get_fe().dofs_per_line; ++d)
                    line->set_mg_dof_index(cell->level(), d, next_free_dof++);
              }

          // for the quads
          if(cell->get_fe().dofs_per_quad > 0)
            for(unsigned int q = 0; q < GeometryInfo<3>::quads_per_cell; ++q)
              {
                typename DoFHandler<dim, spacedim>::quad_iterator quad
                  = cell->quad(q);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if(quad->mg_dof_index(cell->level(), 0)
                   == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for(unsigned int d = 0; d < cell->get_fe().dofs_per_quad; ++d)
                    quad->set_mg_dof_index(cell->level(), d, next_free_dof++);
              }

          // dofs of cell
          if(cell->get_fe().dofs_per_hex > 0)
            for(unsigned int d = 0; d < cell->get_fe().dofs_per_hex; ++d)
              cell->set_mg_dof_index(cell->level(), d, next_free_dof++);

          // note that this cell has been processed
          cell->set_user_flag();

          return next_free_dof;
        }

        // same for the hp::DoFHandler
        template <int spacedim>
        static types::global_dof_index
        distribute_mg_dofs_on_cell(
          const hp::DoFHandler<1, spacedim>& dof_handler,
          const typename hp::DoFHandler<1, spacedim>::active_cell_iterator&
                                  cell,
          types::global_dof_index next_free_dof)
        {
          (void) dof_handler;
          (void) cell;
          (void) next_free_dof;
          return 0;
        }

        template <int spacedim>
        static types::global_dof_index
        distribute_mg_dofs_on_cell(
          const hp::DoFHandler<2, spacedim>& dof_handler,
          const typename hp::DoFHandler<2, spacedim>::active_cell_iterator&
                                  cell,
          types::global_dof_index next_free_dof)
        {
          (void) dof_handler;
          (void) cell;
          (void) next_free_dof;
          return 0;
        }

        template <int spacedim>
        static types::global_dof_index
        distribute_mg_dofs_on_cell(
          const hp::DoFHandler<3, spacedim>& dof_handler,
          const typename hp::DoFHandler<3, spacedim>::active_cell_iterator&
                                  cell,
          types::global_dof_index next_free_dof)
        {
          (void) dof_handler;
          (void) cell;
          (void) next_free_dof;
          return 0;
        }

        template <class DoFHandlerType>
        static types::global_dof_index
        distribute_dofs_on_level(const types::subdomain_id level_subdomain_id,
                                 DoFHandlerType&           dof_handler,
                                 const unsigned int        level)
        {
          const unsigned int dim      = DoFHandlerType::dimension;
          const unsigned int spacedim = DoFHandlerType::space_dimension;

          const dealii::Triangulation<dim, spacedim>& tria
            = dof_handler.get_triangulation();
          Assert(tria.n_levels() > 0, ExcMessage("Empty triangulation"));
          if(level >= tria.n_levels())
            return 0; //this is allowed for multigrid

          // Clear user flags because we will need them. But first we save
          // them and make sure that we restore them later such that at
          // the end of this function the Triangulation will be in the
          // same state as it was at the beginning of this function.
          std::vector<bool> user_flags;
          tria.save_user_flags(user_flags);
          const_cast<dealii::Triangulation<dim, spacedim>&>(tria)
            .clear_user_flags();

          types::global_dof_index next_free_dof = 0;
          typename DoFHandler<dim, spacedim>::level_cell_iterator cell
            = dof_handler.begin(level),
            endc = dof_handler.end(level);

          for(; cell != endc; ++cell)
            if((level_subdomain_id == numbers::invalid_subdomain_id)
               || (cell->level_subdomain_id() == level_subdomain_id))
              next_free_dof
                = Implementation::distribute_mg_dofs_on_cell<dim, spacedim>(
                  cell, next_free_dof, std::integral_constant<int, dim>());

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim, spacedim>&>(tria)
            .load_user_flags(user_flags);

          return next_free_dof;
        }

        /* --------------------- renumber_dofs functionality ---------------- */

        /**
         * The part of the renumber_dofs() functionality that is dimension
         * independent because it renumbers the DoF indices on vertices
         * (which exist for all dimensions).
         *
         * See renumber_dofs() for the meaning of the arguments.
         */
        template <int dim, int spacedim>
        static void
        renumber_vertex_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          DoFHandler<dim, spacedim>&                  dof_handler,
          const bool                                  check_validity)
        {
          // we can not use cell iterators in this function since then
          // we would renumber the dofs on the interface of two cells
          // more than once. Anyway, this way it's not only more
          // correct but also faster; note, however, that dof numbers
          // may be invalid_dof_index, namely when the appropriate
          // vertex/line/etc is unused
          for(std::vector<types::global_dof_index>::iterator i
              = dof_handler.vertex_dofs.begin();
              i != dof_handler.vertex_dofs.end();
              ++i)
            if(*i != numbers::invalid_dof_index)
              *i = (indices.size() == 0) ?
                     (new_numbers[*i]) :
                     (new_numbers[indices.index_within_set(*i)]);
            else if(check_validity)
              // if index is invalid_dof_index: check if this one
              // really is unused
              Assert(dof_handler.get_triangulation().vertex_used(
                       (i - dof_handler.vertex_dofs.begin())
                       / dof_handler.get_fe().dofs_per_vertex)
                       == false,
                     ExcInternalError());
        }

        /**
         * The part of the renumber_dofs() functionality that is dimension
         * independent because it renumbers the DoF indices on cell interiors
         * (which exist for all dimensions).
         *
         * See renumber_dofs() for the meaning of the arguments.
         */
        template <int dim, int spacedim>
        static void
        renumber_cell_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          DoFHandler<dim, spacedim>&                  dof_handler)
        {
          for(unsigned int level = 0; level < dof_handler.levels.size();
              ++level)
            for(std::vector<types::global_dof_index>::iterator i
                = dof_handler.levels[level]->dof_object.dofs.begin();
                i != dof_handler.levels[level]->dof_object.dofs.end();
                ++i)
              if(*i != numbers::invalid_dof_index)
                *i = ((indices.size() == 0) ?
                        new_numbers[*i] :
                        new_numbers[indices.index_within_set(*i)]);
        }

        /**
         * The part of the renumber_dofs() functionality that operates on faces.
         * This part is dimension dependent and so needs to be implemented in
         * three separate specializations of the function.
         *
         * See renumber_dofs() for the meaning of the arguments.
         */
        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index>& /*new_numbers*/,
          const IndexSet& /*indices*/,
          DoFHandler<1, spacedim>& /*dof_handler*/)
        {
          // nothing to do in 1d since there are no separate faces
        }

        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          DoFHandler<2, spacedim>&                    dof_handler)
        {
          // treat dofs on lines
          for(std::vector<types::global_dof_index>::iterator i
              = dof_handler.faces->lines.dofs.begin();
              i != dof_handler.faces->lines.dofs.end();
              ++i)
            if(*i != numbers::invalid_dof_index)
              *i = ((indices.size() == 0) ?
                      new_numbers[*i] :
                      new_numbers[indices.index_within_set(*i)]);
        }

        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          DoFHandler<3, spacedim>&                    dof_handler)
        {
          // treat dofs on lines
          for(std::vector<types::global_dof_index>::iterator i
              = dof_handler.faces->lines.dofs.begin();
              i != dof_handler.faces->lines.dofs.end();
              ++i)
            if(*i != numbers::invalid_dof_index)
              *i = ((indices.size() == 0) ?
                      new_numbers[*i] :
                      new_numbers[indices.index_within_set(*i)]);

          // treat dofs on quads
          for(std::vector<types::global_dof_index>::iterator i
              = dof_handler.faces->quads.dofs.begin();
              i != dof_handler.faces->quads.dofs.end();
              ++i)
            if(*i != numbers::invalid_dof_index)
              *i = ((indices.size() == 0) ?
                      new_numbers[*i] :
                      new_numbers[indices.index_within_set(*i)]);
        }

        template <int dim, int spacedim>
        static void
        renumber_vertex_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          hp::DoFHandler<dim, spacedim>&              dof_handler,
          const bool                                  check_validity)
        {
          for(unsigned int vertex_index = 0;
              vertex_index < dof_handler.get_triangulation().n_vertices();
              ++vertex_index)
            {
              const unsigned int n_active_fe_indices
                = dealii::internal::DoFAccessorImplementation::Implementation::
                  n_active_vertex_fe_indices(dof_handler, vertex_index);

              // if this vertex is unused, then we really ought not to have allocated
              // any space for it, i.e., n_active_fe_indices should be zero, and
              // there is no space to actually store dof indices for this vertex
              if(dof_handler.get_triangulation().vertex_used(vertex_index)
                 == false)
                Assert(n_active_fe_indices == 0, ExcInternalError());

              // otherwise the vertex is used; it may still not hold any dof indices
              // if it is located on an artificial cell and not adjacent to a ghost
              // cell, but in that case there is simply nothing for us to do
              for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                {
                  const unsigned int fe_index
                    = dealii::internal::DoFAccessorImplementation::
                      Implementation::nth_active_vertex_fe_index(
                        dof_handler, vertex_index, f);

                  for(unsigned int d = 0;
                      d < dof_handler.get_fe(fe_index).dofs_per_vertex;
                      ++d)
                    {
                      const types::global_dof_index old_dof_index
                        = dealii::internal::DoFAccessorImplementation::
                          Implementation::get_vertex_dof_index(
                            dof_handler, vertex_index, fe_index, d);

                      // if check_validity was set, then we are to verify that the
                      // previous indices were all valid. this really should be
                      // the case: we allocated space for these vertex dofs,
                      // i.e., at least one adjacent cell has a valid
                      // active_fe_index, so there are DoFs that really live
                      // on this vertex. if check_validity is set, then we
                      // must make sure that they have been set to something
                      // useful
                      if(check_validity)
                        Assert(old_dof_index != numbers::invalid_dof_index,
                               ExcInternalError());

                      if(old_dof_index != numbers::invalid_dof_index)
                        dealii::internal::DoFAccessorImplementation::
                          Implementation::set_vertex_dof_index(
                            dof_handler,
                            vertex_index,
                            fe_index,
                            d,
                            (indices.size() == 0) ?
                              (new_numbers[old_dof_index]) :
                              (new_numbers[indices.index_within_set(
                                old_dof_index)]));
                    }
                }
            }
        }

        template <int dim, int spacedim>
        static void
        renumber_cell_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          hp::DoFHandler<dim, spacedim>&              dof_handler)
        {
          for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator cell
              = dof_handler.begin_active();
              cell != dof_handler.end();
              ++cell)
            if(!cell->is_artificial())
              {
                const unsigned int fe_index = cell->active_fe_index();

                for(unsigned int d = 0;
                    d < dof_handler.get_fe(fe_index)
                          .template n_dofs_per_object<dim>();
                    ++d)
                  {
                    const types::global_dof_index old_dof_index
                      = cell->dof_index(d, fe_index);
                    if(old_dof_index != numbers::invalid_dof_index)
                      cell->set_dof_index(
                        d,
                        (indices.size() == 0) ?
                          (new_numbers[old_dof_index]) :
                          (new_numbers[indices.index_within_set(
                            old_dof_index)]),
                        fe_index);
                  }
              }
        }

        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index>& /*new_numbers*/,
          const IndexSet& /*indices*/,
          hp::DoFHandler<1, spacedim>& /*dof_handler*/)
        {
          // nothing to do in 1d since there are no separate faces -- we've already
          // taken care of this when dealing with the vertices
        }

        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          hp::DoFHandler<2, spacedim>&                dof_handler)
        {
          const unsigned int dim = 2;

          // deal with DoFs on lines
          {
            // save user flags on lines so we can use them to mark lines
            // we've already treated
            std::vector<bool> saved_line_user_flags;
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .save_user_flags_line(saved_line_user_flags);
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .clear_user_flags_line();

            for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator
                  cell
                = dof_handler.begin_active();
                cell != dof_handler.end();
                ++cell)
              if(!cell->is_artificial())
                for(unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell;
                    ++l)
                  if(cell->line(l)->user_flag_set() == false)
                    {
                      const typename hp::DoFHandler<dim,
                                                    spacedim>::line_iterator
                        line
                        = cell->line(l);
                      line->set_user_flag();

                      const unsigned int n_active_fe_indices
                        = line->n_active_fe_indices();

                      for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                        {
                          const unsigned int fe_index
                            = line->nth_active_fe_index(f);

                          for(unsigned int d = 0;
                              d < dof_handler.get_fe(fe_index).dofs_per_line;
                              ++d)
                            {
                              const types::global_dof_index old_dof_index
                                = line->dof_index(d, fe_index);
                              if(old_dof_index != numbers::invalid_dof_index)
                                line->set_dof_index(
                                  d,
                                  (indices.size() == 0) ?
                                    (new_numbers[old_dof_index]) :
                                    (new_numbers[indices.index_within_set(
                                      old_dof_index)]),
                                  fe_index);
                            }
                        }
                    }

            // at the end, restore the user
            // flags for the lines
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .load_user_flags_line(saved_line_user_flags);
          }
        }

        template <int spacedim>
        static void
        renumber_face_dofs(
          const std::vector<types::global_dof_index>& new_numbers,
          const IndexSet&                             indices,
          hp::DoFHandler<3, spacedim>&                dof_handler)
        {
          const unsigned int dim = 3;

          // deal with DoFs on lines
          {
            // save user flags on lines so we can use them to mark lines
            // we've already treated
            std::vector<bool> saved_line_user_flags;
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .save_user_flags_line(saved_line_user_flags);
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .clear_user_flags_line();

            for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator
                  cell
                = dof_handler.begin_active();
                cell != dof_handler.end();
                ++cell)
              if(!cell->is_artificial())
                for(unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell;
                    ++l)
                  if(cell->line(l)->user_flag_set() == false)
                    {
                      const typename hp::DoFHandler<dim,
                                                    spacedim>::line_iterator
                        line
                        = cell->line(l);
                      line->set_user_flag();

                      const unsigned int n_active_fe_indices
                        = line->n_active_fe_indices();

                      for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                        {
                          const unsigned int fe_index
                            = line->nth_active_fe_index(f);

                          for(unsigned int d = 0;
                              d < dof_handler.get_fe(fe_index).dofs_per_line;
                              ++d)
                            {
                              const types::global_dof_index old_dof_index
                                = line->dof_index(d, fe_index);
                              if(old_dof_index != numbers::invalid_dof_index)
                                line->set_dof_index(
                                  d,
                                  (indices.size() == 0) ?
                                    (new_numbers[old_dof_index]) :
                                    (new_numbers[indices.index_within_set(
                                      old_dof_index)]),
                                  fe_index);
                            }
                        }
                    }

            // at the end, restore the user
            // flags for the lines
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .load_user_flags_line(saved_line_user_flags);
          }

          // then deal with dofs on quads
          {
            std::vector<bool> saved_quad_user_flags;
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .save_user_flags_quad(saved_quad_user_flags);
            const_cast<dealii::Triangulation<dim, spacedim>&>(
              dof_handler.get_triangulation())
              .clear_user_flags_quad();

            for(typename hp::DoFHandler<dim, spacedim>::active_cell_iterator
                  cell
                = dof_handler.begin_active();
                cell != dof_handler.end();
                ++cell)
              if(!cell->is_artificial())
                for(unsigned int q = 0; q < GeometryInfo<dim>::quads_per_cell;
                    ++q)
                  if(cell->quad(q)->user_flag_set() == false)
                    {
                      const typename hp::DoFHandler<dim,
                                                    spacedim>::quad_iterator
                        quad
                        = cell->quad(q);
                      quad->set_user_flag();

                      const unsigned int n_active_fe_indices
                        = quad->n_active_fe_indices();

                      for(unsigned int f = 0; f < n_active_fe_indices; ++f)
                        {
                          const unsigned int fe_index
                            = quad->nth_active_fe_index(f);

                          for(unsigned int d = 0;
                              d < dof_handler.get_fe(fe_index).dofs_per_quad;
                              ++d)
                            {
                              const types::global_dof_index old_dof_index
                                = quad->dof_index(d, fe_index);
                              if(old_dof_index != numbers::invalid_dof_index)
                                quad->set_dof_index(
                                  d,
                                  (indices.size() == 0) ?
                                    (new_numbers[old_dof_index]) :
                                    (new_numbers[indices.index_within_set(
                                      old_dof_index)]),
                                  fe_index);
                            }
                        }
                    }

            // at the end, restore the user flags for the quads
            const_cast<dealii::Triangulation<dim, spacedim>&>(
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
        template <class DoFHandlerType>
        static void
        renumber_dofs(const std::vector<types::global_dof_index>& new_numbers,
                      const IndexSet&                             indices,
                      DoFHandlerType&                             dof_handler,
                      const bool check_validity)
        {
          if(DoFHandlerType::dimension == 1)
            Assert(indices == IndexSet(0), ExcNotImplemented());

          // renumber DoF indices on vertices, cells, and faces. this
          // can be done in parallel because the respective functions
          // work on separate data structures
          Threads::TaskGroup<> tasks;
          tasks += Threads::new_task([&]() {
            renumber_vertex_dofs(
              new_numbers, indices, dof_handler, check_validity);
          });
          tasks += Threads::new_task(
            [&]() { renumber_face_dofs(new_numbers, indices, dof_handler); });
          tasks += Threads::new_task(
            [&]() { renumber_cell_dofs(new_numbers, indices, dof_handler); });
          tasks.join_all();

          // update the cache used for cell dof indices
          update_all_active_cell_dof_indices_caches(dof_handler);
        }

        /* --------------------- renumber_mg_dofs functionality ---------------- */

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
          const std::vector<dealii::types::global_dof_index>& new_numbers,
          const IndexSet&                                     indices,
          DoFHandler<dim, spacedim>&                          dof_handler,
          const unsigned int                                  level,
          const bool                                          check_validity)
        {
          Assert(level < dof_handler.get_triangulation().n_levels(),
                 ExcInternalError());

          for(typename std::vector<
                typename DoFHandler<dim, spacedim>::MGVertexDoFs>::iterator i
              = dof_handler.mg_vertex_dofs.begin();
              i != dof_handler.mg_vertex_dofs.end();
              ++i)
            // if the present vertex lives on the current level
            if((i->get_coarsest_level() <= level)
               && (i->get_finest_level() >= level))
              for(unsigned int d = 0; d < dof_handler.get_fe().dofs_per_vertex;
                  ++d)
                {
                  const dealii::types::global_dof_index idx = i->get_index(
                    level, d, dof_handler.get_fe().dofs_per_vertex);

                  if(check_validity)
                    Assert(idx != numbers::invalid_dof_index,
                           ExcInternalError());

                  if(idx != numbers::invalid_dof_index)
                    i->set_index(
                      level,
                      d,
                      dof_handler.get_fe().dofs_per_vertex,
                      (indices.size() == 0) ?
                        (new_numbers[idx]) :
                        (new_numbers[indices.index_within_set(idx)]));
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
          const std::vector<dealii::types::global_dof_index>& new_numbers,
          const IndexSet&                                     indices,
          DoFHandler<dim, spacedim>&                          dof_handler,
          const unsigned int                                  level)
        {
          for(std::vector<types::global_dof_index>::iterator i
              = dof_handler.mg_levels[level]->dof_object.dofs.begin();
              i != dof_handler.mg_levels[level]->dof_object.dofs.end();
              ++i)
            {
              if(*i != numbers::invalid_dof_index)
                {
                  Assert((indices.size() > 0 ? indices.is_element(*i) :
                                               (*i < new_numbers.size())),
                         ExcInternalError());
                  *i = (indices.size() == 0) ?
                         (new_numbers[*i]) :
                         (new_numbers[indices.index_within_set(*i)]);
                }
            }
        }

        /**
         * The part of the renumber_mg_dofs() functionality that operates on faces.
         * This part is dimension dependent and so needs to be implemented in
         * three separate specializations of the function.
         *
         * See renumber_mg_dofs() for the meaning of the arguments.
         */
        template <int spacedim>
        static void
        renumber_face_mg_dofs(
          const std::vector<types::global_dof_index>& /*new_numbers*/,
          const IndexSet& /*indices*/,
          DoFHandler<1, spacedim>& /*dof_handler*/,
          const unsigned int /*level*/,
          const bool /*check_validity*/)
        {
          // nothing to do in 1d because there are no separate faces
        }

        template <int spacedim>
        static void
        renumber_face_mg_dofs(
          const std::vector<dealii::types::global_dof_index>& new_numbers,
          const IndexSet&                                     indices,
          DoFHandler<2, spacedim>&                            dof_handler,
          const unsigned int                                  level,
          const bool                                          check_validity)
        {
          if(dof_handler.get_fe().dofs_per_line > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<2, spacedim>&>(
                dof_handler.get_triangulation())
                .clear_user_flags();

              // flag all lines adjacent to cells of the current
              // level, as those lines logically belong to the same
              // level as the cell, at least for for isotropic
              // refinement
              typename DoFHandler<2, spacedim>::level_cell_iterator cell,
                endc = dof_handler.end(level);
              for(cell = dof_handler.begin(level); cell != endc; ++cell)
                for(unsigned int line = 0;
                    line < GeometryInfo<2>::faces_per_cell;
                    ++line)
                  cell->face(line)->set_user_flag();

              for(typename DoFHandler<2, spacedim>::cell_iterator cell
                  = dof_handler.begin();
                  cell != dof_handler.end();
                  ++cell)
                for(unsigned int l = 0; l < GeometryInfo<2>::lines_per_cell;
                    ++l)
                  if(cell->line(l)->user_flag_set())
                    {
                      for(unsigned int d = 0;
                          d < dof_handler.get_fe().dofs_per_line;
                          ++d)
                        {
                          const dealii::types::global_dof_index idx
                            = cell->line(l)->mg_dof_index(level, d);
                          if(check_validity)
                            Assert(idx != numbers::invalid_dof_index,
                                   ExcInternalError());

                          if(idx != numbers::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index(
                              level,
                              d,
                              ((indices.size() == 0) ?
                                 new_numbers[idx] :
                                 new_numbers[indices.index_within_set(idx)]));
                        }
                      cell->line(l)->clear_user_flag();
                    }
              // finally, restore user flags
              const_cast<dealii::Triangulation<2, spacedim>&>(
                dof_handler.get_triangulation())
                .load_user_flags(user_flags);
            }
        }

        template <int spacedim>
        static void
        renumber_face_mg_dofs(
          const std::vector<dealii::types::global_dof_index>& new_numbers,
          const IndexSet&                                     indices,
          DoFHandler<3, spacedim>&                            dof_handler,
          const unsigned int                                  level,
          const bool                                          check_validity)
        {
          if(dof_handler.get_fe().dofs_per_line > 0
             || dof_handler.get_fe().dofs_per_quad > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<3, spacedim>&>(
                dof_handler.get_triangulation())
                .clear_user_flags();

              // flag all lines adjacent to cells of the current
              // level, as those lines logically belong to the same
              // level as the cell, at least for isotropic refinement
              typename DoFHandler<3, spacedim>::level_cell_iterator cell,
                endc = dof_handler.end(level);
              for(cell = dof_handler.begin(level); cell != endc; ++cell)
                for(unsigned int line = 0;
                    line < GeometryInfo<3>::lines_per_cell;
                    ++line)
                  cell->line(line)->set_user_flag();

              for(typename DoFHandler<3, spacedim>::cell_iterator cell
                  = dof_handler.begin();
                  cell != dof_handler.end();
                  ++cell)
                for(unsigned int l = 0; l < GeometryInfo<3>::lines_per_cell;
                    ++l)
                  if(cell->line(l)->user_flag_set())
                    {
                      for(unsigned int d = 0;
                          d < dof_handler.get_fe().dofs_per_line;
                          ++d)
                        {
                          const dealii::types::global_dof_index idx
                            = cell->line(l)->mg_dof_index(level, d);
                          if(check_validity)
                            Assert(idx != numbers::invalid_dof_index,
                                   ExcInternalError());

                          if(idx != numbers::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index(
                              level,
                              d,
                              ((indices.size() == 0) ?
                                 new_numbers[idx] :
                                 new_numbers[indices.index_within_set(idx)]));
                        }
                      cell->line(l)->clear_user_flag();
                    }

              // flag all quads adjacent to cells of the current level, as
              // those quads logically belong to the same level as the cell,
              // at least for isotropic refinement
              for(cell = dof_handler.begin(level); cell != endc; ++cell)
                for(unsigned int quad = 0;
                    quad < GeometryInfo<3>::quads_per_cell;
                    ++quad)
                  cell->quad(quad)->set_user_flag();

              for(typename DoFHandler<3, spacedim>::cell_iterator cell
                  = dof_handler.begin();
                  cell != dof_handler.end();
                  ++cell)
                for(unsigned int l = 0; l < GeometryInfo<3>::quads_per_cell;
                    ++l)
                  if(cell->quad(l)->user_flag_set())
                    {
                      for(unsigned int d = 0;
                          d < dof_handler.get_fe().dofs_per_quad;
                          ++d)
                        {
                          const dealii::types::global_dof_index idx
                            = cell->quad(l)->mg_dof_index(level, d);
                          if(check_validity)
                            Assert(idx != numbers::invalid_dof_index,
                                   ExcInternalError());

                          if(idx != numbers::invalid_dof_index)
                            cell->quad(l)->set_mg_dof_index(
                              level,
                              d,
                              ((indices.size() == 0) ?
                                 new_numbers[idx] :
                                 new_numbers[indices.index_within_set(idx)]));
                        }
                      cell->quad(l)->clear_user_flag();
                    }

              // finally, restore user flags
              const_cast<dealii::Triangulation<3, spacedim>&>(
                dof_handler.get_triangulation())
                .load_user_flags(user_flags);
            }
        }

        template <int dim, int spacedim>
        static void
        renumber_mg_dofs(
          const std::vector<dealii::types::global_dof_index>& new_numbers,
          const IndexSet&                                     indices,
          DoFHandler<dim, spacedim>&                          dof_handler,
          const unsigned int                                  level,
          const bool                                          check_validity)
        {
          Assert(level < dof_handler.get_triangulation().n_global_levels(),
                 ExcInternalError());

          // renumber DoF indices on vertices, cells, and faces. this
          // can be done in parallel because the respective functions
          // work on separate data structures
          Threads::TaskGroup<> tasks;
          tasks += Threads::new_task([&]() {
            renumber_vertex_mg_dofs(
              new_numbers, indices, dof_handler, level, check_validity);
          });
          tasks += Threads::new_task([&]() {
            renumber_face_mg_dofs(
              new_numbers, indices, dof_handler, level, check_validity);
          });
          tasks += Threads::new_task([&]() {
            renumber_cell_mg_dofs(new_numbers, indices, dof_handler, level);
          });
          tasks.join_all();
        }

        template <int dim, int spacedim>
        static void
        renumber_mg_dofs(
          const std::vector<dealii::types::global_dof_index>& /*new_numbers*/,
          const IndexSet& /*indices*/,
          hp::DoFHandler<dim, spacedim>& /*dof_handler*/,
          const unsigned int /*level*/,
          const bool /*check_validity*/)
        {
          Assert(false, ExcNotImplemented());
        }
      };

      /* --------------------- class Sequential ---------------- */

      template <class DoFHandlerType>
      Sequential<DoFHandlerType>::Sequential(DoFHandlerType& dof_handler)
        : dof_handler(&dof_handler)
      {}

      template <class DoFHandlerType>
      NumberCache
      Sequential<DoFHandlerType>::distribute_dofs() const
      {
        const types::global_dof_index n_dofs = Implementation::distribute_dofs(
          numbers::invalid_subdomain_id, *dof_handler);

        // return a sequential, complete index set
        return NumberCache(n_dofs);
      }

      template <class DoFHandlerType>
      std::vector<NumberCache>
      Sequential<DoFHandlerType>::distribute_mg_dofs() const
      {
        std::vector<bool> user_flags;
        dof_handler->get_triangulation().save_user_flags(user_flags);

        const_cast<dealii::Triangulation<DoFHandlerType::dimension,
                                         DoFHandlerType::space_dimension>&>(
          dof_handler->get_triangulation())
          .clear_user_flags();

        std::vector<NumberCache> number_caches;
        number_caches.reserve(dof_handler->get_triangulation().n_levels());
        for(unsigned int level = 0;
            level < dof_handler->get_triangulation().n_levels();
            ++level)
          {
            // first distribute dofs on this level
            const types::global_dof_index n_level_dofs
              = Implementation::distribute_dofs_on_level(
                numbers::invalid_subdomain_id, *dof_handler, level);

            // then add a complete, sequential index set
            number_caches.emplace_back(NumberCache(n_level_dofs));
          }

        const_cast<dealii::Triangulation<DoFHandlerType::dimension,
                                         DoFHandlerType::space_dimension>&>(
          dof_handler->get_triangulation())
          .load_user_flags(user_flags);

        return number_caches;
      }

      template <class DoFHandlerType>
      NumberCache
      Sequential<DoFHandlerType>::renumber_dofs(
        const std::vector<types::global_dof_index>& new_numbers) const
      {
        Implementation::renumber_dofs(
          new_numbers, IndexSet(0), *dof_handler, true);

        // return a sequential, complete index set. take into account that the
        // number of DoF indices may in fact be smaller than there were before
        // if some previously separately numbered dofs have been identified.
        // this is, for example, what the hp::DoFHandler does: it first
        // enumerates all DoFs on cells independently, and then unifies
        // some located at vertices or faces; this leaves us with fewer
        // DoFs than there were before, so use the largest index as
        // the one to determine the size of the index space
        return NumberCache(
          *std::max_element(new_numbers.begin(), new_numbers.end()) + 1);
      }

      template <class DoFHandlerType>
      NumberCache
      Sequential<DoFHandlerType>::renumber_mg_dofs(
        const unsigned int                          level,
        const std::vector<types::global_dof_index>& new_numbers) const
      {
        Implementation::renumber_mg_dofs(
          new_numbers, IndexSet(0), *dof_handler, level, true);

        // return a sequential, complete index set
        return NumberCache(new_numbers.size());
      }

      /* --------------------- class ParallelShared ---------------- */

      template <class DoFHandlerType>
      ParallelShared<DoFHandlerType>::ParallelShared(
        DoFHandlerType& dof_handler)
        : dof_handler(&dof_handler)
      {}

      namespace
      {
        /**
         * This function is a variation of DoFTools::get_dof_subdomain_association()
         * with the exception that it is (i) specialized for
         * parallel::shared::Triangulation objects, and (ii) does not assume that the
         * internal number cache of the DoFHandler has already been set up. In can,
         * consequently, be called at a point when we are still distributing degrees
         * of freedom.
         */
        template <class DoFHandlerType>
        std::vector<types::subdomain_id>
        get_dof_subdomain_association(const DoFHandlerType&         dof_handler,
                                      const types::global_dof_index n_dofs,
                                      const unsigned int            n_procs)
        {
          (void) n_procs;
          std::vector<types::subdomain_id> subdomain_association(
            n_dofs, numbers::invalid_subdomain_id);
          std::vector<types::global_dof_index> local_dof_indices;
          local_dof_indices.reserve(DoFTools::max_dofs_per_cell(dof_handler));

          // loop over all cells and record which subdomain a DoF belongs to.
          // give to the smaller subdomain_id in case it is on an interface
          typename DoFHandlerType::active_cell_iterator cell
            = dof_handler.begin_active(),
            endc = dof_handler.end();
          for(; cell != endc; ++cell)
            {
              // get the owner of the cell; note that we have made sure above that
              // all cells are either locally owned or ghosts (not artificial), so
              // this call will always yield the true owner
              const types::subdomain_id subdomain_id = cell->subdomain_id();
              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
              local_dof_indices.resize(dofs_per_cell);
              cell->get_dof_indices(local_dof_indices);

              // set subdomain ids. if dofs already have their values set then
              // they must be on partition interfaces. In that case assign them to the
              // processor with the smaller subdomain id.
              for(unsigned int i = 0; i < dofs_per_cell; ++i)
                if(subdomain_association[local_dof_indices[i]]
                   == numbers::invalid_subdomain_id)
                  subdomain_association[local_dof_indices[i]] = subdomain_id;
                else if(subdomain_association[local_dof_indices[i]]
                        > subdomain_id)
                  {
                    subdomain_association[local_dof_indices[i]] = subdomain_id;
                  }
            }

          Assert(std::find(subdomain_association.begin(),
                           subdomain_association.end(),
                           numbers::invalid_subdomain_id)
                   == subdomain_association.end(),
                 ExcInternalError());

          Assert(*std::max_element(subdomain_association.begin(),
                                   subdomain_association.end())
                   < n_procs,
                 ExcInternalError());

          return subdomain_association;
        }

        /**
         * level subdomain association. Similar to the above function only
         * for level meshes. This function assigns boundary dofs in
         * the same way as parallel::distributed::Triangulation (proc with
         * smallest index) instead of the coin flip method above.
         */
        template <class DoFHandlerType>
        std::vector<types::subdomain_id>
        get_dof_level_subdomain_association(
          const DoFHandlerType&         dof_handler,
          const types::global_dof_index n_dofs_on_level,
          const unsigned int            n_procs,
          const unsigned int            level)
        {
          (void) n_procs;
          std::vector<types::subdomain_id> level_subdomain_association(
            n_dofs_on_level, numbers::invalid_subdomain_id);
          std::vector<types::global_dof_index> local_dof_indices;
          local_dof_indices.reserve(DoFTools::max_dofs_per_cell(dof_handler));

          // loop over all cells and record which subdomain a DoF belongs to.
          // interface goes to proccessor with smaller subdomain id
          typename DoFHandlerType::cell_iterator cell
            = dof_handler.begin(level),
            endc = dof_handler.end(level);
          for(; cell != endc; ++cell)
            {
              // get the owner of the cell; note that we have made sure above that
              // all cells are either locally owned or ghosts (not artificial), so
              // this call will always yield the true owner
              const types::subdomain_id level_subdomain_id
                = cell->level_subdomain_id();
              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
              local_dof_indices.resize(dofs_per_cell);
              cell->get_mg_dof_indices(local_dof_indices);

              // set level subdomain ids. if dofs already have their values set then
              // they must be on partition interfaces. In that case assign them to the
              // processor with the smaller subdomain id.
              for(unsigned int i = 0; i < dofs_per_cell; ++i)
                if(level_subdomain_association[local_dof_indices[i]]
                   == numbers::invalid_subdomain_id)
                  level_subdomain_association[local_dof_indices[i]]
                    = level_subdomain_id;
                else if(level_subdomain_association[local_dof_indices[i]]
                        > level_subdomain_id)
                  {
                    level_subdomain_association[local_dof_indices[i]]
                      = level_subdomain_id;
                  }
            }

          Assert(std::find(level_subdomain_association.begin(),
                           level_subdomain_association.end(),
                           numbers::invalid_subdomain_id)
                   == level_subdomain_association.end(),
                 ExcInternalError());

          Assert(*std::max_element(level_subdomain_association.begin(),
                                   level_subdomain_association.end())
                   < n_procs,
                 ExcInternalError());

          return level_subdomain_association;
        }
      } // namespace

      template <class DoFHandlerType>
      NumberCache
      ParallelShared<DoFHandlerType>::distribute_dofs() const
      {
        const unsigned int dim      = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        const parallel::shared::Triangulation<dim, spacedim>* tr
          = (dynamic_cast<
             const parallel::shared::Triangulation<dim, spacedim>*>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        const unsigned int n_procs
          = Utilities::MPI::n_mpi_processes(tr->get_communicator());

        // If the underlying shared::Tria allows artificial cells,
        // then save the current set of subdomain ids, and set
        // subdomain ids to the "true" owner of each cell. we later
        // restore these flags
        std::vector<types::subdomain_id> saved_subdomain_ids;
        if(tr->with_artificial_cells())
          {
            saved_subdomain_ids.resize(tr->n_active_cells());

            typename parallel::shared::Triangulation<dim, spacedim>::
              active_cell_iterator cell
              = this->dof_handler->get_triangulation().begin_active(),
              endc = this->dof_handler->get_triangulation().end();

            const std::vector<types::subdomain_id>& true_subdomain_ids
              = tr->get_true_subdomain_ids_of_cells();

            for(unsigned int index = 0; cell != endc; ++cell, ++index)
              {
                saved_subdomain_ids[index] = cell->subdomain_id();
                cell->set_subdomain_id(true_subdomain_ids[index]);
              }
          }

        // first let the sequential algorithm do its magic. it is going to
        // enumerate DoFs on all cells, regardless of owner
        const types::global_dof_index n_dofs = Implementation::distribute_dofs(
          numbers::invalid_subdomain_id, *this->dof_handler);

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
          n_dofs, numbers::invalid_dof_index);
        {
          // first get the association of each dof with a subdomain and
          // determine the total number of subdomain ids used
          const std::vector<types::subdomain_id> subdomain_association
            = get_dof_subdomain_association(
              *this->dof_handler, n_dofs, n_procs);

          // then renumber the subdomains by first looking at those belonging
          // to subdomain 0, then those of subdomain 1, etc. note that the
          // algorithm is stable, i.e. if two dofs i,j have i<j and belong to
          // the same subdomain, then they will be in this order also after
          // reordering
          types::global_dof_index next_free_index = 0;
          for(types::subdomain_id subdomain = 0; subdomain < n_procs;
              ++subdomain)
            for(types::global_dof_index i = 0; i < n_dofs; ++i)
              if(subdomain_association[i] == subdomain)
                {
                  Assert(new_dof_indices[i] == numbers::invalid_dof_index,
                         ExcInternalError());
                  new_dof_indices[i] = next_free_index;
                  ++next_free_index;
                }

          // we should have numbered all dofs
          Assert(next_free_index == n_dofs, ExcInternalError());
          Assert(std::find(new_dof_indices.begin(),
                           new_dof_indices.end(),
                           numbers::invalid_dof_index)
                   == new_dof_indices.end(),
                 ExcInternalError());
        }
        // finally do the renumbering. we can use the sequential
        // version of the function because we do things on all
        // cells and all cells have their subdomain ids and DoFs
        // correctly set
        Implementation::renumber_dofs(
          new_dof_indices, IndexSet(0), *this->dof_handler, true);

        // update the number cache. for this, we first have to find the subdomain
        // association for each DoF again following renumbering, from which we
        // can then compute the IndexSets of locally owned DoFs for all processors.
        // all other fields then follow from this
        //
        // given the way we enumerate degrees of freedom, the locally owned
        // ranges must all be contiguous and consecutive. this makes filling
        // the IndexSets cheap. an assertion at the top verifies that this
        // assumption is true
        const std::vector<types::subdomain_id> subdomain_association
          = get_dof_subdomain_association(*this->dof_handler, n_dofs, n_procs);

        for(unsigned int i = 1; i < n_dofs; ++i)
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
          while(start_index < n_dofs)
            {
              while((end_index) < n_dofs
                    && (subdomain_association[end_index]
                        == subdomain_association[start_index]))
                ++end_index;

              // we've now identified a range of same indices. set that
              // range in the corresponding IndexSet
              if(end_index > start_index)
                {
                  const unsigned int subdomain_owner
                    = subdomain_association[start_index];
                  locally_owned_dofs_per_processor[subdomain_owner].add_range(
                    start_index, end_index);
                }

              // then move on to thinking about the next range
              start_index = end_index;
            }
        }

        // finally, restore current subdomain ids
        if(tr->with_artificial_cells())
          {
            typename parallel::shared::Triangulation<dim, spacedim>::
              active_cell_iterator cell
              = this->dof_handler->get_triangulation().begin_active(),
              endc = this->dof_handler->get_triangulation().end();

            for(unsigned int index = 0; cell != endc; ++cell, ++index)
              cell->set_subdomain_id(saved_subdomain_ids[index]);
          }

        // return a NumberCache object made up from the sets of locally
        // owned DoFs
        return NumberCache(
          locally_owned_dofs_per_processor,
          this->dof_handler->get_triangulation().locally_owned_subdomain());
      }

      template <class DoFHandlerType>
      std::vector<NumberCache>
      ParallelShared<DoFHandlerType>::distribute_mg_dofs() const
      {
        const unsigned int dim      = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        const parallel::shared::Triangulation<dim, spacedim>* tr
          = (dynamic_cast<
             const parallel::shared::Triangulation<dim, spacedim>*>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        const unsigned int n_procs
          = Utilities::MPI::n_mpi_processes(tr->get_communicator());
        const unsigned int n_levels = tr->n_global_levels();

        std::vector<NumberCache> number_caches;
        number_caches.reserve(n_levels);

        // We create an index set for each level
        for(unsigned int lvl = 0; lvl < n_levels; ++lvl)
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
              typename parallel::shared::Triangulation<dim,
                                                       spacedim>::cell_iterator
                cell
                = this->dof_handler->get_triangulation().begin(lvl),
                endc = this->dof_handler->get_triangulation().end(lvl);

              const std::vector<types::subdomain_id>& true_level_subdomain_ids
                = tr->get_true_level_subdomain_ids_of_cells(lvl);

              for(unsigned int index = 0; cell != endc; ++cell, ++index)
                {
                  saved_level_subdomain_ids[index] = cell->level_subdomain_id();
                  cell->set_level_subdomain_id(true_level_subdomain_ids[index]);
                }
            }

            // Next let the sequential algorithm do its magic. it is going to
            // enumerate DoFs on all cells on the given level, regardless of owner
            const types::global_dof_index n_dofs_on_level
              = Implementation::distribute_dofs_on_level(
                numbers::invalid_subdomain_id, *this->dof_handler, lvl);

            // then re-enumerate them based on their level subdomain association.
            // for this, we first have to identify for each current DoF
            // index which subdomain they belong to. ideally, we would
            // like to call DoFRenumbering::subdomain_wise(), but
            // because the NumberCache of the current DoFHandler is not
            // fully set up yet, we can't quite do that. also, that
            // function has to deal with other kinds of triangulations as
            // well, whereas we here know what kind of triangulation
            // we have and can simplify the code accordingly
            std::vector<types::global_dof_index> new_dof_indices(
              n_dofs_on_level, numbers::invalid_dof_index);
            {
              // first get the association of each dof with a subdomain and
              // determine the total number of subdomain ids used
              const std::vector<types::subdomain_id> level_subdomain_association
                = get_dof_level_subdomain_association(
                  *this->dof_handler, n_dofs_on_level, n_procs, lvl);

              // then renumber the subdomains by first looking at those belonging
              // to subdomain 0, then those of subdomain 1, etc. note that the
              // algorithm is stable, i.e. if two dofs i,j have i<j and belong to
              // the same subdomain, then they will be in this order also after
              // reordering
              types::global_dof_index next_free_index = 0;
              for(types::subdomain_id level_subdomain = 0;
                  level_subdomain < n_procs;
                  ++level_subdomain)
                for(types::global_dof_index i = 0; i < n_dofs_on_level; ++i)
                  if(level_subdomain_association[i] == level_subdomain)
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
                               numbers::invalid_dof_index)
                       == new_dof_indices.end(),
                     ExcInternalError());
            }

            // finally do the renumbering. we can use the sequential
            // version of the function because we do things on all
            // cells and all cells have their subdomain ids and DoFs
            // correctly set
            Implementation::renumber_mg_dofs(
              new_dof_indices, IndexSet(0), *this->dof_handler, lvl, true);

            // update the number cache. for this, we first have to find the level subdomain
            // association for each DoF again following renumbering, from which we
            // can then compute the IndexSets of locally owned DoFs for all processors.
            // all other fields then follow from this
            //
            // given the way we enumerate degrees of freedom, the locally owned
            // ranges must all be contiguous and consecutive. this makes filling
            // the IndexSets cheap. an assertion at the top verifies that this
            // assumption is true
            const std::vector<types::subdomain_id> level_subdomain_association
              = get_dof_level_subdomain_association(
                *this->dof_handler, n_dofs_on_level, n_procs, lvl);

            for(unsigned int i = 1; i < n_dofs_on_level; ++i)
              Assert(level_subdomain_association[i]
                       >= level_subdomain_association[i - 1],
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
              while(start_index < n_dofs_on_level)
                {
                  while((end_index) < n_dofs_on_level
                        && (level_subdomain_association[end_index]
                            == level_subdomain_association[start_index]))
                    ++end_index;

                  // we've now identified a range of same indices. set that
                  // range in the corresponding IndexSet
                  if(end_index > start_index)
                    {
                      const unsigned int level_subdomain_owner
                        = level_subdomain_association[start_index];
                      locally_owned_dofs_per_processor[level_subdomain_owner]
                        .add_range(start_index, end_index);
                    }

                  // then move on to thinking about the next range
                  start_index = end_index;
                }
            }

            // finally, restore current level subdomain ids
            {
              typename parallel::shared::Triangulation<dim,
                                                       spacedim>::cell_iterator
                cell
                = this->dof_handler->get_triangulation().begin(lvl),
                endc = this->dof_handler->get_triangulation().end(lvl);

              for(unsigned int index = 0; cell != endc; ++cell, ++index)
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

      template <class DoFHandlerType>
      NumberCache
      ParallelShared<DoFHandlerType>::renumber_dofs(
        const std::vector<types::global_dof_index>& new_numbers) const
      {
#ifndef DEAL_II_WITH_MPI
        (void) new_numbers;
        Assert(false, ExcNotImplemented());
        return NumberCache();
#else
        const unsigned int dim      = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        // Similar to distribute_dofs() we need to have a special treatment in
        // case artificial cells are present.
        const parallel::shared::Triangulation<dim, spacedim>* tr
          = (dynamic_cast<
             const parallel::shared::Triangulation<dim, spacedim>*>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        typename parallel::shared::Triangulation<dim,
                                                 spacedim>::active_cell_iterator
          cell
          = this->dof_handler->get_triangulation().begin_active(),
          endc = this->dof_handler->get_triangulation().end();
        std::vector<types::subdomain_id> current_subdomain_ids(
          tr->n_active_cells());
        const std::vector<types::subdomain_id>& true_subdomain_ids
          = tr->get_true_subdomain_ids_of_cells();
        if(tr->with_artificial_cells())
          for(unsigned int index = 0; cell != endc; cell++, index++)
            {
              current_subdomain_ids[index] = cell->subdomain_id();
              cell->set_subdomain_id(true_subdomain_ids[index]);
            }

        std::vector<types::global_dof_index> global_gathered_numbers(
          this->dof_handler->n_dofs(), 0);
        // as we call DoFRenumbering::subdomain_wise (*dof_handler) from distribute_dofs(),
        // we need to support sequential-like input.
        // Distributed-like input from, for example, component_wise renumbering is also supported.
        if(new_numbers.size() == this->dof_handler->n_dofs())
          {
            global_gathered_numbers = new_numbers;
          }
        else
          {
            Assert(new_numbers.size()
                     == this->dof_handler->locally_owned_dofs().n_elements(),
                   ExcInternalError());
            const unsigned int n_cpu
              = Utilities::MPI::n_mpi_processes(tr->get_communicator());
            std::vector<types::global_dof_index> gathered_new_numbers(
              this->dof_handler->n_dofs(), 0);
            Assert(Utilities::MPI::this_mpi_process(tr->get_communicator())
                     == this->dof_handler->get_triangulation()
                          .locally_owned_subdomain(),
                   ExcInternalError())

            // gather new numbers among processors into one vector
            {
              std::vector<types::global_dof_index> new_numbers_copy(
                new_numbers);

              // store the number of elements that are to be received from each process
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
              for(unsigned int i = 0; i < n_cpu; i++)
                {
                  displacements[i] = shift;
                  shift += rcounts[i];
                }
              Assert(((int) new_numbers_copy.size())
                       == rcounts[Utilities::MPI::this_mpi_process(
                            tr->get_communicator())],
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

            // put new numbers according to the current locally_owned_dofs_per_processor IndexSets
            types::global_dof_index shift = 0;
            // flag_1 and flag_2 are
            // used to control that there is a
            // one-to-one relation between old and new DoFs.
            std::vector<unsigned int> flag_1(this->dof_handler->n_dofs(), 0);
            std::vector<unsigned int> flag_2(this->dof_handler->n_dofs(), 0);
            for(unsigned int i = 0; i < n_cpu; i++)
              {
                const IndexSet iset
                  = this->dof_handler->locally_owned_dofs_per_processor()[i];
                for(types::global_dof_index ind = 0; ind < iset.n_elements();
                    ind++)
                  {
                    const types::global_dof_index target
                      = iset.nth_index_in_set(ind);
                    const types::global_dof_index value
                      = gathered_new_numbers[shift + ind];
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
        Implementation::renumber_dofs(
          global_gathered_numbers, IndexSet(0), *this->dof_handler, true);

        const NumberCache number_cache(
          DoFTools::locally_owned_dofs_per_subdomain(*this->dof_handler),
          this->dof_handler->get_triangulation().locally_owned_subdomain());

        // restore artificial cells
        cell = tr->begin_active();
        if(tr->with_artificial_cells())
          for(unsigned int index = 0; cell != endc; cell++, index++)
            cell->set_subdomain_id(current_subdomain_ids[index]);

        return number_cache;
#endif
      }

      template <class DoFHandlerType>
      NumberCache
      ParallelShared<DoFHandlerType>::renumber_mg_dofs(
        const unsigned int /*level*/,
        const std::vector<types::global_dof_index>& /*new_numbers*/) const
      {
        // multigrid is not currently implemented for shared triangulations
        Assert(false, ExcNotImplemented());

        return NumberCache();
      }

      /* --------------------- class ParallelDistributed ---------------- */

#ifdef DEAL_II_WITH_P4EST

      namespace
      {
        /**
         * A structure that allows the transfer of DoF indices from one processor
         * to another. It corresponds to a packed buffer that stores a list of
         * cells (in the form of a list of coarse mesh index -- i.e., the tree_index
         * of the cell, and a corresponding list of quadrants within these trees),
         * and a long array of DoF indices.
         *
         * The list of DoF indices stores first the number of indices for the
         * first cell (=tree index and quadrant), then the indices for that cell,
         * then the number of indices of the second cell, then the actual indices
         * of the second cell, etc.
         *
         * The DoF indices array may or may not be used by algorithms using this
         * class.
         */
        template <int dim>
        struct CellDataTransferBuffer
        {
          std::vector<unsigned int> tree_indices;
          std::vector<typename dealii::internal::p4est::types<dim>::quadrant>
                                                       quadrants;
          std::vector<dealii::types::global_dof_index> dof_numbers_and_indices;

          /**
           * Write the data of this object to a stream for the purpose of
           * serialization.
           */
          template <class Archive>
          void
          save(Archive& ar, const unsigned int /*version*/) const
          {
            // we would like to directly serialize the 'quadrants' vector,
            // but the element type is internal to p4est and does not
            // know how to serialize itself. consequently, first copy it over
            // to an array of bytes, and then serialize that
            std::vector<char> quadrants_as_chars(sizeof(quadrants[0])
                                                 * quadrants.size());
            if(quadrants_as_chars.size() > 0)
              {
                Assert(quadrants.data() != nullptr, ExcInternalError());
                std::memcpy(quadrants_as_chars.data(),
                            quadrants.data(),
                            quadrants_as_chars.size());
              }

            // now serialize everything
            ar& quadrants_as_chars& tree_indices& dof_numbers_and_indices;
          }

          /**
           * Read the data of this object from a stream for the purpose of
           * serialization. Throw away the previous content.
           */
          template <class Archive>
          void
          load(Archive& ar, const unsigned int /*version*/)
          {
            // undo the copying trick from the 'save' function
            std::vector<char> quadrants_as_chars;
            ar& quadrants_as_chars& tree_indices& dof_numbers_and_indices;

            if(quadrants_as_chars.size() > 0)
              {
                quadrants.resize(quadrants_as_chars.size()
                                 / sizeof(quadrants[0]));
                std::memcpy(quadrants.data(),
                            quadrants_as_chars.data(),
                            quadrants_as_chars.size());
              }
            else
              quadrants.clear();
          }

          BOOST_SERIALIZATION_SPLIT_MEMBER()

          /**
           * Pack the data that corresponds to this object into a buffer in
           * the form of a vector of chars and return it.
           */
          std::vector<char>
          pack_data() const
          {
            // set up a buffer and then use it as the target of a compressing
            // stream into which we serialize the current object
            std::vector<char> buffer;
            {
#  ifdef DEAL_II_WITH_ZLIB
              boost::iostreams::filtering_ostream out;
              out.push(
                boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(
                  boost::iostreams::gzip::best_compression)));
              out.push(boost::iostreams::back_inserter(buffer));

              boost::archive::binary_oarchive archive(out);

              archive << *this;
              out.flush();
#  else
              std::ostringstream              out;
              boost::archive::binary_oarchive archive(out);
              archive << *this;
              const std::string& s = out.str();
              buffer.reserve(s.size());
              buffer.assign(s.begin(), s.end());
#  endif
            }

            return buffer;
          }

          /**
           * Given a buffer in the form of an array of chars, unpack it and
           * restore the current object to the state that it was when
           * it was packed into said buffer by the pack_data() function.
           */
          void
          unpack_data(const std::vector<char>& buffer)
          {
            std::string decompressed_buffer;

            // first decompress the buffer
            {
#  ifdef DEAL_II_WITH_ZLIB
              boost::iostreams::filtering_ostream decompressing_stream;
              decompressing_stream.push(boost::iostreams::gzip_decompressor());
              decompressing_stream.push(
                boost::iostreams::back_inserter(decompressed_buffer));
              decompressing_stream.write(buffer.data(), buffer.size());
#  else
              decompressed_buffer.assign(buffer.begin(), buffer.end());
#  endif
            }

            // then restore the object from the buffer
            std::istringstream              in(decompressed_buffer);
            boost::archive::binary_iarchive archive(in);

            archive >> *this;
          }
        };

        template <int dim, int spacedim>
        void
        get_mg_dofindices_recursively(
          const parallel::distributed::Triangulation<dim, spacedim>& tria,
          const typename dealii::internal::p4est::types<dim>::quadrant&
            p4est_cell,
          const typename DoFHandler<dim, spacedim>::level_cell_iterator&
            dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant&
                                       quadrant,
          CellDataTransferBuffer<dim>& cell_data_transfer_buffer)
        {
          if(internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              // why would somebody request a cell that is not ours?
              Assert(dealii_cell->level_subdomain_id()
                       == tria.locally_owned_subdomain(),
                     ExcInternalError());

              std::vector<dealii::types::global_dof_index> local_dof_indices(
                dealii_cell->get_fe().dofs_per_cell);
              dealii_cell->get_mg_dof_indices(local_dof_indices);

              cell_data_transfer_buffer.dof_numbers_and_indices.push_back(
                dealii_cell->get_fe().dofs_per_cell);
              cell_data_transfer_buffer.dof_numbers_and_indices.insert(
                cell_data_transfer_buffer.dof_numbers_and_indices.end(),
                local_dof_indices.begin(),
                local_dof_indices.end());
              return; // we are done
            }

          if(!dealii_cell->has_children())
            return;

          if(!internal::p4est::quadrant_is_ancestor<dim>(p4est_cell, quadrant))
            return;

          typename dealii::internal::p4est::types<dim>::quadrant
            p4est_child[GeometryInfo<dim>::max_children_per_cell];
          internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);

          for(unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
              ++c)
            get_mg_dofindices_recursively<dim, spacedim>(
              tria,
              p4est_child[c],
              dealii_cell->child(c),
              quadrant,
              cell_data_transfer_buffer);
        }

        template <int dim, int spacedim>
        void
        find_marked_mg_ghost_cells_recursively(
          const typename parallel::distributed::Triangulation<dim, spacedim>&
                             tria,
          const unsigned int tree_index,
          const typename DoFHandler<dim, spacedim>::level_cell_iterator&
            dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant&
            p4est_cell,
          std::map<dealii::types::subdomain_id, CellDataTransferBuffer<dim>>&
            neighbor_cell_list)
        {
          // recurse...
          if(dealii_cell->has_children())
            {
              typename dealii::internal::p4est::types<dim>::quadrant
                p4est_child[GeometryInfo<dim>::max_children_per_cell];
              internal::p4est::init_quadrant_children<dim>(p4est_cell,
                                                           p4est_child);

              for(unsigned int c = 0;
                  c < GeometryInfo<dim>::max_children_per_cell;
                  ++c)
                find_marked_mg_ghost_cells_recursively<dim, spacedim>(
                  tria,
                  tree_index,
                  dealii_cell->child(c),
                  p4est_child[c],
                  neighbor_cell_list);
            }

          if(dealii_cell->user_flag_set()
             && dealii_cell->level_subdomain_id()
                  != tria.locally_owned_subdomain())
            {
              neighbor_cell_list[dealii_cell->level_subdomain_id()]
                .tree_indices.push_back(tree_index);
              neighbor_cell_list[dealii_cell->level_subdomain_id()]
                .quadrants.push_back(p4est_cell);
            }
        }

        template <int dim, int spacedim>
        void
        set_mg_dofindices_recursively(
          const parallel::distributed::Triangulation<dim, spacedim>& tria,
          const typename dealii::internal::p4est::types<dim>::quadrant&
            p4est_cell,
          const typename DoFHandler<dim, spacedim>::level_cell_iterator&
            dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant&
                                           quadrant,
          dealii::types::global_dof_index* dofs)
        {
          if(internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              Assert(dealii_cell->level_subdomain_id()
                       != dealii::numbers::artificial_subdomain_id,
                     ExcInternalError());

              // update dof indices of cell
              std::vector<dealii::types::global_dof_index> dof_indices(
                dealii_cell->get_fe().dofs_per_cell);
              dealii_cell->get_mg_dof_indices(dof_indices);

              bool complete = true;
              for(unsigned int i = 0; i < dof_indices.size(); ++i)
                if(dofs[i] != numbers::invalid_dof_index)
                  {
                    Assert((dof_indices[i] == (numbers::invalid_dof_index))
                             || (dof_indices[i] == dofs[i]),
                           ExcInternalError());
                    dof_indices[i] = dofs[i];
                  }
                else
                  complete = false;

              if(!complete)
                const_cast<
                  typename DoFHandler<dim, spacedim>::level_cell_iterator&>(
                  dealii_cell)
                  ->set_user_flag();
              else
                const_cast<
                  typename DoFHandler<dim, spacedim>::level_cell_iterator&>(
                  dealii_cell)
                  ->clear_user_flag();

              const_cast<
                typename DoFHandler<dim, spacedim>::level_cell_iterator&>(
                dealii_cell)
                ->set_mg_dof_indices(dof_indices);
              return;
            }

          if(!dealii_cell->has_children())
            return;

          if(!internal::p4est::quadrant_is_ancestor<dim>(p4est_cell, quadrant))
            return;

          typename dealii::internal::p4est::types<dim>::quadrant
            p4est_child[GeometryInfo<dim>::max_children_per_cell];
          internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);

          for(unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
              ++c)
            set_mg_dofindices_recursively<dim, spacedim>(
              tria, p4est_child[c], dealii_cell->child(c), quadrant, dofs);
        }

        template <int dim, int spacedim, class DoFHandlerType>
        void
        communicate_mg_ghost_cells(
          const typename parallel::distributed::Triangulation<dim, spacedim>&
                          tria,
          DoFHandlerType& dof_handler,
          const std::vector<dealii::types::global_dof_index>&
            coarse_cell_to_p4est_tree_permutation,
          const std::vector<dealii::types::global_dof_index>&
            p4est_tree_to_coarse_cell_permutation)
        {
          // build list of cells to request for each neighbor
          std::set<dealii::types::subdomain_id> level_ghost_owners
            = tria.level_ghost_owners();
          typedef std::map<dealii::types::subdomain_id,
                           CellDataTransferBuffer<dim>>
                    cellmap_t;
          cellmap_t neighbor_cell_list;
          for(std::set<dealii::types::subdomain_id>::iterator it
              = level_ghost_owners.begin();
              it != level_ghost_owners.end();
              ++it)
            neighbor_cell_list.insert(
              std::make_pair(*it, CellDataTransferBuffer<dim>()));

          for(typename DoFHandlerType::level_cell_iterator cell
              = dof_handler.begin(0);
              cell != dof_handler.end(0);
              ++cell)
            {
              typename dealii::internal::p4est::types<dim>::quadrant
                p4est_coarse_cell;
              internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

              find_marked_mg_ghost_cells_recursively<dim, spacedim>(
                tria,
                coarse_cell_to_p4est_tree_permutation[cell->index()],
                cell,
                p4est_coarse_cell,
                neighbor_cell_list);
            }
          Assert(level_ghost_owners.size() == neighbor_cell_list.size(),
                 ExcInternalError());

          //* send our requests:
          std::vector<std::vector<char>> sendbuffers(level_ghost_owners.size());
          std::vector<MPI_Request>       requests(level_ghost_owners.size());

          unsigned int idx = 0;
          for(typename cellmap_t::iterator it = neighbor_cell_list.begin();
              it != neighbor_cell_list.end();
              ++it, ++idx)
            {
              // pack all the data into the buffer for this recipient
              // and send it. keep data around till we can make sure
              // that the packet has been received
              sendbuffers[idx] = it->second.pack_data();
              const int ierr   = MPI_Isend(sendbuffers[idx].data(),
                                         sendbuffers[idx].size(),
                                         MPI_BYTE,
                                         it->first,
                                         10101,
                                         tria.get_communicator(),
                                         &requests[idx]);
              AssertThrowMPI(ierr);
            }

          //* receive requests and reply
          std::vector<std::vector<char>> reply_buffers(
            level_ghost_owners.size());
          std::vector<MPI_Request> reply_requests(level_ghost_owners.size());

          for(unsigned int idx = 0; idx < level_ghost_owners.size(); ++idx)
            {
              std::vector<char>           receive;
              CellDataTransferBuffer<dim> cell_data_transfer_buffer;

              MPI_Status status;
              int        len;
              int        ierr = MPI_Probe(
                MPI_ANY_SOURCE, 10101, tria.get_communicator(), &status);
              AssertThrowMPI(ierr);
              ierr = MPI_Get_count(&status, MPI_BYTE, &len);
              AssertThrowMPI(ierr);
              receive.resize(len);

              char* ptr = receive.data();
              ierr      = MPI_Recv(ptr,
                              len,
                              MPI_BYTE,
                              status.MPI_SOURCE,
                              status.MPI_TAG,
                              tria.get_communicator(),
                              &status);
              AssertThrowMPI(ierr);

              cell_data_transfer_buffer.unpack_data(receive);

              // store the dof indices for each cell
              for(unsigned int c = 0;
                  c < cell_data_transfer_buffer.tree_indices.size();
                  ++c)
                {
                  typename DoFHandlerType::level_cell_iterator cell(
                    &dof_handler.get_triangulation(),
                    0,
                    p4est_tree_to_coarse_cell_permutation
                      [cell_data_transfer_buffer.tree_indices[c]],
                    &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant
                    p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  get_mg_dofindices_recursively<dim, spacedim>(
                    tria,
                    p4est_coarse_cell,
                    cell,
                    cell_data_transfer_buffer.quadrants[c],
                    cell_data_transfer_buffer);
                }

              // send reply
              reply_buffers[idx] = cell_data_transfer_buffer.pack_data();
              ierr               = MPI_Isend(&(reply_buffers[idx])[0],
                               reply_buffers[idx].size(),
                               MPI_BYTE,
                               status.MPI_SOURCE,
                               10102,
                               tria.get_communicator(),
                               &reply_requests[idx]);
              AssertThrowMPI(ierr);
            }

          //* finally receive the replies
          for(unsigned int idx = 0; idx < level_ghost_owners.size(); ++idx)
            {
              std::vector<char>           receive;
              CellDataTransferBuffer<dim> cell_data_transfer_buffer;

              MPI_Status status;
              int        len;
              int        ierr = MPI_Probe(
                MPI_ANY_SOURCE, 10102, tria.get_communicator(), &status);
              AssertThrowMPI(ierr);
              ierr = MPI_Get_count(&status, MPI_BYTE, &len);
              AssertThrowMPI(ierr);
              receive.resize(len);

              char* ptr = receive.data();
              ierr      = MPI_Recv(ptr,
                              len,
                              MPI_BYTE,
                              status.MPI_SOURCE,
                              status.MPI_TAG,
                              tria.get_communicator(),
                              &status);
              AssertThrowMPI(ierr);

              cell_data_transfer_buffer.unpack_data(receive);
              if(cell_data_transfer_buffer.tree_indices.size() == 0)
                continue;

              // set the dof indices for each cell
              dealii::types::global_dof_index* dofs
                = cell_data_transfer_buffer.dof_numbers_and_indices.data();
              for(unsigned int c = 0;
                  c < cell_data_transfer_buffer.tree_indices.size();
                  ++c, dofs += 1 + dofs[0])
                {
                  typename DoFHandlerType::level_cell_iterator cell(
                    &tria,
                    0,
                    p4est_tree_to_coarse_cell_permutation
                      [cell_data_transfer_buffer.tree_indices[c]],
                    &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant
                    p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  Assert(cell->get_fe().dofs_per_cell == dofs[0],
                         ExcInternalError());

                  set_mg_dofindices_recursively<dim, spacedim>(
                    tria,
                    p4est_coarse_cell,
                    cell,
                    cell_data_transfer_buffer.quadrants[c],
                    dofs + 1);
                }
            }

          // complete all sends, so that we can safely destroy the
          // buffers.
          if(requests.size() > 0)
            {
              const int ierr = MPI_Waitall(
                requests.size(), requests.data(), MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }
          if(reply_requests.size() > 0)
            {
              const int ierr = MPI_Waitall(reply_requests.size(),
                                           reply_requests.data(),
                                           MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }
        }

        template <int spacedim>
        void
        communicate_mg_ghost_cells(
          const typename parallel::distributed::Triangulation<1, spacedim>&,
          DoFHandler<1, spacedim>&,
          const std::vector<dealii::types::global_dof_index>&,
          const std::vector<dealii::types::global_dof_index>&)
        {
          Assert(false, ExcNotImplemented());
        }

        template <int spacedim>
        void
        communicate_mg_ghost_cells(
          const typename parallel::distributed::Triangulation<1, spacedim>&,
          hp::DoFHandler<1, spacedim>&,
          const std::vector<dealii::types::global_dof_index>&,
          const std::vector<dealii::types::global_dof_index>&)
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
        template <int spacedim>
        void
        communicate_dof_indices_on_marked_cells(
          const DoFHandler<1, spacedim>&,
          const std::map<unsigned int, std::set<dealii::types::subdomain_id>>&,
          const std::vector<dealii::types::global_dof_index>&,
          const std::vector<dealii::types::global_dof_index>&)
        {
          Assert(false, ExcNotImplemented());
        }

        template <int spacedim>
        void
        communicate_dof_indices_on_marked_cells(
          const hp::DoFHandler<1, spacedim>&,
          const std::map<unsigned int, std::set<dealii::types::subdomain_id>>&,
          const std::vector<dealii::types::global_dof_index>&,
          const std::vector<dealii::types::global_dof_index>&)
        {
          Assert(false, ExcNotImplemented());
        }

        template <class DoFHandlerType>
        void
        communicate_dof_indices_on_marked_cells(
          const DoFHandlerType& dof_handler,
          const std::map<unsigned int, std::set<dealii::types::subdomain_id>>&,
          const std::vector<dealii::types::global_dof_index>&,
          const std::vector<dealii::types::global_dof_index>&)
        {
#  ifndef DEAL_II_WITH_MPI
          (void) vertices_with_ghost_neighbors;
          Assert(false, ExcNotImplemented());
#  else
          const unsigned int dim = DoFHandlerType::dimension;
          const unsigned int spacedim = DoFHandlerType::space_dimension;

          // define functions that pack data on cells that are ghost cells
          // somewhere else, and unpack data on cells where we get information
          // from elsewhere
          auto pack
            = [](const typename DoFHandlerType::active_cell_iterator& cell)
            -> boost::optional<std::vector<types::global_dof_index>> {
            Assert(cell->is_locally_owned(), ExcInternalError());

            // first see whether we need to do anything at all on this cell.
            // this is determined by whether the user_flag is set on the
            // cell that indicates that the *complete* set of DoF indices
            // has not been sent
            if(cell->user_flag_set())
              {
                // get dof indices for the current cell
                std::vector<types::global_dof_index> local_dof_indices(
                  cell->get_fe().dofs_per_cell);
                cell->get_dof_indices(local_dof_indices);

                // now see if there are dof indices that were previously
                // unknown. this can only happen in phase 1, and in
                // that case we know that the user flag must have been set
                //
                // in any case, if the cell *is* complete, we do not
                // need to send the data any more in the next phase. indicate
                // this by removing the user flag
                if(std::find(local_dof_indices.begin(),
                             local_dof_indices.end(),
                             numbers::invalid_dof_index)
                   != local_dof_indices.end())
                  {
                    Assert(cell->user_flag_set(), ExcInternalError());
                  }
                else
                  cell->clear_user_flag();

                return local_dof_indices;
              }
            else
              {
                // the fact that the user flag wasn't set means that there is
                // nothing we need to send that hasn't been sent so far.
                // so return an empty array, but also verify that indeed
                // the cell is complete
#    ifdef DEBUG
                std::vector<types::global_dof_index> local_dof_indices(
                  cell->get_fe().dofs_per_cell);
                cell->get_dof_indices(local_dof_indices);

                const bool is_complete = (std::find(local_dof_indices.begin(),
                                                    local_dof_indices.end(),
                                                    numbers::invalid_dof_index)
                                          == local_dof_indices.end());
                Assert(is_complete, ExcInternalError());
#    endif
                return boost::optional<std::vector<types::global_dof_index>>();
              }
          };

          auto unpack =
            [](const typename DoFHandlerType::active_cell_iterator& cell,
               const std::vector<types::global_dof_index>& received_dof_indices)
            -> void {
            // this function should only be called on ghost cells, and
            // on top of that, only on cells that have not been
            // completed -- which we indicate via the user flag.
            // check both
            Assert(cell->is_ghost(), ExcInternalError());
            Assert(cell->user_flag_set(), ExcInternalError());

            // if we just got an incomplete array of DoF indices, then we must
            // be in the first ghost exchange and the user flag must have been
            // set. we tested that already above.
            //
            // if we did get a complete array, then we may be in the first
            // or second ghost exchange, but in any case we need not exchange
            // another time. so delete the user flag
            const bool is_complete = (std::find(received_dof_indices.begin(),
                                                received_dof_indices.end(),
                                                numbers::invalid_dof_index)
                                      == received_dof_indices.end());
            if(is_complete)
              cell->clear_user_flag();

            // in any case, set the DoF indices on this cell. some
            // of the ones we received may still be invalid because
            // the sending processor did not know them yet, so we
            // need to merge the ones we get with those that are
            // already set here and may have already been known. for
            // those that we already know *and* get, they must obviously
            // agree
            //
            // before getting the local dof indices, we need to update the
            // cell dof indices cache because we may have set dof indices
            // on a neighboring ghost cell before this one, which may have
            // affected the dof indices we know about the current cell
            std::vector<types::global_dof_index> local_dof_indices(
              cell->get_fe().dofs_per_cell);
            cell->update_cell_dof_indices_cache();
            cell->get_dof_indices(local_dof_indices);

            for(unsigned int i = 0; i < local_dof_indices.size(); ++i)
              if(local_dof_indices[i] == numbers::invalid_dof_index)
                local_dof_indices[i] = received_dof_indices[i];
              else
                // we already know the dof index. check that there
                // is no conflict
                Assert((received_dof_indices[i] == numbers::invalid_dof_index)
                         || (received_dof_indices[i] == local_dof_indices[i]),
                       ExcInternalError());

            const_cast<typename DoFHandlerType::active_cell_iterator&>(cell)
              ->set_dof_indices(local_dof_indices);
          };

          GridTools::exchange_cell_data_to_ghosts<
            std::vector<types::global_dof_index>,
            DoFHandlerType>(dof_handler, pack, unpack);

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
          if(const auto* triangulation = dynamic_cast<
               const parallel::distributed::Triangulation<dim, spacedim>*>(
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

#endif // DEAL_II_WITH_P4EST

      template <class DoFHandlerType>
      ParallelDistributed<DoFHandlerType>::ParallelDistributed(
        DoFHandlerType& dof_handler)
        : dof_handler(&dof_handler)
      {}

      template <class DoFHandlerType>
      NumberCache
      ParallelDistributed<DoFHandlerType>::distribute_dofs() const
      {
#ifndef DEAL_II_WITH_P4EST
        Assert(false, ExcNotImplemented());
        return NumberCache();
#else
        const unsigned int dim      = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        parallel::distributed::Triangulation<dim, spacedim>* triangulation
          = (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim>*>(
            const_cast<dealii::Triangulation<dim, spacedim>*>(
              &dof_handler->get_triangulation())));
        Assert(triangulation != nullptr, ExcInternalError());

        const unsigned int n_cpus
          = Utilities::MPI::n_mpi_processes(triangulation->get_communicator());

        /*
           The following algorithm has a number of stages that are all documented
           in the paper that describes the parallel::distributed functionality:

           1/ locally enumerate dofs on locally owned cells
           2/ un-numerate those that are on interfaces with ghost
              cells and that we don't own based on the tie-breaking
              criterion; re-enumerate the remaining ones. the
              end result are that we only enumerate locally owned
              DoFs
           3/ shift indices so that each processor has a unique
              range of indices
           4/ for all locally owned cells that are ghost
              cells somewhere else, send our own DoF indices
              to the appropriate set of other processors
         */

        // --------- Phase 1: enumerate dofs on locally owned cells
        const dealii::types::global_dof_index n_initial_local_dofs
          = Implementation::distribute_dofs(
            triangulation->locally_owned_subdomain(), *dof_handler);

        // --------- Phase 2: un-numerate dofs on interfaces to ghost cells
        //                    that we don't own; re-enumerate the remaining ones

        // start with the identity permutation of indices
        std::vector<dealii::types::global_dof_index> renumbering(
          n_initial_local_dofs);
        for(dealii::types::global_dof_index i = 0; i < renumbering.size(); ++i)
          renumbering[i] = i;

        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          for(auto cell : dof_handler->active_cell_iterators())
            if(cell->is_ghost()
               && (cell->subdomain_id()
                   < triangulation->locally_owned_subdomain()))
              {
                // we found a neighboring ghost cell whose subdomain
                // is "stronger" than our own subdomain

                // delete all dofs that live there and that we have
                // previously assigned a number to (i.e. the ones on
                // the interface)
                local_dof_indices.resize(cell->get_fe().dofs_per_cell);
                cell->get_dof_indices(local_dof_indices);
                for(auto& local_dof_index : local_dof_indices)
                  if(local_dof_index != numbers::invalid_dof_index)
                    renumbering[local_dof_index] = numbers::invalid_dof_index;
              }
        }

        // make the remaining indices consecutive
        dealii::types::global_dof_index n_locally_owned_dofs = 0;
        for(auto& new_index : renumbering)
          if(new_index != numbers::invalid_dof_index)
            new_index = n_locally_owned_dofs++;

        // --------- Phase 3: shift indices so that each processor has a unique
        //                    range of indices
        std::vector<dealii::types::global_dof_index>
          n_locally_owned_dofs_per_processor(n_cpus);

        const int ierr
          = MPI_Allgather(&n_locally_owned_dofs,
                          1,
                          DEAL_II_DOF_INDEX_MPI_TYPE,
                          n_locally_owned_dofs_per_processor.data(),
                          1,
                          DEAL_II_DOF_INDEX_MPI_TYPE,
                          triangulation->get_communicator());
        AssertThrowMPI(ierr);

        const dealii::types::global_dof_index my_shift
          = std::accumulate(n_locally_owned_dofs_per_processor.begin(),
                            n_locally_owned_dofs_per_processor.begin()
                              + triangulation->locally_owned_subdomain(),
                            static_cast<dealii::types::global_dof_index>(0));
        for(auto& new_index : renumbering)
          if(new_index != numbers::invalid_dof_index)
            new_index += my_shift;

        // now re-enumerate all dofs to this shifted and condensed
        // numbering form.  we renumber some dofs as invalid, so
        // choose the nocheck-version.
        Implementation::renumber_dofs(
          renumbering, IndexSet(0), *dof_handler, false);

        // now a little bit of housekeeping
        const dealii::types::global_dof_index n_global_dofs
          = std::accumulate(n_locally_owned_dofs_per_processor.begin(),
                            n_locally_owned_dofs_per_processor.end(),
                            dealii::types::global_dof_index(0));

        std::vector<IndexSet> locally_owned_dofs_per_processor(
          n_cpus, IndexSet(n_global_dofs));
        {
          dealii::types::global_dof_index current_shift = 0;
          for(unsigned int i = 0; i < n_cpus; ++i)
            {
              locally_owned_dofs_per_processor[i].add_range(
                current_shift,
                current_shift + n_locally_owned_dofs_per_processor[i]);
              current_shift += n_locally_owned_dofs_per_processor[i];
            }
        }
        NumberCache number_cache(locally_owned_dofs_per_processor,
                                 triangulation->locally_owned_subdomain());
        Assert(
          number_cache
              .locally_owned_dofs_per_processor[triangulation
                                                  ->locally_owned_subdomain()]
              .n_elements()
            == number_cache.n_locally_owned_dofs,
          ExcInternalError());
        Assert(
          !number_cache
              .locally_owned_dofs_per_processor[triangulation
                                                  ->locally_owned_subdomain()]
              .n_elements()
            || number_cache
                   .locally_owned_dofs_per_processor
                     [triangulation->locally_owned_subdomain()]
                   .nth_index_in_set(0)
                 == my_shift,
          ExcInternalError());

        // this ends the phase where we enumerate degrees of freedom on
        // each processor. what is missing is communicating DoF indices
        // on ghost cells

        // --------- Phase 4: for all locally owned cells that are ghost
        //                    cells somewhere else, send our own DoF indices
        //                    to the appropriate set of other processors
        {
          std::vector<bool> user_flags;
          triangulation->save_user_flags(user_flags);
          triangulation->clear_user_flags();

          // figure out which cells are ghost cells on which we have
          // to exchange DoF indices
          const std::map<unsigned int, std::set<dealii::types::subdomain_id>>
            vertices_with_ghost_neighbors
            = triangulation->compute_vertices_with_ghost_neighbors();

          // mark all cells that either have to send data (locally
          // owned cells that are adjacent to ghost neighbors in some
          // way) or receive data (all ghost cells) via the user flags
          for(auto cell : dof_handler->active_cell_iterators())
            if(cell->is_locally_owned())
              {
                for(unsigned int v = 0;
                    v < GeometryInfo<dim>::vertices_per_cell;
                    ++v)
                  if(vertices_with_ghost_neighbors.find(cell->vertex_index(v))
                     != vertices_with_ghost_neighbors.end())
                    {
                      cell->set_user_flag();
                      break;
                    }
              }
            else if(cell->is_ghost())
              cell->set_user_flag();

          // Send and receive cells. After this, only the local cells
          // are marked, that received new data. This has to be
          // communicated in a second communication step.
          //
          // as explained in the 'distributed' paper, this has to be
          // done twice
          communicate_dof_indices_on_marked_cells(
            *dof_handler,
            vertices_with_ghost_neighbors,
            triangulation->coarse_cell_to_p4est_tree_permutation,
            triangulation->p4est_tree_to_coarse_cell_permutation);

          communicate_dof_indices_on_marked_cells(
            *dof_handler,
            vertices_with_ghost_neighbors,
            triangulation->coarse_cell_to_p4est_tree_permutation,
            triangulation->p4est_tree_to_coarse_cell_permutation);

          // at this point, we must have taken care of the data transfer
          // on all cells we had previously marked. verify this
#  ifdef DEBUG
          for(auto cell : dof_handler->active_cell_iterators())
            Assert(cell->user_flag_set() == false, ExcInternalError());
#  endif

          triangulation->load_user_flags(user_flags);
        }

#  ifdef DEBUG
        // check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          for(auto cell : dof_handler->active_cell_iterators())
            if(!cell->is_artificial())
              {
                local_dof_indices.resize(cell->get_fe().dofs_per_cell);
                cell->get_dof_indices(local_dof_indices);
                if(local_dof_indices.end()
                   != std::find(local_dof_indices.begin(),
                                local_dof_indices.end(),
                                numbers::invalid_dof_index))
                  {
                    if(cell->is_ghost())
                      {
                        Assert(
                          false,
                          ExcMessage("A ghost cell ended up with incomplete "
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
#endif   // DEAL_II_WITH_P4EST
      }

      template <class DoFHandlerType>
      std::vector<NumberCache>
      ParallelDistributed<DoFHandlerType>::distribute_mg_dofs() const
      {
#ifndef DEAL_II_WITH_P4EST
        Assert(false, ExcNotImplemented());
        return std::vector<NumberCache>();
#else
        const unsigned int dim      = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        parallel::distributed::Triangulation<dim, spacedim>* triangulation
          = (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim>*>(
            const_cast<dealii::Triangulation<dim, spacedim>*>(
              &dof_handler->get_triangulation())));
        Assert(triangulation != nullptr, ExcInternalError());

        AssertThrow(
          (triangulation->settings
           & parallel::distributed::Triangulation<dim, spacedim>::
               construct_multigrid_hierarchy),
          ExcMessage("Multigrid DoFs can only be distributed on a parallel "
                     "Triangulation if the flag construct_multigrid_hierarchy "
                     "is set in the constructor."));

        const unsigned int n_cpus
          = Utilities::MPI::n_mpi_processes(triangulation->get_communicator());

        // loop over all levels that exist globally (across all
        // processors), even if the current processor does not in fact
        // have any cells on that level or if the local part of the
        // Triangulation has fewer levels. we need to do this because
        // we need to communicate across all processors on all levels
        const unsigned int       n_levels = triangulation->n_global_levels();
        std::vector<NumberCache> number_caches;
        number_caches.reserve(n_levels);
        for(unsigned int level = 0; level < n_levels; ++level)
          {
            NumberCache level_number_cache;

            //* 1. distribute on own subdomain
            const unsigned int n_initial_local_dofs
              = Implementation::distribute_dofs_on_level(
                triangulation->locally_owned_subdomain(), *dof_handler, level);

            //* 2. iterate over ghostcells and kill dofs that are not
            // owned by us
            std::vector<dealii::types::global_dof_index> renumbering(
              n_initial_local_dofs);
            for(dealii::types::global_dof_index i = 0; i < renumbering.size();
                ++i)
              renumbering[i] = i;

            if(level < triangulation->n_levels())
              {
                std::vector<dealii::types::global_dof_index> local_dof_indices;

                typename DoFHandlerType::level_cell_iterator cell
                  = dof_handler->begin(level),
                  endc = dof_handler->end(level);

                for(; cell != endc; ++cell)
                  if(cell->level_subdomain_id()
                       != numbers::artificial_subdomain_id
                     && (cell->level_subdomain_id()
                         < triangulation->locally_owned_subdomain()))
                    {
                      // we found a neighboring ghost cell whose
                      // subdomain is "stronger" than our own
                      // subdomain

                      // delete all dofs that live there and that we
                      // have previously assigned a number to
                      // (i.e. the ones on the interface)
                      local_dof_indices.resize(cell->get_fe().dofs_per_cell);
                      cell->get_mg_dof_indices(local_dof_indices);
                      for(unsigned int i = 0; i < cell->get_fe().dofs_per_cell;
                          ++i)
                        if(local_dof_indices[i] != numbers::invalid_dof_index)
                          renumbering[local_dof_indices[i]]
                            = numbers::invalid_dof_index;
                    }
              }

            // TODO: make this code simpler with the new constructors of NumberCache
            // make indices consecutive
            level_number_cache.n_locally_owned_dofs = 0;
            for(std::vector<dealii::types::global_dof_index>::iterator it
                = renumbering.begin();
                it != renumbering.end();
                ++it)
              if(*it != numbers::invalid_dof_index)
                *it = level_number_cache.n_locally_owned_dofs++;

            //* 3. communicate local dofcount and shift ids to make
            // them unique
            level_number_cache.n_locally_owned_dofs_per_processor.resize(
              n_cpus);

            int ierr = MPI_Allgather(
              &level_number_cache.n_locally_owned_dofs,
              1,
              DEAL_II_DOF_INDEX_MPI_TYPE,
              &level_number_cache.n_locally_owned_dofs_per_processor[0],
              1,
              DEAL_II_DOF_INDEX_MPI_TYPE,
              triangulation->get_communicator());
            AssertThrowMPI(ierr);

            const dealii::types::global_dof_index shift = std::accumulate(
              level_number_cache.n_locally_owned_dofs_per_processor.begin(),
              level_number_cache.n_locally_owned_dofs_per_processor.begin()
                + triangulation->locally_owned_subdomain(),
              static_cast<dealii::types::global_dof_index>(0));
            for(std::vector<dealii::types::global_dof_index>::iterator it
                = renumbering.begin();
                it != renumbering.end();
                ++it)
              if(*it != numbers::invalid_dof_index)
                (*it) += shift;

            // now re-enumerate all dofs to this shifted and condensed
            // numbering form.  we renumber some dofs as invalid, so
            // choose the nocheck-version of the function
            //
            // of course there is nothing for us to renumber if the
            // level we are currently dealing with doesn't even exist
            // within the current triangulation, so skip renumbering
            // in that case
            if(level < triangulation->n_levels())
              Implementation::renumber_mg_dofs(
                renumbering, IndexSet(0), *dof_handler, level, false);

            // now a little bit of housekeeping
            level_number_cache.n_global_dofs = std::accumulate(
              level_number_cache.n_locally_owned_dofs_per_processor.begin(),
              level_number_cache.n_locally_owned_dofs_per_processor.end(),
              static_cast<dealii::types::global_dof_index>(0));

            level_number_cache.locally_owned_dofs
              = IndexSet(level_number_cache.n_global_dofs);
            level_number_cache.locally_owned_dofs.add_range(
              shift, shift + level_number_cache.n_locally_owned_dofs);
            level_number_cache.locally_owned_dofs.compress();

            // fill global_dof_indexsets
            level_number_cache.locally_owned_dofs_per_processor.resize(n_cpus);
            {
              dealii::types::global_dof_index current_shift = 0;
              for(unsigned int i = 0; i < n_cpus; ++i)
                {
                  level_number_cache.locally_owned_dofs_per_processor[i]
                    = IndexSet(level_number_cache.n_global_dofs);
                  level_number_cache.locally_owned_dofs_per_processor[i]
                    .add_range(current_shift,
                               current_shift
                                 + level_number_cache
                                     .n_locally_owned_dofs_per_processor[i]);
                  current_shift
                    += level_number_cache.n_locally_owned_dofs_per_processor[i];
                }
            }
            Assert(level_number_cache
                       .locally_owned_dofs_per_processor
                         [triangulation->locally_owned_subdomain()]
                       .n_elements()
                     == level_number_cache.n_locally_owned_dofs,
                   ExcInternalError());
            Assert(!level_number_cache
                       .locally_owned_dofs_per_processor
                         [triangulation->locally_owned_subdomain()]
                       .n_elements()
                     || level_number_cache
                            .locally_owned_dofs_per_processor
                              [triangulation->locally_owned_subdomain()]
                            .nth_index_in_set(0)
                          == shift,
                   ExcInternalError());

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
            typename DoFHandlerType::level_cell_iterator cell,
              endc = dof_handler->end();
            for(cell = dof_handler->begin(); cell != endc; ++cell)
              if(cell->level_subdomain_id()
                   != dealii::numbers::artificial_subdomain_id
                 && !cell->is_locally_owned_on_level())
                cell->set_user_flag();
          }

          // Phase 1. Request all marked cells from corresponding owners. If we
          // managed to get every DoF, remove the user_flag, otherwise we
          // will request them again in the step below.
          communicate_mg_ghost_cells(
            *triangulation,
            *dof_handler,
            triangulation->coarse_cell_to_p4est_tree_permutation,
            triangulation->p4est_tree_to_coarse_cell_permutation);

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
          communicate_mg_ghost_cells(
            *triangulation,
            *dof_handler,
            triangulation->coarse_cell_to_p4est_tree_permutation,
            triangulation->p4est_tree_to_coarse_cell_permutation);

#  ifdef DEBUG
          // make sure we have removed all flags:
          {
            typename DoFHandlerType::level_cell_iterator cell,
              endc = dof_handler->end();
            for(cell = dof_handler->begin(); cell != endc; ++cell)
              if(cell->level_subdomain_id()
                   != dealii::numbers::artificial_subdomain_id
                 && !cell->is_locally_owned_on_level())
                Assert(cell->user_flag_set() == false, ExcInternalError());
          }
#  endif

          triangulation->load_user_flags(user_flags);
        }

#  ifdef DEBUG
        // check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;
          typename DoFHandlerType::level_cell_iterator cell,
            endc = dof_handler->end();

          for(cell = dof_handler->begin(); cell != endc; ++cell)
            if(cell->level_subdomain_id()
               != dealii::numbers::artificial_subdomain_id)
              {
                local_dof_indices.resize(cell->get_fe().dofs_per_cell);
                cell->get_mg_dof_indices(local_dof_indices);
                if(local_dof_indices.end()
                   != std::find(local_dof_indices.begin(),
                                local_dof_indices.end(),
                                numbers::invalid_dof_index))
                  {
                    Assert(false, ExcMessage("not all DoFs got distributed!"));
                  }
              }
        }
#  endif // DEBUG

        return number_caches;

#endif // DEAL_II_WITH_P4EST
      }

      template <class DoFHandlerType>
      NumberCache
      ParallelDistributed<DoFHandlerType>::renumber_dofs(
        const std::vector<dealii::types::global_dof_index>& new_numbers) const
      {
        (void) new_numbers;

        Assert(new_numbers.size() == dof_handler->n_locally_owned_dofs(),
               ExcInternalError());

#ifndef DEAL_II_WITH_P4EST
        Assert(false, ExcNotImplemented());
        return NumberCache();
#else
        const unsigned int dim      = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        parallel::distributed::Triangulation<dim, spacedim>* triangulation
          = (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim>*>(
            const_cast<dealii::Triangulation<dim, spacedim>*>(
              &dof_handler->get_triangulation())));
        Assert(triangulation != nullptr, ExcInternalError());

        // First figure out the new set of locally owned DoF indices.
        // If we own no DoFs, we still need to go through this function,
        // but we can skip this calculation.
        //
        // The IndexSet::add_indices() function is substantially more
        // efficient if the set of indices is already sorted because
        // it can then insert ranges instead of individual elements.
        // consequently, pre-sort the array of new indices
        IndexSet my_locally_owned_new_dof_indices(dof_handler->n_dofs());
        if(dof_handler->n_locally_owned_dofs() > 0)
          {
            std::vector<dealii::types::global_dof_index> new_numbers_sorted
              = new_numbers;
            std::sort(new_numbers_sorted.begin(), new_numbers_sorted.end());

            my_locally_owned_new_dof_indices.add_indices(
              new_numbers_sorted.begin(), new_numbers_sorted.end());
            my_locally_owned_new_dof_indices.compress();

            Assert(my_locally_owned_new_dof_indices.n_elements()
                     == new_numbers.size(),
                   ExcInternalError());
          }

        // delete all knowledge of DoF indices that are not locally
        // owned. we do so by getting DoF indices on cells, checking
        // whether they are locally owned, if not, setting them to
        // an invalid value, and then setting them again on the current
        // cell
        //
        // DoFs we (i) know about, and (ii) don't own locally must be located
        // either on ghost cells, or on the interface between a locally
        // owned cell and a ghost cell. In any case, it is sufficient
        // to kill them only from the ghost side cell, so loop only over
        // ghost cells
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          for(auto cell : dof_handler->active_cell_iterators())
            if(cell->is_ghost())
              {
                local_dof_indices.resize(cell->get_fe().dofs_per_cell);
                cell->get_dof_indices(local_dof_indices);

                for(unsigned int i = 0; i < cell->get_fe().dofs_per_cell; ++i)
                  // delete a DoF index if it has not already been deleted
                  // (e.g., by visiting a neighboring cell, if it is on the
                  // boundary), and if we don't own it
                  if((local_dof_indices[i] != numbers::invalid_dof_index)
                     && (!dof_handler->locally_owned_dofs().is_element(
                          local_dof_indices[i])))
                    local_dof_indices[i] = numbers::invalid_dof_index;

                cell->set_dof_indices(local_dof_indices);
              }
        }

        // renumber. Skip when there is nothing to do because we own no DoF.
        if(dof_handler->locally_owned_dofs().n_elements() > 0)
          Implementation::renumber_dofs(new_numbers,
                                        dof_handler->locally_owned_dofs(),
                                        *dof_handler,
                                        false);

        // communicate newly assigned DoF indices to other processors
        // and get the same information for our own ghost cells.
        //
        // this is the same as phase 4 in the distribute_dofs() algorithm
        {
          std::vector<bool> user_flags;
          triangulation->save_user_flags(user_flags);
          triangulation->clear_user_flags();

          // mark all own cells for transfer
          for(auto cell : dof_handler->active_cell_iterators())
            if(!cell->is_artificial())
              cell->set_user_flag();

          // figure out which cells are ghost cells on which we have
          // to exchange DoF indices
          const std::map<unsigned int, std::set<dealii::types::subdomain_id>>
            vertices_with_ghost_neighbors
            = triangulation->compute_vertices_with_ghost_neighbors();

          // Send and receive cells. After this, only the local cells
          // are marked, that received new data. This has to be
          // communicated in a second communication step.
          //
          // as explained in the 'distributed' paper, this has to be
          // done twice
          communicate_dof_indices_on_marked_cells(
            *dof_handler,
            vertices_with_ghost_neighbors,
            triangulation->coarse_cell_to_p4est_tree_permutation,
            triangulation->p4est_tree_to_coarse_cell_permutation);

          communicate_dof_indices_on_marked_cells(
            *dof_handler,
            vertices_with_ghost_neighbors,
            triangulation->coarse_cell_to_p4est_tree_permutation,
            triangulation->p4est_tree_to_coarse_cell_permutation);

          triangulation->load_user_flags(user_flags);
        }

        // the last step is to update the NumberCache, including knowing which
        // processor owns which DoF index. this requires communication.
        //
        // this step is substantially more complicated than it is in
        // distribute_dofs() because the IndexSets of locally owned DoFs
        // after renumbering may not be contiguous any more. for
        // distribute_dofs() it was enough to exchange the starting
        // indices for each processor and the global number of DoFs,
        // but here we actually have to serialize the IndexSet
        // objects and shop them across the network.
        const unsigned int n_cpus
          = Utilities::MPI::n_mpi_processes(triangulation->get_communicator());
        std::vector<IndexSet> locally_owned_dofs_per_processor(
          n_cpus, IndexSet(dof_handler->n_dofs()));
        {
          // serialize our own IndexSet
          std::vector<char> my_data;
          {
#  ifdef DEAL_II_WITH_ZLIB

            boost::iostreams::filtering_ostream out;
            out.push(
              boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(
                boost::iostreams::gzip::best_compression)));
            out.push(boost::iostreams::back_inserter(my_data));

            boost::archive::binary_oarchive archive(out);

            archive << my_locally_owned_new_dof_indices;
            out.flush();
#  else
            std::ostringstream              out;
            boost::archive::binary_oarchive archive(out);
            archive << my_locally_owned_new_dof_indices;
            const std::string& s = out.str();
            my_data.reserve(s.size());
            my_data.assign(s.begin(), s.end());
#  endif
          }

          // determine maximum size of IndexSet
          const unsigned int max_size = Utilities::MPI::max(
            my_data.size(), triangulation->get_communicator());

          // as the MPI_Allgather call will be reading max_size elements, and
          // as this may be past the end of my_data, we need to increase the
          // size of the local buffer. This is filled with zeros.
          my_data.resize(max_size);

          std::vector<char> buffer(max_size * n_cpus);
          const int         ierr = MPI_Allgather(my_data.data(),
                                         max_size,
                                         MPI_BYTE,
                                         buffer.data(),
                                         max_size,
                                         MPI_BYTE,
                                         triangulation->get_communicator());
          AssertThrowMPI(ierr);

          for(unsigned int i = 0; i < n_cpus; ++i)
            if(i
               == Utilities::MPI::this_mpi_process(
                    triangulation->get_communicator()))
              locally_owned_dofs_per_processor[i]
                = my_locally_owned_new_dof_indices;
            else
              {
                // copy the data previously received into a stringstream
                // object and then read the IndexSet from it
                std::string decompressed_buffer;

                // first decompress the buffer
                {
#  ifdef DEAL_II_WITH_ZLIB

                  boost::iostreams::filtering_ostream decompressing_stream;
                  decompressing_stream.push(
                    boost::iostreams::gzip_decompressor());
                  decompressing_stream.push(
                    boost::iostreams::back_inserter(decompressed_buffer));

                  decompressing_stream.write(&buffer[i * max_size], max_size);
#  else
                  decompressed_buffer.assign(&buffer[i * max_size], max_size);
#  endif
                }

                // then restore the object from the buffer
                std::istringstream              in(decompressed_buffer);
                boost::archive::binary_iarchive archive(in);

                archive >> locally_owned_dofs_per_processor[i];
              }
        }

        return NumberCache(
          locally_owned_dofs_per_processor,
          Utilities::MPI::this_mpi_process(triangulation->get_communicator()));
#endif
      }

      template <class DoFHandlerType>
      NumberCache
      ParallelDistributed<DoFHandlerType>::renumber_mg_dofs(
        const unsigned int                          level,
        const std::vector<types::global_dof_index>& new_numbers) const
      {
        // we only implement the case where the multigrid numbers are
        // renumbered within the processor's partition, rather than the most
        // general case
        const std::vector<IndexSet>& index_sets
          = dof_handler->locally_owned_mg_dofs_per_processor(level);

        constexpr int dim      = DoFHandlerType::dimension;
        constexpr int spacedim = DoFHandlerType::space_dimension;
        const parallel::Triangulation<dim, spacedim>* tr
          = (dynamic_cast<const parallel::Triangulation<dim, spacedim>*>(
            &this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

#ifdef DEAL_II_WITH_MPI
        const unsigned int my_rank
          = Utilities::MPI::this_mpi_process(tr->get_communicator());

#  ifdef DEBUG
        for(types::global_dof_index i : new_numbers)
          {
            Assert(index_sets[my_rank].is_element(i),
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
        DoFTools::extract_locally_relevant_level_dofs(
          *dof_handler, level, relevant_dofs);
        std::vector<types::global_dof_index> ghosted_new_numbers(
          relevant_dofs.n_elements());
        {
          Utilities::MPI::Partitioner partitioner(
            index_sets[my_rank], relevant_dofs, tr->get_communicator());
          std::vector<types::global_dof_index> temp_array(
            partitioner.n_import_indices());
          const unsigned int       communication_channel = 17;
          std::vector<MPI_Request> requests;
          partitioner.export_to_ghosted_array_start(
            communication_channel,
            make_array_view(new_numbers),
            make_array_view(temp_array),
            ArrayView<types::global_dof_index>(ghosted_new_numbers.data()
                                                 + new_numbers.size(),
                                               partitioner.n_ghost_indices()),
            requests);
          partitioner.export_to_ghosted_array_finish(
            ArrayView<types::global_dof_index>(ghosted_new_numbers.data()
                                                 + new_numbers.size(),
                                               partitioner.n_ghost_indices()),
            requests);

          // we need to fill the indices of the locally owned part into the
          // new numbers array. their right position is somewhere in the
          // middle of the array, so we first copy the ghosted part from
          // smaller ranks to the front, then insert the data in the middle.
          unsigned int n_ghosts_on_smaller_ranks = 0;
          for(std::pair<unsigned int, unsigned int> t :
              partitioner.ghost_targets())
            {
              if(t.first > my_rank)
                break;
              n_ghosts_on_smaller_ranks += t.second;
            }
          if(n_ghosts_on_smaller_ranks > 0)
            {
              Assert(ghosted_new_numbers.data() != nullptr, ExcInternalError());
              std::memmove(ghosted_new_numbers.data(),
                           ghosted_new_numbers.data() + new_numbers.size(),
                           sizeof(types::global_dof_index)
                             * n_ghosts_on_smaller_ranks);
            }
          if(new_numbers.size() > 0)
            {
              Assert(new_numbers.data() != nullptr, ExcInternalError());
              std::memcpy(ghosted_new_numbers.data()
                            + n_ghosts_on_smaller_ranks,
                          new_numbers.data(),
                          sizeof(types::global_dof_index) * new_numbers.size());
            }
        }

        // in case we do not own any of the given level (but only some remote
        // processor), we do not need to call the renumbering
        if(level < this->dof_handler->get_triangulation().n_levels())
          Implementation::renumber_mg_dofs(
            ghosted_new_numbers, relevant_dofs, *dof_handler, level, true);
#else
        (void) new_numbers;
        Assert(false, ExcNotImplemented());
#endif

        return NumberCache(
          index_sets, Utilities::MPI::this_mpi_process(tr->get_communicator()));
      }
    } // namespace Policy
  }   // namespace DoFHandlerImplementation
} // namespace internal

/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler_policy.inst"

DEAL_II_NAMESPACE_CLOSE
