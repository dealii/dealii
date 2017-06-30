// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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


//TODO [TH]: renumber DoFs for multigrid is not done yet

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler_policy.h>
#include <deal.II/fe/fe.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filter/gzip.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <set>
#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DoFHandler
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
        template <typename DoFHandlerType>
        void
        update_all_active_cell_dof_indices_caches (const DoFHandlerType &dof_handler)
        {
          typename DoFHandlerType::active_cell_iterator
          beginc = dof_handler.begin_active(),
          endc   = dof_handler.end();

          auto worker
            = [] (const typename DoFHandlerType::active_cell_iterator &cell,
                  void *,
                  void *)
          {
            if (!cell->is_artificial())
              cell->update_cell_dof_indices_cache ();
          };

          // parallelize filling all of the cell caches. by using
          // WorkStream, we make sure that we only run through the
          // range of iterators once, whereas a parallel_for loop
          // for example has to split the range multiple times,
          // which is expensive because cell iterators are not
          // random access iterators with a cheap operator-
          WorkStream::run (beginc, endc,
                           worker,
                           /* copier */ std::function<void (void *)>(),
                           /* scratch_data */ nullptr,
                           /* copy_data */ nullptr,
                           2*MultithreadInfo::n_threads(),
                           /* chunk_size = */ 32);
        }


        /**
         * Update the cache used for cell dof indices on all (non-artificial)
         * level (multigrid) cells of the given DoFHandler.
         */
        template <typename DoFHandlerType>
        void
        update_all_level_cell_dof_indices_caches (const DoFHandlerType &dof_handler)
        {
          typename DoFHandlerType::level_cell_iterator
          beginc = dof_handler.begin(),
          endc   = dof_handler.end();

          auto worker
            = [] (const typename DoFHandlerType::level_cell_iterator &cell,
                  void *,
                  void *)
          {
            cell->update_cell_dof_indices_cache ();
          };

          // parallelize filling all of the cell caches. by using
          // WorkStream, we make sure that we only run through the
          // range of iterators once, whereas a parallel_for loop
          // for example has to split the range multiple times,
          // which is expensive because cell iterators are not
          // random access iterators with a cheap operator-
          WorkStream::run (beginc, endc,
                           worker,
                           /* copier */ std::function<void (void *)>(),
                           /* scratch_data */ nullptr,
                           /* copy_data */ nullptr,
                           2*MultithreadInfo::n_threads(),
                           /* chunk_size = */ 32);
        }
      }



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
        static
        types::global_dof_index
        distribute_dofs_on_cell (const DoFHandler<1,spacedim>                                &dof_handler,
                                 const typename DoFHandler<1,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index                                      next_free_dof)
        {

          // distribute dofs of vertices
          if (dof_handler.get_fe().dofs_per_vertex > 0)
            for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
              {
                if (cell->vertex_dof_index (v,0) == numbers::invalid_dof_index)
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                    {
                      Assert ((cell->vertex_dof_index (v,d) ==
                               numbers::invalid_dof_index),
                              ExcInternalError());
                      cell->set_vertex_dof_index (v, d, next_free_dof++);
                    }
                else
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                    Assert ((cell->vertex_dof_index (v,d) !=
                             numbers::invalid_dof_index),
                            ExcInternalError());
              }

          // dofs of line
          for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_line; ++d)
            cell->set_dof_index (d, next_free_dof++);

          return next_free_dof;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const DoFHandler<2,spacedim>                                &dof_handler,
                                 const typename DoFHandler<2,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index                                      next_free_dof)
        {
          if (dof_handler.get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
              // check whether dofs for this vertex have been distributed
              // (checking the first dof should be good enough)
              if (cell->vertex_dof_index(vertex, 0) == numbers::invalid_dof_index)
                for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                  cell->set_vertex_dof_index (vertex, d, next_free_dof++);

          // for the four sides
          if (dof_handler.get_fe().dofs_per_line > 0)
            for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
              {
                const typename DoFHandler<2,spacedim>::line_iterator
                line = cell->line(side);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (checking the first dof should be good enough)
                if (line->dof_index(0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_line; ++d)
                    line->set_dof_index (d, next_free_dof++);
              }


          // dofs of quad
          if (dof_handler.get_fe().dofs_per_quad > 0)
            for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_quad; ++d)
              cell->set_dof_index (d, next_free_dof++);

          return next_free_dof;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const DoFHandler<3,spacedim>                                &dof_handler,
                                 const typename DoFHandler<3,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index                                      next_free_dof)
        {
          if (dof_handler.get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
              // check whether dofs for this vertex have been distributed
              // (checking the first dof should be good enough)
              if (cell->vertex_dof_index(vertex, 0) == numbers::invalid_dof_index)
                for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                  cell->set_vertex_dof_index (vertex, d, next_free_dof++);

          // for the lines
          if (dof_handler.get_fe().dofs_per_line > 0)
            for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
              {
                const typename DoFHandler<3,spacedim>::line_iterator
                line = cell->line(l);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (checking the first dof should be good enough)
                if (line->dof_index(0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_line; ++d)
                    line->set_dof_index (d, next_free_dof++);
              }

          // for the quads
          if (dof_handler.get_fe().dofs_per_quad > 0)
            for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
              {
                const typename DoFHandler<3,spacedim>::quad_iterator
                quad = cell->quad(q);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (checking the first dof should be good enough)
                if (quad->dof_index(0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_quad; ++d)
                    quad->set_dof_index (d, next_free_dof++);
              }


          // dofs of hex
          if (dof_handler.get_fe().dofs_per_hex > 0)
            for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_hex; ++d)
              cell->set_dof_index (d, next_free_dof++);

          return next_free_dof;
        }



        // same for the hp::DoFHandler
        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const hp::DoFHandler<1,spacedim>                                &dof_handler,
                                 const typename hp::DoFHandler<1,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index                                          next_free_dof)
        {
          (void)dof_handler;
          (void)cell;
          (void)next_free_dof;
          return 0;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const hp::DoFHandler<2,spacedim>                                &dof_handler,
                                 const typename hp::DoFHandler<2,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index                                          next_free_dof)
        {
          (void)dof_handler;
          (void)cell;
          (void)next_free_dof;
          return 0;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const hp::DoFHandler<3,spacedim>                                &dof_handler,
                                 const typename hp::DoFHandler<3,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index                                          next_free_dof)
        {
          (void)dof_handler;
          (void)cell;
          (void)next_free_dof;
          return 0;
        }



        /**
         * Distribute degrees of freedom on all cells, or on cells with the
         * correct subdomain_id if the corresponding argument is not equal to
         * numbers::invalid_subdomain_id. Return the total number of dofs
         * distributed.
         */
        template <typename DoFHandlerType>
        static
        types::global_dof_index
        distribute_dofs (const types::subdomain_id     subdomain_id,
                         DoFHandlerType               &dof_handler)
        {
          Assert (dof_handler.get_triangulation().n_levels() > 0,
                  ExcMessage("Empty triangulation"));

          types::global_dof_index next_free_dof = 0;
          typename DoFHandlerType::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();

          for (; cell != endc; ++cell)
            if ((subdomain_id == numbers::invalid_subdomain_id)
                ||
                (cell->subdomain_id() == subdomain_id))
              next_free_dof
                = Implementation::distribute_dofs_on_cell (dof_handler,
                                                           cell,
                                                           next_free_dof);

          update_all_active_cell_dof_indices_caches (dof_handler);

          return next_free_dof;
        }


        /**
         * Distribute multilevel dofs on the given cell, with new dofs starting
         * with index @p next_free_dof. Return the next unused index number.
         *
         * Note that unlike for the usual dofs, here all cells and not
         * only active ones are allowed.
         */
        template <int dim, int spacedim>
        static
        types::global_dof_index
        distribute_mg_dofs_on_cell (const typename DoFHandler<dim,spacedim>::level_cell_iterator &cell,
                                    types::global_dof_index   next_free_dof,
                                    const internal::int2type<1> &)
        {
          // distribute dofs of vertices
          if (cell->get_fe().dofs_per_vertex > 0)
            for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
              {
                typename DoFHandler<dim,spacedim>::level_cell_iterator neighbor = cell->neighbor(v);

                if (neighbor.state() == IteratorState::valid)
                  {
                    // has neighbor already been processed?
                    if (neighbor->user_flag_set() &&
                        (neighbor->level() == cell->level()))
                      // copy dofs if the neighbor is on the same level (only then are
                      // mg dofs the same)
                      {
                        if (v==0)
                          for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                            cell->set_mg_vertex_dof_index (cell->level(), 0, d,
                                                           neighbor->mg_vertex_dof_index (cell->level(), 1, d));
                        else
                          for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                            cell->set_mg_vertex_dof_index (cell->level(), 1, d,
                                                           neighbor->mg_vertex_dof_index (cell->level(), 0, d));

                        // next neighbor
                        continue;
                      }
                  }

                // otherwise: create dofs newly
                for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index (cell->level(), v, d, next_free_dof++);
              }

          // dofs of line
          if (cell->get_fe().dofs_per_line > 0)
            for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
              cell->set_mg_dof_index (cell->level(), d, next_free_dof++);

          // note that this cell has been processed
          cell->set_user_flag ();

          return next_free_dof;
        }



        template <int dim, int spacedim>
        static
        types::global_dof_index
        distribute_mg_dofs_on_cell (const typename DoFHandler<dim,spacedim>::level_cell_iterator &cell,
                                    types::global_dof_index   next_free_dof,
                                    const internal::int2type<2> &)
        {
          if (cell->get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
              // check whether dofs for this
              // vertex have been distributed
              // (only check the first dof)
              if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == numbers::invalid_dof_index)
                for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);

          // for the four sides
          if (cell->get_fe().dofs_per_line > 0)
            for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
              {
                typename DoFHandler<dim,spacedim>::line_iterator line = cell->line(side);

                // distribute dofs if necessary: check whether line dof is already
                // numbered (check only first dof)
                if (line->mg_dof_index(cell->level(), 0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
                    line->set_mg_dof_index (cell->level(), d, next_free_dof++);
              }


          // dofs of quad
          if (cell->get_fe().dofs_per_quad > 0)
            for (unsigned int d=0; d<cell->get_fe().dofs_per_quad; ++d)
              cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


          // note that this cell has been processed
          cell->set_user_flag ();

          return next_free_dof;
        }



        template <int dim, int spacedim>
        static
        types::global_dof_index
        distribute_mg_dofs_on_cell (const typename DoFHandler<dim,spacedim>::level_cell_iterator &cell,
                                    types::global_dof_index   next_free_dof,
                                    const internal::int2type<3> &)
        {
          if (cell->get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
              // check whether dofs for this vertex have been distributed
              // (only check the first dof)
              if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == numbers::invalid_dof_index)
                for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);

          // for the lines
          if (cell->get_fe().dofs_per_line > 0)
            for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
              {
                typename DoFHandler<dim,spacedim>::line_iterator line = cell->line(l);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if (line->mg_dof_index(cell->level(), 0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
                    line->set_mg_dof_index (cell->level(), d, next_free_dof++);
              }

          // for the quads
          if (cell->get_fe().dofs_per_quad > 0)
            for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
              {
                typename DoFHandler<dim,spacedim>::quad_iterator quad = cell->quad(q);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if (quad->mg_dof_index(cell->level(), 0) == numbers::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<cell->get_fe().dofs_per_quad; ++d)
                    quad->set_mg_dof_index (cell->level(), d, next_free_dof++);
              }


          // dofs of cell
          if (cell->get_fe().dofs_per_hex > 0)
            for (unsigned int d=0; d<cell->get_fe().dofs_per_hex; ++d)
              cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


          // note that this cell has been processed
          cell->set_user_flag ();

          return next_free_dof;
        }



        // same for the hp::DoFHandler
        template <int spacedim>
        static
        types::global_dof_index
        distribute_mg_dofs_on_cell (const hp::DoFHandler<1,spacedim>                                &dof_handler,
                                    const typename hp::DoFHandler<1,spacedim>::active_cell_iterator &cell,
                                    types::global_dof_index                                          next_free_dof)
        {
          (void)dof_handler;
          (void)cell;
          (void)next_free_dof;
          return 0;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_mg_dofs_on_cell (const hp::DoFHandler<2,spacedim>                                &dof_handler,
                                    const typename hp::DoFHandler<2,spacedim>::active_cell_iterator &cell,
                                    types::global_dof_index                                          next_free_dof)
        {
          (void)dof_handler;
          (void)cell;
          (void)next_free_dof;
          return 0;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_mg_dofs_on_cell (const hp::DoFHandler<3,spacedim>                                &dof_handler,
                                    const typename hp::DoFHandler<3,spacedim>::active_cell_iterator &cell,
                                    types::global_dof_index                                          next_free_dof)
        {
          (void)dof_handler;
          (void)cell;
          (void)next_free_dof;
          return 0;
        }



        template <typename DoFHandlerType>
        static
        types::global_dof_index
        distribute_dofs_on_level (const types::subdomain_id level_subdomain_id,
                                  DoFHandlerType           &dof_handler,
                                  const unsigned int        level)
        {
          const unsigned int dim      = DoFHandlerType::dimension;
          const unsigned int spacedim = DoFHandlerType::space_dimension;

          const dealii::Triangulation<dim,spacedim> &tria
            = dof_handler.get_triangulation();
          Assert (tria.n_levels() > 0, ExcMessage("Empty triangulation"));
          if (level>=tria.n_levels())
            return 0; //this is allowed for multigrid

          // Clear user flags because we will need them. But first we save
          // them and make sure that we restore them later such that at
          // the end of this function the Triangulation will be in the
          // same state as it was at the beginning of this function.
          std::vector<bool> user_flags;
          tria.save_user_flags(user_flags);
          const_cast<dealii::Triangulation<dim,spacedim> &>(tria).clear_user_flags ();

          types::global_dof_index next_free_dof = 0;
          typename DoFHandler<dim,spacedim>::level_cell_iterator
          cell = dof_handler.begin(level),
          endc = dof_handler.end(level);

          for (; cell != endc; ++cell)
            if ((level_subdomain_id == numbers::invalid_subdomain_id)
                ||
                (cell->level_subdomain_id() == level_subdomain_id))
              next_free_dof
                = Implementation::distribute_mg_dofs_on_cell<dim,spacedim> (cell, next_free_dof,
                    internal::int2type<dim>());

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim,spacedim> &>(tria).load_user_flags(user_flags);

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
        static
        void
        renumber_vertex_dofs (const std::vector<types::global_dof_index> &new_numbers,
                              const IndexSet                             &indices,
                              DoFHandler<dim,spacedim>                   &dof_handler,
                              const bool                                  check_validity)
        {
          // we can not use cell iterators in this function since then
          // we would renumber the dofs on the interface of two cells
          // more than once. Anyway, this way it's not only more
          // correct but also faster; note, however, that dof numbers
          // may be invalid_dof_index, namely when the appropriate
          // vertex/line/etc is unused
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.vertex_dofs.begin();
               i!=dof_handler.vertex_dofs.end(); ++i)
            if (*i != numbers::invalid_dof_index)
              *i = (indices.size() == 0)?
                   (new_numbers[*i]) :
                   (new_numbers[indices.index_within_set(*i)]);
            else if (check_validity)
              // if index is invalid_dof_index: check if this one
              // really is unused
              Assert (dof_handler.get_triangulation()
                      .vertex_used((i-dof_handler.vertex_dofs.begin()) /
                                   dof_handler.selected_fe->dofs_per_vertex)
                      == false,
                      ExcInternalError ());
        }



        /**
         * The part of the renumber_dofs() functionality that is dimension
         * independent because it renumbers the DoF indices on cell interiors
         * (which exist for all dimensions).
         *
         * See renumber_dofs() for the meaning of the arguments.
         */
        template <int dim, int spacedim>
        static
        void
        renumber_cell_dofs (const std::vector<types::global_dof_index> &new_numbers,
                            const IndexSet                             &indices,
                            DoFHandler<dim,spacedim>                   &dof_handler)
        {
          for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
            for (std::vector<types::global_dof_index>::iterator
                 i=dof_handler.levels[level]->dof_object.dofs.begin();
                 i!=dof_handler.levels[level]->dof_object.dofs.end(); ++i)
              if (*i != numbers::invalid_dof_index)
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
        static
        void
        renumber_face_dofs (const std::vector<types::global_dof_index> &/*new_numbers*/,
                            const IndexSet                             &/*indices*/,
                            DoFHandler<1,spacedim>                     &/*dof_handler*/)
        {
          // nothing to do in 1d since there are no separate faces
        }



        template <int spacedim>
        static
        void
        renumber_face_dofs (const std::vector<types::global_dof_index> &new_numbers,
                            const IndexSet                             &indices,
                            DoFHandler<2,spacedim>                     &dof_handler)
        {
          // treat dofs on lines
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.faces->lines.dofs.begin();
               i!=dof_handler.faces->lines.dofs.end(); ++i)
            if (*i != numbers::invalid_dof_index)
              *i = ((indices.size() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);
        }



        template <int spacedim>
        static
        void
        renumber_face_dofs (const std::vector<types::global_dof_index> &new_numbers,
                            const IndexSet                             &indices,
                            DoFHandler<3,spacedim>                     &dof_handler)
        {
          // treat dofs on lines
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.faces->lines.dofs.begin();
               i!=dof_handler.faces->lines.dofs.end(); ++i)
            if (*i != numbers::invalid_dof_index)
              *i = ((indices.size() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);

          // treat dofs on quads
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.faces->quads.dofs.begin();
               i!=dof_handler.faces->quads.dofs.end(); ++i)
            if (*i != numbers::invalid_dof_index)
              *i = ((indices.size() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);
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
        template <int dim, int spacedim>
        static
        void
        renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
                       const IndexSet                             &indices,
                       DoFHandler<dim,spacedim>                   &dof_handler,
                       const bool                                  check_validity)
        {
          if (dim == 1)
            Assert (indices == IndexSet(0), ExcNotImplemented());

          // renumber DoF indices on vertices, cells, and faces. this
          // can be done in parallel because the respective functions
          // work on separate data structures
          Threads::TaskGroup<> tasks;
          tasks += Threads::new_task ([&]()
          {
            renumber_vertex_dofs (new_numbers, indices, dof_handler, check_validity);
          });
          tasks += Threads::new_task ([&]()
          {
            renumber_face_dofs (new_numbers, indices, dof_handler);
          });
          tasks += Threads::new_task ([&]()
          {
            renumber_cell_dofs (new_numbers, indices, dof_handler);
          });
          tasks.join_all ();

          // update the cache used for cell dof indices
          update_all_level_cell_dof_indices_caches(dof_handler);
        }



        template <int dim, int spacedim>
        static
        void
        renumber_dofs (const std::vector<types::global_dof_index> &/*new_numbers*/,
                       const IndexSet                             &/*indices*/,
                       hp::DoFHandler<dim,spacedim>               &/*dof_handler*/,
                       const bool                                  /*check_validity*/)
        {
          Assert (false, ExcNotImplemented());
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
        static
        void
        renumber_vertex_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                                 const IndexSet                                     &indices,
                                 DoFHandler<dim,spacedim>                           &dof_handler,
                                 const unsigned int                                  level,
                                 const bool                                          check_validity)
        {
          Assert (level<dof_handler.get_triangulation().n_levels(),
                  ExcInternalError());

          for (typename std::vector<typename DoFHandler<dim,spacedim>::MGVertexDoFs>::iterator
               i=dof_handler.mg_vertex_dofs.begin();
               i!=dof_handler.mg_vertex_dofs.end();
               ++i)
            // if the present vertex lives on the current level
            if ((i->get_coarsest_level() <= level) &&
                (i->get_finest_level() >= level))
              for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                {
                  const dealii::types::global_dof_index idx = i->get_index (level, d);
                  if (check_validity)
                    Assert(idx != numbers::invalid_dof_index, ExcInternalError ());

                  if (idx != numbers::invalid_dof_index)
                    i->set_index (level, d,
                                  (indices.size() == 0)?
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
        static
        void
        renumber_cell_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                               const IndexSet                                     &indices,
                               DoFHandler<dim,spacedim>                           &dof_handler,
                               const unsigned int                                  level)
        {
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.mg_levels[level]->dof_object.dofs.begin();
               i!=dof_handler.mg_levels[level]->dof_object.dofs.end();
               ++i)
            {
              if (*i != numbers::invalid_dof_index)
                {
                  Assert(*i<new_numbers.size(), ExcInternalError());
                  *i = (indices.size() == 0)?
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
        static
        void
        renumber_face_mg_dofs (const std::vector<types::global_dof_index> &/*new_numbers*/,
                               const IndexSet                             &/*indices*/,
                               DoFHandler<1,spacedim>                     &/*dof_handler*/,
                               const unsigned int                          /*level*/,
                               const bool                                  /*check_validity*/)
        {
          // nothing to do in 1d because there are no separate faces
        }



        template <int spacedim>
        static
        void
        renumber_face_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                               const IndexSet                                     &indices,
                               DoFHandler<2,spacedim>                             &dof_handler,
                               const unsigned int                                  level,
                               const bool                                          check_validity)
        {
          if (dof_handler.get_fe().dofs_per_line > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<2,spacedim> &>(dof_handler.get_triangulation()).clear_user_flags ();

              // flag all lines adjacent to cells of the current
              // level, as those lines logically belong to the same
              // level as the cell, at least for for isotropic
              // refinement
              typename DoFHandler<2,spacedim>::level_cell_iterator cell,
                       endc = dof_handler.end(level);
              for (cell = dof_handler.begin(level); cell != endc; ++cell)
                for (unsigned int line=0; line < GeometryInfo<2>::faces_per_cell; ++line)
                  cell->face(line)->set_user_flag();

              for (typename DoFHandler<2,spacedim>::cell_iterator cell = dof_handler.begin();
                   cell != dof_handler.end(); ++cell)
                for (unsigned int l=0; l<GeometryInfo<2>::lines_per_cell; ++l)
                  if (cell->line(l)->user_flag_set())
                    {
                      for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_line; ++d)
                        {
                          const dealii::types::global_dof_index idx = cell->line(l)->mg_dof_index(level, d);
                          if (check_validity)
                            Assert(idx != numbers::invalid_dof_index, ExcInternalError ());

                          if (idx != numbers::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index (level, d, ((indices.size() == 0) ?
                                                                        new_numbers[idx] :
                                                                        new_numbers[indices.index_within_set(idx)]));
                        }
                      cell->line(l)->clear_user_flag();
                    }
              // finally, restore user flags
              const_cast<dealii::Triangulation<2,spacedim> &>(dof_handler.get_triangulation()).load_user_flags (user_flags);
            }
        }



        template <int spacedim>
        static
        void
        renumber_face_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                               const IndexSet                                     &indices,
                               DoFHandler<3,spacedim>                             &dof_handler,
                               const unsigned int                                  level,
                               const bool                                          check_validity)
        {
          if (dof_handler.get_fe().dofs_per_line > 0 ||
              dof_handler.get_fe().dofs_per_quad > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<3,spacedim> &>(dof_handler.get_triangulation()).clear_user_flags ();

              // flag all lines adjacent to cells of the current
              // level, as those lines logically belong to the same
              // level as the cell, at least for isotropic refinement
              typename DoFHandler<3,spacedim>::level_cell_iterator cell,
                       endc = dof_handler.end(level);
              for (cell = dof_handler.begin(level); cell != endc; ++cell)
                for (unsigned int line=0; line < GeometryInfo<3>::lines_per_cell; ++line)
                  cell->line(line)->set_user_flag();

              for (typename DoFHandler<3,spacedim>::cell_iterator cell = dof_handler.begin();
                   cell != dof_handler.end(); ++cell)
                for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
                  if (cell->line(l)->user_flag_set())
                    {
                      for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_line; ++d)
                        {
                          const dealii::types::global_dof_index idx = cell->line(l)->mg_dof_index(level, d);
                          if (check_validity)
                            Assert(idx != numbers::invalid_dof_index, ExcInternalError ());

                          if (idx != numbers::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index (level, d, ((indices.size() == 0) ?
                                                                        new_numbers[idx] :
                                                                        new_numbers[indices.index_within_set(idx)]));
                        }
                      cell->line(l)->clear_user_flag();
                    }

              // flag all quads adjacent to cells of the current level, as
              // those quads logically belong to the same level as the cell,
              // at least for isotropic refinement
              for (cell = dof_handler.begin(level); cell != endc; ++cell)
                for (unsigned int quad=0; quad < GeometryInfo<3>::quads_per_cell; ++quad)
                  cell->quad(quad)->set_user_flag();

              for (typename DoFHandler<3,spacedim>::cell_iterator cell = dof_handler.begin();
                   cell != dof_handler.end(); ++cell)
                for (unsigned int l=0; l<GeometryInfo<3>::quads_per_cell; ++l)
                  if (cell->quad(l)->user_flag_set())
                    {
                      for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_quad; ++d)
                        {
                          const dealii::types::global_dof_index idx = cell->quad(l)->mg_dof_index(level, d);
                          if (check_validity)
                            Assert(idx != numbers::invalid_dof_index, ExcInternalError ());

                          if (idx != numbers::invalid_dof_index)
                            cell->quad(l)->set_mg_dof_index (level, d, ((indices.size() == 0) ?
                                                                        new_numbers[idx] :
                                                                        new_numbers[indices.index_within_set(idx)]));
                        }
                      cell->quad(l)->clear_user_flag();
                    }

              // finally, restore user flags
              const_cast<dealii::Triangulation<3,spacedim> &>(dof_handler.get_triangulation()).load_user_flags (user_flags);
            }
        }



        template <int dim, int spacedim>
        static
        void
        renumber_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                          const IndexSet                                     &indices,
                          DoFHandler<dim,spacedim>                           &dof_handler,
                          const unsigned int                                  level,
                          const bool                                          check_validity)
        {
          Assert (level<dof_handler.get_triangulation().n_levels(),
                  ExcInternalError());

          // renumber DoF indices on vertices, cells, and faces. this
          // can be done in parallel because the respective functions
          // work on separate data structures
          Threads::TaskGroup<> tasks;
          tasks += Threads::new_task ([&]()
          {
            renumber_vertex_mg_dofs (new_numbers, indices, dof_handler, level, check_validity);
          });
          tasks += Threads::new_task ([&]()
          {
            renumber_face_mg_dofs (new_numbers, indices, dof_handler, level, check_validity);
          });
          tasks += Threads::new_task ([&]()
          {
            renumber_cell_mg_dofs (new_numbers, indices, dof_handler, level);
          });
          tasks.join_all ();
        }



        template <int dim, int spacedim>
        static
        void
        renumber_mg_dofs (const std::vector<dealii::types::global_dof_index> &/*new_numbers*/,
                          const IndexSet                                     &/*indices*/,
                          hp::DoFHandler<dim,spacedim>                       &/*dof_handler*/,
                          const unsigned int                                  /*level*/,
                          const bool                                          /*check_validity*/)
        {
          Assert (false, ExcNotImplemented());
        }

      };



      /* --------------------- class Sequential ---------------- */



      template <typename DoFHandlerType>
      Sequential<DoFHandlerType>::
      Sequential (DoFHandlerType &dof_handler)
        :
        dof_handler (&dof_handler)
      {}



      template <typename DoFHandlerType>
      NumberCache
      Sequential<DoFHandlerType>::
      distribute_dofs () const
      {
        const types::global_dof_index n_dofs =
          Implementation::distribute_dofs (numbers::invalid_subdomain_id,
                                           *dof_handler);

        // return a sequential, complete index set
        return NumberCache (n_dofs);
      }



      template <typename DoFHandlerType>
      std::vector<NumberCache>
      Sequential<DoFHandlerType>::
      distribute_mg_dofs () const
      {
        std::vector<bool> user_flags;
        dof_handler->get_triangulation().save_user_flags (user_flags);

        const_cast<dealii::Triangulation<DoFHandlerType::dimension, DoFHandlerType::space_dimension>&>
        (dof_handler->get_triangulation()).clear_user_flags ();

        std::vector<NumberCache> number_caches;
        number_caches.reserve (dof_handler->get_triangulation().n_levels());
        for (unsigned int level = 0; level < dof_handler->get_triangulation().n_levels(); ++level)
          {
            // first distribute dofs on this level
            const types::global_dof_index n_level_dofs
              = Implementation::distribute_dofs_on_level(numbers::invalid_subdomain_id,
                                                         *dof_handler, level);

            // then add a complete, sequential index set
            number_caches.emplace_back (NumberCache(n_level_dofs));
          }

        const_cast<dealii::Triangulation<DoFHandlerType::dimension, DoFHandlerType::space_dimension>&>
        (dof_handler->get_triangulation()).load_user_flags (user_flags);

        return number_caches;
      }



      template <typename DoFHandlerType>
      NumberCache
      Sequential<DoFHandlerType>::
      renumber_dofs (const std::vector<types::global_dof_index> &new_numbers) const
      {
        Implementation::renumber_dofs (new_numbers, IndexSet(0),
                                       *dof_handler, true);

        // return a sequential, complete index set
        return NumberCache (new_numbers.size());
      }


      /* --------------------- class ParallelShared ---------------- */


      template <int dim, int spacedim>
      ParallelShared<dim,spacedim>::
      ParallelShared (dealii::DoFHandler<dim,spacedim> &dof_handler)
        :
        dof_handler (&dof_handler)
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
        get_dof_subdomain_association (const DoFHandlerType          &dof_handler,
                                       const types::global_dof_index  n_dofs,
                                       const unsigned int             n_procs)
        {
          (void)n_procs;
          std::vector<types::subdomain_id> subdomain_association (n_dofs,
                                                                  numbers::invalid_subdomain_id);
          std::vector<types::global_dof_index> local_dof_indices;
          local_dof_indices.reserve (DoFTools::max_dofs_per_cell(dof_handler));

          // pseudo-randomly assign variables which lie on the interface between
          // subdomains to each of the two or more
          bool coin_flip = true;

          // loop over all cells and record which subdomain a DoF belongs to.
          // toss a coin in case it is on an interface
          typename DoFHandlerType::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell!=endc; ++cell)
            {
              // get the owner of the cell; note that we have made sure above that
              // all cells are either locally owned or ghosts (not artificial), so
              // this call will always yield the true owner
              const types::subdomain_id subdomain_id = cell->subdomain_id();
              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
              local_dof_indices.resize (dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);

              // set subdomain ids. if dofs already have their values set then
              // they must be on partition interfaces. in that case randomly
              // assign them to either the previous association or the current
              // one, where we take "random" to be "once this way once that way"
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                if (subdomain_association[local_dof_indices[i]] ==
                    numbers::invalid_subdomain_id)
                  subdomain_association[local_dof_indices[i]] = subdomain_id;
                else
                  {
                    if (coin_flip == true)
                      subdomain_association[local_dof_indices[i]] = subdomain_id;
                    coin_flip = !coin_flip;
                  }
            }

          Assert (std::find (subdomain_association.begin(),
                             subdomain_association.end(),
                             numbers::invalid_subdomain_id)
                  == subdomain_association.end(),
                  ExcInternalError());

          Assert (*std::max_element (subdomain_association.begin(),
                                     subdomain_association.end())
                  < n_procs,
                  ExcInternalError());

          return subdomain_association;
        }
      }



      template <int dim, int spacedim>
      NumberCache
      ParallelShared<dim,spacedim>::
      distribute_dofs () const
      {
        const parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>*> (&this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        const unsigned int n_procs = Utilities::MPI::n_mpi_processes(tr->get_communicator());

        // If the underlying shared::Tria allows artificial cells,
        // then save the current set of subdomain ids, and set
        // subdomain ids to the "true" owner of each cell. we later
        // restore these flags
        std::vector<types::subdomain_id> saved_subdomain_ids;
        if (tr->with_artificial_cells())
          {
            saved_subdomain_ids.resize (tr->n_active_cells());

            typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator
            cell = this->dof_handler->get_triangulation().begin_active(),
            endc = this->dof_handler->get_triangulation().end();

            const std::vector<types::subdomain_id> &true_subdomain_ids
              = tr->get_true_subdomain_ids_of_cells();

            for (unsigned int index=0; cell != endc; ++cell, ++index)
              {
                saved_subdomain_ids[index] = cell->subdomain_id();
                cell->set_subdomain_id(true_subdomain_ids[index]);
              }
          }

        // first let the sequential algorithm do its magic. it is going to
        // enumerate DoFs on all cells, regardless of owner
        const types::global_dof_index n_dofs =
          Implementation::distribute_dofs (numbers::invalid_subdomain_id,
                                           *this->dof_handler);

        // then re-enumerate them based on their subdomain association.
        // for this, we first have to identify for each current DoF
        // index which subdomain they belong to. ideally, we would
        // like to call DoFRenumbering::subdomain_wise(), but
        // because the NumberCache of the current DoFHandler is not
        // fully set up yet, we can't quite do that. also, that
        // function has to deal with other kinds of triangulations as
        // well, whereas we here know what kind of triangulation
        // we have and can simplify the code accordingly
        std::vector<types::global_dof_index> new_dof_indices (n_dofs,
                                                              numbers::invalid_dof_index);
        {
          // first get the association of each dof with a subdomain and
          // determine the total number of subdomain ids used
          const std::vector<types::subdomain_id> subdomain_association
            = get_dof_subdomain_association (*this->dof_handler, n_dofs, n_procs);

          // then renumber the subdomains by first looking at those belonging
          // to subdomain 0, then those of subdomain 1, etc. note that the
          // algorithm is stable, i.e. if two dofs i,j have i<j and belong to
          // the same subdomain, then they will be in this order also after
          // reordering
          types::global_dof_index next_free_index = 0;
          for (types::subdomain_id subdomain=0; subdomain<n_procs; ++subdomain)
            for (types::global_dof_index i=0; i<n_dofs; ++i)
              if (subdomain_association[i] == subdomain)
                {
                  Assert (new_dof_indices[i] == numbers::invalid_dof_index,
                          ExcInternalError());
                  new_dof_indices[i] = next_free_index;
                  ++next_free_index;
                }

          // we should have numbered all dofs
          Assert (next_free_index == n_dofs, ExcInternalError());
          Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
                             numbers::invalid_dof_index)
                  == new_dof_indices.end(),
                  ExcInternalError());
        }
        // finally do the renumbering. we can use the sequential
        // version of the function because we do things on all
        // cells and all cells have their subdomain ids and DoFs
        // correctly set
        Implementation::renumber_dofs (new_dof_indices, IndexSet(0),
                                       *this->dof_handler, true);

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
          = get_dof_subdomain_association (*this->dof_handler, n_dofs, n_procs);

        for (unsigned int i=1; i<n_dofs; ++i)
          Assert (subdomain_association[i] >= subdomain_association[i-1],
                  ExcInternalError());

        std::vector<IndexSet> locally_owned_dofs_per_processor (n_procs,
                                                                IndexSet(n_dofs));
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
                     (subdomain_association[end_index] == subdomain_association[start_index]))
                ++end_index;

              // we've now identified a range of same indices. set that
              // range in the corresponding IndexSet
              if (end_index > start_index)
                {
                  const unsigned int subdomain_owner = subdomain_association[start_index];
                  locally_owned_dofs_per_processor[subdomain_owner]
                  .add_range (start_index, end_index);
                }

              // then move on to thinking about the next range
              start_index = end_index;
            }
        }

        // finally, restore current subdomain ids
        if (tr->with_artificial_cells())
          {
            typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator
            cell = this->dof_handler->get_triangulation().begin_active(),
            endc = this->dof_handler->get_triangulation().end();

            for (unsigned int index=0; cell != endc; ++cell, ++index)
              cell->set_subdomain_id(saved_subdomain_ids[index]);
          }

        // return a NumberCache object made up from the sets of locally
        // owned DoFs
        return NumberCache(locally_owned_dofs_per_processor,
                           this->dof_handler->get_triangulation().locally_owned_subdomain ());
      }



      template <int dim, int spacedim>
      std::vector<NumberCache>
      ParallelShared<dim,spacedim>::
      distribute_mg_dofs () const
      {
        // this is not currently implemented; the algorithm should work
        // as above, though: first call the sequential numbering
        // algorithm, then re-enumerate subdomain-wise
        Assert(false,ExcNotImplemented());

        return std::vector<NumberCache>();
      }



      template <int dim, int spacedim>
      NumberCache
      ParallelShared<dim,spacedim>::
      renumber_dofs (const std::vector<types::global_dof_index> &new_numbers) const
      {

#ifndef DEAL_II_WITH_MPI
        (void)new_numbers;
        Assert (false, ExcNotImplemented());
        return NumberCache();
#else
        // Similar to distribute_dofs() we need to have a special treatment in
        // case artificial cells are present.
        const parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>*> (&this->dof_handler->get_triangulation()));
        Assert(tr != nullptr, ExcInternalError());

        typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator
        cell = this->dof_handler->get_triangulation().begin_active(),
        endc = this->dof_handler->get_triangulation().end();
        std::vector<types::subdomain_id> current_subdomain_ids(tr->n_active_cells());
        const std::vector<types::subdomain_id> &true_subdomain_ids = tr->get_true_subdomain_ids_of_cells();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            {
              current_subdomain_ids[index] = cell->subdomain_id();
              cell->set_subdomain_id(true_subdomain_ids[index]);
            }

        std::vector<types::global_dof_index> global_gathered_numbers (this->dof_handler->n_dofs (), 0);
        // as we call DoFRenumbering::subdomain_wise (*dof_handler) from distribute_dofs(),
        // we need to support sequential-like input.
        // Distributed-like input from, for example, component_wise renumbering is also supported.
        if (new_numbers.size () == this->dof_handler->n_dofs ())
          {
            global_gathered_numbers = new_numbers;
          }
        else
          {
            Assert(new_numbers.size() == this->dof_handler->locally_owned_dofs().n_elements(),
                   ExcInternalError());
            const unsigned int n_cpu = Utilities::MPI::n_mpi_processes (tr->get_communicator ());
            std::vector<types::global_dof_index> gathered_new_numbers (this->dof_handler->n_dofs (), 0);
            Assert(Utilities::MPI::this_mpi_process (tr->get_communicator ()) ==
                   this->dof_handler->get_triangulation().locally_owned_subdomain (),
                   ExcInternalError())

            // gather new numbers among processors into one vector
            {
              std::vector<types::global_dof_index> new_numbers_copy (new_numbers);

              // store the number of elements that are to be received from each process
              std::vector<int> rcounts(n_cpu);

              types::global_dof_index shift = 0;
              // set rcounts based on new_numbers:
              int cur_count = new_numbers_copy.size ();
              int ierr = MPI_Allgather (&cur_count,  1, MPI_INT,
                                        &rcounts[0], 1, MPI_INT,
                                        tr->get_communicator ());
              AssertThrowMPI(ierr);

              // compute the displacements (relative to recvbuf)
              // at which to place the incoming data from process i
              std::vector<int> displacements(n_cpu);
              for (unsigned int i = 0; i < n_cpu; i++)
                {
                  displacements[i]  = shift;
                  shift            += rcounts[i];
                }
              Assert(((int)new_numbers_copy.size()) ==
                     rcounts[Utilities::MPI::this_mpi_process (tr->get_communicator ())],
                     ExcInternalError());
              ierr = MPI_Allgatherv (&new_numbers_copy[0],     new_numbers_copy.size (),
                                     DEAL_II_DOF_INDEX_MPI_TYPE,
                                     &gathered_new_numbers[0], &rcounts[0],
                                     &displacements[0],
                                     DEAL_II_DOF_INDEX_MPI_TYPE,
                                     tr->get_communicator ());
              AssertThrowMPI(ierr);
            }

            // put new numbers according to the current locally_owned_dofs_per_processor IndexSets
            types::global_dof_index shift = 0;
            // flag_1 and flag_2 are
            // used to control that there is a
            // one-to-one relation between old and new DoFs.
            std::vector<unsigned int> flag_1 (this->dof_handler->n_dofs (), 0);
            std::vector<unsigned int> flag_2 (this->dof_handler->n_dofs (), 0);
            for (unsigned int i = 0; i < n_cpu; i++)
              {
                const IndexSet iset =
                  this->dof_handler->locally_owned_dofs_per_processor()[i];
                for (types::global_dof_index ind = 0;
                     ind < iset.n_elements (); ind++)
                  {
                    const types::global_dof_index target = iset.nth_index_in_set (ind);
                    const types::global_dof_index value  = gathered_new_numbers[shift + ind];
                    Assert(target < this->dof_handler->n_dofs(), ExcInternalError());
                    Assert(value  < this->dof_handler->n_dofs(), ExcInternalError());
                    global_gathered_numbers[target] = value;
                    flag_1[target]++;
                    flag_2[value]++;
                  }
                shift += iset.n_elements ();
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
        Implementation::renumber_dofs (global_gathered_numbers, IndexSet(0),
                                       *this->dof_handler, true);

        const NumberCache number_cache (DoFTools::locally_owned_dofs_per_subdomain (*this->dof_handler),
                                        this->dof_handler->get_triangulation().locally_owned_subdomain ());

        // restore artificial cells
        cell = tr->begin_active();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            cell->set_subdomain_id(current_subdomain_ids[index]);

        return number_cache;
#endif
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
          std::vector<unsigned int>                                           tree_index;
          std::vector<typename dealii::internal::p4est::types<dim>::quadrant> quadrants;
          std::vector<dealii::types::global_dof_index>                        dof_numbers_and_indices;


          unsigned int bytes_for_buffer () const
          {
            return (sizeof(unsigned int)*2 +
                    tree_index.size() * sizeof(unsigned int) +
                    quadrants.size() * sizeof(typename dealii::internal::p4est
                                              ::types<dim>::quadrant) +
                    dof_numbers_and_indices.size() * sizeof(dealii::types::global_dof_index));
          }


          /**
           * Write the data of this object to a stream for the purpose of
           * serialization.
           */
          template <class Archive>
          void save (Archive &ar,
                     const unsigned int /*version*/) const
          {
            // we would like to directly serialize the 'quadrants' vector,
            // but the element type is internal to p4est and does not
            // know how to serialize itself. consequently, first copy it over
            // to an array of bytes, and then serialize that
            std::vector<char> quadrants_as_chars (sizeof(quadrants[0]) * quadrants.size());
            std::memcpy(&quadrants_as_chars[0],
                        &quadrants[0],
                        quadrants_as_chars.size());

            // now serialize everything
            ar &quadrants_as_chars
            &tree_index
            &dof_numbers_and_indices;
          }

          /**
           * Read the data of this object from a stream for the purpose of
           * serialization. Throw away the previous content.
           */
          template <class Archive>
          void load (Archive &ar,
                     const unsigned int /*version*/)
          {
            // undo the copying trick from the 'save' function
            std::vector<char> quadrants_as_chars;
            ar &quadrants_as_chars
            &tree_index
            &dof_numbers_and_indices;

            quadrants.resize (quadrants_as_chars.size() / sizeof(quadrants[0]));
            std::memcpy(&quadrants[0],
                        &quadrants_as_chars[0],
                        quadrants_as_chars.size());
          }

          BOOST_SERIALIZATION_SPLIT_MEMBER()


          /**
           * Pack the data that corresponds to this object into a buffer in
           * the form of a vector of chars and return it.
           */
          std::vector<char>
          pack_data () const
          {
            // set up a buffer and then use it as the target of a compressing
            // stream into which we serialize the current object
            std::vector<char> buffer;
            {
              boost::iostreams::filtering_ostream out;
              out.push(boost::iostreams::gzip_compressor
                       (boost::iostreams::gzip_params
                        (boost::iostreams::gzip::best_compression)));
              out.push(boost::iostreams::back_inserter(buffer));

              boost::archive::binary_oarchive archive(out);

              archive << *this;
              out.flush();
            }

            return buffer;
          }


          /**
           * Given a buffer in the form of an array of chars, unpack it and
           * restore the current object to the state that it was when
           * it was packed into said buffer by the pack_data() function.
           */
          void unpack_data (const std::vector<char> &buffer)
          {
            std::string decompressed_buffer;

            // first decompress the buffer
            {
              boost::iostreams::filtering_ostream decompressing_stream;
              decompressing_stream.push(boost::iostreams::gzip_decompressor());
              decompressing_stream.push(boost::iostreams::back_inserter(decompressed_buffer));

              for (const auto p : buffer)
                decompressing_stream << p;
            }

            // then restore the object from the buffer
            std::istringstream in(decompressed_buffer);
            boost::archive::binary_iarchive archive(in);

            archive >> *this;
          }
        };



        template <int dim, int spacedim>
        void
        fill_dofindices_recursively (const typename parallel::distributed::Triangulation<dim,spacedim> &tria,
                                     const unsigned int tree_index,
                                     const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
                                     const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
                                     const std::map<unsigned int, std::set<dealii::types::subdomain_id> > &vertices_with_ghost_neighbors,
                                     std::map<dealii::types::subdomain_id, CellDataTransferBuffer<dim>> &needs_to_get_cell)
        {
          // see if we have to recurse...
          if (dealii_cell->has_children())
            {
              typename dealii::internal::p4est::types<dim>::quadrant
              p4est_child[GeometryInfo<dim>::max_children_per_cell];
              internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);


              for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                fill_dofindices_recursively<dim,spacedim>(tria,
                                                          tree_index,
                                                          dealii_cell->child(c),
                                                          p4est_child[c],
                                                          vertices_with_ghost_neighbors,
                                                          needs_to_get_cell);
              return;
            }

          // we're at a leaf cell. see if the cell is flagged as
          // interesting. note that we have only flagged our own cells
          // before
          if (dealii_cell->user_flag_set() && !dealii_cell->is_ghost())
            {
              Assert (!dealii_cell->is_artificial(), ExcInternalError());

              // check each vertex if it is interesting and push
              // dofindices if yes
              std::set<dealii::types::subdomain_id> send_to;
              for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
                {
                  const std::map<unsigned int, std::set<dealii::types::subdomain_id> >::const_iterator
                  neighbor_subdomains_of_vertex
                    = vertices_with_ghost_neighbors.find (dealii_cell->vertex_index(v));

                  if (neighbor_subdomains_of_vertex ==
                      vertices_with_ghost_neighbors.end())
                    continue;

                  Assert(neighbor_subdomains_of_vertex->second.size()!=0,
                         ExcInternalError());

                  send_to.insert(neighbor_subdomains_of_vertex->second.begin(),
                                 neighbor_subdomains_of_vertex->second.end());
                }

              if (send_to.size() > 0)
                {
                  // this cell's dof_indices need to be sent to
                  // someone
                  std::vector<dealii::types::global_dof_index>
                  local_dof_indices (dealii_cell->get_fe().dofs_per_cell);
                  dealii_cell->get_dof_indices (local_dof_indices);

                  for (std::set<dealii::types::subdomain_id>::iterator it=send_to.begin();
                       it!=send_to.end(); ++it)
                    {
                      const dealii::types::subdomain_id subdomain = *it;

                      // get an iterator to what needs to be sent to
                      // that subdomain (if already exists), or create
                      // such an object
                      typename std::map<dealii::types::subdomain_id, CellDataTransferBuffer<dim>>::iterator
                          p
                          = needs_to_get_cell.insert (std::make_pair(subdomain,
                                                                     CellDataTransferBuffer<dim>()))
                            .first;

                      p->second.tree_index.push_back(tree_index);
                      p->second.quadrants.push_back(p4est_cell);

                      p->second.dof_numbers_and_indices.push_back(dealii_cell->get_fe().dofs_per_cell);
                      p->second.dof_numbers_and_indices.insert(p->second.dof_numbers_and_indices.end(),
                                                               local_dof_indices.begin(),
                                                               local_dof_indices.end());

                    }
                }
            }
        }

        template <int dim, int spacedim>
        void
        get_mg_dofindices_recursively (
          const parallel::distributed::Triangulation<dim,spacedim> &tria,
          const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
          const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant &quadrant,
          CellDataTransferBuffer<dim> &cell_data_transfer_buffer)
        {
          if (internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              // why would somebody request a cell that is not ours?
              Assert(dealii_cell->level_subdomain_id()==tria.locally_owned_subdomain(), ExcInternalError());


              std::vector<dealii::types::global_dof_index>
              local_dof_indices (dealii_cell->get_fe().dofs_per_cell);
              dealii_cell->get_mg_dof_indices (local_dof_indices);

              cell_data_transfer_buffer.dof_numbers_and_indices.push_back(dealii_cell->get_fe().dofs_per_cell);
              cell_data_transfer_buffer.dof_numbers_and_indices.insert(cell_data_transfer_buffer.dof_numbers_and_indices.end(),
                                                                       local_dof_indices.begin(),
                                                                       local_dof_indices.end());
              return; // we are done
            }

          if (! dealii_cell->has_children())
            return;

          if (! internal::p4est::quadrant_is_ancestor<dim> (p4est_cell, quadrant))
            return;

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
          internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            get_mg_dofindices_recursively<dim,spacedim> (tria, p4est_child[c],
                                                         dealii_cell->child(c),
                                                         quadrant, cell_data_transfer_buffer);
        }


        template <int dim, int spacedim>
        void
        find_marked_mg_ghost_cells_recursively(const typename parallel::distributed::Triangulation<dim,spacedim> &tria,
                                               const unsigned int tree_index,
                                               const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
                                               const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
                                               std::map<dealii::types::subdomain_id, CellDataTransferBuffer<dim>> &neighbor_cell_list)
        {
          // recurse...
          if (dealii_cell->has_children())
            {
              typename dealii::internal::p4est::types<dim>::quadrant
              p4est_child[GeometryInfo<dim>::max_children_per_cell];
              internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);


              for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                find_marked_mg_ghost_cells_recursively<dim,spacedim>(tria,
                                                                     tree_index,
                                                                     dealii_cell->child(c),
                                                                     p4est_child[c],
                                                                     neighbor_cell_list);
            }

          if (dealii_cell->user_flag_set() && dealii_cell->level_subdomain_id() != tria.locally_owned_subdomain())
            {
              neighbor_cell_list[dealii_cell->level_subdomain_id()].tree_index.push_back(tree_index);
              neighbor_cell_list[dealii_cell->level_subdomain_id()].quadrants.push_back(p4est_cell);
            }
        }


        template <int dim, int spacedim>
        void
        set_mg_dofindices_recursively (
          const parallel::distributed::Triangulation<dim,spacedim> &tria,
          const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
          const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant &quadrant,
          dealii::types::global_dof_index *dofs)
        {
          if (internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              Assert(dealii_cell->level_subdomain_id()!=dealii::numbers::artificial_subdomain_id, ExcInternalError());

              // update dof indices of cell
              std::vector<dealii::types::global_dof_index>
              dof_indices (dealii_cell->get_fe().dofs_per_cell);
              dealii_cell->get_mg_dof_indices(dof_indices);

              bool complete = true;
              for (unsigned int i=0; i<dof_indices.size(); ++i)
                if (dofs[i] != numbers::invalid_dof_index)
                  {
                    Assert((dof_indices[i] ==
                            (numbers::invalid_dof_index))
                           ||
                           (dof_indices[i]==dofs[i]),
                           ExcInternalError());
                    dof_indices[i]=dofs[i];
                  }
                else
                  complete=false;

              if (!complete)
                const_cast
                <typename DoFHandler<dim,spacedim>::level_cell_iterator &>
                (dealii_cell)->set_user_flag();
              else
                const_cast
                <typename DoFHandler<dim,spacedim>::level_cell_iterator &>
                (dealii_cell)->clear_user_flag();

              const_cast
              <typename DoFHandler<dim,spacedim>::level_cell_iterator &>
              (dealii_cell)->set_mg_dof_indices(dof_indices);
              return;
            }

          if (! dealii_cell->has_children())
            return;

          if (! internal::p4est::quadrant_is_ancestor<dim> (p4est_cell, quadrant))
            return;

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
          internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            set_mg_dofindices_recursively<dim,spacedim> (tria, p4est_child[c],
                                                         dealii_cell->child(c),
                                                         quadrant, dofs);

        }



        template <int dim, int spacedim>
        void
        communicate_mg_ghost_cells(const typename parallel::distributed::Triangulation<dim,spacedim> &tria,
                                   DoFHandler<dim,spacedim> &dof_handler,
                                   const std::vector<dealii::types::global_dof_index> &coarse_cell_to_p4est_tree_permutation,
                                   const std::vector<dealii::types::global_dof_index> &p4est_tree_to_coarse_cell_permutation)
        {
          // build list of cells to request for each neighbor
          std::set<dealii::types::subdomain_id> level_ghost_owners = tria.level_ghost_owners();
          typedef std::map<dealii::types::subdomain_id, CellDataTransferBuffer<dim>> cellmap_t;
          cellmap_t neighbor_cell_list;
          for (std::set<dealii::types::subdomain_id>::iterator it = level_ghost_owners.begin();
               it != level_ghost_owners.end();
               ++it)
            neighbor_cell_list.insert(std::make_pair(*it, CellDataTransferBuffer<dim>()));

          for (typename DoFHandler<dim,spacedim>::level_cell_iterator
               cell = dof_handler.begin(0);
               cell != dof_handler.end(0);
               ++cell)
            {
              typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
              internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

              find_marked_mg_ghost_cells_recursively<dim,spacedim>
              (tria,
               coarse_cell_to_p4est_tree_permutation[cell->index()],
               cell,
               p4est_coarse_cell,
               neighbor_cell_list);
            }
          Assert(level_ghost_owners.size() == neighbor_cell_list.size(), ExcInternalError());

          //* send our requests:
          std::vector<std::vector<char> > sendbuffers (level_ghost_owners.size());
          std::vector<MPI_Request> requests (level_ghost_owners.size());

          unsigned int idx=0;
          for (typename cellmap_t::iterator it = neighbor_cell_list.begin();
               it!=neighbor_cell_list.end();
               ++it, ++idx)
            {
              // pack all the data into the buffer for this recipient
              // and send it. keep data around till we can make sure
              // that the packet has been received
              sendbuffers[idx] = it->second.pack_data ();
              const int ierr = MPI_Isend(sendbuffers[idx].data(), sendbuffers[idx].size(),
                                         MPI_BYTE, it->first,
                                         1100101, tria.get_communicator(), &requests[idx]);
              AssertThrowMPI(ierr);
            }

          //* receive requests and reply
          std::vector<std::vector<char> > reply_buffers (level_ghost_owners.size());
          std::vector<MPI_Request> reply_requests (level_ghost_owners.size());

          for (unsigned int idx=0; idx<level_ghost_owners.size(); ++idx)
            {
              std::vector<char> receive;
              CellDataTransferBuffer<dim> cell_data_transfer_buffer;

              MPI_Status status;
              int len;
              int ierr = MPI_Probe(MPI_ANY_SOURCE, 1100101, tria.get_communicator(), &status);
              AssertThrowMPI(ierr);
              ierr = MPI_Get_count(&status, MPI_BYTE, &len);
              AssertThrowMPI(ierr);
              receive.resize(len);

              char *ptr = &receive[0];
              ierr = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                              tria.get_communicator(), &status);
              AssertThrowMPI(ierr);

              cell_data_transfer_buffer.unpack_data(receive);

              // store the dof indices for each cell
              for (unsigned int c=0; c<cell_data_transfer_buffer.tree_index.size(); ++c)
                {
                  typename DoFHandler<dim,spacedim>::level_cell_iterator
                  cell (&dof_handler.get_triangulation(),
                        0,
                        p4est_tree_to_coarse_cell_permutation[cell_data_transfer_buffer.tree_index[c]],
                        &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  get_mg_dofindices_recursively<dim,spacedim> (tria,
                                                               p4est_coarse_cell,
                                                               cell,
                                                               cell_data_transfer_buffer.quadrants[c],
                                                               cell_data_transfer_buffer);
                }

              // send reply
              reply_buffers[idx] = cell_data_transfer_buffer.pack_data();
              ierr = MPI_Isend(&(reply_buffers[idx])[0], reply_buffers[idx].size(),
                               MPI_BYTE, status.MPI_SOURCE,
                               1100102, tria.get_communicator(), &reply_requests[idx]);
              AssertThrowMPI(ierr);
            }

          //* finally receive the replies
          for (unsigned int idx=0; idx<level_ghost_owners.size(); ++idx)
            {
              std::vector<char> receive;
              CellDataTransferBuffer<dim> cell_data_transfer_buffer;

              MPI_Status status;
              int len;
              int ierr = MPI_Probe(MPI_ANY_SOURCE, 1100102, tria.get_communicator(), &status);
              AssertThrowMPI(ierr);
              ierr = MPI_Get_count(&status, MPI_BYTE, &len);
              AssertThrowMPI(ierr);
              receive.resize(len);

              char *ptr = &receive[0];
              ierr = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                              tria.get_communicator(), &status);
              AssertThrowMPI(ierr);

              cell_data_transfer_buffer.unpack_data(receive);
              if (cell_data_transfer_buffer.tree_index.size()==0)
                continue;

              // set the dof indices for each cell
              dealii::types::global_dof_index *dofs = cell_data_transfer_buffer.dof_numbers_and_indices.data();
              for (unsigned int c=0; c<cell_data_transfer_buffer.tree_index.size(); ++c, dofs+=1+dofs[0])
                {
                  typename DoFHandler<dim,spacedim>::level_cell_iterator
                  cell (&tria,
                        0,
                        p4est_tree_to_coarse_cell_permutation[cell_data_transfer_buffer.tree_index[c]],
                        &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  Assert(cell->get_fe().dofs_per_cell==dofs[0], ExcInternalError());

                  set_mg_dofindices_recursively<dim,spacedim> (tria,
                                                               p4est_coarse_cell,
                                                               cell,
                                                               cell_data_transfer_buffer.quadrants[c],
                                                               dofs+1);
                }
            }

          // complete all sends, so that we can safely destroy the
          // buffers.
          if (requests.size() > 0)
            {
              const int ierr = MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }
          if (reply_requests.size() > 0)
            {
              const int ierr = MPI_Waitall(reply_requests.size(), &reply_requests[0], MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }
        }



        template <int spacedim>
        void
        communicate_mg_ghost_cells(const typename parallel::distributed::Triangulation<1,spacedim> &,
                                   DoFHandler<1,spacedim> &,
                                   const std::vector<dealii::types::global_dof_index> &,
                                   const std::vector<dealii::types::global_dof_index> &)
        {
          Assert (false, ExcNotImplemented());
        }




        template <int dim, int spacedim>
        void
        set_dofindices_recursively (
          const parallel::distributed::Triangulation<dim,spacedim> &tria,
          const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
          const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant &quadrant,
          dealii::types::global_dof_index *dofs)
        {
          if (internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              Assert(!dealii_cell->has_children(), ExcInternalError());
              Assert(dealii_cell->is_ghost(), ExcInternalError());

              // update dof indices of cell
              std::vector<dealii::types::global_dof_index>
              dof_indices (dealii_cell->get_fe().dofs_per_cell);
              dealii_cell->update_cell_dof_indices_cache();
              dealii_cell->get_dof_indices(dof_indices);

              bool complete = true;
              for (unsigned int i=0; i<dof_indices.size(); ++i)
                if (dofs[i] != numbers::invalid_dof_index)
                  {
                    Assert((dof_indices[i] ==
                            (numbers::invalid_dof_index))
                           ||
                           (dof_indices[i]==dofs[i]),
                           ExcInternalError());
                    dof_indices[i]=dofs[i];
                  }
                else
                  complete=false;

              if (!complete)
                const_cast
                <typename DoFHandler<dim,spacedim>::level_cell_iterator &>
                (dealii_cell)->set_user_flag();
              else
                const_cast
                <typename DoFHandler<dim,spacedim>::level_cell_iterator &>
                (dealii_cell)->clear_user_flag();

              const_cast
              <typename DoFHandler<dim,spacedim>::level_cell_iterator &>
              (dealii_cell)->set_dof_indices(dof_indices);

              return;
            }

          if (! dealii_cell->has_children())
            return;

          if (! internal::p4est::quadrant_is_ancestor<dim> (p4est_cell, quadrant))
            return;

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
          internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            set_dofindices_recursively<dim,spacedim> (tria, p4est_child[c],
                                                      dealii_cell->child(c),
                                                      quadrant, dofs);
        }



        template <int spacedim>
        void
        communicate_dof_indices_on_marked_cells
        (const DoFHandler<1,spacedim> &,
         const std::map<unsigned int, std::set<dealii::types::subdomain_id> > &,
         const std::vector<dealii::types::global_dof_index> &,
         const std::vector<dealii::types::global_dof_index> &)
        {
          Assert (false, ExcNotImplemented());
        }



        template <int dim, int spacedim>
        void
        communicate_dof_indices_on_marked_cells
        (const DoFHandler<dim,spacedim> &dof_handler,
         const std::map<unsigned int, std::set<dealii::types::subdomain_id> > &vertices_with_ghost_neighbors,
         const std::vector<dealii::types::global_dof_index> &coarse_cell_to_p4est_tree_permutation,
         const std::vector<dealii::types::global_dof_index> &p4est_tree_to_coarse_cell_permutation)
        {
#ifndef DEAL_II_WITH_P4EST
          (void)vertices_with_ghost_neighbors;
          Assert (false, ExcNotImplemented());
#else

          const parallel::distributed::Triangulation< dim, spacedim > *tr
            = (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
               (&dof_handler.get_triangulation()));
          Assert (tr != nullptr, ExcInternalError());

          // now collect cells and their dof_indices for the
          // interested neighbors
          typedef
          std::map<dealii::types::subdomain_id, CellDataTransferBuffer<dim>>
              cellmap_t;
          cellmap_t needs_to_get_cells;

          for (typename DoFHandler<dim,spacedim>::level_cell_iterator
               cell = dof_handler.begin(0);
               cell != dof_handler.end(0);
               ++cell)
            {
              typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
              internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

              fill_dofindices_recursively<dim,spacedim>
              (*tr,
               coarse_cell_to_p4est_tree_permutation[cell->index()],
               cell,
               p4est_coarse_cell,
               vertices_with_ghost_neighbors,
               needs_to_get_cells);
            }


          // sending
          std::vector<std::vector<char> > sendbuffers (needs_to_get_cells.size());
          std::vector<MPI_Request> requests (needs_to_get_cells.size());

          unsigned int                              idx    = 0;
          std::vector<std::vector<char> >::iterator buffer = sendbuffers.begin();
          for (typename cellmap_t::iterator it=needs_to_get_cells.begin();
               it!=needs_to_get_cells.end();
               ++it, ++buffer, ++idx)
            {
              const unsigned int num_cells = it->second.tree_index.size();
              (void)num_cells;

              Assert(num_cells==it->second.quadrants.size(), ExcInternalError());
              Assert(num_cells>0, ExcInternalError());

              // pack all the data into the buffer for this recipient
              // and send it. keep data around till we can make sure
              // that the packet has been received
              *buffer = it->second.pack_data ();
              const int ierr = MPI_Isend(&(*buffer)[0], buffer->size(),
                                         MPI_BYTE, it->first,
                                         123, tr->get_communicator(), &requests[idx]);
              AssertThrowMPI(ierr);
            }


          // mark all of our own cells that miss some dof_data and collect
          // the neighbors that are going to send stuff to us
          std::set<dealii::types::subdomain_id> senders;
          {
            std::vector<dealii::types::global_dof_index> local_dof_indices;
            typename DoFHandler<dim,spacedim>::active_cell_iterator
            cell, endc = dof_handler.end();

            for (cell = dof_handler.begin_active(); cell != endc; ++cell)
              if (!cell->is_artificial())
                {
                  if (cell->is_ghost())
                    {
                      if (cell->user_flag_set())
                        senders.insert(cell->subdomain_id());
                    }
                  else
                    {
                      local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                      cell->get_dof_indices (local_dof_indices);
                      if (local_dof_indices.end() !=
                          std::find (local_dof_indices.begin(),
                                     local_dof_indices.end(),
                                     numbers::invalid_dof_index))
                        cell->set_user_flag();
                      else
                        cell->clear_user_flag();
                    }
                }
          }


          //* 5. receive ghostcelldata
          std::vector<char> receive;
          for (unsigned int i=0; i<senders.size(); ++i)
            {
              MPI_Status status;
              int len;
              int ierr = MPI_Probe(MPI_ANY_SOURCE, 123, tr->get_communicator(), &status);
              AssertThrowMPI(ierr);
              ierr = MPI_Get_count(&status, MPI_BYTE, &len);
              AssertThrowMPI(ierr);
              receive.resize(len);

              char *ptr = &receive[0];
              ierr = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                              tr->get_communicator(), &status);
              AssertThrowMPI(ierr);

              CellDataTransferBuffer<dim> cell_data_transfer_buffer;
              cell_data_transfer_buffer.unpack_data(receive);
              unsigned int cells = cell_data_transfer_buffer.tree_index.size();
              dealii::types::global_dof_index *dofs = cell_data_transfer_buffer.dof_numbers_and_indices.data();

              // the dofs pointer contains for each cell the number of
              // dofs on that cell (dofs[0]) followed by the dof
              // indices itself.
              for (unsigned int c=0; c<cells; ++c, dofs+=1+dofs[0])
                {
                  typename DoFHandler<dim,spacedim>::level_cell_iterator
                  cell (&dof_handler.get_triangulation(),
                        0,
                        p4est_tree_to_coarse_cell_permutation[cell_data_transfer_buffer.tree_index[c]],
                        &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  Assert(cell->get_fe().dofs_per_cell==dofs[0], ExcInternalError());

                  set_dofindices_recursively<dim,spacedim> (*tr,
                                                            p4est_coarse_cell,
                                                            cell,
                                                            cell_data_transfer_buffer.quadrants[c],
                                                            (dofs+1));
                }
            }

          // complete all sends, so that we can safely destroy the
          // buffers.
          if (requests.size() > 0)
            {
              const int ierr = MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }


#ifdef DEBUG
          {
            // check all msgs got sent and received
            unsigned int sum_send=0;
            unsigned int sum_recv=0;
            unsigned int sent=needs_to_get_cells.size();
            unsigned int recv=senders.size();

            int ierr = MPI_Allreduce(&sent, &sum_send, 1, MPI_UNSIGNED, MPI_SUM, tr->get_communicator());
            AssertThrowMPI(ierr);
            ierr = MPI_Allreduce(&recv, &sum_recv, 1, MPI_UNSIGNED, MPI_SUM, tr->get_communicator());
            AssertThrowMPI(ierr);
            Assert(sum_send==sum_recv, ExcInternalError());
          }
#endif

          // update dof index caches on all cells
          update_all_active_cell_dof_indices_caches (dof_handler);

          // have a barrier so that sends between two calls to this
          // function are not mixed up.
          //
          // this is necessary because above we just see if there are
          // messages and then receive them, without discriminating
          // where they come from and whether they were sent in phase
          // 1 or 2. the need for a global communication step like
          // this barrier could be avoided by receiving messages
          // specifically from those processors from which we expect
          // messages, and by using different tags for phase 1 and 2
          const int ierr = MPI_Barrier(tr->get_communicator());
          AssertThrowMPI(ierr);
#endif
        }







      }

#endif // DEAL_II_WITH_P4EST



      template <int dim, int spacedim>
      ParallelDistributed<dim,spacedim>::
      ParallelDistributed (dealii::DoFHandler<dim,spacedim> &dof_handler)
        :
        dof_handler (&dof_handler)
      {}






      template <int dim, int spacedim>
      NumberCache
      ParallelDistributed<dim, spacedim>::
      distribute_dofs () const
      {
        NumberCache number_cache;

#ifndef DEAL_II_WITH_P4EST
        Assert (false, ExcNotImplemented());
#else

        parallel::distributed::Triangulation< dim, spacedim > *tr
          = (dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (const_cast<dealii::Triangulation< dim, spacedim >*>
              (&dof_handler->get_triangulation())));
        Assert (tr != nullptr, ExcInternalError());

        const unsigned int
        n_cpus = Utilities::MPI::n_mpi_processes (tr->get_communicator());

        //* 1. distribute on own subdomain
        const dealii::types::global_dof_index n_initial_local_dofs =
          Implementation::distribute_dofs (tr->locally_owned_subdomain(),
                                           *dof_handler);

        //* 2. iterate over ghostcells and kill dofs that are not
        // owned by us
        std::vector<dealii::types::global_dof_index> renumbering(n_initial_local_dofs);
        for (unsigned int i=0; i<renumbering.size(); ++i)
          renumbering[i] = i;

        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          typename DoFHandler<dim,spacedim>::active_cell_iterator
          cell = dof_handler->begin_active(),
          endc = dof_handler->end();

          for (; cell != endc; ++cell)
            if (cell->is_ghost() &&
                (cell->subdomain_id() < tr->locally_owned_subdomain()))
              {
                // we found a neighboring ghost cell whose subdomain
                // is "stronger" than our own subdomain

                // delete all dofs that live there and that we have
                // previously assigned a number to (i.e. the ones on
                // the interface)
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                  if (local_dof_indices[i] != numbers::invalid_dof_index)
                    renumbering[local_dof_indices[i]]
                      = numbers::invalid_dof_index;
              }
        }


        // make indices consecutive
        number_cache.n_locally_owned_dofs = 0;
        for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
             it!=renumbering.end(); ++it)
          if (*it != numbers::invalid_dof_index)
            *it = number_cache.n_locally_owned_dofs++;

        //* 3. communicate local dofcount and shift ids to make them
        // unique
        number_cache.n_locally_owned_dofs_per_processor.resize(n_cpus);

        const int ierr = MPI_Allgather ( &number_cache.n_locally_owned_dofs,
                                         1, DEAL_II_DOF_INDEX_MPI_TYPE,
                                         &number_cache.n_locally_owned_dofs_per_processor[0],
                                         1, DEAL_II_DOF_INDEX_MPI_TYPE,
                                         tr->get_communicator());
        AssertThrowMPI(ierr);

        const dealii::types::global_dof_index
        shift = std::accumulate (number_cache
                                 .n_locally_owned_dofs_per_processor.begin(),
                                 number_cache
                                 .n_locally_owned_dofs_per_processor.begin()
                                 + tr->locally_owned_subdomain(),
                                 static_cast<dealii::types::global_dof_index>(0));
        for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
             it!=renumbering.end(); ++it)
          if (*it != numbers::invalid_dof_index)
            (*it) += shift;

        // now re-enumerate all dofs to this shifted and condensed
        // numbering form.  we renumber some dofs as invalid, so
        // choose the nocheck-version.
        Implementation::renumber_dofs (renumbering, IndexSet(0),
                                       *dof_handler, false);

        // now a little bit of housekeeping
        number_cache.n_global_dofs
          = std::accumulate (number_cache
                             .n_locally_owned_dofs_per_processor.begin(),
                             number_cache
                             .n_locally_owned_dofs_per_processor.end(),
                             static_cast<dealii::types::global_dof_index>(0));

        number_cache.locally_owned_dofs = IndexSet(number_cache.n_global_dofs);
        number_cache.locally_owned_dofs
        .add_range(shift,
                   shift+number_cache.n_locally_owned_dofs);
        number_cache.locally_owned_dofs.compress();

        // fill global_dof_indexsets
        number_cache.locally_owned_dofs_per_processor.resize(n_cpus);
        {
          dealii::types::global_dof_index lshift = 0;
          for (unsigned int i=0; i<n_cpus; ++i)
            {
              number_cache.locally_owned_dofs_per_processor[i]
                = IndexSet(number_cache.n_global_dofs);
              number_cache.locally_owned_dofs_per_processor[i]
              .add_range(lshift,
                         lshift +
                         number_cache.n_locally_owned_dofs_per_processor[i]);
              lshift += number_cache.n_locally_owned_dofs_per_processor[i];
            }
        }
        Assert(number_cache.locally_owned_dofs_per_processor
               [tr->locally_owned_subdomain()].n_elements()
               ==
               number_cache.n_locally_owned_dofs,
               ExcInternalError());
        Assert(!number_cache.locally_owned_dofs_per_processor
               [tr->locally_owned_subdomain()].n_elements()
               ||
               number_cache.locally_owned_dofs_per_processor
               [tr->locally_owned_subdomain()].nth_index_in_set(0)
               == shift,
               ExcInternalError());

        //* 4. send dofids of cells that are ghostcells on other
        // machines

        std::vector<bool> user_flags;
        tr->save_user_flags(user_flags);
        tr->clear_user_flags ();

        // mark all own cells for transfer
        for (typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler->begin_active();
             cell != dof_handler->end(); ++cell)
          if (!cell->is_artificial())
            cell->set_user_flag();

        // add each ghostcell's subdomain to the vertex and keep track
        // of interesting neighbors
        std::map<unsigned int, std::set<dealii::types::subdomain_id> >
        vertices_with_ghost_neighbors;

        tr->fill_vertices_with_ghost_neighbors (vertices_with_ghost_neighbors);


        // Send and receive cells. After this, only the local cells
        // are marked, that received new data. This has to be
        // communicated in a second communication step.
        communicate_dof_indices_on_marked_cells (*dof_handler,
                                                 vertices_with_ghost_neighbors,
                                                 tr->coarse_cell_to_p4est_tree_permutation,
                                                 tr->p4est_tree_to_coarse_cell_permutation);

        communicate_dof_indices_on_marked_cells (*dof_handler,
                                                 vertices_with_ghost_neighbors,
                                                 tr->coarse_cell_to_p4est_tree_permutation,
                                                 tr->p4est_tree_to_coarse_cell_permutation);

        tr->load_user_flags(user_flags);

#ifdef DEBUG
        // check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          for (typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler->begin_active();
               cell != dof_handler->end(); ++cell)
            if (!cell->is_artificial())
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                if (local_dof_indices.end() !=
                    std::find (local_dof_indices.begin(),
                               local_dof_indices.end(),
                               numbers::invalid_dof_index))
                  {
                    if (cell->is_ghost())
                      {
                        Assert(false, ExcMessage ("Not a ghost cell"));
                      }
                    else
                      {
                        Assert(false, ExcMessage ("Not one of our own cells"));
                      }
                  }
              }
        }
#endif // DEBUG
#endif // DEAL_II_WITH_P4EST

        return number_cache;
      }



      template <int dim, int spacedim>
      std::vector<NumberCache>
      ParallelDistributed<dim, spacedim>::
      distribute_mg_dofs () const
      {
#ifndef DEAL_II_WITH_P4EST
        Assert (false, ExcNotImplemented());
        return std::vector<NumberCache>();
#else

        parallel::distributed::Triangulation< dim, spacedim > *tr
          = (dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (const_cast<dealii::Triangulation< dim, spacedim >*>
              (&dof_handler->get_triangulation())));
        Assert (tr != nullptr, ExcInternalError());

        AssertThrow(
          (tr->settings & parallel::distributed::Triangulation< dim, spacedim >::construct_multigrid_hierarchy),
          ExcMessage("Multigrid DoFs can only be distributed on a parallel "
                     "Triangulation if the flag construct_multigrid_hierarchy "
                     "is set in the constructor."));


        const unsigned int
        n_cpus = Utilities::MPI::n_mpi_processes (tr->get_communicator());

        // loop over all levels that exist globally (across all
        // processors), even if the current processor does not in fact
        // have any cells on that level or if the local part of the
        // Triangulation has fewer levels. we need to do this because
        // we need to communicate across all processors on all levels
        const unsigned int n_levels = tr->n_global_levels();
        std::vector<NumberCache> number_caches;
        number_caches.reserve(n_levels);
        for (unsigned int level = 0; level < n_levels; ++level)
          {
            NumberCache level_number_cache;

            //* 1. distribute on own subdomain
            const unsigned int n_initial_local_dofs =
              Implementation::distribute_dofs_on_level(tr->locally_owned_subdomain(),
                                                       *dof_handler,
                                                       level);

            //* 2. iterate over ghostcells and kill dofs that are not
            // owned by us
            std::vector<dealii::types::global_dof_index> renumbering(n_initial_local_dofs);
            for (dealii::types::global_dof_index i=0; i<renumbering.size(); ++i)
              renumbering[i] = i;

            if (level < tr->n_levels())
              {
                std::vector<dealii::types::global_dof_index> local_dof_indices;

                typename DoFHandler<dim,spacedim>::level_cell_iterator
                cell = dof_handler->begin(level),
                endc = dof_handler->end(level);

                for (; cell != endc; ++cell)
                  if (cell->level_subdomain_id()!=numbers::artificial_subdomain_id &&
                      (cell->level_subdomain_id() < tr->locally_owned_subdomain()))
                    {
                      // we found a neighboring ghost cell whose
                      // subdomain is "stronger" than our own
                      // subdomain

                      // delete all dofs that live there and that we
                      // have previously assigned a number to
                      // (i.e. the ones on the interface)
                      local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                      cell->get_mg_dof_indices (local_dof_indices);
                      for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                        if (local_dof_indices[i] != numbers::invalid_dof_index)
                          renumbering[local_dof_indices[i]]
                            = numbers::invalid_dof_index;
                    }
              }

// TODO: make this code simpler with the new constructors of NumberCache
            // make indices consecutive
            level_number_cache.n_locally_owned_dofs = 0;
            for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
                 it!=renumbering.end(); ++it)
              if (*it != numbers::invalid_dof_index)
                *it = level_number_cache.n_locally_owned_dofs++;

            //* 3. communicate local dofcount and shift ids to make
            // them unique
            level_number_cache.n_locally_owned_dofs_per_processor.resize(n_cpus);

            int ierr = MPI_Allgather ( &level_number_cache.n_locally_owned_dofs,
                                       1, DEAL_II_DOF_INDEX_MPI_TYPE,
                                       &level_number_cache.n_locally_owned_dofs_per_processor[0],
                                       1, DEAL_II_DOF_INDEX_MPI_TYPE,
                                       tr->get_communicator());
            AssertThrowMPI(ierr);

            const dealii::types::global_dof_index
            shift = std::accumulate (level_number_cache
                                     .n_locally_owned_dofs_per_processor.begin(),
                                     level_number_cache
                                     .n_locally_owned_dofs_per_processor.begin()
                                     + tr->locally_owned_subdomain(),
                                     static_cast<dealii::types::global_dof_index>(0));
            for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
                 it!=renumbering.end(); ++it)
              if (*it != numbers::invalid_dof_index)
                (*it) += shift;

            // now re-enumerate all dofs to this shifted and condensed
            // numbering form.  we renumber some dofs as invalid, so
            // choose the nocheck-version of the function
            //
            // of course there is nothing for us to renumber if the
            // level we are currently dealing with doesn't even exist
            // within the current triangulation, so skip renumbering
            // in that case
            if (level < tr->n_levels())
              Implementation::renumber_mg_dofs (renumbering, IndexSet(0),
                                                *dof_handler, level,
                                                false);

            // now a little bit of housekeeping
            level_number_cache.n_global_dofs
              = std::accumulate (level_number_cache
                                 .n_locally_owned_dofs_per_processor.begin(),
                                 level_number_cache
                                 .n_locally_owned_dofs_per_processor.end(),
                                 static_cast<dealii::types::global_dof_index>(0));

            level_number_cache.locally_owned_dofs = IndexSet(level_number_cache.n_global_dofs);
            level_number_cache.locally_owned_dofs
            .add_range(shift,
                       shift+level_number_cache.n_locally_owned_dofs);
            level_number_cache.locally_owned_dofs.compress();

            // fill global_dof_indexsets
            level_number_cache.locally_owned_dofs_per_processor.resize(n_cpus);
            {
              dealii::types::global_dof_index lshift = 0;
              for (unsigned int i=0; i<n_cpus; ++i)
                {
                  level_number_cache.locally_owned_dofs_per_processor[i]
                    = IndexSet(level_number_cache.n_global_dofs);
                  level_number_cache.locally_owned_dofs_per_processor[i]
                  .add_range(lshift,
                             lshift +
                             level_number_cache.n_locally_owned_dofs_per_processor[i]);
                  lshift += level_number_cache.n_locally_owned_dofs_per_processor[i];
                }
            }
            Assert(level_number_cache.locally_owned_dofs_per_processor
                   [tr->locally_owned_subdomain()].n_elements()
                   ==
                   level_number_cache.n_locally_owned_dofs,
                   ExcInternalError());
            Assert(!level_number_cache.locally_owned_dofs_per_processor
                   [tr->locally_owned_subdomain()].n_elements()
                   ||
                   level_number_cache.locally_owned_dofs_per_processor
                   [tr->locally_owned_subdomain()].nth_index_in_set(0)
                   == shift,
                   ExcInternalError());

            number_caches.emplace_back (level_number_cache);
          }


        //* communicate ghost DoFs
        // We mark all ghost cells by setting the user_flag and then request
        // these cells from the corresponding owners. As this information
        // can be incomplete,
        {
          std::vector<bool> user_flags;
          tr->save_user_flags(user_flags);
          tr->clear_user_flags ();

          // mark all ghost cells for transfer
          {
            typename DoFHandler<dim,spacedim>::level_cell_iterator
            cell, endc = dof_handler->end();
            for (cell = dof_handler->begin(); cell != endc; ++cell)
              if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id
                  && !cell->is_locally_owned_on_level())
                cell->set_user_flag();
          }

          // Phase 1. Request all marked cells from corresponding owners. If we
          // managed to get every DoF, remove the user_flag, otherwise we
          // will request them again in the step below.
          communicate_mg_ghost_cells(*tr,
                                     *dof_handler,
                                     tr->coarse_cell_to_p4est_tree_permutation,
                                     tr->p4est_tree_to_coarse_cell_permutation);

          // This barrier is crucial so that messages between phases
          // 1&2 don't mix.
          const int ierr = MPI_Barrier(tr->get_communicator());
          AssertThrowMPI(ierr);

          // Phase 2, only request the cells that were not completed
          // in Phase 1.
          communicate_mg_ghost_cells(*tr,
                                     *dof_handler,
                                     tr->coarse_cell_to_p4est_tree_permutation,
                                     tr->p4est_tree_to_coarse_cell_permutation);

#ifdef DEBUG
          // make sure we have removed all flags:
          {
            typename DoFHandler<dim,spacedim>::level_cell_iterator
            cell, endc = dof_handler->end();
            for (cell = dof_handler->begin(); cell != endc; ++cell)
              if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id
                  && !cell->is_locally_owned_on_level())
                Assert(cell->user_flag_set()==false, ExcInternalError());
          }
#endif

          tr->load_user_flags(user_flags);
        }




#ifdef DEBUG
        // check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;
          typename DoFHandler<dim,spacedim>::level_cell_iterator
          cell, endc = dof_handler->end();

          for (cell = dof_handler->begin(); cell != endc; ++cell)
            if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id)
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_mg_dof_indices (local_dof_indices);
                if (local_dof_indices.end() !=
                    std::find (local_dof_indices.begin(),
                               local_dof_indices.end(),
                               numbers::invalid_dof_index))
                  {
                    Assert(false, ExcMessage ("not all DoFs got distributed!"));
                  }
              }
        }
#endif // DEBUG

        return number_caches;

#endif // DEAL_II_WITH_P4EST
      }


      template <int dim, int spacedim>
      NumberCache
      ParallelDistributed<dim, spacedim>::
      renumber_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers) const
      {
        (void)new_numbers;

        Assert (new_numbers.size() == dof_handler->locally_owned_dofs().n_elements(),
                ExcInternalError());

        NumberCache number_cache;

#ifndef DEAL_II_WITH_P4EST
        Assert (false, ExcNotImplemented());
#else


        // calculate new IndexSet. First try to find out if the new indices
        // are contiguous blocks. This avoids inserting each index
        // individually into the IndexSet, which is slow.  If we own no DoFs,
        // we still need to go through this function, but we can skip this
        // calculation.

        number_cache.locally_owned_dofs = IndexSet (dof_handler->n_dofs());
        if (dof_handler->locally_owned_dofs().n_elements()>0)
          {
            std::vector<dealii::types::global_dof_index> new_numbers_sorted (new_numbers);
            std::sort(new_numbers_sorted.begin(), new_numbers_sorted.end());
            std::vector<dealii::types::global_dof_index>::const_iterator it = new_numbers_sorted.begin();
            const unsigned int n_blocks = dof_handler->get_fe().n_blocks();
            std::vector<std::pair<dealii::types::global_dof_index,unsigned int> > block_indices(n_blocks);
            block_indices[0].first = *it++;
            block_indices[0].second = 1;
            unsigned int current_block = 0, n_filled_blocks = 1;
            for ( ; it != new_numbers_sorted.end(); ++it)
              {
                bool done = false;

                // search from the current block onwards whether the next
                // index is shifted by one from the previous one.
                for (unsigned int i=0; i<n_filled_blocks; ++i)
                  if (*it == block_indices[current_block].first
                      +block_indices[current_block].second)
                    {
                      block_indices[current_block].second++;
                      done = true;
                      break;
                    }
                  else
                    {
                      if (current_block == n_filled_blocks-1)
                        current_block = 0;
                      else
                        ++current_block;
                    }

                // could not find any contiguous range: need to add a new
                // block if possible. Abort otherwise, which will add all
                // elements individually to the IndexSet.
                if (done == false)
                  {
                    if (n_filled_blocks < n_blocks)
                      {
                        block_indices[n_filled_blocks].first = *it;
                        block_indices[n_filled_blocks].second = 1;
                        current_block = n_filled_blocks;
                        ++n_filled_blocks;
                      }
                    else
                      break;
                  }
              }

            // check whether all indices could be assigned to blocks. If yes,
            // we can add the block ranges to the IndexSet, otherwise we need
            // to go through the indices once again and add each element
            // individually
            unsigned int sum = 0;
            for (unsigned int i=0; i<n_filled_blocks; ++i)
              sum += block_indices[i].second;
            if (sum == new_numbers.size())
              for (unsigned int i=0; i<n_filled_blocks; ++i)
                number_cache.locally_owned_dofs.add_range (block_indices[i].first,
                                                           block_indices[i].first+
                                                           block_indices[i].second);
            else
              number_cache.locally_owned_dofs.add_indices(new_numbers_sorted.begin(),
                                                          new_numbers_sorted.end());
          }


        number_cache.locally_owned_dofs.compress();
        Assert (number_cache.locally_owned_dofs.n_elements() == new_numbers.size(),
                ExcInternalError());
        // also check with the number of locally owned degrees of freedom that
        // the DoFHandler object still stores
        Assert (number_cache.locally_owned_dofs.n_elements() ==
                dof_handler->n_locally_owned_dofs(),
                ExcInternalError());

        // then also set this number in our own copy
        number_cache.n_locally_owned_dofs = dof_handler->n_locally_owned_dofs();

        // mark not locally active DoFs as invalid
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          typename DoFHandler<dim,spacedim>::active_cell_iterator
          cell = dof_handler->begin_active(),
          endc = dof_handler->end();

          for (; cell != endc; ++cell)
            if (!cell->is_artificial())
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                  {
                    if (local_dof_indices[i] == numbers::invalid_dof_index)
                      continue;

                    if (!dof_handler->locally_owned_dofs().is_element(local_dof_indices[i]))
                      {
                        //this DoF is not owned by us, so set it to invalid.
                        local_dof_indices[i]
                          = numbers::invalid_dof_index;
                      }
                  }

                cell->set_dof_indices (local_dof_indices);
              }
        }


        // renumber. Skip when there is nothing to do because we own no DoF.
        if (dof_handler->locally_owned_dofs().n_elements() > 0)
          Implementation::renumber_dofs (new_numbers,
                                         dof_handler->locally_owned_dofs(),
                                         *dof_handler,
                                         false);

        // communication
        {
          parallel::distributed::Triangulation< dim, spacedim > *tr
            = (dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
               (const_cast<dealii::Triangulation< dim, spacedim >*>
                (&dof_handler->get_triangulation())));
          Assert (tr != nullptr, ExcInternalError());

          std::vector<bool> user_flags;
          tr->save_user_flags(user_flags);
          tr->clear_user_flags ();

          // mark all own cells for transfer
          typename DoFHandler<dim,spacedim>::active_cell_iterator
          cell, endc = dof_handler->end();
          for (cell = dof_handler->begin_active(); cell != endc; ++cell)
            if (!cell->is_artificial())
              cell->set_user_flag();

          // add each ghostcells' subdomain to the vertex and keep track of
          // interesting neighbors
          std::map<unsigned int, std::set<dealii::types::subdomain_id> >
          vertices_with_ghost_neighbors;

          tr->fill_vertices_with_ghost_neighbors (vertices_with_ghost_neighbors);

          // Send and receive cells. After this, only the local cells are
          // marked, that received new data. This has to be communicated in a
          // second communication step.
          communicate_dof_indices_on_marked_cells (*dof_handler,
                                                   vertices_with_ghost_neighbors,
                                                   tr->coarse_cell_to_p4est_tree_permutation,
                                                   tr->p4est_tree_to_coarse_cell_permutation);

          communicate_dof_indices_on_marked_cells (*dof_handler,
                                                   vertices_with_ghost_neighbors,
                                                   tr->coarse_cell_to_p4est_tree_permutation,
                                                   tr->p4est_tree_to_coarse_cell_permutation);


          // * Create global_dof_indexsets by transferring our own owned_dofs
          // to every other machine.
          const unsigned int n_cpus = Utilities::MPI::n_mpi_processes (tr->get_communicator());

          // Serialize our IndexSet and determine size.
          std::ostringstream oss;
          number_cache.locally_owned_dofs.block_write(oss);
          std::string oss_str=oss.str();
          std::vector<char> my_data(oss_str.begin(), oss_str.end());
          unsigned int my_size = oss_str.size();

          // determine maximum size of IndexSet
          const unsigned int max_size
            = Utilities::MPI::max (my_size, tr->get_communicator());

          // as we are reading past the end, we need to increase the size of
          // the local buffer. This is filled with zeros.
          my_data.resize(max_size);

          std::vector<char> buffer(max_size*n_cpus);
          const int ierr = MPI_Allgather(&my_data[0], max_size, MPI_BYTE,
                                         &buffer[0], max_size, MPI_BYTE,
                                         tr->get_communicator());
          AssertThrowMPI(ierr);

          number_cache.locally_owned_dofs_per_processor.resize (n_cpus);
          number_cache.n_locally_owned_dofs_per_processor.resize (n_cpus);
          for (unsigned int i=0; i<n_cpus; ++i)
            {
              std::stringstream strstr;
              strstr.write(&buffer[i*max_size], max_size);
              // This does not read the whole buffer, when the size is smaller
              // than max_size. Therefore we need to create a new stringstream
              // in each iteration (resetting would be fine too).
              number_cache.locally_owned_dofs_per_processor[i]
              .block_read(strstr);
              number_cache.n_locally_owned_dofs_per_processor[i]
                = number_cache.locally_owned_dofs_per_processor[i].n_elements();
            }

          number_cache.n_global_dofs
            = std::accumulate (number_cache
                               .n_locally_owned_dofs_per_processor.begin(),
                               number_cache
                               .n_locally_owned_dofs_per_processor.end(),
                               static_cast<dealii::types::global_dof_index>(0));

          tr->load_user_flags(user_flags);
        }
#endif

        return number_cache;
      }
    }
  }
}




/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler_policy.inst"


DEAL_II_NAMESPACE_CLOSE
