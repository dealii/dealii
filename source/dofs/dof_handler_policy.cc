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
#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler_policy.h>
#include <deal.II/fe/fe.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

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

      struct Implementation
      {

        /* -------------- distribute_dofs functionality ------------- */

        /**
         * Distribute dofs on the given cell,
         * with new dofs starting with index
         * @p next_free_dof. Return the next
         * unused index number.
         *
         * This function is excluded from the
         * @p distribute_dofs function since
         * it can not be implemented dimension
         * independent.
         */
        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const DoFHandler<1,spacedim> &dof_handler,
                                 const typename DoFHandler<1,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index next_free_dof)
        {

          // distribute dofs of vertices
          if (dof_handler.get_fe().dofs_per_vertex > 0)
            for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
              {
                if (cell->vertex_dof_index (v,0) ==
                    DoFHandler<1,spacedim>::invalid_dof_index)
                  for (unsigned int d=0;
                       d<dof_handler.get_fe().dofs_per_vertex; ++d)
                    {
                      Assert ((cell->vertex_dof_index (v,d) ==
                               DoFHandler<1,spacedim>::invalid_dof_index),
                              ExcInternalError());
                      cell->set_vertex_dof_index (v, d, next_free_dof++);
                    }
                else
                  for (unsigned int d=0;
                       d<dof_handler.get_fe().dofs_per_vertex; ++d)
                    Assert ((cell->vertex_dof_index (v,d) !=
                             DoFHandler<1,spacedim>::invalid_dof_index),
                            ExcInternalError());
              }

          // dofs of line
          for (unsigned int d=0;
               d<dof_handler.get_fe().dofs_per_line; ++d)
            cell->set_dof_index (d, next_free_dof++);

          // note that this cell has been
          // processed
          cell->set_user_flag ();

          return next_free_dof;
        }



        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const DoFHandler<2,spacedim> &dof_handler,
                                 const typename DoFHandler<2,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index next_free_dof)
        {
          if (dof_handler.get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
              // check whether dofs for this
              // vertex have been distributed
              // (only check the first dof)
              if (cell->vertex_dof_index(vertex, 0) == DoFHandler<2,spacedim>::invalid_dof_index)
                for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                  cell->set_vertex_dof_index (vertex, d, next_free_dof++);

          // for the four sides
          if (dof_handler.get_fe().dofs_per_line > 0)
            for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
              {
                const typename DoFHandler<2,spacedim>::line_iterator
                line = cell->line(side);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if (line->dof_index(0) == DoFHandler<2,spacedim>::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_line; ++d)
                    line->set_dof_index (d, next_free_dof++);
              }


          // dofs of quad
          if (dof_handler.get_fe().dofs_per_quad > 0)
            for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_quad; ++d)
              cell->set_dof_index (d, next_free_dof++);


          // note that this cell has been processed
          cell->set_user_flag ();

          return next_free_dof;
        }


        template <int spacedim>
        static
        types::global_dof_index
        distribute_dofs_on_cell (const DoFHandler<3,spacedim> &dof_handler,
                                 const typename DoFHandler<3,spacedim>::active_cell_iterator &cell,
                                 types::global_dof_index next_free_dof)
        {
          if (dof_handler.get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
              // check whether dofs for this
              // vertex have been distributed
              // (only check the first dof)
              if (cell->vertex_dof_index(vertex, 0) == DoFHandler<3,spacedim>::invalid_dof_index)
                for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                  cell->set_vertex_dof_index (vertex, d, next_free_dof++);

          // for the lines
          if (dof_handler.get_fe().dofs_per_line > 0)
            for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
              {
                const typename DoFHandler<3,spacedim>::line_iterator
                line = cell->line(l);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if (line->dof_index(0) == DoFHandler<3,spacedim>::invalid_dof_index)
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

                // distribute dofs if necessary:
                // check whether quad dof is already
                // numbered (check only first dof)
                if (quad->dof_index(0) == DoFHandler<3,spacedim>::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_quad; ++d)
                    quad->set_dof_index (d, next_free_dof++);
              }


          // dofs of hex
          if (dof_handler.get_fe().dofs_per_hex > 0)
            for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_hex; ++d)
              cell->set_dof_index (d, next_free_dof++);


          // note that this cell has been
          // processed
          cell->set_user_flag ();

          return next_free_dof;
        }


        /**
         * Distribute degrees of freedom on all cells, or on cells with the
         * correct subdomain_id if the corresponding argument is not equal to
         * numbers::invalid_subdomain_id. Return the total number of dofs
         * distributed.
         */
        template <int dim, int spacedim>
        static
        types::global_dof_index
        distribute_dofs (const types::global_dof_index offset,
                         const types::subdomain_id subdomain_id,
                         DoFHandler<dim,spacedim> &dof_handler)
        {
          const dealii::Triangulation<dim,spacedim> &tria
            = dof_handler.get_triangulation();
          Assert (tria.n_levels() > 0, ExcMessage("Empty triangulation"));

          // Clear user flags because we will need them. But first we save
          // them and make sure that we restore them later such that at the
          // end of this function the Triangulation will be in the same state
          // as it was at the beginning of this function.
          std::vector<bool> user_flags;
          tria.save_user_flags(user_flags);
          const_cast<dealii::Triangulation<dim,spacedim> &>(tria).clear_user_flags ();

          types::global_dof_index next_free_dof = offset;
          typename DoFHandler<dim,spacedim>::active_cell_iterator
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

          // update the cache used for cell dof indices
          for (cell = dof_handler.begin_active(); cell != endc; ++cell)
            if (!cell->is_artificial())
              cell->update_cell_dof_indices_cache ();

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim,spacedim> &>(tria).load_user_flags(user_flags);

          return next_free_dof;
        }


        /**
         * Distribute dofs on the given
         * cell, with new dofs starting
         * with index
         * @p next_free_dof. Return the
         * next unused index number.
         *
         * This function is excluded from
         * the @p distribute_dofs
         * function since it can not be
         * implemented dimension
         * independent.
         *
         * Note that unlike for the usual
         * dofs, here all cells and not
         * only active ones are allowed.
         */

        // These three function
        // have an unused
        // DoFHandler object as
        // their first
        // argument. Without it,
        // the file was not
        // compileable under gcc
        // 4.4.5 (Debian).
        template <int spacedim>
        static
        unsigned int
        distribute_mg_dofs_on_cell (const DoFHandler<1,spacedim> &,
                                    typename DoFHandler<1,spacedim>::level_cell_iterator &cell,
                                    unsigned int   next_free_dof)
        {
          const unsigned int dim = 1;

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
                      // copy dofs if the neighbor is on
                      // the same level (only then are
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
                      };
                  };

                // otherwise: create dofs newly
                for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index (cell->level(), v, d, next_free_dof++);
              };

          // dofs of line
          if (cell->get_fe().dofs_per_line > 0)
            for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
              cell->set_mg_dof_index (cell->level(), d, next_free_dof++);

          // note that this cell has been processed
          cell->set_user_flag ();

          return next_free_dof;
        }


        template <int spacedim>
        static
        unsigned int
        distribute_mg_dofs_on_cell (const DoFHandler<2,spacedim> &,
                                    typename DoFHandler<2,spacedim>::level_cell_iterator &cell,
                                    unsigned int   next_free_dof)
        {
          const unsigned int dim = 2;
          if (cell->get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
              // check whether dofs for this
              // vertex have been distributed
              // (only check the first dof)
              if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == DoFHandler<2>::invalid_dof_index)
                for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
                  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);

          // for the four sides
          if (cell->get_fe().dofs_per_line > 0)
            for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
              {
                typename DoFHandler<dim,spacedim>::line_iterator line = cell->line(side);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if (line->mg_dof_index(cell->level(), 0) == DoFHandler<2>::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
                    line->set_mg_dof_index (cell->level(), d, next_free_dof++);
              };


          // dofs of quad
          if (cell->get_fe().dofs_per_quad > 0)
            for (unsigned int d=0; d<cell->get_fe().dofs_per_quad; ++d)
              cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


          // note that this cell has been processed
          cell->set_user_flag ();

          return next_free_dof;
        }


        template <int spacedim>
        static
        unsigned int
        distribute_mg_dofs_on_cell (const DoFHandler<3,spacedim> &,
                                    typename DoFHandler<3,spacedim>::level_cell_iterator &cell,
                                    unsigned int   next_free_dof)
        {
          const unsigned int dim = 3;
          if (cell->get_fe().dofs_per_vertex > 0)
            // number dofs on vertices
            for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
              // check whether dofs for this
              // vertex have been distributed
              // (only check the first dof)
              if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == DoFHandler<3>::invalid_dof_index)
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
                if (line->mg_dof_index(cell->level(), 0) == DoFHandler<3>::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
                    line->set_mg_dof_index (cell->level(), d, next_free_dof++);
              };

          // for the quads
          if (cell->get_fe().dofs_per_quad > 0)
            for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
              {
                typename DoFHandler<dim,spacedim>::quad_iterator quad = cell->quad(q);

                // distribute dofs if necessary:
                // check whether line dof is already
                // numbered (check only first dof)
                if (quad->mg_dof_index(cell->level(), 0) == DoFHandler<3>::invalid_dof_index)
                  // if not: distribute dofs
                  for (unsigned int d=0; d<cell->get_fe().dofs_per_quad; ++d)
                    quad->set_mg_dof_index (cell->level(), d, next_free_dof++);
              };


          // dofs of cell
          if (cell->get_fe().dofs_per_hex > 0)
            for (unsigned int d=0; d<cell->get_fe().dofs_per_hex; ++d)
              cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


          // note that this cell has
          // been processed
          cell->set_user_flag ();

          return next_free_dof;
        }


        template <int dim, int spacedim>
        static
        unsigned int
        distribute_dofs_on_level (const unsigned int        offset,
                                  const types::subdomain_id level_subdomain_id,
                                  DoFHandler<dim,spacedim> &dof_handler,
                                  const unsigned int level)
        {
          const dealii::Triangulation<dim,spacedim> &tria
            = dof_handler.get_triangulation();
          Assert (tria.n_levels() > 0, ExcMessage("Empty triangulation"));
          if (level>=tria.n_levels())
            return 0; //this is allowed for multigrid

          // Clear user flags because we will
          // need them. But first we save
          // them and make sure that we
          // restore them later such that at
          // the end of this function the
          // Triangulation will be in the
          // same state as it was at the
          // beginning of this function.
          std::vector<bool> user_flags;
          tria.save_user_flags(user_flags);
          const_cast<dealii::Triangulation<dim,spacedim> &>(tria).clear_user_flags ();

          unsigned int next_free_dof = offset;
          typename DoFHandler<dim,spacedim>::level_cell_iterator
          cell = dof_handler.begin(level),
          endc = dof_handler.end(level);

          for (; cell != endc; ++cell)
            if ((level_subdomain_id == numbers::invalid_subdomain_id)
                ||
                (cell->level_subdomain_id() == level_subdomain_id))
              next_free_dof
                = Implementation::distribute_mg_dofs_on_cell (dof_handler, cell, next_free_dof);

//                                               // update the cache used
//                                               // for cell dof indices
//              for (typename DoFHandler<dim,spacedim>::level_cell_iterator
//                     cell = dof_handler.begin(); cell != dof_handler.end(); ++cell)
//                if (cell->subdomain_id() != numbers::artificial_subdomain_id)
//                  cell->update_cell_dof_indices_cache ();

          // finally restore the user flags
          const_cast<dealii::Triangulation<dim,spacedim> &>(tria).load_user_flags(user_flags);

          return next_free_dof;
        }


        /* --------------------- renumber_dofs functionality ---------------- */


        /**
         * Implementation of the
         * general template of same
         * name.
         *
         * If the second argument
         * has any elements set,
         * elements of the then the
         * vector of new numbers do
         * not relate to the old
         * DoF number but instead
         * to the index of the old
         * DoF number within the
         * set of locally owned
         * DoFs.
         */
        template <int spacedim>
        static
        void
        renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
                       const IndexSet &,
                       DoFHandler<1,spacedim>          &dof_handler,
                       const bool check_validity)
        {
          // note that we can not use cell
          // iterators in this function since
          // then we would renumber the dofs on
          // the interface of two cells more
          // than once. Anyway, this way it's
          // not only more correct but also
          // faster; note, however, that dof
          // numbers may be invalid_dof_index,
          // namely when the appropriate
          // vertex/line/etc is unused
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.vertex_dofs.begin();
               i!=dof_handler.vertex_dofs.end(); ++i)
            if (*i != DoFHandler<1,spacedim>::invalid_dof_index)
              *i = new_numbers[*i];
            else if (check_validity)
              // if index is
              // invalid_dof_index:
              // check if this one
              // really is unused
              Assert (dof_handler.get_triangulation()
                      .vertex_used((i-dof_handler.vertex_dofs.begin()) /
                                   dof_handler.selected_fe->dofs_per_vertex)
                      == false,
                      ExcInternalError ());

          for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
            for (std::vector<types::global_dof_index>::iterator
                 i=dof_handler.levels[level]->dof_object.dofs.begin();
                 i!=dof_handler.levels[level]->dof_object.dofs.end(); ++i)
              if (*i != DoFHandler<1,spacedim>::invalid_dof_index)
                *i = new_numbers[*i];

          // update the cache
          // used for cell dof
          // indices
          for (typename DoFHandler<1,spacedim>::level_cell_iterator
               cell = dof_handler.begin();
               cell != dof_handler.end(); ++cell)
            cell->update_cell_dof_indices_cache ();
        }

        template <int spacedim>
        static
        void
        renumber_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                          const IndexSet &indices,
                          DoFHandler<1,spacedim>          &dof_handler,
                          const unsigned int level,
                          const bool check_validity)
        {
          for (typename std::vector<typename DoFHandler<1,spacedim>::MGVertexDoFs>::iterator
               i=dof_handler.mg_vertex_dofs.begin();
               i!=dof_handler.mg_vertex_dofs.end();
               ++i)
            // if the present vertex lives on
            // the current level
            if ((i->get_coarsest_level() <= level) &&
                (i->get_finest_level() >= level))
              for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                {
                  dealii::types::global_dof_index idx = i->get_index (level, d);
                  if (idx != DoFHandler<1>::invalid_dof_index)
                    i->set_index (level, d,
                                  (indices.n_elements() == 0)?
                                  (new_numbers[idx]) :
                                  (new_numbers[indices.index_within_set(idx)]));

                  if (check_validity)
                    Assert(idx != DoFHandler<1>::invalid_dof_index, ExcInternalError ());
                }


          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.mg_levels[level]->dof_object.dofs.begin();
               i!=dof_handler.mg_levels[level]->dof_object.dofs.end();
               ++i)
            {
              if (*i != DoFHandler<1>::invalid_dof_index)
                {
                  Assert(*i<new_numbers.size(), ExcInternalError());
                  *i = (indices.n_elements() == 0)?
                       (new_numbers[*i]) :
                       (new_numbers[indices.index_within_set(*i)]);
                }
            }
        }


        template <int spacedim>
        static
        void
        renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
                       const IndexSet &indices,
                       DoFHandler<2,spacedim> &dof_handler,
                       const bool check_validity)
        {
          // note that we can not use cell
          // iterators in this function since
          // then we would renumber the dofs on
          // the interface of two cells more
          // than once. Anyway, this way it's
          // not only more correct but also
          // faster; note, however, that dof
          // numbers may be invalid_dof_index,
          // namely when the appropriate
          // vertex/line/etc is unused
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.vertex_dofs.begin();
               i!=dof_handler.vertex_dofs.end(); ++i)
            if (*i != DoFHandler<2,spacedim>::invalid_dof_index)
              *i = (indices.n_elements() == 0)?
                   (new_numbers[*i]) :
                   (new_numbers[indices.index_within_set(*i)]);
            else if (check_validity)
              // if index is invalid_dof_index:
              // check if this one really is
              // unused
              Assert (dof_handler.get_triangulation()
                      .vertex_used((i-dof_handler.vertex_dofs.begin()) /
                                   dof_handler.selected_fe->dofs_per_vertex)
                      == false,
                      ExcInternalError ());

          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.faces->lines.dofs.begin();
               i!=dof_handler.faces->lines.dofs.end(); ++i)
            if (*i != DoFHandler<2,spacedim>::invalid_dof_index)
              *i = ((indices.n_elements() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);

          for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
            {
              for (std::vector<types::global_dof_index>::iterator
                   i=dof_handler.levels[level]->dof_object.dofs.begin();
                   i!=dof_handler.levels[level]->dof_object.dofs.end(); ++i)
                if (*i != DoFHandler<2,spacedim>::invalid_dof_index)
                  *i = ((indices.n_elements() == 0) ?
                        new_numbers[*i] :
                        new_numbers[indices.index_within_set(*i)]);
            }

          // update the cache
          // used for cell dof
          // indices
          for (typename DoFHandler<2,spacedim>::level_cell_iterator
               cell = dof_handler.begin();
               cell != dof_handler.end(); ++cell)
            cell->update_cell_dof_indices_cache ();
        }

        template <int spacedim>
        static
        void
        renumber_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                          const IndexSet &indices,
                          DoFHandler<2,spacedim>          &dof_handler,
                          const unsigned int level,
                          const bool check_validity)
        {
          if (level>=dof_handler.get_triangulation().n_levels())
            return;
          for (typename std::vector<typename DoFHandler<2,spacedim>::MGVertexDoFs>::iterator i=dof_handler.mg_vertex_dofs.begin();
               i!=dof_handler.mg_vertex_dofs.end(); ++i)
            // if the present vertex lives on
            // the present level
            if ((i->get_coarsest_level() <= level) &&
                (i->get_finest_level() >= level))
              for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                {
                  dealii::types::global_dof_index idx =i->get_index (level, d/*,                    dof_handler.get_fe().dofs_per_vertex*/);
                  if (idx != DoFHandler<1>::invalid_dof_index)
                    i->set_index (level, d/*, dof_handler.get_fe().dofs_per_vertex*/,
                                  ((indices.n_elements() == 0) ?
                                   new_numbers[idx] :
                                   new_numbers[indices.index_within_set(idx)]));

                  if (check_validity)
                    Assert(idx != DoFHandler<2>::invalid_dof_index, ExcInternalError ());
                }

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
                          dealii::types::global_dof_index idx = cell->line(l)->mg_dof_index(level, d);
                          if (idx != DoFHandler<1>::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index (level, d, ((indices.n_elements() == 0) ?
                                                                        new_numbers[idx] :
                                                                        new_numbers[indices.index_within_set(idx)]));
                          if (check_validity)
                            Assert(idx != DoFHandler<2>::invalid_dof_index, ExcInternalError ());
                        }
                      cell->line(l)->clear_user_flag();
                    }
              // finally, restore user flags
              const_cast<dealii::Triangulation<2,spacedim> &>(dof_handler.get_triangulation()).load_user_flags (user_flags);
            }

          for (std::vector<types::global_dof_index>::iterator i=dof_handler.mg_levels[level]->dof_object.dofs.begin();
               i!=dof_handler.mg_levels[level]->dof_object.dofs.end(); ++i)
            {
              if (*i != DoFHandler<2>::invalid_dof_index)
                {
                  Assert(*i<new_numbers.size(), ExcInternalError());
                  *i = ((indices.n_elements() == 0) ?
                        new_numbers[*i] :
                        new_numbers[indices.index_within_set(*i)]);
                }
            }

        }


        template <int spacedim>
        static
        void
        renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
                       const IndexSet &indices,
                       DoFHandler<3,spacedim>          &dof_handler,
                       const bool check_validity)
        {
          // note that we can not use cell
          // iterators in this function since
          // then we would renumber the dofs on
          // the interface of two cells more
          // than once. Anyway, this way it's
          // not only more correct but also
          // faster; note, however, that dof
          // numbers may be invalid_dof_index,
          // namely when the appropriate
          // vertex/line/etc is unused
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.vertex_dofs.begin();
               i!=dof_handler.vertex_dofs.end(); ++i)
            if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
              *i = ((indices.n_elements() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);
            else if (check_validity)
              // if index is invalid_dof_index:
              // check if this one really is
              // unused
              Assert (dof_handler.get_triangulation()
                      .vertex_used((i-dof_handler.vertex_dofs.begin()) /
                                   dof_handler.selected_fe->dofs_per_vertex)
                      == false,
                      ExcInternalError ());

          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.faces->lines.dofs.begin();
               i!=dof_handler.faces->lines.dofs.end(); ++i)
            if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
              *i = ((indices.n_elements() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);
          for (std::vector<types::global_dof_index>::iterator
               i=dof_handler.faces->quads.dofs.begin();
               i!=dof_handler.faces->quads.dofs.end(); ++i)
            if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
              *i = ((indices.n_elements() == 0) ?
                    new_numbers[*i] :
                    new_numbers[indices.index_within_set(*i)]);

          for (unsigned int level=0; level<dof_handler.levels.size(); ++level)
            {
              for (std::vector<types::global_dof_index>::iterator
                   i=dof_handler.levels[level]->dof_object.dofs.begin();
                   i!=dof_handler.levels[level]->dof_object.dofs.end(); ++i)
                if (*i != DoFHandler<3,spacedim>::invalid_dof_index)
                  *i = ((indices.n_elements() == 0) ?
                        new_numbers[*i] :
                        new_numbers[indices.index_within_set(*i)]);
            }

          // update the cache
          // used for cell dof
          // indices
          for (typename DoFHandler<3,spacedim>::level_cell_iterator
               cell = dof_handler.begin();
               cell != dof_handler.end(); ++cell)
            cell->update_cell_dof_indices_cache ();
        }

        template <int spacedim>
        static
        void
        renumber_mg_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                          const IndexSet &indices,
                          DoFHandler<3,spacedim> &dof_handler,
                          const unsigned int level,
                          const bool check_validity)
        {
          if (level>=dof_handler.get_triangulation().n_levels())
            return;
          for (typename std::vector<typename DoFHandler<3,spacedim>::MGVertexDoFs>::iterator i=dof_handler.mg_vertex_dofs.begin();
               i!=dof_handler.mg_vertex_dofs.end(); ++i)
            // if the present vertex lives on
            // the present level
            if ((i->get_coarsest_level() <= level) &&
                (i->get_finest_level() >= level))
              for (unsigned int d=0; d<dof_handler.get_fe().dofs_per_vertex; ++d)
                {
                  dealii::types::global_dof_index idx =i->get_index (level, d);
                  if (idx != DoFHandler<3>::invalid_dof_index)
                    i->set_index (level, d,
                                  ((indices.n_elements() == 0) ?
                                   new_numbers[idx] :
                                   new_numbers[indices.index_within_set(idx)]));

                  if (check_validity)
                    Assert(idx != DoFHandler<3>::invalid_dof_index, ExcInternalError ());
                }

          if (dof_handler.get_fe().dofs_per_line > 0 ||
              dof_handler.get_fe().dofs_per_quad > 0)
            {
              // save user flags as they will be modified
              std::vector<bool> user_flags;
              dof_handler.get_triangulation().save_user_flags(user_flags);
              const_cast<dealii::Triangulation<3,spacedim> &>(dof_handler.get_triangulation()).clear_user_flags ();

              // flag all lines adjacent to cells of the current level, as
              // those lines logically belong to the same level as the cell,
              // at least for isotropic refinement
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
                          dealii::types::global_dof_index idx = cell->line(l)->mg_dof_index(level, d);
                          if (idx != DoFHandler<3>::invalid_dof_index)
                            cell->line(l)->set_mg_dof_index (level, d, ((indices.n_elements() == 0) ?
                                                                        new_numbers[idx] :
                                                                        new_numbers[indices.index_within_set(idx)]));
                          if (check_validity)
                            Assert(idx != DoFHandler<3>::invalid_dof_index, ExcInternalError ());
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
                          dealii::types::global_dof_index idx = cell->quad(l)->mg_dof_index(level, d);
                          if (idx != DoFHandler<1>::invalid_dof_index)
                            cell->quad(l)->set_mg_dof_index (level, d, ((indices.n_elements() == 0) ?
                                                                        new_numbers[idx] :
                                                                        new_numbers[indices.index_within_set(idx)]));
                          if (check_validity)
                            Assert(idx != DoFHandler<3>::invalid_dof_index, ExcInternalError ());
                        }
                      cell->quad(l)->clear_user_flag();
                    }

              // finally, restore user flags
              const_cast<dealii::Triangulation<3,spacedim> &>(dof_handler.get_triangulation()).load_user_flags (user_flags);
            }

          for (std::vector<types::global_dof_index>::iterator i=dof_handler.mg_levels[level]->dof_object.dofs.begin();
               i!=dof_handler.mg_levels[level]->dof_object.dofs.end(); ++i)
            {
              if (*i != DoFHandler<3>::invalid_dof_index)
                {
                  Assert(*i<new_numbers.size(), ExcInternalError());
                  *i = ((indices.n_elements() == 0) ?
                        new_numbers[*i] :
                        new_numbers[indices.index_within_set(*i)]);
                }
            }
        }


      };



      /* --------------------- class PolicyBase ---------------- */

      template <int dim, int spacedim>
      PolicyBase<dim,spacedim>::~PolicyBase ()
      {}


      /* --------------------- class Sequential ---------------- */


      template <int dim, int spacedim>
      void
      Sequential<dim,spacedim>::
      distribute_dofs (DoFHandler<dim,spacedim> &dof_handler,
                       NumberCache &number_cache_current ) const
      {
        const types::global_dof_index n_dofs =
          Implementation::distribute_dofs (0,
                                           numbers::invalid_subdomain_id,
                                           dof_handler);

        // now set the elements of the
        // number cache appropriately
        NumberCache number_cache;
        number_cache.n_global_dofs        = n_dofs;
        number_cache.n_locally_owned_dofs = number_cache.n_global_dofs;

        number_cache.locally_owned_dofs
          = IndexSet (number_cache.n_global_dofs);
        number_cache.locally_owned_dofs.add_range (0,
                                                   number_cache.n_global_dofs);
        number_cache.locally_owned_dofs.compress();

        number_cache.n_locally_owned_dofs_per_processor
          = std::vector<types::global_dof_index> (1,
                                                  number_cache.n_global_dofs);

        number_cache.locally_owned_dofs_per_processor
          = std::vector<IndexSet> (1,
                                   number_cache.locally_owned_dofs);
        number_cache_current = number_cache;
      }


      template <int dim, int spacedim>
      void
      Sequential<dim,spacedim>::
      distribute_mg_dofs (DoFHandler<dim,spacedim> &dof_handler,
                          std::vector<NumberCache> &number_caches) const
      {
        std::vector<bool> user_flags;

        dof_handler.get_triangulation().save_user_flags (user_flags);
        const_cast<dealii::Triangulation<dim, spacedim>&>(dof_handler.get_triangulation()).clear_user_flags ();

        for (unsigned int level = 0; level < dof_handler.get_triangulation().n_levels(); ++level)
          {
            types::global_dof_index next_free_dof = Implementation::distribute_dofs_on_level(0, numbers::invalid_subdomain_id, dof_handler, level);

            number_caches[level].n_global_dofs = next_free_dof;
            number_caches[level].n_locally_owned_dofs = next_free_dof;
            number_caches[level].locally_owned_dofs = complete_index_set(next_free_dof);
            number_caches[level].locally_owned_dofs_per_processor.resize(1);
            number_caches[level].locally_owned_dofs_per_processor[0] = complete_index_set(next_free_dof);
            number_caches[level].n_locally_owned_dofs_per_processor.resize(1);
            number_caches[level].n_locally_owned_dofs_per_processor[0] = next_free_dof;
          }
        const_cast<dealii::Triangulation<dim, spacedim>&>(dof_handler.get_triangulation()).load_user_flags (user_flags);
      }

      template <int dim, int spacedim>
      void
      Sequential<dim,spacedim>::
      renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
                     dealii::DoFHandler<dim,spacedim> &dof_handler,
                     NumberCache &number_cache_current) const
      {
        Implementation::renumber_dofs (new_numbers, IndexSet(0),
                                       dof_handler, true);

        // in the sequential case,
        // the number cache should
        // not have changed but we
        // have to set the elements
        // of the structure
        // appropriately anyway
        NumberCache number_cache;
        number_cache.n_global_dofs        = dof_handler.n_dofs();
        number_cache.n_locally_owned_dofs = number_cache.n_global_dofs;

        number_cache.locally_owned_dofs
          = IndexSet (number_cache.n_global_dofs);
        number_cache.locally_owned_dofs.add_range (0,
                                                   number_cache.n_global_dofs);
        number_cache.locally_owned_dofs.compress();

        number_cache.n_locally_owned_dofs_per_processor
          = std::vector<types::global_dof_index> (1,
                                                  number_cache.n_global_dofs);

        number_cache.locally_owned_dofs_per_processor
          = std::vector<IndexSet> (1,
                                   number_cache.locally_owned_dofs);
        number_cache_current = number_cache;
      }

      /* --------------------- class ParallelShared ---------------- */

      template <int dim, int spacedim>
      void
      ParallelShared<dim,spacedim>::
      distribute_dofs (DoFHandler<dim,spacedim> &dof_handler,
                       NumberCache &number_cache) const
      {
        // If the underlying shared::Tria allows artifical cells, we need to do
        // some tricks here to make Sequential algorithms play nicely.
        // Namely, we first restore original partition (without artificial cells)
        // and then turn artificial cells on at the end of this function.
        const parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>*> (&dof_handler.get_triangulation()));
        Assert(tr != 0, ExcInternalError());
        typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator
        cell = dof_handler.get_triangulation().begin_active(),
        endc = dof_handler.get_triangulation().end();
        std::vector<types::subdomain_id> current_subdomain_ids(tr->n_active_cells());
        const std::vector<types::subdomain_id> &true_subdomain_ids = tr->get_true_subdomain_ids_of_cells();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            {
              current_subdomain_ids[index] = cell->subdomain_id();
              cell->set_subdomain_id(true_subdomain_ids[index]);
            }

        Sequential<dim,spacedim>::distribute_dofs (dof_handler,number_cache);
        DoFRenumbering::subdomain_wise (dof_handler);
        // dofrenumbering will reset subdomains, this is ugly but we need to do it again:
        cell = tr->begin_active();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            cell->set_subdomain_id(true_subdomain_ids[index]);

        number_cache.locally_owned_dofs_per_processor = DoFTools::locally_owned_dofs_per_subdomain (dof_handler);
        number_cache.locally_owned_dofs = number_cache.locally_owned_dofs_per_processor[dof_handler.get_triangulation().locally_owned_subdomain()];
        number_cache.n_locally_owned_dofs_per_processor.resize (number_cache.locally_owned_dofs_per_processor.size());
        for (unsigned int i = 0; i < number_cache.n_locally_owned_dofs_per_processor.size(); i++)
          number_cache.n_locally_owned_dofs_per_processor[i] = number_cache.locally_owned_dofs_per_processor[i].n_elements();
        number_cache.n_locally_owned_dofs = number_cache.n_locally_owned_dofs_per_processor[dof_handler.get_triangulation().locally_owned_subdomain()];

        // restore current subdomain ids
        cell = tr->begin_active();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            cell->set_subdomain_id(current_subdomain_ids[index]);
      }

      template <int dim, int spacedim>
      void
      ParallelShared<dim,spacedim>::
      distribute_mg_dofs (DoFHandler<dim,spacedim> &dof_handler,
                          std::vector<NumberCache> &number_caches) const
      {
        // first, call the sequential function to distribute dofs
        Sequential<dim,spacedim>:: distribute_mg_dofs (dof_handler, number_caches);
        // now we need to update the number cache.
        // This part is not yet implemented.
        AssertThrow(false,ExcNotImplemented());
      }

      template <int dim, int spacedim>
      void
      ParallelShared<dim,spacedim>::
      renumber_dofs (const std::vector<types::global_dof_index> &new_numbers,
                     dealii::DoFHandler<dim,spacedim> &dof_handler,
                     NumberCache &number_cache) const
      {

#ifndef DEAL_II_WITH_MPI
        (void)new_numbers;
        (void)dof_handler;
        (void)number_cache;
        Assert (false, ExcNotImplemented());
#else
        // Similar to distribute_dofs() we need to have a special treatment in
        // case artificial cells are present.
        const parallel::shared::Triangulation<dim, spacedim> *tr =
          (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>*> (&dof_handler.get_triangulation()));
        Assert(tr != 0, ExcInternalError());
        typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator
        cell = dof_handler.get_triangulation().begin_active(),
        endc = dof_handler.get_triangulation().end();
        std::vector<types::subdomain_id> current_subdomain_ids(tr->n_active_cells());
        const std::vector<types::subdomain_id> &true_subdomain_ids = tr->get_true_subdomain_ids_of_cells();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            {
              current_subdomain_ids[index] = cell->subdomain_id();
              cell->set_subdomain_id(true_subdomain_ids[index]);
            }

        std::vector<types::global_dof_index> global_gathered_numbers (dof_handler.n_dofs (), 0);
        // as we call DoFRenumbering::subdomain_wise (dof_handler) from distribute_dofs(),
        // we need to support sequential-like input.
        // Distributed-like input from, for example, component_wise renumbering is also supported.
        if (new_numbers.size () == dof_handler.n_dofs ())
          {
            global_gathered_numbers = new_numbers;
          }
        else
          {
            Assert(new_numbers.size() == dof_handler.locally_owned_dofs().n_elements(),
                   ExcInternalError());
            const unsigned int n_cpu = Utilities::MPI::n_mpi_processes (tr->get_communicator ());
            std::vector<types::global_dof_index> gathered_new_numbers (dof_handler.n_dofs (), 0);
            Assert(Utilities::MPI::this_mpi_process (tr->get_communicator ()) ==
                   dof_handler.get_triangulation().locally_owned_subdomain (),
                   ExcInternalError())

            //gather new numbers among processors into one vector
            {
              std::vector<types::global_dof_index> new_numbers_copy (new_numbers);
              // displs:
              // Entry i specifies the displacement (relative to recvbuf )
              // at which to place the incoming data from process i
              // rcounts:
              // containing the number of elements that are to be received from each process
              std::vector<int> displs(n_cpu),
                  rcounts(n_cpu);
              types::global_dof_index shift = 0;
              //set rcounts based on new_numbers:
              int cur_count = new_numbers_copy.size ();
              MPI_Allgather (&cur_count,  1, MPI_INT,
                             &rcounts[0], 1, MPI_INT,
                             tr->get_communicator ());

              for (unsigned int i = 0; i < n_cpu; i++)
                {
                  displs[i]  = shift;
                  shift     += rcounts[i];
                }
              Assert(((int)new_numbers_copy.size()) ==
                     rcounts[Utilities::MPI::this_mpi_process (tr->get_communicator ())],
                     ExcInternalError());
              MPI_Allgatherv (&new_numbers_copy[0],     new_numbers_copy.size (),
                              DEAL_II_DOF_INDEX_MPI_TYPE,
                              &gathered_new_numbers[0], &rcounts[0],
                              &displs[0],
                              DEAL_II_DOF_INDEX_MPI_TYPE,
                              tr->get_communicator ());
            }

            // put new numbers according to the current locally_owned_dofs_per_processor IndexSets
            types::global_dof_index shift = 0;
            // flag_1 and flag_2 are
            // used to control that there is a
            // one-to-one relation between old and new DoFs.
            std::vector<unsigned int> flag_1 (dof_handler.n_dofs (), 0),
                flag_2 (dof_handler.n_dofs (), 0);
            for (unsigned int i = 0; i < n_cpu; i++)
              {
                const IndexSet &iset =
                  number_cache.locally_owned_dofs_per_processor[i];
                for (types::global_dof_index ind = 0;
                     ind < iset.n_elements (); ind++)
                  {
                    const types::global_dof_index target = iset.nth_index_in_set (ind);
                    const types::global_dof_index value  = gathered_new_numbers[shift + ind];
                    Assert(target < dof_handler.n_dofs(), ExcInternalError());
                    Assert(value  < dof_handler.n_dofs(), ExcInternalError());
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
        Sequential<dim, spacedim>::renumber_dofs (global_gathered_numbers, dof_handler, number_cache);
        // correct number_cache:
        number_cache.locally_owned_dofs_per_processor =
          DoFTools::locally_owned_dofs_per_subdomain (dof_handler);
        number_cache.locally_owned_dofs =
          number_cache.locally_owned_dofs_per_processor[dof_handler.get_triangulation().locally_owned_subdomain ()];
        // sequential renumbering returns a vector of size 1 here,
        // correct this:
        number_cache.n_locally_owned_dofs_per_processor.resize(number_cache.locally_owned_dofs_per_processor.size());
        for (unsigned int i = 0;
             i < number_cache.n_locally_owned_dofs_per_processor.size (); i++)
          number_cache.n_locally_owned_dofs_per_processor[i] = number_cache.locally_owned_dofs_per_processor[i].n_elements ();

        number_cache.n_locally_owned_dofs =
          number_cache.n_locally_owned_dofs_per_processor[dof_handler.get_triangulation().locally_owned_subdomain ()];

        // restore artificial cells
        cell = tr->begin_active();
        if (tr->with_artificial_cells())
          for (unsigned int index=0; cell != endc; cell++, index++)
            cell->set_subdomain_id(current_subdomain_ids[index]);
#endif
      }

      /* --------------------- class ParallelDistributed ---------------- */

#ifdef DEAL_II_WITH_P4EST

      namespace
      {
        template <int dim>
        struct types
        {

          /**
           * A list of tree+quadrant and optionally their dof indices. This
           * data structure supports (de)serialization into an array of bytes
           * to allow transfer over MPI.
           *
           * If the dofs array exists, it has the format num_dofindices of
           * quadrant 0, followed by
           * num_dofindices indices,
           * num_dofindices of quadrant 1,
           * ...
           */
          struct cellinfo
          {
            std::vector<unsigned int> tree_index;
            std::vector<typename dealii::internal::p4est::types<dim>::quadrant> quadrants;
            std::vector<dealii::types::global_dof_index> dofs;

            unsigned int bytes_for_buffer () const
            {
              return (sizeof(unsigned int)*2 +
                      tree_index.size() * sizeof(unsigned int) +
                      quadrants.size() * sizeof(typename dealii::internal::p4est
                                                ::types<dim>::quadrant) +
                      dofs.size() * sizeof(dealii::types::global_dof_index));
            }

            void pack_data (std::vector<char> &buffer) const
            {
              buffer.resize(bytes_for_buffer());

              char *ptr = &buffer[0];

              const unsigned int num_cells = tree_index.size();
              const unsigned int num_dofs = dofs.size();
              std::memcpy(ptr, &num_cells, sizeof(unsigned int));
              ptr += sizeof(unsigned int);
              std::memcpy(ptr, &num_dofs, sizeof(unsigned int));
              ptr += sizeof(unsigned int);

              std::memcpy(ptr,
                          &tree_index[0],
                          num_cells*sizeof(unsigned int));
              ptr += num_cells*sizeof(unsigned int);

              std::memcpy(ptr,
                          &quadrants[0],
                          num_cells * sizeof(typename dealii::internal::p4est::
                                             types<dim>::quadrant));
              ptr += num_cells*sizeof(typename dealii::internal::p4est::types<dim>::
                                      quadrant);

              if (num_dofs>0)
                std::memcpy(ptr,
                            &dofs[0],
                            dofs.size() * sizeof(dealii::types::global_dof_index));
              ptr += dofs.size() * sizeof(dealii::types::global_dof_index);

              Assert (ptr == &buffer[0]+buffer.size(),
                      ExcInternalError());
            }

            void unpack_data(const std::vector<char> &buffer)
            {
              const char *ptr = &buffer[0];
              unsigned int num_cells;
              unsigned int num_dofs;
              memcpy(&num_cells, ptr, sizeof(unsigned int));
              ptr+=sizeof(unsigned int);
              memcpy(&num_dofs, ptr, sizeof(unsigned int));
              ptr+=sizeof(unsigned int);

              tree_index.resize(num_cells);
              std::memcpy(&tree_index[0],
                          ptr,
                          num_cells*sizeof(unsigned int));
              ptr += num_cells*sizeof(unsigned int);

              quadrants.resize(num_cells);
              std::memcpy(&quadrants[0],
                          ptr,
                          num_cells * sizeof(typename dealii::internal::p4est::
                                             types<dim>::quadrant));
              ptr += num_cells*sizeof(typename dealii::internal::p4est::types<dim>::
                                      quadrant);

              dofs.resize(num_dofs);
              if (num_dofs>0)
                std::memcpy(&dofs[0],
                            ptr,
                            dofs.size() * sizeof(dealii::types::global_dof_index));
              ptr += dofs.size() * sizeof(dealii::types::global_dof_index);

              Assert (ptr == &buffer[0]+buffer.size(),
                      ExcInternalError());
            }

          };

        };



        template <int dim, int spacedim>
        void
        fill_dofindices_recursively (const typename parallel::distributed::Triangulation<dim,spacedim> &tria,
                                     const unsigned int tree_index,
                                     const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
                                     const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
                                     const std::map<unsigned int, std::set<dealii::types::subdomain_id> > &vertices_with_ghost_neighbors,
                                     std::map<dealii::types::subdomain_id, typename types<dim>::cellinfo> &needs_to_get_cell)
        {
          // see if we have to
          // recurse...
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

          // we're at a leaf cell. see if
          // the cell is flagged as
          // interesting. note that we
          // have only flagged our own
          // cells before
          if (dealii_cell->user_flag_set() && !dealii_cell->is_ghost())
            {
              Assert (!dealii_cell->is_artificial(), ExcInternalError());

              // check each vertex if
              // it is interesting and
              // push dofindices if yes
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
                  // this cell's dof_indices
                  // need to be sent to
                  // someone
                  std::vector<dealii::types::global_dof_index>
                  local_dof_indices (dealii_cell->get_fe().dofs_per_cell);
                  dealii_cell->get_dof_indices (local_dof_indices);

                  for (std::set<dealii::types::subdomain_id>::iterator it=send_to.begin();
                       it!=send_to.end(); ++it)
                    {
                      const dealii::types::subdomain_id subdomain = *it;

                      // get an iterator
                      // to what needs to
                      // be sent to that
                      // subdomain (if
                      // already exists),
                      // or create such
                      // an object
                      typename std::map<dealii::types::subdomain_id, typename types<dim>::cellinfo>::iterator
                      p
                        = needs_to_get_cell.insert (std::make_pair(subdomain,
                                                                   typename types<dim>::cellinfo()))
                          .first;

                      p->second.tree_index.push_back(tree_index);
                      p->second.quadrants.push_back(p4est_cell);

                      p->second.dofs.push_back(dealii_cell->get_fe().dofs_per_cell);
                      p->second.dofs.insert(p->second.dofs.end(),
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
          typename types<dim>::cellinfo &cell_info)
        {
          if (internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              // why would somebody request a cell that is not ours?
              Assert(dealii_cell->level_subdomain_id()==tria.locally_owned_subdomain(), ExcInternalError());


              std::vector<dealii::types::global_dof_index>
              local_dof_indices (dealii_cell->get_fe().dofs_per_cell);
              dealii_cell->get_mg_dof_indices (local_dof_indices);

              cell_info.dofs.push_back(dealii_cell->get_fe().dofs_per_cell);
              cell_info.dofs.insert(cell_info.dofs.end(),
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
                                                         quadrant, cell_info);
        }


        template <int dim, int spacedim>
        void
        find_marked_mg_ghost_cells_recursively(const typename parallel::distributed::Triangulation<dim,spacedim> &tria,
                                               const unsigned int tree_index,
                                               const typename DoFHandler<dim,spacedim>::level_cell_iterator &dealii_cell,
                                               const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
                                               std::map<dealii::types::subdomain_id, typename types<dim>::cellinfo> &neighbor_cell_list)
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
              typename types<dim>::cellinfo &cell_info = neighbor_cell_list[dealii_cell->level_subdomain_id()];

              cell_info.tree_index.push_back(tree_index);
              cell_info.quadrants.push_back(p4est_cell);
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
                if (dofs[i] != DoFHandler<dim,spacedim>::invalid_dof_index)
                  {
                    Assert((dof_indices[i] ==
                            (DoFHandler<dim,spacedim>::invalid_dof_index))
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
          std::set<unsigned int> level_ghost_owners = tria.level_ghost_owners();
          typedef std::map<dealii::types::subdomain_id, typename types<dim>::cellinfo> cellmap_t;
          cellmap_t neighbor_cell_list;
          for (std::set<unsigned int>::iterator it = level_ghost_owners.begin();
               it != level_ghost_owners.end();
               ++it)
            neighbor_cell_list.insert(std::make_pair(*it, typename types<dim>::cellinfo()));

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
              // pack all the data into
              // the buffer for this
              // recipient and send
              // it. keep data around
              // till we can make sure
              // that the packet has been
              // received
              it->second.pack_data (sendbuffers[idx]);
              MPI_Isend(sendbuffers[idx].data(), sendbuffers[idx].size(),
                        MPI_BYTE, it->first,
                        1100101, tria.get_communicator(), &requests[idx]);
            }

          //* receive requests and reply
          std::vector<std::vector<char> > reply_buffers (level_ghost_owners.size());
          std::vector<MPI_Request> reply_requests (level_ghost_owners.size());

          for (unsigned int idx=0; idx<level_ghost_owners.size(); ++idx)
            {
              std::vector<char> receive;
              typename types<dim>::cellinfo cellinfo;

              MPI_Status status;
              int len;
              MPI_Probe(MPI_ANY_SOURCE, 1100101, tria.get_communicator(), &status);
              MPI_Get_count(&status, MPI_BYTE, &len);
              receive.resize(len);

              char *ptr = &receive[0];
              MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                       tria.get_communicator(), &status);

              cellinfo.unpack_data(receive);

              // store the dof indices for each cell
              for (unsigned int c=0; c<cellinfo.tree_index.size(); ++c)
                {
                  typename DoFHandler<dim,spacedim>::level_cell_iterator
                  cell (&dof_handler.get_triangulation(),
                        0,
                        p4est_tree_to_coarse_cell_permutation[cellinfo.tree_index[c]],
                        &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  get_mg_dofindices_recursively<dim,spacedim> (tria,
                                                               p4est_coarse_cell,
                                                               cell,
                                                               cellinfo.quadrants[c],
                                                               cellinfo);
                }

              //send reply
              cellinfo.pack_data(reply_buffers[idx]);
              MPI_Isend(&(reply_buffers[idx])[0], reply_buffers[idx].size(),
                        MPI_BYTE, status.MPI_SOURCE,
                        1100102, tria.get_communicator(), &reply_requests[idx]);
            }

          // * finally receive the replies
          for (unsigned int idx=0; idx<level_ghost_owners.size(); ++idx)
            {
              std::vector<char> receive;
              typename types<dim>::cellinfo cellinfo;

              MPI_Status status;
              int len;
              MPI_Probe(MPI_ANY_SOURCE, 1100102, tria.get_communicator(), &status);
              MPI_Get_count(&status, MPI_BYTE, &len);
              receive.resize(len);

              char *ptr = &receive[0];
              MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                       tria.get_communicator(), &status);

              cellinfo.unpack_data(receive);
              if (cellinfo.tree_index.size()==0)
                continue;

              // set the dof indices for each cell
              dealii::types::global_dof_index *dofs = cellinfo.dofs.data();
              for (unsigned int c=0; c<cellinfo.tree_index.size(); ++c, dofs+=1+dofs[0])
                {
                  typename DoFHandler<dim,spacedim>::level_cell_iterator
                  cell (&tria,
                        0,
                        p4est_tree_to_coarse_cell_permutation[cellinfo.tree_index[c]],
                        &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  Assert(cell->get_fe().dofs_per_cell==dofs[0], ExcInternalError());

                  set_mg_dofindices_recursively<dim,spacedim> (tria,
                                                               p4est_coarse_cell,
                                                               cell,
                                                               cellinfo.quadrants[c],
                                                               dofs+1);
                }
            }

          // complete all sends, so that we can
          // safely destroy the buffers.
          if (requests.size() > 0)
            MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
          if (reply_requests.size() > 0)
            MPI_Waitall(reply_requests.size(), &reply_requests[0], MPI_STATUSES_IGNORE);

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
                if (dofs[i] != DoFHandler<dim,spacedim>::invalid_dof_index)
                  {
                    Assert((dof_indices[i] ==
                            (DoFHandler<dim,spacedim>::invalid_dof_index))
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
          Assert (tr != 0, ExcInternalError());

          // now collect cells and their
          // dof_indices for the
          // interested neighbors
          typedef
          std::map<dealii::types::subdomain_id, typename types<dim>::cellinfo>
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


          //sending
          std::vector<std::vector<char> > sendbuffers (needs_to_get_cells.size());
          std::vector<std::vector<char> >::iterator buffer = sendbuffers.begin();
          std::vector<MPI_Request> requests (needs_to_get_cells.size());

          unsigned int idx=0;

          for (typename cellmap_t::iterator it=needs_to_get_cells.begin();
               it!=needs_to_get_cells.end();
               ++it, ++buffer, ++idx)
            {
              const unsigned int num_cells = it->second.tree_index.size();
              (void)num_cells;

              Assert(num_cells==it->second.quadrants.size(), ExcInternalError());
              Assert(num_cells>0, ExcInternalError());

              // pack all the data into
              // the buffer for this
              // recipient and send
              // it. keep data around
              // till we can make sure
              // that the packet has been
              // received
              it->second.pack_data (*buffer);
              MPI_Isend(&(*buffer)[0], buffer->size(),
                        MPI_BYTE, it->first,
                        123, tr->get_communicator(), &requests[idx]);
            }


          // mark all own cells, that miss some
          // dof_data and collect the neighbors
          // that are going to send stuff to us
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
                                     DoFHandler<dim,spacedim>::invalid_dof_index))
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
              MPI_Probe(MPI_ANY_SOURCE, 123, tr->get_communicator(), &status);
              MPI_Get_count(&status, MPI_BYTE, &len);
              receive.resize(len);

              char *ptr = &receive[0];
              MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                       tr->get_communicator(), &status);

              typename types<dim>::cellinfo cellinfo;
              cellinfo.unpack_data(receive);
              unsigned int cells = cellinfo.tree_index.size();
              dealii::types::global_dof_index *dofs = cellinfo.dofs.data();

              // the dofs pointer contains for each cell the number of dofs
              // on that cell (dofs[0]) followed by the dof indices itself.
              for (unsigned int c=0; c<cells; ++c, dofs+=1+dofs[0])
                {
                  typename DoFHandler<dim,spacedim>::level_cell_iterator
                  cell (&dof_handler.get_triangulation(),
                        0,
                        p4est_tree_to_coarse_cell_permutation[cellinfo.tree_index[c]],
                        &dof_handler);

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  Assert(cell->get_fe().dofs_per_cell==dofs[0], ExcInternalError());

                  set_dofindices_recursively<dim,spacedim> (*tr,
                                                            p4est_coarse_cell,
                                                            cell,
                                                            cellinfo.quadrants[c],
                                                            (dofs+1));
                }
            }

          // complete all sends, so that we can
          // safely destroy the buffers.
          if (requests.size() > 0)
            MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);


#ifdef DEBUG
          {
            //check all msgs got sent and received
            unsigned int sum_send=0;
            unsigned int sum_recv=0;
            unsigned int sent=needs_to_get_cells.size();
            unsigned int recv=senders.size();

            MPI_Allreduce(&sent, &sum_send, 1, MPI_UNSIGNED, MPI_SUM, tr->get_communicator());
            MPI_Allreduce(&recv, &sum_recv, 1, MPI_UNSIGNED, MPI_SUM, tr->get_communicator());
            Assert(sum_send==sum_recv, ExcInternalError());
          }
#endif

          //update dofindices
          {
            typename DoFHandler<dim,spacedim>::active_cell_iterator
            cell, endc = dof_handler.end();

            for (cell = dof_handler.begin_active(); cell != endc; ++cell)
              if (!cell->is_artificial())
                cell->update_cell_dof_indices_cache();
          }

          // important, so that sends between two
          // calls to this function are not mixed
          // up.
          //
          // this is necessary because above we
          // just see if there are messages and
          // then receive them, without
          // discriminating where they come from
          // and whether they were sent in phase
          // 1 or 2. the need for a global
          // communication step like this barrier
          // could be avoided by receiving
          // messages specifically from those
          // processors from which we expect
          // messages, and by using different
          // tags for phase 1 and 2
          MPI_Barrier(tr->get_communicator());
#endif
        }







      }

#endif // DEAL_II_WITH_P4EST



      template <int dim, int spacedim>
      void
      ParallelDistributed<dim, spacedim>::
      distribute_dofs (DoFHandler<dim,spacedim> &dof_handler,
                       NumberCache &number_cache_current) const
      {
        NumberCache number_cache;

#ifndef DEAL_II_WITH_P4EST
        (void)dof_handler;
        Assert (false, ExcNotImplemented());
#else

        parallel::distributed::Triangulation< dim, spacedim > *tr
          = (dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (const_cast<dealii::Triangulation< dim, spacedim >*>
              (&dof_handler.get_triangulation())));
        Assert (tr != 0, ExcInternalError());

        const unsigned int
        n_cpus = Utilities::MPI::n_mpi_processes (tr->get_communicator());

        //* 1. distribute on own
        //* subdomain
        const dealii::types::global_dof_index n_initial_local_dofs =
          Implementation::distribute_dofs (0, tr->locally_owned_subdomain(),
                                           dof_handler);

        //* 2. iterate over ghostcells and
        //kill dofs that are not owned
        //by us
        std::vector<dealii::types::global_dof_index> renumbering(n_initial_local_dofs);
        for (unsigned int i=0; i<renumbering.size(); ++i)
          renumbering[i] = i;

        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          typename DoFHandler<dim,spacedim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();

          for (; cell != endc; ++cell)
            if (cell->is_ghost() &&
                (cell->subdomain_id() < tr->locally_owned_subdomain()))
              {
                // we found a
                // neighboring ghost
                // cell whose subdomain
                // is "stronger" than
                // our own subdomain

                // delete all dofs that
                // live there and that
                // we have previously
                // assigned a number to
                // (i.e. the ones on
                // the interface)
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                  if (local_dof_indices[i] != DoFHandler<dim,spacedim>::invalid_dof_index)
                    renumbering[local_dof_indices[i]]
                      = DoFHandler<dim,spacedim>::invalid_dof_index;
              }
        }


        // make indices consecutive
        number_cache.n_locally_owned_dofs = 0;
        for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
             it!=renumbering.end(); ++it)
          if (*it != DoFHandler<dim,spacedim>::invalid_dof_index)
            *it = number_cache.n_locally_owned_dofs++;

        //* 3. communicate local dofcount and
        //shift ids to make them unique
        number_cache.n_locally_owned_dofs_per_processor.resize(n_cpus);

        MPI_Allgather ( &number_cache.n_locally_owned_dofs,
                        1, DEAL_II_DOF_INDEX_MPI_TYPE,
                        &number_cache.n_locally_owned_dofs_per_processor[0],
                        1, DEAL_II_DOF_INDEX_MPI_TYPE,
                        tr->get_communicator());

        const dealii::types::global_dof_index
        shift = std::accumulate (number_cache
                                 .n_locally_owned_dofs_per_processor.begin(),
                                 number_cache
                                 .n_locally_owned_dofs_per_processor.begin()
                                 + tr->locally_owned_subdomain(),
                                 static_cast<dealii::types::global_dof_index>(0));
        for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
             it!=renumbering.end(); ++it)
          if (*it != DoFHandler<dim,spacedim>::invalid_dof_index)
            (*it) += shift;

        // now re-enumerate all dofs to
        // this shifted and condensed
        // numbering form.  we renumber
        // some dofs as invalid, so
        // choose the nocheck-version.
        Implementation::renumber_dofs (renumbering, IndexSet(0),
                                       dof_handler, false);

        // now a little bit of
        // housekeeping
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

        //* 4. send dofids of cells that are
        //ghostcells on other machines

        std::vector<bool> user_flags;
        tr->save_user_flags(user_flags);
        tr->clear_user_flags ();

        //mark all own cells for transfer
        for (typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active();
             cell != dof_handler.end(); ++cell)
          if (!cell->is_artificial())
            cell->set_user_flag();

        // add each ghostcells'
        // subdomain to the vertex and
        // keep track of interesting
        // neighbors
        std::map<unsigned int, std::set<dealii::types::subdomain_id> >
        vertices_with_ghost_neighbors;

        tr->fill_vertices_with_ghost_neighbors (vertices_with_ghost_neighbors);


        /* Send and receive cells. After this,
           only the local cells are marked,
           that received new data. This has to
           be communicated in a second
           communication step. */
        communicate_dof_indices_on_marked_cells (dof_handler,
                                                 vertices_with_ghost_neighbors,
                                                 tr->coarse_cell_to_p4est_tree_permutation,
                                                 tr->p4est_tree_to_coarse_cell_permutation);

        communicate_dof_indices_on_marked_cells (dof_handler,
                                                 vertices_with_ghost_neighbors,
                                                 tr->coarse_cell_to_p4est_tree_permutation,
                                                 tr->p4est_tree_to_coarse_cell_permutation);

        tr->load_user_flags(user_flags);

#ifdef DEBUG
        //check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          for (typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active();
               cell != dof_handler.end(); ++cell)
            if (!cell->is_artificial())
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                if (local_dof_indices.end() !=
                    std::find (local_dof_indices.begin(),
                               local_dof_indices.end(),
                               DoFHandler<dim,spacedim>::invalid_dof_index))
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

        number_cache_current = number_cache;
      }



      template <int dim, int spacedim>
      void
      ParallelDistributed<dim, spacedim>::
      distribute_mg_dofs (DoFHandler<dim,spacedim> &dof_handler,
                          std::vector<NumberCache> &number_caches) const
      {
#ifndef DEAL_II_WITH_P4EST
        (void)dof_handler;
        (void)number_caches;
        Assert (false, ExcNotImplemented());
#else

        parallel::distributed::Triangulation< dim, spacedim > *tr
          = (dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (const_cast<dealii::Triangulation< dim, spacedim >*>
              (&dof_handler.get_triangulation())));
        Assert (tr != 0, ExcInternalError());

        AssertThrow(
          (tr->settings &  parallel::distributed::Triangulation< dim, spacedim >::construct_multigrid_hierarchy),
          ExcMessage("Multigrid DoFs can only be distributed on a parallel Triangulation if the flag construct_multigrid_hierarchy is set in the constructor."));


        const unsigned int
        n_cpus = Utilities::MPI::n_mpi_processes (tr->get_communicator());

        unsigned int n_levels = Utilities::MPI::max(dof_handler.get_triangulation().n_levels(), tr->get_communicator());

        for (unsigned int level = 0; level < n_levels; ++level)
          {
            NumberCache &number_cache = number_caches[level];

            //* 1. distribute on own
            //* subdomain
            const unsigned int n_initial_local_dofs =
              Implementation::distribute_dofs_on_level(0, tr->locally_owned_subdomain(), dof_handler, level);

            //* 2. iterate over ghostcells and
            //kill dofs that are not owned
            //by us
            std::vector<dealii::types::global_dof_index> renumbering(n_initial_local_dofs);
            for (dealii::types::global_dof_index i=0; i<renumbering.size(); ++i)
              renumbering[i] = i;

            if (level<tr->n_levels())
              {
                std::vector<dealii::types::global_dof_index> local_dof_indices;

                typename DoFHandler<dim,spacedim>::level_cell_iterator
                cell = dof_handler.begin(level),
                endc = dof_handler.end(level);

                for (; cell != endc; ++cell)
                  if (cell->level_subdomain_id()!=numbers::artificial_subdomain_id &&
                      (cell->level_subdomain_id() < tr->locally_owned_subdomain()))
                    {
                      // we found a
                      // neighboring ghost
                      // cell whose subdomain
                      // is "stronger" than
                      // our own subdomain

                      // delete all dofs that
                      // live there and that
                      // we have previously
                      // assigned a number to
                      // (i.e. the ones on
                      // the interface)
                      local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                      cell->get_mg_dof_indices (local_dof_indices);
                      for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                        if (local_dof_indices[i] != DoFHandler<dim,spacedim>::invalid_dof_index)
                          renumbering[local_dof_indices[i]]
                            = DoFHandler<dim,spacedim>::invalid_dof_index;
                    }
              }


            // make indices consecutive
            number_cache.n_locally_owned_dofs = 0;
            for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
                 it!=renumbering.end(); ++it)
              if (*it != DoFHandler<dim,spacedim>::invalid_dof_index)
                *it = number_cache.n_locally_owned_dofs++;

            //* 3. communicate local dofcount and
            //shift ids to make them unique
            number_cache.n_locally_owned_dofs_per_processor.resize(n_cpus);

            MPI_Allgather ( &number_cache.n_locally_owned_dofs,
                            1, DEAL_II_DOF_INDEX_MPI_TYPE,
                            &number_cache.n_locally_owned_dofs_per_processor[0],
                            1, DEAL_II_DOF_INDEX_MPI_TYPE,
                            tr->get_communicator());

            const dealii::types::global_dof_index
            shift = std::accumulate (number_cache
                                     .n_locally_owned_dofs_per_processor.begin(),
                                     number_cache
                                     .n_locally_owned_dofs_per_processor.begin()
                                     + tr->locally_owned_subdomain(),
                                     static_cast<dealii::types::global_dof_index>(0));
            for (std::vector<dealii::types::global_dof_index>::iterator it=renumbering.begin();
                 it!=renumbering.end(); ++it)
              if (*it != DoFHandler<dim,spacedim>::invalid_dof_index)
                (*it) += shift;

            // now re-enumerate all dofs to
            // this shifted and condensed
            // numbering form.  we renumber
            // some dofs as invalid, so
            // choose the nocheck-version.
            Implementation::renumber_mg_dofs (renumbering, IndexSet(0),
                                              dof_handler, level, false);

            // now a little bit of
            // housekeeping
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
            cell, endc = dof_handler.end();
            for (cell = dof_handler.begin(); cell != endc; ++cell)
              if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id
                  && !cell->is_locally_owned_on_level())
                cell->set_user_flag();
          }

          // Phase 1. Request all marked cells from corresponding owners. If we
          // managed to get every DoF, remove the user_flag, otherwise we
          // will request them again in the step below.
          communicate_mg_ghost_cells(*tr,
                                     dof_handler,
                                     tr->coarse_cell_to_p4est_tree_permutation,
                                     tr->p4est_tree_to_coarse_cell_permutation);

          // This barrier is crucial so that messages between phase 1&2 don't
          // mix.
          MPI_Barrier(tr->get_communicator());

          // Phase 2, only request the cells that were not completed in Phase
          // 1.
          communicate_mg_ghost_cells(*tr,
                                     dof_handler,
                                     tr->coarse_cell_to_p4est_tree_permutation,
                                     tr->p4est_tree_to_coarse_cell_permutation);

#ifdef DEBUG
          // make sure we have removed all flags:
          {
            typename DoFHandler<dim,spacedim>::level_cell_iterator
            cell, endc = dof_handler.end();
            for (cell = dof_handler.begin(); cell != endc; ++cell)
              if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id
                  && !cell->is_locally_owned_on_level())
                Assert(cell->user_flag_set()==false, ExcInternalError());
          }
#endif

          tr->load_user_flags(user_flags);
        }




#ifdef DEBUG
        //check that we are really done
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;
          typename DoFHandler<dim,spacedim>::level_cell_iterator
          cell, endc = dof_handler.end();

          for (cell = dof_handler.begin(); cell != endc; ++cell)
            if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id)
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_mg_dof_indices (local_dof_indices);
                if (local_dof_indices.end() !=
                    std::find (local_dof_indices.begin(),
                               local_dof_indices.end(),
                               DoFHandler<dim,spacedim>::invalid_dof_index))
                  {
                    Assert(false, ExcMessage ("not all DoFs got distributed!"));
                  }
              }
        }
#endif // DEBUG

#endif // DEAL_II_WITH_P4EST
      }


      template <int dim, int spacedim>
      void
      ParallelDistributed<dim, spacedim>::
      renumber_dofs (const std::vector<dealii::types::global_dof_index> &new_numbers,
                     dealii::DoFHandler<dim,spacedim> &dof_handler,
                     NumberCache &number_cache_current) const
      {
        (void)new_numbers;
        (void)dof_handler;

        Assert (new_numbers.size() == dof_handler.locally_owned_dofs().n_elements(),
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

        number_cache.locally_owned_dofs = IndexSet (dof_handler.n_dofs());
        if (dof_handler.locally_owned_dofs().n_elements()>0)
          {
            std::vector<dealii::types::global_dof_index> new_numbers_sorted (new_numbers);
            std::sort(new_numbers_sorted.begin(), new_numbers_sorted.end());
            std::vector<dealii::types::global_dof_index>::const_iterator it = new_numbers_sorted.begin();
            const unsigned int n_blocks = dof_handler.get_fe().n_blocks();
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
                dof_handler.n_locally_owned_dofs(),
                ExcInternalError());

        // then also set this number in our own copy
        number_cache.n_locally_owned_dofs = dof_handler.n_locally_owned_dofs();

        // mark not locally active DoFs as invalid
        {
          std::vector<dealii::types::global_dof_index> local_dof_indices;

          typename DoFHandler<dim,spacedim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();

          for (; cell != endc; ++cell)
            if (!cell->is_artificial())
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                  {
                    if (local_dof_indices[i] == DoFHandler<dim,spacedim>::invalid_dof_index)
                      continue;

                    if (!dof_handler.locally_owned_dofs().is_element(local_dof_indices[i]))
                      {
                        //this DoF is not owned by us, so set it to invalid.
                        local_dof_indices[i]
                          = DoFHandler<dim,spacedim>::invalid_dof_index;
                      }
                  }

                cell->set_dof_indices (local_dof_indices);
              }
        }


        // renumber. Skip when there is nothing to do because we own no DoF.
        if (dof_handler.locally_owned_dofs().n_elements() > 0)
          Implementation::renumber_dofs (new_numbers,
                                         dof_handler.locally_owned_dofs(),
                                         dof_handler,
                                         false);

        // communication
        {
          parallel::distributed::Triangulation< dim, spacedim > *tr
            = (dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
               (const_cast<dealii::Triangulation< dim, spacedim >*>
                (&dof_handler.get_triangulation())));
          Assert (tr != 0, ExcInternalError());

          std::vector<bool> user_flags;
          tr->save_user_flags(user_flags);
          tr->clear_user_flags ();

          //mark all own cells for transfer
          typename DoFHandler<dim,spacedim>::active_cell_iterator
          cell, endc = dof_handler.end();
          for (cell = dof_handler.begin_active(); cell != endc; ++cell)
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
          communicate_dof_indices_on_marked_cells (dof_handler,
                                                   vertices_with_ghost_neighbors,
                                                   tr->coarse_cell_to_p4est_tree_permutation,
                                                   tr->p4est_tree_to_coarse_cell_permutation);

          communicate_dof_indices_on_marked_cells (dof_handler,
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
          MPI_Allgather(&my_data[0], max_size, MPI_BYTE,
                        &buffer[0], max_size, MPI_BYTE,
                        tr->get_communicator());

          number_cache.locally_owned_dofs_per_processor.resize (n_cpus);
          number_cache.n_locally_owned_dofs_per_processor.resize (n_cpus);
          for (unsigned int i=0; i<n_cpus; ++i)
            {
              std::stringstream strstr;
              strstr.write(&buffer[i*max_size],max_size);
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

        number_cache_current = number_cache;
      }
    }
  }
}




/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler_policy.inst"


DEAL_II_NAMESPACE_CLOSE
