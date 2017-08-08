// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/dof_level.h>
#include <deal.II/hp/dof_faces.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler_policy.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_levels.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <set>
#include <algorithm>
#include <functional>

DEAL_II_NAMESPACE_OPEN

// The following is necessary for compilation under Visual Studio which is unable to correctly
// distinguish between dealii::DoFHandler and dealii::hp::DoFHandler.
// Plus it makes code in dof_handler.cc easier to read.
// Requires C++11 support which is in Visual Studio 2013 and newer.
#if defined(_MSC_VER) && (_MSC_VER >= 1800)
template <int dim, int spacedim> using HpDoFHandler = ::dealii::hp::DoFHandler<dim, spacedim>;
#else
// When using older Visual Studio or a different compiler just fall back.
#define HpDoFHandler DoFHandler
#endif

namespace parallel
{
  namespace distributed
  {
    template <int, int> class Triangulation;
  }
}




namespace internal
{
  namespace hp
  {
    namespace DoFHandler
    {
      // access class dealii::hp::DoFHandler instead of namespace
      // internal::hp::DoFHandler, etc
      using dealii::hp::DoFHandler;

      /**
       * A class with the same purpose as the similarly named class of the
       * Triangulation class. See there for more information.
       */
      struct Implementation
      {
        /**
         * Do that part of reserving space that pertains to vertices,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static
        void
        reserve_space_vertices (DoFHandler<dim,spacedim> &dof_handler)
        {
          // The final step is allocating memory is to set up vertex
          // dof information. since vertices are sequentially
          // numbered, what we do first is to set up an array in which
          // we record whether a vertex is associated with any of the
          // given fe's, by setting a bit. in a later step, we then
          // actually allocate memory for the required dofs
          std::vector<std::vector<bool> >
          vertex_fe_association (dof_handler.finite_elements->size(),
                                 std::vector<bool> (dof_handler.tria->n_vertices(), false));

          for (typename HpDoFHandler<dim,spacedim>::active_cell_iterator
               cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
              vertex_fe_association[cell->active_fe_index()][cell->vertex_index(v)]
                = true;

          // in debug mode, make sure that each vertex is associated
          // with at least one fe (note that except for unused
          // vertices, all vertices are actually active)
#ifdef DEBUG
          for (unsigned int v=0; v<dof_handler.tria->n_vertices(); ++v)
            if (dof_handler.tria->vertex_used(v) == true)
              {
                unsigned int fe=0;
                for (; fe<dof_handler.finite_elements->size(); ++fe)
                  if (vertex_fe_association[fe][v] == true)
                    break;
                Assert (fe != dof_handler.finite_elements->size(), ExcInternalError());
              }
#endif

          // next count how much memory we actually need. for each
          // vertex, we need one slot per fe to store the fe_index,
          // plus dofs_per_vertex for this fe. in addition, we need
          // one slot as the end marker for the fe_indices. at the
          // same time already fill the vertex_dofs_offsets field
          dof_handler.vertex_dofs_offsets.resize (dof_handler.tria->n_vertices(),
                                                  numbers::invalid_dof_index);

          unsigned int vertex_slots_needed = 0;
          for (unsigned int v=0; v<dof_handler.tria->n_vertices(); ++v)
            if (dof_handler.tria->vertex_used(v) == true)
              {
                dof_handler.vertex_dofs_offsets[v] = vertex_slots_needed;

                for (unsigned int fe=0; fe<dof_handler.finite_elements->size(); ++fe)
                  if (vertex_fe_association[fe][v] == true)
                    vertex_slots_needed += (*dof_handler.finite_elements)[fe].dofs_per_vertex + 1;
                ++vertex_slots_needed;
              }

          // now allocate the space we have determined we need, and
          // set up the linked lists for each of the vertices
          dof_handler.vertex_dofs.resize (vertex_slots_needed,
                                          numbers::invalid_dof_index);
          for (unsigned int v=0; v<dof_handler.tria->n_vertices(); ++v)
            if (dof_handler.tria->vertex_used(v) == true)
              {
                types::global_dof_index pointer = dof_handler.vertex_dofs_offsets[v];
                for (unsigned int fe=0; fe<dof_handler.finite_elements->size(); ++fe)
                  if (vertex_fe_association[fe][v] == true)
                    {
                      // if this vertex uses this fe, then set the
                      // fe_index and move the pointer ahead
                      dof_handler.vertex_dofs[pointer] = fe;
                      pointer += (*dof_handler.finite_elements)[fe].dofs_per_vertex + 1;
                    }
                // finally place the end marker
                dof_handler.vertex_dofs[pointer] = numbers::invalid_dof_index;
              }
        }



        /**
         * Reserve enough space in the <tt>levels[]</tt> objects to
         * store the numbers of the degrees of freedom needed for the
         * given element. The given element is that one which was
         * selected when calling @p distribute_dofs the last time.
         */
        template <int spacedim>
        static
        void
        reserve_space (DoFHandler<1,spacedim> &dof_handler)
        {
          const unsigned int dim = 1;

          typedef DoFHandler<dim,spacedim> BaseClass;

          Assert (dof_handler.finite_elements != nullptr,
                  typename BaseClass::ExcNoFESelected());
          Assert (dof_handler.finite_elements->size() > 0,
                  typename BaseClass::ExcNoFESelected());
          Assert (dof_handler.tria->n_levels() > 0,
                  typename
                  BaseClass::ExcInvalidTriangulation());
          Assert (dof_handler.tria->n_levels() == dof_handler.levels.size (),
                  ExcInternalError ());

          // Release all space except the active_fe_indices field
          // which we have to backup before
          {
            std::vector<std::vector<DoFLevel::active_fe_index_type> >
            active_fe_backup(dof_handler.levels.size ());
            for (unsigned int level = 0; level<dof_handler.levels.size (); ++level)
              std::swap (dof_handler.levels[level]->active_fe_indices,
                         active_fe_backup[level]);

            // delete all levels and set them up newly, since vectors
            // are troublesome if you want to change their size
            dof_handler.clear_space ();

            for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
              {
                dof_handler.levels.emplace_back (new internal::hp::DoFLevel);
                std::swap (active_fe_backup[level],
                           dof_handler.levels[level]->active_fe_indices);
              }
          }

          // LINE (CELL) DOFs

          // count how much space we need on each level for the cell
          // dofs and set the dof_*_offsets data. initially set the
          // latter to an invalid index, and only later set it to
          // something reasonable for active dof_handler.cells
          //
          // note that for dof_handler.cells, the situation is simpler
          // than for other (lower dimensional) objects since exactly
          // one finite element is used for it
          for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
            {
              dof_handler.levels[level]->dof_offsets
                = std::vector<DoFLevel::offset_type> (
                    dof_handler.tria->n_raw_lines(level),
                    (DoFLevel::offset_type)(-1));
              dof_handler.levels[level]->cell_cache_offsets
                = std::vector<DoFLevel::offset_type> (
                    dof_handler.tria->n_raw_lines(level),
                    (DoFLevel::offset_type)(-1));

              types::global_dof_index next_free_dof = 0;
              types::global_dof_index cache_size = 0;
              typename HpDoFHandler<dim,spacedim>::active_cell_iterator
              cell=dof_handler.begin_active(level),
              endc=dof_handler.end_active(level);
              for (; cell!=endc; ++cell)
                if (!cell->has_children())
                  {
                    dof_handler.levels[level]->dof_offsets[cell->index()] = next_free_dof;
                    next_free_dof += cell->get_fe().dofs_per_line;

                    dof_handler.levels[level]->cell_cache_offsets[cell->index()] = cache_size;
                    cache_size += cell->get_fe().dofs_per_cell;
                  }

              dof_handler.levels[level]->dof_indices
                = std::vector<types::global_dof_index> (next_free_dof,
                                                        numbers::invalid_dof_index);
              dof_handler.levels[level]->cell_dof_indices_cache
                = std::vector<types::global_dof_index> (cache_size,
                                                        numbers::invalid_dof_index);
            }

          // safety check: make sure that the number of DoFs we
          // allocated is actually correct (above we have also set the
          // dof_*_offsets field, so we couldn't use this simpler
          // algorithm)
#ifdef DEBUG
          for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
            {
              types::global_dof_index counter = 0;
              typename HpDoFHandler<dim,spacedim>::active_cell_iterator
              cell=dof_handler.begin_active(level),
              endc=dof_handler.end_active(level);
              for (; cell!=endc; ++cell)
                if (!cell->has_children())
                  counter += cell->get_fe().dofs_per_line;

              Assert (dof_handler.levels[level]->dof_indices.size() == counter,
                      ExcInternalError());
              Assert (static_cast<unsigned int>
                      (std::count (dof_handler.levels[level]->dof_offsets.begin(),
                                   dof_handler.levels[level]->dof_offsets.end(),
                                   (DoFLevel::offset_type)(-1)))
                      ==
                      dof_handler.tria->n_raw_lines(level) - dof_handler.tria->n_active_lines(level),
                      ExcInternalError());
            }
#endif


          // VERTEX DOFS
          reserve_space_vertices (dof_handler);
        }


        template <int spacedim>
        static
        void
        reserve_space (DoFHandler<2,spacedim> &dof_handler)
        {
          const unsigned int dim = 2;

          typedef DoFHandler<dim,spacedim> BaseClass;

          Assert (dof_handler.finite_elements != nullptr,
                  typename BaseClass::ExcNoFESelected());
          Assert (dof_handler.finite_elements->size() > 0,
                  typename BaseClass::ExcNoFESelected());
          Assert (dof_handler.tria->n_levels() > 0,
                  typename BaseClass::ExcInvalidTriangulation());
          Assert (dof_handler.tria->n_levels() == dof_handler.levels.size (),
                  ExcInternalError ());

          // Release all space except the active_fe_indices field
          // which we have to backup before
          {
            std::vector<std::vector<DoFLevel::active_fe_index_type> >
            active_fe_backup(dof_handler.levels.size ());
            for (unsigned int level = 0; level<dof_handler.levels.size (); ++level)
              std::swap (dof_handler.levels[level]->active_fe_indices,
                         active_fe_backup[level]);

            // delete all levels and set them up newly, since vectors
            // are troublesome if you want to change their size
            dof_handler.clear_space ();

            for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
              {
                dof_handler.levels.emplace_back (new internal::hp::DoFLevel);
                std::swap (active_fe_backup[level],
                           dof_handler.levels[level]->active_fe_indices);
              }
            dof_handler.faces.reset (new internal::hp::DoFIndicesOnFaces<2>);
          }

          // QUAD (CELL) DOFs

          // count how much space we need on each level for the cell
          // dofs and set the dof_*_offsets data. initially set the
          // latter to an invalid index, and only later set it to
          // something reasonable for active dof_handler.cells
          //
          // note that for dof_handler.cells, the situation is simpler
          // than for other (lower dimensional) objects since exactly
          // one finite element is used for it
          for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
            {
              dof_handler.levels[level]->dof_offsets
                = std::vector<DoFLevel::offset_type> (
                    dof_handler.tria->n_raw_quads(level),
                    (DoFLevel::offset_type)(-1));
              dof_handler.levels[level]->cell_cache_offsets
                = std::vector<DoFLevel::offset_type> (
                    dof_handler.tria->n_raw_quads(level),
                    (DoFLevel::offset_type)(-1));

              types::global_dof_index next_free_dof = 0;
              types::global_dof_index cache_size = 0;
              typename HpDoFHandler<dim,spacedim>::active_cell_iterator
              cell=dof_handler.begin_active(level),
              endc=dof_handler.end_active(level);
              for (; cell!=endc; ++cell)
                if (!cell->has_children())
                  {
                    dof_handler.levels[level]->dof_offsets[cell->index()] = next_free_dof;
                    next_free_dof += cell->get_fe().dofs_per_quad;

                    dof_handler.levels[level]->cell_cache_offsets[cell->index()] = cache_size;
                    cache_size += cell->get_fe().dofs_per_cell;
                  }

              dof_handler.levels[level]->dof_indices
                = std::vector<types::global_dof_index> (next_free_dof,
                                                        numbers::invalid_dof_index);
              dof_handler.levels[level]->cell_dof_indices_cache
                = std::vector<types::global_dof_index> (cache_size,
                                                        numbers::invalid_dof_index);
            }

          // safety check: make sure that the number of DoFs we
          // allocated is actually correct (above we have also set the
          // dof_*_offsets field, so we couldn't use this simpler
          // algorithm)
#ifdef DEBUG
          for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
            {
              types::global_dof_index counter = 0;
              typename HpDoFHandler<dim,spacedim>::active_cell_iterator
              cell=dof_handler.begin_active(level),
              endc=dof_handler.end_active(level);
              for (; cell!=endc; ++cell)
                if (!cell->has_children())
                  counter += cell->get_fe().dofs_per_quad;

              Assert (dof_handler.levels[level]->dof_indices.size() == counter,
                      ExcInternalError());
              Assert (static_cast<unsigned int>
                      (std::count (dof_handler.levels[level]->dof_offsets.begin(),
                                   dof_handler.levels[level]->dof_offsets.end(),
                                   (DoFLevel::offset_type)(-1)))
                      ==
                      dof_handler.tria->n_raw_quads(level) - dof_handler.tria->n_active_quads(level),
                      ExcInternalError());
            }
#endif


          // LINE DOFS
          //
          // same here: count line dofs, then allocate as much space
          // as we need and prime the linked list for lines (see the
          // description in hp::DoFLevel) with the indices we will
          // need. note that our task is more complicated since two
          // adjacent dof_handler.cells may have different
          // active_fe_indices, in which case we need to allocate
          // *two* sets of line dofs for the same line
          //
          // the way we do things is that we loop over all active
          // dof_handler.cells (these are the ones that have DoFs only
          // anyway) and all their dof_handler.faces. We note in the
          // user flags whether we have previously visited a face and
          // if so skip it (consequently, we have to save and later
          // restore the line flags)
          {
            std::vector<bool> saved_line_user_flags;
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .save_user_flags_line (saved_line_user_flags);
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .clear_user_flags_line ();

            // an array to hold how many slots (see the hp::DoFLevel
            // class) we will have to store on each level
            unsigned int n_line_slots = 0;

            for (typename HpDoFHandler<dim,spacedim>::active_cell_iterator
                 cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                if (! cell->face(face)->user_flag_set())
                  {
                    // ok, face has not been visited. so we need to
                    // allocate space for it. let's see how much we
                    // need: we need one set if a) there is no
                    // neighbor behind this face, or b) the neighbor
                    // is either coarser or finer than we are, or c)
                    // the neighbor is neither coarser nor finer, but
                    // has happens to have the same active_fe_index:
                    if (cell->at_boundary(face)
                        ||
                        cell->face(face)->has_children()
                        ||
                        cell->neighbor_is_coarser(face)
                        ||
                        (!cell->at_boundary(face)
                         &&
                         (cell->active_fe_index() == cell->neighbor(face)->active_fe_index())))
                      // ok, one set of dofs. that makes one index, 1
                      // times dofs_per_line dofs, and one stop index
                      n_line_slots
                      += (*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_line + 2;

                    // otherwise we do indeed need two sets, i.e. two
                    // indices, two sets of dofs, and one stop index:
                    else
                      n_line_slots
                      += ((*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_line
                          +
                          (*dof_handler.finite_elements)[cell->neighbor(face)->active_fe_index()]
                          .dofs_per_line
                          +
                          3);

                    // mark this face as visited
                    cell->face(face)->set_user_flag ();
                  }

            // now that we know how many line dofs we will have to
            // have on each level, allocate the memory. note that we
            // allocate offsets for all lines, though only the active
            // ones will have a non-invalid value later on
            dof_handler.faces->lines.dof_offsets
              = std::vector<unsigned int> (dof_handler.tria->n_raw_lines(),
                                           (unsigned int)(-1));
            dof_handler.faces->lines.dofs
              = std::vector<types::global_dof_index> (n_line_slots,
                                                      numbers::invalid_dof_index);

            // with the memory now allocated, loop over the
            // dof_handler.cells again and prime the _offset values as
            // well as the fe_index fields
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .clear_user_flags_line ();

            unsigned int next_free_line_slot = 0;

            for (typename HpDoFHandler<dim,spacedim>::active_cell_iterator
                 cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                if (! cell->face(face)->user_flag_set())
                  {
                    // same decision tree as before
                    if (cell->at_boundary(face)
                        ||
                        cell->face(face)->has_children()
                        ||
                        cell->neighbor_is_coarser(face)
                        ||
                        (!cell->at_boundary(face)
                         &&
                         (cell->active_fe_index() == cell->neighbor(face)->active_fe_index())))
                      {
                        dof_handler.faces
                        ->lines.dof_offsets[cell->face(face)->index()]
                          = next_free_line_slot;

                        // set first slot for this line to
                        // active_fe_index of this face
                        dof_handler.faces
                        ->lines.dofs[next_free_line_slot]
                          = cell->active_fe_index();

                        // the next dofs_per_line indices remain unset
                        // for the moment (i.e. at invalid_dof_index).
                        // following this comes the stop index, which
                        // also is invalid_dof_index and therefore
                        // does not have to be explicitly set

                        // finally, mark those slots as used
                        next_free_line_slot
                        += (*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_line + 2;
                      }
                    else
                      {
                        dof_handler.faces
                        ->lines.dof_offsets[cell->face(face)->index()]
                          = next_free_line_slot;

                        // set first slot for this line to
                        // active_fe_index of this face
                        dof_handler.faces
                        ->lines.dofs[next_free_line_slot]
                          = cell->active_fe_index();

                        // the next dofs_per_line indices remain unset
                        // for the moment (i.e. at invalid_dof_index).
                        //
                        // then comes the fe_index for the neighboring
                        // cell:
                        dof_handler.faces
                        ->lines.dofs[next_free_line_slot
                                     +
                                     (*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_line
                                     +
                                     1]
                          = cell->neighbor(face)->active_fe_index();
                        // then again a set of dofs that we need not
                        // set right now
                        //
                        // following this comes the stop index, which
                        // also is invalid_dof_index and therefore
                        // does not have to be explicitly set

                        // finally, mark those slots as used
                        next_free_line_slot
                        += ((*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_line
                            +
                            (*dof_handler.finite_elements)[cell->neighbor(face)->active_fe_index()]
                            .dofs_per_line
                            +
                            3);
                      }

                    // mark this face as visited
                    cell->face(face)->set_user_flag ();
                  }

            // we should have moved the cursor for each level to the
            // total number of dofs on that level. check that
            Assert (next_free_line_slot == n_line_slots,
                    ExcInternalError());

            // at the end, restore the user flags for the lines
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .load_user_flags_line (saved_line_user_flags);
          }


          // VERTEX DOFS
          reserve_space_vertices (dof_handler);
        }


        template <int spacedim>
        static
        void
        reserve_space (DoFHandler<3,spacedim> &dof_handler)
        {
          const unsigned int dim = 3;

          typedef DoFHandler<dim,spacedim> BaseClass;

          Assert (dof_handler.finite_elements != nullptr,
                  typename BaseClass::ExcNoFESelected());
          Assert (dof_handler.finite_elements->size() > 0,
                  typename BaseClass::ExcNoFESelected());
          Assert (dof_handler.tria->n_levels() > 0,
                  typename BaseClass::ExcInvalidTriangulation());
          Assert (dof_handler.tria->n_levels() == dof_handler.levels.size (),
                  ExcInternalError ());

          // Release all space except the active_fe_indices field
          // which we have to backup before
          {
            std::vector<std::vector<DoFLevel::active_fe_index_type> >
            active_fe_backup(dof_handler.levels.size ());
            for (unsigned int level = 0; level<dof_handler.levels.size (); ++level)
              std::swap (dof_handler.levels[level]->active_fe_indices,
                         active_fe_backup[level]);

            // delete all levels and set them up newly, since vectors
            // are troublesome if you want to change their size
            dof_handler.clear_space ();

            for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
              {
                dof_handler.levels.emplace_back (new internal::hp::DoFLevel);
                std::swap (active_fe_backup[level],
                           dof_handler.levels[level]->active_fe_indices);
              }
            dof_handler.faces.reset (new internal::hp::DoFIndicesOnFaces<3>);
          }

          // HEX (CELL) DOFs

          // count how much space we need on each level for the cell
          // dofs and set the dof_*_offsets data. initially set the
          // latter to an invalid index, and only later set it to
          // something reasonable for active dof_handler.cells
          //
          // note that for dof_handler.cells, the situation is simpler
          // than for other (lower dimensional) objects since exactly
          // one finite element is used for it
          for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
            {
              dof_handler.levels[level]->dof_offsets
                = std::vector<DoFLevel::offset_type> (
                    dof_handler.tria->n_raw_hexs(level),
                    (DoFLevel::offset_type)(-1));
              dof_handler.levels[level]->cell_cache_offsets
                = std::vector<DoFLevel::offset_type> (
                    dof_handler.tria->n_raw_hexs(level),
                    (DoFLevel::offset_type)(-1));

              types::global_dof_index next_free_dof = 0;
              types::global_dof_index cache_size = 0;
              typename HpDoFHandler<dim,spacedim>::active_cell_iterator
              cell=dof_handler.begin_active(level),
              endc=dof_handler.end_active(level);
              for (; cell!=endc; ++cell)
                if (!cell->has_children())
                  {
                    dof_handler.levels[level]->dof_offsets[cell->index()] = next_free_dof;
                    next_free_dof += cell->get_fe().dofs_per_hex;

                    dof_handler.levels[level]->cell_cache_offsets[cell->index()] = cache_size;
                    cache_size += cell->get_fe().dofs_per_cell;
                  }

              dof_handler.levels[level]->dof_indices
                = std::vector<types::global_dof_index> (next_free_dof,
                                                        numbers::invalid_dof_index);
              dof_handler.levels[level]->cell_dof_indices_cache
                = std::vector<types::global_dof_index> (cache_size,
                                                        numbers::invalid_dof_index);
            }

          // safety check: make sure that the number of DoFs we
          // allocated is actually correct (above we have also set the
          // dof_*_offsets field, so we couldn't use this simpler
          // algorithm)
#ifdef DEBUG
          for (unsigned int level=0; level<dof_handler.tria->n_levels(); ++level)
            {
              types::global_dof_index counter = 0;
              typename HpDoFHandler<dim,spacedim>::active_cell_iterator
              cell=dof_handler.begin_active(level),
              endc=dof_handler.end_active(level);
              for (; cell!=endc; ++cell)
                if (!cell->has_children())
                  counter += cell->get_fe().dofs_per_hex;

              Assert (dof_handler.levels[level]->dof_indices.size() == counter,
                      ExcInternalError());
              Assert (static_cast<unsigned int>
                      (std::count (dof_handler.levels[level]->dof_offsets.begin(),
                                   dof_handler.levels[level]->dof_offsets.end(),
                                   (DoFLevel::offset_type)(-1)))
                      ==
                      dof_handler.tria->n_raw_hexs(level) - dof_handler.tria->n_active_hexs(level),
                      ExcInternalError());
            }
#endif


          // QUAD DOFS
          //
          // same here: count quad dofs, then allocate as much space
          // as we need and prime the linked list for quad (see the
          // description in hp::DoFLevel) with the indices we will
          // need. note that our task is more complicated since two
          // adjacent dof_handler.cells may have different
          // active_fe_indices, in which case we need to allocate
          // *two* sets of line dofs for the same line
          //
          // the way we do things is that we loop over all active
          // dof_handler.cells (these are the ones that have DoFs only
          // anyway) and all their dof_handler.faces. We note in the
          // user flags whether we have previously visited a face and
          // if so skip it (consequently, we have to save and later
          // restore the line flags)
          {
            std::vector<bool> saved_quad_user_flags;
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .save_user_flags_quad (saved_quad_user_flags);
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .clear_user_flags_quad ();

            // examine, how how many slots (see the hp::DoFLevel
            // class) we will have to store
            unsigned int n_quad_slots = 0;

            for (typename HpDoFHandler<dim,spacedim>::active_cell_iterator
                 cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                if (! cell->face(face)->user_flag_set())
                  {
                    // ok, face has not been visited. so we need to
                    // allocate space for it. let's see how much we
                    // need: we need one set if a) there is no
                    // neighbor behind this face, or b) the neighbor
                    // is not on the same level or further refined, or
                    // c) the neighbor is on the same level, but
                    // happens to have the same active_fe_index:
                    if (cell->at_boundary(face)
                        ||
                        cell->face(face)->has_children()
                        ||
                        cell->neighbor_is_coarser(face)
                        ||
                        (!cell->at_boundary(face)
                         &&
                         (cell->active_fe_index() == cell->neighbor(face)->active_fe_index())))
                      // ok, one set of dofs. that makes one index, 1
                      // times dofs_per_quad dofs, and one stop index
                      n_quad_slots
                      += (*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_quad + 2;

                    // otherwise we do indeed need two sets, i.e. two
                    // indices, two sets of dofs, and one stop index:
                    else
                      n_quad_slots
                      += ((*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_quad
                          +
                          (*dof_handler.finite_elements)[cell->neighbor(face)->active_fe_index()]
                          .dofs_per_quad
                          +
                          3);

                    // mark this face as visited
                    cell->face(face)->set_user_flag ();
                  }

            // now that we know how many quad dofs we will have to
            // have, allocate the memory. note that we allocate
            // offsets for all quads, though only the active ones will
            // have a non-invalid value later on
            if (true)
              {
                dof_handler.faces->quads.dof_offsets
                  = std::vector<unsigned int>
                    (dof_handler.tria->n_raw_quads(),
                     (unsigned int)(-1));
                dof_handler.faces->quads.dofs
                  = std::vector<types::global_dof_index> (n_quad_slots,
                                                          numbers::invalid_dof_index);
              }

            // with the memory now allocated, loop over the
            // dof_handler.cells again and prime the _offset values as
            // well as the fe_index fields
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .clear_user_flags_quad ();

            unsigned int next_free_quad_slot = 0;

            for (typename HpDoFHandler<dim,spacedim>::active_cell_iterator
                 cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                if (! cell->face(face)->user_flag_set())
                  {
                    // same decision tree as before
                    if (cell->at_boundary(face)
                        ||
                        cell->face(face)->has_children()
                        ||
                        cell->neighbor_is_coarser(face)
                        ||
                        (!cell->at_boundary(face)
                         &&
                         (cell->active_fe_index() == cell->neighbor(face)->active_fe_index())))
                      {
                        dof_handler.faces
                        ->quads.dof_offsets[cell->face(face)->index()]
                          = next_free_quad_slot;

                        // set first slot for this quad to
                        // active_fe_index of this face
                        dof_handler.faces
                        ->quads.dofs[next_free_quad_slot]
                          = cell->active_fe_index();

                        // the next dofs_per_quad indices remain unset
                        // for the moment (i.e. at invalid_dof_index).
                        // following this comes the stop index, which
                        // also is invalid_dof_index and therefore
                        // does not have to be explicitly set

                        // finally, mark those slots as used
                        next_free_quad_slot
                        += (*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_quad + 2;
                      }
                    else
                      {
                        dof_handler.faces
                        ->quads.dof_offsets[cell->face(face)->index()]
                          = next_free_quad_slot;

                        // set first slot for this quad to
                        // active_fe_index of this face
                        dof_handler.faces
                        ->quads.dofs[next_free_quad_slot]
                          = cell->active_fe_index();

                        // the next dofs_per_quad indices remain unset
                        // for the moment (i.e. at invalid_dof_index).
                        //
                        // then comes the fe_index for the neighboring
                        // cell:
                        dof_handler.faces
                        ->quads.dofs[next_free_quad_slot
                                     +
                                     (*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_quad
                                     +
                                     1]
                          = cell->neighbor(face)->active_fe_index();
                        // then again a set of dofs that we need not
                        // set right now
                        //
                        // following this comes the stop index, which
                        // also is invalid_dof_index and therefore
                        // does not have to be explicitly set

                        // finally, mark those slots as used
                        next_free_quad_slot
                        += ((*dof_handler.finite_elements)[cell->active_fe_index()].dofs_per_quad
                            +
                            (*dof_handler.finite_elements)[cell->neighbor(face)->active_fe_index()]
                            .dofs_per_quad
                            +
                            3);
                      }

                    // mark this face as visited
                    cell->face(face)->set_user_flag ();
                  }

            // we should have moved the cursor to the total number of
            // dofs. check that
            Assert (next_free_quad_slot == n_quad_slots,
                    ExcInternalError());

            // at the end, restore the user flags for the quads
            const_cast<dealii::Triangulation<dim,spacedim>&>(*dof_handler.tria)
            .load_user_flags_quad (saved_quad_user_flags);
          }


          // LINE DOFS

          // the situation here is pretty much like with vertices:
          // there can be an arbitrary number of finite elements
          // associated with each line.
          //
          // the algorithm we use is somewhat similar to what we do in
          // reserve_space_vertices()
          if (true)
            {
              // what we do first is to set up an array in which we
              // record whether a line is associated with any of the
              // given fe's, by setting a bit. in a later step, we
              // then actually allocate memory for the required dofs
              std::vector<std::vector<bool> >
              line_fe_association (dof_handler.finite_elements->size(),
                                   std::vector<bool> (dof_handler.tria->n_raw_lines(),
                                                      false));

              for (typename HpDoFHandler<dim,spacedim>::active_cell_iterator
                   cell=dof_handler.begin_active();
                   cell!=dof_handler.end(); ++cell)
                for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
                  line_fe_association[cell->active_fe_index()][cell->line_index(l)]
                    = true;

              // first check which of the lines is used at all,
              // i.e. is associated with a finite element. we do this
              // since not all lines may actually be used, in which
              // case we do not have to allocate any memory at all
              std::vector<bool> line_is_used (dof_handler.tria->n_raw_lines(), false);
              for (unsigned int line=0; line<dof_handler.tria->n_raw_lines(); ++line)
                for (unsigned int fe=0; fe<dof_handler.finite_elements->size(); ++fe)
                  if (line_fe_association[fe][line] == true)
                    {
                      line_is_used[line] = true;
                      break;
                    }

              // next count how much memory we actually need. for each
              // line, we need one slot per fe to store the fe_index,
              // plus dofs_per_line for this fe. in addition, we need
              // one slot as the end marker for the fe_indices. at the
              // same time already fill the line_dofs_offsets field
              dof_handler.faces->lines.dof_offsets
              .resize (dof_handler.tria->n_raw_lines(),
                       numbers::invalid_unsigned_int);

              unsigned int line_slots_needed = 0;
              for (unsigned int line=0; line<dof_handler.tria->n_raw_lines(); ++line)
                if (line_is_used[line] == true)
                  {
                    dof_handler.faces->lines.dof_offsets[line] = line_slots_needed;

                    for (unsigned int fe=0; fe<dof_handler.finite_elements->size(); ++fe)
                      if (line_fe_association[fe][line] == true)
                        line_slots_needed += (*dof_handler.finite_elements)[fe].dofs_per_line + 1;
                    ++line_slots_needed;
                  }

              // now allocate the space we have determined we need,
              // and set up the linked lists for each of the lines
              dof_handler.faces->lines.dofs.resize (line_slots_needed,
                                                    numbers::invalid_dof_index);
              for (unsigned int line=0; line<dof_handler.tria->n_raw_lines(); ++line)
                if (line_is_used[line] == true)
                  {
                    unsigned int pointer = dof_handler.faces->lines.dof_offsets[line];
                    for (unsigned int fe=0; fe<dof_handler.finite_elements->size(); ++fe)
                      if (line_fe_association[fe][line] == true)
                        {
                          // if this line uses this fe, then set the
                          // fe_index and move the pointer ahead
                          dof_handler.faces->lines.dofs[pointer] = fe;
                          pointer += (*dof_handler.finite_elements)[fe].dofs_per_line + 1;
                        }
                    // finally place the end marker
                    dof_handler.faces->lines.dofs[pointer] = numbers::invalid_dof_index;
                  }
            }



          // VERTEX DOFS
          reserve_space_vertices (dof_handler);
        }


        /**
         * Implement the function of same name in the mother class.
         */
        template <int spacedim>
        static
        unsigned int
        max_couplings_between_dofs (const DoFHandler<1,spacedim> &dof_handler)
        {
          return std::min(static_cast<types::global_dof_index> (3*
                                                                dof_handler.finite_elements->max_dofs_per_vertex() +
                                                                2*dof_handler.finite_elements->max_dofs_per_line()),
                          dof_handler.n_dofs());
        }



        template <int spacedim>
        static
        unsigned int
        max_couplings_between_dofs (const DoFHandler<2,spacedim> &dof_handler)
        {
          // get these numbers by drawing pictures and counting...
          // example:
          //   |     |     |
          // --x-----x--x--X--
          //   |     |  |  |
          //   |     x--x--x
          //   |     |  |  |
          // --x--x--*--x--x--
          //   |  |  |     |
          //   x--x--x     |
          //   |  |  |     |
          // --X--x--x-----x--
          //   |     |     |
          // x = vertices connected with center vertex *;
          //   = total of 19
          //
          // (the X vertices are connected with * if the vertices
          // adjacent to X are hanging nodes) count lines -> 28 (don't
          // forget to count mother and children separately!)
          types::global_dof_index max_couplings;
          switch (dof_handler.tria->max_adjacent_cells())
            {
            case 4:
              max_couplings=19*dof_handler.finite_elements->max_dofs_per_vertex() +
                            28*dof_handler.finite_elements->max_dofs_per_line() +
                            8*dof_handler.finite_elements->max_dofs_per_quad();
              break;
            case 5:
              max_couplings=21*dof_handler.finite_elements->max_dofs_per_vertex() +
                            31*dof_handler.finite_elements->max_dofs_per_line() +
                            9*dof_handler.finite_elements->max_dofs_per_quad();
              break;
            case 6:
              max_couplings=28*dof_handler.finite_elements->max_dofs_per_vertex() +
                            42*dof_handler.finite_elements->max_dofs_per_line() +
                            12*dof_handler.finite_elements->max_dofs_per_quad();
              break;
            case 7:
              max_couplings=30*dof_handler.finite_elements->max_dofs_per_vertex() +
                            45*dof_handler.finite_elements->max_dofs_per_line() +
                            13*dof_handler.finite_elements->max_dofs_per_quad();
              break;
            case 8:
              max_couplings=37*dof_handler.finite_elements->max_dofs_per_vertex() +
                            56*dof_handler.finite_elements->max_dofs_per_line() +
                            16*dof_handler.finite_elements->max_dofs_per_quad();
              break;
            default:
              Assert (false, ExcNotImplemented());
              max_couplings=0;
            };
          return std::min(max_couplings,dof_handler.n_dofs());
        }


        template <int spacedim>
        static
        unsigned int
        max_couplings_between_dofs (const DoFHandler<3,spacedim> &dof_handler)
        {
//TODO:[?] Invent significantly better estimates than the ones in this function
          // doing the same thing here is a rather complicated thing,
          // compared to the 2d case, since it is hard to draw
          // pictures with several refined hexahedra :-) so I
          // presently only give a coarse estimate for the case that
          // at most 8 hexes meet at each vertex
          //
          // can anyone give better estimate here?
          const unsigned int max_adjacent_cells = dof_handler.tria->max_adjacent_cells();

          types::global_dof_index max_couplings;
          if (max_adjacent_cells <= 8)
            max_couplings=7*7*7*dof_handler.finite_elements->max_dofs_per_vertex() +
                          7*6*7*3*dof_handler.finite_elements->max_dofs_per_line() +
                          9*4*7*3*dof_handler.finite_elements->max_dofs_per_quad() +
                          27*dof_handler.finite_elements->max_dofs_per_hex();
          else
            {
              Assert (false, ExcNotImplemented());
              max_couplings=0;
            }

          return std::min(max_couplings,dof_handler.n_dofs());
        }
      };
    }
  }
}


namespace hp
{
  template <int dim, int spacedim>
  const unsigned int DoFHandler<dim,spacedim>::dimension;

  template <int dim, int spacedim>
  const types::global_dof_index DoFHandler<dim,spacedim>::invalid_dof_index;

  template <int dim, int spacedim>
  const unsigned int DoFHandler<dim,spacedim>::default_fe_index;



  template <int dim, int spacedim>
  DoFHandler<dim,spacedim>::DoFHandler (const Triangulation<dim,spacedim> &tria)
    :
    tria(&tria, typeid(*this).name()),
    faces (nullptr)
  {
    Assert ((dynamic_cast<const parallel::distributed::Triangulation< dim, spacedim >*>
             (&tria)
             == nullptr),
            ExcMessage ("The given triangulation is parallel distributed but "
                        "this class does not currently support this."));

    // decide whether we need a sequential or a parallel distributed policy
    if (dynamic_cast<const parallel::shared::Triangulation< dim, spacedim>*> (&*this->tria) != nullptr)
      policy.reset (new internal::DoFHandler::Policy::ParallelShared<DoFHandler<dim,spacedim> >(*this));
    else if (dynamic_cast<const parallel::distributed::Triangulation< dim, spacedim >*> (&*this->tria) != nullptr)
      //policy.reset (new internal::DoFHandler::Policy::ParallelDistributed<DoFHandler<dim,spacedim> >(*this));
      {
        Assert (false, ExcNotImplemented());
      }
    else
      policy.reset (new internal::DoFHandler::Policy::Sequential<DoFHandler<dim,spacedim> >(*this));

    create_active_fe_table ();

    tria_listeners.push_back
    (tria.signals.pre_refinement
     .connect (std::bind (&DoFHandler<dim,spacedim>::pre_refinement_action,
                          std::ref(*this))));
    tria_listeners.push_back
    (tria.signals.post_refinement
     .connect (std::bind (&DoFHandler<dim,spacedim>::post_refinement_action,
                          std::ref(*this))));
    tria_listeners.push_back
    (tria.signals.create
     .connect (std::bind (&DoFHandler<dim,spacedim>::post_refinement_action,
                          std::ref(*this))));
  }


  template <int dim, int spacedim>
  DoFHandler<dim,spacedim>::~DoFHandler ()
  {
    // unsubscribe as a listener to refinement of the underlying
    // triangulation
    for (unsigned int i=0; i<tria_listeners.size(); ++i)
      tria_listeners[i].disconnect ();
    tria_listeners.clear ();

    // ...and release allocated memory
    clear ();
  }


  /*------------------------ Cell iterator functions ------------------------*/

  template <int dim, int spacedim>
  typename DoFHandler<dim,spacedim>::cell_iterator
  DoFHandler<dim, spacedim>::begin(const unsigned int level) const
  {
    return cell_iterator (*this->get_triangulation().begin(level),
                          this);
  }



  template <int dim, int spacedim>
  typename DoFHandler<dim,spacedim>::active_cell_iterator
  DoFHandler<dim,spacedim>::begin_active (const unsigned int level) const
  {
    // level is checked in begin
    cell_iterator i = begin (level);
    if (i.state() != IteratorState::valid)
      return i;
    while (i->has_children())
      if ((++i).state() != IteratorState::valid)
        return i;
    return i;
  }



  template <int dim, int spacedim>
  typename DoFHandler<dim,spacedim>::cell_iterator
  DoFHandler<dim,spacedim>::end () const
  {
    return cell_iterator (&this->get_triangulation(),
                          -1,
                          -1,
                          this);
  }


  template <int dim, int spacedim>
  typename DoFHandler<dim,spacedim>::cell_iterator
  DoFHandler<dim,spacedim>::end (const unsigned int level) const
  {
    return (level == this->get_triangulation().n_levels()-1 ?
            end() :
            begin (level+1));
  }


  template <int dim, int spacedim>
  typename DoFHandler<dim, spacedim>::active_cell_iterator
  DoFHandler<dim, spacedim>::end_active (const unsigned int level) const
  {
    return (level == this->get_triangulation().n_levels()-1 ?
            active_cell_iterator(end()) :
            begin_active (level+1));
  }



  template <int dim, int spacedim>
  IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>
  DoFHandler<dim, spacedim>::cell_iterators () const
  {
    return
      IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>
      (begin(), end());
  }


  template <int dim, int spacedim>
  IteratorRange<typename DoFHandler<dim, spacedim>::active_cell_iterator>
  DoFHandler<dim, spacedim>::active_cell_iterators () const
  {
    return
      IteratorRange<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      (begin_active(), end());
  }



  template <int dim, int spacedim>
  IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>
  DoFHandler<dim, spacedim>::cell_iterators_on_level (const unsigned int level) const
  {
    return
      IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>
      (begin(level), end(level));
  }



  template <int dim, int spacedim>
  IteratorRange<typename DoFHandler<dim, spacedim>::active_cell_iterator>
  DoFHandler<dim, spacedim>::active_cell_iterators_on_level (const unsigned int level) const
  {
    return
      IteratorRange<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      (begin_active(level), end_active(level));
  }




//------------------------------------------------------------------


  template <int dim, int spacedim>
  types::global_dof_index DoFHandler<dim,spacedim>::n_boundary_dofs () const
  {
    Assert (finite_elements != nullptr, ExcNoFESelected());

    std::set<types::global_dof_index> boundary_dofs;
    std::vector<types::global_dof_index> dofs_on_face;
    dofs_on_face.reserve (this->get_fe ().max_dofs_per_face());

    // loop over all faces to check whether they are at a
    // boundary. note that we need not take special care of single
    // lines in 3d (using @p{cell->has_boundary_lines}), since we do
    // not support boundaries of dimension dim-2, and so every
    // boundary line is also part of a boundary face.
    typename HpDoFHandler<dim,spacedim>::active_cell_iterator cell = this->begin_active (),
                                                              endc = this->end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);

            cell->face(f)->get_dof_indices (dofs_on_face,
                                            cell->active_fe_index());
            for (unsigned int i=0; i<dofs_per_face; ++i)
              boundary_dofs.insert(dofs_on_face[i]);
          }
    return boundary_dofs.size();
  }



  template <int dim, int spacedim>
  types::global_dof_index
  DoFHandler<dim,spacedim>::n_boundary_dofs (const std::set<types::boundary_id> &boundary_ids) const
  {
    Assert (finite_elements != nullptr, ExcNoFESelected());
    Assert (boundary_ids.find (numbers::internal_face_boundary_id) == boundary_ids.end(),
            ExcInvalidBoundaryIndicator());

    // same as above, but with additional checks for set of boundary
    // indicators
    std::set<types::global_dof_index> boundary_dofs;
    std::vector<types::global_dof_index> dofs_on_face;
    dofs_on_face.reserve (this->get_fe ().max_dofs_per_face());

    typename HpDoFHandler<dim,spacedim>::active_cell_iterator cell = this->begin_active (),
                                                              endc = this->end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f) &&
            (boundary_ids.find(cell->face(f)->boundary_id()) !=
             boundary_ids.end()))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);

            cell->face(f)->get_dof_indices (dofs_on_face,
                                            cell->active_fe_index());
            for (unsigned int i=0; i<dofs_per_face; ++i)
              boundary_dofs.insert(dofs_on_face[i]);
          }
    return boundary_dofs.size();
  }



  template <>
  types::global_dof_index DoFHandler<2,3>::n_boundary_dofs () const
  {
    Assert(false,ExcNotImplemented());
    return 0;
  }



  template <>
  template <typename number>
  types::global_dof_index DoFHandler<2,3>::n_boundary_dofs (const std::map<types::boundary_id, const Function<3,number>*> &) const
  {
    Assert(false,ExcNotImplemented());
    return 0;
  }



  template <>
  types::global_dof_index DoFHandler<2,3>::n_boundary_dofs (const std::set<types::boundary_id> &) const
  {
    Assert(false,ExcNotImplemented());
    return 0;
  }



  template <int dim, int spacedim>
  std::size_t
  DoFHandler<dim,spacedim>::memory_consumption () const
  {
    std::size_t mem = (MemoryConsumption::memory_consumption (tria) +
                       MemoryConsumption::memory_consumption (finite_elements) +
                       MemoryConsumption::memory_consumption (tria) +
                       MemoryConsumption::memory_consumption (levels) +
                       MemoryConsumption::memory_consumption (*faces) +
                       MemoryConsumption::memory_consumption (number_cache) +
                       MemoryConsumption::memory_consumption (vertex_dofs) +
                       MemoryConsumption::memory_consumption (vertex_dofs_offsets) +
                       MemoryConsumption::memory_consumption (has_children));
    for (unsigned int i=0; i<levels.size(); ++i)
      mem += MemoryConsumption::memory_consumption (*levels[i]);
    mem += MemoryConsumption::memory_consumption (*faces);

    return mem;
  }






  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::set_active_fe_indices (const std::vector<unsigned int> &active_fe_indices)
  {
    Assert(active_fe_indices.size()==get_triangulation().n_active_cells(),
           ExcDimensionMismatch(active_fe_indices.size(), get_triangulation().n_active_cells()));

    create_active_fe_table ();
    // we could set the values directly, since they are stored as
    // protected data of this object, but for simplicity we use the
    // cell-wise access. this way we also have to pass some debug-mode
    // tests which we would have to duplicate ourselves otherwise
    active_cell_iterator cell=begin_active(),
                         endc=end();
    for (unsigned int i=0; cell!=endc; ++cell, ++i)
      cell->set_active_fe_index(active_fe_indices[i]);
  }



  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::get_active_fe_indices (std::vector<unsigned int> &active_fe_indices) const
  {
    active_fe_indices.resize(get_triangulation().n_active_cells());

    // we could try to extract the values directly, since they are
    // stored as protected data of this object, but for simplicity we
    // use the cell-wise access.
    active_cell_iterator cell=begin_active(),
                         endc=end();
    for (unsigned int i=0; cell!=endc; ++cell, ++i)
      active_fe_indices[i]=cell->active_fe_index();
  }



  namespace
  {
    /**
     * Given a hp::DoFHandler object, make sure that the active_fe_indices that
     * a user has set for locally owned cells are communicated to all ghost
     * cells as well.
     */
    template <int dim, int spacedim>
    void communicate_active_fe_indices (hp::DoFHandler<dim,spacedim> &dof_handler)
    {
      if (const parallel::shared::Triangulation<dim, spacedim> *tr =
            dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>*> (&dof_handler.get_triangulation()))
        {
          // we have a shared triangulation. in this case, every processor
          // knows about all cells, but every processor only has knowledge
          // about the active_fe_index on the cells it owns.
          //
          // we can create a complete set of active_fe_indices by letting
          // every processor create a vector of indices for all cells,
          // filling only those on the cells it owns and setting the indices
          // on the other cells to zero. then we add all of these vectors
          // up, and because every vector entry has exactly one processor
          // that owns it, the sum is correct
          std::vector<unsigned int> active_fe_indices (tr->n_active_cells(),
                                                       (unsigned int)0);
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              active_fe_indices[cell->active_cell_index()] = cell->active_fe_index();

          Utilities::MPI::sum (active_fe_indices,
                               tr->get_communicator(),
                               active_fe_indices);

          // now go back and fill the active_fe_index on ghost cells
          for (auto cell : dof_handler.active_cell_iterators())
            if (cell->is_ghost())
              cell->set_active_fe_index(active_fe_indices[cell->active_cell_index()]);
        }
      else if (const parallel::distributed::Triangulation<dim, spacedim> *tr =
                 dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>*> (&dof_handler.get_triangulation()))
        {
          Assert (false, ExcNotImplemented());
        }
      else
        {
          // a sequential triangulation. there is nothing we need to do here
          Assert ((dynamic_cast<const parallel::Triangulation<dim, spacedim>*> (&dof_handler.get_triangulation())
                   == nullptr),
                  ExcInternalError());
        }
    }
  }






  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::distribute_dofs (const hp::FECollection<dim,spacedim> &ff)
  {
    Assert (tria->n_levels() > 0, ExcInvalidTriangulation());

    finite_elements = &ff;

    // at the beginning, make sure every processor knows the
    // active_fe_indices on both its own cells and all ghost cells
    communicate_active_fe_indices (*this);

    // This call ensures that the active_fe_indices vectors are
    // initialized correctly.
    create_active_fe_table ();

    // up front make sure that the fe collection is large enough to
    // cover all fe indices presently in use on the mesh
    for (active_cell_iterator cell = begin_active(); cell != end(); ++cell)
      Assert (cell->active_fe_index() < finite_elements->size(),
              ExcInvalidFEIndex (cell->active_fe_index(),
                                 finite_elements->size()));


    // then allocate space for all the other tables
    dealii::internal::hp::DoFHandler::Implementation::reserve_space (*this);

    // Clear user flags because we will need them. But first we save
    // them and make sure that we restore them later such that at the
    // end of this function the Triangulation will be in the same
    // state as it was at the beginning of this function.
    std::vector<bool> user_flags;
    tria->save_user_flags(user_flags);
    const_cast<Triangulation<dim,spacedim> &>(*tria).clear_user_flags ();


    /////////////////////////////////

    // Now for the real work:
    number_cache = policy->distribute_dofs ();

    /////////////////////////////////


    Assert ((dynamic_cast<const parallel::distributed::Triangulation< dim, spacedim >*>
             (&this->get_triangulation())
             == nullptr),
            ExcNotImplemented());

    {
      Threads::TaskGroup<> tg;
      for (int level=levels.size()-1; level>=0; --level)
        tg += Threads::new_task (&dealii::internal::hp::DoFLevel::compress_data<dim,spacedim>,
                                 *levels[level], *finite_elements);
      tg.join_all ();
    }

    // finally restore the user flags
    const_cast<Triangulation<dim,spacedim> &>(*tria).load_user_flags(user_flags);
  }



  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::clear ()
  {
    // release lock to old fe
    finite_elements = nullptr;

    // release memory
    clear_space ();
  }



  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::renumber_dofs (const std::vector<types::global_dof_index> &new_numbers)
  {
    Assert(levels.size()>0, ExcMessage("You need to distribute DoFs before you can renumber them."));

    Assert (new_numbers.size() == n_locally_owned_dofs(),
            ExcRenumberingIncomplete());

#ifdef DEBUG
    // assert that the new indices are
    // consecutively numbered if we are
    // working on a single
    // processor. this doesn't need to
    // hold in the case of a parallel
    // mesh since we map the interval
    // [0...n_dofs()) into itself but
    // only globally, not on each
    // processor
    if (n_locally_owned_dofs() == n_dofs())
      {
        std::vector<types::global_dof_index> tmp(new_numbers);
        std::sort (tmp.begin(), tmp.end());
        std::vector<types::global_dof_index>::const_iterator p = tmp.begin();
        types::global_dof_index i = 0;
        for (; p!=tmp.end(); ++p, ++i)
          Assert (*p == i, ExcNewNumbersNotConsecutive(i));
      }
    else
      for (types::global_dof_index i=0; i<new_numbers.size(); ++i)
        Assert (new_numbers[i] < n_dofs(),
                ExcMessage ("New DoF index is not less than the total number of dofs."));
#endif

    // uncompress the internal storage scheme of dofs on cells so that
    // we can access dofs in turns. uncompress in parallel, starting
    // with the most expensive levels (the highest ones)
    {
      Threads::TaskGroup<> tg;
      for (int level=levels.size()-1; level>=0; --level)
        tg += Threads::new_task (&dealii::internal::hp::DoFLevel::uncompress_data<dim,spacedim>,
                                 *levels[level], *finite_elements);
      tg.join_all ();
    }

    // do the renumbering
    number_cache = policy->renumber_dofs(new_numbers);

    // now re-compress the dof indices
    {
      Threads::TaskGroup<> tg;
      for (int level=levels.size()-1; level>=0; --level)
        tg += Threads::new_task (&dealii::internal::hp::DoFLevel::compress_data<dim,spacedim>,
                                 *levels[level], *finite_elements);
      tg.join_all ();
    }
  }



  template <int dim, int spacedim>
  unsigned int
  DoFHandler<dim, spacedim>::max_couplings_between_dofs () const
  {
    Assert (finite_elements != nullptr, ExcNoFESelected());
    return dealii::internal::hp::DoFHandler::Implementation::max_couplings_between_dofs (*this);
  }



  template <int dim, int spacedim>
  unsigned int
  DoFHandler<dim,spacedim>::max_couplings_between_boundary_dofs () const
  {
    Assert (finite_elements != nullptr, ExcNoFESelected());

    switch (dim)
      {
      case 1:
        return finite_elements->max_dofs_per_vertex();
      case 2:
        return (3*finite_elements->max_dofs_per_vertex()
                +
                2*finite_elements->max_dofs_per_line());
      case 3:
        // we need to take refinement of one boundary face into
        // consideration here; in fact, this function returns what
        // #max_coupling_between_dofs<2> returns
        //
        // we assume here, that only four faces meet at the boundary;
        // this assumption is not justified and needs to be fixed some
        // time. fortunately, omitting it for now does no harm since
        // the matrix will cry foul if its requirements are not
        // satisfied
        return (19*finite_elements->max_dofs_per_vertex() +
                28*finite_elements->max_dofs_per_line() +
                8*finite_elements->max_dofs_per_quad());
      default:
        Assert (false, ExcNotImplemented());
        return 0;
      }
  }



  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::create_active_fe_table ()
  {
    // Create sufficiently many hp::DoFLevels.
    while (levels.size () < tria->n_levels ())
      levels.emplace_back (new dealii::internal::hp::DoFLevel);

    // then make sure that on each level we have the appropriate size
    // of active_fe_indices; preset them to zero, i.e. the default FE
    for (unsigned int level=0; level<levels.size(); ++level)
      {
        if (levels[level]->active_fe_indices.size () == 0)
          levels[level]->active_fe_indices.resize (tria->n_raw_cells(level),
                                                   0);
        else
          {
            // Either the active_fe_indices have size zero because
            // they were just created, or the correct size. Other
            // sizes indicate that something went wrong.
            Assert (levels[level]->active_fe_indices.size () ==
                    tria->n_raw_cells(level),
                    ExcInternalError ());
          }

        // it may be that the previous table was compressed; in that
        // case, restore the correct active_fe_index. the fact that
        // this no longer matches the indices in the table is of no
        // importance because the current function is called at a
        // point where we have to recreate the dof_indices tables in
        // the levels anyway
        levels[level]->normalize_active_fe_indices ();
      }
  }


  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::pre_refinement_action ()
  {
    create_active_fe_table ();

    // Remember if the cells already have children. That will make the
    // transfer of the active_fe_index to the finer levels easier.
    Assert (has_children.size () == 0, ExcInternalError ());
    for (unsigned int i=0; i<levels.size(); ++i)
      {
        const unsigned int cells_on_level = tria->n_raw_cells(i);
        std::unique_ptr<std::vector<bool> > has_children_level(new std::vector<bool> (cells_on_level));

        // Check for each cell, if it has children. in 1d, we don't
        // store refinement cases, so use the 'children' vector
        // instead
        if (dim == 1)
          std::transform (tria->levels[i]->cells.children.begin (),
                          tria->levels[i]->cells.children.end (),
                          has_children_level->begin (),
                          std::bind(std::not_equal_to<int>(),
                                    std::placeholders::_1,
                                    -1));
        else
          std::transform (tria->levels[i]->cells.refinement_cases.begin (),
                          tria->levels[i]->cells.refinement_cases.end (),
                          has_children_level->begin (),
                          std::bind (std::not_equal_to<unsigned char>(),
                                     std::placeholders::_1,
                                     static_cast<unsigned char>(RefinementCase<dim>::no_refinement)));

        has_children.emplace_back (std::move(has_children_level));
      }
  }



  template <int dim, int spacedim>
  void
  DoFHandler<dim,spacedim>::post_refinement_action ()
  {
    Assert (has_children.size () == levels.size (), ExcInternalError ());

    // Normally only one level is added, but if this Triangulation
    // is created by copy_triangulation, it can be more than one level.
    while (levels.size () < tria->n_levels ())
      levels.emplace_back (new dealii::internal::hp::DoFLevel);

    // Coarsening can lead to the loss of levels. Hence remove them.
    while (levels.size () > tria->n_levels ())
      {
        // drop the last element. that also releases the memory pointed to
        levels.pop_back ();
      }

    Assert(levels.size () == tria->n_levels (), ExcInternalError());

    // Resize active_fe_indices vectors. use zero indicator to extend
    for (unsigned int i=0; i<levels.size(); ++i)
      levels[i]->active_fe_indices.resize (tria->n_raw_cells(i), 0);

    // if a finite element collection has already been set, then
    // actually try to set active_fe_indices for child cells of
    // refined cells to the active_fe_index of the mother cell. if no
    // finite element collection has been assigned yet, then all
    // indicators are zero anyway, and there is no point trying to set
    // anything (besides, we would trip over an assertion in
    // set_active_fe_index)
    if (finite_elements != nullptr)
      {
        cell_iterator cell = begin(),
                      endc = end ();
        for (; cell != endc; ++cell)
          {
            // Look if the cell got children during refinement by
            // checking whether it has children now but didn't have
            // children before refinement (the has_children array is
            // set in pre-refinement action)
            //
            // Note: Although one level is added to the DoFHandler
            // levels, when the triangulation got one, for the buffer
            // has_children this new level is not required, because
            // the cells on the finest level never have
            // children. Hence cell->has_children () will always
            // return false on that level, which would cause shortcut
            // evaluation of the following expression. Thus an index
            // error in has_children should never occur.
            if (cell->has_children () &&
                !(*has_children [cell->level ()])[cell->index ()])
              {
                // Set active_fe_index in children to the same value
                // as in the parent cell. we can't access the
                // active_fe_index in the parent cell any more through
                // cell->active_fe_index() since that function is not
                // allowed for inactive cells, but we can access this
                // information from the DoFLevels directly
                for (unsigned int i = 0; i < cell->n_children(); ++i)
                  cell->child (i)->set_active_fe_index
                  (levels[cell->level()]->active_fe_index (cell->index()));
              }
          }
      }

    // Free buffer objects
    has_children.clear ();
  }


  template <int dim, int spacedim>
  template <int structdim>
  types::global_dof_index
  DoFHandler<dim,spacedim>::get_dof_index (const unsigned int,
                                           const unsigned int,
                                           const unsigned int,
                                           const unsigned int) const
  {
    Assert (false, ExcNotImplemented());
    return numbers::invalid_dof_index;
  }


  template <int dim, int spacedim>
  template <int structdim>
  void
  DoFHandler<dim,spacedim>::set_dof_index (const unsigned int,
                                           const unsigned int,
                                           const unsigned int,
                                           const unsigned int,
                                           const types::global_dof_index) const
  {
    Assert (false, ExcNotImplemented());
  }


  template <int dim, int spacedim>
  void DoFHandler<dim,spacedim>::clear_space ()
  {
    levels.clear ();
    faces.reset ();

    {
      std::vector<types::global_dof_index> tmp;
      std::swap (vertex_dofs, tmp);
    }

    {
      std::vector<types::global_dof_index> tmp;
      std::swap (vertex_dofs_offsets, tmp);
    }
  }
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler.inst"


DEAL_II_NAMESPACE_CLOSE
