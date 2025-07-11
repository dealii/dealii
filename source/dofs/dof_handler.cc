// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/distributed/cell_data_transfer.templates.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler_policy.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_levels.h>

#include <algorithm>
#include <memory>
#include <set>
#include <unordered_set>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
const types::fe_index DoFHandler<dim, spacedim>::default_fe_index;
#endif

namespace internal
{
  template <int dim, int spacedim>
  std::string
  policy_to_string(const dealii::internal::DoFHandlerImplementation::Policy::
                     PolicyBase<dim, spacedim> &policy)
  {
    std::string policy_name;
    if (dynamic_cast<const typename dealii::internal::DoFHandlerImplementation::
                       Policy::Sequential<dim, spacedim> *>(&policy))
      policy_name = "Policy::Sequential<";
    else if (dynamic_cast<
               const typename dealii::internal::DoFHandlerImplementation::
                 Policy::ParallelDistributed<dim, spacedim> *>(&policy))
      policy_name = "Policy::ParallelDistributed<";
    else if (dynamic_cast<
               const typename dealii::internal::DoFHandlerImplementation::
                 Policy::ParallelShared<dim, spacedim> *>(&policy))
      policy_name = "Policy::ParallelShared<";
    else
      AssertThrow(false, ExcNotImplemented());
    policy_name += Utilities::int_to_string(dim) + "," +
                   Utilities::int_to_string(spacedim) + ">";
    return policy_name;
  }


  namespace DoFHandlerImplementation
  {
    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Implement the function of same name in
       * the mother class.
       */
      template <int spacedim>
      static unsigned int
      max_couplings_between_dofs(const DoFHandler<1, spacedim> &dof_handler)
      {
        return std::min(static_cast<types::global_dof_index>(
                          3 * dof_handler.fe_collection.max_dofs_per_vertex() +
                          2 * dof_handler.fe_collection.max_dofs_per_line()),
                        dof_handler.n_dofs());
      }

      template <int spacedim>
      static unsigned int
      max_couplings_between_dofs(const DoFHandler<2, spacedim> &dof_handler)
      {
        // get these numbers by drawing pictures
        // and counting...
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
        // (the X vertices are connected with * if
        // the vertices adjacent to X are hanging
        // nodes)
        // count lines -> 28 (don't forget to count
        // mother and children separately!)
        types::global_dof_index max_couplings;
        switch (dof_handler.tria->max_adjacent_cells())
          {
            case 4:
              max_couplings =
                19 * dof_handler.fe_collection.max_dofs_per_vertex() +
                28 * dof_handler.fe_collection.max_dofs_per_line() +
                8 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 5:
              max_couplings =
                21 * dof_handler.fe_collection.max_dofs_per_vertex() +
                31 * dof_handler.fe_collection.max_dofs_per_line() +
                9 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 6:
              max_couplings =
                28 * dof_handler.fe_collection.max_dofs_per_vertex() +
                42 * dof_handler.fe_collection.max_dofs_per_line() +
                12 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 7:
              max_couplings =
                30 * dof_handler.fe_collection.max_dofs_per_vertex() +
                45 * dof_handler.fe_collection.max_dofs_per_line() +
                13 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 8:
              max_couplings =
                37 * dof_handler.fe_collection.max_dofs_per_vertex() +
                56 * dof_handler.fe_collection.max_dofs_per_line() +
                16 * dof_handler.fe_collection.max_dofs_per_quad();
              break;

            // the following numbers are not based on actual counting but by
            // extrapolating the number sequences from the previous ones (for
            // example, for n_dofs_per_vertex(), the sequence above is 19, 21,
            // 28, 30, 37, and is continued as follows):
            case 9:
              max_couplings =
                39 * dof_handler.fe_collection.max_dofs_per_vertex() +
                59 * dof_handler.fe_collection.max_dofs_per_line() +
                17 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 10:
              max_couplings =
                46 * dof_handler.fe_collection.max_dofs_per_vertex() +
                70 * dof_handler.fe_collection.max_dofs_per_line() +
                20 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 11:
              max_couplings =
                48 * dof_handler.fe_collection.max_dofs_per_vertex() +
                73 * dof_handler.fe_collection.max_dofs_per_line() +
                21 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 12:
              max_couplings =
                55 * dof_handler.fe_collection.max_dofs_per_vertex() +
                84 * dof_handler.fe_collection.max_dofs_per_line() +
                24 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 13:
              max_couplings =
                57 * dof_handler.fe_collection.max_dofs_per_vertex() +
                87 * dof_handler.fe_collection.max_dofs_per_line() +
                25 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 14:
              max_couplings =
                63 * dof_handler.fe_collection.max_dofs_per_vertex() +
                98 * dof_handler.fe_collection.max_dofs_per_line() +
                28 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 15:
              max_couplings =
                65 * dof_handler.fe_collection.max_dofs_per_vertex() +
                103 * dof_handler.fe_collection.max_dofs_per_line() +
                29 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 16:
              max_couplings =
                72 * dof_handler.fe_collection.max_dofs_per_vertex() +
                114 * dof_handler.fe_collection.max_dofs_per_line() +
                32 * dof_handler.fe_collection.max_dofs_per_quad();
              break;

            default:
              DEAL_II_NOT_IMPLEMENTED();
              max_couplings = 0;
          }
        return std::min(max_couplings, dof_handler.n_dofs());
      }

      template <int spacedim>
      static unsigned int
      max_couplings_between_dofs(const DoFHandler<3, spacedim> &dof_handler)
      {
        // TODO:[?] Invent significantly better estimates than the ones in this
        // function

        // doing the same thing here is a rather complicated thing, compared
        // to the 2d case, since it is hard to draw pictures with several
        // refined hexahedra :-) so I presently only give a coarse
        // estimate for the case that at most 8 hexes meet at each vertex
        //
        // can anyone give better estimate here?
        const unsigned int max_adjacent_cells =
          dof_handler.tria->max_adjacent_cells();

        types::global_dof_index max_couplings;
        if (max_adjacent_cells <= 8)
          max_couplings =
            7 * 7 * 7 * dof_handler.fe_collection.max_dofs_per_vertex() +
            7 * 6 * 7 * 3 * dof_handler.fe_collection.max_dofs_per_line() +
            9 * 4 * 7 * 3 * dof_handler.fe_collection.max_dofs_per_quad() +
            27 * dof_handler.fe_collection.max_dofs_per_hex();
        else
          {
            DEAL_II_NOT_IMPLEMENTED();
            max_couplings = 0;
          }

        return std::min(max_couplings, dof_handler.n_dofs());
      }

      /**
       * Do that part of reserving space that pertains to releasing
       * the previously used memory.
       */
      template <int dim, int spacedim>
      static void
      reset_to_empty_objects(DoFHandler<dim, spacedim> &dof_handler)
      {
        dof_handler.object_dof_indices.clear();
        dof_handler.object_dof_indices.resize(dof_handler.tria->n_levels());
        dof_handler.object_dof_indices.shrink_to_fit();

        dof_handler.object_dof_ptr.clear();
        dof_handler.object_dof_ptr.resize(dof_handler.tria->n_levels());
        dof_handler.object_dof_ptr.shrink_to_fit();
      }

      /**
       * Reserve space for non-artificial cells.
       */
      template <int dim, int spacedim>
      static void
      reserve_cells(DoFHandler<dim, spacedim> &dof_handler,
                    const unsigned int         n_inner_dofs_per_cell)
      {
        for (unsigned int i = 0; i < dof_handler.tria->n_levels(); ++i)
          {
            dof_handler.object_dof_ptr[i][dim].assign(
              dof_handler.tria->n_raw_cells(i) + 1, 0);

            for (const auto &cell :
                 dof_handler.tria->cell_iterators_on_level(i))
              if (cell->is_active() && !cell->is_artificial())
                dof_handler.object_dof_ptr[i][dim][cell->index() + 1] =
                  n_inner_dofs_per_cell;

            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i); ++j)
              dof_handler.object_dof_ptr[i][dim][j + 1] +=
                dof_handler.object_dof_ptr[i][dim][j];

            dof_handler.object_dof_indices[i][dim].resize(
              dof_handler.object_dof_ptr[i][dim].back(),
              numbers::invalid_dof_index);
          }
      }

      /**
       * Reserve space for @p structdim-dimensional objects connected to
       * non-artificial cells.
       */
      template <int dim, int spacedim, typename T>
      static void
      reserve_subentities(DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int         structdim,
                          const unsigned int         n_raw_entities,
                          const T                   &cell_process)
      {
        if (dof_handler.tria->n_cells() == 0)
          return;

        dof_handler.object_dof_ptr[0][structdim].assign(n_raw_entities + 1, -1);
        // determine for each entity the number of dofs
        for (const auto &cell : dof_handler.tria->cell_iterators())
          if (cell->is_active() && !cell->is_artificial())
            cell_process(
              cell,
              [&](const unsigned int n_dofs_per_entity,
                  const unsigned int index) {
                auto &n_dofs_per_entity_target =
                  dof_handler.object_dof_ptr[0][structdim][index + 1];

                // make sure that either the entity has not been visited or
                // the entity has the same number of dofs assigned
                Assert((n_dofs_per_entity_target ==
                          static_cast<
                            typename DoFHandler<dim, spacedim>::offset_type>(
                            -1) ||
                        n_dofs_per_entity_target == n_dofs_per_entity),
                       ExcNotImplemented());

                n_dofs_per_entity_target = n_dofs_per_entity;
              });

        // convert the absolute numbers to CRS
        dof_handler.object_dof_ptr[0][structdim][0] = 0;
        for (unsigned int i = 1; i < n_raw_entities + 1; ++i)
          {
            if (dof_handler.object_dof_ptr[0][structdim][i] ==
                static_cast<typename DoFHandler<dim, spacedim>::offset_type>(
                  -1))
              dof_handler.object_dof_ptr[0][structdim][i] =
                dof_handler.object_dof_ptr[0][structdim][i - 1];
            else
              dof_handler.object_dof_ptr[0][structdim][i] +=
                dof_handler.object_dof_ptr[0][structdim][i - 1];
          }

        // allocate memory for indices
        dof_handler.object_dof_indices[0][structdim].resize(
          dof_handler.object_dof_ptr[0][structdim].back(),
          numbers::invalid_dof_index);
      }

      /**
       * Reserve enough space in the <tt>levels[]</tt> objects to store the
       * numbers of the degrees of freedom needed for the given element. The
       * given element is that one which was selected when calling
       * @p distribute_dofs the last time.
       */
      template <int dim, int spacedim>
      static void
      reserve_space(DoFHandler<dim, spacedim> &dof_handler)
      {
        reset_to_empty_objects(dof_handler);

        const auto &fe = dof_handler.get_fe();

        // cell
        reserve_cells(dof_handler,
                      dim == 1 ? fe.n_dofs_per_line() :
                                 (dim == 2 ? fe.n_dofs_per_quad(0) :
                                             fe.n_dofs_per_hex()));

        // vertices
        reserve_subentities(dof_handler,
                            0,
                            dof_handler.tria->n_vertices(),
                            [&](const auto &cell, const auto &process) {
                              for (const auto vertex_index :
                                   cell->vertex_indices())
                                process(fe.n_dofs_per_vertex(),
                                        cell->vertex_index(vertex_index));
                            });

        // lines
        if (dim == 2 || dim == 3)
          reserve_subentities(dof_handler,
                              1,
                              dof_handler.tria->n_raw_lines(),
                              [&](const auto &cell, const auto &process) {
                                for (const auto line_index :
                                     cell->line_indices())
                                  process(fe.n_dofs_per_line(),
                                          cell->line(line_index)->index());
                              });

        // quads
        if (dim == 3)
          reserve_subentities(dof_handler,
                              2,
                              dof_handler.tria->n_raw_quads(),
                              [&](const auto &cell, const auto &process) {
                                for (const auto face_index :
                                     cell->face_indices())
                                  process(fe.n_dofs_per_quad(face_index),
                                          cell->face(face_index)->index());
                              });
      }

      template <int spacedim>
      static void
      reserve_space_mg(DoFHandler<1, spacedim> &dof_handler)
      {
        Assert(dof_handler.get_triangulation().n_levels() > 0,
               ExcMessage("Invalid triangulation"));
        dof_handler.clear_mg_space();

        const dealii::Triangulation<1, spacedim> &tria =
          dof_handler.get_triangulation();
        const unsigned int dofs_per_line =
          dof_handler.get_fe().n_dofs_per_line();
        const unsigned int n_levels = tria.n_levels();

        for (unsigned int i = 0; i < n_levels; ++i)
          {
            dof_handler.mg_levels.emplace_back(
              new internal::DoFHandlerImplementation::DoFLevel<1>);
            dof_handler.mg_levels.back()->dof_object.dofs =
              std::vector<types::global_dof_index>(tria.n_raw_lines(i) *
                                                     dofs_per_line,
                                                   numbers::invalid_dof_index);
          }

        const unsigned int n_vertices = tria.n_vertices();

        dof_handler.mg_vertex_dofs.resize(n_vertices);

        std::vector<unsigned int> max_level(n_vertices, 0);
        std::vector<unsigned int> min_level(n_vertices, n_levels);

        for (typename dealii::Triangulation<1, spacedim>::cell_iterator cell =
               tria.begin();
             cell != tria.end();
             ++cell)
          {
            const unsigned int level = cell->level();

            for (const auto vertex : cell->vertex_indices())
              {
                const unsigned int vertex_index = cell->vertex_index(vertex);

                if (min_level[vertex_index] > level)
                  min_level[vertex_index] = level;

                if (max_level[vertex_index] < level)
                  max_level[vertex_index] = level;
              }
          }

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          if (tria.vertex_used(vertex))
            {
              Assert(min_level[vertex] < n_levels, ExcInternalError());
              Assert(max_level[vertex] >= min_level[vertex],
                     ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(
                min_level[vertex],
                max_level[vertex],
                dof_handler.get_fe().n_dofs_per_vertex());
            }

          else
            {
              Assert(min_level[vertex] == n_levels, ExcInternalError());
              Assert(max_level[vertex] == 0, ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(1, 0, 0);
            }
      }

      template <int spacedim>
      static void
      reserve_space_mg(DoFHandler<2, spacedim> &dof_handler)
      {
        Assert(dof_handler.get_triangulation().n_levels() > 0,
               ExcMessage("Invalid triangulation"));
        dof_handler.clear_mg_space();

        const dealii::FiniteElement<2, spacedim> &fe = dof_handler.get_fe();
        const dealii::Triangulation<2, spacedim> &tria =
          dof_handler.get_triangulation();
        const unsigned int n_levels = tria.n_levels();

        for (unsigned int i = 0; i < n_levels; ++i)
          {
            dof_handler.mg_levels.emplace_back(
              std::make_unique<
                internal::DoFHandlerImplementation::DoFLevel<2>>());
            dof_handler.mg_levels.back()->dof_object.dofs =
              std::vector<types::global_dof_index>(
                tria.n_raw_quads(i) *
                  fe.n_dofs_per_quad(0 /*note: in 2d there is only one quad*/),
                numbers::invalid_dof_index);
          }

        dof_handler.mg_faces =
          std::make_unique<internal::DoFHandlerImplementation::DoFFaces<2>>();
        dof_handler.mg_faces->lines.dofs =
          std::vector<types::global_dof_index>(tria.n_raw_lines() *
                                                 fe.n_dofs_per_line(),
                                               numbers::invalid_dof_index);

        const unsigned int n_vertices = tria.n_vertices();

        dof_handler.mg_vertex_dofs.resize(n_vertices);

        std::vector<unsigned int> max_level(n_vertices, 0);
        std::vector<unsigned int> min_level(n_vertices, n_levels);

        for (typename dealii::Triangulation<2, spacedim>::cell_iterator cell =
               tria.begin();
             cell != tria.end();
             ++cell)
          {
            const unsigned int level = cell->level();

            for (const auto vertex : cell->vertex_indices())
              {
                const unsigned int vertex_index = cell->vertex_index(vertex);

                if (min_level[vertex_index] > level)
                  min_level[vertex_index] = level;

                if (max_level[vertex_index] < level)
                  max_level[vertex_index] = level;
              }
          }

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          if (tria.vertex_used(vertex))
            {
              Assert(min_level[vertex] < n_levels, ExcInternalError());
              Assert(max_level[vertex] >= min_level[vertex],
                     ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(min_level[vertex],
                                                      max_level[vertex],
                                                      fe.n_dofs_per_vertex());
            }

          else
            {
              Assert(min_level[vertex] == n_levels, ExcInternalError());
              Assert(max_level[vertex] == 0, ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(1, 0, 0);
            }
      }

      template <int spacedim>
      static void
      reserve_space_mg(DoFHandler<3, spacedim> &dof_handler)
      {
        Assert(dof_handler.get_triangulation().n_levels() > 0,
               ExcMessage("Invalid triangulation"));
        dof_handler.clear_mg_space();

        const dealii::FiniteElement<3, spacedim> &fe = dof_handler.get_fe();
        const dealii::Triangulation<3, spacedim> &tria =
          dof_handler.get_triangulation();
        const unsigned int n_levels = tria.n_levels();

        for (unsigned int i = 0; i < n_levels; ++i)
          {
            dof_handler.mg_levels.emplace_back(
              std::make_unique<
                internal::DoFHandlerImplementation::DoFLevel<3>>());
            dof_handler.mg_levels.back()->dof_object.dofs =
              std::vector<types::global_dof_index>(tria.n_raw_hexs(i) *
                                                     fe.n_dofs_per_hex(),
                                                   numbers::invalid_dof_index);
          }

        dof_handler.mg_faces =
          std::make_unique<internal::DoFHandlerImplementation::DoFFaces<3>>();
        dof_handler.mg_faces->lines.dofs =
          std::vector<types::global_dof_index>(tria.n_raw_lines() *
                                                 fe.n_dofs_per_line(),
                                               numbers::invalid_dof_index);

        // TODO: the implementation makes the assumption that all faces have the
        // same number of dofs
        AssertDimension(fe.n_unique_faces(), 1);
        dof_handler.mg_faces->quads.dofs = std::vector<types::global_dof_index>(
          tria.n_raw_quads() * fe.n_dofs_per_quad(0 /*=face_no*/),
          numbers::invalid_dof_index);

        const unsigned int n_vertices = tria.n_vertices();

        dof_handler.mg_vertex_dofs.resize(n_vertices);

        std::vector<unsigned int> max_level(n_vertices, 0);
        std::vector<unsigned int> min_level(n_vertices, n_levels);

        for (typename dealii::Triangulation<3, spacedim>::cell_iterator cell =
               tria.begin();
             cell != tria.end();
             ++cell)
          {
            const unsigned int level = cell->level();

            for (const auto vertex : cell->vertex_indices())
              {
                const unsigned int vertex_index = cell->vertex_index(vertex);

                if (min_level[vertex_index] > level)
                  min_level[vertex_index] = level;

                if (max_level[vertex_index] < level)
                  max_level[vertex_index] = level;
              }
          }

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          if (tria.vertex_used(vertex))
            {
              Assert(min_level[vertex] < n_levels, ExcInternalError());
              Assert(max_level[vertex] >= min_level[vertex],
                     ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(min_level[vertex],
                                                      max_level[vertex],
                                                      fe.n_dofs_per_vertex());
            }

          else
            {
              Assert(min_level[vertex] == n_levels, ExcInternalError());
              Assert(max_level[vertex] == 0, ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(1, 0, 0);
            }
      }
    };
  } // namespace DoFHandlerImplementation



  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      /**
       * A class with the same purpose as the similarly named class of the
       * Triangulation class. See there for more information.
       */
      struct Implementation
      {
        /**
         * No future FE indices should have been assigned when partitioning a
         * triangulation, since they are only available locally and will not be
         * communicated.
         */
        template <int dim, int spacedim>
        static void
        ensure_absence_of_future_fe_indices(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          (void)dof_handler;
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              Assert(
                !cell->future_fe_index_set(),
                ExcMessage(
                  "There shouldn't be any cells flagged for p-adaptation when partitioning."));
        }



        /**
         * Do that part of reserving space that pertains to vertices,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_vertices(DoFHandler<dim, spacedim> &dof_handler)
        {
          // The final step in all of the reserve_space() functions is to set
          // up vertex dof information. since vertices are sequentially
          // numbered, what we do first is to set up an array in which
          // we record whether a vertex is associated with any of the
          // given fe's, by setting a bit. in a later step, we then
          // actually allocate memory for the required dofs
          //
          // in the following, we only need to consider vertices that are
          // adjacent to either a locally owned or a ghost cell; we never
          // store anything on vertices that are only surrounded by
          // artificial cells. so figure out that subset of vertices
          // first
          std::vector<bool> locally_used_vertices(
            dof_handler.tria->n_vertices(), false);
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              for (const auto v : cell->vertex_indices())
                locally_used_vertices[cell->vertex_index(v)] = true;

          std::vector<std::vector<bool>> vertex_fe_association(
            dof_handler.fe_collection.size(),
            std::vector<bool>(dof_handler.tria->n_vertices(), false));

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              for (const auto v : cell->vertex_indices())
                vertex_fe_association[cell->active_fe_index()]
                                     [cell->vertex_index(v)] = true;

          // in debug mode, make sure that each vertex is associated
          // with at least one FE (note that except for unused
          // vertices, all vertices are actually active). this is of
          // course only true for vertices that are part of either
          // ghost or locally owned cells
          if constexpr (running_in_debug_mode())
            {
              for (unsigned int v = 0; v < dof_handler.tria->n_vertices(); ++v)
                if (locally_used_vertices[v] == true)
                  if (dof_handler.tria->vertex_used(v) == true)
                    {
                      unsigned int fe = 0;
                      for (; fe < dof_handler.fe_collection.size(); ++fe)
                        if (vertex_fe_association[fe][v] == true)
                          break;
                      Assert(fe != dof_handler.fe_collection.size(),
                             ExcInternalError());
                    }
            }

          const unsigned int d = 0;
          const unsigned int l = 0;

          dof_handler.hp_object_fe_ptr[d].clear();
          dof_handler.hp_object_fe_indices[d].clear();
          dof_handler.object_dof_ptr[l][d].clear();
          dof_handler.object_dof_indices[l][d].clear();

          dof_handler.hp_object_fe_ptr[d].reserve(
            dof_handler.tria->n_vertices() + 1);

          unsigned int vertex_slots_needed = 0;
          unsigned int fe_slots_needed     = 0;

          for (unsigned int v = 0; v < dof_handler.tria->n_vertices(); ++v)
            {
              dof_handler.hp_object_fe_ptr[d].push_back(fe_slots_needed);

              if (dof_handler.tria->vertex_used(v) && locally_used_vertices[v])
                {
                  for (unsigned int fe = 0;
                       fe < dof_handler.fe_collection.size();
                       ++fe)
                    if (vertex_fe_association[fe][v] == true)
                      {
                        ++fe_slots_needed;
                        vertex_slots_needed +=
                          dof_handler.get_fe(fe).n_dofs_per_vertex();
                      }
                }
            }

          dof_handler.hp_object_fe_ptr[d].push_back(fe_slots_needed);

          dof_handler.hp_object_fe_indices[d].reserve(fe_slots_needed);
          dof_handler.object_dof_ptr[l][d].reserve(fe_slots_needed + 1);

          dof_handler.object_dof_indices[l][d].reserve(vertex_slots_needed);

          for (unsigned int v = 0; v < dof_handler.tria->n_vertices(); ++v)
            if (dof_handler.tria->vertex_used(v) && locally_used_vertices[v])
              {
                for (unsigned int fe = 0; fe < dof_handler.fe_collection.size();
                     ++fe)
                  if (vertex_fe_association[fe][v] == true)
                    {
                      dof_handler.hp_object_fe_indices[d].push_back(fe);
                      dof_handler.object_dof_ptr[l][d].push_back(
                        dof_handler.object_dof_indices[l][d].size());

                      for (unsigned int i = 0;
                           i < dof_handler.get_fe(fe).n_dofs_per_vertex();
                           i++)
                        dof_handler.object_dof_indices[l][d].push_back(
                          numbers::invalid_dof_index);
                    }
              }


          dof_handler.object_dof_ptr[l][d].push_back(
            dof_handler.object_dof_indices[l][d].size());

          AssertDimension(vertex_slots_needed,
                          dof_handler.object_dof_indices[l][d].size());
          AssertDimension(fe_slots_needed,
                          dof_handler.hp_object_fe_indices[d].size());
          AssertDimension(fe_slots_needed + 1,
                          dof_handler.object_dof_ptr[l][d].size());
          AssertDimension(dof_handler.tria->n_vertices() + 1,
                          dof_handler.hp_object_fe_ptr[d].size());

          dof_handler.object_dof_indices[l][d].assign(
            vertex_slots_needed, numbers::invalid_dof_index);
        }



        /**
         * Do that part of reserving space that pertains to cells,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_cells(DoFHandler<dim, spacedim> &dof_handler)
        {
          (void)dof_handler;
          // count how much space we need on each level for the cell
          // dofs and set the dof_*_offsets data. initially set the
          // latter to an invalid index, and only later set it to
          // something reasonable for active dof_handler.cells
          //
          // note that for dof_handler.cells, the situation is simpler
          // than for other (lower dimensional) objects since exactly
          // one finite element is used for it
          for (unsigned int level = 0; level < dof_handler.tria->n_levels();
               ++level)
            {
              dof_handler.object_dof_ptr[level][dim] =
                std::vector<typename DoFHandler<dim, spacedim>::offset_type>(
                  dof_handler.tria->n_raw_cells(level),
                  static_cast<typename DoFHandler<dim, spacedim>::offset_type>(
                    -1));

              types::global_dof_index next_free_dof = 0;
              for (auto cell :
                   dof_handler.active_cell_iterators_on_level(level))
                if (cell->is_active() && !cell->is_artificial())
                  {
                    dof_handler.object_dof_ptr[level][dim][cell->index()] =
                      next_free_dof;
                    next_free_dof +=
                      cell->get_fe().template n_dofs_per_object<dim>();
                  }

              dof_handler.object_dof_indices[level][dim] =
                std::vector<types::global_dof_index>(
                  next_free_dof, numbers::invalid_dof_index);
            }
        }



        /**
         * Do that part of reserving space that pertains to faces,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_faces(DoFHandler<dim, spacedim> &dof_handler)
        {
          // FACE DOFS
          //
          // Count face dofs, then allocate as much space
          // as we need and prime the linked list for faces (see the
          // description in hp::DoFLevel) with the indices we will
          // need. Note that our task is more complicated than for the
          // cell case above since two adjacent cells may have different
          // active FE indices, in which case we need to allocate
          // *two* sets of face dofs for the same face. But they don't
          // *have* to be different, and so we need to prepare for this
          // as well.
          //
          // The way we do things is that we loop over all active cells (these
          // are the only ones that have DoFs anyway) and all their faces. We
          // note in the vector face_touched whether we have previously
          // visited a face and if so skip it
          {
            std::vector<bool> face_touched(dim == 2 ?
                                             dof_handler.tria->n_raw_lines() :
                                             dof_handler.tria->n_raw_quads());

            const unsigned int d = dim - 1;
            const unsigned int l = 0;

            dof_handler.hp_object_fe_ptr[d].clear();
            dof_handler.hp_object_fe_indices[d].clear();
            dof_handler.object_dof_ptr[l][d].clear();
            dof_handler.object_dof_indices[l][d].clear();

            dof_handler.hp_object_fe_ptr[d].resize(
              dof_handler.tria->n_raw_faces() + 1);

            // An array to hold how many slots (see the hp::DoFLevel
            // class) we will have to store on each level
            unsigned int n_face_slots = 0;

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const auto face : cell->face_indices())
                  if (!face_touched[cell->face(face)->index()])
                    {
                      unsigned int fe_slots_needed = 0;

                      if (cell->at_boundary(face) ||
                          cell->face(face)->has_children() ||
                          cell->neighbor_is_coarser(face) ||
                          (!cell->at_boundary(face) &&
                           cell->neighbor(face)->is_artificial()) ||
                          (!cell->at_boundary(face) &&
                           !cell->neighbor(face)->is_artificial() &&
                           (cell->active_fe_index() ==
                            cell->neighbor(face)->active_fe_index())))
                        {
                          fe_slots_needed = 1;
                          n_face_slots +=
                            dof_handler.get_fe(cell->active_fe_index())
                              .template n_dofs_per_object<dim - 1>(face);
                        }
                      else
                        {
                          fe_slots_needed = 2;
                          n_face_slots +=
                            dof_handler.get_fe(cell->active_fe_index())
                              .template n_dofs_per_object<dim - 1>(face) +
                            dof_handler
                              .get_fe(cell->neighbor(face)->active_fe_index())
                              .template n_dofs_per_object<dim - 1>(
                                cell->neighbor_face_no(face));
                        }

                      // mark this face as visited
                      face_touched[cell->face(face)->index()] = true;

                      dof_handler
                        .hp_object_fe_ptr[d][cell->face(face)->index() + 1] =
                        fe_slots_needed;
                    }

            for (unsigned int i = 1; i < dof_handler.hp_object_fe_ptr[d].size();
                 i++)
              dof_handler.hp_object_fe_ptr[d][i] +=
                dof_handler.hp_object_fe_ptr[d][i - 1];


            dof_handler.hp_object_fe_indices[d].resize(
              dof_handler.hp_object_fe_ptr[d].back());
            dof_handler.object_dof_ptr[l][d].resize(
              dof_handler.hp_object_fe_ptr[d].back() + 1);

            dof_handler.object_dof_indices[l][d].reserve(n_face_slots);


            // With the memory now allocated, loop over the
            // dof_handler cells again and prime the _offset values as
            // well as the fe_index fields
            face_touched = std::vector<bool>(face_touched.size());

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const auto face : cell->face_indices())
                  if (!face_touched[cell->face(face)->index()])
                    {
                      // Same decision tree as before
                      if (cell->at_boundary(face) ||
                          cell->face(face)->has_children() ||
                          cell->neighbor_is_coarser(face) ||
                          (!cell->at_boundary(face) &&
                           cell->neighbor(face)->is_artificial()) ||
                          (!cell->at_boundary(face) &&
                           !cell->neighbor(face)->is_artificial() &&
                           (cell->active_fe_index() ==
                            cell->neighbor(face)->active_fe_index())))
                        {
                          const types::fe_index fe = cell->active_fe_index();
                          const unsigned int    n_dofs =
                            dof_handler.get_fe(fe)
                              .template n_dofs_per_object<dim - 1>(face);
                          const unsigned int offset =
                            dof_handler
                              .hp_object_fe_ptr[d][cell->face(face)->index()];

                          dof_handler.hp_object_fe_indices[d][offset]  = fe;
                          dof_handler.object_dof_ptr[l][d][offset + 1] = n_dofs;

                          for (unsigned int i = 0; i < n_dofs; ++i)
                            dof_handler.object_dof_indices[l][d].push_back(
                              numbers::invalid_dof_index);
                        }
                      else
                        {
                          types::fe_index fe_1      = cell->active_fe_index();
                          unsigned int    face_no_1 = face;
                          types::fe_index fe_2 =
                            cell->neighbor(face)->active_fe_index();
                          unsigned int face_no_2 = cell->neighbor_face_no(face);

                          if (fe_2 < fe_1)
                            {
                              std::swap(fe_1, fe_2);
                              std::swap(face_no_1, face_no_2);
                            }

                          const unsigned int n_dofs_1 =
                            dof_handler.get_fe(fe_1)
                              .template n_dofs_per_object<dim - 1>(face_no_1);

                          const unsigned int n_dofs_2 =
                            dof_handler.get_fe(fe_2)
                              .template n_dofs_per_object<dim - 1>(face_no_2);

                          const unsigned int offset =
                            dof_handler
                              .hp_object_fe_ptr[d][cell->face(face)->index()];

                          dof_handler.hp_object_fe_indices[d].push_back(
                            cell->active_fe_index());
                          dof_handler.object_dof_ptr[l][d].push_back(
                            dof_handler.object_dof_indices[l][d].size());

                          dof_handler.hp_object_fe_indices[d][offset + 0] =
                            fe_1;
                          dof_handler.hp_object_fe_indices[d][offset + 1] =
                            fe_2;
                          dof_handler.object_dof_ptr[l][d][offset + 1] =
                            n_dofs_1;
                          dof_handler.object_dof_ptr[l][d][offset + 2] =
                            n_dofs_2;


                          for (unsigned int i = 0; i < n_dofs_1 + n_dofs_2; ++i)
                            dof_handler.object_dof_indices[l][d].push_back(
                              numbers::invalid_dof_index);
                        }

                      // mark this face as visited
                      face_touched[cell->face(face)->index()] = true;
                    }

            for (unsigned int i = 1;
                 i < dof_handler.object_dof_ptr[l][d].size();
                 i++)
              dof_handler.object_dof_ptr[l][d][i] +=
                dof_handler.object_dof_ptr[l][d][i - 1];
          }
        }



        /**
         * Reserve enough space in the <tt>levels[]</tt> objects to
         * store the numbers of the degrees of freedom needed for the
         * given element. The given element is that one which was
         * selected when calling @p distribute_dofs the last time.
         */
        template <int spacedim>
        static void
        reserve_space(DoFHandler<1, spacedim> &dof_handler)
        {
          Assert(dof_handler.fe_collection.size() > 0,
                 (typename DoFHandler<1, spacedim>::ExcNoFESelected()));
          Assert(dof_handler.tria->n_levels() > 0,
                 ExcMessage("The current Triangulation must not be empty."));
          Assert(dof_handler.tria->n_levels() ==
                   dof_handler.hp_cell_future_fe_indices.size(),
                 ExcInternalError());

          internal::DoFHandlerImplementation::Implementation::
            reset_to_empty_objects(dof_handler);

          Threads::TaskGroup<> tasks;
          tasks +=
            Threads::new_task(&reserve_space_cells<1, spacedim>, dof_handler);
          tasks += Threads::new_task(&reserve_space_vertices<1, spacedim>,
                                     dof_handler);
          tasks.join_all();
        }



        template <int spacedim>
        static void
        reserve_space(DoFHandler<2, spacedim> &dof_handler)
        {
          Assert(dof_handler.fe_collection.size() > 0,
                 (typename DoFHandler<1, spacedim>::ExcNoFESelected()));
          Assert(dof_handler.tria->n_levels() > 0,
                 ExcMessage("The current Triangulation must not be empty."));
          Assert(dof_handler.tria->n_levels() ==
                   dof_handler.hp_cell_future_fe_indices.size(),
                 ExcInternalError());

          internal::DoFHandlerImplementation::Implementation::
            reset_to_empty_objects(dof_handler);

          Threads::TaskGroup<> tasks;
          tasks +=
            Threads::new_task(&reserve_space_cells<2, spacedim>, dof_handler);
          tasks +=
            Threads::new_task(&reserve_space_faces<2, spacedim>, dof_handler);
          tasks += Threads::new_task(&reserve_space_vertices<2, spacedim>,
                                     dof_handler);
          tasks.join_all();
        }



        template <int spacedim>
        static void
        reserve_space(DoFHandler<3, spacedim> &dof_handler)
        {
          Assert(dof_handler.fe_collection.size() > 0,
                 (typename DoFHandler<1, spacedim>::ExcNoFESelected()));
          Assert(dof_handler.tria->n_levels() > 0,
                 ExcMessage("The current Triangulation must not be empty."));
          Assert(dof_handler.tria->n_levels() ==
                   dof_handler.hp_cell_future_fe_indices.size(),
                 ExcInternalError());

          internal::DoFHandlerImplementation::Implementation::
            reset_to_empty_objects(dof_handler);

          Threads::TaskGroup<> tasks;
          tasks +=
            Threads::new_task(&reserve_space_cells<3, spacedim>, dof_handler);
          tasks +=
            Threads::new_task(&reserve_space_faces<3, spacedim>, dof_handler);
          tasks += Threads::new_task(&reserve_space_vertices<3, spacedim>,
                                     dof_handler);

          // While the tasks above are running, we can turn to line dofs

          // the situation here is pretty much like with vertices:
          // there can be an arbitrary number of finite elements
          // associated with each line.
          //
          // the algorithm we use is somewhat similar to what we do in
          // reserve_space_vertices()
          {
            // what we do first is to set up an array in which we
            // record whether a line is associated with any of the
            // given fe's, by setting a bit. in a later step, we
            // then actually allocate memory for the required dofs
            std::vector<std::vector<bool>> line_fe_association(
              dof_handler.fe_collection.size(),
              std::vector<bool>(dof_handler.tria->n_raw_lines(), false));

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const auto l : cell->line_indices())
                  line_fe_association[cell->active_fe_index()]
                                     [cell->line_index(l)] = true;

            // first check which of the lines is used at all,
            // i.e. is associated with a finite element. we do this
            // since not all lines may actually be used, in which
            // case we do not have to allocate any memory at all
            std::vector<bool> line_is_used(dof_handler.tria->n_raw_lines(),
                                           false);
            for (unsigned int line = 0; line < dof_handler.tria->n_raw_lines();
                 ++line)
              for (unsigned int fe = 0; fe < dof_handler.fe_collection.size();
                   ++fe)
                if (line_fe_association[fe][line] == true)
                  {
                    line_is_used[line] = true;
                    break;
                  }



            const unsigned int d = 1;
            const unsigned int l = 0;

            dof_handler.hp_object_fe_ptr[d].clear();
            dof_handler.hp_object_fe_indices[d].clear();
            dof_handler.object_dof_ptr[l][d].clear();
            dof_handler.object_dof_indices[l][d].clear();

            dof_handler.hp_object_fe_ptr[d].reserve(
              dof_handler.tria->n_raw_lines() + 1);

            unsigned int line_slots_needed = 0;
            unsigned int fe_slots_needed   = 0;

            for (unsigned int line = 0; line < dof_handler.tria->n_raw_lines();
                 ++line)
              {
                dof_handler.hp_object_fe_ptr[d].push_back(fe_slots_needed);

                if (line_is_used[line] == true)
                  {
                    for (unsigned int fe = 0;
                         fe < dof_handler.fe_collection.size();
                         ++fe)
                      if (line_fe_association[fe][line] == true)
                        {
                          ++fe_slots_needed;
                          line_slots_needed +=
                            dof_handler.get_fe(fe).n_dofs_per_line();
                        }
                  }
              }

            dof_handler.hp_object_fe_ptr[d].push_back(fe_slots_needed);

            // make sure that all entries have been set
            AssertDimension(dof_handler.hp_object_fe_ptr[d].size(),
                            dof_handler.tria->n_raw_lines() + 1);

            dof_handler.hp_object_fe_indices[d].reserve(fe_slots_needed);
            dof_handler.object_dof_ptr[l][d].reserve(fe_slots_needed + 1);

            dof_handler.object_dof_indices[l][d].reserve(line_slots_needed);

            for (unsigned int line = 0; line < dof_handler.tria->n_raw_lines();
                 ++line)
              if (line_is_used[line] == true)
                {
                  for (unsigned int fe = 0;
                       fe < dof_handler.fe_collection.size();
                       ++fe)
                    if (line_fe_association[fe][line] == true)
                      {
                        dof_handler.hp_object_fe_indices[d].push_back(fe);
                        dof_handler.object_dof_ptr[l][d].push_back(
                          dof_handler.object_dof_indices[l][d].size());

                        for (unsigned int i = 0;
                             i < dof_handler.get_fe(fe).n_dofs_per_line();
                             i++)
                          dof_handler.object_dof_indices[l][d].push_back(
                            numbers::invalid_dof_index);
                      }
                }

            dof_handler.object_dof_ptr[l][d].push_back(
              dof_handler.object_dof_indices[l][d].size());

            // make sure that all entries have been set
            AssertDimension(dof_handler.hp_object_fe_indices[d].size(),
                            fe_slots_needed);
            AssertDimension(dof_handler.object_dof_ptr[l][d].size(),
                            fe_slots_needed + 1);
            AssertDimension(dof_handler.object_dof_indices[l][d].size(),
                            line_slots_needed);
          }

          // Ensure that everything is done at this point.
          tasks.join_all();
        }



        /**
         * Given a DoFHandler object in hp-mode, make sure that the
         * active FE indices that a user has set for locally owned cells are
         * communicated to all other relevant cells as well.
         *
         * For parallel::shared::Triangulation objects,
         * this information is distributed on both ghost and artificial cells.
         *
         * In case a parallel::distributed::Triangulation is used,
         * indices are communicated only to ghost cells.
         */
        template <int dim, int spacedim>
        static void
        communicate_active_fe_indices(DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          if (const dealii::parallel::shared::Triangulation<dim, spacedim> *tr =
                dynamic_cast<
                  const dealii::parallel::shared::Triangulation<dim, spacedim>
                    *>(&dof_handler.get_triangulation()))
            {
              // we have a shared triangulation. in this case, every processor
              // knows about all cells, but every processor only has knowledge
              // about the active FE index on the cells it owns.
              //
              // we can create a complete set of active FE indices by letting
              // every processor create a vector of indices for all cells,
              // filling only those on the cells it owns and setting the indices
              // on the other cells to zero. then we add all of these vectors
              // up, and because every vector entry has exactly one processor
              // that owns it, the sum is correct
              std::vector<types::fe_index> active_fe_indices(
                tr->n_active_cells(), 0u);
              for (const auto &cell : dof_handler.active_cell_iterators())
                if (cell->is_locally_owned())
                  active_fe_indices[cell->active_cell_index()] =
                    cell->active_fe_index();

              Utilities::MPI::sum(active_fe_indices,
                                  tr->get_mpi_communicator(),
                                  active_fe_indices);

              // now go back and fill the active FE index on all other
              // cells. we would like to call cell->set_active_fe_index(),
              // but that function does not allow setting these indices on
              // non-locally_owned cells. so we have to work around the
              // issue a little bit by accessing the underlying data
              // structures directly
              for (const auto &cell : dof_handler.active_cell_iterators())
                if (!cell->is_locally_owned())
                  dof_handler
                    .hp_cell_active_fe_indices[cell->level()][cell->index()] =
                    active_fe_indices[cell->active_cell_index()];
            }
          else if (const dealii::parallel::
                     DistributedTriangulationBase<dim, spacedim> *tr =
                       dynamic_cast<
                         const dealii::parallel::
                           DistributedTriangulationBase<dim, spacedim> *>(
                         &dof_handler.get_triangulation()))
            {
              // For completely distributed meshes, use the function that is
              // able to move data from locally owned cells on one processor to
              // the corresponding ghost cells on others. To this end, we need
              // to have functions that can pack and unpack the data we want to
              // transport -- namely, the single unsigned int active_fe_index
              // objects
              auto pack =
                [](
                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                    &cell) -> types::fe_index {
                return cell->active_fe_index();
              };

              auto unpack =
                [&dof_handler](
                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                                       &cell,
                  const types::fe_index active_fe_index) -> void {
                // we would like to say
                //   cell->set_active_fe_index(active_fe_index);
                // but this is not allowed on cells that are not
                // locally owned, and we are on a ghost cell
                dof_handler
                  .hp_cell_active_fe_indices[cell->level()][cell->index()] =
                  active_fe_index;
              };

              GridTools::exchange_cell_data_to_ghosts<
                types::fe_index,
                DoFHandler<dim, spacedim>>(dof_handler, pack, unpack);
            }
          else
            {
              // a sequential triangulation. there is nothing we need to do here
              Assert(
                (dynamic_cast<
                   const dealii::parallel::TriangulationBase<dim, spacedim> *>(
                   &dof_handler.get_triangulation()) == nullptr),
                ExcInternalError());
            }
        }



        /**
         * Same as above, but for future FE indices.
         *
         * Given a DoFHandler object in hp-mode, make sure that the
         * future FE indices that a user has set for locally owned cells are
         * communicated to all other relevant cells as well.
         *
         * For parallel::shared::Triangulation objects,
         * this information is distributed on both ghost and artificial cells.
         *
         * In case a parallel::distributed::Triangulation is used,
         * indices are communicated only to ghost cells.
         */
        template <int dim, int spacedim>
        static void
        communicate_future_fe_indices(DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(
            dof_handler.hp_capability_enabled == true,
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          if (const dealii::parallel::shared::Triangulation<dim, spacedim> *tr =
                dynamic_cast<
                  const dealii::parallel::shared::Triangulation<dim, spacedim>
                    *>(&dof_handler.get_triangulation()))
            {
              std::vector<types::fe_index> future_fe_indices(
                tr->n_active_cells(), 0u);
              for (const auto &cell : dof_handler.active_cell_iterators() |
                                        IteratorFilters::LocallyOwnedCell())
                future_fe_indices[cell->active_cell_index()] =
                  dof_handler
                    .hp_cell_future_fe_indices[cell->level()][cell->index()];

              Utilities::MPI::sum(future_fe_indices,
                                  tr->get_mpi_communicator(),
                                  future_fe_indices);

              for (const auto &cell : dof_handler.active_cell_iterators())
                if (!cell->is_locally_owned())
                  dof_handler
                    .hp_cell_future_fe_indices[cell->level()][cell->index()] =
                    future_fe_indices[cell->active_cell_index()];
            }
          else if (const dealii::parallel::
                     DistributedTriangulationBase<dim, spacedim> *tr =
                       dynamic_cast<
                         const dealii::parallel::
                           DistributedTriangulationBase<dim, spacedim> *>(
                         &dof_handler.get_triangulation()))
            {
              auto pack =
                [&dof_handler](
                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                    &cell) -> types::fe_index {
                return dof_handler
                  .hp_cell_future_fe_indices[cell->level()][cell->index()];
              };

              auto unpack =
                [&dof_handler](
                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                                       &cell,
                  const types::fe_index future_fe_index) -> void {
                dof_handler
                  .hp_cell_future_fe_indices[cell->level()][cell->index()] =
                  future_fe_index;
              };

              GridTools::exchange_cell_data_to_ghosts<
                types::fe_index,
                DoFHandler<dim, spacedim>>(dof_handler, pack, unpack);
            }
          else
            {
              Assert(
                (dynamic_cast<
                   const dealii::parallel::TriangulationBase<dim, spacedim> *>(
                   &dof_handler.get_triangulation()) == nullptr),
                ExcInternalError());
            }
        }



        /**
         * Collect all finite element indices on cells that will be affected by
         * future refinement and coarsening. Further, prepare those indices to
         * be distributed on on the updated triangulation later.
         *
         * On cells to be refined, the active FE index will be inherited to
         * their children and thus will be stored as such.
         *
         * On cells to be coarsened, we choose the finite element on the parent
         * cell from those assigned to their children to be the one that is
         * dominated by all children. If none was found, we pick the most
         * dominant element in the whole collection that is dominated by all
         * children. See documentation of
         * hp::FECollection::find_dominated_fe_extended() for further
         * information.
         *
         * On cells intended for p-refinement or p-coarsening, those
         * active FE indices will be determined by the corresponding flags that
         * have been set on the relevant cells.
         */
        template <int dim, int spacedim>
        static void
        collect_fe_indices_on_cells_to_be_refined(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          const auto &fe_transfer = dof_handler.active_fe_index_transfer;

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                if (cell->refine_flag_set())
                  {
                    // Store the active FE index of each cell that will be
                    // refined to and distribute it later on its children.
                    // Pick their future index if flagged for p-refinement.
                    fe_transfer->refined_cells_fe_index.insert(
                      {cell, cell->future_fe_index()});
                  }
                else if (cell->coarsen_flag_set())
                  {
                    // From all cells that will be coarsened, determine their
                    // parent and calculate its proper active FE index, so that
                    // it can be set after refinement. But first, check if that
                    // particular cell has a parent at all.
                    Assert(cell->level() > 0, ExcInternalError());
                    const auto &parent = cell->parent();

                    // Check if the active FE index for the current cell has
                    // been determined already.
                    if (fe_transfer->coarsened_cells_fe_index.find(parent) ==
                        fe_transfer->coarsened_cells_fe_index.end())
                      {
                        // Find a suitable active FE index for the parent cell
                        // based on the 'least dominant finite element' of its
                        // children. Consider the childrens' hypothetical future
                        // index when they have been flagged for p-refinement.
                        if constexpr (library_build_mode ==
                                      LibraryBuildMode::debug)
                          {
                            for (const auto &child : parent->child_iterators())
                              Assert(child->is_active() &&
                                       child->coarsen_flag_set(),
                                     typename dealii::Triangulation<
                                       dim>::ExcInconsistentCoarseningFlags());
                          }

                        const types::fe_index fe_index = dealii::internal::hp::
                          DoFHandlerImplementation::Implementation::
                            dominated_future_fe_on_children<dim, spacedim>(
                              parent);

                        fe_transfer->coarsened_cells_fe_index.insert(
                          {parent, fe_index});
                      }
                  }
                else
                  {
                    // No h-refinement is scheduled for this cell.
                    // However, it may have p-refinement indicators, so we
                    // choose a new active FE index based on its flags.
                    if (cell->future_fe_index_set() == true)
                      fe_transfer->persisting_cells_fe_index.insert(
                        {cell, cell->future_fe_index()});
                  }
              }
        }



        /**
         * Distribute active finite element indices that have been previously
         * prepared in collect_fe_indices_on_cells_to_be_refined().
         */
        template <int dim, int spacedim>
        static void
        distribute_fe_indices_on_refined_cells(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          const auto &fe_transfer = dof_handler.active_fe_index_transfer;

          // Set active FE indices on persisting cells.
          for (const auto &persist : fe_transfer->persisting_cells_fe_index)
            {
              const auto &cell = persist.first;

              if (cell->is_locally_owned())
                {
                  Assert(cell->is_active(), ExcInternalError());
                  cell->set_active_fe_index(persist.second);
                }
            }

          // Distribute active FE indices from all refined cells on their
          // respective children.
          for (const auto &refine : fe_transfer->refined_cells_fe_index)
            {
              const auto &parent = refine.first;

              for (const auto &child : parent->child_iterators())
                if (child->is_locally_owned())
                  {
                    Assert(child->is_active(), ExcInternalError());
                    child->set_active_fe_index(refine.second);
                  }
            }

          // Set active FE indices on coarsened cells that have been determined
          // before the actual coarsening happened.
          for (const auto &coarsen : fe_transfer->coarsened_cells_fe_index)
            {
              const auto &cell = coarsen.first;

              if (cell->is_locally_owned())
                {
                  Assert(cell->is_active(), ExcInternalError());
                  cell->set_active_fe_index(coarsen.second);
                }
            }
        }


        /**
         * Coarsening strategy for the CellDataTransfer object responsible for
         * transferring the active FE index of each cell on
         * parallel::distributed::Triangulation objects that have been refined.
         *
         * A finite element index needs to be determined for the (not yet
         * active) parent cell from its (still active) children.  Out of the set
         * of elements previously assigned to the former children, we choose the
         * one dominated by all children for the parent cell.
         */
        template <int dim, int spacedim>
        static types::fe_index
        determine_fe_from_children(
          const typename Triangulation<dim, spacedim>::cell_iterator &,
          const std::vector<types::fe_index>            &children_fe_indices,
          const dealii::hp::FECollection<dim, spacedim> &fe_collection)
        {
          Assert(!children_fe_indices.empty(), ExcInternalError());

          // convert vector to set
          // TODO: Change set to types::fe_index
          const std::set<unsigned int> children_fe_indices_set(
            children_fe_indices.begin(), children_fe_indices.end());

          const types::fe_index dominated_fe_index =
            fe_collection.find_dominated_fe_extended(children_fe_indices_set,
                                                     /*codim=*/0);

          Assert(dominated_fe_index != numbers::invalid_fe_index,
                 ExcNoDominatedFiniteElementOnChildren());

          return dominated_fe_index;
        }


        /**
         * Return the index of the finite element from the entire
         * hp::FECollection that is dominated by those assigned as future finite
         * elements to the children of @p parent.
         *
         * See documentation in the header file for more information.
         */
        template <int dim, int spacedim>
        static types::fe_index
        dominated_future_fe_on_children(
          const typename DoFHandler<dim, spacedim>::cell_iterator &parent)
        {
          Assert(
            !parent->is_active(),
            ExcMessage(
              "You ask for information on children of this cell which is only "
              "available for active cells. This cell has no children."));

          const auto &dof_handler = parent->get_dof_handler();
          Assert(
            dof_handler.has_hp_capabilities(),
            (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

          // TODO: Change set to types::fe_index
          std::set<unsigned int> future_fe_indices_children;
          for (const auto &child : parent->child_iterators())
            {
              Assert(
                child->is_active(),
                ExcMessage(
                  "You ask for information on children of this cell which is only "
                  "available for active cells. One of its children is not active."));

              // Ghost siblings might occur on parallel Triangulation
              // objects. The public interface does not allow to access future
              // FE indices on ghost cells. However, we need this information
              // here and thus call the internal function that does not check
              // for cell ownership. This requires that future FE indices have
              // been communicated prior to calling this function.
              const types::fe_index future_fe_index_child =
                dealii::internal::DoFCellAccessorImplementation::
                  Implementation::future_fe_index<dim, spacedim, false>(*child);

              future_fe_indices_children.insert(future_fe_index_child);
            }
          Assert(!future_fe_indices_children.empty(), ExcInternalError());

          const types::fe_index future_fe_index =
            dof_handler.fe_collection.find_dominated_fe_extended(
              future_fe_indices_children,
              /*codim=*/0);

          Assert(future_fe_index != numbers::invalid_fe_index,
                 ExcNoDominatedFiniteElementOnChildren());

          return future_fe_index;
        }
      };



      /**
       * Public wrapper of Implementation::communicate_future_fe_indices().
       */
      template <int dim, int spacedim>
      void
      communicate_future_fe_indices(DoFHandler<dim, spacedim> &dof_handler)
      {
        Implementation::communicate_future_fe_indices<dim, spacedim>(
          dof_handler);
      }



      /**
       * Public wrapper of Implementation::dominated_future_fe_on_children().
       */
      template <int dim, int spacedim>
      unsigned int
      dominated_future_fe_on_children(
        const typename DoFHandler<dim, spacedim>::cell_iterator &parent)
      {
        return Implementation::dominated_future_fe_on_children<dim, spacedim>(
          parent);
      }
    } // namespace DoFHandlerImplementation
  }   // namespace hp
} // namespace internal

#ifndef DOXYGEN

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
DoFHandler<dim, spacedim>::DoFHandler()
  : hp_capability_enabled(true)
  , tria(nullptr, typeid(*this).name())
  , mg_faces(nullptr)
{}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
DoFHandler<dim, spacedim>::DoFHandler(const Triangulation<dim, spacedim> &tria)
  : DoFHandler()
{
  reinit(tria);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
DoFHandler<dim, spacedim>::~DoFHandler()
{
  // unsubscribe all attachments to signals of the underlying triangulation
  for (auto &connection : this->tria_listeners)
    connection.disconnect();
  this->tria_listeners.clear();

  for (auto &connection : this->tria_listeners_for_transfer)
    connection.disconnect();
  this->tria_listeners_for_transfer.clear();

  // release allocated memory
  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  DoFHandler<dim, spacedim>::clear();

  // also release the policy. this needs to happen before the
  // current object disappears because the policy objects
  // store references to the DoFhandler object they work on
  this->policy.reset();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::reinit(const Triangulation<dim, spacedim> &tria)
{
  //
  // call destructor
  //
  // remove association with old triangulation
  for (auto &connection : this->tria_listeners)
    connection.disconnect();
  this->tria_listeners.clear();

  for (auto &connection : this->tria_listeners_for_transfer)
    connection.disconnect();
  this->tria_listeners_for_transfer.clear();

  // release allocated memory and policy
  DoFHandler<dim, spacedim>::clear();
  this->policy.reset();

  // reset the finite element collection
  this->fe_collection = hp::FECollection<dim, spacedim>();

  //
  // call constructor
  //
  // establish connection to new triangulation
  this->tria = &tria;
  this->setup_policy();

  // start in hp-mode and let distribute_dofs toggle it if necessary
  hp_capability_enabled = true;
  this->connect_to_triangulation_signals();
  this->create_active_fe_table();
}

#endif
/*------------------------ Cell iterator functions ------------------------*/
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::cell_iterator
  DoFHandler<dim, spacedim>::begin(const unsigned int level) const
{
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().begin(level);
  if (cell == this->get_triangulation().end(level))
    return end(level);
  return cell_iterator(*cell, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::active_cell_iterator
  DoFHandler<dim, spacedim>::begin_active(const unsigned int level) const
{
  // level is checked in begin
  cell_iterator i = begin(level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::cell_iterator
  DoFHandler<dim, spacedim>::end() const
{
  return cell_iterator(&this->get_triangulation(), -1, -1, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::cell_iterator
  DoFHandler<dim, spacedim>::end(const unsigned int level) const
{
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().end(level);
  if (cell.state() != IteratorState::valid)
    return end();
  return cell_iterator(*cell, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::active_cell_iterator
  DoFHandler<dim, spacedim>::end_active(const unsigned int level) const
{
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().end_active(level);
  if (cell.state() != IteratorState::valid)
    return active_cell_iterator(end());
  return active_cell_iterator(*cell, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::level_cell_iterator
  DoFHandler<dim, spacedim>::begin_mg(const unsigned int level) const
{
  Assert(this->has_level_dofs(),
         ExcMessage("You can only iterate over mg "
                    "levels if mg dofs got distributed."));
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().begin(level);
  if (cell == this->get_triangulation().end(level))
    return end_mg(level);
  return level_cell_iterator(*cell, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::level_cell_iterator
  DoFHandler<dim, spacedim>::end_mg(const unsigned int level) const
{
  Assert(this->has_level_dofs(),
         ExcMessage("You can only iterate over mg "
                    "levels if mg dofs got distributed."));
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().end(level);
  if (cell.state() != IteratorState::valid)
    return end();
  return level_cell_iterator(*cell, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
typename DoFHandler<dim, spacedim>::level_cell_iterator
  DoFHandler<dim, spacedim>::end_mg() const
{
  return level_cell_iterator(&this->get_triangulation(), -1, -1, this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator> DoFHandler<
  dim,
  spacedim>::cell_iterators() const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>(
    begin(), end());
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
IteratorRange<typename DoFHandler<dim, spacedim>::
                active_cell_iterator> DoFHandler<dim, spacedim>::
  active_cell_iterators() const
{
  return IteratorRange<
    typename DoFHandler<dim, spacedim>::active_cell_iterator>(begin_active(),
                                                              end());
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
IteratorRange<
  typename DoFHandler<dim, spacedim>::
    level_cell_iterator> DoFHandler<dim, spacedim>::mg_cell_iterators() const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::level_cell_iterator>(
    begin_mg(), end_mg());
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator> DoFHandler<
  dim,
  spacedim>::cell_iterators_on_level(const unsigned int level) const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>(
    begin(level), end(level));
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
IteratorRange<typename DoFHandler<dim, spacedim>::
                active_cell_iterator> DoFHandler<dim, spacedim>::
  active_cell_iterators_on_level(const unsigned int level) const
{
  return IteratorRange<
    typename DoFHandler<dim, spacedim>::active_cell_iterator>(
    begin_active(level), end_active(level));
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
IteratorRange<typename DoFHandler<dim, spacedim>::
                level_cell_iterator> DoFHandler<dim, spacedim>::
  mg_cell_iterators_on_level(const unsigned int level) const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::level_cell_iterator>(
    begin_mg(level), end_mg(level));
}



//---------------------------------------------------------------------------



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
types::global_dof_index DoFHandler<dim, spacedim>::n_boundary_dofs() const
{
  Assert(!(dim == 2 && spacedim == 3) || hp_capability_enabled == false,
         ExcNotImplementedWithHP());

  Assert(this->fe_collection.size() > 0, ExcNoFESelected());

  std::unordered_set<types::global_dof_index> boundary_dofs;
  std::vector<types::global_dof_index>        dofs_on_face;
  dofs_on_face.reserve(this->get_fe_collection().max_dofs_per_face());

  const IndexSet &owned_dofs = locally_owned_dofs();

  // loop over all faces to check whether they are at a
  // boundary. note that we need not take special care of single
  // lines in 3d (using @p{cell->has_boundary_lines}), since we do
  // not support boundaries of dimension dim-2, and so every
  // boundary line is also part of a boundary face.
  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned() && cell->at_boundary())
      {
        for (const auto iface : cell->face_indices())
          {
            const auto face = cell->face(iface);
            if (face->at_boundary())
              {
                const unsigned int dofs_per_face =
                  cell->get_fe().n_dofs_per_face(iface);
                dofs_on_face.resize(dofs_per_face);

                face->get_dof_indices(dofs_on_face, cell->active_fe_index());
                for (unsigned int i = 0; i < dofs_per_face; ++i)
                  {
                    const unsigned int global_idof_index = dofs_on_face[i];
                    if (owned_dofs.is_element(global_idof_index))
                      {
                        boundary_dofs.insert(global_idof_index);
                      }
                  }
              }
          }
      }
  return boundary_dofs.size();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
types::global_dof_index DoFHandler<dim, spacedim>::n_boundary_dofs(
  const std::set<types::boundary_id> &boundary_ids) const
{
  Assert(!(dim == 2 && spacedim == 3) || hp_capability_enabled == false,
         ExcNotImplementedWithHP());

  Assert(this->fe_collection.size() > 0, ExcNoFESelected());
  Assert(boundary_ids.find(numbers::internal_face_boundary_id) ==
           boundary_ids.end(),
         ExcInvalidBoundaryIndicator());

  // same as above, but with additional checks for set of boundary
  // indicators
  std::unordered_set<types::global_dof_index> boundary_dofs;
  std::vector<types::global_dof_index>        dofs_on_face;
  dofs_on_face.reserve(this->get_fe_collection().max_dofs_per_face());

  const IndexSet &owned_dofs = locally_owned_dofs();

  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned() && cell->at_boundary())
      {
        for (const auto iface : cell->face_indices())
          {
            const auto         face        = cell->face(iface);
            const unsigned int boundary_id = face->boundary_id();
            if (face->at_boundary() &&
                (boundary_ids.find(boundary_id) != boundary_ids.end()))
              {
                const unsigned int dofs_per_face =
                  cell->get_fe().n_dofs_per_face(iface);
                dofs_on_face.resize(dofs_per_face);

                face->get_dof_indices(dofs_on_face, cell->active_fe_index());
                for (unsigned int i = 0; i < dofs_per_face; ++i)
                  {
                    const unsigned int global_idof_index = dofs_on_face[i];
                    if (owned_dofs.is_element(global_idof_index))
                      {
                        boundary_dofs.insert(global_idof_index);
                      }
                  }
              }
          }
      }
  return boundary_dofs.size();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
std::size_t DoFHandler<dim, spacedim>::memory_consumption() const
{
  std::size_t mem = MemoryConsumption::memory_consumption(this->tria) +
                    MemoryConsumption::memory_consumption(this->fe_collection) +
                    MemoryConsumption::memory_consumption(this->number_cache);

  mem += MemoryConsumption::memory_consumption(object_dof_indices) +
         MemoryConsumption::memory_consumption(object_dof_ptr) +
         MemoryConsumption::memory_consumption(hp_object_fe_indices) +
         MemoryConsumption::memory_consumption(hp_object_fe_ptr) +
         MemoryConsumption::memory_consumption(hp_cell_active_fe_indices) +
         MemoryConsumption::memory_consumption(hp_cell_future_fe_indices);


  if (hp_capability_enabled)
    {
      // nothing to add
    }
  else
    {
      // collect size of multigrid data structures

      mem += MemoryConsumption::memory_consumption(this->block_info_object);

      for (unsigned int level = 0; level < this->mg_levels.size(); ++level)
        mem += this->mg_levels[level]->memory_consumption();

      if (this->mg_faces != nullptr)
        mem += MemoryConsumption::memory_consumption(*this->mg_faces);

      for (unsigned int i = 0; i < this->mg_vertex_dofs.size(); ++i)
        mem += sizeof(MGVertexDoFs) +
               (1 + this->mg_vertex_dofs[i].get_finest_level() -
                this->mg_vertex_dofs[i].get_coarsest_level()) *
                 sizeof(types::global_dof_index);
    }

  return mem;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::distribute_dofs(
  const FiniteElement<dim, spacedim> &fe)
{
  this->distribute_dofs(hp::FECollection<dim, spacedim>(fe));
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::distribute_dofs(
  const hp::FECollection<dim, spacedim> &ff)
{
  Assert(this->tria != nullptr,
         ExcMessage(
           "You need to set the Triangulation in the DoFHandler using reinit() "
           "or in the constructor before you can distribute DoFs."));
  Assert(this->tria->n_levels() > 0,
         ExcMessage("The Triangulation you are using is empty!"));

  // verify size of provided FE collection
  Assert(ff.size() > 0, ExcMessage("The given hp::FECollection is empty!"));
  Assert((ff.size() <= std::numeric_limits<types::fe_index>::max()) &&
           (ff.size() != numbers::invalid_fe_index),
         ExcMessage("The given hp::FECollection contains more finite elements "
                    "than the DoFHandler can cover with active FE indices."));

  if constexpr (running_in_debug_mode())
    {
      // make sure that the provided FE collection is large enough to
      // cover all FE indices presently in use on the mesh
      if ((hp_cell_active_fe_indices.size() > 0) &&
          (hp_cell_future_fe_indices.size() > 0))
        {
          Assert(hp_capability_enabled, ExcInternalError());

          for (const auto &cell : this->active_cell_iterators() |
                                    IteratorFilters::LocallyOwnedCell())
            {
              Assert(cell->active_fe_index() < ff.size(),
                     ExcInvalidFEIndex(cell->active_fe_index(), ff.size()));
              Assert(cell->future_fe_index() < ff.size(),
                     ExcInvalidFEIndex(cell->future_fe_index(), ff.size()));
            }
        }
    }

  //
  // register the new finite element collection
  //
  // don't create a new object if the one we have is identical
  if (this->fe_collection != ff)
    {
      this->fe_collection = hp::FECollection<dim, spacedim>(ff);

      const bool contains_multiple_fes = (this->fe_collection.size() > 1);

      // disable hp-mode if only a single finite element has been registered
      if (hp_capability_enabled && !contains_multiple_fes)
        {
          hp_capability_enabled = false;

          // unsubscribe connections to signals that are only relevant for
          // hp-mode, since we only have a single element here
          for (auto &connection : this->tria_listeners_for_transfer)
            connection.disconnect();
          this->tria_listeners_for_transfer.clear();

          // release active and future finite element tables
          this->hp_cell_active_fe_indices.clear();
          this->hp_cell_active_fe_indices.shrink_to_fit();
          this->hp_cell_future_fe_indices.clear();
          this->hp_cell_future_fe_indices.shrink_to_fit();
        }

      // re-enabling hp-mode is not permitted since the active and future FE
      // tables are no longer available
      AssertThrow(
        hp_capability_enabled || !contains_multiple_fes,
        ExcMessage(
          "You cannot re-enable hp-capabilities after you registered a single "
          "finite element. Please call reinit() or create a new DoFHandler "
          "object instead."));
    }

  //
  // enumerate all degrees of freedom
  //
  if (hp_capability_enabled)
    {
      // make sure every processor knows the active FE indices
      // on both its own cells and all ghost cells
      internal::hp::DoFHandlerImplementation::Implementation::
        communicate_active_fe_indices(*this);
    }

  {
    // We would like to enumerate all dofs for shared::Triangulations. If an
    // underlying shared::Tria allows artificial cells, we need to restore the
    // true cell owners temporarily.
    // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
    // the current set of subdomain ids, set subdomain ids to the "true" owner
    // of each cell upon construction of the TemporarilyRestoreSubdomainIds
    // object, and later restore these flags when it is destroyed.
    const internal::parallel::shared::TemporarilyRestoreSubdomainIds<dim,
                                                                     spacedim>
      subdomain_modifier(this->get_triangulation());

    // Adjust size of levels to the triangulation. Note that we still have to
    // allocate space for all degrees of freedom on this mesh (including ghost
    // and cells that are entirely stored on different processors), though we
    // may not assign numbers to some of them (i.e. they will remain at
    // invalid_dof_index). We need to allocate the space because we will want
    // to be able to query the dof_indices on each cell, and simply be told
    // that we don't know them on some cell (i.e. get back invalid_dof_index)
    if (hp_capability_enabled)
      internal::hp::DoFHandlerImplementation::Implementation::reserve_space(
        *this);
    else
      internal::DoFHandlerImplementation::Implementation::reserve_space(*this);
  }

  // hand the actual work over to the policy
  this->number_cache = this->policy->distribute_dofs();

  // do some housekeeping: compress indices
  // if(hp_capability_enabled)
  //   {
  //     Threads::TaskGroup<> tg;
  //     for (int level = this->levels_hp.size() - 1; level >= 0; --level)
  //       tg += Threads::new_task(
  //         &dealii::internal::hp::DoFLevel::compress_data<dim, spacedim>,
  //         *this->levels_hp[level],
  //         this->fe_collection);
  //     tg.join_all();
  //   }

  // Initialize the block info object only if this is a sequential
  // triangulation. It doesn't work correctly yet if it is parallel and has not
  // yet been implemented for hp-mode.
  if (!hp_capability_enabled &&
      dynamic_cast<const parallel::DistributedTriangulationBase<dim, spacedim>
                     *>(&*this->tria) == nullptr)
    this->block_info_object.initialize(*this, false, true);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::distribute_mg_dofs()
{
  AssertThrow(hp_capability_enabled == false, ExcNotImplementedWithHP());

  Assert(
    this->object_dof_indices.size() > 0,
    ExcMessage(
      "Distribute active DoFs using distribute_dofs() before calling distribute_mg_dofs()."));

  Assert(
    ((this->tria->get_mesh_smoothing() &
      Triangulation<dim, spacedim>::limit_level_difference_at_vertices) !=
     Triangulation<dim, spacedim>::none),
    ExcMessage(
      "The mesh smoothing requirement 'limit_level_difference_at_vertices' has to be set for using multigrid!"));

  this->clear_mg_space();

  internal::DoFHandlerImplementation::Implementation::reserve_space_mg(*this);
  this->mg_number_cache = this->policy->distribute_mg_dofs();

  // initialize the block info object only if this is a sequential
  // triangulation. it doesn't work correctly yet if it is parallel
  if (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &*this->tria) == nullptr)
    this->block_info_object.initialize(*this, true, false);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::initialize_local_block_info()
{
  AssertThrow(hp_capability_enabled == false, ExcNotImplementedWithHP());

  this->block_info_object.initialize_local(*this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::setup_policy()
{
  // decide whether we need a sequential or a parallel distributed policy
  if (dynamic_cast<const dealii::parallel::shared::Triangulation<dim, spacedim>
                     *>(&this->get_triangulation()) != nullptr)
    this->policy = std::make_unique<internal::DoFHandlerImplementation::Policy::
                                      ParallelShared<dim, spacedim>>(*this);
  else if (dynamic_cast<
             const dealii::parallel::DistributedTriangulationBase<dim, spacedim>
               *>(&this->get_triangulation()) == nullptr)
    this->policy = std::make_unique<
      internal::DoFHandlerImplementation::Policy::Sequential<dim, spacedim>>(
      *this);
  else
    this->policy =
      std::make_unique<internal::DoFHandlerImplementation::Policy::
                         ParallelDistributed<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::clear()
{
  // release memory
  this->clear_space();
  this->clear_mg_space();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::clear_space()
{
  object_dof_indices.clear();

  object_dof_ptr.clear();

  this->number_cache.clear();

  this->hp_cell_active_fe_indices.clear();
  this->hp_cell_future_fe_indices.clear();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::clear_mg_space()
{
  this->mg_levels.clear();
  this->mg_faces.reset();

  std::vector<MGVertexDoFs> tmp;

  std::swap(this->mg_vertex_dofs, tmp);

  this->mg_number_cache.clear();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::renumber_dofs(
  const std::vector<types::global_dof_index> &new_numbers)
{
  if (hp_capability_enabled)
    {
      Assert(this->hp_cell_future_fe_indices.size() > 0,
             ExcMessage(
               "You need to distribute DoFs before you can renumber them."));

      AssertDimension(new_numbers.size(), this->n_locally_owned_dofs());

      if constexpr (running_in_debug_mode())
        {
          // assert that the new indices are consecutively numbered if we are
          // working on a single processor. this doesn't need to
          // hold in the case of a parallel mesh since we map the interval
          // [0...n_dofs()) into itself but only globally, not on each processor
          if (this->n_locally_owned_dofs() == this->n_dofs())
            {
              std::vector<types::global_dof_index> tmp(new_numbers);
              std::sort(tmp.begin(), tmp.end());
              std::vector<types::global_dof_index>::const_iterator p =
                tmp.begin();
              types::global_dof_index i = 0;
              for (; p != tmp.end(); ++p, ++i)
                Assert(*p == i, ExcNewNumbersNotConsecutive(i));
            }
          else
            for (const auto new_number : new_numbers)
              Assert(
                new_number < this->n_dofs(),
                ExcMessage(
                  "New DoF index is not less than the total number of dofs."));
        }

      // uncompress the internal storage scheme of dofs on cells so that
      // we can access dofs in turns. uncompress in parallel, starting
      // with the most expensive levels (the highest ones)
      //{
      //  Threads::TaskGroup<> tg;
      //  for (int level = this->levels_hp.size() - 1; level >= 0; --level)
      //    tg += Threads::new_task(
      //      &dealii::internal::hp::DoFLevel::uncompress_data<dim, spacedim>,
      //      *this->levels_hp[level],
      //      this->fe_collection);
      //  tg.join_all();
      //}

      // do the renumbering
      this->number_cache = this->policy->renumber_dofs(new_numbers);

      // now re-compress the dof indices
      //{
      //  Threads::TaskGroup<> tg;
      //  for (int level = this->levels_hp.size() - 1; level >= 0; --level)
      //    tg += Threads::new_task(
      //      &dealii::internal::hp::DoFLevel::compress_data<dim, spacedim>,
      //     *this->levels_hp[level],
      //      this->fe_collection);
      //  tg.join_all();
      //}
    }
  else
    {
      Assert(this->object_dof_indices.size() > 0,
             ExcMessage(
               "You need to distribute DoFs before you can renumber them."));

      if constexpr (running_in_debug_mode())
        {
          if (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>
                             *>(&*this->tria) != nullptr)
            {
              Assert(new_numbers.size() == this->n_dofs() ||
                       new_numbers.size() == this->n_locally_owned_dofs(),
                     ExcMessage("Incorrect size of the input array."));
            }
          else if (dynamic_cast<
                     const parallel::DistributedTriangulationBase<dim, spacedim>
                       *>(&*this->tria) != nullptr)
            {
              AssertDimension(new_numbers.size(), this->n_locally_owned_dofs());
            }
          else
            {
              AssertDimension(new_numbers.size(), this->n_dofs());
            }

          // assert that the new indices are consecutively numbered if we are
          // working on a single processor. this doesn't need to
          // hold in the case of a parallel mesh since we map the interval
          // [0...n_dofs()) into itself but only globally, not on each processor
          if (this->n_locally_owned_dofs() == this->n_dofs())
            {
              std::vector<types::global_dof_index> tmp(new_numbers);
              std::sort(tmp.begin(), tmp.end());
              std::vector<types::global_dof_index>::const_iterator p =
                tmp.begin();
              types::global_dof_index i = 0;
              for (; p != tmp.end(); ++p, ++i)
                Assert(*p == i, ExcNewNumbersNotConsecutive(i));
            }
          else
            for (const auto new_number : new_numbers)
              Assert(
                new_number < this->n_dofs(),
                ExcMessage(
                  "New DoF index is not less than the total number of dofs."));
        }

      this->number_cache = this->policy->renumber_dofs(new_numbers);
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::renumber_dofs(
  const unsigned int                          level,
  const std::vector<types::global_dof_index> &new_numbers)
{
  AssertThrow(hp_capability_enabled == false, ExcNotImplementedWithHP());

  Assert(
    this->mg_levels.size() > 0 && this->object_dof_indices.size() > 0,
    ExcMessage(
      "You need to distribute active and level DoFs before you can renumber level DoFs."));
  AssertIndexRange(level, this->get_triangulation().n_global_levels());
  AssertDimension(new_numbers.size(),
                  this->locally_owned_mg_dofs(level).n_elements());

  if constexpr (running_in_debug_mode())
    {
      // assert that the new indices are consecutively numbered if we are
      // working on a single processor. this doesn't need to hold in the case of
      // a parallel mesh since we map the interval [0...n_dofs(level)) into
      // itself but only globally, not on each processor
      if (this->n_locally_owned_dofs() == this->n_dofs())
        {
          std::vector<types::global_dof_index> tmp(new_numbers);
          std::sort(tmp.begin(), tmp.end());
          std::vector<types::global_dof_index>::const_iterator p = tmp.begin();
          types::global_dof_index                              i = 0;
          for (; p != tmp.end(); ++p, ++i)
            Assert(*p == i, ExcNewNumbersNotConsecutive(i));
        }
      else
        for (const auto new_number : new_numbers)
          Assert(new_number < this->n_dofs(level),
                 ExcMessage(
                   "New DoF index is not less than the total number of dofs."));
    }

  this->mg_number_cache[level] =
    this->policy->renumber_mg_dofs(level, new_numbers);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
unsigned int DoFHandler<dim, spacedim>::max_couplings_between_boundary_dofs()
  const
{
  Assert(this->fe_collection.size() > 0, ExcNoFESelected());

  switch (dim)
    {
      case 1:
        return this->fe_collection.max_dofs_per_vertex();
      case 2:
        return (3 * this->fe_collection.max_dofs_per_vertex() +
                2 * this->fe_collection.max_dofs_per_line());
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
        return (19 * this->fe_collection.max_dofs_per_vertex() +
                28 * this->fe_collection.max_dofs_per_line() +
                8 * this->fe_collection.max_dofs_per_quad());
      default:
        DEAL_II_NOT_IMPLEMENTED();
        return 0;
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
unsigned int DoFHandler<dim, spacedim>::max_couplings_between_dofs() const
{
  Assert(this->fe_collection.size() > 0, ExcNoFESelected());
  return internal::DoFHandlerImplementation::Implementation::
    max_couplings_between_dofs(*this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::set_active_fe_indices(
  const std::vector<types::fe_index> &active_fe_indices)
{
  Assert(active_fe_indices.size() == this->get_triangulation().n_active_cells(),
         ExcDimensionMismatch(active_fe_indices.size(),
                              this->get_triangulation().n_active_cells()));

  this->create_active_fe_table();
  // we could set the values directly, since they are stored as
  // protected data of this object, but for simplicity we use the
  // cell-wise access. this way we also have to pass some debug-mode
  // tests which we would have to duplicate ourselves otherwise
  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_active_fe_index(active_fe_indices[cell->active_cell_index()]);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::set_active_fe_indices(
  const std::vector<unsigned int> &active_fe_indices)
{
  set_active_fe_indices(std::vector<types::fe_index>(active_fe_indices.begin(),
                                                     active_fe_indices.end()));
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
std::vector<types::fe_index> DoFHandler<dim, spacedim>::get_active_fe_indices()
  const
{
  std::vector<types::fe_index> active_fe_indices(
    this->get_triangulation().n_active_cells(), numbers::invalid_fe_index);

  // we could try to extract the values directly, since they are
  // stored as protected data of this object, but for simplicity we
  // use the cell-wise access.
  for (const auto &cell : this->active_cell_iterators())
    if (!cell->is_artificial())
      active_fe_indices[cell->active_cell_index()] = cell->active_fe_index();

  return active_fe_indices;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::get_active_fe_indices(
  std::vector<unsigned int> &active_fe_indices) const
{
  const std::vector<types::fe_index> indices = get_active_fe_indices();

  active_fe_indices.assign(indices.begin(), indices.end());
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::set_future_fe_indices(
  const std::vector<types::fe_index> &future_fe_indices)
{
  Assert(future_fe_indices.size() == this->get_triangulation().n_active_cells(),
         ExcDimensionMismatch(future_fe_indices.size(),
                              this->get_triangulation().n_active_cells()));

  this->create_active_fe_table();
  // we could set the values directly, since they are stored as
  // protected data of this object, but for simplicity we use the
  // cell-wise access. this way we also have to pass some debug-mode
  // tests which we would have to duplicate ourselves otherwise
  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned() &&
        future_fe_indices[cell->active_cell_index()] !=
          numbers::invalid_fe_index)
      cell->set_future_fe_index(future_fe_indices[cell->active_cell_index()]);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
std::vector<types::fe_index> DoFHandler<dim, spacedim>::get_future_fe_indices()
  const
{
  std::vector<types::fe_index> future_fe_indices(
    this->get_triangulation().n_active_cells(), numbers::invalid_fe_index);

  // we could try to extract the values directly, since they are
  // stored as protected data of this object, but for simplicity we
  // use the cell-wise access.
  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned() && cell->future_fe_index_set())
      future_fe_indices[cell->active_cell_index()] = cell->future_fe_index();

  return future_fe_indices;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::connect_to_triangulation_signals()
{
  // make sure this is called during initialization in hp-mode
  Assert(hp_capability_enabled, ExcOnlyAvailableWithHP());

  // connect functions to signals of the underlying triangulation
  this->tria_listeners.push_back(this->tria->signals.create.connect(
    [this]() { this->reinit(*(this->tria)); }));
  this->tria_listeners.push_back(
    this->tria->signals.clear.connect([this]() { this->clear(); }));

  // attach corresponding callback functions dealing with the transfer of
  // active FE indices depending on the type of triangulation
  if (dynamic_cast<
        const dealii::parallel::fullydistributed::Triangulation<dim, spacedim>
          *>(&this->get_triangulation()))
    {
      // no transfer of active FE indices for this class
    }
  else if (dynamic_cast<
             const dealii::parallel::distributed::Triangulation<dim, spacedim>
               *>(&this->get_triangulation()))
    {
      // repartitioning signals
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.pre_distributed_repartition.connect([this]() {
          internal::hp::DoFHandlerImplementation::Implementation::
            ensure_absence_of_future_fe_indices<dim, spacedim>(*this);
        }));
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.pre_distributed_repartition.connect(
          [this]() { this->pre_distributed_transfer_action(); }));
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_distributed_repartition.connect(
          [this]() { this->post_distributed_transfer_action(); }));

      // refinement signals
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_p4est_refinement.connect(
          [this]() { this->pre_distributed_transfer_action(); }));
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_distributed_refinement.connect(
          [this]() { this->post_distributed_transfer_action(); }));

      // serialization signals
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_distributed_save.connect(
          [this]() { this->active_fe_index_transfer.reset(); }));
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_distributed_load.connect(
          [this]() { this->update_active_fe_table(); }));
    }
  else if (dynamic_cast<
             const dealii::parallel::shared::Triangulation<dim, spacedim> *>(
             &this->get_triangulation()) != nullptr)
    {
      // partitioning signals
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.pre_partition.connect([this]() {
          internal::hp::DoFHandlerImplementation::Implementation::
            ensure_absence_of_future_fe_indices(*this);
        }));

      // refinement signals
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.pre_refinement.connect(
          [this]() { this->pre_transfer_action(); }));
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_refinement.connect(
          [this]() { this->post_transfer_action(); }));
    }
  else
    {
      // refinement signals
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.pre_refinement.connect(
          [this]() { this->pre_transfer_action(); }));
      this->tria_listeners_for_transfer.push_back(
        this->tria->signals.post_refinement.connect(
          [this]() { this->post_transfer_action(); }));
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::create_active_fe_table()
{
  AssertThrow(hp_capability_enabled == true, ExcOnlyAvailableWithHP());


  // Create sufficiently many hp::DoFLevels.
  //  while (this->levels_hp.size() < this->tria->n_levels())
  //    this->levels_hp.emplace_back(new dealii::internal::hp::DoFLevel);

  this->hp_cell_active_fe_indices.resize(this->tria->n_levels());
  this->hp_cell_future_fe_indices.resize(this->tria->n_levels());

  // then make sure that on each level we have the appropriate size
  // of active FE indices; preset them to zero, i.e. the default FE
  for (unsigned int level = 0; level < this->hp_cell_future_fe_indices.size();
       ++level)
    {
      if (this->hp_cell_active_fe_indices[level].empty() &&
          this->hp_cell_future_fe_indices[level].empty())
        {
          this->hp_cell_active_fe_indices[level].resize(
            this->tria->n_raw_cells(level), 0);
          this->hp_cell_future_fe_indices[level].resize(
            this->tria->n_raw_cells(level), numbers::invalid_fe_index);
        }
      else
        {
          // Either the active FE indices have size zero because
          // they were just created, or the correct size. Other
          // sizes indicate that something went wrong.
          Assert(this->hp_cell_active_fe_indices[level].size() ==
                     this->tria->n_raw_cells(level) &&
                   this->hp_cell_future_fe_indices[level].size() ==
                     this->tria->n_raw_cells(level),
                 ExcInternalError());
        }

      // it may be that the previous table was compressed; in that
      // case, restore the correct active FE index. the fact that
      // this no longer matches the indices in the table is of no
      // importance because the current function is called at a
      // point where we have to recreate the dof_indices tables in
      // the levels anyway
      // this->levels_hp[level]->normalize_active_fe_indices();
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::update_active_fe_table()
{
  //  // Normally only one level is added, but if this Triangulation
  //  // is created by copy_triangulation, it can be more than one level.
  //  while (this->levels_hp.size() < this->tria->n_levels())
  //    this->levels_hp.emplace_back(new dealii::internal::hp::DoFLevel);
  //
  //  // Coarsening can lead to the loss of levels. Hence remove them.
  //  while (this->levels_hp.size() > this->tria->n_levels())
  //    {
  //      // drop the last element. that also releases the memory pointed to
  //      this->levels_hp.pop_back();
  //    }

  this->hp_cell_active_fe_indices.resize(this->tria->n_levels());
  this->hp_cell_active_fe_indices.shrink_to_fit();

  this->hp_cell_future_fe_indices.resize(this->tria->n_levels());
  this->hp_cell_future_fe_indices.shrink_to_fit();

  for (unsigned int i = 0; i < this->hp_cell_future_fe_indices.size(); ++i)
    {
      // Resize active FE indices vectors. Use zero indicator to extend.
      this->hp_cell_active_fe_indices[i].resize(this->tria->n_raw_cells(i), 0);

      // Resize future FE indices vectors. Make sure that all
      // future FE indices have been cleared after refinement happened.
      //
      // We have used future FE indices to update all active FE indices
      // before refinement happened, thus we are safe to clear them now.
      this->hp_cell_future_fe_indices[i].assign(this->tria->n_raw_cells(i),
                                                numbers::invalid_fe_index);
    }
}


template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::pre_transfer_action()
{
  Assert(this->active_fe_index_transfer == nullptr, ExcInternalError());

  this->active_fe_index_transfer = std::make_unique<ActiveFEIndexTransfer>();

  internal::hp::DoFHandlerImplementation::Implementation::
    communicate_future_fe_indices(*this);

  dealii::internal::hp::DoFHandlerImplementation::Implementation::
    collect_fe_indices_on_cells_to_be_refined(*this);
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::pre_distributed_transfer_action()
{
#  ifndef DEAL_II_WITH_P4EST
  Assert(false,
         ExcMessage(
           "You are attempting to use a functionality that is only available "
           "if deal.II was configured to use p4est, but cmake did not find a "
           "valid p4est library."));
#  else
  // the implementation below requires a p:d:T currently
  Assert(
    (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
       &this->get_triangulation()) != nullptr),
    ExcNotImplemented());

  Assert(active_fe_index_transfer == nullptr, ExcInternalError());

  active_fe_index_transfer = std::make_unique<ActiveFEIndexTransfer>();

  // If we work on a p::d::Triangulation, we have to transfer all
  // active FE indices since ownership of cells may change. We will
  // use our p::d::CellDataTransfer member to achieve this. Further,
  // we prepare the values in such a way that they will correspond to
  // the active FE indices on the new mesh.

  // Gather all current future FE indices.
  active_fe_index_transfer->active_fe_indices.resize(
    get_triangulation().n_active_cells(), numbers::invalid_fe_index);

  // Collect future FE indices on locally owned and ghost cells.
  // The public interface does not allow to access future FE indices
  // on ghost cells. However, we need this information here and thus
  // call the internal function that does not check for cell ownership.
  internal::hp::DoFHandlerImplementation::Implementation::
    communicate_future_fe_indices(*this);

  for (const auto &cell : active_cell_iterators())
    if (cell->is_artificial() == false)
      active_fe_index_transfer->active_fe_indices[cell->active_cell_index()] =
        dealii::internal::DoFCellAccessorImplementation::Implementation::
          future_fe_index<dim, spacedim, false>(*cell);

  // Create transfer object and attach to it.
  const auto *distributed_tria =
    dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
      &this->get_triangulation());

  active_fe_index_transfer->cell_data_transfer = std::make_unique<
    parallel::distributed::
      CellDataTransfer<dim, spacedim, std::vector<types::fe_index>>>(
    *distributed_tria,
    /*transfer_variable_size_data=*/false,
    /*refinement_strategy=*/
    &dealii::AdaptationStrategies::Refinement::
      preserve<dim, spacedim, types::fe_index>,
    /*coarsening_strategy=*/
    [this](const typename Triangulation<dim, spacedim>::cell_iterator &parent,
           const std::vector<types::fe_index> &children_fe_indices)
      -> types::fe_index {
      return dealii::internal::hp::DoFHandlerImplementation::Implementation::
        determine_fe_from_children<dim, spacedim>(parent,
                                                  children_fe_indices,
                                                  fe_collection);
    });

  active_fe_index_transfer->cell_data_transfer
    ->prepare_for_coarsening_and_refinement(
      active_fe_index_transfer->active_fe_indices);
#  endif
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::post_transfer_action()
{
  update_active_fe_table();

  Assert(this->active_fe_index_transfer != nullptr, ExcInternalError());

  dealii::internal::hp::DoFHandlerImplementation::Implementation::
    distribute_fe_indices_on_refined_cells(*this);

  // We have to distribute the information about active FE indices
  // of all cells (including the artificial ones) on all processors,
  // if a parallel::shared::Triangulation has been used.
  dealii::internal::hp::DoFHandlerImplementation::Implementation::
    communicate_active_fe_indices(*this);

  // Free memory.
  this->active_fe_index_transfer.reset();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::post_distributed_transfer_action()
{
#  ifndef DEAL_II_WITH_P4EST
  DEAL_II_ASSERT_UNREACHABLE();
#  else
  update_active_fe_table();

  Assert(this->active_fe_index_transfer != nullptr, ExcInternalError());

  // Unpack active FE indices.
  this->active_fe_index_transfer->active_fe_indices.resize(
    this->get_triangulation().n_active_cells(), numbers::invalid_fe_index);
  this->active_fe_index_transfer->cell_data_transfer->unpack(
    this->active_fe_index_transfer->active_fe_indices);

  // Update all locally owned active FE indices.
  this->set_active_fe_indices(
    this->active_fe_index_transfer->active_fe_indices);

  // Update active FE indices on ghost cells.
  dealii::internal::hp::DoFHandlerImplementation::Implementation::
    communicate_active_fe_indices(*this);

  // Free memory.
  this->active_fe_index_transfer.reset();
#  endif
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::prepare_for_serialization_of_active_fe_indices()
{
#  ifndef DEAL_II_WITH_P4EST
  Assert(false,
         ExcMessage(
           "You are attempting to use a functionality that is only available "
           "if deal.II was configured to use p4est, but cmake did not find a "
           "valid p4est library."));
#  else
  // the implementation below requires a p:d:T currently
  Assert(
    (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
       &this->get_triangulation()) != nullptr),
    ExcNotImplemented());

  Assert(active_fe_index_transfer == nullptr, ExcInternalError());

  active_fe_index_transfer = std::make_unique<ActiveFEIndexTransfer>();

  // Create transfer object and attach to it.
  const auto *distributed_tria =
    dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
      &this->get_triangulation());

  active_fe_index_transfer->cell_data_transfer = std::make_unique<
    parallel::distributed::
      CellDataTransfer<dim, spacedim, std::vector<types::fe_index>>>(
    *distributed_tria,
    /*transfer_variable_size_data=*/false,
    /*refinement_strategy=*/
    &dealii::AdaptationStrategies::Refinement::
      preserve<dim, spacedim, types::fe_index>,
    /*coarsening_strategy=*/
    [this](const typename Triangulation<dim, spacedim>::cell_iterator &parent,
           const std::vector<types::fe_index> &children_fe_indices)
      -> types::fe_index {
      return dealii::internal::hp::DoFHandlerImplementation::Implementation::
        determine_fe_from_children<dim, spacedim>(parent,
                                                  children_fe_indices,
                                                  fe_collection);
    });

  // If we work on a p::d::Triangulation, we have to transfer all
  // active FE indices since ownership of cells may change.

  // Gather all current active FE indices
  active_fe_index_transfer->active_fe_indices = get_active_fe_indices();

  // Attach to transfer object
  active_fe_index_transfer->cell_data_transfer->prepare_for_serialization(
    active_fe_index_transfer->active_fe_indices);
#  endif
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::deserialize_active_fe_indices()
{
#  ifndef DEAL_II_WITH_P4EST
  Assert(false,
         ExcMessage(
           "You are attempting to use a functionality that is only available "
           "if deal.II was configured to use p4est, but cmake did not find a "
           "valid p4est library."));
#  else
  // the implementation below requires a p:d:T currently
  Assert(
    (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
       &this->get_triangulation()) != nullptr),
    ExcNotImplemented());

  Assert(active_fe_index_transfer == nullptr, ExcInternalError());

  active_fe_index_transfer = std::make_unique<ActiveFEIndexTransfer>();

  // Create transfer object and attach to it.
  const auto *distributed_tria =
    dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
      &this->get_triangulation());

  active_fe_index_transfer->cell_data_transfer = std::make_unique<
    parallel::distributed::
      CellDataTransfer<dim, spacedim, std::vector<types::fe_index>>>(
    *distributed_tria,
    /*transfer_variable_size_data=*/false,
    /*refinement_strategy=*/
    &dealii::AdaptationStrategies::Refinement::
      preserve<dim, spacedim, types::fe_index>,
    /*coarsening_strategy=*/
    [this](const typename Triangulation<dim, spacedim>::cell_iterator &parent,
           const std::vector<types::fe_index> &children_fe_indices)
      -> types::fe_index {
      return dealii::internal::hp::DoFHandlerImplementation::Implementation::
        determine_fe_from_children<dim, spacedim>(parent,
                                                  children_fe_indices,
                                                  fe_collection);
    });

  // Unpack active FE indices.
  active_fe_index_transfer->active_fe_indices.resize(
    get_triangulation().n_active_cells(), numbers::invalid_fe_index);
  active_fe_index_transfer->cell_data_transfer->deserialize(
    active_fe_index_transfer->active_fe_indices);

  // Update all locally owned active FE indices.
  set_active_fe_indices(active_fe_index_transfer->active_fe_indices);

  // Update active FE indices on ghost cells.
  dealii::internal::hp::DoFHandlerImplementation::Implementation::
    communicate_active_fe_indices(*this);

  // Free memory.
  active_fe_index_transfer.reset();
#  endif
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
DoFHandler<dim, spacedim>::MGVertexDoFs::MGVertexDoFs()
  : coarsest_level(numbers::invalid_unsigned_int)
  , finest_level(0)
{}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void DoFHandler<dim, spacedim>::MGVertexDoFs::init(
  const unsigned int cl,
  const unsigned int fl,
  const unsigned int dofs_per_vertex)
{
  coarsest_level = cl;
  finest_level   = fl;

  if (coarsest_level <= finest_level)
    {
      const unsigned int n_levels  = finest_level - coarsest_level + 1;
      const unsigned int n_indices = n_levels * dofs_per_vertex;

      indices = std::make_unique<types::global_dof_index[]>(n_indices);
      std::fill(indices.get(),
                indices.get() + n_indices,
                numbers::invalid_dof_index);
    }
  else
    indices.reset();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
unsigned int DoFHandler<dim, spacedim>::MGVertexDoFs::get_coarsest_level() const
{
  return coarsest_level;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
unsigned int DoFHandler<dim, spacedim>::MGVertexDoFs::get_finest_level() const
{
  return finest_level;
}
#endif
/*-------------- Explicit Instantiations -------------------------------*/
#include "dofs/dof_handler.inst"



DEAL_II_NAMESPACE_CLOSE
