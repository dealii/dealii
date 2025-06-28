// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria.h>

#ifdef DEAL_II_WITH_P4EST
#  include <p4est.h>
#  include <p8est.h>
#  include <sc_containers.h>

// Below, we will use the P4EST_QUADRANT_INIT and P8EST_QUADRANT_INIT
// function-like macros. If we are building the library based on
// header files, we get these from the <p4est.h> and <p8est.h> header
// inclusions. But if we build a C++20 module, we only import
// declarations, not preprocessor macros. As a consequence, let us
// duplicate these macros here, hoping that at some point, the p4est
// library folks add regular functions that can do the job.
#  ifndef P4EST_QUADRANT_INIT
#    define P4EST_QUADRANT_INIT(q) \
      ((void)std::memset((q), -1, sizeof(p4est_quadrant_t)))
#  endif

#  ifndef P8EST_QUADRANT_INIT
#    define P8EST_QUADRANT_INIT(q) \
      ((void)std::memset((q), -1, sizeof(p8est_quadrant_t)))
#  endif

#endif


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_P4EST

namespace internal
{
  namespace p4est
  {
    namespace
    {
      template <int dim, int spacedim>
      typename dealii::Triangulation<dim, spacedim>::cell_iterator
      cell_from_quad(
        const dealii::parallel::distributed::Triangulation<dim, spacedim>
          *triangulation,
        const typename dealii::internal::p4est::types<dim>::topidx    treeidx,
        const typename dealii::internal::p4est::types<dim>::quadrant &quad)
      {
        int                             i, l = quad.level;
        dealii::types::global_dof_index dealii_index =
          triangulation->get_p4est_tree_to_coarse_cell_permutation()[treeidx];

        for (i = 0; i < l; ++i)
          {
            typename dealii::Triangulation<dim, spacedim>::cell_iterator cell(
              triangulation, i, dealii_index);
            const int child_id =
              dealii::internal::p4est::functions<dim>::quadrant_ancestor_id(
                &quad, i + 1);
            Assert(cell->has_children(),
                   ExcMessage("p4est quadrant does not correspond to a cell!"));
            dealii_index = cell->child_index(child_id);
          }

        typename dealii::Triangulation<dim, spacedim>::cell_iterator out_cell(
          triangulation, l, dealii_index);

        return out_cell;
      }

      /**
       * This is the callback data structure used to fill
       * vertices_with_ghost_neighbors via the p4est_iterate tool
       */
      template <int dim, int spacedim>
      struct FindGhosts
      {
        const typename dealii::parallel::distributed::Triangulation<dim,
                                                                    spacedim>
                   *triangulation;
        sc_array_t *subids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id>>
          *vertices_with_ghost_neighbors;
      };


      /** At a corner (vertex), determine if any of the neighboring cells are
       * ghosts.  If there are, find out their subdomain ids, and if this is a
       * local vertex, then add these subdomain ids to the map
       * vertices_with_ghost_neighbors of that index
       */
      template <int dim, int spacedim>
      void
      find_ghosts_corner(
        typename dealii::internal::p4est::iter<dim>::corner_info *info,
        void                                                     *user_data)
      {
        int   i, j;
        int   nsides = info->sides.elem_count;
        auto *sides  = reinterpret_cast<
          typename dealii::internal::p4est::iter<dim>::corner_side *>(
          info->sides.array);
        FindGhosts<dim, spacedim> *fg =
          static_cast<FindGhosts<dim, spacedim> *>(user_data);
        sc_array_t *subids = fg->subids;
        const dealii::parallel::distributed::Triangulation<dim, spacedim>
                                    *triangulation = fg->triangulation;
        int                          nsubs;
        dealii::types::subdomain_id *subdomain_ids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id>>
          *vertices_with_ghost_neighbors = fg->vertices_with_ghost_neighbors;

        subids->elem_count = 0;
        for (i = 0; i < nsides; ++i)
          {
            if (sides[i].is_ghost)
              {
                typename dealii::parallel::distributed::
                  Triangulation<dim, spacedim>::cell_iterator cell =
                    cell_from_quad(triangulation,
                                   sides[i].treeid,
                                   *(sides[i].quad));
                Assert(cell->is_ghost(),
                       ExcMessage("ghost quad did not find ghost cell"));
                dealii::types::subdomain_id *subid =
                  static_cast<dealii::types::subdomain_id *>(
                    sc_array_push(subids));
                *subid = cell->subdomain_id();
              }
          }

        if (!subids->elem_count)
          {
            return;
          }

        nsubs = static_cast<int>(subids->elem_count);
        subdomain_ids =
          reinterpret_cast<dealii::types::subdomain_id *>(subids->array);

        for (i = 0; i < nsides; ++i)
          {
            if (!sides[i].is_ghost)
              {
                typename dealii::parallel::distributed::
                  Triangulation<dim, spacedim>::cell_iterator cell =
                    cell_from_quad(triangulation,
                                   sides[i].treeid,
                                   *(sides[i].quad));

                Assert(!cell->is_ghost(),
                       ExcMessage("local quad found ghost cell"));

                for (j = 0; j < nsubs; ++j)
                  {
                    (*vertices_with_ghost_neighbors)[cell->vertex_index(
                                                       sides[i].corner)]
                      .insert(subdomain_ids[j]);
                  }
              }
          }

        subids->elem_count = 0;
      }

      /** Similar to find_ghosts_corner, but for the hanging vertex in the
       * middle of an edge
       */
      template <int dim, int spacedim>
      void
      find_ghosts_edge(
        typename dealii::internal::p4est::iter<dim>::edge_info *info,
        void                                                   *user_data)
      {
        int   i, j, k;
        int   nsides = info->sides.elem_count;
        auto *sides  = reinterpret_cast<
          typename dealii::internal::p4est::iter<dim>::edge_side *>(
          info->sides.array);
        auto       *fg = static_cast<FindGhosts<dim, spacedim> *>(user_data);
        sc_array_t *subids = fg->subids;
        const dealii::parallel::distributed::Triangulation<dim, spacedim>
                                    *triangulation = fg->triangulation;
        int                          nsubs;
        dealii::types::subdomain_id *subdomain_ids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id>>
          *vertices_with_ghost_neighbors = fg->vertices_with_ghost_neighbors;

        subids->elem_count = 0;
        for (i = 0; i < nsides; ++i)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < 2; ++j)
                  {
                    if (sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::
                          Triangulation<dim, spacedim>::cell_iterator cell =
                            cell_from_quad(triangulation,
                                           sides[i].treeid,
                                           *(sides[i].is.hanging.quad[j]));
                        dealii::types::subdomain_id *subid =
                          static_cast<dealii::types::subdomain_id *>(
                            sc_array_push(subids));
                        *subid = cell->subdomain_id();
                      }
                  }
              }
          }

        if (!subids->elem_count)
          {
            return;
          }

        nsubs = static_cast<int>(subids->elem_count);
        subdomain_ids =
          reinterpret_cast<dealii::types::subdomain_id *>(subids->array);

        for (i = 0; i < nsides; ++i)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < 2; ++j)
                  {
                    if (!sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::
                          Triangulation<dim, spacedim>::cell_iterator cell =
                            cell_from_quad(triangulation,
                                           sides[i].treeid,
                                           *(sides[i].is.hanging.quad[j]));

                        for (k = 0; k < nsubs; ++k)
                          {
                            (*vertices_with_ghost_neighbors)
                              [cell->vertex_index(
                                 p8est_edge_corners[sides[i].edge][1 ^ j])]
                                .insert(subdomain_ids[k]);
                          }
                      }
                  }
              }
          }

        subids->elem_count = 0;
      }

      /** Similar to find_ghosts_corner, but for the hanging vertex in the
       * middle of a face
       */
      template <int dim, int spacedim>
      void
      find_ghosts_face(
        typename dealii::internal::p4est::iter<dim>::face_info *info,
        void                                                   *user_data)
      {
        int   i, j, k;
        int   nsides = info->sides.elem_count;
        auto *sides  = reinterpret_cast<
          typename dealii::internal::p4est::iter<dim>::face_side *>(
          info->sides.array);
        FindGhosts<dim, spacedim> *fg =
          static_cast<FindGhosts<dim, spacedim> *>(user_data);
        sc_array_t *subids = fg->subids;
        const dealii::parallel::distributed::Triangulation<dim, spacedim>
                                    *triangulation = fg->triangulation;
        int                          nsubs;
        dealii::types::subdomain_id *subdomain_ids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id>>
           *vertices_with_ghost_neighbors = fg->vertices_with_ghost_neighbors;
        int limit                         = (dim == 2) ? 2 : 4;

        subids->elem_count = 0;
        for (i = 0; i < nsides; ++i)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < limit; ++j)
                  {
                    if (sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::
                          Triangulation<dim, spacedim>::cell_iterator cell =
                            cell_from_quad(triangulation,
                                           sides[i].treeid,
                                           *(sides[i].is.hanging.quad[j]));
                        dealii::types::subdomain_id *subid =
                          static_cast<dealii::types::subdomain_id *>(
                            sc_array_push(subids));
                        *subid = cell->subdomain_id();
                      }
                  }
              }
          }

        if (!subids->elem_count)
          {
            return;
          }

        nsubs = static_cast<int>(subids->elem_count);
        subdomain_ids =
          reinterpret_cast<dealii::types::subdomain_id *>(subids->array);

        for (i = 0; i < nsides; ++i)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < limit; ++j)
                  {
                    if (!sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::
                          Triangulation<dim, spacedim>::cell_iterator cell =
                            cell_from_quad(triangulation,
                                           sides[i].treeid,
                                           *(sides[i].is.hanging.quad[j]));

                        for (k = 0; k < nsubs; ++k)
                          {
                            if (dim == 2)
                              {
                                (*vertices_with_ghost_neighbors)
                                  [cell->vertex_index(
                                     p4est_face_corners[sides[i].face]
                                                       [(limit - 1) ^ j])]
                                    .insert(subdomain_ids[k]);
                              }
                            else
                              {
                                (*vertices_with_ghost_neighbors)
                                  [cell->vertex_index(
                                     p8est_face_corners[sides[i].face]
                                                       [(limit - 1) ^ j])]
                                    .insert(subdomain_ids[k]);
                              }
                          }
                      }
                  }
              }
          }

        subids->elem_count = 0;
      }
    } // namespace


    int (&functions<2>::quadrant_compare)(const void *v1, const void *v2) =
      p4est_quadrant_compare;

    void (&functions<2>::quadrant_childrenv)(const types<2>::quadrant *q,
                                             types<2>::quadrant        c[]) =
      p4est_quadrant_childrenv;

    int (&functions<2>::quadrant_overlaps_tree)(types<2>::tree           *tree,
                                                const types<2>::quadrant *q) =
      p4est_quadrant_overlaps_tree;

    void (&functions<2>::quadrant_set_morton)(types<2>::quadrant *quadrant,
                                              int                 level,
                                              std::uint64_t       id) =
      p4est_quadrant_set_morton;

    void
    functions<2>::quadrant_init(types<2>::quadrant &q)
    {
      P4EST_QUADRANT_INIT(&q);
    }

    int (&functions<2>::quadrant_is_equal)(const types<2>::quadrant *q1,
                                           const types<2>::quadrant *q2) =
      p4est_quadrant_is_equal;

    int (&functions<2>::quadrant_is_sibling)(const types<2>::quadrant *q1,
                                             const types<2>::quadrant *q2) =
      p4est_quadrant_is_sibling;

    int (&functions<2>::quadrant_is_ancestor)(const types<2>::quadrant *q1,
                                              const types<2>::quadrant *q2) =
      p4est_quadrant_is_ancestor;

    int (&functions<2>::quadrant_ancestor_id)(const types<2>::quadrant *q,
                                              int                       level) =
      p4est_quadrant_ancestor_id;

    int (&functions<2>::comm_find_owner)(types<2>::forest         *p4est,
                                         const types<2>::locidx    which_tree,
                                         const types<2>::quadrant *q,
                                         const int                 guess) =
      p4est_comm_find_owner;

    types<2>::connectivity *(&functions<2>::connectivity_new)(
      types<2>::topidx num_vertices,
      types<2>::topidx num_trees,
      types<2>::topidx num_corners,
      types<2>::topidx num_vtt) = p4est_connectivity_new;

    types<2>::connectivity *(&functions<2>::connectivity_new_copy)(
      types<2>::topidx        num_vertices,
      types<2>::topidx        num_trees,
      types<2>::topidx        num_corners,
      const double           *vertices,
      const types<2>::topidx *ttv,
      const types<2>::topidx *ttt,
      const int8_t           *ttf,
      const types<2>::topidx *ttc,
      const types<2>::topidx *coff,
      const types<2>::topidx *ctt,
      const int8_t           *ctc) = p4est_connectivity_new_copy;

    void (&functions<2>::connectivity_join_faces)(types<2>::connectivity *conn,
                                                  types<2>::topidx tree_left,
                                                  types<2>::topidx tree_right,
                                                  int              face_left,
                                                  int              face_right,
                                                  int orientation) =
      p4est_connectivity_join_faces;

    void (&functions<2>::connectivity_destroy)(
      p4est_connectivity_t *connectivity) = p4est_connectivity_destroy;

    types<2>::forest *(&functions<2>::new_forest)(
      MPI_Comm                mpicomm,
      types<2>::connectivity *connectivity,
      types<2>::locidx        min_quadrants,
      int                     min_level,
      int                     fill_uniform,
      std::size_t             data_size,
      p4est_init_t            init_fn,
      void                   *user_pointer) = p4est_new_ext;

    types<2>::forest *(&functions<2>::copy_forest)(types<2>::forest *input,
                                                   int copy_data) = p4est_copy;

    void (&functions<2>::destroy)(types<2>::forest *p4est) = p4est_destroy;

    void (&functions<2>::refine)(types<2>::forest *p4est,
                                 int               refine_recursive,
                                 p4est_refine_t    refine_fn,
                                 p4est_init_t      init_fn) = p4est_refine;

    void (&functions<2>::coarsen)(types<2>::forest *p4est,
                                  int               coarsen_recursive,
                                  p4est_coarsen_t   coarsen_fn,
                                  p4est_init_t      init_fn) = p4est_coarsen;

    void (&functions<2>::balance)(types<2>::forest      *p4est,
                                  types<2>::balance_type btype,
                                  p4est_init_t init_fn) = p4est_balance;

    types<2>::gloidx (&functions<2>::partition)(types<2>::forest *p4est,
                                                int partition_for_coarsening,
                                                p4est_weight_t weight_fn) =
      p4est_partition_ext;

    void (&functions<2>::save)(const char       *filename,
                               types<2>::forest *p4est,
                               int               save_data) = p4est_save;

    types<2>::forest *(&functions<2>::load_ext)(
      const char              *filename,
      MPI_Comm                 mpicomm,
      std::size_t              data_size,
      int                      load_data,
      int                      autopartition,
      int                      broadcasthead,
      void                    *user_pointer,
      types<2>::connectivity **p4est) = p4est_load_ext;

    int (&functions<2>::connectivity_save)(
      const char             *filename,
      types<2>::connectivity *connectivity) = p4est_connectivity_save;

    int (&functions<2>::connectivity_is_valid)(
      types<2>::connectivity *connectivity) = p4est_connectivity_is_valid;

    types<2>::connectivity *(&functions<2>::connectivity_load)(
      const char  *filename,
      std::size_t *length) = p4est_connectivity_load;

    unsigned int (&functions<2>::checksum)(types<2>::forest *p4est) =
      p4est_checksum;

    void (&functions<2>::vtk_write_file)(types<2>::forest *p4est,
                                         p4est_geometry_t *,
                                         const char *baseName) =
      p4est_vtk_write_file;

    types<2>::ghost *(&functions<2>::ghost_new)(types<2>::forest      *p4est,
                                                types<2>::balance_type btype) =
      p4est_ghost_new;

    void (&functions<2>::ghost_destroy)(types<2>::ghost *ghost) =
      p4est_ghost_destroy;

    void (&functions<2>::reset_data)(types<2>::forest *p4est,
                                     std::size_t       data_size,
                                     p4est_init_t      init_fn,
                                     void *user_pointer) = p4est_reset_data;

    std::size_t (&functions<2>::forest_memory_used)(types<2>::forest *p4est) =
      p4est_memory_used;

    std::size_t (&functions<2>::connectivity_memory_used)(
      types<2>::connectivity *p4est) = p4est_connectivity_memory_used;

    constexpr unsigned int functions<2>::max_level;

    void (&functions<2>::transfer_fixed)(const types<2>::gloidx *dest_gfq,
                                         const types<2>::gloidx *src_gfq,
                                         MPI_Comm                mpicomm,
                                         int                     tag,
                                         void                   *dest_data,
                                         const void             *src_data,
                                         std::size_t             data_size) =
      p4est_transfer_fixed;

    types<2>::transfer_context *(&functions<2>::transfer_fixed_begin)(
      const types<2>::gloidx *dest_gfq,
      const types<2>::gloidx *src_gfq,
      MPI_Comm                mpicomm,
      int                     tag,
      void                   *dest_data,
      const void             *src_data,
      std::size_t             data_size) = p4est_transfer_fixed_begin;

    void (&functions<2>::transfer_fixed_end)(types<2>::transfer_context *tc) =
      p4est_transfer_fixed_end;

    void (&functions<2>::transfer_custom)(const types<2>::gloidx *dest_gfq,
                                          const types<2>::gloidx *src_gfq,
                                          MPI_Comm                mpicomm,
                                          int                     tag,
                                          void                   *dest_data,
                                          const int              *dest_sizes,
                                          const void             *src_data,
                                          const int              *src_sizes) =
      p4est_transfer_custom;

    types<2>::transfer_context *(&functions<2>::transfer_custom_begin)(
      const types<2>::gloidx *dest_gfq,
      const types<2>::gloidx *src_gfq,
      MPI_Comm                mpicomm,
      int                     tag,
      void                   *dest_data,
      const int              *dest_sizes,
      const void             *src_data,
      const int              *src_sizes) = p4est_transfer_custom_begin;

    void (&functions<2>::transfer_custom_end)(types<2>::transfer_context *tc) =
      p4est_transfer_custom_end;

    void (&functions<2>::search_partition)(
      types<2>::forest                   *p4est,
      int                                 call_post,
      types<2>::search_partition_callback quadrant_fn,
      types<2>::search_partition_callback point_fn,
      sc_array_t                         *points) = p4est_search_partition;

    void (&functions<2>::quadrant_coord_to_vertex)(
      types<2>::connectivity  *connectivity,
      types<2>::topidx         treeid,
      types<2>::quadrant_coord x,
      types<2>::quadrant_coord y,
      double                   vxyz[3]) = p4est_qcoord_to_vertex;

    int (&functions<3>::quadrant_compare)(const void *v1, const void *v2) =
      p8est_quadrant_compare;

    void (&functions<3>::quadrant_childrenv)(const types<3>::quadrant *q,
                                             types<3>::quadrant        c[]) =
      p8est_quadrant_childrenv;

    int (&functions<3>::quadrant_overlaps_tree)(types<3>::tree           *tree,
                                                const types<3>::quadrant *q) =
      p8est_quadrant_overlaps_tree;

    void (&functions<3>::quadrant_set_morton)(types<3>::quadrant *quadrant,
                                              int                 level,
                                              std::uint64_t       id) =
      p8est_quadrant_set_morton;

    void
    functions<3>::quadrant_init(types<3>::quadrant &q)
    {
      P8EST_QUADRANT_INIT(&q);
    }

    int (&functions<3>::quadrant_is_equal)(const types<3>::quadrant *q1,
                                           const types<3>::quadrant *q2) =
      p8est_quadrant_is_equal;

    int (&functions<3>::quadrant_is_sibling)(const types<3>::quadrant *q1,
                                             const types<3>::quadrant *q2) =
      p8est_quadrant_is_sibling;

    int (&functions<3>::quadrant_is_ancestor)(const types<3>::quadrant *q1,
                                              const types<3>::quadrant *q2) =
      p8est_quadrant_is_ancestor;

    int (&functions<3>::quadrant_ancestor_id)(const types<3>::quadrant *q,
                                              int                       level) =
      p8est_quadrant_ancestor_id;

    int (&functions<3>::comm_find_owner)(types<3>::forest         *p4est,
                                         const types<3>::locidx    which_tree,
                                         const types<3>::quadrant *q,
                                         const int                 guess) =
      p8est_comm_find_owner;

    types<3>::connectivity *(&functions<3>::connectivity_new)(
      types<3>::topidx num_vertices,
      types<3>::topidx num_trees,
      types<3>::topidx num_edges,
      types<3>::topidx num_ett,
      types<3>::topidx num_corners,
      types<3>::topidx num_ctt) = p8est_connectivity_new;

    types<3>::connectivity *(&functions<3>::connectivity_new_copy)(
      types<3>::topidx        num_vertices,
      types<3>::topidx        num_trees,
      types<3>::topidx        num_edges,
      types<3>::topidx        num_corners,
      const double           *vertices,
      const types<3>::topidx *ttv,
      const types<3>::topidx *ttt,
      const int8_t           *ttf,
      const types<3>::topidx *tte,
      const types<3>::topidx *eoff,
      const types<3>::topidx *ett,
      const int8_t           *ete,
      const types<3>::topidx *ttc,
      const types<3>::topidx *coff,
      const types<3>::topidx *ctt,
      const int8_t           *ctc) = p8est_connectivity_new_copy;

    void (&functions<3>::connectivity_destroy)(
      p8est_connectivity_t *connectivity) = p8est_connectivity_destroy;

    void (&functions<3>::connectivity_join_faces)(types<3>::connectivity *conn,
                                                  types<3>::topidx tree_left,
                                                  types<3>::topidx tree_right,
                                                  int              face_left,
                                                  int              face_right,
                                                  int orientation) =
      p8est_connectivity_join_faces;

    types<3>::forest *(&functions<3>::new_forest)(
      MPI_Comm                mpicomm,
      types<3>::connectivity *connectivity,
      types<3>::locidx        min_quadrants,
      int                     min_level,
      int                     fill_uniform,
      std::size_t             data_size,
      p8est_init_t            init_fn,
      void                   *user_pointer) = p8est_new_ext;

    types<3>::forest *(&functions<3>::copy_forest)(types<3>::forest *input,
                                                   int copy_data) = p8est_copy;

    void (&functions<3>::destroy)(types<3>::forest *p8est) = p8est_destroy;

    void (&functions<3>::refine)(types<3>::forest *p8est,
                                 int               refine_recursive,
                                 p8est_refine_t    refine_fn,
                                 p8est_init_t      init_fn) = p8est_refine;

    void (&functions<3>::coarsen)(types<3>::forest *p8est,
                                  int               coarsen_recursive,
                                  p8est_coarsen_t   coarsen_fn,
                                  p8est_init_t      init_fn) = p8est_coarsen;

    void (&functions<3>::balance)(types<3>::forest      *p8est,
                                  types<3>::balance_type btype,
                                  p8est_init_t init_fn) = p8est_balance;

    types<3>::gloidx (&functions<3>::partition)(types<3>::forest *p8est,
                                                int partition_for_coarsening,
                                                p8est_weight_t weight_fn) =
      p8est_partition_ext;

    void (&functions<3>::save)(const char       *filename,
                               types<3>::forest *p4est,
                               int               save_data) = p8est_save;

    types<3>::forest *(&functions<3>::load_ext)(
      const char              *filename,
      MPI_Comm                 mpicomm,
      std::size_t              data_size,
      int                      load_data,
      int                      autopartition,
      int                      broadcasthead,
      void                    *user_pointer,
      types<3>::connectivity **p4est) = p8est_load_ext;

    int (&functions<3>::connectivity_save)(
      const char             *filename,
      types<3>::connectivity *connectivity) = p8est_connectivity_save;

    int (&functions<3>::connectivity_is_valid)(
      types<3>::connectivity *connectivity) = p8est_connectivity_is_valid;

    types<3>::connectivity *(&functions<3>::connectivity_load)(
      const char  *filename,
      std::size_t *length) = p8est_connectivity_load;

    unsigned int (&functions<3>::checksum)(types<3>::forest *p8est) =
      p8est_checksum;

    void (&functions<3>::vtk_write_file)(types<3>::forest *p8est,
                                         p8est_geometry_t *,
                                         const char *baseName) =
      p8est_vtk_write_file;

    types<3>::ghost *(&functions<3>::ghost_new)(types<3>::forest      *p4est,
                                                types<3>::balance_type btype) =
      p8est_ghost_new;

    void (&functions<3>::ghost_destroy)(types<3>::ghost *ghost) =
      p8est_ghost_destroy;

    void (&functions<3>::reset_data)(types<3>::forest *p4est,
                                     std::size_t       data_size,
                                     p8est_init_t      init_fn,
                                     void *user_pointer) = p8est_reset_data;

    std::size_t (&functions<3>::forest_memory_used)(types<3>::forest *p4est) =
      p8est_memory_used;

    std::size_t (&functions<3>::connectivity_memory_used)(
      types<3>::connectivity *p4est) = p8est_connectivity_memory_used;

    constexpr unsigned int functions<3>::max_level;

    void (&functions<3>::transfer_fixed)(const types<3>::gloidx *dest_gfq,
                                         const types<3>::gloidx *src_gfq,
                                         MPI_Comm                mpicomm,
                                         int                     tag,
                                         void                   *dest_data,
                                         const void             *src_data,
                                         std::size_t             data_size) =
      p8est_transfer_fixed;

    types<3>::transfer_context *(&functions<3>::transfer_fixed_begin)(
      const types<3>::gloidx *dest_gfq,
      const types<3>::gloidx *src_gfq,
      MPI_Comm                mpicomm,
      int                     tag,
      void                   *dest_data,
      const void             *src_data,
      std::size_t             data_size) = p8est_transfer_fixed_begin;

    void (&functions<3>::transfer_fixed_end)(types<3>::transfer_context *tc) =
      p8est_transfer_fixed_end;

    void (&functions<3>::transfer_custom)(const types<3>::gloidx *dest_gfq,
                                          const types<3>::gloidx *src_gfq,
                                          MPI_Comm                mpicomm,
                                          int                     tag,
                                          void                   *dest_data,
                                          const int              *dest_sizes,
                                          const void             *src_data,
                                          const int              *src_sizes) =
      p8est_transfer_custom;

    types<3>::transfer_context *(&functions<3>::transfer_custom_begin)(
      const types<3>::gloidx *dest_gfq,
      const types<3>::gloidx *src_gfq,
      MPI_Comm                mpicomm,
      int                     tag,
      void                   *dest_data,
      const int              *dest_sizes,
      const void             *src_data,
      const int              *src_sizes) = p8est_transfer_custom_begin;

    void (&functions<3>::transfer_custom_end)(types<3>::transfer_context *tc) =
      p8est_transfer_custom_end;

    void (&functions<3>::search_partition)(
      types<3>::forest                   *p4est,
      int                                 call_post,
      types<3>::search_partition_callback quadrant_fn,
      types<3>::search_partition_callback point_fn,
      sc_array_t                         *points) = p8est_search_partition;

    void (&functions<3>::quadrant_coord_to_vertex)(
      types<3>::connectivity  *connectivity,
      types<3>::topidx         treeid,
      types<3>::quadrant_coord x,
      types<3>::quadrant_coord y,
      types<3>::quadrant_coord z,
      double                   vxyz[3]) = p8est_qcoord_to_vertex;

    template <int dim>
    void
    init_quadrant_children(
      const typename types<dim>::quadrant &p4est_cell,
      typename types<dim>::quadrant (
        &p4est_children)[dealii::GeometryInfo<dim>::max_children_per_cell])
    {
      for (unsigned int c = 0;
           c < dealii::GeometryInfo<dim>::max_children_per_cell;
           ++c)
        functions<dim>::quadrant_init(p4est_children[c]);

      functions<dim>::quadrant_childrenv(&p4est_cell, p4est_children);
    }

    template <int dim>
    void
    init_coarse_quadrant(typename types<dim>::quadrant &quad)
    {
      functions<dim>::quadrant_init(quad);
      functions<dim>::quadrant_set_morton(&quad,
                                          /*level=*/0,
                                          /*index=*/0);
    }

    template <int dim>
    bool
    quadrant_is_equal(const typename types<dim>::quadrant &q1,
                      const typename types<dim>::quadrant &q2)
    {
      return functions<dim>::quadrant_is_equal(&q1, &q2);
    }



    template <int dim>
    bool
    quadrant_is_ancestor(const typename types<dim>::quadrant &q1,
                         const typename types<dim>::quadrant &q2)
    {
      return functions<dim>::quadrant_is_ancestor(&q1, &q2);
    }

    template <int dim>
    bool
    tree_exists_locally(const typename types<dim>::forest *parallel_forest,
                        const typename types<dim>::topidx  coarse_grid_cell)
    {
      Assert(coarse_grid_cell < parallel_forest->connectivity->num_trees,
             ExcInternalError());
      return ((coarse_grid_cell >= parallel_forest->first_local_tree) &&
              (coarse_grid_cell <= parallel_forest->last_local_tree));
    }



    // template specializations

    template <>
    typename types<2>::connectivity *
    copy_connectivity<2>(const typename types<2>::connectivity *connectivity)
    {
      return functions<2>::connectivity_new_copy(
        connectivity->num_vertices,
        connectivity->num_trees,
        connectivity->num_corners,
        connectivity->vertices,
        connectivity->tree_to_vertex,
        connectivity->tree_to_tree,
        connectivity->tree_to_face,
        connectivity->tree_to_corner,
        connectivity->ctt_offset,
        connectivity->corner_to_tree,
        connectivity->corner_to_corner);
    }

    template <>
    typename types<3>::connectivity *
    copy_connectivity<3>(const typename types<3>::connectivity *connectivity)
    {
      return functions<3>::connectivity_new_copy(
        connectivity->num_vertices,
        connectivity->num_trees,
        connectivity->num_edges,
        connectivity->num_corners,
        connectivity->vertices,
        connectivity->tree_to_vertex,
        connectivity->tree_to_tree,
        connectivity->tree_to_face,
        connectivity->tree_to_edge,
        connectivity->ett_offset,
        connectivity->edge_to_tree,
        connectivity->edge_to_edge,
        connectivity->tree_to_corner,
        connectivity->ctt_offset,
        connectivity->corner_to_tree,
        connectivity->corner_to_corner);
    }



    template <>
    bool
    quadrant_is_equal<1>(const typename types<1>::quadrant &q1,
                         const typename types<1>::quadrant &q2)
    {
      return q1 == q2;
    }



    template <>
    bool
    quadrant_is_ancestor<1>(const types<1>::quadrant &q1,
                            const types<1>::quadrant &q2)
    {
      // determine level of quadrants
      const int level_1 = (q1 << types<1>::max_n_child_indices_bits) >>
                          types<1>::max_n_child_indices_bits;
      const int level_2 = (q2 << types<1>::max_n_child_indices_bits) >>
                          types<1>::max_n_child_indices_bits;

      // q1 can be an ancestor of q2 if q1's level is smaller
      if (level_1 >= level_2)
        return false;

      // extract path of quadrants up to level of possible ancestor q1
      const int truncated_id_1 = (q1 >> (types<1>::n_bits - 1 - level_1))
                                 << (types<1>::n_bits - 1 - level_1);
      const int truncated_id_2 = (q2 >> (types<1>::n_bits - 1 - level_1))
                                 << (types<1>::n_bits - 1 - level_1);

      // compare paths
      return truncated_id_1 == truncated_id_2;
    }



    template <>
    void
    init_quadrant_children<1>(
      const typename types<1>::quadrant &q,
      typename types<1>::quadrant (
        &p4est_children)[dealii::GeometryInfo<1>::max_children_per_cell])
    {
      // determine the current level of quadrant
      const int level_parent = (q << types<1>::max_n_child_indices_bits) >>
                               types<1>::max_n_child_indices_bits;
      const int level_child = level_parent + 1;

      // left child: only n_child_indices has to be incremented
      p4est_children[0] = (q + 1);

      // right child: increment and set a bit to 1 indicating that it is a right
      // child
      p4est_children[1] = (q + 1) | (1 << (types<1>::n_bits - 1 - level_child));
    }



    template <>
    void
    init_coarse_quadrant<1>(typename types<1>::quadrant &quad)
    {
      quad = 0;
    }

  } // namespace p4est
} // namespace internal

#endif // DEAL_II_WITH_P4EST

/*-------------- Explicit Instantiations -------------------------------*/
#include "distributed/p4est_wrappers.inst"


DEAL_II_NAMESPACE_CLOSE
