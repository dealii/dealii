// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_p4est_wrappers_h
#define dealii_p4est_wrappers_h

#include <deal.II/base/geometry_info.h>

#ifdef DEAL_II_WITH_P4EST
#  include <p4est_bits.h>
#  include <p4est_communication.h>
#  include <p4est_extended.h>
#  include <p4est_ghost.h>
#  include <p4est_iterate.h>
#  include <p4est_vtk.h>
#  include <p8est_bits.h>
#  include <p8est_communication.h>
#  include <p8est_extended.h>
#  include <p8est_ghost.h>
#  include <p8est_iterate.h>
#  include <p8est_vtk.h>

#  include <map>
#  include <set>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    class Triangulation;
  }
} // namespace parallel


namespace internal
{
  namespace p4est
  {
    /**
     * A structure whose explicit specializations contain typedefs to the
     * relevant p4est_* and p8est_* types. Using this structure, for example
     * by saying <tt>types<dim>::connectivity</tt> we can write code in a
     * dimension independent way, either referring to p4est_connectivity_t or
     * p8est_connectivity_t, depending on template argument.
     */
    template <int>
    struct types;

    template <>
    struct types<2>
    {
      typedef p4est_connectivity_t connectivity;
      typedef p4est_t              forest;
      typedef p4est_tree_t         tree;
      typedef p4est_quadrant_t     quadrant;
      typedef p4est_topidx_t       topidx;
      typedef p4est_locidx_t       locidx;
#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      typedef p4est_connect_type_t balance_type;
#  else
      typedef p4est_balance_type_t balance_type;
#  endif
      typedef p4est_ghost_t ghost;
    };

    template <>
    struct types<3>
    {
      typedef p8est_connectivity_t connectivity;
      typedef p8est_t              forest;
      typedef p8est_tree_t         tree;
      typedef p8est_quadrant_t     quadrant;
      typedef p4est_topidx_t       topidx;
      typedef p4est_locidx_t       locidx;
#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      typedef p8est_connect_type_t balance_type;
#  else
      typedef p8est_balance_type_t balance_type;
#  endif
      typedef p8est_ghost_t ghost;
    };



    /**
     * A structure whose explicit specializations contain pointers to the
     * relevant p4est_* and p8est_* functions. Using this structure, for
     * example by saying functions<dim>::quadrant_compare, we can write code
     * in a dimension independent way, either calling p4est_quadrant_compare
     * or p8est_quadrant_compare, depending on template argument.
     */
    template <int dim>
    struct functions;

    template <>
    struct functions<2>
    {
      static int (&quadrant_compare)(const void *v1, const void *v2);

      static void (&quadrant_childrenv)(const types<2>::quadrant *q,
                                        types<2>::quadrant        c[]);

      static int (&quadrant_overlaps_tree)(types<2>::tree *          tree,
                                           const types<2>::quadrant *q);

      static void (&quadrant_set_morton)(types<2>::quadrant *quadrant,
                                         int                 level,
                                         uint64_t            id);

      static int (&quadrant_is_equal)(const types<2>::quadrant *q1,
                                      const types<2>::quadrant *q2);

      static int (&quadrant_is_sibling)(const types<2>::quadrant *q1,
                                        const types<2>::quadrant *q2);

      static int (&quadrant_is_ancestor)(const types<2>::quadrant *q1,
                                         const types<2>::quadrant *q2);

      static int (&quadrant_ancestor_id)(const types<2>::quadrant *q,
                                         int                       level);

      static int (&comm_find_owner)(types<2>::forest *        p4est,
                                    const types<2>::locidx    which_tree,
                                    const types<2>::quadrant *q,
                                    const int                 guess);

      static types<2>::connectivity *(&connectivity_new)(
        types<2>::topidx num_vertices,
        types<2>::topidx num_trees,
        types<2>::topidx num_corners,
        types<2>::topidx num_vtt);
      static void (&connectivity_join_faces)(types<2>::connectivity *conn,
                                             types<2>::topidx        tree_left,
                                             types<2>::topidx        tree_right,
                                             int                     face_left,
                                             int                     face_right,
                                             int orientation);



      static void (&connectivity_destroy)(p4est_connectivity_t *connectivity);

      static types<2>::forest *(&new_forest)(
        MPI_Comm                mpicomm,
        types<2>::connectivity *connectivity,
        types<2>::locidx        min_quadrants,
        int                     min_level,
        int                     fill_uniform,
        size_t                  data_size,
        p4est_init_t            init_fn,
        void *                  user_pointer);

      static void (&destroy)(types<2>::forest *p4est);

      static void (&refine)(types<2>::forest *p4est,
                            int               refine_recursive,
                            p4est_refine_t    refine_fn,
                            p4est_init_t      init_fn);

      static void (&coarsen)(types<2>::forest *p4est,
                             int               coarsen_recursive,
                             p4est_coarsen_t   coarsen_fn,
                             p4est_init_t      init_fn);

      static void (&balance)(types<2>::forest *     p4est,
                             types<2>::balance_type btype,
                             p4est_init_t           init_fn);

#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static p4est_gloidx_t (&partition)(types<2>::forest *p4est,
                                         int partition_for_coarsening,
                                         p4est_weight_t weight_fn);
#  else
      static void (&partition)(types<2>::forest *p4est,
                               int               partition_for_coarsening,
                               p4est_weight_t    weight_fn);
#  endif

      static void (&save)(const char *      filename,
                          types<2>::forest *p4est,
                          int               save_data);

#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static types<2>::forest *(&load_ext)(const char *filename,
                                           MPI_Comm    mpicomm,
                                           size_t      data_size,
                                           int         load_data,
                                           int         autopartition,
                                           int         broadcasthead,
                                           void *      user_pointer,
                                           types<2>::connectivity **p4est);
#  else
      static types<2>::forest *(&load)(const char *             filename,
                                       MPI_Comm                 mpicomm,
                                       size_t                   data_size,
                                       int                      load_data,
                                       void *                   user_pointer,
                                       types<2>::connectivity **p4est);
#  endif

#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static int (&connectivity_save)(const char *            filename,
                                      types<2>::connectivity *connectivity);
#  else
      static void (&connectivity_save)(const char *            filename,
                                       types<2>::connectivity *connectivity);
#  endif

      static int (&connectivity_is_valid)(types<2>::connectivity *connectivity);

#  if DEAL_II_P4EST_VERSION_GTE(1, 0, 0, 0)
      static types<2>::connectivity *(&connectivity_load)(const char *filename,
                                                          size_t *    length);
#  elif DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static types<2>::connectivity *(
        &connectivity_load)(const char *filename, long unsigned *length);
#  else
      static types<2>::connectivity *(&connectivity_load)(const char *filename,
                                                          long *      length);
#  endif

      static unsigned int (&checksum)(types<2>::forest *p4est);

      static void (&vtk_write_file)(types<2>::forest *p4est,
                                    p4est_geometry_t *,
                                    const char *baseName);

      static types<2>::ghost *(&ghost_new)(types<2>::forest *     p4est,
                                           types<2>::balance_type btype);

      static void (&ghost_destroy)(types<2>::ghost *ghost);

      static void (&reset_data)(types<2>::forest *p4est,
                                size_t            data_size,
                                p4est_init_t      init_fn,
                                void *            user_pointer);

      static size_t (&forest_memory_used)(types<2>::forest *p4est);

      static size_t (&connectivity_memory_used)(types<2>::connectivity *p4est);

      template <int spacedim>
      static void
        iterate(dealii::internal::p4est::types<2>::forest *parallel_forest,
                dealii::internal::p4est::types<2>::ghost * parallel_ghost,
                void *                                     user_data);

      static const unsigned int max_level = P4EST_MAXLEVEL;
    };


    template <>
    struct functions<3>
    {
      static int (&quadrant_compare)(const void *v1, const void *v2);

      static void (&quadrant_childrenv)(const types<3>::quadrant *q,
                                        types<3>::quadrant        c[]);

      static int (&quadrant_overlaps_tree)(types<3>::tree *          tree,
                                           const types<3>::quadrant *q);

      static void (&quadrant_set_morton)(types<3>::quadrant *quadrant,
                                         int                 level,
                                         uint64_t            id);

      static int (&quadrant_is_equal)(const types<3>::quadrant *q1,
                                      const types<3>::quadrant *q2);

      static int (&quadrant_is_sibling)(const types<3>::quadrant *q1,
                                        const types<3>::quadrant *q2);

      static int (&quadrant_is_ancestor)(const types<3>::quadrant *q1,
                                         const types<3>::quadrant *q2);
      static int (&quadrant_ancestor_id)(const types<3>::quadrant *q,
                                         int                       level);

      static int (&comm_find_owner)(types<3>::forest *        p4est,
                                    const types<3>::locidx    which_tree,
                                    const types<3>::quadrant *q,
                                    const int                 guess);

      static types<3>::connectivity *(&connectivity_new)(
        types<3>::topidx num_vertices,
        types<3>::topidx num_trees,
        types<3>::topidx num_edges,
        types<3>::topidx num_ett,
        types<3>::topidx num_corners,
        types<3>::topidx num_ctt);

      static void (&connectivity_join_faces)(types<3>::connectivity *conn,
                                             types<3>::topidx        tree_left,
                                             types<3>::topidx        tree_right,
                                             int                     face_left,
                                             int                     face_right,
                                             int orientation);

      static void (&connectivity_destroy)(p8est_connectivity_t *connectivity);

      static types<3>::forest *(&new_forest)(
        MPI_Comm                mpicomm,
        types<3>::connectivity *connectivity,
        types<3>::locidx        min_quadrants,
        int                     min_level,
        int                     fill_uniform,
        size_t                  data_size,
        p8est_init_t            init_fn,
        void *                  user_pointer);

      static void (&destroy)(types<3>::forest *p8est);

      static void (&refine)(types<3>::forest *p8est,
                            int               refine_recursive,
                            p8est_refine_t    refine_fn,
                            p8est_init_t      init_fn);

      static void (&coarsen)(types<3>::forest *p8est,
                             int               coarsen_recursive,
                             p8est_coarsen_t   coarsen_fn,
                             p8est_init_t      init_fn);

      static void (&balance)(types<3>::forest *     p8est,
                             types<3>::balance_type btype,
                             p8est_init_t           init_fn);

#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static p4est_gloidx_t (&partition)(types<3>::forest *p8est,
                                         int partition_for_coarsening,
                                         p8est_weight_t weight_fn);
#  else
      static void (&partition)(types<3>::forest *p8est,
                               int               partition_for_coarsening,
                               p8est_weight_t    weight_fn);
#  endif

      static void (&save)(const char *      filename,
                          types<3>::forest *p4est,
                          int               save_data);

#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static types<3>::forest *(&load_ext)(const char *filename,
                                           MPI_Comm    mpicomm,
                                           std::size_t data_size,
                                           int         load_data,
                                           int         autopartition,
                                           int         broadcasthead,
                                           void *      user_pointer,
                                           types<3>::connectivity **p4est);
#  else
      static types<3>::forest *(&load)(const char *             filename,
                                       MPI_Comm                 mpicomm,
                                       std::size_t              data_size,
                                       int                      load_data,
                                       void *                   user_pointer,
                                       types<3>::connectivity **p4est);
#  endif

#  if DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static int (&connectivity_save)(const char *            filename,
                                      types<3>::connectivity *connectivity);
#  else
      static void (&connectivity_save)(const char *            filename,
                                       types<3>::connectivity *connectivity);
#  endif

      static int (&connectivity_is_valid)(types<3>::connectivity *connectivity);

#  if DEAL_II_P4EST_VERSION_GTE(1, 0, 0, 0)
      static types<3>::connectivity *(&connectivity_load)(const char *filename,
                                                          size_t *    length);
#  elif DEAL_II_P4EST_VERSION_GTE(0, 3, 4, 3)
      static types<3>::connectivity *(
        &connectivity_load)(const char *filename, long unsigned *length);
#  else
      static types<3>::connectivity *(&connectivity_load)(const char *filename,
                                                          long *      length);
#  endif

      static unsigned int (&checksum)(types<3>::forest *p8est);

      static void (&vtk_write_file)(types<3>::forest *p8est,
                                    p8est_geometry_t *,
                                    const char *baseName);
      static types<3>::ghost *(&ghost_new)(types<3>::forest *     p4est,
                                           types<3>::balance_type btype);

      static void (&ghost_destroy)(types<3>::ghost *ghost);

      static void (&reset_data)(types<3>::forest *p4est,
                                size_t            data_size,
                                p8est_init_t      init_fn,
                                void *            user_pointer);

      static size_t (&forest_memory_used)(types<3>::forest *p4est);

      static size_t (&connectivity_memory_used)(types<3>::connectivity *p4est);



      static const unsigned int max_level = P8EST_MAXLEVEL;
    };



    /**
     * This struct templatizes the p4est iterate structs and function
     * prototypes, which are used to execute callback functions for faces,
     * edges, and corners that require local neighborhood information, i.e.
     * the neighboring cells
     */
    template <int dim>
    struct iter;

    template <>
    struct iter<2>
    {
      typedef p4est_iter_corner_info_t corner_info;
      typedef p4est_iter_corner_side_t corner_side;
      typedef p4est_iter_corner_t      corner_iter;
      typedef p4est_iter_face_info_t   face_info;
      typedef p4est_iter_face_side_t   face_side;
      typedef p4est_iter_face_t        face_iter;
    };

    template <>
    struct iter<3>
    {
      typedef p8est_iter_corner_info_t corner_info;
      typedef p8est_iter_corner_side_t corner_side;
      typedef p8est_iter_corner_t      corner_iter;
      typedef p8est_iter_edge_info_t   edge_info;
      typedef p8est_iter_edge_side_t   edge_side;
      typedef p8est_iter_edge_t        edge_iter;
      typedef p8est_iter_face_info_t   face_info;
      typedef p8est_iter_face_side_t   face_side;
      typedef p8est_iter_face_t        face_iter;
    };



    /**
     * Initialize the GeometryInfo<dim>::max_children_per_cell children of the
     * cell p4est_cell.
     */
    template <int dim>
    void
    init_quadrant_children(
      const typename types<dim>::quadrant &p4est_cell,
      typename types<dim>::quadrant (
        &p4est_children)[dealii::GeometryInfo<dim>::max_children_per_cell]);



    /**
     * Initialize quadrant to represent a coarse cell.
     */
    template <int dim>
    void
    init_coarse_quadrant(typename types<dim>::quadrant &quad);



    /**
     * Return whether q1 and q2 are equal
     */
    template <int dim>
    bool
    quadrant_is_equal(const typename types<dim>::quadrant &q1,
                      const typename types<dim>::quadrant &q2);



    /**
     * Return whether q1 is an ancestor of q2
     */
    template <int dim>
    bool
    quadrant_is_ancestor(const typename types<dim>::quadrant &q1,
                         const typename types<dim>::quadrant &q2);



    /**
     * Return whether the children of a coarse cell are stored locally
     */
    template <int dim>
    bool
    tree_exists_locally(const typename types<dim>::forest *parallel_forest,
                        const typename types<dim>::topidx  coarse_grid_cell);



    /**
     * Compute the ghost neighbors surrounding each vertex by querying p4est
     */
    template <int dim, int spacedim>
    std::map<unsigned int, std::set<dealii::types::subdomain_id>>
    compute_vertices_with_ghost_neighbors(
      const dealii::parallel::distributed::Triangulation<dim, spacedim> &tria,
      typename dealii::internal::p4est::types<dim>::forest *parallel_forest,
      typename dealii::internal::p4est::types<dim>::ghost * parallel_ghost);

  } // namespace p4est
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_P4EST

#endif // dealii_p4est_wrappers_h
