// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/tria.h>

#ifdef DEAL_II_WITH_P4EST
#  include <p4est_bits.h>
#  include <p4est_extended.h>
#  include <p4est_vtk.h>
#  include <p4est_ghost.h>
#  include <p4est_communication.h>
#  include <p4est_iterate.h>

#  include <p8est_bits.h>
#  include <p8est_extended.h>
#  include <p8est_vtk.h>
#  include <p8est_ghost.h>
#  include <p8est_communication.h>
#  include <p8est_iterate.h>
#endif

#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>


DEAL_II_NAMESPACE_OPEN


#ifdef DEAL_II_WITH_P4EST

namespace internal
{
  namespace p4est
  {
    /**
     * A structure whose explicit specializations contain pointers to the
     * relevant p4est_* and p8est_* functions. Using this structure, for
     * example by saying functions<dim>::quadrant_compare, we can write code
     * in a dimension independent way, either calling p4est_quadrant_compare
     * or p8est_quadrant_compare, depending on template argument.
     */
    template <int dim> struct functions;

    template <> struct functions<2>
    {
      static
      int (&quadrant_compare) (const void *v1, const void *v2);

      static
      void (&quadrant_childrenv) (const types<2>::quadrant *q,
                                  types<2>::quadrant c[]);

      static
      int (&quadrant_overlaps_tree) (types<2>::tree *tree,
                                     const types<2>::quadrant *q);

      static
      void (&quadrant_set_morton) (types<2>::quadrant *quadrant,
                                   int level,
                                   uint64_t id);

      static
      int (&quadrant_is_equal) (const types<2>::quadrant *q1,
                                const types<2>::quadrant *q2);

      static
      int (&quadrant_is_sibling) (const types<2>::quadrant *q1,
                                  const types<2>::quadrant *q2);

      static
      int (&quadrant_is_ancestor) (const types<2>::quadrant *q1,
                                   const types<2>::quadrant *q2);

      static
      int (&quadrant_ancestor_id) (const types<2>::quadrant *q,
                                   int level);

      static
      int (&comm_find_owner) (types<2>::forest *p4est,
                              const types<2>::locidx which_tree,
                              const types<2>::quadrant *q,
                              const int guess);

      static
      types<2>::connectivity *(&connectivity_new) (types<2>::topidx num_vertices,
                                                   types<2>::topidx num_trees,
                                                   types<2>::topidx num_corners,
                                                   types<2>::topidx num_vtt);

      static
      void (&connectivity_join_faces) (types<2>::connectivity *conn,
                                       types<2>::topidx tree_left,
                                       types<2>::topidx tree_right,
                                       int face_left,
                                       int face_right,
                                       int orientation);



      static
      void (&connectivity_destroy) (p4est_connectivity_t *connectivity);

      static
      types<2>::forest *(&new_forest) (MPI_Comm mpicomm,
                                       types<2>::connectivity *connectivity,
                                       types<2>::locidx min_quadrants,
                                       int min_level,
                                       int fill_uniform,
                                       size_t data_size,
                                       p4est_init_t init_fn,
                                       void *user_pointer);

      static
      void (&destroy) (types<2>::forest *p4est);

      static
      void (&refine) (types<2>::forest *p4est,
                      int refine_recursive,
                      p4est_refine_t refine_fn,
                      p4est_init_t init_fn);

      static
      void (&coarsen) (types<2>::forest *p4est,
                       int coarsen_recursive,
                       p4est_coarsen_t coarsen_fn,
                       p4est_init_t init_fn);

      static
      void (&balance) (types<2>::forest *p4est,
                       types<2>::balance_type btype,
                       p4est_init_t init_fn);

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      p4est_gloidx_t (&partition) (types<2>::forest *p4est,
                                   int partition_for_coarsening,
                                   p4est_weight_t weight_fn);
#else
      static
      void (&partition) (types<2>::forest *p4est,
                         int partition_for_coarsening,
                         p4est_weight_t weight_fn);
#endif

      static
      void (&save) (const char *filename,
                    types<2>::forest *p4est,
                    int save_data);

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      types<2>::forest *(&load_ext) (const char *filename,
                                     MPI_Comm mpicomm,
                                     size_t data_size,
                                     int load_data,
                                     int autopartition,
                                     int broadcasthead,
                                     void *user_pointer,
                                     types<2>::connectivity **p4est);
#else
      static
      types<2>::forest *(&load) (const char *filename,
                                 MPI_Comm mpicomm,
                                 size_t data_size,
                                 int load_data,
                                 void *user_pointer,
                                 types<2>::connectivity **p4est);
#endif

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      int (&connectivity_save) (const char *filename,
                                types<2>::connectivity *connectivity);
#else
      static
      void (&connectivity_save) (const char *filename,
                                 types<2>::connectivity *connectivity);
#endif

      static
      int (&connectivity_is_valid) (types<2>::connectivity *connectivity);

#if DEAL_II_P4EST_VERSION_GTE(1,0,0,0)
      static
      types<2>::connectivity *(&connectivity_load) (const char *filename,
                                                    size_t *length);
#elif DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      types<2>::connectivity *(&connectivity_load) (const char *filename,
                                                    long unsigned *length);
#else
      static
      types<2>::connectivity *(&connectivity_load) (const char *filename,
                                                    long *length);
#endif

      static
      unsigned int (&checksum) (types<2>::forest *p4est);

      static
      void (&vtk_write_file) (types<2>::forest *p4est,
                              p4est_geometry_t *,
                              const char *baseName);

      static
      types<2>::ghost *(&ghost_new) (types<2>::forest *p4est,
                                     types<2>::balance_type btype);

      static
      void (&ghost_destroy) (types<2>::ghost *ghost);

      static
      void (&reset_data) (types<2>::forest *p4est,
                          size_t data_size,
                          p4est_init_t init_fn,
                          void *user_pointer);

      static
      size_t (&forest_memory_used) (types<2>::forest *p4est);

      static
      size_t (&connectivity_memory_used) (types<2>::connectivity *p4est);

      static const unsigned max_level;
    };

    int (&functions<2>::quadrant_compare) (const void *v1, const void *v2)
      = p4est_quadrant_compare;

    void (&functions<2>::quadrant_childrenv) (const types<2>::quadrant *q,
                                              types<2>::quadrant c[])
      = p4est_quadrant_childrenv;

    int (&functions<2>::quadrant_overlaps_tree) (types<2>::tree *tree,
                                                 const types<2>::quadrant *q)
      = p4est_quadrant_overlaps_tree;

    void (&functions<2>::quadrant_set_morton) (types<2>::quadrant *quadrant,
                                               int level,
                                               uint64_t id)
      = p4est_quadrant_set_morton;

    int (&functions<2>::quadrant_is_equal) (const types<2>::quadrant *q1,
                                            const types<2>::quadrant *q2)
      = p4est_quadrant_is_equal;

    int (&functions<2>::quadrant_is_sibling) (const types<2>::quadrant *q1,
                                              const types<2>::quadrant *q2)
      = p4est_quadrant_is_sibling;

    int (&functions<2>::quadrant_is_ancestor) (const types<2>::quadrant *q1,
                                               const types<2>::quadrant *q2)
      = p4est_quadrant_is_ancestor;

    int (&functions<2>::quadrant_ancestor_id) (const types<2>::quadrant *q,
                                               int level)
      = p4est_quadrant_ancestor_id;

    int (&functions<2>::comm_find_owner) (types<2>::forest *p4est,
                                          const types<2>::locidx which_tree,
                                          const types<2>::quadrant *q,
                                          const int guess)
      = p4est_comm_find_owner;

    types<2>::connectivity *(&functions<2>::connectivity_new) (types<2>::topidx num_vertices,
                                                               types<2>::topidx num_trees,
                                                               types<2>::topidx num_corners,
                                                               types<2>::topidx num_vtt)
      = p4est_connectivity_new;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,1)
    void (&functions<2>::connectivity_join_faces) (types<2>::connectivity *conn,
                                                   types<2>::topidx tree_left,
                                                   types<2>::topidx tree_right,
                                                   int face_left,
                                                   int face_right,
                                                   int orientation)
      = p4est_connectivity_join_faces;
#endif

    void (&functions<2>::connectivity_destroy) (p4est_connectivity_t *connectivity)
      = p4est_connectivity_destroy;

    types<2>::forest *(&functions<2>::new_forest) (MPI_Comm mpicomm,
                                                   types<2>::connectivity *connectivity,
                                                   types<2>::locidx min_quadrants,
                                                   int min_level,
                                                   int fill_uniform,
                                                   size_t data_size,
                                                   p4est_init_t init_fn,
                                                   void *user_pointer)
      = p4est_new_ext;

    void (&functions<2>::destroy) (types<2>::forest *p4est)
      = p4est_destroy;

    void (&functions<2>::refine) (types<2>::forest *p4est,
                                  int refine_recursive,
                                  p4est_refine_t refine_fn,
                                  p4est_init_t init_fn)
      = p4est_refine;

    void (&functions<2>::coarsen) (types<2>::forest *p4est,
                                   int coarsen_recursive,
                                   p4est_coarsen_t coarsen_fn,
                                   p4est_init_t init_fn)
      = p4est_coarsen;

    void (&functions<2>::balance) (types<2>::forest *p4est,
                                   types<2>::balance_type btype,
                                   p4est_init_t init_fn)
      = p4est_balance;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    p4est_gloidx_t (&functions<2>::partition) (types<2>::forest *p4est,
                                               int partition_for_coarsening,
                                               p4est_weight_t weight_fn)
      = p4est_partition_ext;
#else
    void (&functions<2>::partition) (types<2>::forest *p4est,
                                     int partition_for_coarsening,
                                     p4est_weight_t weight_fn)
      = p4est_partition_ext;
#endif

    void (&functions<2>::save) (const char *filename,
                                types<2>::forest *p4est,
                                int save_data)
      = p4est_save;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    types<2>::forest *
    (&functions<2>::load_ext) (const char *filename,
                               MPI_Comm mpicomm,
                               std::size_t data_size,
                               int load_data,
                               int autopartition,
                               int broadcasthead,
                               void *user_pointer,
                               types<2>::connectivity **p4est)
      = p4est_load_ext;
#else
    types<2>::forest *
    (&functions<2>::load) (const char *filename,
                           MPI_Comm mpicomm,
                           std::size_t data_size,
                           int load_data,
                           void *user_pointer,
                           types<2>::connectivity **p4est)
      = p4est_load;
#endif

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    int (&functions<2>::connectivity_save) (const char *filename,
                                            types<2>::connectivity *connectivity)
      = p4est_connectivity_save;
#else
    void (&functions<2>::connectivity_save) (const char *filename,
                                             types<2>::connectivity *connectivity)
      = p4est_connectivity_save;
#endif

    int (&functions<2>::connectivity_is_valid) (types<2>::connectivity
                                                *connectivity)
      = p4est_connectivity_is_valid;

#if DEAL_II_P4EST_VERSION_GTE(1,0,0,0)
    types<2>::connectivity *
    (&functions<2>::connectivity_load) (const char *filename,
                                        size_t *length)
      = p4est_connectivity_load;
#elif DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    types<2>::connectivity *
    (&functions<2>::connectivity_load) (const char *filename,
                                        long unsigned *length)
      = p4est_connectivity_load;
#else
    types<2>::connectivity *
    (&functions<2>::connectivity_load) (const char *filename,
                                        long *length)
      = p4est_connectivity_load;
#endif

    unsigned int (&functions<2>::checksum) (types<2>::forest *p4est)
      = p4est_checksum;

    void (&functions<2>::vtk_write_file) (types<2>::forest *p4est,
                                          p4est_geometry_t *,
                                          const char *baseName)
      = p4est_vtk_write_file;

    types<2>::ghost *(&functions<2>::ghost_new) (types<2>::forest *p4est,
                                                 types<2>::balance_type btype)
      = p4est_ghost_new;

    void (&functions<2>::ghost_destroy) (types<2>::ghost *ghost)
      = p4est_ghost_destroy;

    void (&functions<2>::reset_data) (types<2>::forest *p4est,
                                      size_t data_size,
                                      p4est_init_t init_fn,
                                      void *user_pointer)
      = p4est_reset_data;

    size_t (&functions<2>::forest_memory_used) (types<2>::forest *p4est)
      = p4est_memory_used;

    size_t (&functions<2>::connectivity_memory_used) (types<2>::connectivity *p4est)
      = p4est_connectivity_memory_used;

    const unsigned int functions<2>::max_level = P4EST_MAXLEVEL;

    template <> struct functions<3>
    {
      static
      int (&quadrant_compare) (const void *v1, const void *v2);

      static
      void (&quadrant_childrenv) (const types<3>::quadrant *q,
                                  types<3>::quadrant c[]);

      static
      int (&quadrant_overlaps_tree) (types<3>::tree *tree,
                                     const types<3>::quadrant *q);

      static
      void (&quadrant_set_morton) (types<3>::quadrant *quadrant,
                                   int level,
                                   uint64_t id);

      static
      int (&quadrant_is_equal) (const types<3>::quadrant *q1,
                                const types<3>::quadrant *q2);

      static
      int (&quadrant_is_sibling) (const types<3>::quadrant *q1,
                                  const types<3>::quadrant *q2);

      static
      int (&quadrant_is_ancestor) (const types<3>::quadrant *q1,
                                   const types<3>::quadrant *q2);

      static
      int (&quadrant_ancestor_id) (const types<3>::quadrant *q,
                                   int level);

      static
      int (&comm_find_owner) (types<3>::forest *p4est,
                              const types<3>::locidx which_tree,
                              const types<3>::quadrant *q,
                              const int guess);

      static
      types<3>::connectivity *(&connectivity_new) (types<3>::topidx num_vertices,
                                                   types<3>::topidx num_trees,
                                                   types<3>::topidx num_edges,
                                                   types<3>::topidx num_ett,
                                                   types<3>::topidx num_corners,
                                                   types<3>::topidx num_ctt);

      static
      void (&connectivity_join_faces) (types<3>::connectivity *conn,
                                       types<3>::topidx tree_left,
                                       types<3>::topidx tree_right,
                                       int face_left,
                                       int face_right,
                                       int orientation);

      static
      void (&connectivity_destroy) (p8est_connectivity_t *connectivity);

      static
      types<3>::forest *(&new_forest) (MPI_Comm mpicomm,
                                       types<3>::connectivity *connectivity,
                                       types<3>::locidx min_quadrants,
                                       int min_level,
                                       int fill_uniform,
                                       size_t data_size,
                                       p8est_init_t init_fn,
                                       void *user_pointer);

      static
      void (&destroy) (types<3>::forest *p8est);

      static
      void (&refine) (types<3>::forest *p8est,
                      int refine_recursive,
                      p8est_refine_t refine_fn,
                      p8est_init_t init_fn);

      static
      void (&coarsen) (types<3>::forest *p8est,
                       int coarsen_recursive,
                       p8est_coarsen_t coarsen_fn,
                       p8est_init_t init_fn);

      static
      void (&balance) (types<3>::forest *p8est,
                       types<3>::balance_type btype,
                       p8est_init_t init_fn);

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      p4est_gloidx_t (&partition) (types<3>::forest *p8est,
                                   int partition_for_coarsening,
                                   p8est_weight_t weight_fn);
#else
      static
      void (&partition) (types<3>::forest *p8est,
                         int partition_for_coarsening,
                         p8est_weight_t weight_fn);
#endif

      static
      void (&save) (const char *filename,
                    types<3>::forest *p4est,
                    int save_data);

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      types<3>::forest *(&load_ext) (const char *filename,
                                     MPI_Comm mpicomm,
                                     std::size_t data_size,
                                     int load_data,
                                     int autopartition,
                                     int broadcasthead,
                                     void *user_pointer,
                                     types<3>::connectivity **p4est);
#else
      static
      types<3>::forest *(&load) (const char *filename,
                                 MPI_Comm mpicomm,
                                 std::size_t data_size,
                                 int load_data,
                                 void *user_pointer,
                                 types<3>::connectivity **p4est);
#endif

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      int (&connectivity_save) (const char *filename,
                                types<3>::connectivity *connectivity);
#else
      static
      void (&connectivity_save) (const char *filename,
                                 types<3>::connectivity *connectivity);
#endif

      static
      int (&connectivity_is_valid) (types<3>::connectivity *connectivity);

#if DEAL_II_P4EST_VERSION_GTE(1,0,0,0)
      static
      types<3>::connectivity *(&connectivity_load) (const char *filename,
                                                    size_t *length);
#elif DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      static
      types<3>::connectivity *(&connectivity_load) (const char *filename,
                                                    long unsigned *length);
#else
      static
      types<3>::connectivity *(&connectivity_load) (const char *filename,
                                                    long *length);
#endif

      static
      unsigned int (&checksum) (types<3>::forest *p8est);

      static
      void (&vtk_write_file) (types<3>::forest *p8est,
                              p8est_geometry_t *,
                              const char *baseName);

      static
      types<3>::ghost *(&ghost_new) (types<3>::forest *p4est,
                                     types<3>::balance_type btype);

      static
      void (&ghost_destroy) (types<3>::ghost *ghost);

      static
      void (&reset_data) (types<3>::forest *p4est,
                          size_t data_size,
                          p8est_init_t init_fn,
                          void *user_pointer);

      static
      size_t (&forest_memory_used) (types<3>::forest *p4est);

      static
      size_t (&connectivity_memory_used) (types<3>::connectivity *p4est);

      static const unsigned max_level;
    };


    int (&functions<3>::quadrant_compare) (const void *v1, const void *v2)
      = p8est_quadrant_compare;

    void (&functions<3>::quadrant_childrenv) (const types<3>::quadrant *q,
                                              types<3>::quadrant c[])
      = p8est_quadrant_childrenv;

    int (&functions<3>::quadrant_overlaps_tree) (types<3>::tree *tree,
                                                 const types<3>::quadrant *q)
      = p8est_quadrant_overlaps_tree;

    void (&functions<3>::quadrant_set_morton) (types<3>::quadrant *quadrant,
                                               int level,
                                               uint64_t id)
      = p8est_quadrant_set_morton;

    int (&functions<3>::quadrant_is_equal) (const types<3>::quadrant *q1,
                                            const types<3>::quadrant *q2)
      = p8est_quadrant_is_equal;

    int (&functions<3>::quadrant_is_sibling) (const types<3>::quadrant *q1,
                                              const types<3>::quadrant *q2)
      = p8est_quadrant_is_sibling;

    int (&functions<3>::quadrant_is_ancestor) (const types<3>::quadrant *q1,
                                               const types<3>::quadrant *q2)
      = p8est_quadrant_is_ancestor;

    int (&functions<3>::quadrant_ancestor_id) (const types<3>::quadrant *q,
                                               int level)
      = p8est_quadrant_ancestor_id;

    int (&functions<3>::comm_find_owner) (types<3>::forest *p4est,
                                          const types<3>::locidx which_tree,
                                          const types<3>::quadrant *q,
                                          const int guess)
      = p8est_comm_find_owner;

    types<3>::connectivity *(&functions<3>::connectivity_new) (types<3>::topidx num_vertices,
                                                               types<3>::topidx num_trees,
                                                               types<3>::topidx num_edges,
                                                               types<3>::topidx num_ett,
                                                               types<3>::topidx num_corners,
                                                               types<3>::topidx num_ctt)
      = p8est_connectivity_new;

    void (&functions<3>::connectivity_destroy) (p8est_connectivity_t *connectivity)
      = p8est_connectivity_destroy;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,1)
    void (&functions<3>::connectivity_join_faces) (types<3>::connectivity *conn,
                                                   types<3>::topidx tree_left,
                                                   types<3>::topidx tree_right,
                                                   int face_left,
                                                   int face_right,
                                                   int orientation)
      = p8est_connectivity_join_faces;
#endif

    types<3>::forest *(&functions<3>::new_forest) (MPI_Comm mpicomm,
                                                   types<3>::connectivity *connectivity,
                                                   types<3>::locidx min_quadrants,
                                                   int min_level,
                                                   int fill_uniform,
                                                   size_t data_size,
                                                   p8est_init_t init_fn,
                                                   void *user_pointer)
      = p8est_new_ext;

    void (&functions<3>::destroy) (types<3>::forest *p8est)
      = p8est_destroy;

    void (&functions<3>::refine) (types<3>::forest *p8est,
                                  int refine_recursive,
                                  p8est_refine_t refine_fn,
                                  p8est_init_t init_fn)
      = p8est_refine;

    void (&functions<3>::coarsen) (types<3>::forest *p8est,
                                   int coarsen_recursive,
                                   p8est_coarsen_t coarsen_fn,
                                   p8est_init_t init_fn)
      = p8est_coarsen;

    void (&functions<3>::balance) (types<3>::forest *p8est,
                                   types<3>::balance_type btype,
                                   p8est_init_t init_fn)
      = p8est_balance;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    p4est_gloidx_t (&functions<3>::partition) (types<3>::forest *p8est,
                                               int partition_for_coarsening,
                                               p8est_weight_t weight_fn)
      = p8est_partition_ext;
#else
    void (&functions<3>::partition) (types<3>::forest *p8est,
                                     int partition_for_coarsening,
                                     p8est_weight_t weight_fn)
      = p8est_partition_ext;
#endif

    void (&functions<3>::save) (const char *filename,
                                types<3>::forest *p4est,
                                int save_data)
      = p8est_save;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    types<3>::forest *
    (&functions<3>::load_ext) (const char *filename,
                               MPI_Comm mpicomm,
                               std::size_t data_size,
                               int load_data,
                               int autopartition,
                               int broadcasthead,
                               void *user_pointer,
                               types<3>::connectivity **p4est)
      = p8est_load_ext;
#else
    types<3>::forest *
    (&functions<3>::load) (const char *filename,
                           MPI_Comm mpicomm,
                           std::size_t data_size,
                           int load_data,
                           void *user_pointer,
                           types<3>::connectivity **p4est)
      = p8est_load;
#endif

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    int (&functions<3>::connectivity_save) (const char *filename,
                                            types<3>::connectivity *connectivity)
      = p8est_connectivity_save;
#else
    void (&functions<3>::connectivity_save) (const char *filename,
                                             types<3>::connectivity *connectivity)
      = p8est_connectivity_save;
#endif

    int (&functions<3>::connectivity_is_valid) (types<3>::connectivity
                                                *connectivity)
      = p8est_connectivity_is_valid;

#if DEAL_II_P4EST_VERSION_GTE(1,0,0,0)
    types<3>::connectivity *
    (&functions<3>::connectivity_load) (const char *filename,
                                        size_t *length)
      = p8est_connectivity_load;
#elif DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
    types<3>::connectivity *
    (&functions<3>::connectivity_load) (const char *filename,
                                        long unsigned *length)
      = p8est_connectivity_load;
#else
    types<3>::connectivity *
    (&functions<3>::connectivity_load) (const char *filename,
                                        long *length)
      = p8est_connectivity_load;
#endif

    unsigned int (&functions<3>::checksum) (types<3>::forest *p8est)
      = p8est_checksum;

    void (&functions<3>::vtk_write_file) (types<3>::forest *p8est,
                                          p8est_geometry_t *,
                                          const char *baseName)
      = p8est_vtk_write_file;

    types<3>::ghost *(&functions<3>::ghost_new) (types<3>::forest *p4est,
                                                 types<3>::balance_type btype)
      = p8est_ghost_new;

    void (&functions<3>::ghost_destroy) (types<3>::ghost *ghost)
      = p8est_ghost_destroy;

    void (&functions<3>::reset_data) (types<3>::forest *p4est,
                                      size_t data_size,
                                      p8est_init_t init_fn,
                                      void *user_pointer)
      = p8est_reset_data;

    size_t (&functions<3>::forest_memory_used) (types<3>::forest *p4est)
      = p8est_memory_used;

    size_t (&functions<3>::connectivity_memory_used) (types<3>::connectivity *p4est)
      = p8est_connectivity_memory_used;

    const unsigned int functions<3>::max_level = P8EST_MAXLEVEL;


    template <int dim>
    void
    init_quadrant_children
    (const typename types<dim>::quadrant &p4est_cell,
     typename types<dim>::quadrant (&p4est_children)[GeometryInfo<dim>::max_children_per_cell])
    {

      for (unsigned int c=0;
           c<GeometryInfo<dim>::max_children_per_cell; ++c)
        switch (dim)
          {
          case 2:
            P4EST_QUADRANT_INIT(&p4est_children[c]);
            break;
          case 3:
            P8EST_QUADRANT_INIT(&p4est_children[c]);
            break;
          default:
            Assert (false, ExcNotImplemented());
          }


      functions<dim>::quadrant_childrenv (&p4est_cell,
                                          p4est_children);

    }


    template <int dim>
    void
    init_coarse_quadrant(typename types<dim>::quadrant &quad)
    {
      switch (dim)
        {
        case 2:
          P4EST_QUADRANT_INIT(&quad);
          break;
        case 3:
          P8EST_QUADRANT_INIT(&quad);
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
      functions<dim>::quadrant_set_morton (&quad,
                                           /*level=*/0,
                                           /*index=*/0);
    }


    template <int dim>
    bool
    quadrant_is_equal (const typename types<dim>::quadrant &q1,
                       const typename types<dim>::quadrant &q2)
    {
      return functions<dim>::quadrant_is_equal(&q1, &q2);
    }



    template <int dim>
    bool
    quadrant_is_ancestor (const typename types<dim>::quadrant &q1,
                          const typename types<dim>::quadrant &q2)
    {
      return functions<dim>::quadrant_is_ancestor(&q1, &q2);
    }

    /**
     * This struct templatizes the p4est iterate structs and function
     * prototypes, which are used to execute callback functions for faces,
     * edges, and corners that require local neighborhood information, i.e.
     * the neighboring cells */
    template <int dim> struct iter;

    template <> struct iter<2>
    {
      typedef p4est_iter_corner_info_t corner_info;
      typedef p4est_iter_corner_side_t corner_side;
      typedef p4est_iter_corner_t      corner_iter;
      typedef p4est_iter_face_info_t face_info;
      typedef p4est_iter_face_side_t face_side;
      typedef p4est_iter_face_t      face_iter;
    };

    template <> struct iter<3>
    {
      typedef p8est_iter_corner_info_t corner_info;
      typedef p8est_iter_corner_side_t corner_side;
      typedef p8est_iter_corner_t      corner_iter;
      typedef p8est_iter_edge_info_t edge_info;
      typedef p8est_iter_edge_side_t edge_side;
      typedef p8est_iter_edge_t      edge_iter;
      typedef p8est_iter_face_info_t face_info;
      typedef p8est_iter_face_side_t face_side;
      typedef p8est_iter_face_t      face_iter;
    };

  }
}


namespace
{
  template <int dim, int spacedim>
  void
  get_vertex_to_cell_mappings (const Triangulation<dim,spacedim> &triangulation,
                               std::vector<unsigned int> &vertex_touch_count,
                               std::vector<std::list<
                               std::pair<typename Triangulation<dim,spacedim>::active_cell_iterator,unsigned int> > >
                               &vertex_to_cell)
  {
    vertex_touch_count.resize (triangulation.n_vertices());
    vertex_to_cell.resize (triangulation.n_vertices());

    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          ++vertex_touch_count[cell->vertex_index(v)];
          vertex_to_cell[cell->vertex_index(v)]
          .push_back (std::make_pair (cell, v));
        }
  }



  template <int dim, int spacedim>
  void
  get_edge_to_cell_mappings (const Triangulation<dim,spacedim> &triangulation,
                             std::vector<unsigned int> &edge_touch_count,
                             std::vector<std::list<
                             std::pair<typename Triangulation<dim,spacedim>::active_cell_iterator,unsigned int> > >
                             &edge_to_cell)
  {
    Assert (triangulation.n_levels() == 1, ExcInternalError());

    edge_touch_count.resize (triangulation.n_active_lines());
    edge_to_cell.resize (triangulation.n_active_lines());

    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        {
          ++edge_touch_count[cell->line(l)->index()];
          edge_to_cell[cell->line(l)->index()]
          .push_back (std::make_pair (cell, l));
        }
  }



  /**
   * Set all vertex and cell related information in the p4est connectivity
   * structure.
   */
  template <int dim, int spacedim>
  void
  set_vertex_and_cell_info (const Triangulation<dim,spacedim> &triangulation,
                            const std::vector<unsigned int> &vertex_touch_count,
                            const std::vector<std::list<
                            std::pair<typename Triangulation<dim,spacedim>::active_cell_iterator,unsigned int> > >
                            &vertex_to_cell,
                            const std::vector<types::global_dof_index> &coarse_cell_to_p4est_tree_permutation,
                            const bool set_vertex_info,
                            typename internal::p4est::types<dim>::connectivity *connectivity)
  {
    // copy the vertices into the connectivity structure. the triangulation
    // exports the array of vertices, but some of the entries are sometimes
    // unused; this shouldn't be the case for a newly created triangulation,
    // but make sure
    //
    // note that p4est stores coordinates as a triplet of values even in 2d
    Assert (triangulation.get_used_vertices().size() ==
            triangulation.get_vertices().size(),
            ExcInternalError());
    Assert (std::find (triangulation.get_used_vertices().begin(),
                       triangulation.get_used_vertices().end(),
                       false)
            == triangulation.get_used_vertices().end(),
            ExcInternalError());
    if (set_vertex_info == true)
      for (unsigned int v=0; v<triangulation.n_vertices(); ++v)
        {
          connectivity->vertices[3*v  ] = triangulation.get_vertices()[v][0];
          connectivity->vertices[3*v+1] = triangulation.get_vertices()[v][1];
          connectivity->vertices[3*v+2] = (spacedim == 2 ?
                                           0
                                           :
                                           triangulation.get_vertices()[v][2]);
        }

    // next store the tree_to_vertex indices (each tree is here only a single
    // cell in the coarse mesh). p4est requires vertex numbering in clockwise
    // orientation
    //
    // while we're at it, also copy the neighborship information between cells
    typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    for (; cell != endc; ++cell)
      {
        const unsigned int
        index = coarse_cell_to_p4est_tree_permutation[cell->index()];

        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            if (set_vertex_info == true)
              connectivity->tree_to_vertex[index*GeometryInfo<dim>::vertices_per_cell+v] = cell->vertex_index(v);
            connectivity->tree_to_corner[index*GeometryInfo<dim>::vertices_per_cell+v] = cell->vertex_index(v);
          }

        // neighborship information. if a cell is at a boundary, then enter
        // the index of the cell itself here
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary() == false)
            connectivity->tree_to_tree[index*GeometryInfo<dim>::faces_per_cell + f]
              = coarse_cell_to_p4est_tree_permutation[cell->neighbor(f)->index()];
          else
            connectivity->tree_to_tree[index*GeometryInfo<dim>::faces_per_cell + f]
              = coarse_cell_to_p4est_tree_permutation[cell->index()];

        // fill tree_to_face, which is essentially neighbor_to_neighbor;
        // however, we have to remap the resulting face number as well
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary() == false)
            {
              switch (dim)
                {
                case 2:
                {
                  connectivity->tree_to_face[index*GeometryInfo<dim>::faces_per_cell + f]
                    = cell->neighbor_of_neighbor (f);
                  break;
                }

                case 3:
                {
                  /*
                   * The values for tree_to_face are in 0..23 where ttf % 6
                   * gives the face number and ttf / 4 the face orientation
                   * code.  The orientation is determined as follows.  Let
                   * my_face and other_face be the two face numbers of the
                   * connecting trees in 0..5.  Then the first face vertex of
                   * the lower of my_face and other_face connects to a face
                   * vertex numbered 0..3 in the higher of my_face and
                   * other_face.  The face orientation is defined as this
                   * number.  If my_face == other_face, treating either of
                   * both faces as the lower one leads to the same result.
                   */

                  connectivity->tree_to_face[index*6 + f]
                    = cell->neighbor_of_neighbor (f);

                  unsigned int face_idx_list[2] =
                  {f, cell->neighbor_of_neighbor (f)};
                  typename Triangulation<dim>::active_cell_iterator
                  cell_list[2] = {cell, cell->neighbor(f)};
                  unsigned int smaller_idx = 0;

                  if (f>cell->neighbor_of_neighbor (f))
                    smaller_idx = 1;

                  unsigned int larger_idx = (smaller_idx+1) % 2;
                  //smaller = *_list[smaller_idx]
                  //larger = *_list[larger_idx]

                  unsigned int v = 0;

                  // global vertex index of vertex 0 on face of cell with
                  // smaller local face index
                  unsigned int g_idx =
                    cell_list[smaller_idx]->vertex_index(
                      GeometryInfo<dim>::face_to_cell_vertices(
                        face_idx_list[smaller_idx],
                        0,
                        cell_list[smaller_idx]->face_orientation(face_idx_list[smaller_idx]),
                        cell_list[smaller_idx]->face_flip(face_idx_list[smaller_idx]),
                        cell_list[smaller_idx]->face_rotation(face_idx_list[smaller_idx]))
                    );

                  // loop over vertices on face from other cell and compare
                  // global vertex numbers
                  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
                    {
                      unsigned int idx
                        =
                          cell_list[larger_idx]->vertex_index(
                            GeometryInfo<dim>::face_to_cell_vertices(
                              face_idx_list[larger_idx],
                              i)
                          );

                      if (idx==g_idx)
                        {
                          v = i;
                          break;
                        }
                    }

                  connectivity->tree_to_face[index*6 + f] += 6*v;
                  break;
                }

                default:
                  Assert (false, ExcNotImplemented());
                }
            }
          else
            connectivity->tree_to_face[index*GeometryInfo<dim>::faces_per_cell + f] = f;
      }

    // now fill the vertex information
    connectivity->ctt_offset[0] = 0;
    std::partial_sum (vertex_touch_count.begin(),
                      vertex_touch_count.end(),
                      &connectivity->ctt_offset[1]);

    const typename internal::p4est::types<dim>::locidx
    num_vtt = std::accumulate (vertex_touch_count.begin(),
                               vertex_touch_count.end(),
                               0);
    (void)num_vtt;
    Assert (connectivity->ctt_offset[triangulation.n_vertices()] ==
            num_vtt,
            ExcInternalError());

    for (unsigned int v=0; v<triangulation.n_vertices(); ++v)
      {
        Assert (vertex_to_cell[v].size() == vertex_touch_count[v],
                ExcInternalError());

        typename std::list<std::pair
        <typename Triangulation<dim,spacedim>::active_cell_iterator,
        unsigned int> >::const_iterator
        p = vertex_to_cell[v].begin();
        for (unsigned int c=0; c<vertex_touch_count[v]; ++c, ++p)
          {
            connectivity->corner_to_tree[connectivity->ctt_offset[v]+c]
              = coarse_cell_to_p4est_tree_permutation[p->first->index()];
            connectivity->corner_to_corner[connectivity->ctt_offset[v]+c]
              = p->second;
          }
      }
  }



  template <int dim, int spacedim>
  bool
  tree_exists_locally (const typename internal::p4est::types<dim>::forest *parallel_forest,
                       const typename internal::p4est::types<dim>::topidx coarse_grid_cell)
  {
    Assert (coarse_grid_cell < parallel_forest->connectivity->num_trees,
            ExcInternalError());
    return ((coarse_grid_cell >= parallel_forest->first_local_tree)
            &&
            (coarse_grid_cell <= parallel_forest->last_local_tree));
  }


  template <int dim, int spacedim>
  void
  delete_all_children_and_self (const typename Triangulation<dim,spacedim>::cell_iterator &cell)
  {
    if (cell->has_children())
      for (unsigned int c=0; c<cell->n_children(); ++c)
        delete_all_children_and_self<dim,spacedim> (cell->child(c));
    else
      cell->set_coarsen_flag ();
  }



  template <int dim, int spacedim>
  void
  delete_all_children (const typename Triangulation<dim,spacedim>::cell_iterator &cell)
  {
    if (cell->has_children())
      for (unsigned int c=0; c<cell->n_children(); ++c)
        delete_all_children_and_self<dim,spacedim> (cell->child(c));
  }


  template <int dim, int spacedim>
  void
  determine_level_subdomain_id_recursively (const typename internal::p4est::types<dim>::tree     &tree,
                                            const typename internal::p4est::types<dim>::locidx &tree_index,
                                            const typename Triangulation<dim,spacedim>::cell_iterator     &dealii_cell,
                                            const typename internal::p4est::types<dim>::quadrant &p4est_cell,
                                            typename internal::p4est::types<dim>::forest   &forest,
                                            const types::subdomain_id                           my_subdomain,
                                            const std::vector<std::vector<bool> > &marked_vertices)
  {
    if (dealii_cell->level_subdomain_id()==numbers::artificial_subdomain_id)
      {
        //important: only assign the level_subdomain_id if it is a ghost cell
        // even though we could fill in all.
        bool used = false;
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            if (marked_vertices[dealii_cell->level()][dealii_cell->vertex_index(v)])
              {
                used = true;
                break;
              }
          }

        // Special case: if this cell is active we might be a ghost neighbor
        // to a locally owned cell across a vertex that is finer.
        // Example (M= my, O=dealii_cell, owned by somebody else):
        //         *------*
        //         |      |
        //         |  O   |
        //         |      |
        // *---*---*------*
        // | M | M |
        // *---*---*
        // |   | M |
        // *---*---*
        if (!used && dealii_cell->active() && dealii_cell->is_artificial()==false
            && dealii_cell->level()+1<static_cast<int>(marked_vertices.size()))
          {
            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
              {
                if (marked_vertices[dealii_cell->level()+1][dealii_cell->vertex_index(v)])
                  {
                    used = true;
                    break;
                  }
              }
          }

        // Like above, but now the other way around
        if (!used && dealii_cell->active() && dealii_cell->is_artificial()==false
            && dealii_cell->level()>0)
          {
            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
              {
                if (marked_vertices[dealii_cell->level()-1][dealii_cell->vertex_index(v)])
                  {
                    used = true;
                    break;
                  }
              }
          }

        if (used)
          {
            int owner = internal::p4est::functions<dim>::comm_find_owner (&forest,
                        tree_index,
                        &p4est_cell,
                        my_subdomain);
            Assert((owner!=-2) && (owner!=-1), ExcMessage("p4est should know the owner."));
            dealii_cell->set_level_subdomain_id(owner);
          }

      }

    if (dealii_cell->has_children ())
      {
        typename internal::p4est::types<dim>::quadrant
        p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c=0;
             c<GeometryInfo<dim>::max_children_per_cell; ++c)
          switch (dim)
            {
            case 2:
              P4EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            case 3:
              P8EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            default:
              Assert (false, ExcNotImplemented());
            }


        internal::p4est::functions<dim>::
        quadrant_childrenv (&p4est_cell,
                            p4est_child);

        for (unsigned int c=0;
             c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            determine_level_subdomain_id_recursively <dim,spacedim> (tree,tree_index,
                                                                     dealii_cell->child(c),
                                                                     p4est_child[c],
                                                                     forest,
                                                                     my_subdomain,
                                                                     marked_vertices);
          }
      }
  }


  template <int dim, int spacedim>
  void
  match_tree_recursively (const typename internal::p4est::types<dim>::tree     &tree,
                          const typename Triangulation<dim,spacedim>::cell_iterator     &dealii_cell,
                          const typename internal::p4est::types<dim>::quadrant &p4est_cell,
                          const typename internal::p4est::types<dim>::forest   &forest,
                          const types::subdomain_id                           my_subdomain)
  {
    // check if this cell exists in the local p4est cell
    if (sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                         &p4est_cell,
                         internal::p4est::functions<dim>::quadrant_compare)
        != -1)
      {
        // yes, cell found in local part of p4est
        delete_all_children<dim,spacedim> (dealii_cell);
        if (!dealii_cell->has_children())
          dealii_cell->set_subdomain_id(my_subdomain);
      }
    else
      {
        // no, cell not found in local part of p4est. this means that the
        // local part is more refined than the current cell. if this cell has
        // no children of its own, we need to refine it, and if it does
        // already have children then loop over all children and see if they
        // are locally available as well
        if (dealii_cell->has_children () == false)
          dealii_cell->set_refine_flag ();
        else
          {
            typename internal::p4est::types<dim>::quadrant
            p4est_child[GeometryInfo<dim>::max_children_per_cell];
            for (unsigned int c=0;
                 c<GeometryInfo<dim>::max_children_per_cell; ++c)
              switch (dim)
                {
                case 2:
                  P4EST_QUADRANT_INIT(&p4est_child[c]);
                  break;
                case 3:
                  P8EST_QUADRANT_INIT(&p4est_child[c]);
                  break;
                default:
                  Assert (false, ExcNotImplemented());
                }


            internal::p4est::functions<dim>::
            quadrant_childrenv (&p4est_cell,
                                p4est_child);

            for (unsigned int c=0;
                 c<GeometryInfo<dim>::max_children_per_cell; ++c)
              if (internal::p4est::functions<dim>::
                  quadrant_overlaps_tree (const_cast<typename internal::p4est::types<dim>::tree *>(&tree),
                                          &p4est_child[c])
                  == false)
                {
                  // no, this child is locally not available in the p4est.
                  // delete all its children but, because this may not be
                  // successfull, make sure to mark all children recursively
                  // as not local.
                  delete_all_children<dim,spacedim> (dealii_cell->child(c));
                  dealii_cell->child(c)
                  ->recursively_set_subdomain_id(numbers::artificial_subdomain_id);
                }
              else
                {
                  // at least some part of the tree rooted in this child is
                  // locally available
                  match_tree_recursively<dim,spacedim> (tree,
                                                        dealii_cell->child(c),
                                                        p4est_child[c],
                                                        forest,
                                                        my_subdomain);
                }
          }
      }
  }


  template <int dim, int spacedim>
  void
  match_quadrant (const dealii::Triangulation<dim,spacedim> *tria,
                  unsigned int dealii_index,
                  typename internal::p4est::types<dim>::quadrant &ghost_quadrant,
                  unsigned int ghost_owner)
  {
    int i, child_id;
    int l = ghost_quadrant.level;

    for (i = 0; i < l; i++)
      {
        typename Triangulation<dim,spacedim>::cell_iterator cell (tria, i, dealii_index);
        if (cell->has_children () == false)
          {
            cell->clear_coarsen_flag();
            cell->set_refine_flag ();
            return;
          }

        child_id = internal::p4est::functions<dim>::quadrant_ancestor_id (&ghost_quadrant, i + 1);
        dealii_index = cell->child_index(child_id);
      }

    typename Triangulation<dim,spacedim>::cell_iterator cell (tria, l, dealii_index);
    if (cell->has_children())
      delete_all_children<dim,spacedim> (cell);
    else
      {
        cell->clear_coarsen_flag();
        cell->set_subdomain_id(ghost_owner);
      }
  }



  template <int dim, int spacedim>
  void
  attach_mesh_data_recursively (const typename internal::p4est::types<dim>::tree &tree,
                                const typename Triangulation<dim,spacedim>::cell_iterator &dealii_cell,
                                const typename internal::p4est::types<dim>::quadrant &p4est_cell,
                                const typename std::list<std::pair<unsigned int, typename std_cxx11::function<
                                void(typename parallel::distributed::Triangulation<dim,spacedim>::cell_iterator,
                                     typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus,
                                     void *)
                                > > > &attached_data_pack_callbacks)
  {
    typedef std::list<std::pair<unsigned int, typename std_cxx11::function<
    void(typename parallel::distributed::Triangulation<dim,spacedim>::cell_iterator,
         typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus,
         void *)
    > > > callback_list_t;

    int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                               &p4est_cell,
                               internal::p4est::functions<dim>::quadrant_compare);

    if (idx == -1 && (internal::p4est::functions<dim>::
                      quadrant_overlaps_tree (const_cast<typename internal::p4est::types<dim>::tree *>(&tree),
                                              &p4est_cell)
                      == false))
      return; //this quadrant and none of its children belongs to us.

    bool p4est_has_children = (idx == -1);

    if (p4est_has_children && dealii_cell->has_children())
      {
        //recurse further
        typename internal::p4est::types<dim>::quadrant
        p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          switch (dim)
            {
            case 2:
              P4EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            case 3:
              P8EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            default:
              Assert (false, ExcNotImplemented());
            }

        internal::p4est::functions<dim>::
        quadrant_childrenv (&p4est_cell, p4est_child);

        for (unsigned int c=0;
             c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            attach_mesh_data_recursively<dim,spacedim> (tree,
                                                        dealii_cell->child(c),
                                                        p4est_child[c],
                                                        attached_data_pack_callbacks);
          }
      }
    else if (!p4est_has_children && !dealii_cell->has_children())
      {
        //this active cell didn't change
        typename internal::p4est::types<dim>::quadrant *q;
        q = static_cast<typename internal::p4est::types<dim>::quadrant *> (
              sc_array_index (const_cast<sc_array_t *>(&tree.quadrants), idx)
            );
        *static_cast<typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus *>(q->p.user_data) = parallel::distributed::Triangulation<dim,spacedim>::CELL_PERSIST;

        for (typename callback_list_t::const_iterator it = attached_data_pack_callbacks.begin();
             it != attached_data_pack_callbacks.end();
             ++it)
          {
            void *ptr = static_cast<char *>(q->p.user_data) + (*it).first; //add offset
            ((*it).second)(dealii_cell,
                           parallel::distributed::Triangulation<dim,spacedim>::CELL_PERSIST,
                           ptr);
          }
      }
    else if (p4est_has_children)
      {
        //this cell got refined

        //attach to the first child, because we can only attach to active
        // quadrants
        typename internal::p4est::types<dim>::quadrant
        p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          switch (dim)
            {
            case 2:
              P4EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            case 3:
              P8EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            default:
              Assert (false, ExcNotImplemented());
            }

        internal::p4est::functions<dim>::
        quadrant_childrenv (&p4est_cell, p4est_child);
        int child0_idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                                          &p4est_child[0],
                                          internal::p4est::functions<dim>::quadrant_compare);
        Assert(child0_idx != -1, ExcMessage("the first child should exist as an active quadrant!"));

        typename internal::p4est::types<dim>::quadrant *q;
        q = static_cast<typename internal::p4est::types<dim>::quadrant *> (
              sc_array_index (const_cast<sc_array_t *>(&tree.quadrants), child0_idx)
            );
        *static_cast<typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus *>(q->p.user_data) = parallel::distributed::Triangulation<dim,spacedim>::CELL_REFINE;

        for (typename callback_list_t::const_iterator it = attached_data_pack_callbacks.begin();
             it != attached_data_pack_callbacks.end();
             ++it)
          {
            void *ptr = static_cast<char *>(q->p.user_data) + (*it).first; //add offset

            ((*it).second)(dealii_cell,
                           parallel::distributed::Triangulation<dim,spacedim>::CELL_REFINE,
                           ptr);
          }

        //mark other children as invalid, so that unpack only happens once
        for (unsigned int i=1; i<GeometryInfo<dim>::max_children_per_cell; ++i)
          {
            int child_idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                                             &p4est_child[i],
                                             internal::p4est::functions<dim>::quadrant_compare);
            q = static_cast<typename internal::p4est::types<dim>::quadrant *> (
                  sc_array_index (const_cast<sc_array_t *>(&tree.quadrants), child_idx)
                );
            *static_cast<typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus *>(q->p.user_data) = parallel::distributed::Triangulation<dim,spacedim>::CELL_INVALID;
          }


      }
    else
      {
        //its children got coarsened into this cell
        typename internal::p4est::types<dim>::quadrant *q;
        q = static_cast<typename internal::p4est::types<dim>::quadrant *> (
              sc_array_index (const_cast<sc_array_t *>(&tree.quadrants), idx)
            );
        *static_cast<typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus *>(q->p.user_data) = parallel::distributed::Triangulation<dim,spacedim>::CELL_COARSEN;

        for (typename callback_list_t::const_iterator it = attached_data_pack_callbacks.begin();
             it != attached_data_pack_callbacks.end();
             ++it)
          {
            void *ptr = static_cast<char *>(q->p.user_data) + (*it).first; //add offset
            ((*it).second)(dealii_cell,
                           parallel::distributed::Triangulation<dim,spacedim>::CELL_COARSEN,
                           ptr);
          }
      }
  }

  template <int dim, int spacedim>
  void
  get_cell_weights_recursively (const typename internal::p4est::types<dim>::tree &tree,
                                const typename Triangulation<dim,spacedim>::cell_iterator &dealii_cell,
                                const typename internal::p4est::types<dim>::quadrant &p4est_cell,
                                const typename Triangulation<dim,spacedim>::Signals &signals,
                                std::vector<unsigned int> &weight)
  {
    const int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                                     &p4est_cell,
                                     internal::p4est::functions<dim>::quadrant_compare);

    if (idx == -1 && (internal::p4est::functions<dim>::
                      quadrant_overlaps_tree (const_cast<typename internal::p4est::types<dim>::tree *>(&tree),
                                              &p4est_cell)
                      == false))
      return; // This quadrant and none of its children belongs to us.

    const bool p4est_has_children = (idx == -1);

    if (p4est_has_children && dealii_cell->has_children())
      {
        //recurse further
        typename internal::p4est::types<dim>::quadrant
        p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          switch (dim)
            {
            case 2:
              P4EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            case 3:
              P8EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            default:
              Assert (false, ExcNotImplemented());
            }

        internal::p4est::functions<dim>::
        quadrant_childrenv (&p4est_cell, p4est_child);

        for (unsigned int c=0;
             c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            get_cell_weights_recursively<dim,spacedim> (tree,
                                                        dealii_cell->child(c),
                                                        p4est_child[c],
                                                        signals,
                                                        weight);
          }
      }
    else if (!p4est_has_children && !dealii_cell->has_children())
      {
        // This active cell didn't change
        weight.push_back(1000);
        weight.back() += signals.cell_weight(dealii_cell,
                                             parallel::distributed::Triangulation<dim,spacedim>::CELL_PERSIST);
      }
    else if (p4est_has_children)
      {
        // This cell will be refined
        unsigned int parent_weight(1000);
        parent_weight += signals.cell_weight(dealii_cell,
                                             parallel::distributed::Triangulation<dim,spacedim>::CELL_REFINE);

        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            // We assign the weight of the parent cell equally to all children
            weight.push_back(parent_weight);
          }
      }
    else
      {
        // This cell's children will be coarsened into this cell
        weight.push_back(1000);
        weight.back() += signals.cell_weight(dealii_cell,
                                             parallel::distributed::Triangulation<dim,spacedim>::CELL_COARSEN);
      }
  }


  template <int dim, int spacedim>
  void
  post_mesh_data_recursively (const typename internal::p4est::types<dim>::tree &tree,
                              const typename Triangulation<dim,spacedim>::cell_iterator &dealii_cell,
                              const typename Triangulation<dim,spacedim>::cell_iterator &parent_cell,
                              const typename internal::p4est::types<dim>::quadrant &p4est_cell,
                              const unsigned int offset,
                              const typename std_cxx11::function<
                              void(typename parallel::distributed::Triangulation<dim,spacedim>::cell_iterator, typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus, void *)
                              > &unpack_callback)
  {
    int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                               &p4est_cell,
                               internal::p4est::functions<dim>::quadrant_compare);
    if (idx == -1 && (internal::p4est::functions<dim>::
                      quadrant_overlaps_tree (const_cast<typename internal::p4est::types<dim>::tree *>(&tree),
                                              &p4est_cell)
                      == false))
      // this quadrant and none of its children belong to us.
      return;


    const bool p4est_has_children = (idx == -1);
    if (p4est_has_children)
      {
        Assert(dealii_cell->has_children(), ExcInternalError());

        //recurse further
        typename internal::p4est::types<dim>::quadrant
        p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          switch (dim)
            {
            case 2:
              P4EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            case 3:
              P8EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            default:
              Assert (false, ExcNotImplemented());
            }

        internal::p4est::functions<dim>::
        quadrant_childrenv (&p4est_cell, p4est_child);

        for (unsigned int c=0;
             c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            post_mesh_data_recursively<dim,spacedim> (tree,
                                                      dealii_cell->child(c),
                                                      dealii_cell,
                                                      p4est_child[c],
                                                      offset,
                                                      unpack_callback);
          }
      }
    else
      {
        Assert(! dealii_cell->has_children(), ExcInternalError());

        typename internal::p4est::types<dim>::quadrant *q;
        q = static_cast<typename internal::p4est::types<dim>::quadrant *> (
              sc_array_index (const_cast<sc_array_t *>(&tree.quadrants), idx)
            );

        void *ptr = static_cast<char *>(q->p.user_data) + offset;
        typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus
        status = * static_cast<
                 typename parallel::distributed::Triangulation<dim,spacedim>::CellStatus *
                 >(q->p.user_data);
        switch (status)
          {
          case parallel::distributed::Triangulation<dim,spacedim>::CELL_PERSIST:
          {
            unpack_callback(dealii_cell, status, ptr);
            break;
          }
          case parallel::distributed::Triangulation<dim,spacedim>::CELL_REFINE:
          {
            unpack_callback(parent_cell, status, ptr);
            break;
          }
          case parallel::distributed::Triangulation<dim,spacedim>::CELL_COARSEN:
          {
            unpack_callback(dealii_cell, status, ptr);
            break;
          }
          case parallel::distributed::Triangulation<dim,spacedim>::CELL_INVALID:
          {
            break;
          }
          default:
            AssertThrow (false, ExcInternalError());
          }
      }
  }



  /**
   * A data structure that we use to store which cells (indicated by
   * internal::p4est::types<dim>::quadrant objects) shall be refined and which
   * shall be coarsened.
   */
  template <int dim, int spacedim>
  class RefineAndCoarsenList
  {
  public:
    RefineAndCoarsenList (const Triangulation<dim,spacedim> &triangulation,
                          const std::vector<types::global_dof_index> &p4est_tree_to_coarse_cell_permutation,
                          const types::subdomain_id                   my_subdomain);

    /**
     * A callback function that we pass to the p4est data structures when a
     * forest is to be refined. The p4est functions call it back with a tree
     * (the index of the tree that grows out of a given coarse cell) and a
     * refinement path from that coarse cell to a terminal/leaf cell. The
     * function returns whether the corresponding cell in the deal.II
     * triangulation has the refined flag set.
     */
    static
    int
    refine_callback (typename internal::p4est::types<dim>::forest *forest,
                     typename internal::p4est::types<dim>::topidx  coarse_cell_index,
                     typename internal::p4est::types<dim>::quadrant *quadrant);

    /**
     * Same as the refine_callback function, but return whether all four of
     * the given children of a non-terminal cell are to be coarsened away.
     */
    static
    int
    coarsen_callback (typename internal::p4est::types<dim>::forest *forest,
                      typename internal::p4est::types<dim>::topidx  coarse_cell_index,
                      typename internal::p4est::types<dim>::quadrant *children[]);

    bool pointers_are_at_end () const;

  private:
    std::vector<typename internal::p4est::types<dim>::quadrant> refine_list;
    typename std::vector<typename internal::p4est::types<dim>::quadrant>::const_iterator current_refine_pointer;

    std::vector<typename internal::p4est::types<dim>::quadrant> coarsen_list;
    typename std::vector<typename internal::p4est::types<dim>::quadrant>::const_iterator current_coarsen_pointer;

    void build_lists (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                      const typename internal::p4est::types<dim>::quadrant &p4est_cell,
                      const unsigned int myid);
  };



  template <int dim, int spacedim>
  bool
  RefineAndCoarsenList<dim,spacedim>::
  pointers_are_at_end () const
  {
    return ((current_refine_pointer == refine_list.end())
            &&
            (current_coarsen_pointer == coarsen_list.end()));
  }



  template <int dim, int spacedim>
  RefineAndCoarsenList<dim,spacedim>::
  RefineAndCoarsenList (const Triangulation<dim,spacedim>            &triangulation,
                        const std::vector<types::global_dof_index>   &p4est_tree_to_coarse_cell_permutation,
                        const types::subdomain_id                    my_subdomain)
  {
    // count how many flags are set and allocate that much memory
    unsigned int n_refine_flags  = 0,
                 n_coarsen_flags = 0;
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
        //skip cells that are not local
        if (cell->subdomain_id() != my_subdomain)
          continue;

        if (cell->refine_flag_set())
          ++n_refine_flags;
        else if (cell->coarsen_flag_set())
          ++n_coarsen_flags;
      }

    refine_list.reserve (n_refine_flags);
    coarsen_list.reserve (n_coarsen_flags);


    // now build the lists of cells that are flagged. note that p4est will
    // traverse its cells in the order in which trees appear in the
    // forest. this order is not the same as the order of coarse cells in the
    // deal.II Triangulation because we have translated everything by the
    // coarse_cell_to_p4est_tree_permutation permutation. in order to make
    // sure that the output array is already in the correct order, traverse
    // our coarse cells in the same order in which p4est will:
    for (unsigned int c=0; c<triangulation.n_cells(0); ++c)
      {
        unsigned int coarse_cell_index =
          p4est_tree_to_coarse_cell_permutation[c];

        const typename Triangulation<dim,spacedim>::cell_iterator
        cell (&triangulation, 0, coarse_cell_index);

        typename internal::p4est::types<dim>::quadrant p4est_cell;
        internal::p4est::functions<dim>::
        quadrant_set_morton (&p4est_cell,
                             /*level=*/0,
                             /*index=*/0);
        p4est_cell.p.which_tree = c;
        build_lists (cell, p4est_cell, my_subdomain);
      }


    Assert(refine_list.size() == n_refine_flags,
           ExcInternalError());
    Assert(coarsen_list.size() == n_coarsen_flags,
           ExcInternalError());

    // make sure that our ordering in fact worked
    for (unsigned int i=1; i<refine_list.size(); ++i)
      Assert (refine_list[i].p.which_tree >=
              refine_list[i-1].p.which_tree,
              ExcInternalError());
    for (unsigned int i=1; i<coarsen_list.size(); ++i)
      Assert (coarsen_list[i].p.which_tree >=
              coarsen_list[i-1].p.which_tree,
              ExcInternalError());

    current_refine_pointer  = refine_list.begin();
    current_coarsen_pointer = coarsen_list.begin();
  }



  template <int dim, int spacedim>
  void
  RefineAndCoarsenList<dim,spacedim>::
  build_lists (const typename Triangulation<dim,spacedim>::cell_iterator     &cell,
               const typename internal::p4est::types<dim>::quadrant &p4est_cell,
               const types::subdomain_id my_subdomain)
  {
    if (!cell->has_children())
      {
        if (cell->subdomain_id() == my_subdomain)
          {
            if (cell->refine_flag_set())
              refine_list.push_back (p4est_cell);
            else if (cell->coarsen_flag_set())
              coarsen_list.push_back (p4est_cell);
          }
      }
    else
      {
        typename internal::p4est::types<dim>::quadrant
        p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          switch (dim)
            {
            case 2:
              P4EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            case 3:
              P8EST_QUADRANT_INIT(&p4est_child[c]);
              break;
            default:
              Assert (false, ExcNotImplemented());
            }
        internal::p4est::functions<dim>::
        quadrant_childrenv (&p4est_cell,
                            p4est_child);
        for (unsigned int c=0;
             c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            p4est_child[c].p.which_tree = p4est_cell.p.which_tree;
            build_lists (cell->child(c),
                         p4est_child[c],
                         my_subdomain);
          }
      }
  }


  template <int dim, int spacedim>
  int
  RefineAndCoarsenList<dim,spacedim>::
  refine_callback (typename internal::p4est::types<dim>::forest *forest,
                   typename internal::p4est::types<dim>::topidx  coarse_cell_index,
                   typename internal::p4est::types<dim>::quadrant *quadrant)
  {
    RefineAndCoarsenList<dim,spacedim> *this_object
      = reinterpret_cast<RefineAndCoarsenList<dim,spacedim>*>(forest->user_pointer);

    // if there are no more cells in our list the current cell can't be
    // flagged for refinement
    if (this_object->current_refine_pointer == this_object->refine_list.end())
      return false;

    Assert (coarse_cell_index <=
            this_object->current_refine_pointer->p.which_tree,
            ExcInternalError());

    // if p4est hasn't yet reached the tree of the next flagged cell the
    // current cell can't be flagged for refinement
    if (coarse_cell_index <
        this_object->current_refine_pointer->p.which_tree)
      return false;

    // now we're in the right tree in the forest
    Assert (coarse_cell_index <=
            this_object->current_refine_pointer->p.which_tree,
            ExcInternalError());

    // make sure that the p4est loop over cells hasn't gotten ahead of our own
    // pointer
    Assert (internal::p4est::functions<dim>::
            quadrant_compare (quadrant,
                              &*this_object->current_refine_pointer) <= 0,
            ExcInternalError());

    // now, if the p4est cell is one in the list, it is supposed to be refined
    if (internal::p4est::functions<dim>::
        quadrant_is_equal (quadrant, &*this_object->current_refine_pointer))
      {
        ++this_object->current_refine_pointer;
        return true;
      }

    // p4est cell is not in list
    return false;
  }



  template <int dim, int spacedim>
  int
  RefineAndCoarsenList<dim,spacedim>::
  coarsen_callback (typename internal::p4est::types<dim>::forest *forest,
                    typename internal::p4est::types<dim>::topidx  coarse_cell_index,
                    typename internal::p4est::types<dim>::quadrant *children[])
  {
    RefineAndCoarsenList<dim,spacedim> *this_object
      = reinterpret_cast<RefineAndCoarsenList<dim,spacedim>*>(forest->user_pointer);

    // if there are no more cells in our list the current cell can't be
    // flagged for coarsening
    if (this_object->current_coarsen_pointer ==
        this_object->coarsen_list.end())
      return false;

    Assert (coarse_cell_index <=
            this_object->current_coarsen_pointer->p.which_tree,
            ExcInternalError());

    // if p4est hasn't yet reached the tree of the next flagged cell the
    // current cell can't be flagged for coarsening
    if (coarse_cell_index <
        this_object->current_coarsen_pointer->p.which_tree)
      return false;

    // now we're in the right tree in the forest
    Assert (coarse_cell_index <=
            this_object->current_coarsen_pointer->p.which_tree,
            ExcInternalError());

    // make sure that the p4est loop over cells hasn't gotten ahead of our own
    // pointer
    Assert (internal::p4est::functions<dim>::
            quadrant_compare (children[0],
                              &*this_object->current_coarsen_pointer) <= 0,
            ExcInternalError());

    // now, if the p4est cell is one in the list, it is supposed to be
    // coarsened
    if (internal::p4est::functions<dim>::
        quadrant_is_equal (children[0],
                           &*this_object->current_coarsen_pointer))
      {
        // move current pointer one up
        ++this_object->current_coarsen_pointer;

        // note that the next 3 cells in our list need to correspond to the
        // other siblings of the cell we have just found
        for (unsigned int c=1; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          {
            Assert (internal::p4est::functions<dim>::
                    quadrant_is_equal (children[c],
                                       &*this_object->current_coarsen_pointer),
                    ExcInternalError());
            ++this_object->current_coarsen_pointer;
          }

        return true;
      }

    // p4est cell is not in list
    return false;
  }



  /**
   * A data structure that we use to store the weights of all cells to
   * be used upon partitioning. The class stores them in the order in
   * which p4est will encounter cells, not in the order in which
   * deal.II walks over them.
   */
  template <int dim, int spacedim>
  class PartitionWeights
  {
  public:
    /**
     * This constructor assumes the cell_weights are already sorted in the
     * order that p4est will encounter the cells, and they do not contain
     * ghost cells or artificial cells.
     */
    PartitionWeights (const std::vector<unsigned int> &cell_weights);

    /**
     * A callback function that we pass to the p4est data structures when a
     * forest is to be partitioned. The p4est functions call it back with a tree
     * (the index of the tree that grows out of a given coarse cell) and a
     * refinement path from that coarse cell to a terminal/leaf cell. The
     * function returns the weight of the cell.
     */
    static
    int
    cell_weight (typename internal::p4est::types<dim>::forest *forest,
                 typename internal::p4est::types<dim>::topidx  coarse_cell_index,
                 typename internal::p4est::types<dim>::quadrant *quadrant);

  private:
    std::vector<unsigned int> cell_weights_list;
    std::vector<unsigned int>::const_iterator current_pointer;
  };


  template <int dim, int spacedim>
  PartitionWeights<dim,spacedim>::
  PartitionWeights (const std::vector<unsigned int> &cell_weights)
    :
    cell_weights_list(cell_weights)
  {
    // set the current pointer to the first element of the list, given that
    // we will walk through it sequentially
    current_pointer  = cell_weights_list.begin();
  }


  template <int dim, int spacedim>
  int
  PartitionWeights<dim,spacedim>::
  cell_weight (typename internal::p4est::types<dim>::forest *forest,
               typename internal::p4est::types<dim>::topidx,
               typename internal::p4est::types<dim>::quadrant *)
  {
    // the function gets two additional arguments, but we don't need them
    // since we know in which order p4est will walk through the cells
    // and have already built our weight lists in this order

    PartitionWeights<dim,spacedim> *this_object
      = reinterpret_cast<PartitionWeights<dim,spacedim>*>(forest->user_pointer);

    Assert (this_object->current_pointer >= this_object->cell_weights_list.begin(),
            ExcInternalError());
    Assert (this_object->current_pointer < this_object->cell_weights_list.end(),
            ExcInternalError());

    // get the weight, increment the pointer, and return the weight
    return *this_object->current_pointer++;
  }
}


// initialize p4est
namespace internal
{
  namespace p4est
  {
    struct InitFinalize
    {
    private:
      struct Singleton
      {
        Singleton ()
        {
          // ensure that the initialization code is run only once, even if we
          // link with 1d, 2d, and 3d libraries
          static bool initialized = false;

          if (initialized == false)
            {
              sc_init (MPI_COMM_WORLD,
                       0, 0, 0, SC_LP_SILENT);
              p4est_init (0, SC_LP_SILENT);

              initialized = true;
            }
        }

        ~Singleton ()
        {
          // same here
          static bool deinitialized = false;

          if (deinitialized == false)
            {
              // p4est has no p4est_finalize function
              sc_finalize ();

              deinitialized = true;
            }
        }
      };

    public:
      // do run the initialization code, at least the first time around we get
      // to this function
      static void do_initialize ()
      {
        static Singleton singleton;
      }
    };
  }
}


namespace parallel
{
  namespace distributed
  {

    /* ---------------------- class Triangulation<dim,spacedim> ------------------------------ */


    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::
    Triangulation (MPI_Comm mpi_communicator,
                   const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid,
                   const Settings settings_)
      :
      // do not check for distorted cells
      dealii::parallel::Triangulation<dim,spacedim>
      (mpi_communicator,
       smooth_grid,
       false),
      settings(settings_),
      triangulation_has_content (false),
      connectivity (0),
      parallel_forest (0),
      refinement_in_progress (false),
      attached_data_size(0),
      n_attached_datas(0),
      n_attached_deserialize(0)
    {
      // initialize p4est. do this in a separate function since it has to
      // happen only once, even if we have triangulation objects for several
      // different space dimensions
      dealii::internal::p4est::InitFinalize::do_initialize ();

      parallel_ghost = 0;
    }



    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::~Triangulation ()
    {
      clear ();

      Assert (triangulation_has_content == false,
              ExcInternalError());
      Assert (connectivity == 0,    ExcInternalError());
      Assert (parallel_forest == 0, ExcInternalError());
      Assert (refinement_in_progress == false, ExcInternalError());
    }




    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    create_triangulation (const std::vector<Point<spacedim> >    &vertices,
                          const std::vector<CellData<dim> > &cells,
                          const SubCellData                 &subcelldata)
    {
      try
        {
          dealii::Triangulation<dim,spacedim>::
          create_triangulation (vertices, cells, subcelldata);
        }
      catch (const typename dealii::Triangulation<dim,spacedim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

      // note that now we have some content in the p4est objects and call the
      // functions that do the actual work (which are dimension dependent, so
      // separate)
      triangulation_has_content = true;

      setup_coarse_cell_to_p4est_tree_permutation ();

      copy_new_triangulation_to_p4est (dealii::internal::int2type<dim>());

      try
        {
          copy_local_forest_to_triangulation ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

      this->update_number_cache ();
      this->update_periodic_face_map();
    }


    // This anonymous namespace contains utility for
    // the function Triangulation::communicate_locally_moved_vertices
    namespace CommunicateLocallyMovedVertices
    {
      namespace
      {
        /**
          * A list of tree+quadrant and their vertex indices.
          * The bool vector describes which vertices are of interest
          * and should be set on the receiving processes.
          */
        template <int dim, int spacedim>
        struct CellInfo
        {
          // store all the tree_indices we send/receive consecutively (n_cells entries)
          std::vector<unsigned int> tree_index;
          // store all the quadrants we send/receive consecutively (n_cells entries)
          std::vector<typename dealii::internal::p4est::types<dim>::quadrant> quadrants;
          // store for each cell the number of vertices we send/receive
          // and then the vertex indices (for each cell: n_vertices+1 entries)
          std::vector<unsigned int> vertex_indices;
          // store for each cell the vertices we send/receive
          // (for each cell n_vertices entries)
          std::vector<dealii::Point<spacedim> > vertices;
          // for receiving and unpacking data we need to store pointers to the
          // first vertex and vertex_index on each cell additionally
          // both vectors have as many entries as there are cells
          std::vector<unsigned int * > first_vertex_indices;
          std::vector<dealii::Point<spacedim>* > first_vertices;

          unsigned int bytes_for_buffer () const
          {
            return (sizeof(unsigned int) +
                    tree_index.size() * sizeof(unsigned int) +
                    quadrants.size() * sizeof(typename dealii::internal::p4est
                                              ::types<dim>::quadrant) +
                    vertices.size() * sizeof(dealii::Point<spacedim>)) +
                   vertex_indices.size() * sizeof(unsigned int);
          }

          void pack_data (std::vector<char> &buffer) const
          {
            buffer.resize(bytes_for_buffer());

            char *ptr = &buffer[0];

            const unsigned int num_cells = tree_index.size();
            std::memcpy(ptr, &num_cells, sizeof(unsigned int));
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

            std::memcpy(ptr,
                        &vertex_indices[0],
                        vertex_indices.size() * sizeof(unsigned int));
            ptr += vertex_indices.size() * sizeof(unsigned int);

            std::memcpy(ptr,
                        &vertices[0],
                        vertices.size() * sizeof(dealii::Point<spacedim>));
            ptr += vertices.size() * sizeof(dealii::Point<spacedim>);

            Assert (ptr == &buffer[0]+buffer.size(),
                    ExcInternalError());

          }

          void unpack_data (const std::vector<char> &buffer)
          {
            const char *ptr = &buffer[0];
            unsigned int cells;
            memcpy(&cells, ptr, sizeof(unsigned int));
            ptr += sizeof(unsigned int);

            tree_index.resize(cells);
            memcpy(&tree_index[0],ptr,sizeof(unsigned int)*cells);
            ptr += sizeof(unsigned int)*cells;

            quadrants.resize(cells);
            memcpy(&quadrants[0],ptr,
                   sizeof(typename dealii::internal::p4est::types<dim>::quadrant)*cells);
            ptr += sizeof(typename dealii::internal::p4est::types<dim>::quadrant)*cells;

            vertex_indices.clear();
            first_vertex_indices.resize(cells);
            std::vector<unsigned int> n_vertices_on_cell(cells);
            std::vector<unsigned int> first_indices (cells);
            for (unsigned int c=0; c<cells; ++c)
              {
                // The first 'vertex index' is the number of vertices.
                // Additionally, we need to store the pointer to this
                // vertex index with respect to the std::vector
                const unsigned int *const vertex_index
                  = reinterpret_cast<const unsigned int *const>(ptr);
                first_indices[c] = vertex_indices.size();
                vertex_indices.push_back(*vertex_index);
                n_vertices_on_cell[c] = *vertex_index;
                ptr += sizeof(unsigned int);
                // Now copy all the 'real' vertex_indices
                vertex_indices.resize(vertex_indices.size() + n_vertices_on_cell[c]);
                memcpy(&vertex_indices[vertex_indices.size() - n_vertices_on_cell[c]],
                       ptr, n_vertices_on_cell[c]*sizeof(unsigned int));
                ptr += n_vertices_on_cell[c]*sizeof(unsigned int);
              }
            for (unsigned int c=0; c<cells; ++c)
              first_vertex_indices[c] = &vertex_indices[first_indices[c]];

            vertices.clear();
            first_vertices.resize(cells);
            for (unsigned int c=0; c<cells; ++c)
              {
                // We need to store a pointer to the first vertex.
                const dealii::Point<spacedim> *const vertex
                  = reinterpret_cast<const dealii::Point<spacedim> * const>(ptr);
                first_indices[c] = vertices.size();
                vertices.push_back(*vertex);
                ptr += sizeof(dealii::Point<spacedim>);
                vertices.resize(vertices.size() + n_vertices_on_cell[c]-1);
                memcpy(&vertices[vertices.size() - (n_vertices_on_cell[c]-1)],
                       ptr, (n_vertices_on_cell[c]-1)*sizeof(dealii::Point<spacedim>));
                ptr += (n_vertices_on_cell[c]-1)*sizeof(dealii::Point<spacedim>);
              }
            for (unsigned int c=0; c<cells; ++c)
              first_vertices[c] = &vertices[first_indices[c]];

            Assert (ptr == &buffer[0]+buffer.size(),
                    ExcInternalError());
          }
        };



        // This function is responsible for gathering the information
        // we want to send to each process.
        // For each dealii cell on the coarsest level the corresponding
        // p4est_cell has to be provided when calling this function.
        // By recursing through all children we consider each active cell.
        // vertices_with_ghost_neighbors tells us which vertices
        // are in the ghost layer and for which processes they might
        // be interesting.
        // Whether a vertex has actually been updated locally is
        // stored in vertex_locally_moved. Only those are considered
        // for sending.
        // The gathered information is saved into needs_to_get_cell.
        template <int dim, int spacedim>
        void
        fill_vertices_recursively (const typename parallel::distributed::Triangulation<dim,spacedim> &tria,
                                   const unsigned int tree_index,
                                   const typename Triangulation<dim,spacedim>::cell_iterator &dealii_cell,
                                   const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
                                   const std::map<unsigned int, std::set<dealii::types::subdomain_id> > &vertices_with_ghost_neighbors,
                                   const std::vector<bool> &vertex_locally_moved,
                                   std::map<dealii::types::subdomain_id, CellInfo<dim, spacedim> > &needs_to_get_cell)
        {
          // see if we have to
          // recurse...
          if (dealii_cell->has_children())
            {
              typename dealii::internal::p4est::types<dim>::quadrant
              p4est_child[GeometryInfo<dim>::max_children_per_cell];
              dealii::internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);


              for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                fill_vertices_recursively<dim,spacedim>(tria,
                                                        tree_index,
                                                        dealii_cell->child(c),
                                                        p4est_child[c],
                                                        vertices_with_ghost_neighbors,
                                                        vertex_locally_moved,
                                                        needs_to_get_cell);
              return;
            }

          // We're at a leaf cell. If the cell is locally owned, we may
          // have to send its vertices to other processors if any of
          // its vertices is adjacent to a ghost cell and has been moved
          //
          // If one of the vertices of the cell is interesting,
          // send all moved vertices of the cell to all processors
          // adjacent to all cells adjacent to this vertex
          if (dealii_cell->is_locally_owned())
            {
              std::set<dealii::types::subdomain_id> send_to;
              for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
                {
                  const std::map<unsigned int, std::set<dealii::types::subdomain_id> >::const_iterator
                  neighbor_subdomains_of_vertex
                    = vertices_with_ghost_neighbors.find (dealii_cell->vertex_index(v));

                  if (neighbor_subdomains_of_vertex
                      != vertices_with_ghost_neighbors.end())
                    {
                      Assert(neighbor_subdomains_of_vertex->second.size()!=0,
                             ExcInternalError());
                      send_to.insert(neighbor_subdomains_of_vertex->second.begin(),
                                     neighbor_subdomains_of_vertex->second.end());
                    }
                }

              if (send_to.size() > 0)
                {
                  std::vector<unsigned int> vertex_indices;
                  std::vector<dealii::Point<spacedim> > local_vertices;
                  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
                    if (vertex_locally_moved[dealii_cell->vertex_index(v)])
                      {
                        vertex_indices.push_back(v);
                        local_vertices.push_back(dealii_cell->vertex(v));
                      }

                  if (vertex_indices.size()>0)
                    for (std::set<dealii::types::subdomain_id>::iterator it=send_to.begin();
                         it!=send_to.end(); ++it)
                      {
                        const dealii::types::subdomain_id subdomain = *it;

                        // get an iterator to what needs to be sent to that
                        // subdomain (if already exists), or create such an object
                        const typename std::map<dealii::types::subdomain_id, CellInfo<dim, spacedim> >::iterator
                        p
                          = needs_to_get_cell.insert (std::make_pair(subdomain,
                                                                     CellInfo<dim,spacedim>()))
                            .first;

                        p->second.tree_index.push_back(tree_index);
                        p->second.quadrants.push_back(p4est_cell);

                        p->second.vertex_indices.push_back(vertex_indices.size());
                        p->second.vertex_indices.insert(p->second.vertex_indices.end(),
                                                        vertex_indices.begin(),
                                                        vertex_indices.end());

                        p->second.vertices.insert(p->second.vertices.end(),
                                                  local_vertices.begin(),
                                                  local_vertices.end());
                      }
                }
            }
        }



        // After the cell data has been received this function is responsible
        // for moving the vertices in the corresponding ghost layer locally.
        // As in fill_vertices_recursively for each dealii cell on the
        // coarsest level the corresponding p4est_cell has to be provided
        // when calling this function. By recursing through through all
        // children we consider each active cell.
        // Additionally, we need to give a pointer to the first vertex indices
        // and vertices. Since the first information saved in vertex_indices
        // is the number of vertices this all the information we need.
        template <int dim, int spacedim>
        void
        set_vertices_recursively (
          const parallel::distributed::Triangulation<dim,spacedim> &tria,
          const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell,
          const typename Triangulation<dim,spacedim>::cell_iterator &dealii_cell,
          const typename dealii::internal::p4est::types<dim>::quadrant &quadrant,
          const dealii::Point<spacedim> *const vertices,
          const unsigned int *const vertex_indices)
        {
          if (dealii::internal::p4est::quadrant_is_equal<dim>(p4est_cell, quadrant))
            {
              Assert(!dealii_cell->is_artificial(), ExcInternalError());
              Assert(!dealii_cell->has_children(), ExcInternalError());
              Assert(!dealii_cell->is_locally_owned(), ExcInternalError());

              const unsigned int n_vertices = vertex_indices[0];

              // update dof indices of cell
              for (unsigned int i=0; i<n_vertices; ++i)
                dealii_cell->vertex(vertex_indices[i+1]) = vertices[i];

              return;
            }

          if (! dealii_cell->has_children())
            return;

          if (! dealii::internal::p4est::quadrant_is_ancestor<dim> (p4est_cell, quadrant))
            return;

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
          dealii::internal::p4est::init_quadrant_children<dim>(p4est_cell, p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            set_vertices_recursively<dim,spacedim> (tria, p4est_child[c],
                                                    dealii_cell->child(c),
                                                    quadrant, vertices,
                                                    vertex_indices);
        }
      }
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::clear ()
    {
      triangulation_has_content = false;

      if (parallel_ghost != 0)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy (parallel_ghost);
          parallel_ghost = 0;
        }

      if (parallel_forest != 0)
        {
          dealii::internal::p4est::functions<dim>::destroy (parallel_forest);
          parallel_forest = 0;
        }

      if (connectivity != 0)
        {
          dealii::internal::p4est::functions<dim>::connectivity_destroy (connectivity);
          connectivity = 0;
        }

      coarse_cell_to_p4est_tree_permutation.resize (0);
      p4est_tree_to_coarse_cell_permutation.resize (0);

      dealii::Triangulation<dim,spacedim>::clear ();

      this->update_number_cache ();
    }

    template <int dim, int spacedim>
    bool
    Triangulation<dim,spacedim>::has_hanging_nodes () const
    {
      if (this->n_global_levels()<=1)
        return false; // can not have hanging nodes without refined cells

      // if there are any active cells with level less than n_global_levels()-1, then
      // there is obviously also one with level n_global_levels()-1, and
      // consequently there must be a hanging node somewhere.
      //
      // The problem is that we cannot just ask for the first active cell, but
      // instead need to filter over locally owned cells.
      bool have_coarser_cell = false;
      for (typename Triangulation<dim, spacedim>::active_cell_iterator cell = this->begin_active(this->n_global_levels()-2);
           cell != this->end(this->n_global_levels()-2);
           ++cell)
        if (cell->is_locally_owned())
          {
            have_coarser_cell = true;
            break;
          }

      // return true if at least one process has a coarser cell
      return 0<Utilities::MPI::max(have_coarser_cell?1:0, this->mpi_communicator);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::setup_coarse_cell_to_p4est_tree_permutation ()
    {
      DynamicSparsityPattern cell_connectivity;
      GridTools::get_vertex_connectivity_of_cells (*this, cell_connectivity);
      coarse_cell_to_p4est_tree_permutation.resize (this->n_cells(0));
      SparsityTools::
      reorder_hierarchical (cell_connectivity,
                            coarse_cell_to_p4est_tree_permutation);

      p4est_tree_to_coarse_cell_permutation
        = Utilities::invert_permutation (coarse_cell_to_p4est_tree_permutation);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::write_mesh_vtk (const char *file_basename) const
    {
      Assert (parallel_forest != 0,
              ExcMessage ("Can't produce output when no forest is created yet."));
      dealii::internal::p4est::functions<dim>::
      vtk_write_file (parallel_forest, 0, file_basename);
    }


    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    save(const char *filename) const
    {
      Assert(n_attached_deserialize==0,
             ExcMessage ("not all SolutionTransfer's got deserialized after the last load()"));
      int real_data_size = 0;
      if (attached_data_size>0)
        real_data_size = attached_data_size+sizeof(CellStatus);

      Assert(this->n_cells()>0, ExcMessage("Can not save() an empty Triangulation."));

      if (this->my_subdomain==0)
        {
          std::string fname=std::string(filename)+".info";
          std::ofstream f(fname.c_str());
          f << "version nproc attached_bytes n_attached_objs n_coarse_cells" << std::endl
            << 2 << " "
            << Utilities::MPI::n_mpi_processes (this->mpi_communicator) << " "
            << real_data_size << " "
            << attached_data_pack_callbacks.size() << " "
            << this->n_cells(0)
            << std::endl;
        }

      if (attached_data_size>0)
        {
          const_cast<dealii::parallel::distributed::Triangulation<dim, spacedim>*>(this)
          ->attach_mesh_data();
        }

      dealii::internal::p4est::functions<dim>::save(filename, parallel_forest, attached_data_size>0);

      dealii::parallel::distributed::Triangulation<dim, spacedim> *tria
        = const_cast<dealii::parallel::distributed::Triangulation<dim, spacedim>*>(this);

      tria->n_attached_datas = 0;
      tria->attached_data_size = 0;
      tria->attached_data_pack_callbacks.clear();

      // and release the data
      void *userptr = parallel_forest->user_pointer;
      dealii::internal::p4est::functions<dim>::reset_data (parallel_forest, 0, NULL, NULL);
      parallel_forest->user_pointer = userptr;
    }


    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    load (const char *filename,
          const bool autopartition)
    {
      Assert(this->n_cells()>0, ExcMessage("load() only works if the Triangulation already contains a coarse mesh!"));
      Assert(this->n_levels()==1, ExcMessage("Triangulation may only contain coarse cells when calling load()."));


      if (parallel_ghost != 0)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy (parallel_ghost);
          parallel_ghost = 0;
        }
      dealii::internal::p4est::functions<dim>::destroy (parallel_forest);
      parallel_forest = 0;
      dealii::internal::p4est::functions<dim>::connectivity_destroy (connectivity);
      connectivity = 0;

      unsigned int version, numcpus, attached_size, attached_count, n_coarse_cells;
      {
        std::string fname=std::string(filename)+".info";
        std::ifstream f(fname.c_str());
        std::string firstline;
        getline(f, firstline); //skip first line
        f >> version >> numcpus >> attached_size >> attached_count >> n_coarse_cells;
      }

      Assert(version == 2, ExcMessage("Incompatible version found in .info file."));
      Assert(this->n_cells(0) == n_coarse_cells, ExcMessage("Number of coarse cells differ!"));
#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
#else
      AssertThrow(numcpus <= Utilities::MPI::n_mpi_processes (this->mpi_communicator),
                  ExcMessage("parallel::distributed::Triangulation::load() only supports loading "
                             "saved data with a greater or equal number of processes than were used to "
                             "save() when using p4est 0.3.4.2."));
#endif

      attached_data_size = 0;
      n_attached_datas = 0;
      n_attached_deserialize = attached_count;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      parallel_forest = dealii::internal::p4est::functions<dim>::load_ext (
                          filename, this->mpi_communicator,
                          attached_size, attached_size>0,
                          autopartition, 0,
                          this,
                          &connectivity);
#else
      (void)autopartition;
      parallel_forest = dealii::internal::p4est::functions<dim>::load (
                          filename, this->mpi_communicator,
                          attached_size, attached_size>0,
                          this,
                          &connectivity);
#endif
      if (numcpus != Utilities::MPI::n_mpi_processes (this->mpi_communicator))
        // We are changing the number of CPUs so we need to repartition.
        // Note that p4est actually distributes the cells between the changed
        // number of CPUs and so everything works without this call, but
        // this command changes the distribution for some reason, so we
        // will leave it in here.
        repartition();

      try
        {
          copy_local_forest_to_triangulation ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying
          // triangulation should not
          // be checking for
          // distorted cells
          AssertThrow (false, ExcInternalError());
        }

      this->update_number_cache ();
      this->update_periodic_face_map();
    }



    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim,spacedim>::get_checksum () const
    {
      Assert (parallel_forest != 0,
              ExcMessage ("Can't produce a check sum when no forest is created yet."));
      return dealii::internal::p4est::functions<dim>::checksum (parallel_forest);
    }

    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    update_number_cache ()
    {
      parallel::Triangulation<dim,spacedim>::update_number_cache();

      if (this->n_levels() == 0)
        return;

      if (settings & construct_multigrid_hierarchy)
        {
          // find level ghost owners
          for (typename Triangulation<dim,spacedim>::cell_iterator
               cell = this->begin();
               cell != this->end();
               ++cell)
            if (cell->level_subdomain_id() != numbers::artificial_subdomain_id
                && cell->level_subdomain_id() != this->locally_owned_subdomain())
              this->number_cache.level_ghost_owners.insert(cell->level_subdomain_id());

#ifdef DEBUG
          // Check that level_ghost_owners is symmetric by sending a message
          // to everyone
          {

            MPI_Barrier(this->mpi_communicator);

            // important: preallocate to avoid (re)allocation:
            std::vector<MPI_Request> requests (this->number_cache.level_ghost_owners.size());
            int dummy = 0;
            unsigned int req_counter = 0;

            for (std::set<unsigned int>::iterator it = this->number_cache.level_ghost_owners.begin();
                 it != this->number_cache.level_ghost_owners.end();
                 ++it, ++req_counter)
              {
                MPI_Isend(&dummy, 1, MPI_INT,
                          *it, 9001, this->mpi_communicator,
                          &requests[req_counter]);
              }

            for (std::set<unsigned int>::iterator it = this->number_cache.level_ghost_owners.begin();
                 it != this->number_cache.level_ghost_owners.end();
                 ++it)
              {
                int dummy;
                MPI_Recv(&dummy, 1, MPI_INT,
                         *it, 9001, this->mpi_communicator,
                         MPI_STATUS_IGNORE);
              }

            if (requests.size() > 0)
              MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

            MPI_Barrier(this->mpi_communicator);
          }
#endif

          Assert(this->number_cache.level_ghost_owners.size() < Utilities::MPI::n_mpi_processes(this->mpi_communicator), ExcInternalError());
        }
    }


    template <int dim, int spacedim>
    typename dealii::internal::p4est::types<dim>::tree *
    Triangulation<dim,spacedim>::
    init_tree(const int dealii_coarse_cell_index) const
    {
      const unsigned int tree_index
        = coarse_cell_to_p4est_tree_permutation[dealii_coarse_cell_index];
      typename dealii::internal::p4est::types<dim>::tree *tree
        = static_cast<typename dealii::internal::p4est::types<dim>::tree *>
          (sc_array_index (parallel_forest->trees,
                           tree_index));

      return tree;
    }



    template <>
    void
    Triangulation<2,2>::copy_new_triangulation_to_p4est (dealii::internal::int2type<2>)
    {
      const unsigned int dim = 2, spacedim = 2;
      Assert (this->n_cells(0) > 0, ExcInternalError());
      Assert (this->n_levels() == 1, ExcInternalError());

      // data structures that counts how many cells touch each vertex
      // (vertex_touch_count), and which cells touch a given vertex (together
      // with the local numbering of that vertex within the cells that touch
      // it)
      std::vector<unsigned int> vertex_touch_count;
      std::vector<
      std::list<
      std::pair<Triangulation<dim,spacedim>::active_cell_iterator,
          unsigned int> > >
          vertex_to_cell;
      get_vertex_to_cell_mappings (*this,
                                   vertex_touch_count,
                                   vertex_to_cell);
      const dealii::internal::p4est::types<2>::locidx
      num_vtt = std::accumulate (vertex_touch_count.begin(),
                                 vertex_touch_count.end(),
                                 0);

      // now create a connectivity object with the right sizes for all
      // arrays. set vertex information only in debug mode (saves a few bytes
      // in optimized mode)
      const bool set_vertex_info
#ifdef DEBUG
        = true
#else
        = false
#endif
          ;

      connectivity
        = dealii::internal::p4est::functions<2>::
          connectivity_new ((set_vertex_info == true ? this->n_vertices() : 0),
                            this->n_cells(0),
                            this->n_vertices(),
                            num_vtt);

      set_vertex_and_cell_info (*this,
                                vertex_touch_count,
                                vertex_to_cell,
                                coarse_cell_to_p4est_tree_permutation,
                                set_vertex_info,
                                connectivity);

      Assert (p4est_connectivity_is_valid (connectivity) == 1,
              ExcInternalError());

      // now create a forest out of the connectivity data structure
      parallel_forest
        = dealii::internal::p4est::functions<2>::
          new_forest (this->mpi_communicator,
                      connectivity,
                      /* minimum initial number of quadrants per tree */ 0,
                      /* minimum level of upfront refinement */ 0,
                      /* use uniform upfront refinement */ 1,
                      /* user_data_size = */ 0,
                      /* user_data_constructor = */ NULL,
                      /* user_pointer */ this);
    }



    // TODO: This is a verbatim copy of the 2,2 case. However, we can't just
    // specialize the dim template argument, but let spacedim open
    template <>
    void
    Triangulation<2,3>::copy_new_triangulation_to_p4est (dealii::internal::int2type<2>)
    {
      const unsigned int dim = 2, spacedim = 3;
      Assert (this->n_cells(0) > 0, ExcInternalError());
      Assert (this->n_levels() == 1, ExcInternalError());

      // data structures that counts how many cells touch each vertex
      // (vertex_touch_count), and which cells touch a given vertex (together
      // with the local numbering of that vertex within the cells that touch
      // it)
      std::vector<unsigned int> vertex_touch_count;
      std::vector<
      std::list<
      std::pair<Triangulation<dim,spacedim>::active_cell_iterator,
          unsigned int> > >
          vertex_to_cell;
      get_vertex_to_cell_mappings (*this,
                                   vertex_touch_count,
                                   vertex_to_cell);
      const dealii::internal::p4est::types<2>::locidx
      num_vtt = std::accumulate (vertex_touch_count.begin(),
                                 vertex_touch_count.end(),
                                 0);

      // now create a connectivity object with the right sizes for all
      // arrays. set vertex information only in debug mode (saves a few bytes
      // in optimized mode)
      const bool set_vertex_info
#ifdef DEBUG
        = true
#else
        = false
#endif
          ;

      connectivity
        = dealii::internal::p4est::functions<2>::
          connectivity_new ((set_vertex_info == true ? this->n_vertices() : 0),
                            this->n_cells(0),
                            this->n_vertices(),
                            num_vtt);

      set_vertex_and_cell_info (*this,
                                vertex_touch_count,
                                vertex_to_cell,
                                coarse_cell_to_p4est_tree_permutation,
                                set_vertex_info,
                                connectivity);

      Assert (p4est_connectivity_is_valid (connectivity) == 1,
              ExcInternalError());

      // now create a forest out of the connectivity data structure
      parallel_forest
        = dealii::internal::p4est::functions<2>::
          new_forest (this->mpi_communicator,
                      connectivity,
                      /* minimum initial number of quadrants per tree */ 0,
                      /* minimum level of upfront refinement */ 0,
                      /* use uniform upfront refinement */ 1,
                      /* user_data_size = */ 0,
                      /* user_data_constructor = */ NULL,
                      /* user_pointer */ this);
    }



    template <>
    void
    Triangulation<3,3>::copy_new_triangulation_to_p4est (dealii::internal::int2type<3>)
    {
      const int dim = 3, spacedim = 3;
      Assert (this->n_cells(0) > 0, ExcInternalError());
      Assert (this->n_levels() == 1, ExcInternalError());

      // data structures that counts how many cells touch each vertex
      // (vertex_touch_count), and which cells touch a given vertex (together
      // with the local numbering of that vertex within the cells that touch
      // it)
      std::vector<unsigned int> vertex_touch_count;
      std::vector<
      std::list<
      std::pair<Triangulation<3>::active_cell_iterator,
          unsigned int> > >
          vertex_to_cell;
      get_vertex_to_cell_mappings (*this,
                                   vertex_touch_count,
                                   vertex_to_cell);
      const dealii::internal::p4est::types<2>::locidx
      num_vtt = std::accumulate (vertex_touch_count.begin(),
                                 vertex_touch_count.end(),
                                 0);

      std::vector<unsigned int> edge_touch_count;
      std::vector<
      std::list<
      std::pair<Triangulation<3>::active_cell_iterator,
          unsigned int> > >
          edge_to_cell;
      get_edge_to_cell_mappings (*this,
                                 edge_touch_count,
                                 edge_to_cell);
      const dealii::internal::p4est::types<2>::locidx
      num_ett = std::accumulate (edge_touch_count.begin(),
                                 edge_touch_count.end(),
                                 0);

      // now create a connectivity object with the right sizes for all arrays
      const bool set_vertex_info
#ifdef DEBUG
        = true
#else
        = false
#endif
          ;

      connectivity
        = dealii::internal::p4est::functions<3>::
          connectivity_new ((set_vertex_info == true ? this->n_vertices() : 0),
                            this->n_cells(0),
                            this->n_active_lines(),
                            num_ett,
                            this->n_vertices(),
                            num_vtt);

      set_vertex_and_cell_info (*this,
                                vertex_touch_count,
                                vertex_to_cell,
                                coarse_cell_to_p4est_tree_permutation,
                                set_vertex_info,
                                connectivity);

      // next to tree-to-edge
      // data. note that in p4est lines
      // are ordered as follows
      //      *---3---*        *---3---*
      //     /|       |       /       /|
      //    6 |       11     6       7 11
      //   /  10      |     /       /  |
      //  *   |       |    *---2---*   |
      //  |   *---1---*    |       |   *
      //  |  /       /     |       9  /
      //  8 4       5      8       | 5
      //  |/       /       |       |/
      //  *---0---*        *---0---*
      // whereas in deal.II they are like this:
      //      *---7---*        *---7---*
      //     /|       |       /       /|
      //    4 |       11     4       5 11
      //   /  10      |     /       /  |
      //  *   |       |    *---6---*   |
      //  |   *---3---*    |       |   *
      //  |  /       /     |       9  /
      //  8 0       1      8       | 1
      //  |/       /       |       |/
      //  *---2---*        *---2---*

      const unsigned int deal_to_p4est_line_index[12]
        = { 4, 5, 0, 1,  6, 7, 2, 3, 8, 9, 10, 11 } ;

      for (Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end(); ++cell)
        {
          const unsigned int
          index = coarse_cell_to_p4est_tree_permutation[cell->index()];
          for (unsigned int e=0; e<GeometryInfo<3>::lines_per_cell; ++e)
            connectivity->tree_to_edge[index*GeometryInfo<3>::lines_per_cell+
                                       deal_to_p4est_line_index[e]]
              = cell->line(e)->index();
        }

      // now also set edge-to-tree
      // information
      connectivity->ett_offset[0] = 0;
      std::partial_sum (edge_touch_count.begin(),
                        edge_touch_count.end(),
                        &connectivity->ett_offset[1]);

      Assert (connectivity->ett_offset[this->n_active_lines()] ==
              num_ett,
              ExcInternalError());

      for (unsigned int v=0; v<this->n_active_lines(); ++v)
        {
          Assert (edge_to_cell[v].size() == edge_touch_count[v],
                  ExcInternalError());

          std::list<std::pair
          <Triangulation<dim,spacedim>::active_cell_iterator,
          unsigned int> >::const_iterator
          p = edge_to_cell[v].begin();
          for (unsigned int c=0; c<edge_touch_count[v]; ++c, ++p)
            {
              connectivity->edge_to_tree[connectivity->ett_offset[v]+c]
                = coarse_cell_to_p4est_tree_permutation[p->first->index()];
              connectivity->edge_to_edge[connectivity->ett_offset[v]+c]
                = deal_to_p4est_line_index[p->second];
            }
        }

      Assert (p8est_connectivity_is_valid (connectivity) == 1,
              ExcInternalError());

      // now create a forest out of the connectivity data structure
      parallel_forest
        = dealii::internal::p4est::functions<3>::
          new_forest (this->mpi_communicator,
                      connectivity,
                      /* minimum initial number of quadrants per tree */ 0,
                      /* minimum level of upfront refinement */ 0,
                      /* use uniform upfront refinement */ 1,
                      /* user_data_size = */ 0,
                      /* user_data_constructor = */ NULL,
                      /* user_pointer */ this);
    }



    namespace
    {
      // ensures the 2:1 mesh balance for periodic boundary conditions in the
      // artificial cell layer (the active cells are taken care of by p4est)
      template <int dim, int spacedim>
      bool enforce_mesh_balance_over_periodic_boundaries
      (Triangulation<dim,spacedim> &tria)
      {
        if (tria.get_periodic_face_map().size()==0)
          return false;

        std::vector<bool> flags_before[2];
        tria.save_coarsen_flags (flags_before[0]);
        tria.save_refine_flags (flags_before[1]);

        std::vector<unsigned int> topological_vertex_numbering(tria.n_vertices());
        for (unsigned int i=0; i<topological_vertex_numbering.size(); ++i)
          topological_vertex_numbering[i] = i;
        // combine vertices that have different locations (and thus, different
        // vertex_index) but represent the same topological entity over periodic
        // boundaries. The vector topological_vertex_numbering contains a linear
        // map from 0 to n_vertices at input and at output relates periodic
        // vertices with only one vertex index. The output is used to always
        // identify the same vertex according to the periodicity, e.g. when
        // finding the maximum cell level around a vertex.
        //
        // Example: On a 3D cell with vertices numbered from 0 to 7 and periodic
        // boundary conditions in x direction, the vector
        // topological_vertex_numbering will contain the numbers
        // {0,0,2,2,4,4,6,6} (because the vertex pairs {0,1}, {2,3}, {4,5},
        // {6,7} belong together, respectively). If periodicity is set in x and
        // z direction, the output is {0,0,2,2,0,0,2,2}, and if periodicity is
        // in all directions, the output is simply {0,0,0,0,0,0,0,0}.
        typedef typename Triangulation<dim, spacedim>::cell_iterator cell_iterator;
        typename std::map<std::pair<cell_iterator, unsigned int>,
                 std::pair<std::pair<cell_iterator,unsigned int>, std::bitset<3> > >::const_iterator it;
        for (it = tria.get_periodic_face_map().begin(); it!= tria.get_periodic_face_map().end(); ++it)
          {
            const cell_iterator &cell_1 = it->first.first;
            const unsigned int face_no_1 = it->first.second;
            const cell_iterator &cell_2 = it->second.first.first;
            const unsigned int face_no_2 = it->second.first.second;
            const std::bitset<3> face_orientation = it->second.second;

            if (cell_1->level() == cell_2->level())
              {
                for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
                  {
                    // take possible non-standard orientation of face on cell[0] into
                    // account
                    const unsigned int vface0 =
                      GeometryInfo<dim>::standard_to_real_face_vertex(v,face_orientation[0],
                                                                      face_orientation[1],
                                                                      face_orientation[2]);
                    const unsigned int vi0 = topological_vertex_numbering[cell_1->face(face_no_1)->vertex_index(vface0)];
                    const unsigned int vi1 = topological_vertex_numbering[cell_2->face(face_no_2)->vertex_index(v)];
                    const unsigned int min_index = std::min(vi0, vi1);
                    topological_vertex_numbering[cell_1->face(face_no_1)->vertex_index(vface0)]
                      = topological_vertex_numbering[cell_2->face(face_no_2)->vertex_index(v)]
                        = min_index;
                  }
              }
          }

        // There must not be any chains!
        for (unsigned int i=0; i<topological_vertex_numbering.size(); ++i)
          {
            const unsigned int j = topological_vertex_numbering[i];
            if (j != i)
              Assert(topological_vertex_numbering[j] == j,
                     ExcInternalError());
          }


        // this code is replicated from grid/tria.cc but using an indirection
        // for periodic boundary conditions
        bool continue_iterating = true;
        std::vector<int> vertex_level(tria.n_vertices());
        while (continue_iterating)
          {
            // store highest level one of the cells adjacent to a vertex
            // belongs to
            std::fill (vertex_level.begin(), vertex_level.end(), 0);
            typename Triangulation<dim,spacedim>::active_cell_iterator
            cell = tria.begin_active(), endc = tria.end();
            for (; cell!=endc; ++cell)
              {
                if (cell->refine_flag_set())
                  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                       ++vertex)
                    vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]]
                      = std::max (vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]],
                                  cell->level()+1);
                else if (!cell->coarsen_flag_set())
                  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                       ++vertex)
                    vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]]
                      = std::max (vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]],
                                  cell->level());
                else
                  {
                    // if coarsen flag is set then tentatively assume
                    // that the cell will be coarsened. this isn't
                    // always true (the coarsen flag could be removed
                    // again) and so we may make an error here. we try
                    // to correct this by iterating over the entire
                    // process until we are converged
                    Assert (cell->coarsen_flag_set(), ExcInternalError());
                    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                         ++vertex)
                      vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]]
                        = std::max (vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]],
                                    cell->level()-1);
                  }
              }

            continue_iterating = false;

            // loop over all cells in reverse order. do so because we
            // can then update the vertex levels on the adjacent
            // vertices and maybe already flag additional cells in this
            // loop
            //
            // note that not only may we have to add additional
            // refinement flags, but we will also have to remove
            // coarsening flags on cells adjacent to vertices that will
            // see refinement
            for (cell=tria.last_active(); cell != endc; --cell)
              if (cell->refine_flag_set() == false)
                {
                  for (unsigned int vertex=0;
                       vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
                    if (vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]] >=
                        cell->level()+1)
                      {
                        // remove coarsen flag...
                        cell->clear_coarsen_flag();

                        // ...and if necessary also refine the current
                        // cell, at the same time updating the level
                        // information about vertices
                        if (vertex_level[topological_vertex_numbering[cell->vertex_index(vertex)]] >
                            cell->level()+1)
                          {
                            cell->set_refine_flag();
                            continue_iterating = true;

                            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell;
                                 ++v)
                              vertex_level[topological_vertex_numbering[cell->vertex_index(v)]]
                                = std::max (vertex_level[topological_vertex_numbering[cell->vertex_index(v)]],
                                            cell->level()+1);
                          }

                        // continue and see whether we may, for example,
                        // go into the inner 'if' above based on a
                        // different vertex
                      }
                }

            // clear coarsen flag if not all children were marked
            for (typename Triangulation<dim,spacedim>::cell_iterator cell = tria.begin();
                 cell!=tria.end(); ++cell)
              {
                // nothing to do if we are already on the finest level
                if (cell->active())
                  continue;

                const unsigned int n_children=cell->n_children();
                unsigned int flagged_children=0;
                for (unsigned int child=0; child<n_children; ++child)
                  if (cell->child(child)->active() &&
                      cell->child(child)->coarsen_flag_set())
                    ++flagged_children;

                // if not all children were flagged for coarsening, remove
                // coarsen flags
                if (flagged_children < n_children)
                  for (unsigned int child=0; child<n_children; ++child)
                    if (cell->child(child)->active())
                      cell->child(child)->clear_coarsen_flag();
              }
          }
        std::vector<bool> flags_after[2];
        tria.save_coarsen_flags (flags_after[0]);
        tria.save_refine_flags (flags_after[1]);
        return ((flags_before[0] != flags_after[0]) ||
                (flags_before[1] != flags_after[1]));
      }
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim,spacedim>::prepare_coarsening_and_refinement()
    {
      std::vector<bool> flags_before[2];
      this->save_coarsen_flags (flags_before[0]);
      this->save_refine_flags (flags_before[1]);

      bool mesh_changed = false;
      do
        {
          this->dealii::Triangulation<dim,spacedim>::prepare_coarsening_and_refinement();
          this->update_periodic_face_map();
          // enforce 2:1 mesh balance over periodic boundaries
          if (this->smooth_grid &
              dealii::Triangulation<dim,spacedim>::limit_level_difference_at_vertices)
            mesh_changed = enforce_mesh_balance_over_periodic_boundaries(*this);
        }
      while (mesh_changed);

      // check if any of the refinement flags were changed during this
      // function and return that value
      std::vector<bool> flags_after[2];
      this->save_coarsen_flags (flags_after[0]);
      this->save_refine_flags (flags_after[1]);
      return ((flags_before[0] != flags_after[0]) ||
              (flags_before[1] != flags_after[1]));
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::copy_local_forest_to_triangulation ()
    {
      // disable mesh smoothing for recreating the deal.II triangulation,
      // otherwise we might not be able to reproduce the p4est mesh
      // exactly. We restore the original smoothing at the end of this
      // function. Note that the smoothing flag is used in the normal
      // refinement process.
      typename Triangulation<dim,spacedim>::MeshSmoothing
      save_smooth = this->smooth_grid;

      // We will refine manually to match the p4est further down, which
      // obeys a level difference of 2 at each vertex (see the balance call
      // to p4est). We can disable this here so we store fewer artificial
      // cells (in some cases). For geometric multigrid it turns out that
      // we will miss level cells at shared vertices if we ignore this.
      // See tests/mpi/mg_06.
      if (settings & construct_multigrid_hierarchy)
        this->smooth_grid = dealii::Triangulation<dim,spacedim>::limit_level_difference_at_vertices;
      else
        this->smooth_grid = dealii::Triangulation<dim,spacedim>::none;

      bool mesh_changed = false;

      // remove all deal.II refinements. Note that we could skip this and
      // start from our current state, because the algorithm later coarsens as
      // necessary. This has the advantage of being faster when large parts
      // of the local partition changes (likely) and gives a deterministic
      // ordering of the cells (useful for snapshot/resume).
      // TODO: is there a more efficient way to do this?
      if (settings & mesh_reconstruction_after_repartitioning)
        while (this->begin_active()->level() > 0)
          {
            for (typename Triangulation<dim, spacedim>::active_cell_iterator
                 cell = this->begin_active();
                 cell != this->end();
                 ++cell)
              {
                cell->set_coarsen_flag();
              }

            this->prepare_coarsening_and_refinement();
            const bool saved_refinement_in_progress = refinement_in_progress;
            refinement_in_progress = true;

            try
              {
                this->execute_coarsening_and_refinement();
              }
            catch (const typename Triangulation<dim, spacedim>::DistortedCellList &)
              {
                // the underlying triangulation should not be checking for
                // distorted cells
                AssertThrow (false, ExcInternalError());
              }

            refinement_in_progress = saved_refinement_in_progress;
          }


      // query p4est for the ghost cells
      if (parallel_ghost != 0)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy (parallel_ghost);
          parallel_ghost = 0;
        }
      parallel_ghost = dealii::internal::p4est::functions<dim>::ghost_new (parallel_forest,
                       (dim == 2
                        ?
                        typename dealii::internal::p4est::types<dim>::
                        balance_type(P4EST_CONNECT_CORNER)
                        :
                        typename dealii::internal::p4est::types<dim>::
                        balance_type(P8EST_CONNECT_CORNER)));

      Assert (parallel_ghost, ExcInternalError());


      // set all cells to artificial. we will later set it to the correct
      // subdomain in match_tree_recursively
      for (typename Triangulation<dim,spacedim>::cell_iterator
           cell = this->begin(0);
           cell != this->end(0);
           ++cell)
        cell->recursively_set_subdomain_id(numbers::artificial_subdomain_id);

      do
        {
          for (typename Triangulation<dim,spacedim>::cell_iterator
               cell = this->begin(0);
               cell != this->end(0);
               ++cell)
            {
              // if this processor stores no part of the forest that comes out
              // of this coarse grid cell, then we need to delete all children
              // of this cell (the coarse grid cell remains)
              if (tree_exists_locally<dim,spacedim>(parallel_forest,
                                                    coarse_cell_to_p4est_tree_permutation[cell->index()])
                  == false)
                {
                  delete_all_children<dim,spacedim> (cell);
                  if (!cell->has_children())
                    cell->set_subdomain_id (numbers::artificial_subdomain_id);
                }

              else
                {
                  // this processor stores at least a part of the tree that
                  // comes out of this cell.

                  typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
                  typename dealii::internal::p4est::types<dim>::tree *tree =
                    init_tree(cell->index());

                  dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

                  match_tree_recursively<dim,spacedim> (*tree, cell,
                                                        p4est_coarse_cell,
                                                        *parallel_forest,
                                                        this->my_subdomain);
                }
            }

          // check mesh for ghostcells, refine as necessary.  iterate over
          // every ghostquadrant, find corresponding deal coarsecell and
          // recurse.
          typename dealii::internal::p4est::types<dim>::quadrant *quadr;
          unsigned int ghost_owner=0;
          typename dealii::internal::p4est::types<dim>::topidx ghost_tree=0;

          for (unsigned int g_idx=0; g_idx<parallel_ghost->ghosts.elem_count; ++g_idx)
            {
              while (g_idx >= (unsigned int)parallel_ghost->proc_offsets[ghost_owner+1])
                ++ghost_owner;
              while (g_idx >= (unsigned int)parallel_ghost->tree_offsets[ghost_tree+1])
                ++ghost_tree;

              quadr = static_cast<typename dealii::internal::p4est::types<dim>::quadrant *>
                      ( sc_array_index(&parallel_ghost->ghosts, g_idx) );

              unsigned int coarse_cell_index =
                p4est_tree_to_coarse_cell_permutation[ghost_tree];

              match_quadrant<dim,spacedim> (this, coarse_cell_index, *quadr, ghost_owner);
            }

          // fix all the flags to make sure we have a consistent mesh
          this->prepare_coarsening_and_refinement ();

          // see if any flags are still set
          mesh_changed = false;
          for (typename Triangulation<dim,spacedim>::active_cell_iterator
               cell = this->begin_active();
               cell != this->end();
               ++cell)
            if (cell->refine_flag_set() || cell->coarsen_flag_set())
              {
                mesh_changed = true;
                break;
              }

          // actually do the refinement but prevent the refinement hook below
          // from taking over
          const bool saved_refinement_in_progress = refinement_in_progress;
          refinement_in_progress = true;

          try
            {
              this->execute_coarsening_and_refinement();
            }
          catch (const typename Triangulation<dim,spacedim>::DistortedCellList &)
            {
              // the underlying triangulation should not be checking for
              // distorted cells
              AssertThrow (false, ExcInternalError());
            }

          refinement_in_progress = saved_refinement_in_progress;
        }
      while (mesh_changed);

#ifdef DEBUG
      // check if correct number of ghosts is created
      unsigned int num_ghosts = 0;

      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end();
           ++cell)
        {
          if (cell->subdomain_id() != this->my_subdomain
              &&
              cell->subdomain_id() != numbers::artificial_subdomain_id)
            ++num_ghosts;
        }

      Assert( num_ghosts == parallel_ghost->ghosts.elem_count, ExcInternalError());
#endif



      // fill level_subdomain_ids for geometric multigrid
      // the level ownership of a cell is defined as the owner if the cell is active or as the owner of child(0)
      // we need this information for all our ancestors and the same-level neighbors of our own cells (=level ghosts)
      if (settings & construct_multigrid_hierarchy)
        {
          // step 1: We set our own ids all the way down and all the others to
          // -1. Note that we do not fill other cells we could figure out the
          // same way, because we might accidentally set an id for a cell that
          // is not a ghost cell.
          for (unsigned int lvl=this->n_levels(); lvl>0; )
            {
              --lvl;
              typename Triangulation<dim,spacedim>::cell_iterator cell, endc = this->end(lvl);
              for (cell = this->begin(lvl); cell!=endc; ++cell)
                {
                  if ((!cell->has_children() && cell->subdomain_id()==this->locally_owned_subdomain())
                      || (cell->has_children() && cell->child(0)->level_subdomain_id()==this->locally_owned_subdomain()))
                    cell->set_level_subdomain_id(this->locally_owned_subdomain());
                  else
                    {
                      //not our cell
                      cell->set_level_subdomain_id(numbers::artificial_subdomain_id);
                    }
                }
            }

          //step 2: make sure all the neighbors to our level_cells exist. Need
          //to look up in p4est...
          std::vector<std::vector<bool> > marked_vertices(this->n_levels());
          for (unsigned int lvl=0; lvl < this->n_levels(); ++lvl)
            marked_vertices[lvl] = mark_locally_active_vertices_on_level(lvl);

          for (typename Triangulation<dim,spacedim>::cell_iterator cell = this->begin(0); cell!=this->end(0); ++cell)
            {
              typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
              const unsigned int tree_index
                = coarse_cell_to_p4est_tree_permutation[cell->index()];
              typename dealii::internal::p4est::types<dim>::tree *tree =
                init_tree(cell->index());

              dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

              determine_level_subdomain_id_recursively<dim,spacedim> (*tree, tree_index, cell,
                                                                      p4est_coarse_cell,
                                                                      *parallel_forest,
                                                                      this->my_subdomain,
                                                                      marked_vertices);
            }

          //step 3: make sure we have the parent of our level cells
          for (unsigned int lvl=this->n_levels(); lvl>0;)
            {
              --lvl;
              typename Triangulation<dim,spacedim>::cell_iterator cell, endc = this->end(lvl);
              for (cell = this->begin(lvl); cell!=endc; ++cell)
                {
                  if (cell->has_children())
                    for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                      {
                        if (cell->child(c)->level_subdomain_id()==this->locally_owned_subdomain())
                          {
                            //at least one of the children belongs to us, so
                            //make sure we set the level subdomain id
                            types::subdomain_id mark = numbers::artificial_subdomain_id;
                            mark = cell->child(0)->level_subdomain_id();
                            Assert(mark != numbers::artificial_subdomain_id, ExcInternalError()); //we should know the child(0)
                            cell->set_level_subdomain_id(mark);
                            break;
                          }
                      }
                }
            }

        }



      // check that our local copy has exactly as many cells as the p4est
      // original (at least if we are on only one processor); for parallel
      // computations, we want to check that we have at least as many as p4est
      // stores locally (in the future we should check that we have exactly as
      // many non-artificial cells as parallel_forest->local_num_quadrants)
      {
        const unsigned int total_local_cells = this->n_active_cells();
        (void)total_local_cells;

        if (Utilities::MPI::n_mpi_processes (this->mpi_communicator) == 1)
          Assert (static_cast<unsigned int>(parallel_forest->local_num_quadrants) ==
                  total_local_cells,
                  ExcInternalError())
          else
            Assert (static_cast<unsigned int>(parallel_forest->local_num_quadrants) <=
                    total_local_cells,
                    ExcInternalError());

        // count the number of owned, active cells and compare with p4est.
        unsigned int n_owned = 0;
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
             cell = this->begin_active();
             cell != this->end(); ++cell)
          {
            if (cell->subdomain_id() == this->my_subdomain)
              ++n_owned;
          }

        Assert(static_cast<unsigned int>(parallel_forest->local_num_quadrants) ==
               n_owned, ExcInternalError());

      }

      this->smooth_grid = save_smooth;
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_coarsening_and_refinement ()
    {
      // first make sure that recursive calls are handled correctly
      if (refinement_in_progress == true)
        {
          dealii::Triangulation<dim,spacedim>::execute_coarsening_and_refinement ();
          return;
        }

      // do not allow anisotropic refinement
#ifdef DEBUG
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end(); ++cell)
        if (cell->is_locally_owned() && cell->refine_flag_set())
          Assert (cell->refine_flag_set() ==
                  RefinementPossibilities<dim>::isotropic_refinement,
                  ExcMessage ("This class does not support anisotropic refinement"));
#endif


      // safety check: p4est has an upper limit on the level of a cell
      if (this->n_levels()==dealii::internal::p4est::functions<dim>::max_level)
        {
          for (typename Triangulation<dim,spacedim>::active_cell_iterator
               cell = this->begin_active(dealii::internal::p4est::functions<dim>::max_level-1);
               cell != this->end(dealii::internal::p4est::functions<dim>::max_level-1); ++cell)
            {
              AssertThrow(!(cell->refine_flag_set()),
                          ExcMessage("Fatal Error: maximum refinement level of p4est reached."));
            }
        }

      // now do the work we're supposed to do when we are in charge
      refinement_in_progress = true;
      this->prepare_coarsening_and_refinement ();

      // make sure all flags are cleared on cells we don't own, since nothing
      // good can come of that if they are still around
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end(); ++cell)
        if (cell->is_ghost() || cell->is_artificial())
          {
            cell->clear_refine_flag ();
            cell->clear_coarsen_flag ();
          }


      // count how many cells will be refined and coarsened, and allocate that
      // much memory
      RefineAndCoarsenList<dim,spacedim>
      refine_and_coarsen_list (*this,
                               p4est_tree_to_coarse_cell_permutation,
                               this->my_subdomain);

      // copy refine and coarsen flags into p4est and execute the refinement
      // and coarsening. this uses the refine_and_coarsen_list just built,
      // which is communicated to the callback functions through
      // p4est's user_pointer object
      Assert (parallel_forest->user_pointer == this,
              ExcInternalError());
      parallel_forest->user_pointer = &refine_and_coarsen_list;

      if (parallel_ghost != 0)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy (parallel_ghost);
          parallel_ghost = 0;
        }
      dealii::internal::p4est::functions<dim>::
      refine (parallel_forest, /* refine_recursive */ false,
              &RefineAndCoarsenList<dim,spacedim>::refine_callback,
              /*init_callback=*/NULL);
      dealii::internal::p4est::functions<dim>::
      coarsen (parallel_forest, /* coarsen_recursive */ false,
               &RefineAndCoarsenList<dim,spacedim>::coarsen_callback,
               /*init_callback=*/NULL);

      // make sure all cells in the lists have been consumed
      Assert (refine_and_coarsen_list.pointers_are_at_end(),
              ExcInternalError());

      // reset the pointer
      parallel_forest->user_pointer = this;

      // enforce 2:1 hanging node condition
      dealii::internal::p4est::functions<dim>::
      balance (parallel_forest,
               /* face and corner balance */
               (dim == 2
                ?
                typename dealii::internal::p4est::types<dim>::
                balance_type(P4EST_CONNECT_FULL)
                :
                typename dealii::internal::p4est::types<dim>::
                balance_type(P8EST_CONNECT_FULL)),
               /*init_callback=*/NULL);

      // before repartitioning the mesh let others attach mesh related info
      // (such as SolutionTransfer data) to the p4est
      attach_mesh_data();

      if (!(settings & no_automatic_repartitioning))
        {
          // partition the new mesh between all processors. If cell weights have
          // not been given balance the number of cells.
          if (this->signals.cell_weight.num_slots() == 0)
            dealii::internal::p4est::functions<dim>::
            partition (parallel_forest,
                       /* prepare coarsening */ 1,
                       /* weight_callback */ NULL);
          else
            {
              // get cell weights for a weighted repartitioning.
              const std::vector<unsigned int> cell_weights = get_cell_weights();

              PartitionWeights<dim,spacedim> partition_weights (cell_weights);

              // attach (temporarily) a pointer to the cell weights through p4est's
              // user_pointer object
              Assert (parallel_forest->user_pointer == this,
                      ExcInternalError());
              parallel_forest->user_pointer = &partition_weights;

              dealii::internal::p4est::functions<dim>::
              partition (parallel_forest,
                         /* prepare coarsening */ 1,
                         /* weight_callback */ &PartitionWeights<dim,spacedim>::cell_weight);

              // reset the user pointer to its previous state
              parallel_forest->user_pointer = this;
            }
        }

      // finally copy back from local part of tree to deal.II
      // triangulation. before doing so, make sure there are no refine or
      // coarsen flags pending
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end(); ++cell)
        {
          cell->clear_refine_flag();
          cell->clear_coarsen_flag();
        }

      try
        {
          copy_local_forest_to_triangulation ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

#ifdef DEBUG
      // Check that we know the level subdomain ids of all our neighbors. This
      // also involves coarser cells that share a vertex if they are active.
      //
      // Example (M= my, O=other):
      //         *------*
      //         |      |
      //         |  O   |
      //         |      |
      // *---*---*------*
      // | M | M |
      // *---*---*
      // |   | M |
      // *---*---*
      //  ^- the parent can be owned by somebody else, so O is not a neighbor
      // one level coarser
      if (settings & construct_multigrid_hierarchy)
        {
          for (unsigned int lvl=0; lvl<this->n_global_levels(); ++lvl)
            {
              std::vector<bool> active_verts = this->mark_locally_active_vertices_on_level(lvl);

              const unsigned int maybe_coarser_lvl = (lvl>0) ? (lvl-1) : lvl;
              typename Triangulation<dim, spacedim>::cell_iterator cell = this->begin(maybe_coarser_lvl),
                                                                   endc = this->end(lvl);
              for (; cell != endc; ++cell)
                if (cell->level() == static_cast<int>(lvl) || cell->active())
                  {
                    const bool is_level_artificial =
                      (cell->level_subdomain_id() == numbers::artificial_subdomain_id);
                    bool need_to_know = false;
                    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                         ++vertex)
                      if (active_verts[cell->vertex_index(vertex)])
                        {
                          need_to_know = true;
                          break;
                        }

                    Assert(!need_to_know || !is_level_artificial,
                           ExcMessage("Internal error: the owner of cell"
                                      + cell->id().to_string()
                                      + " is unknown even though it is needed for geometric multigrid."));
                  }
            }
        }
#endif

      refinement_in_progress = false;
      this->update_number_cache ();
      this->update_periodic_face_map();
    }

    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::repartition ()
    {

#ifdef DEBUG
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end(); ++cell)
        if (cell->is_locally_owned())
          Assert (
            !cell->refine_flag_set() && !cell->coarsen_flag_set(),
            ExcMessage ("Error: There shouldn't be any cells flagged for coarsening/refinement when calling repartition()."));
#endif

      refinement_in_progress = true;

      // before repartitioning the mesh let others attach mesh related info
      // (such as SolutionTransfer data) to the p4est
      attach_mesh_data();

      if (this->signals.cell_weight.num_slots() == 0)
        {
          // no cell weights given -- call p4est's 'partition' without a
          // callback for cell weights
          dealii::internal::p4est::functions<dim>::
          partition (parallel_forest,
                     /* prepare coarsening */ 1,
                     /* weight_callback */ NULL);
        }
      else
        {
          // get cell weights for a weighted repartitioning.
          const std::vector<unsigned int> cell_weights = get_cell_weights();

          PartitionWeights<dim,spacedim> partition_weights (cell_weights);

          // attach (temporarily) a pointer to the cell weights through p4est's
          // user_pointer object
          Assert (parallel_forest->user_pointer == this,
                  ExcInternalError());
          parallel_forest->user_pointer = &partition_weights;

          dealii::internal::p4est::functions<dim>::
          partition (parallel_forest,
                     /* prepare coarsening */ 1,
                     /* weight_callback */ &PartitionWeights<dim,spacedim>::cell_weight);

          // reset the user pointer to its previous state
          parallel_forest->user_pointer = this;
        }

      try
        {
          copy_local_forest_to_triangulation ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

      refinement_in_progress = false;

      // update how many cells, edges, etc, we store locally
      this->update_number_cache ();
      this->update_periodic_face_map();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    communicate_locally_moved_vertices (const std::vector<bool> &vertex_locally_moved)
    {
      Assert (vertex_locally_moved.size() == this->n_vertices(),
              ExcDimensionMismatch(vertex_locally_moved.size(),
                                   this->n_vertices()));
#ifdef DEBUG
      {
        const std::vector<bool> locally_owned_vertices
          = GridTools::get_locally_owned_vertices (*this);
        for (unsigned int i=0; i<locally_owned_vertices.size(); ++i)
          Assert ((vertex_locally_moved[i] == false)
                  ||
                  (locally_owned_vertices[i] == true),
                  ExcMessage ("The vertex_locally_moved argument must not "
                              "contain vertices that are not locally owned"));
      }
#endif

      std::map<unsigned int, std::set<dealii::types::subdomain_id> >
      vertices_with_ghost_neighbors;

      // First find out which process should receive which vertices.
      // these are specifically the ones that sit on ghost cells and,
      // among these, the ones that we own locally
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell=this->begin_active(); cell!=this->end();
           ++cell)
        if (cell->is_ghost())
          for (unsigned int vertex_no=0;
               vertex_no<GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
            {
              const unsigned int process_local_vertex_no = cell->vertex_index(vertex_no);
              vertices_with_ghost_neighbors[process_local_vertex_no].insert
              (cell->subdomain_id());
            }

      // now collect cells and their vertices
      // for the interested neighbors
      typedef
      std::map<dealii::types::subdomain_id, CommunicateLocallyMovedVertices::CellInfo<dim,spacedim> > cellmap_t;
      cellmap_t needs_to_get_cells;

      for (typename Triangulation<dim,spacedim>::cell_iterator
           cell = this->begin(0);
           cell != this->end(0);
           ++cell)
        {
          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

          CommunicateLocallyMovedVertices::fill_vertices_recursively<dim,spacedim>
          (*this,
           this->get_coarse_cell_to_p4est_tree_permutation()[cell->index()],
           cell,
           p4est_coarse_cell,
           vertices_with_ghost_neighbors,
           vertex_locally_moved,
           needs_to_get_cells);
        }

      // sending
      std::vector<std::vector<char> > sendbuffers (needs_to_get_cells.size());
      std::vector<std::vector<char> >::iterator buffer = sendbuffers.begin();
      std::vector<MPI_Request> requests (needs_to_get_cells.size());
      std::vector<unsigned int> destinations;

      unsigned int idx=0;

      for (typename cellmap_t::iterator it=needs_to_get_cells.begin();
           it!=needs_to_get_cells.end();
           ++it, ++buffer, ++idx)
        {
          const unsigned int num_cells = it->second.tree_index.size();
          (void)num_cells;
          destinations.push_back(it->first);

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
                    123, this->get_communicator(), &requests[idx]);
        }

      Assert(destinations.size()==needs_to_get_cells.size(), ExcInternalError());

      // collect the neighbors
      // that are going to send stuff to us
      const std::vector<unsigned int> senders
        = Utilities::MPI::compute_point_to_point_communication_pattern
          (this->get_communicator(), destinations);

      // receive ghostcelldata
      std::vector<char> receive;
      CommunicateLocallyMovedVertices::CellInfo<dim,spacedim> cellinfo;
      for (unsigned int i=0; i<senders.size(); ++i)
        {
          MPI_Status status;
          int len;
          MPI_Probe(MPI_ANY_SOURCE, 123, this->get_communicator(), &status);
          MPI_Get_count(&status, MPI_BYTE, &len);
          receive.resize(len);

          char *ptr = &receive[0];
          MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                   this->get_communicator(), &status);

          cellinfo.unpack_data(receive);
          const unsigned int cells = cellinfo.tree_index.size();
          for (unsigned int c=0; c<cells; ++c)
            {
              typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator
              cell (this,
                    0,
                    this->get_p4est_tree_to_coarse_cell_permutation()[cellinfo.tree_index[c]]);

              typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
              dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

              CommunicateLocallyMovedVertices::set_vertices_recursively<dim,spacedim> (*this,
                  p4est_coarse_cell,
                  cell,
                  cellinfo.quadrants[c],
                  cellinfo.first_vertices[c],
                  cellinfo.first_vertex_indices[c]);
            }
        }

      // complete all sends, so that we can
      // safely destroy the buffers.
      if (requests.size() > 0)
        MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

      //check all msgs got sent and received
      Assert(Utilities::MPI::sum(needs_to_get_cells.size(), this->get_communicator())
             == Utilities::MPI::sum(senders.size(), this->get_communicator()),
             ExcInternalError());
    }

    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim,spacedim>::
    register_data_attach (const std::size_t size,
                          const std_cxx11::function<void(const cell_iterator &,
                                                         const CellStatus,
                                                         void *)> &pack_callback)
    {
      Assert(size>0, ExcMessage("register_data_attach(), size==0"));
      Assert(attached_data_pack_callbacks.size()==n_attached_datas,
             ExcMessage("register_data_attach(), not all data has been unpacked last time?"));

      unsigned int offset = attached_data_size+sizeof(CellStatus);
      ++n_attached_datas;
      attached_data_size+=size;
      attached_data_pack_callbacks.push_back(
        std::pair<unsigned int, pack_callback_t> (offset, pack_callback)
      );
      return offset;
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    notify_ready_to_unpack (const unsigned int offset,
                            const std_cxx11::function<void (const cell_iterator &,
                                                            const CellStatus,
                                                            const void *)> &unpack_callback)
    {
      Assert (offset >= sizeof(CellStatus),
              ExcMessage ("invalid offset in notify_ready_to_unpack()"));
      Assert (offset < sizeof(CellStatus)+attached_data_size,
              ExcMessage ("invalid offset in notify_ready_to_unpack()"));
      Assert (n_attached_datas > 0, ExcMessage ("notify_ready_to_unpack() called too often"));

      // Recurse over p4est and hand the caller the data back
      for (typename Triangulation<dim, spacedim>::cell_iterator
           cell = this->begin (0);
           cell != this->end (0);
           ++cell)
        {
          //skip coarse cells, that are not ours
          if (tree_exists_locally<dim, spacedim> (parallel_forest,
                                                  coarse_cell_to_p4est_tree_permutation[cell->index() ])
              == false)
            continue;

          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          typename dealii::internal::p4est::types<dim>::tree *tree =
            init_tree (cell->index());

          dealii::internal::p4est::init_coarse_quadrant<dim> (p4est_coarse_cell);

          // parent_cell is not correct here, but is only used in a refined
          // cell
          post_mesh_data_recursively<dim, spacedim> (*tree,
                                                     cell,
                                                     cell,
                                                     p4est_coarse_cell,
                                                     offset,
                                                     unpack_callback);
        }

      --n_attached_datas;
      if (n_attached_deserialize > 0)
        {
          --n_attached_deserialize;
          attached_data_pack_callbacks.pop_front();
        }

      // important: only remove data if we are not in the deserialization
      // process. There, each SolutionTransfer registers and unpacks before
      // the next one does this, so n_attached_datas is only 1 here.  This
      // would destroy the saved data before the second SolutionTransfer can
      // get it. This created a bug that is documented in
      // tests/mpi/p4est_save_03 with more than one SolutionTransfer.
      if (!n_attached_datas && n_attached_deserialize == 0)
        {
          // everybody got his data, time for cleanup!
          attached_data_size = 0;
          attached_data_pack_callbacks.clear();

          // and release the data
          void *userptr = parallel_forest->user_pointer;
          dealii::internal::p4est::functions<dim>::reset_data (parallel_forest, 0, NULL, NULL);
          parallel_forest->user_pointer = userptr;
        }
    }


    template <int dim, int spacedim>
    const std::vector<types::global_dof_index> &
    Triangulation<dim, spacedim>::get_p4est_tree_to_coarse_cell_permutation() const
    {
      return p4est_tree_to_coarse_cell_permutation;
    }



    template <int dim, int spacedim>
    const std::vector<types::global_dof_index> &
    Triangulation<dim, spacedim>::get_coarse_cell_to_p4est_tree_permutation() const
    {
      return coarse_cell_to_p4est_tree_permutation;
    }



    namespace
    {
      /**
       * This is the callback data structure used to fill
       * vertices_with_ghost_neighbors via the p4est_iterate tool
       */
      template <int dim, int spacedim>
      struct FindGhosts
      {
        typename dealii::parallel::distributed::Triangulation<dim,spacedim> *triangulation;
        sc_array_t *subids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id> >
        *vertices_with_ghost_neighbors;
      };

      /** At a corner (vertex), determine if any of the neighboring cells are
       * ghosts.  If there are, find out their subdomain ids, and if this is a
       * local vertex, then add these subdomain ids to the map
       * vertices_with_ghost_neighbors of that index
       */
      template <int dim, int spacedim>
      void
      find_ghosts_corner
      (typename dealii::internal::p4est::iter<dim>::corner_info *info,
       void *user_data)
      {
        int i, j;
        int nsides = info->sides.elem_count;
        typename dealii::internal::p4est::iter<dim>::corner_side *sides =
          (typename dealii::internal::p4est::iter<dim>::corner_side *)
          (info->sides.array);
        FindGhosts<dim,spacedim> *fg = static_cast<FindGhosts<dim,spacedim> *>(user_data);
        sc_array_t *subids = fg->subids;
        typename dealii::parallel::distributed::Triangulation<dim,spacedim> *triangulation = fg->triangulation;
        int nsubs;
        dealii::types::subdomain_id *subdomain_ids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id> >
        *vertices_with_ghost_neighbors = fg->vertices_with_ghost_neighbors;

        subids->elem_count = 0;
        for (i = 0; i < nsides; i++)
          {
            if (sides[i].is_ghost)
              {
                typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell = cell_from_quad (triangulation, sides[i].treeid, *(sides[i].quad));
                Assert (cell->is_ghost(), ExcMessage ("ghost quad did not find ghost cell"));
                dealii::types::subdomain_id *subid =
                  static_cast<dealii::types::subdomain_id *>(sc_array_push (subids));
                *subid = cell->subdomain_id();
              }
          }

        if (!subids->elem_count)
          {
            return;
          }

        nsubs = (int) subids->elem_count;
        subdomain_ids = (dealii::types::subdomain_id *) (subids->array);

        for (i = 0; i < nsides; i++)
          {
            if (!sides[i].is_ghost)
              {
                typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell = cell_from_quad (triangulation, sides[i].treeid, *(sides[i].quad));

                Assert (!cell->is_ghost(), ExcMessage ("local quad found ghost cell"));

                for (j = 0; j < nsubs; j++)
                  {
                    (*vertices_with_ghost_neighbors)[cell->vertex_index(sides[i].corner)]
                    .insert (subdomain_ids[j]);
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
      find_ghosts_edge
      (typename dealii::internal::p4est::iter<dim>::edge_info *info,
       void *user_data)
      {
        int i, j, k;
        int nsides = info->sides.elem_count;
        typename dealii::internal::p4est::iter<dim>::edge_side *sides =
          (typename dealii::internal::p4est::iter<dim>::edge_side *)
          (info->sides.array);
        FindGhosts<dim,spacedim> *fg = static_cast<FindGhosts<dim,spacedim> *>(user_data);
        sc_array_t *subids = fg->subids;
        typename dealii::parallel::distributed::Triangulation<dim,spacedim> *triangulation = fg->triangulation;
        int nsubs;
        dealii::types::subdomain_id *subdomain_ids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id> >
        *vertices_with_ghost_neighbors = fg->vertices_with_ghost_neighbors;

        subids->elem_count = 0;
        for (i = 0; i < nsides; i++)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < 2; j++)
                  {
                    if (sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell = cell_from_quad (triangulation, sides[i].treeid, *(sides[i].is.hanging.quad[j]));
                        dealii::types::subdomain_id *subid =
                          static_cast<dealii::types::subdomain_id *>(sc_array_push (subids));
                        *subid = cell->subdomain_id();
                      }
                  }
              }
          }

        if (!subids->elem_count)
          {
            return;
          }

        nsubs = (int) subids->elem_count;
        subdomain_ids = (dealii::types::subdomain_id *) (subids->array);

        for (i = 0; i < nsides; i++)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < 2; j++)
                  {
                    if (!sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell = cell_from_quad (triangulation, sides[i].treeid, *(sides[i].is.hanging.quad[j]));

                        for (k = 0; k < nsubs; k++)
                          {
                            (*vertices_with_ghost_neighbors)[cell->vertex_index(p8est_edge_corners[sides[i].edge][1^j])]
                            .insert (subdomain_ids[k]);
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
      find_ghosts_face
      (typename dealii::internal::p4est::iter<dim>::face_info *info,
       void *user_data)
      {
        int i, j, k;
        int nsides = info->sides.elem_count;
        typename dealii::internal::p4est::iter<dim>::face_side *sides =
          (typename dealii::internal::p4est::iter<dim>::face_side *)
          (info->sides.array);
        FindGhosts<dim,spacedim> *fg = static_cast<FindGhosts<dim,spacedim> *>(user_data);
        sc_array_t *subids = fg->subids;
        typename dealii::parallel::distributed::Triangulation<dim,spacedim> *triangulation = fg->triangulation;
        int nsubs;
        dealii::types::subdomain_id *subdomain_ids;
        std::map<unsigned int, std::set<dealii::types::subdomain_id> >
        *vertices_with_ghost_neighbors = fg->vertices_with_ghost_neighbors;
        int limit = (dim == 2) ? 2 : 4;

        subids->elem_count = 0;
        for (i = 0; i < nsides; i++)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < limit; j++)
                  {
                    if (sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell = cell_from_quad (triangulation, sides[i].treeid, *(sides[i].is.hanging.quad[j]));
                        dealii::types::subdomain_id *subid =
                          static_cast<dealii::types::subdomain_id *>(sc_array_push (subids));
                        *subid = cell->subdomain_id();
                      }
                  }
              }
          }

        if (!subids->elem_count)
          {
            return;
          }

        nsubs = (int) subids->elem_count;
        subdomain_ids = (dealii::types::subdomain_id *) (subids->array);

        for (i = 0; i < nsides; i++)
          {
            if (sides[i].is_hanging)
              {
                for (j = 0; j < limit; j++)
                  {
                    if (!sides[i].is.hanging.is_ghost[j])
                      {
                        typename dealii::parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell = cell_from_quad (triangulation, sides[i].treeid, *(sides[i].is.hanging.quad[j]));

                        for (k = 0; k < nsubs; k++)
                          {
                            if (dim == 2)
                              {
                                (*vertices_with_ghost_neighbors)[cell->vertex_index(p4est_face_corners[sides[i].face][(limit - 1)^j])]
                                .insert (subdomain_ids[k]);
                              }
                            else
                              {
                                (*vertices_with_ghost_neighbors)[cell->vertex_index(p8est_face_corners[sides[i].face][(limit - 1)^j])]
                                .insert (subdomain_ids[k]);
                              }
                          }
                      }
                  }
              }
          }

        subids->elem_count = 0;
      }
    }



    /**
     * Determine the neighboring subdomains that are adjacent to each vertex.
     * This is achieved via the p4est_iterate/p8est_iterate tool
     */
    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    fill_vertices_with_ghost_neighbors
    (std::map<unsigned int, std::set<dealii::types::subdomain_id> >
     &vertices_with_ghost_neighbors)
    {
      Assert (dim>1, ExcNotImplemented());

      FindGhosts<dim,spacedim> fg;
      fg.subids = sc_array_new (sizeof (dealii::types::subdomain_id));
      fg.triangulation = this;
      fg.vertices_with_ghost_neighbors = &vertices_with_ghost_neighbors;

      // switch between functions. to make the compiler happy, we need to cast
      // the first two arguments to the type p[48]est_iterate wants to see. this
      // cast is the identity cast in each of the two branches, so it is safe.
      switch (dim)
        {
        case 2:
          p4est_iterate (reinterpret_cast<dealii::internal::p4est::types<2>::forest *>(this->parallel_forest),
                         reinterpret_cast<dealii::internal::p4est::types<2>::ghost *>(this->parallel_ghost),
                         static_cast<void *>(&fg),
                         NULL, find_ghosts_face<2,spacedim>, find_ghosts_corner<2,spacedim>);
          break;

        case 3:
          p8est_iterate (reinterpret_cast<dealii::internal::p4est::types<3>::forest *>(this->parallel_forest),
                         reinterpret_cast<dealii::internal::p4est::types<3>::ghost *>(this->parallel_ghost),
                         static_cast<void *>(&fg),
                         NULL, find_ghosts_face<3,spacedim>, find_ghosts_edge<3,spacedim>, find_ghosts_corner<3,spacedim>);
          break;

        default:
          Assert (false, ExcNotImplemented());
        }

      sc_array_destroy (fg.subids);
    }



    /**
     * Determine the neighboring subdomains that are adjacent to each vertex
     * on the given multigrid level
     */
    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    fill_level_vertices_with_ghost_neighbors
    (const int level,
     std::map<unsigned int, std::set<dealii::types::subdomain_id> >
     &vertices_with_ghost_neighbors)
    {
      const std::vector<bool> locally_active_vertices =
        mark_locally_active_vertices_on_level(level);
      cell_iterator cell = this->begin(level),
                    endc = this->end(level);
      for ( ; cell != endc; ++cell)
        if (cell->level_subdomain_id() != dealii::numbers::artificial_subdomain_id
            && cell->level_subdomain_id() != this->locally_owned_subdomain())
          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
            if (locally_active_vertices[cell->vertex_index(v)])
              vertices_with_ghost_neighbors[cell->vertex_index(v)]
              .insert (cell->level_subdomain_id());

      //now for the vertices on periodic faces
      typename std::map<std::pair<cell_iterator, unsigned int>,
               std::pair<std::pair<cell_iterator,unsigned int>, std::bitset<3> > >::const_iterator it;

      for (it = this->get_periodic_face_map().begin(); it!= this->get_periodic_face_map().end(); ++it)
        {
          const cell_iterator &cell_1 = it->first.first;
          const unsigned int face_no_1 = it->first.second;
          const cell_iterator &cell_2 = it->second.first.first;
          const unsigned int face_no_2 = it->second.first.second;
          const std::bitset<3> face_orientation = it->second.second;

          if (cell_1->level() == level &&
              cell_2->level() == level)
            {
              for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
                {
                  // take possible non-standard orientation of faces into account
                  const unsigned int vface0 =
                    GeometryInfo<dim>::standard_to_real_face_vertex(v,face_orientation[0],
                                                                    face_orientation[1],
                                                                    face_orientation[2]);
                  const unsigned int idx0 = cell_1->face(face_no_1)->vertex_index(vface0);
                  const unsigned int idx1 = cell_2->face(face_no_2)->vertex_index(v);
                  if (vertices_with_ghost_neighbors.find(idx0) != vertices_with_ghost_neighbors.end())
                    vertices_with_ghost_neighbors[idx1].insert(vertices_with_ghost_neighbors[idx0].begin(),
                                                               vertices_with_ghost_neighbors[idx0].end());
                  if (vertices_with_ghost_neighbors.find(idx1) != vertices_with_ghost_neighbors.end())
                    vertices_with_ghost_neighbors[idx0].insert(vertices_with_ghost_neighbors[idx1].begin(),
                                                               vertices_with_ghost_neighbors[idx1].end());
                }
            }
        }
    }



    template<int dim, int spacedim>
    std::vector<bool>
    Triangulation<dim,spacedim>
    ::mark_locally_active_vertices_on_level (const int level) const
    {
      Assert (dim>1, ExcNotImplemented());

      std::vector<bool> marked_vertices(this->n_vertices(), false);
      cell_iterator cell = this->begin(level),
                    endc = this->end(level);
      for ( ; cell != endc; ++cell)
        if (cell->level_subdomain_id() == this->locally_owned_subdomain())
          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
            marked_vertices[cell->vertex_index(v)] = true;

      /**
       * ensure that if one of the two vertices on a periodic face is marked
       * as active (i.e., belonging to an owned level cell), also the other
       * one is active
       */
      typename std::map<std::pair<cell_iterator, unsigned int>,
               std::pair<std::pair<cell_iterator,unsigned int>, std::bitset<3> > >::const_iterator it;

      for (it = this->get_periodic_face_map().begin(); it!= this->get_periodic_face_map().end(); ++it)
        {
          const cell_iterator &cell_1 = it->first.first;
          const unsigned int face_no_1 = it->first.second;
          const cell_iterator &cell_2 = it->second.first.first;
          const unsigned int face_no_2 = it->second.first.second;
          const std::bitset<3> &face_orientation = it->second.second;

          if (cell_1->level() == level &&
              cell_2->level() == level)
            {
              for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
                {
                  // take possible non-standard orientation of faces into account
                  const unsigned int vface0 =
                    GeometryInfo<dim>::standard_to_real_face_vertex(v,face_orientation[0],
                                                                    face_orientation[1],
                                                                    face_orientation[2]);
                  if (marked_vertices[cell_1->face(face_no_1)->vertex_index(vface0)] ||
                      marked_vertices[cell_2->face(face_no_2)->vertex_index(v)])
                    marked_vertices[cell_1->face(face_no_1)->vertex_index(vface0)]
                      = marked_vertices[cell_2->face(face_no_2)->vertex_index(v)]
                        = true;
                }
            }
        }

      return marked_vertices;
    }



    template<int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::add_periodicity
    (const std::vector<GridTools::PeriodicFacePair<cell_iterator> > &
     periodicity_vector)
    {
#if DEAL_II_P4EST_VERSION_GTE(0,3,4,1)
      Assert (triangulation_has_content == true,
              ExcMessage ("The triangulation is empty!"));
      Assert (this->n_levels() == 1,
              ExcMessage ("The triangulation is refined!"));

      typedef std::vector<GridTools::PeriodicFacePair<cell_iterator> >
      FaceVector;
      typename FaceVector::const_iterator it, periodic_end;
      it = periodicity_vector.begin();
      periodic_end = periodicity_vector.end();

      for (; it<periodic_end; ++it)
        {
          const cell_iterator first_cell = it->cell[0];
          const cell_iterator second_cell = it->cell[1];
          const unsigned int face_left = it->face_idx[0];
          const unsigned int face_right = it->face_idx[1];

          //respective cells of the matching faces in p4est
          const unsigned int tree_left
            = coarse_cell_to_p4est_tree_permutation[std::distance(this->begin(),
                                                                  first_cell)];
          const unsigned int tree_right
            = coarse_cell_to_p4est_tree_permutation[std::distance(this->begin(),
                                                                  second_cell)];

          // p4est wants to know which corner the first corner on
          // the face with the lower id is mapped to on the face with
          // with the higher id. For d==2 there are only two possibilities
          // that are determined by it->orientation[1].
          // For d==3 we have to use GridTools::OrientationLookupTable.
          // The result is given below.

          unsigned int p4est_orientation = 0;
          if (dim==2)
            p4est_orientation = it->orientation[1];
          else
            {
              const unsigned int face_idx_list[] = {face_left, face_right};
              const cell_iterator cell_list[] = {first_cell, second_cell};
              unsigned int lower_idx, higher_idx;
              if (face_left<=face_right)
                {
                  higher_idx = 1;
                  lower_idx = 0;
                }
              else
                {
                  higher_idx = 0;
                  lower_idx = 1;
                }

              // get the cell index of the first index on the face with the lower id
              unsigned int first_p4est_idx_on_cell = p8est_face_corners[face_idx_list[lower_idx]][0];
              unsigned int first_dealii_idx_on_face = numbers::invalid_unsigned_int;
              for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
                {
                  const unsigned int first_dealii_idx_on_cell
                    =  GeometryInfo<dim>::face_to_cell_vertices
                       (face_idx_list[lower_idx], i,
                        cell_list[lower_idx]->face_orientation(face_idx_list[lower_idx]),
                        cell_list[lower_idx]->face_flip(face_idx_list[lower_idx]),
                        cell_list[lower_idx]->face_rotation(face_idx_list[lower_idx]));
                  if (first_p4est_idx_on_cell == first_dealii_idx_on_cell)
                    {
                      first_dealii_idx_on_face = i;
                      break;
                    }
                }
              Assert( first_dealii_idx_on_face != numbers::invalid_unsigned_int, ExcInternalError());
              // Now map dealii_idx_on_face according to the orientation
              const unsigned int left_to_right [8][4] = {{0,2,1,3},{0,1,2,3},{3,1,2,0},{3,2,1,0},
                {2,3,0,1},{1,3,0,2},{1,0,3,2},{2,0,3,1}
              };
              const unsigned int right_to_left [8][4] = {{0,2,1,3},{0,1,2,3},{3,1,2,0},{3,2,1,0},
                {2,3,0,1},{2,0,3,1},{1,0,3,2},{1,3,0,2}
              };
              const unsigned int second_dealii_idx_on_face
                = lower_idx==0?left_to_right[it->orientation.to_ulong()][first_dealii_idx_on_face]:
                  right_to_left[it->orientation.to_ulong()][first_dealii_idx_on_face];
              const unsigned int second_dealii_idx_on_cell
                = GeometryInfo<dim>::face_to_cell_vertices
                  (face_idx_list[higher_idx], second_dealii_idx_on_face,
                   cell_list[higher_idx]->face_orientation(face_idx_list[higher_idx]),
                   cell_list[higher_idx]->face_flip(face_idx_list[higher_idx]),
                   cell_list[higher_idx]->face_rotation(face_idx_list[higher_idx]));
              //map back to p4est
              const unsigned int second_p4est_idx_on_face
                = p8est_corner_face_corners[second_dealii_idx_on_cell][face_idx_list[higher_idx]];
              p4est_orientation = second_p4est_idx_on_face;
            }

          dealii::internal::p4est::functions<dim>::
          connectivity_join_faces (connectivity,
                                   tree_left,
                                   tree_right,
                                   face_left,
                                   face_right,
                                   p4est_orientation);
        }


      Assert(dealii::internal::p4est::functions<dim>::connectivity_is_valid
             (connectivity) == 1, ExcInternalError());

      // now create a forest out of the connectivity data structure
      dealii::internal::p4est::functions<dim>::destroy (parallel_forest);
      parallel_forest
        = dealii::internal::p4est::functions<dim>::
          new_forest (this->mpi_communicator,
                      connectivity,
                      /* minimum initial number of quadrants per tree */ 0,
                      /* minimum level of upfront refinement */ 0,
                      /* use uniform upfront refinement */ 1,
                      /* user_data_size = */ 0,
                      /* user_data_constructor = */ NULL,
                      /* user_pointer */ this);


      try
        {
          copy_local_forest_to_triangulation ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

      //finally call the base class for storing the periodicity information
      dealii::Triangulation<dim, spacedim>::add_periodicity(periodicity_vector);
#else
      Assert(false, ExcMessage ("Need p4est version >= 0.3.4.1!"));
#endif
    }



    template <int dim, int spacedim>
    std::size_t
    Triangulation<dim,spacedim>::memory_consumption () const
    {
      std::size_t mem=
        this->dealii::parallel::Triangulation<dim,spacedim>::memory_consumption()
        + MemoryConsumption::memory_consumption(triangulation_has_content)
        + MemoryConsumption::memory_consumption(connectivity)
        + MemoryConsumption::memory_consumption(parallel_forest)
        + MemoryConsumption::memory_consumption(refinement_in_progress)
        + MemoryConsumption::memory_consumption(attached_data_size)
        + MemoryConsumption::memory_consumption(n_attached_datas)
//      + MemoryConsumption::memory_consumption(attached_data_pack_callbacks) //TODO[TH]: how?
        + MemoryConsumption::memory_consumption(coarse_cell_to_p4est_tree_permutation)
        + MemoryConsumption::memory_consumption(p4est_tree_to_coarse_cell_permutation)
        + memory_consumption_p4est();

      return mem;
    }



    template <int dim, int spacedim>
    std::size_t
    Triangulation<dim,spacedim>::memory_consumption_p4est () const
    {
      return dealii::internal::p4est::functions<dim>::forest_memory_used(parallel_forest)
             + dealii::internal::p4est::functions<dim>::connectivity_memory_used(connectivity);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    copy_triangulation (const dealii::Triangulation<dim, spacedim> &old_tria)
    {
      try
        {
          dealii::Triangulation<dim,spacedim>::
          copy_triangulation (old_tria);
        }
      catch (const typename dealii::Triangulation<dim,spacedim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

      // note that now we have some content in the p4est objects and call the
      // functions that do the actual work (which are dimension dependent, so
      // separate)
      triangulation_has_content = true;

      Assert (old_tria.n_levels() == 1,
              ExcMessage ("Parallel distributed triangulations can only be copied, "
                          "if they are not refined!"));

      if (const dealii::parallel::distributed::Triangulation<dim,spacedim> *
          old_tria_x = dynamic_cast<const dealii::parallel::distributed::Triangulation<dim,spacedim> *>(&old_tria))
        {
          Assert (!old_tria_x->refinement_in_progress,
                  ExcMessage ("Parallel distributed triangulations can only "
                              "be copied, if no refinement is in progress!"));

          // duplicate MPI communicator, stored in the base class
          dealii::parallel::Triangulation<dim,spacedim>::copy_triangulation (old_tria);

          coarse_cell_to_p4est_tree_permutation = old_tria_x->coarse_cell_to_p4est_tree_permutation;
          p4est_tree_to_coarse_cell_permutation = old_tria_x->p4est_tree_to_coarse_cell_permutation;
          attached_data_size = old_tria_x->attached_data_size;
          n_attached_datas   = old_tria_x->n_attached_datas;

          settings           = old_tria_x->settings;
        }
      else
        {
          setup_coarse_cell_to_p4est_tree_permutation ();
        }

      copy_new_triangulation_to_p4est (dealii::internal::int2type<dim>());

      try
        {
          copy_local_forest_to_triangulation ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }

      this->update_number_cache ();
      this->update_periodic_face_map();
    }


    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::
    attach_mesh_data()
    {
      // determine size of memory in bytes to attach to each cell. This needs
      // to be constant because of p4est.
      if (attached_data_size==0)
        {
          Assert(n_attached_datas==0, ExcInternalError());

          //nothing to do
          return;
        }

      // realloc user_data in p4est
      void *userptr = parallel_forest->user_pointer;
      dealii::internal::p4est::functions<dim>::reset_data (parallel_forest,
                                                           attached_data_size+sizeof(CellStatus),
                                                           NULL, NULL);
      parallel_forest->user_pointer = userptr;


      // Recurse over p4est and Triangulation
      // to find refined/coarsened/kept
      // cells. Then query and attach the data.
      for (typename Triangulation<dim,spacedim>::cell_iterator
           cell = this->begin(0);
           cell != this->end(0);
           ++cell)
        {
          //skip coarse cells, that are not ours
          if (tree_exists_locally<dim,spacedim>(parallel_forest,
                                                coarse_cell_to_p4est_tree_permutation[cell->index()])
              == false)
            continue;

          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          typename dealii::internal::p4est::types<dim>::tree *tree =
            init_tree(cell->index());

          dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

          attach_mesh_data_recursively<dim,spacedim>(*tree,
                                                     cell,
                                                     p4est_coarse_cell,
                                                     attached_data_pack_callbacks);
        }
    }

    template <int dim, int spacedim>
    std::vector<unsigned int>
    Triangulation<dim,spacedim>::
    get_cell_weights()
    {
      // Allocate the space for the weights. In fact we do not know yet, how
      // many cells we own after the refinement (only p4est knows that
      // at this point). We simply reserve n_active_cells space and if many
      // more cells are refined than coarsened than additional reallocation
      // will be done inside get_cell_weights_recursively.
      std::vector<unsigned int> weights;
      weights.reserve(this->n_active_cells());

      // Recurse over p4est and Triangulation
      // to find refined/coarsened/kept
      // cells. Then append cell_weight.
      // Note that we need to follow the p4est ordering
      // instead of the deal.II ordering to get the cell_weights
      // in the same order p4est will encounter them during repartitioning.
      for (unsigned int c=0; c<this->n_cells(0); ++c)
        {
          // skip coarse cells, that are not ours
          if (tree_exists_locally<dim,spacedim>(parallel_forest,c) == false)
            continue;

          const unsigned int coarse_cell_index =
            p4est_tree_to_coarse_cell_permutation[c];

          const typename Triangulation<dim,spacedim>::cell_iterator
          dealii_coarse_cell (this, 0, coarse_cell_index);

          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          dealii::internal::p4est::functions<dim>::
          quadrant_set_morton (&p4est_coarse_cell,
                               /*level=*/0,
                               /*index=*/0);
          p4est_coarse_cell.p.which_tree = c;

          const typename dealii::internal::p4est::types<dim>::tree *tree =
            init_tree(coarse_cell_index);

          get_cell_weights_recursively<dim,spacedim>(*tree,
                                                     dealii_coarse_cell,
                                                     p4est_coarse_cell,
                                                     this->signals,
                                                     weights);
        }

      return weights;
    }

    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::cell_iterator
    cell_from_quad
    (typename dealii::parallel::distributed::Triangulation<dim,spacedim> *triangulation,
     typename dealii::internal::p4est::types<dim>::topidx treeidx,
     typename dealii::internal::p4est::types<dim>::quadrant &quad)
    {
      int i, l = quad.level;
      int child_id;
      types::global_dof_index dealii_index =
        triangulation->get_p4est_tree_to_coarse_cell_permutation()[treeidx];

      for (i = 0; i < l; i++)
        {
          typename dealii::Triangulation<dim,spacedim>::cell_iterator cell (triangulation, i, dealii_index);
          child_id = dealii::internal::p4est::functions<dim>::quadrant_ancestor_id (&quad, i + 1);
          Assert (cell->has_children (), ExcMessage ("p4est quadrant does not correspond to a cell!"));
          dealii_index = cell->child_index(child_id);
        }

      typename dealii::Triangulation<dim,spacedim>::cell_iterator out_cell (triangulation, l, dealii_index);

      return out_cell;
    }



    template <int spacedim>
    Triangulation<1,spacedim>::Triangulation (MPI_Comm)
      :
      dealii::parallel::Triangulation<1,spacedim>(MPI_COMM_WORLD,
                                                  typename dealii::Triangulation<1,spacedim>::MeshSmoothing(),
                                                  false)
    {
      Assert (false, ExcNotImplemented());
    }


    template <int spacedim>
    Triangulation<1,spacedim>::~Triangulation ()
    {
      Assert (false, ExcNotImplemented());
    }



    template <int spacedim>
    void
    Triangulation<1,spacedim>::communicate_locally_moved_vertices
    (const std::vector<bool> &/*vertex_locally_moved*/)
    {
      Assert (false, ExcNotImplemented());
    }


    template <int spacedim>
    const std::vector<types::global_dof_index> &
    Triangulation<1,spacedim>::get_p4est_tree_to_coarse_cell_permutation() const
    {
      static std::vector<types::global_dof_index> a;
      return a;
    }

    template <int spacedim>
    void
    Triangulation<1,spacedim>::
    fill_vertices_with_ghost_neighbors
    (std::map<unsigned int, std::set<dealii::types::subdomain_id> >
     &/*vertices_with_ghost_neighbors*/)
    {
      Assert (false, ExcNotImplemented());
    }

    template <int spacedim>
    void
    Triangulation<1,spacedim>::
    fill_level_vertices_with_ghost_neighbors
    (const unsigned int /*level*/,
     std::map<unsigned int, std::set<dealii::types::subdomain_id> >
     &/*vertices_with_ghost_neighbors*/)
    {
      Assert (false, ExcNotImplemented());
    }

    template <int spacedim>
    std::vector<bool>
    Triangulation<1,spacedim>::
    mark_locally_active_vertices_on_level (const unsigned int) const
    {
      Assert (false, ExcNotImplemented());
      return std::vector<bool>();
    }
  }
}


#else // DEAL_II_WITH_P4EST

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::Triangulation ()
      :
      dealii::parallel::Triangulation<dim,spacedim>(MPI_COMM_SELF)
    {
      Assert (false, ExcNotImplemented());
    }
  }
}

#endif // DEAL_II_WITH_P4EST




/*-------------- Explicit Instantiations -------------------------------*/
#include "tria.inst"


DEAL_II_NAMESPACE_CLOSE
