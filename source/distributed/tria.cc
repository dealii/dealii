// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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
#include <deal.II/lac/sparsity_pattern.h>
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

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
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

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
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

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
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

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
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
                   * gives the face number and ttf / 6 the face orientation
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
      return; //this quadrant and none of it's childs belongs to us.

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

        //mark other childs as invalid, so that unpack only happens once
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
        //it's children got coarsened into this cell
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
      // this quadrant and none of it's children belong to us.
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
                          const types::subdomain_id                   my_subdomain,
                          typename internal::p4est::types<dim>::forest &forest);

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

    typename internal::p4est::types<dim>::forest forest;

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
                        const types::subdomain_id                    my_subdomain,
                        typename internal::p4est::types<dim>::forest &forest)
    :
    forest(forest)
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
    Triangulation<dim,spacedim>::NumberCache::NumberCache()
      :
      n_global_active_cells(0),
      n_global_levels(0)
    {}



    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::
    Triangulation (MPI_Comm mpi_communicator,
                   const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid,
                   const Settings settings_)
      :
      // do not check for distorted cells
      dealii::Triangulation<dim,spacedim>
      (smooth_grid,
       false),
      mpi_communicator (Utilities::MPI::
                        duplicate_communicator(mpi_communicator)),
      settings(settings_),
      my_subdomain (Utilities::MPI::this_mpi_process (this->mpi_communicator)),
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

      number_cache.n_locally_owned_active_cells
      .resize (Utilities::MPI::n_mpi_processes (mpi_communicator));
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

      // get rid of the unique communicator used here again
      MPI_Comm_free (&mpi_communicator);
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

      update_number_cache ();
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

      update_number_cache ();
    }

    template <int dim, int spacedim>
    bool
    Triangulation<dim,spacedim>::has_hanging_nodes () const
    {
      // if there are any active cells with level less than n_global_levels()-1, then
      // there is obviously also one with level n_global_levels()-1, and
      // consequently there must be a hanging node somewhere.
      //
      // the problem is that we cannot just ask for the first active cell, but
      // instead need to filter over locally owned cells
      bool res_local = false;
      for (typename Triangulation<dim, spacedim>::active_cell_iterator cell = this->begin_active();
           (cell != this->end()) && (cell->level() < (int)(n_global_levels()-1));
           cell++)
        if (cell->is_locally_owned())
          {
            res_local = true;
            break;
          }

      // reduce over MPI
      bool res;
      MPI_Allreduce(&res_local, &res, 1, MPI::BOOL, MPI_LOR, mpi_communicator);
      return res;
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::setup_coarse_cell_to_p4est_tree_permutation ()
    {
      SparsityPattern cell_connectivity;
      GridTools::get_face_connectivity_of_cells (*this, cell_connectivity);
      coarse_cell_to_p4est_tree_permutation.resize (this->n_cells(0));
      SparsityTools::
      reorder_Cuthill_McKee (cell_connectivity,
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

      if (my_subdomain==0)
        {
          std::string fname=std::string(filename)+".info";
          std::ofstream f(fname.c_str());
          f << "version nproc attached_bytes n_attached_objs n_coarse_cells" << std::endl
            << 2 << " "
            << Utilities::MPI::n_mpi_processes (mpi_communicator) << " "
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
      AssertThrow(numcpus <= Utilities::MPI::n_mpi_processes (mpi_communicator),
                  ExcMessage("parallel::distributed::Triangulation::load() only supports loading "
                             "saved data with a greater or equal number of processes than were used to "
                             "save() when using p4est 0.3.4.2."));
#endif

      attached_data_size = 0;
      n_attached_datas = 0;
      n_attached_deserialize = attached_count;

#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      parallel_forest = dealii::internal::p4est::functions<dim>::load_ext (
                          filename, mpi_communicator,
                          attached_size, attached_size>0,
                          autopartition, 0,
                          this,
                          &connectivity);
#else
      parallel_forest = dealii::internal::p4est::functions<dim>::load (
                          filename, mpi_communicator,
                          attached_size, attached_size>0,
                          this,
                          &connectivity);
#endif
      if (numcpus != Utilities::MPI::n_mpi_processes (mpi_communicator))
        {
          // We are changing the number of CPUs so we need to repartition.
          // Note that p4est actually distributes the cells between the changed
          // number of CPUs and so everything works without this call, but
          // this command changes the distribution for some reason, so we
          // will leave it in here.
          dealii::internal::p4est::functions<dim>::
          partition (parallel_forest,
                     /* prepare coarsening */ 1,
                     /* weight_callback */ NULL);


        }

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

      update_number_cache ();
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
          new_forest (mpi_communicator,
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
          new_forest (mpi_communicator,
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
          new_forest (mpi_communicator,
                      connectivity,
                      /* minimum initial number of quadrants per tree */ 0,
                      /* minimum level of upfront refinement */ 0,
                      /* use uniform upfront refinement */ 1,
                      /* user_data_size = */ 0,
                      /* user_data_constructor = */ NULL,
                      /* user_pointer */ this);
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
                                                        my_subdomain);
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
          if (cell->subdomain_id() != my_subdomain
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
          for (unsigned int lvl=this->n_levels(); lvl>0;)
            {
              --lvl;
              for (typename Triangulation<dim,spacedim>::cell_iterator cell = this->begin(lvl); cell!=this->end(lvl); ++cell)
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

          for (unsigned int lvl=0; lvl<this->n_levels(); ++lvl)
            {
              marked_vertices[lvl].resize(this->n_vertices(), false);

              for (typename dealii::Triangulation<dim,spacedim>::cell_iterator
                   cell = this->begin(lvl);
                   cell != this->end(lvl); ++cell)
                if (cell->level_subdomain_id() == this->locally_owned_subdomain())
                  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
                    marked_vertices[lvl][cell->vertex_index(v)] = true;
            }



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
                                                                      my_subdomain,
                                                                      marked_vertices);
            }

          //step 3: make sure we have the parent of our level cells
          for (unsigned int lvl=this->n_levels(); lvl>0;)
            {
              --lvl;
              for (typename Triangulation<dim,spacedim>::cell_iterator cell = this->begin(lvl); cell!=this->end(lvl); ++cell)
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

          //step 4: Special case: on each level we need all the face neighbors
          // of our own level cells these are normally on the same level,
          // unless the neighbor is active and coarser. It can end up on a
          // different processor. Luckily, the level_subdomain_id can be
          // figured out without communication, because the cell is active
          // (and so level_subdomain_id=subdomain_id):
          for (typename Triangulation<dim,spacedim>::cell_iterator cell = this->begin(); cell!=this->end(); ++cell)
            {
              if (cell->level_subdomain_id()!=this->locally_owned_subdomain())
                continue;

              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                  if (cell->face(f)->at_boundary())
                    continue;
                  if (cell->neighbor(f)->level() < cell->level()
                      &&
                      cell->neighbor(f)->level_subdomain_id()!=this->locally_owned_subdomain())
                    {
                      Assert(cell->neighbor(f)->active(), ExcInternalError());
                      Assert(cell->neighbor(f)->subdomain_id() != numbers::artificial_subdomain_id, ExcInternalError());
                      Assert(cell->neighbor(f)->level_subdomain_id() == numbers::artificial_subdomain_id
                             || cell->neighbor(f)->level_subdomain_id() == cell->neighbor(f)->subdomain_id(), ExcInternalError());
                      cell->neighbor(f)->set_level_subdomain_id(cell->neighbor(f)->subdomain_id());
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

        if (Utilities::MPI::n_mpi_processes (mpi_communicator) == 1)
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
            if (cell->subdomain_id() == my_subdomain)
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
                               my_subdomain,
                               *parallel_forest);

      // copy refine and coarsen flags into p4est and execute the refinement
      // and coarsening. this uses the refine_and_coarsen_list just built,
      // which is communicated to the callback functions through the
      // user_pointer
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

      // partition the new mesh between all processors
      dealii::internal::p4est::functions<dim>::
      partition (parallel_forest,
                 /* prepare coarsening */ 1,
                 /* weight_callback */ NULL);

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


      refinement_in_progress = false;

      update_number_cache ();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::update_number_cache ()
    {
      Assert (number_cache.n_locally_owned_active_cells.size()
              ==
              Utilities::MPI::n_mpi_processes (mpi_communicator),
              ExcInternalError());

      std::fill (number_cache.n_locally_owned_active_cells.begin(),
                 number_cache.n_locally_owned_active_cells.end(),
                 0);

      if (this->n_levels() > 0)
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
             cell = this->begin_active();
             cell != this->end(); ++cell)
          if (cell->subdomain_id() == my_subdomain)
            ++number_cache.n_locally_owned_active_cells[my_subdomain];

      unsigned int send_value
        = number_cache.n_locally_owned_active_cells[my_subdomain];
      MPI_Allgather (&send_value,
                     1,
                     MPI_UNSIGNED,
                     &number_cache.n_locally_owned_active_cells[0],
                     1,
                     MPI_UNSIGNED,
                     mpi_communicator);

      number_cache.n_global_active_cells
        = std::accumulate (number_cache.n_locally_owned_active_cells.begin(),
                           number_cache.n_locally_owned_active_cells.end(),
                           /* ensure sum is computed with correct data type:*/
                           static_cast<types::global_dof_index>(0));
      number_cache.n_global_levels = Utilities::MPI::max(this->n_levels(), mpi_communicator);
    }



    template <int dim, int spacedim>
    types::subdomain_id
    Triangulation<dim,spacedim>::locally_owned_subdomain () const
    {
      Assert (dim > 1, ExcNotImplemented());
      return my_subdomain;
    }



    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim,spacedim>::n_locally_owned_active_cells () const
    {
      return number_cache.n_locally_owned_active_cells[my_subdomain];
    }



    template <int dim, int spacedim>
    types::global_dof_index
    Triangulation<dim,spacedim>::n_global_active_cells () const
    {
      return number_cache.n_global_active_cells;
    }



    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim,spacedim>::n_global_levels () const
    {
      return number_cache.n_global_levels;
    }



    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    Triangulation<dim,spacedim>::n_locally_owned_active_cells_per_processor () const
    {
      return number_cache.n_locally_owned_active_cells;
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

    namespace
    {
      /**
       * This is the callback data structure used to fill
       * vertices_with_ghost_neighbors via the p4est_iterate tool
       */
      template <int dim, int spacedim>
      struct find_ghosts
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
        struct find_ghosts<dim,spacedim> *fg = static_cast<struct find_ghosts<dim,spacedim> *>(user_data);
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
        struct find_ghosts<dim,spacedim> *fg = static_cast<struct find_ghosts<dim,spacedim> *>(user_data);
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
        struct find_ghosts<dim,spacedim> *fg = static_cast<struct find_ghosts<dim,spacedim> *>(user_data);
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

      struct find_ghosts<dim,spacedim> fg;

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


    template <int dim, int spacedim>
    MPI_Comm
    Triangulation<dim,spacedim>::get_communicator () const
    {
      return mpi_communicator;
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

          //TODO Add support for non default orientation.
          Assert(it->orientation == 1,
                 ExcMessage("Found a face match with non standard orientation. "
                            "This function is only suitable for meshes with "
                            "cells in default orientation"));

          dealii::internal::p4est::functions<dim>::
          connectivity_join_faces (connectivity,
                                   tree_left,
                                   tree_right,
                                   face_left,
                                   face_right,
                                   /* orientation */ 0);

          /* The orientation parameter above describes the difference between
           * the cell on the left and the cell on the right would number of the
           * corners of the face.  In the periodic domains most users will want,
           * this orientation will be 0, i.e. the two cells would number the
           * corners the same way.  More exotic periodicity, like moebius strips
           * or converting an unstructured quad/hex mesh into a periodic domain,
           * are not supported right now, and undefined behavior will occur if
           * users try to make them periodic.  This may be addressed at a later
           * date.
           */
        }


      Assert(dealii::internal::p4est::functions<dim>::connectivity_is_valid
             (connectivity) == 1, ExcInternalError());

      // now create a forest out of the connectivity data structure
      dealii::internal::p4est::functions<dim>::destroy (parallel_forest);
      parallel_forest
        = dealii::internal::p4est::functions<dim>::
          new_forest (mpi_communicator,
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

#else
      Assert(false, ExcMessage ("Need p4est version >= 0.3.4.1!"));
#endif
    }



    template <int dim, int spacedim>
    std::size_t
    Triangulation<dim,spacedim>::memory_consumption () const
    {
      std::size_t mem=
        this->dealii::Triangulation<dim,spacedim>::memory_consumption()
        + MemoryConsumption::memory_consumption(mpi_communicator)
        + MemoryConsumption::memory_consumption(my_subdomain)
        + MemoryConsumption::memory_consumption(triangulation_has_content)
        + MemoryConsumption::memory_consumption(number_cache.n_locally_owned_active_cells)
        + MemoryConsumption::memory_consumption(number_cache.n_global_active_cells)
        + MemoryConsumption::memory_consumption(number_cache.n_global_levels)
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

          mpi_communicator = Utilities::MPI::duplicate_communicator (old_tria_x->get_communicator ());

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

      update_number_cache ();
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
    typename dealii::Triangulation<dim,spacedim>::cell_iterator
    cell_from_quad
    (typename dealii::parallel::distributed::Triangulation<dim,spacedim> *triangulation,
     typename dealii::internal::p4est::types<dim>::topidx treeidx,
     typename dealii::internal::p4est::types<dim>::quadrant &quad)
    {
      int i, l = quad.level;
      int child_id;
      const std::vector<types::global_dof_index> perm = triangulation->get_p4est_tree_to_coarse_cell_permutation ();
      types::global_dof_index dealii_index = perm[treeidx];

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
    {
      Assert (false, ExcNotImplemented());
    }


    template <int spacedim>
    Triangulation<1,spacedim>::~Triangulation ()
    {
      Assert (false, ExcNotImplemented());
    }



    template <int spacedim>
    types::subdomain_id
    Triangulation<1,spacedim>::locally_owned_subdomain () const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }


    template <int spacedim>
    types::global_dof_index
    Triangulation<1,spacedim>::n_global_active_cells () const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }


    template <int spacedim>
    unsigned int
    Triangulation<1,spacedim>::n_global_levels () const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }


    template <int spacedim>
    MPI_Comm
    Triangulation<1,spacedim>::get_communicator () const
    {
      return MPI_COMM_WORLD;
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
     &vertices_with_ghost_neighbors)
    {
      Assert (false, ExcNotImplemented());
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
    {
      Assert (false, ExcNotImplemented());
    }


    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::~Triangulation ()
    {
      Assert (false, ExcNotImplemented());
    }



    template <int dim, int spacedim>
    types::subdomain_id
    Triangulation<dim,spacedim>::locally_owned_subdomain () const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }


#ifdef DEAL_II_WITH_MPI
    template <int dim, int spacedim>
    MPI_Comm
    Triangulation<dim,spacedim>::get_communicator () const
    {
      Assert (false, ExcNotImplemented());
      return MPI_COMM_WORLD;
    }
#endif
  }
}

#endif // DEAL_II_WITH_P4EST




/*-------------- Explicit Instantiations -------------------------------*/
#include "tria.inst"


DEAL_II_NAMESPACE_CLOSE
