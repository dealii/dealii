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

#ifndef dealii_p4est_wrappers_h
#define dealii_p4est_wrappers_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>

#ifdef DEAL_II_WITH_P4EST

#  include <deal.II/base/mpi_stub.h>

#  include <p4est_bits.h>
#  include <p4est_communication.h>
#  include <p4est_extended.h>
#  include <p4est_ghost.h>
#  include <p4est_iterate.h>
#  include <p4est_search.h>
#  include <p4est_vtk.h>
#  include <p8est_bits.h>
#  include <p8est_communication.h>
#  include <p8est_extended.h>
#  include <p8est_ghost.h>
#  include <p8est_iterate.h>
#  include <p8est_search.h>
#  include <p8est_vtk.h>

#  include <limits>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#  ifndef DOXYGEN
namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    class Triangulation;
  }
} // namespace parallel
#  endif

namespace internal
{
  namespace p4est
  {
    /**
     * A structure whose explicit specializations contain alias to the
     * relevant p4est_* and p8est_* types. Using this structure, for example
     * by saying <tt>types<dim>::connectivity</tt> we can write code in a
     * dimension independent way, either referring to p4est_connectivity_t or
     * p8est_connectivity_t, depending on template argument.
     */
    template <int>
    struct types;

    // these struct mimics p4est for 1d
    template <>
    struct types<1>
    {
      // id of a quadrant is an integeger
      using quadrant = int;

      // maximum number of children
      static const int max_n_child_indices_bits = 27;

      // number of bits the data type of id has
      static const int n_bits = std::numeric_limits<quadrant>::digits;
    };

    template <>
    struct types<2>
    {
      using connectivity              = p4est_connectivity_t;
      using forest                    = p4est_t;
      using tree                      = p4est_tree_t;
      using quadrant                  = p4est_quadrant_t;
      using quadrant_coord            = p4est_qcoord_t;
      using topidx                    = p4est_topidx_t;
      using locidx                    = p4est_locidx_t;
      using gloidx                    = p4est_gloidx_t;
      using balance_type              = p4est_connect_type_t;
      using ghost                     = p4est_ghost_t;
      using transfer_context          = p4est_transfer_context_t;
      using search_partition_callback = p4est_search_partition_t;
    };

    template <>
    struct types<3>
    {
      using connectivity              = p8est_connectivity_t;
      using forest                    = p8est_t;
      using tree                      = p8est_tree_t;
      using quadrant                  = p8est_quadrant_t;
      using quadrant_coord            = p4est_qcoord_t;
      using topidx                    = p4est_topidx_t;
      using locidx                    = p4est_locidx_t;
      using gloidx                    = p4est_gloidx_t;
      using balance_type              = p8est_connect_type_t;
      using ghost                     = p8est_ghost_t;
      using transfer_context          = p8est_transfer_context_t;
      using search_partition_callback = p8est_search_partition_t;
    };



    /**
     * A structure whose explicit specializations represent the
     * relevant p4est_* and p8est_* functions. Using this structure, for
     * example by saying functions<dim>::quadrant_compare(...), we can write
     * code in a dimension independent way, either calling
     * p4est_quadrant_compare or p8est_quadrant_compare, depending on template
     * argument.
     *
     * In most cases, the members of this class are simply pointers to
     * p4est_* or p8est_* functions. In one case, it's simply a static
     * member function that dispatches to things p4est chooses to
     * implement via a macro.
     */
    template <int dim>
    struct functions;

    template <>
    struct functions<2>
    {
      static int (&quadrant_compare)(const void *v1, const void *v2);

      static void (&quadrant_childrenv)(const types<2>::quadrant *q,
                                        types<2>::quadrant        c[]);

      static int (&quadrant_overlaps_tree)(types<2>::tree           *tree,
                                           const types<2>::quadrant *q);

      static void (&quadrant_set_morton)(types<2>::quadrant *quadrant,
                                         int                 level,
                                         std::uint64_t       id);

      static void
      quadrant_init(types<2>::quadrant &q);

      static int (&quadrant_is_equal)(const types<2>::quadrant *q1,
                                      const types<2>::quadrant *q2);

      static int (&quadrant_is_sibling)(const types<2>::quadrant *q1,
                                        const types<2>::quadrant *q2);

      static int (&quadrant_is_ancestor)(const types<2>::quadrant *q1,
                                         const types<2>::quadrant *q2);

      static int (&quadrant_ancestor_id)(const types<2>::quadrant *q,
                                         int                       level);

      static int (&comm_find_owner)(types<2>::forest         *p4est,
                                    const types<2>::locidx    which_tree,
                                    const types<2>::quadrant *q,
                                    const int                 guess);

      static types<2>::connectivity *(&connectivity_new)(
        types<2>::topidx num_vertices,
        types<2>::topidx num_trees,
        types<2>::topidx num_corners,
        types<2>::topidx num_vtt);

      static types<2>::connectivity *(&connectivity_new_copy)(
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
        const int8_t           *ctc);

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
        std::size_t             data_size,
        p4est_init_t            init_fn,
        void                   *user_pointer);

      static types<2>::forest *(&copy_forest)(types<2>::forest *input,
                                              int               copy_data);

      static void (&destroy)(types<2>::forest *p4est);

      static void (&refine)(types<2>::forest *p4est,
                            int               refine_recursive,
                            p4est_refine_t    refine_fn,
                            p4est_init_t      init_fn);

      static void (&coarsen)(types<2>::forest *p4est,
                             int               coarsen_recursive,
                             p4est_coarsen_t   coarsen_fn,
                             p4est_init_t      init_fn);

      static void (&balance)(types<2>::forest      *p4est,
                             types<2>::balance_type btype,
                             p4est_init_t           init_fn);

      static types<2>::gloidx (&partition)(types<2>::forest *p4est,
                                           int partition_for_coarsening,
                                           p4est_weight_t weight_fn);

      static void (&save)(const char       *filename,
                          types<2>::forest *p4est,
                          int               save_data);

      static types<2>::forest *(&load_ext)(const char *filename,
                                           MPI_Comm    mpicomm,
                                           std::size_t data_size,
                                           int         load_data,
                                           int         autopartition,
                                           int         broadcasthead,
                                           void       *user_pointer,
                                           types<2>::connectivity **p4est);

      static int (&connectivity_save)(const char             *filename,
                                      types<2>::connectivity *connectivity);

      static int (&connectivity_is_valid)(types<2>::connectivity *connectivity);

      static types<2>::connectivity *(&connectivity_load)(const char  *filename,
                                                          std::size_t *length);

      static unsigned int (&checksum)(types<2>::forest *p4est);

      static void (&vtk_write_file)(types<2>::forest *p4est,
                                    p4est_geometry_t *,
                                    const char *baseName);

      static types<2>::ghost *(&ghost_new)(types<2>::forest      *p4est,
                                           types<2>::balance_type btype);

      static void (&ghost_destroy)(types<2>::ghost *ghost);

      static void (&reset_data)(types<2>::forest *p4est,
                                std::size_t       data_size,
                                p4est_init_t      init_fn,
                                void             *user_pointer);

      static std::size_t (&forest_memory_used)(types<2>::forest *p4est);

      static std::size_t (&connectivity_memory_used)(
        types<2>::connectivity *p4est);

      template <int spacedim>
      static void
      iterate(dealii::internal::p4est::types<2>::forest *parallel_forest,
              dealii::internal::p4est::types<2>::ghost  *parallel_ghost,
              void                                      *user_data);

      static constexpr unsigned int max_level = P4EST_MAXLEVEL;

      static void (&transfer_fixed)(const types<2>::gloidx *dest_gfq,
                                    const types<2>::gloidx *src_gfq,
                                    MPI_Comm                mpicomm,
                                    int                     tag,
                                    void                   *dest_data,
                                    const void             *src_data,
                                    std::size_t             data_size);

      static types<2>::transfer_context *(&transfer_fixed_begin)(
        const types<2>::gloidx *dest_gfq,
        const types<2>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void                   *dest_data,
        const void             *src_data,
        std::size_t             data_size);

      static void (&transfer_fixed_end)(types<2>::transfer_context *tc);

      static void (&transfer_custom)(const types<2>::gloidx *dest_gfq,
                                     const types<2>::gloidx *src_gfq,
                                     MPI_Comm                mpicomm,
                                     int                     tag,
                                     void                   *dest_data,
                                     const int              *dest_sizes,
                                     const void             *src_data,
                                     const int              *src_sizes);

      static types<2>::transfer_context *(&transfer_custom_begin)(
        const types<2>::gloidx *dest_gfq,
        const types<2>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void                   *dest_data,
        const int              *dest_sizes,
        const void             *src_data,
        const int              *src_sizes);

      static void (&transfer_custom_end)(types<2>::transfer_context *tc);

      static void (&search_partition)(
        types<2>::forest                   *forest,
        int                                 call_post,
        types<2>::search_partition_callback quadrant_fn,
        types<2>::search_partition_callback point_fn,
        sc_array_t                         *points);

      static void (&quadrant_coord_to_vertex)(
        types<2>::connectivity  *connectivity,
        types<2>::topidx         treeid,
        types<2>::quadrant_coord x,
        types<2>::quadrant_coord y,
        double                   vxyz[3]);
    };


    template <>
    struct functions<3>
    {
      static int (&quadrant_compare)(const void *v1, const void *v2);

      static void (&quadrant_childrenv)(const types<3>::quadrant *q,
                                        types<3>::quadrant        c[]);

      static int (&quadrant_overlaps_tree)(types<3>::tree           *tree,
                                           const types<3>::quadrant *q);

      static void (&quadrant_set_morton)(types<3>::quadrant *quadrant,
                                         int                 level,
                                         std::uint64_t       id);

      static void
      quadrant_init(types<3>::quadrant &q);

      static int (&quadrant_is_equal)(const types<3>::quadrant *q1,
                                      const types<3>::quadrant *q2);

      static int (&quadrant_is_sibling)(const types<3>::quadrant *q1,
                                        const types<3>::quadrant *q2);

      static int (&quadrant_is_ancestor)(const types<3>::quadrant *q1,
                                         const types<3>::quadrant *q2);
      static int (&quadrant_ancestor_id)(const types<3>::quadrant *q,
                                         int                       level);

      static int (&comm_find_owner)(types<3>::forest         *p4est,
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

      static types<3>::connectivity *(&connectivity_new_copy)(
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
        const int8_t           *ctc);

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
        std::size_t             data_size,
        p8est_init_t            init_fn,
        void                   *user_pointer);

      static types<3>::forest *(&copy_forest)(types<3>::forest *input,
                                              int               copy_data);

      static void (&destroy)(types<3>::forest *p8est);

      static void (&refine)(types<3>::forest *p8est,
                            int               refine_recursive,
                            p8est_refine_t    refine_fn,
                            p8est_init_t      init_fn);

      static void (&coarsen)(types<3>::forest *p8est,
                             int               coarsen_recursive,
                             p8est_coarsen_t   coarsen_fn,
                             p8est_init_t      init_fn);

      static void (&balance)(types<3>::forest      *p8est,
                             types<3>::balance_type btype,
                             p8est_init_t           init_fn);

      static types<3>::gloidx (&partition)(types<3>::forest *p8est,
                                           int partition_for_coarsening,
                                           p8est_weight_t weight_fn);

      static void (&save)(const char       *filename,
                          types<3>::forest *p4est,
                          int               save_data);

      static types<3>::forest *(&load_ext)(const char *filename,
                                           MPI_Comm    mpicomm,
                                           std::size_t data_size,
                                           int         load_data,
                                           int         autopartition,
                                           int         broadcasthead,
                                           void       *user_pointer,
                                           types<3>::connectivity **p4est);

      static int (&connectivity_save)(const char             *filename,
                                      types<3>::connectivity *connectivity);

      static int (&connectivity_is_valid)(types<3>::connectivity *connectivity);

      static types<3>::connectivity *(&connectivity_load)(const char  *filename,
                                                          std::size_t *length);

      static unsigned int (&checksum)(types<3>::forest *p8est);

      static void (&vtk_write_file)(types<3>::forest *p8est,
                                    p8est_geometry_t *,
                                    const char *baseName);
      static types<3>::ghost *(&ghost_new)(types<3>::forest      *p4est,
                                           types<3>::balance_type btype);

      static void (&ghost_destroy)(types<3>::ghost *ghost);

      static void (&reset_data)(types<3>::forest *p4est,
                                std::size_t       data_size,
                                p8est_init_t      init_fn,
                                void             *user_pointer);

      static std::size_t (&forest_memory_used)(types<3>::forest *p4est);

      static std::size_t (&connectivity_memory_used)(
        types<3>::connectivity *p4est);

      static constexpr unsigned int max_level = P8EST_MAXLEVEL;

      static void (&transfer_fixed)(const types<3>::gloidx *dest_gfq,
                                    const types<3>::gloidx *src_gfq,
                                    MPI_Comm                mpicomm,
                                    int                     tag,
                                    void                   *dest_data,
                                    const void             *src_data,
                                    std::size_t             data_size);

      static types<3>::transfer_context *(&transfer_fixed_begin)(
        const types<3>::gloidx *dest_gfq,
        const types<3>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void                   *dest_data,
        const void             *src_data,
        std::size_t             data_size);

      static void (&transfer_fixed_end)(types<3>::transfer_context *tc);

      static void (&transfer_custom)(const types<3>::gloidx *dest_gfq,
                                     const types<3>::gloidx *src_gfq,
                                     MPI_Comm                mpicomm,
                                     int                     tag,
                                     void                   *dest_data,
                                     const int              *dest_sizes,
                                     const void             *src_data,
                                     const int              *src_sizes);

      static types<3>::transfer_context *(&transfer_custom_begin)(
        const types<3>::gloidx *dest_gfq,
        const types<3>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void                   *dest_data,
        const int              *dest_sizes,
        const void             *src_data,
        const int              *src_sizes);

      static void (&transfer_custom_end)(types<3>::transfer_context *tc);

      static void (&search_partition)(
        types<3>::forest                   *forest,
        int                                 call_post,
        types<3>::search_partition_callback quadrant_fn,
        types<3>::search_partition_callback point_fn,
        sc_array_t                         *points);

      static void (&quadrant_coord_to_vertex)(
        types<3>::connectivity  *connectivity,
        types<3>::topidx         treeid,
        types<3>::quadrant_coord x,
        types<3>::quadrant_coord y,
        types<3>::quadrant_coord z,
        double                   vxyz[3]);
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
      using corner_info = p4est_iter_corner_info_t;
      using corner_side = p4est_iter_corner_side_t;
      using corner_iter = p4est_iter_corner_t;
      using face_info   = p4est_iter_face_info_t;
      using face_side   = p4est_iter_face_side_t;
      using face_iter   = p4est_iter_face_t;
    };

    template <>
    struct iter<3>
    {
      using corner_info = p8est_iter_corner_info_t;
      using corner_side = p8est_iter_corner_side_t;
      using corner_iter = p8est_iter_corner_t;
      using edge_info   = p8est_iter_edge_info_t;
      using edge_side   = p8est_iter_edge_side_t;
      using edge_iter   = p8est_iter_edge_t;
      using face_info   = p8est_iter_face_info_t;
      using face_side   = p8est_iter_face_side_t;
      using face_iter   = p8est_iter_face_t;
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
     * Deep copy a p4est connectivity object.
     */
    template <int dim>
    typename types<dim>::connectivity *
    copy_connectivity(const typename types<dim>::connectivity *connectivity);

#  ifndef DOXYGEN
    template <>
    typename types<2>::connectivity *
    copy_connectivity<2>(const typename types<2>::connectivity *connectivity);

    template <>
    typename types<3>::connectivity *
    copy_connectivity<3>(const typename types<3>::connectivity *connectivity);
#  endif
  } // namespace p4est
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_P4EST

#endif // dealii_p4est_wrappers_h
