#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria.h>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_P4EST

namespace internal
{
  namespace p4est
  {
    namespace
    {


    }


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



    const unsigned int functions<2>::max_level;




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

    const unsigned int functions<3>::max_level;

    template <int dim>
    void iterate(typename dealii::internal::p4est::types<dim>::forest *parallel_forest,
                 typename dealii::internal::p4est::types<dim>::ghost *parallel_ghost,
                 std::function<void (typename dealii::internal::p4est::iter<dim>::face_info *info)> face_callback,
                 std::function<void (typename dealii::internal::p4est::iter<dim>::edge_info *info)> edge_callback,
                 std::function<void (typename dealii::internal::p4est::iter<dim>::corner_info *info)> corner_callback)
    {
      Assert(false, ExcNotImplemented());
    }

    template <>
    void iterate<2>(typename dealii::internal::p4est::types<2>::forest *parallel_forest,
                    typename dealii::internal::p4est::types<2>::ghost *parallel_ghost,
                    std::function<void (typename dealii::internal::p4est::iter<2>::face_info *info)> face_callback,
                    std::function<void (typename dealii::internal::p4est::iter<3>::edge_info *info)> /*edge_callback*/,
                    std::function<void (typename dealii::internal::p4est::iter<2>::corner_info *info)> corner_callback)
    {
      struct user_data_t
      {
        std::function<void (typename dealii::internal::p4est::iter<2>::face_info *info)> face_callback;
        std::function<void (typename dealii::internal::p4est::iter<2>::corner_info *info)> corner_callback;

      } user_data;
      user_data.face_callback = face_callback;
      user_data.corner_callback = corner_callback;

      p4est_iterate (parallel_forest, parallel_ghost, &user_data,
                     nullptr,
                     [] (typename dealii::internal::p4est::iter<2>::face_info *info, void *user_data)
      {
        auto my_user_data = static_cast<user_data_t *> (user_data);
        my_user_data->face_callback(info);
      },
      [] (typename dealii::internal::p4est::iter<2>::corner_info *info, void *user_data)
      {
        auto my_user_data = static_cast<user_data_t *> (user_data);
        my_user_data->corner_callback(info);
      });
    }

    template <>
    void iterate<3>(typename dealii::internal::p4est::types<3>::forest *parallel_forest,
                    typename dealii::internal::p4est::types<3>::ghost *parallel_ghost,
                    std::function<void (typename dealii::internal::p4est::iter<3>::face_info *info)> face_callback,
                    std::function<void (typename dealii::internal::p4est::iter<3>::edge_info *info)> edge_callback,
                    std::function<void (typename dealii::internal::p4est::iter<3>::corner_info *info)> corner_callback)
    {
      struct user_data_t
      {
        std::function<void (typename dealii::internal::p4est::iter<3>::face_info *info)> face_callback;
        std::function<void (typename dealii::internal::p4est::iter<3>::edge_info *info)> edge_callback;
        std::function<void (typename dealii::internal::p4est::iter<3>::corner_info *info)> corner_callback;

      } user_data;
      user_data.face_callback = face_callback;
      user_data.edge_callback = edge_callback;
      user_data.corner_callback = corner_callback;

      p8est_iterate (parallel_forest, parallel_ghost, &user_data,
                     nullptr,
                     [] (typename dealii::internal::p4est::iter<3>::face_info *info, void *user_data)
      {
        auto my_user_data = static_cast<user_data_t *> (user_data);
        my_user_data->face_callback(info);
      },
      [] (typename dealii::internal::p4est::iter<3>::edge_info *info, void *user_data)
      {
        auto my_user_data = static_cast<user_data_t *> (user_data);
        my_user_data->edge_callback(info);
      },
      [] (typename dealii::internal::p4est::iter<3>::corner_info *info, void *user_data)
      {
        auto my_user_data = static_cast<user_data_t *> (user_data);
        my_user_data->corner_callback(info);
      });

    }



    template <int dim>
    void
    init_quadrant_children
    (const typename types<dim>::quadrant &p4est_cell,
     typename types<dim>::quadrant (&p4est_children)[dealii::GeometryInfo<dim>::max_children_per_cell])
    {

      for (unsigned int c=0;
           c<dealii::GeometryInfo<dim>::max_children_per_cell; ++c)
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

    template <int dim>
    bool
    tree_exists_locally (const typename types<dim>::forest *parallel_forest,
                         const typename types<dim>::topidx coarse_grid_cell)
    {
      Assert (coarse_grid_cell < parallel_forest->connectivity->num_trees,
              ExcInternalError());
      return ((coarse_grid_cell >= parallel_forest->first_local_tree)
              &&
              (coarse_grid_cell <= parallel_forest->last_local_tree));
    }
  }
}

#endif //DEAL_II_WITH_P4EST

/*-------------- Explicit Instantiations -------------------------------*/
#include "p4est_wrappers.inst"


DEAL_II_NAMESPACE_CLOSE
