// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>


DEAL_II_NAMESPACE_OPEN


#ifdef DEAL_II_WITH_P4EST

namespace
{
  template <int dim, int spacedim>
  void
  get_vertex_to_cell_mappings(
    const Triangulation<dim, spacedim> &triangulation,
    std::vector<unsigned int> &         vertex_touch_count,
    std::vector<std::list<
      std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                unsigned int>>> &       vertex_to_cell)
  {
    vertex_touch_count.resize(triangulation.n_vertices());
    vertex_to_cell.resize(triangulation.n_vertices());

    for (const auto &cell : triangulation.active_cell_iterators())
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        {
          ++vertex_touch_count[cell->vertex_index(v)];
          vertex_to_cell[cell->vertex_index(v)].emplace_back(cell, v);
        }
  }



  template <int dim, int spacedim>
  void
  get_edge_to_cell_mappings(
    const Triangulation<dim, spacedim> &triangulation,
    std::vector<unsigned int> &         edge_touch_count,
    std::vector<std::list<
      std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                unsigned int>>> &       edge_to_cell)
  {
    Assert(triangulation.n_levels() == 1, ExcInternalError());

    edge_touch_count.resize(triangulation.n_active_lines());
    edge_to_cell.resize(triangulation.n_active_lines());

    for (const auto &cell : triangulation.active_cell_iterators())
      for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
        {
          ++edge_touch_count[cell->line(l)->index()];
          edge_to_cell[cell->line(l)->index()].emplace_back(cell, l);
        }
  }



  /**
   * Set all vertex and cell related information in the p4est connectivity
   * structure.
   */
  template <int dim, int spacedim>
  void
  set_vertex_and_cell_info(
    const Triangulation<dim, spacedim> &triangulation,
    const std::vector<unsigned int> &   vertex_touch_count,
    const std::vector<std::list<
      std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                unsigned int>>> &       vertex_to_cell,
    const std::vector<types::global_dof_index>
      &        coarse_cell_to_p4est_tree_permutation,
    const bool set_vertex_info,
    typename internal::p4est::types<dim>::connectivity *connectivity)
  {
    // copy the vertices into the connectivity structure. the triangulation
    // exports the array of vertices, but some of the entries are sometimes
    // unused; this shouldn't be the case for a newly created triangulation,
    // but make sure
    //
    // note that p4est stores coordinates as a triplet of values even in 2d
    Assert(triangulation.get_used_vertices().size() ==
             triangulation.get_vertices().size(),
           ExcInternalError());
    Assert(std::find(triangulation.get_used_vertices().begin(),
                     triangulation.get_used_vertices().end(),
                     false) == triangulation.get_used_vertices().end(),
           ExcInternalError());
    if (set_vertex_info == true)
      for (unsigned int v = 0; v < triangulation.n_vertices(); ++v)
        {
          connectivity->vertices[3 * v] = triangulation.get_vertices()[v][0];
          connectivity->vertices[3 * v + 1] =
            triangulation.get_vertices()[v][1];
          connectivity->vertices[3 * v + 2] =
            (spacedim == 2 ? 0 : triangulation.get_vertices()[v][2]);
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
        const unsigned int index =
          coarse_cell_to_p4est_tree_permutation[cell->index()];

        for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
          {
            if (set_vertex_info == true)
              connectivity
                ->tree_to_vertex[index * GeometryInfo<dim>::vertices_per_cell +
                                 v] = cell->vertex_index(v);
            connectivity
              ->tree_to_corner[index * GeometryInfo<dim>::vertices_per_cell +
                               v] = cell->vertex_index(v);
          }

        // neighborship information. if a cell is at a boundary, then enter
        // the index of the cell itself here
        for (auto f : GeometryInfo<dim>::face_indices())
          if (cell->face(f)->at_boundary() == false)
            connectivity
              ->tree_to_tree[index * GeometryInfo<dim>::faces_per_cell + f] =
              coarse_cell_to_p4est_tree_permutation[cell->neighbor(f)->index()];
          else
            connectivity
              ->tree_to_tree[index * GeometryInfo<dim>::faces_per_cell + f] =
              coarse_cell_to_p4est_tree_permutation[cell->index()];

        // fill tree_to_face, which is essentially neighbor_to_neighbor;
        // however, we have to remap the resulting face number as well
        for (auto f : GeometryInfo<dim>::face_indices())
          if (cell->face(f)->at_boundary() == false)
            {
              switch (dim)
                {
                  case 2:
                    {
                      connectivity->tree_to_face
                        [index * GeometryInfo<dim>::faces_per_cell + f] =
                        cell->neighbor_of_neighbor(f);
                      break;
                    }

                  case 3:
                    {
                      /*
                       * The values for tree_to_face are in 0..23 where ttf % 6
                       * gives the face number and ttf / 4 the face orientation
                       * code.  The orientation is determined as follows.  Let
                       * my_face and other_face be the two face numbers of the
                       * connecting trees in 0..5.  Then the first face vertex
                       * of the lower of my_face and other_face connects to a
                       * face vertex numbered 0..3 in the higher of my_face and
                       * other_face.  The face orientation is defined as this
                       * number.  If my_face == other_face, treating either of
                       * both faces as the lower one leads to the same result.
                       */

                      connectivity->tree_to_face[index * 6 + f] =
                        cell->neighbor_of_neighbor(f);

                      unsigned int face_idx_list[2] = {
                        f, cell->neighbor_of_neighbor(f)};
                      typename Triangulation<dim>::active_cell_iterator
                                   cell_list[2] = {cell, cell->neighbor(f)};
                      unsigned int smaller_idx  = 0;

                      if (f > cell->neighbor_of_neighbor(f))
                        smaller_idx = 1;

                      unsigned int larger_idx = (smaller_idx + 1) % 2;
                      // smaller = *_list[smaller_idx]
                      // larger = *_list[larger_idx]

                      unsigned int v = 0;

                      // global vertex index of vertex 0 on face of cell with
                      // smaller local face index
                      unsigned int g_idx = cell_list[smaller_idx]->vertex_index(
                        GeometryInfo<dim>::face_to_cell_vertices(
                          face_idx_list[smaller_idx],
                          0,
                          cell_list[smaller_idx]->face_orientation(
                            face_idx_list[smaller_idx]),
                          cell_list[smaller_idx]->face_flip(
                            face_idx_list[smaller_idx]),
                          cell_list[smaller_idx]->face_rotation(
                            face_idx_list[smaller_idx])));

                      // loop over vertices on face from other cell and compare
                      // global vertex numbers
                      for (unsigned int i = 0;
                           i < GeometryInfo<dim>::vertices_per_face;
                           ++i)
                        {
                          unsigned int idx =
                            cell_list[larger_idx]->vertex_index(
                              GeometryInfo<dim>::face_to_cell_vertices(
                                face_idx_list[larger_idx], i));

                          if (idx == g_idx)
                            {
                              v = i;
                              break;
                            }
                        }

                      connectivity->tree_to_face[index * 6 + f] += 6 * v;
                      break;
                    }

                  default:
                    Assert(false, ExcNotImplemented());
                }
            }
          else
            connectivity
              ->tree_to_face[index * GeometryInfo<dim>::faces_per_cell + f] = f;
      }

    // now fill the vertex information
    connectivity->ctt_offset[0] = 0;
    std::partial_sum(vertex_touch_count.begin(),
                     vertex_touch_count.end(),
                     &connectivity->ctt_offset[1]);

    const typename internal::p4est::types<dim>::locidx num_vtt =
      std::accumulate(vertex_touch_count.begin(), vertex_touch_count.end(), 0u);
    (void)num_vtt;
    Assert(connectivity->ctt_offset[triangulation.n_vertices()] == num_vtt,
           ExcInternalError());

    for (unsigned int v = 0; v < triangulation.n_vertices(); ++v)
      {
        Assert(vertex_to_cell[v].size() == vertex_touch_count[v],
               ExcInternalError());

        typename std::list<
          std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                    unsigned int>>::const_iterator p =
          vertex_to_cell[v].begin();
        for (unsigned int c = 0; c < vertex_touch_count[v]; ++c, ++p)
          {
            connectivity->corner_to_tree[connectivity->ctt_offset[v] + c] =
              coarse_cell_to_p4est_tree_permutation[p->first->index()];
            connectivity->corner_to_corner[connectivity->ctt_offset[v] + c] =
              p->second;
          }
      }
  }



  template <int dim, int spacedim>
  bool
  tree_exists_locally(
    const typename internal::p4est::types<dim>::forest *parallel_forest,
    const typename internal::p4est::types<dim>::topidx  coarse_grid_cell)
  {
    Assert(coarse_grid_cell < parallel_forest->connectivity->num_trees,
           ExcInternalError());
    return ((coarse_grid_cell >= parallel_forest->first_local_tree) &&
            (coarse_grid_cell <= parallel_forest->last_local_tree));
  }


  template <int dim, int spacedim>
  void
  delete_all_children_and_self(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  {
    if (cell->has_children())
      for (unsigned int c = 0; c < cell->n_children(); ++c)
        delete_all_children_and_self<dim, spacedim>(cell->child(c));
    else
      cell->set_coarsen_flag();
  }



  template <int dim, int spacedim>
  void
  delete_all_children(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  {
    if (cell->has_children())
      for (unsigned int c = 0; c < cell->n_children(); ++c)
        delete_all_children_and_self<dim, spacedim>(cell->child(c));
  }


  template <int dim, int spacedim>
  void
  determine_level_subdomain_id_recursively(
    const typename internal::p4est::types<dim>::tree &          tree,
    const typename internal::p4est::types<dim>::locidx &        tree_index,
    const typename Triangulation<dim, spacedim>::cell_iterator &dealii_cell,
    const typename internal::p4est::types<dim>::quadrant &      p4est_cell,
    typename internal::p4est::types<dim>::forest &              forest,
    const types::subdomain_id                                   my_subdomain,
    const std::vector<std::vector<bool>> &                      marked_vertices)
  {
    if (dealii_cell->level_subdomain_id() == numbers::artificial_subdomain_id)
      {
        // important: only assign the level_subdomain_id if it is a ghost cell
        // even though we could fill in all.
        bool used = false;
        for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
          {
            if (marked_vertices[dealii_cell->level()]
                               [dealii_cell->vertex_index(v)])
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
        if (!used && dealii_cell->is_active() &&
            dealii_cell->is_artificial() == false &&
            dealii_cell->level() + 1 < static_cast<int>(marked_vertices.size()))
          {
            for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
              {
                if (marked_vertices[dealii_cell->level() + 1]
                                   [dealii_cell->vertex_index(v)])
                  {
                    used = true;
                    break;
                  }
              }
          }

        // Like above, but now the other way around
        if (!used && dealii_cell->is_active() &&
            dealii_cell->is_artificial() == false && dealii_cell->level() > 0)
          {
            for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
              {
                if (marked_vertices[dealii_cell->level() - 1]
                                   [dealii_cell->vertex_index(v)])
                  {
                    used = true;
                    break;
                  }
              }
          }

        if (used)
          {
            int owner = internal::p4est::functions<dim>::comm_find_owner(
              &forest, tree_index, &p4est_cell, my_subdomain);
            Assert((owner != -2) && (owner != -1),
                   ExcMessage("p4est should know the owner."));
            dealii_cell->set_level_subdomain_id(owner);
          }
      }

    if (dealii_cell->has_children())
      {
        typename internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          switch (dim)
            {
              case 2:
                P4EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              case 3:
                P8EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }


        internal::p4est::functions<dim>::quadrant_childrenv(&p4est_cell,
                                                            p4est_child);

        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          {
            determine_level_subdomain_id_recursively<dim, spacedim>(
              tree,
              tree_index,
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
  match_tree_recursively(
    const typename internal::p4est::types<dim>::tree &          tree,
    const typename Triangulation<dim, spacedim>::cell_iterator &dealii_cell,
    const typename internal::p4est::types<dim>::quadrant &      p4est_cell,
    const typename internal::p4est::types<dim>::forest &        forest,
    const types::subdomain_id                                   my_subdomain)
  {
    // check if this cell exists in the local p4est cell
    if (sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                         &p4est_cell,
                         internal::p4est::functions<dim>::quadrant_compare) !=
        -1)
      {
        // yes, cell found in local part of p4est
        delete_all_children<dim, spacedim>(dealii_cell);
        if (dealii_cell->is_active())
          dealii_cell->set_subdomain_id(my_subdomain);
      }
    else
      {
        // no, cell not found in local part of p4est. this means that the
        // local part is more refined than the current cell. if this cell has
        // no children of its own, we need to refine it, and if it does
        // already have children then loop over all children and see if they
        // are locally available as well
        if (dealii_cell->is_active())
          dealii_cell->set_refine_flag();
        else
          {
            typename internal::p4est::types<dim>::quadrant
              p4est_child[GeometryInfo<dim>::max_children_per_cell];
            for (unsigned int c = 0;
                 c < GeometryInfo<dim>::max_children_per_cell;
                 ++c)
              switch (dim)
                {
                  case 2:
                    P4EST_QUADRANT_INIT(&p4est_child[c]);
                    break;
                  case 3:
                    P8EST_QUADRANT_INIT(&p4est_child[c]);
                    break;
                  default:
                    Assert(false, ExcNotImplemented());
                }


            internal::p4est::functions<dim>::quadrant_childrenv(&p4est_cell,
                                                                p4est_child);

            for (unsigned int c = 0;
                 c < GeometryInfo<dim>::max_children_per_cell;
                 ++c)
              if (internal::p4est::functions<dim>::quadrant_overlaps_tree(
                    const_cast<typename internal::p4est::types<dim>::tree *>(
                      &tree),
                    &p4est_child[c]) == false)
                {
                  // no, this child is locally not available in the p4est.
                  // delete all its children but, because this may not be
                  // successful, make sure to mark all children recursively
                  // as not local.
                  delete_all_children<dim, spacedim>(dealii_cell->child(c));
                  dealii_cell->child(c)->recursively_set_subdomain_id(
                    numbers::artificial_subdomain_id);
                }
              else
                {
                  // at least some part of the tree rooted in this child is
                  // locally available
                  match_tree_recursively<dim, spacedim>(tree,
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
  match_quadrant(
    const dealii::Triangulation<dim, spacedim> *          tria,
    unsigned int                                          dealii_index,
    const typename internal::p4est::types<dim>::quadrant &ghost_quadrant,
    types::subdomain_id                                   ghost_owner)
  {
    const int l = ghost_quadrant.level;

    for (int i = 0; i < l; ++i)
      {
        typename Triangulation<dim, spacedim>::cell_iterator cell(tria,
                                                                  i,
                                                                  dealii_index);
        if (cell->is_active())
          {
            cell->clear_coarsen_flag();
            cell->set_refine_flag();
            return;
          }

        const int child_id =
          internal::p4est::functions<dim>::quadrant_ancestor_id(&ghost_quadrant,
                                                                i + 1);
        dealii_index = cell->child_index(child_id);
      }

    typename Triangulation<dim, spacedim>::cell_iterator cell(tria,
                                                              l,
                                                              dealii_index);
    if (cell->has_children())
      delete_all_children<dim, spacedim>(cell);
    else
      {
        cell->clear_coarsen_flag();
        cell->set_subdomain_id(ghost_owner);
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
    RefineAndCoarsenList(const Triangulation<dim, spacedim> &triangulation,
                         const std::vector<types::global_dof_index>
                           &p4est_tree_to_coarse_cell_permutation,
                         const types::subdomain_id my_subdomain);

    /**
     * A callback function that we pass to the p4est data structures when a
     * forest is to be refined. The p4est functions call it back with a tree
     * (the index of the tree that grows out of a given coarse cell) and a
     * refinement path from that coarse cell to a terminal/leaf cell. The
     * function returns whether the corresponding cell in the deal.II
     * triangulation has the refined flag set.
     */
    static int
    refine_callback(
      typename internal::p4est::types<dim>::forest *  forest,
      typename internal::p4est::types<dim>::topidx    coarse_cell_index,
      typename internal::p4est::types<dim>::quadrant *quadrant);

    /**
     * Same as the refine_callback function, but return whether all four of
     * the given children of a non-terminal cell are to be coarsened away.
     */
    static int
    coarsen_callback(
      typename internal::p4est::types<dim>::forest *  forest,
      typename internal::p4est::types<dim>::topidx    coarse_cell_index,
      typename internal::p4est::types<dim>::quadrant *children[]);

    bool
    pointers_are_at_end() const;

  private:
    std::vector<typename internal::p4est::types<dim>::quadrant> refine_list;
    typename std::vector<typename internal::p4est::types<dim>::quadrant>::
      const_iterator current_refine_pointer;

    std::vector<typename internal::p4est::types<dim>::quadrant> coarsen_list;
    typename std::vector<typename internal::p4est::types<dim>::quadrant>::
      const_iterator current_coarsen_pointer;

    void
    build_lists(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename internal::p4est::types<dim>::quadrant &      p4est_cell,
      const types::subdomain_id                                   myid);
  };



  template <int dim, int spacedim>
  bool
  RefineAndCoarsenList<dim, spacedim>::pointers_are_at_end() const
  {
    return ((current_refine_pointer == refine_list.end()) &&
            (current_coarsen_pointer == coarsen_list.end()));
  }



  template <int dim, int spacedim>
  RefineAndCoarsenList<dim, spacedim>::RefineAndCoarsenList(
    const Triangulation<dim, spacedim> &triangulation,
    const std::vector<types::global_dof_index>
      &                       p4est_tree_to_coarse_cell_permutation,
    const types::subdomain_id my_subdomain)
  {
    // count how many flags are set and allocate that much memory
    unsigned int n_refine_flags = 0, n_coarsen_flags = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        // skip cells that are not local
        if (cell->subdomain_id() != my_subdomain)
          continue;

        if (cell->refine_flag_set())
          ++n_refine_flags;
        else if (cell->coarsen_flag_set())
          ++n_coarsen_flags;
      }

    refine_list.reserve(n_refine_flags);
    coarsen_list.reserve(n_coarsen_flags);


    // now build the lists of cells that are flagged. note that p4est will
    // traverse its cells in the order in which trees appear in the
    // forest. this order is not the same as the order of coarse cells in the
    // deal.II Triangulation because we have translated everything by the
    // coarse_cell_to_p4est_tree_permutation permutation. in order to make
    // sure that the output array is already in the correct order, traverse
    // our coarse cells in the same order in which p4est will:
    for (unsigned int c = 0; c < triangulation.n_cells(0); ++c)
      {
        unsigned int coarse_cell_index =
          p4est_tree_to_coarse_cell_permutation[c];

        const typename Triangulation<dim, spacedim>::cell_iterator cell(
          &triangulation, 0, coarse_cell_index);

        typename internal::p4est::types<dim>::quadrant p4est_cell;
        internal::p4est::functions<dim>::quadrant_set_morton(&p4est_cell,
                                                             /*level=*/0,
                                                             /*index=*/0);
        p4est_cell.p.which_tree = c;
        build_lists(cell, p4est_cell, my_subdomain);
      }


    Assert(refine_list.size() == n_refine_flags, ExcInternalError());
    Assert(coarsen_list.size() == n_coarsen_flags, ExcInternalError());

    // make sure that our ordering in fact worked
    for (unsigned int i = 1; i < refine_list.size(); ++i)
      Assert(refine_list[i].p.which_tree >= refine_list[i - 1].p.which_tree,
             ExcInternalError());
    for (unsigned int i = 1; i < coarsen_list.size(); ++i)
      Assert(coarsen_list[i].p.which_tree >= coarsen_list[i - 1].p.which_tree,
             ExcInternalError());

    current_refine_pointer  = refine_list.begin();
    current_coarsen_pointer = coarsen_list.begin();
  }



  template <int dim, int spacedim>
  void
  RefineAndCoarsenList<dim, spacedim>::build_lists(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const typename internal::p4est::types<dim>::quadrant &      p4est_cell,
    const types::subdomain_id                                   my_subdomain)
  {
    if (cell->is_active())
      {
        if (cell->subdomain_id() == my_subdomain)
          {
            if (cell->refine_flag_set())
              refine_list.push_back(p4est_cell);
            else if (cell->coarsen_flag_set())
              coarsen_list.push_back(p4est_cell);
          }
      }
    else
      {
        typename internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          switch (dim)
            {
              case 2:
                P4EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              case 3:
                P8EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }
        internal::p4est::functions<dim>::quadrant_childrenv(&p4est_cell,
                                                            p4est_child);
        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          {
            p4est_child[c].p.which_tree = p4est_cell.p.which_tree;
            build_lists(cell->child(c), p4est_child[c], my_subdomain);
          }
      }
  }


  template <int dim, int spacedim>
  int
  RefineAndCoarsenList<dim, spacedim>::refine_callback(
    typename internal::p4est::types<dim>::forest *  forest,
    typename internal::p4est::types<dim>::topidx    coarse_cell_index,
    typename internal::p4est::types<dim>::quadrant *quadrant)
  {
    RefineAndCoarsenList<dim, spacedim> *this_object =
      reinterpret_cast<RefineAndCoarsenList<dim, spacedim> *>(
        forest->user_pointer);

    // if there are no more cells in our list the current cell can't be
    // flagged for refinement
    if (this_object->current_refine_pointer == this_object->refine_list.end())
      return false;

    Assert(coarse_cell_index <=
             this_object->current_refine_pointer->p.which_tree,
           ExcInternalError());

    // if p4est hasn't yet reached the tree of the next flagged cell the
    // current cell can't be flagged for refinement
    if (coarse_cell_index < this_object->current_refine_pointer->p.which_tree)
      return false;

    // now we're in the right tree in the forest
    Assert(coarse_cell_index <=
             this_object->current_refine_pointer->p.which_tree,
           ExcInternalError());

    // make sure that the p4est loop over cells hasn't gotten ahead of our own
    // pointer
    Assert(internal::p4est::functions<dim>::quadrant_compare(
             quadrant, &*this_object->current_refine_pointer) <= 0,
           ExcInternalError());

    // now, if the p4est cell is one in the list, it is supposed to be refined
    if (internal::p4est::functions<dim>::quadrant_is_equal(
          quadrant, &*this_object->current_refine_pointer))
      {
        ++this_object->current_refine_pointer;
        return true;
      }

    // p4est cell is not in list
    return false;
  }



  template <int dim, int spacedim>
  int
  RefineAndCoarsenList<dim, spacedim>::coarsen_callback(
    typename internal::p4est::types<dim>::forest *  forest,
    typename internal::p4est::types<dim>::topidx    coarse_cell_index,
    typename internal::p4est::types<dim>::quadrant *children[])
  {
    RefineAndCoarsenList<dim, spacedim> *this_object =
      reinterpret_cast<RefineAndCoarsenList<dim, spacedim> *>(
        forest->user_pointer);

    // if there are no more cells in our list the current cell can't be
    // flagged for coarsening
    if (this_object->current_coarsen_pointer == this_object->coarsen_list.end())
      return false;

    Assert(coarse_cell_index <=
             this_object->current_coarsen_pointer->p.which_tree,
           ExcInternalError());

    // if p4est hasn't yet reached the tree of the next flagged cell the
    // current cell can't be flagged for coarsening
    if (coarse_cell_index < this_object->current_coarsen_pointer->p.which_tree)
      return false;

    // now we're in the right tree in the forest
    Assert(coarse_cell_index <=
             this_object->current_coarsen_pointer->p.which_tree,
           ExcInternalError());

    // make sure that the p4est loop over cells hasn't gotten ahead of our own
    // pointer
    Assert(internal::p4est::functions<dim>::quadrant_compare(
             children[0], &*this_object->current_coarsen_pointer) <= 0,
           ExcInternalError());

    // now, if the p4est cell is one in the list, it is supposed to be
    // coarsened
    if (internal::p4est::functions<dim>::quadrant_is_equal(
          children[0], &*this_object->current_coarsen_pointer))
      {
        // move current pointer one up
        ++this_object->current_coarsen_pointer;

        // note that the next 3 cells in our list need to correspond to the
        // other siblings of the cell we have just found
        for (unsigned int c = 1; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          {
            Assert(internal::p4est::functions<dim>::quadrant_is_equal(
                     children[c], &*this_object->current_coarsen_pointer),
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
    explicit PartitionWeights(const std::vector<unsigned int> &cell_weights);

    /**
     * A callback function that we pass to the p4est data structures when a
     * forest is to be partitioned. The p4est functions call it back with a tree
     * (the index of the tree that grows out of a given coarse cell) and a
     * refinement path from that coarse cell to a terminal/leaf cell. The
     * function returns the weight of the cell.
     */
    static int
    cell_weight(typename internal::p4est::types<dim>::forest *forest,
                typename internal::p4est::types<dim>::topidx  coarse_cell_index,
                typename internal::p4est::types<dim>::quadrant *quadrant);

  private:
    std::vector<unsigned int>                 cell_weights_list;
    std::vector<unsigned int>::const_iterator current_pointer;
  };


  template <int dim, int spacedim>
  PartitionWeights<dim, spacedim>::PartitionWeights(
    const std::vector<unsigned int> &cell_weights)
    : cell_weights_list(cell_weights)
  {
    // set the current pointer to the first element of the list, given that
    // we will walk through it sequentially
    current_pointer = cell_weights_list.begin();
  }


  template <int dim, int spacedim>
  int
  PartitionWeights<dim, spacedim>::cell_weight(
    typename internal::p4est::types<dim>::forest *forest,
    typename internal::p4est::types<dim>::topidx,
    typename internal::p4est::types<dim>::quadrant *)
  {
    // the function gets two additional arguments, but we don't need them
    // since we know in which order p4est will walk through the cells
    // and have already built our weight lists in this order

    PartitionWeights<dim, spacedim> *this_object =
      reinterpret_cast<PartitionWeights<dim, spacedim> *>(forest->user_pointer);

    Assert(this_object->current_pointer >=
             this_object->cell_weights_list.begin(),
           ExcInternalError());
    Assert(this_object->current_pointer < this_object->cell_weights_list.end(),
           ExcInternalError());

    // get the weight, increment the pointer, and return the weight
    return *this_object->current_pointer++;
  }

  template <int dim, int spacedim>
  using cell_relation_t = typename std::pair<
    typename dealii::Triangulation<dim, spacedim>::cell_iterator,
    typename dealii::Triangulation<dim, spacedim>::CellStatus>;

  /**
   * Adds a pair of a @p dealii_cell and its @p status
   * to the vector containing all relations @p cell_rel.
   * The pair will be inserted in the position corresponding to the one
   * of the p4est quadrant in the underlying p4est sc_array. The position
   * will be determined from @p idx, which is the position of the quadrant
   * in its corresponding @p tree. The p4est quadrant will be deduced from
   * the @p tree by @p idx.
   */
  template <int dim, int spacedim>
  inline void
  add_single_cell_relation(
    std::vector<cell_relation_t<dim, spacedim>> &               cell_rel,
    const typename dealii::internal::p4est::types<dim>::tree &  tree,
    const unsigned int                                          idx,
    const typename Triangulation<dim, spacedim>::cell_iterator &dealii_cell,
    const typename Triangulation<dim, spacedim>::CellStatus     status)
  {
    const unsigned int local_quadrant_index = tree.quadrants_offset + idx;

    // check if we will be writing into valid memory
    Assert(local_quadrant_index < cell_rel.size(), ExcInternalError());

    // store relation
    cell_rel[local_quadrant_index] = std::make_pair(dealii_cell, status);
  }



  /**
   * This is the recursive part of the member function
   * update_cell_relations().
   *
   * Find the relation between the @p p4est_cell and the @p dealii_cell in the
   * corresponding @p tree. Depending on the CellStatus relation between the two,
   * a new entry will either be inserted in @p cell_rel or the recursion
   * will be continued.
   */
  template <int dim, int spacedim>
  void
  update_cell_relations_recursively(
    std::vector<cell_relation_t<dim, spacedim>> &                 cell_rel,
    const typename dealii::internal::p4est::types<dim>::tree &    tree,
    const typename Triangulation<dim, spacedim>::cell_iterator &  dealii_cell,
    const typename dealii::internal::p4est::types<dim>::quadrant &p4est_cell)
  {
    // find index of p4est_cell in the quadrants array of the corresponding tree
    const int idx = sc_array_bsearch(
      const_cast<sc_array_t *>(&tree.quadrants),
      &p4est_cell,
      dealii::internal::p4est::functions<dim>::quadrant_compare);
    if (idx == -1 &&
        (dealii::internal::p4est::functions<dim>::quadrant_overlaps_tree(
           const_cast<typename dealii::internal::p4est::types<dim>::tree *>(
             &tree),
           &p4est_cell) == false))
      // this quadrant and none of its children belong to us.
      return;

    // recurse further if both p4est and dealii still have children
    const bool p4est_has_children = (idx == -1);
    if (p4est_has_children && dealii_cell->has_children())
      {
        // recurse further
        typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          switch (dim)
            {
              case 2:
                P4EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              case 3:
                P8EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }

        dealii::internal::p4est::functions<dim>::quadrant_childrenv(
          &p4est_cell, p4est_child);

        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          {
            update_cell_relations_recursively<dim, spacedim>(
              cell_rel, tree, dealii_cell->child(c), p4est_child[c]);
          }
      }
    else if (!p4est_has_children && !dealii_cell->has_children())
      {
        // this active cell didn't change
        // save pair into corresponding position
        add_single_cell_relation<dim, spacedim>(
          cell_rel,
          tree,
          idx,
          dealii_cell,
          Triangulation<dim, spacedim>::CELL_PERSIST);
      }
    else if (p4est_has_children) // based on the conditions above, we know that
                                 // dealii_cell has no children
      {
        // this cell got refined in p4est, but the dealii_cell has not yet been
        // refined

        // this quadrant is not active
        // generate its children, and store information in those
        typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];
        for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
             ++c)
          switch (dim)
            {
              case 2:
                P4EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              case 3:
                P8EST_QUADRANT_INIT(&p4est_child[c]);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }

        dealii::internal::p4est::functions<dim>::quadrant_childrenv(
          &p4est_cell, p4est_child);

        // mark first child with CELL_REFINE and the remaining children with
        // CELL_INVALID, but associate them all with the parent cell unpack
        // algorithm will be called only on CELL_REFINE flagged quadrant
        int                                               child_idx;
        typename Triangulation<dim, spacedim>::CellStatus cell_status;
        for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
             ++i)
          {
            child_idx = sc_array_bsearch(
              const_cast<sc_array_t *>(&tree.quadrants),
              &p4est_child[i],
              dealii::internal::p4est::functions<dim>::quadrant_compare);

            cell_status = (i == 0) ? Triangulation<dim, spacedim>::CELL_REFINE :
                                     Triangulation<dim, spacedim>::CELL_INVALID;

            add_single_cell_relation<dim, spacedim>(
              cell_rel, tree, child_idx, dealii_cell, cell_status);
          }
      }
    else // based on the conditions above, we know that p4est_cell has no
         // children, and the dealii_cell does
      {
        // its children got coarsened into this cell in p4est,
        // but the dealii_cell still has its children
        add_single_cell_relation<dim, spacedim>(
          cell_rel,
          tree,
          idx,
          dealii_cell,
          Triangulation<dim, spacedim>::CELL_COARSEN);
      }
  }
} // namespace



namespace parallel
{
  namespace distributed
  {
    /*----------------- class Triangulation<dim,spacedim> ---------------\*/
    template <int dim, int spacedim>
    Triangulation<dim, spacedim>::Triangulation(
      const MPI_Comm &mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                     smooth_grid,
      const Settings settings)
      : // Do not check for distorted cells.
        // For multigrid, we need limit_level_difference_at_vertices
        // to make sure the transfer operators only need to consider two levels.
      dealii::parallel::DistributedTriangulationBase<dim, spacedim>(
        mpi_communicator,
        (settings & construct_multigrid_hierarchy) ?
          static_cast<
            typename dealii::Triangulation<dim, spacedim>::MeshSmoothing>(
            smooth_grid |
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices) :
          smooth_grid,
        false)
      , settings(settings)
      , triangulation_has_content(false)
      , connectivity(nullptr)
      , parallel_forest(nullptr)
    {
      parallel_ghost = nullptr;
    }



    template <int dim, int spacedim>
    Triangulation<dim, spacedim>::~Triangulation()
    {
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      try
        {
          Triangulation<dim, spacedim>::clear();
        }
      catch (...)
        {}

      AssertNothrow(triangulation_has_content == false, ExcInternalError());
      AssertNothrow(connectivity == nullptr, ExcInternalError());
      AssertNothrow(parallel_forest == nullptr, ExcInternalError());
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::create_triangulation(
      const std::vector<Point<spacedim>> &vertices,
      const std::vector<CellData<dim>> &  cells,
      const SubCellData &                 subcelldata)
    {
      try
        {
          dealii::Triangulation<dim, spacedim>::create_triangulation(
            vertices, cells, subcelldata);
        }
      catch (
        const typename dealii::Triangulation<dim, spacedim>::DistortedCellList
          &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      Assert(
        this->all_reference_cells_are_hyper_cube(),
        ExcMessage(
          "The class parallel::distributed::Triangulation only supports meshes "
          "consisting only of hypercube-like cells."));

      // note that now we have some content in the p4est objects and call the
      // functions that do the actual work (which are dimension dependent, so
      // separate)
      triangulation_has_content = true;

      setup_coarse_cell_to_p4est_tree_permutation();

      copy_new_triangulation_to_p4est(std::integral_constant<int, dim>());

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      this->update_periodic_face_map();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::create_triangulation(
      const TriangulationDescription::Description<dim, spacedim>
        &construction_data)
    {
      (void)construction_data;

      Assert(false, ExcInternalError());
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::clear()
    {
      triangulation_has_content = false;

      if (parallel_ghost != nullptr)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy(
            parallel_ghost);
          parallel_ghost = nullptr;
        }

      if (parallel_forest != nullptr)
        {
          dealii::internal::p4est::functions<dim>::destroy(parallel_forest);
          parallel_forest = nullptr;
        }

      if (connectivity != nullptr)
        {
          dealii::internal::p4est::functions<dim>::connectivity_destroy(
            connectivity);
          connectivity = nullptr;
        }

      coarse_cell_to_p4est_tree_permutation.resize(0);
      p4est_tree_to_coarse_cell_permutation.resize(0);

      dealii::parallel::DistributedTriangulationBase<dim, spacedim>::clear();

      this->update_number_cache();
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      return settings &
             Triangulation<dim, spacedim>::construct_multigrid_hierarchy;
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_transfer(
      const typename dealii::internal::p4est::types<dim>::forest
        *parallel_forest,
      const typename dealii::internal::p4est::types<dim>::gloidx
        *previous_global_first_quadrant)
    {
      Assert(this->data_transfer.sizes_fixed_cumulative.size() > 0,
             ExcMessage("No data has been packed!"));

      // Resize memory according to the data that we will receive.
      this->data_transfer.dest_data_fixed.resize(
        parallel_forest->local_num_quadrants *
        this->data_transfer.sizes_fixed_cumulative.back());

      // Execute non-blocking fixed size transfer.
      typename dealii::internal::p4est::types<dim>::transfer_context
        *tf_context;
      tf_context =
        dealii::internal::p4est::functions<dim>::transfer_fixed_begin(
          parallel_forest->global_first_quadrant,
          previous_global_first_quadrant,
          parallel_forest->mpicomm,
          0,
          this->data_transfer.dest_data_fixed.data(),
          this->data_transfer.src_data_fixed.data(),
          this->data_transfer.sizes_fixed_cumulative.back());

      if (this->data_transfer.variable_size_data_stored)
        {
          // Resize memory according to the data that we will receive.
          this->data_transfer.dest_sizes_variable.resize(
            parallel_forest->local_num_quadrants);

          // Execute fixed size transfer of data sizes for variable size
          // transfer.
          dealii::internal::p4est::functions<dim>::transfer_fixed(
            parallel_forest->global_first_quadrant,
            previous_global_first_quadrant,
            parallel_forest->mpicomm,
            1,
            this->data_transfer.dest_sizes_variable.data(),
            this->data_transfer.src_sizes_variable.data(),
            sizeof(unsigned int));
        }

      dealii::internal::p4est::functions<dim>::transfer_fixed_end(tf_context);

      // Release memory of previously packed data.
      this->data_transfer.src_data_fixed.clear();
      this->data_transfer.src_data_fixed.shrink_to_fit();

      if (this->data_transfer.variable_size_data_stored)
        {
          // Resize memory according to the data that we will receive.
          this->data_transfer.dest_data_variable.resize(
            std::accumulate(this->data_transfer.dest_sizes_variable.begin(),
                            this->data_transfer.dest_sizes_variable.end(),
                            std::vector<int>::size_type(0)));

#  if DEAL_II_P4EST_VERSION_GTE(2, 0, 65, 0)
#  else
          // ----- WORKAROUND -----
          // An assertion in p4est prevents us from sending/receiving no data
          // at all, which is mandatory if one of our processes does not own
          // any quadrant. This bypasses the assertion from being triggered.
          //   - see: https://github.com/cburstedde/p4est/issues/48
          if (this->data_transfer.src_sizes_variable.size() == 0)
            this->data_transfer.src_sizes_variable.resize(1);
          if (this->data_transfer.dest_sizes_variable.size() == 0)
            this->data_transfer.dest_sizes_variable.resize(1);
#  endif

          // Execute variable size transfer.
          dealii::internal::p4est::functions<dim>::transfer_custom(
            parallel_forest->global_first_quadrant,
            previous_global_first_quadrant,
            parallel_forest->mpicomm,
            1,
            this->data_transfer.dest_data_variable.data(),
            this->data_transfer.dest_sizes_variable.data(),
            this->data_transfer.src_data_variable.data(),
            this->data_transfer.src_sizes_variable.data());

          // Release memory of previously packed data.
          this->data_transfer.src_sizes_variable.clear();
          this->data_transfer.src_sizes_variable.shrink_to_fit();
          this->data_transfer.src_data_variable.clear();
          this->data_transfer.src_data_variable.shrink_to_fit();
        }
    }


    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::has_hanging_nodes() const
    {
      if (this->n_global_levels() <= 1)
        return false; // can not have hanging nodes without refined cells

      // if there are any active cells with level less than n_global_levels()-1,
      // then there is obviously also one with level n_global_levels()-1, and
      // consequently there must be a hanging node somewhere.
      //
      // The problem is that we cannot just ask for the first active cell, but
      // instead need to filter over locally owned cells.
      const bool have_coarser_cell =
        std::any_of(this->begin_active(this->n_global_levels() - 2),
                    this->end_active(this->n_global_levels() - 2),
                    [](const CellAccessor<dim, spacedim> &cell) {
                      return cell.is_locally_owned();
                    });

      // return true if at least one process has a coarser cell
      return 0 < Utilities::MPI::max(have_coarser_cell ? 1 : 0,
                                     this->mpi_communicator);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::setup_coarse_cell_to_p4est_tree_permutation()
    {
      DynamicSparsityPattern cell_connectivity;
      dealii::GridTools::get_vertex_connectivity_of_cells(*this,
                                                          cell_connectivity);
      coarse_cell_to_p4est_tree_permutation.resize(this->n_cells(0));
      SparsityTools::reorder_hierarchical(
        cell_connectivity, coarse_cell_to_p4est_tree_permutation);

      p4est_tree_to_coarse_cell_permutation =
        Utilities::invert_permutation(coarse_cell_to_p4est_tree_permutation);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::write_mesh_vtk(
      const std::string &file_basename) const
    {
      Assert(parallel_forest != nullptr,
             ExcMessage("Can't produce output when no forest is created yet."));
      dealii::internal::p4est::functions<dim>::vtk_write_file(
        parallel_forest, nullptr, file_basename.c_str());
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::save(const std::string &filename) const
    {
      Assert(
        this->cell_attached_data.n_attached_deserialize == 0,
        ExcMessage(
          "Not all SolutionTransfer objects have been deserialized after the last call to load()."));
      Assert(this->n_cells() > 0,
             ExcMessage("Can not save() an empty Triangulation."));

      const int myrank =
        Utilities::MPI::this_mpi_process(this->mpi_communicator);

      // signal that serialization is going to happen
      this->signals.pre_distributed_save();

      if (this->my_subdomain == 0)
        {
          std::string   fname = std::string(filename) + ".info";
          std::ofstream f(fname.c_str());
          f << "version nproc n_attached_fixed_size_objs n_attached_variable_size_objs n_coarse_cells"
            << std::endl
            << 4 << " "
            << Utilities::MPI::n_mpi_processes(this->mpi_communicator) << " "
            << this->cell_attached_data.pack_callbacks_fixed.size() << " "
            << this->cell_attached_data.pack_callbacks_variable.size() << " "
            << this->n_cells(0) << std::endl;
        }

      // each cell should have been flagged `CELL_PERSIST`
      for (const auto &cell_rel : this->local_cell_relations)
        {
          (void)cell_rel;
          Assert(
            (cell_rel.second == // cell_status
             parallel::distributed::Triangulation<dim, spacedim>::CELL_PERSIST),
            ExcInternalError());
        }

      // Save cell attached data.
      this->save_attached_data(parallel_forest->global_first_quadrant[myrank],
                               parallel_forest->global_num_quadrants,
                               filename);

      dealii::internal::p4est::functions<dim>::save(filename.c_str(),
                                                    parallel_forest,
                                                    false);

      // signal that serialization has finished
      this->signals.post_distributed_save();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::load(const std::string &filename,
                                       const bool         autopartition)
    {
      Assert(
        this->n_cells() > 0,
        ExcMessage(
          "load() only works if the Triangulation already contains a coarse mesh!"));
      Assert(
        this->n_levels() == 1,
        ExcMessage(
          "Triangulation may only contain coarse cells when calling load()."));

      const int myrank =
        Utilities::MPI::this_mpi_process(this->mpi_communicator);

      // signal that de-serialization is going to happen
      this->signals.pre_distributed_load();

      if (parallel_ghost != nullptr)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy(
            parallel_ghost);
          parallel_ghost = nullptr;
        }
      dealii::internal::p4est::functions<dim>::destroy(parallel_forest);
      parallel_forest = nullptr;
      dealii::internal::p4est::functions<dim>::connectivity_destroy(
        connectivity);
      connectivity = nullptr;

      unsigned int version, numcpus, attached_count_fixed,
        attached_count_variable, n_coarse_cells;
      {
        std::string   fname = std::string(filename) + ".info";
        std::ifstream f(fname.c_str());
        AssertThrow(f, ExcIO());
        std::string firstline;
        getline(f, firstline); // skip first line
        f >> version >> numcpus >> attached_count_fixed >>
          attached_count_variable >> n_coarse_cells;
      }

      AssertThrow(version == 4,
                  ExcMessage("Incompatible version found in .info file."));
      Assert(this->n_cells(0) == n_coarse_cells,
             ExcMessage("Number of coarse cells differ!"));

      // clear all of the callback data, as explained in the documentation of
      // register_data_attach()
      this->cell_attached_data.n_attached_data_sets = 0;
      this->cell_attached_data.n_attached_deserialize =
        attached_count_fixed + attached_count_variable;

      parallel_forest = dealii::internal::p4est::functions<dim>::load_ext(
        filename.c_str(),
        this->mpi_communicator,
        0,
        false,
        autopartition,
        0,
        this,
        &connectivity);

      if (numcpus != Utilities::MPI::n_mpi_processes(this->mpi_communicator))
        {
          // We are changing the number of CPUs so we need to repartition.
          // Note that p4est actually distributes the cells between the changed
          // number of CPUs and so everything works without this call, but
          // this command changes the distribution for some reason, so we
          // will leave it in here.
          if (this->signals.cell_weight.num_slots() == 0)
            {
              // no cell weights given -- call p4est's 'partition' without a
              // callback for cell weights
              dealii::internal::p4est::functions<dim>::partition(
                parallel_forest,
                /* prepare coarsening */ 1,
                /* weight_callback */ nullptr);
            }
          else
            {
              // get cell weights for a weighted repartitioning.
              const std::vector<unsigned int> cell_weights = get_cell_weights();

              PartitionWeights<dim, spacedim> partition_weights(cell_weights);

              // attach (temporarily) a pointer to the cell weights through
              // p4est's user_pointer object
              Assert(parallel_forest->user_pointer == this, ExcInternalError());
              parallel_forest->user_pointer = &partition_weights;

              dealii::internal::p4est::functions<dim>::partition(
                parallel_forest,
                /* prepare coarsening */ 1,
                /* weight_callback */
                &PartitionWeights<dim, spacedim>::cell_weight);

              // reset the user pointer to its previous state
              parallel_forest->user_pointer = this;
            }
        }

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying
          // triangulation should not
          // be checking for
          // distorted cells
          Assert(false, ExcInternalError());
        }

      // Load attached cell data, if any was stored.
      this->load_attached_data(parallel_forest->global_first_quadrant[myrank],
                               parallel_forest->global_num_quadrants,
                               parallel_forest->local_num_quadrants,
                               filename,
                               attached_count_fixed,
                               attached_count_variable);

      // signal that de-serialization is finished
      this->signals.post_distributed_load();

      this->update_periodic_face_map();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::load(
      const typename dealii::internal::p4est::types<dim>::forest *forest)
    {
      Assert(this->n_cells() > 0,
             ExcMessage(
               "load() only works if the Triangulation already contains "
               "a coarse mesh!"));
      Assert(this->n_cells() == forest->trees->elem_count,
             ExcMessage(
               "Coarse mesh of the Triangulation does not match the one "
               "of the provided forest!"));

      // clear the old forest
      if (parallel_ghost != nullptr)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy(
            parallel_ghost);
          parallel_ghost = nullptr;
        }
      dealii::internal::p4est::functions<dim>::destroy(parallel_forest);
      parallel_forest = nullptr;

      // note: we can keep the connectivity, since the coarse grid does not
      // change

      // create deep copy of the new forest
      typename dealii::internal::p4est::types<dim>::forest *temp =
        const_cast<typename dealii::internal::p4est::types<dim>::forest *>(
          forest);
      parallel_forest =
        dealii::internal::p4est::functions<dim>::copy_forest(temp, false);
      parallel_forest->connectivity = connectivity;
      parallel_forest->user_pointer = this;

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      this->update_periodic_face_map();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim, spacedim>::get_checksum() const
    {
      Assert(parallel_forest != nullptr,
             ExcMessage(
               "Can't produce a check sum when no forest is created yet."));
      return dealii::internal::p4est::functions<dim>::checksum(parallel_forest);
    }



    template <int dim, int spacedim>
    const typename dealii::internal::p4est::types<dim>::forest *
    Triangulation<dim, spacedim>::get_p4est() const
    {
      Assert(parallel_forest != nullptr,
             ExcMessage("The forest has not been allocated yet."));
      return parallel_forest;
    }



    template <int dim, int spacedim>
    typename dealii::internal::p4est::types<dim>::tree *
    Triangulation<dim, spacedim>::init_tree(
      const int dealii_coarse_cell_index) const
    {
      const unsigned int tree_index =
        coarse_cell_to_p4est_tree_permutation[dealii_coarse_cell_index];
      typename dealii::internal::p4est::types<dim>::tree *tree =
        static_cast<typename dealii::internal::p4est::types<dim>::tree *>(
          sc_array_index(parallel_forest->trees, tree_index));

      return tree;
    }



    // Note: this has been added here to prevent that these functions
    // appear in the Doxygen documentation of dealii::Triangulation
#  ifndef DOXYGEN

    template <>
    void
    Triangulation<2, 2>::copy_new_triangulation_to_p4est(
      std::integral_constant<int, 2>)
    {
      const unsigned int dim = 2, spacedim = 2;
      Assert(this->n_cells(0) > 0, ExcInternalError());
      Assert(this->n_levels() == 1, ExcInternalError());

      // data structures that counts how many cells touch each vertex
      // (vertex_touch_count), and which cells touch a given vertex (together
      // with the local numbering of that vertex within the cells that touch
      // it)
      std::vector<unsigned int> vertex_touch_count;
      std::vector<
        std::list<std::pair<Triangulation<dim, spacedim>::active_cell_iterator,
                            unsigned int>>>
        vertex_to_cell;
      get_vertex_to_cell_mappings(*this, vertex_touch_count, vertex_to_cell);
      const dealii::internal::p4est::types<2>::locidx num_vtt =
        std::accumulate(vertex_touch_count.begin(),
                        vertex_touch_count.end(),
                        0u);

      // now create a connectivity object with the right sizes for all
      // arrays. set vertex information only in debug mode (saves a few bytes
      // in optimized mode)
      const bool set_vertex_info
#    ifdef DEBUG
        = true
#    else
        = false
#    endif
        ;

      connectivity = dealii::internal::p4est::functions<2>::connectivity_new(
        (set_vertex_info == true ? this->n_vertices() : 0),
        this->n_cells(0),
        this->n_vertices(),
        num_vtt);

      set_vertex_and_cell_info(*this,
                               vertex_touch_count,
                               vertex_to_cell,
                               coarse_cell_to_p4est_tree_permutation,
                               set_vertex_info,
                               connectivity);

      Assert(p4est_connectivity_is_valid(connectivity) == 1,
             ExcInternalError());

      // now create a forest out of the connectivity data structure
      parallel_forest = dealii::internal::p4est::functions<2>::new_forest(
        this->mpi_communicator,
        connectivity,
        /* minimum initial number of quadrants per tree */ 0,
        /* minimum level of upfront refinement */ 0,
        /* use uniform upfront refinement */ 1,
        /* user_data_size = */ 0,
        /* user_data_constructor = */ nullptr,
        /* user_pointer */ this);
    }



    // TODO: This is a verbatim copy of the 2,2 case. However, we can't just
    // specialize the dim template argument, but let spacedim open
    template <>
    void
    Triangulation<2, 3>::copy_new_triangulation_to_p4est(
      std::integral_constant<int, 2>)
    {
      const unsigned int dim = 2, spacedim = 3;
      Assert(this->n_cells(0) > 0, ExcInternalError());
      Assert(this->n_levels() == 1, ExcInternalError());

      // data structures that counts how many cells touch each vertex
      // (vertex_touch_count), and which cells touch a given vertex (together
      // with the local numbering of that vertex within the cells that touch
      // it)
      std::vector<unsigned int> vertex_touch_count;
      std::vector<
        std::list<std::pair<Triangulation<dim, spacedim>::active_cell_iterator,
                            unsigned int>>>
        vertex_to_cell;
      get_vertex_to_cell_mappings(*this, vertex_touch_count, vertex_to_cell);
      const dealii::internal::p4est::types<2>::locidx num_vtt =
        std::accumulate(vertex_touch_count.begin(),
                        vertex_touch_count.end(),
                        0u);

      // now create a connectivity object with the right sizes for all
      // arrays. set vertex information only in debug mode (saves a few bytes
      // in optimized mode)
      const bool set_vertex_info
#    ifdef DEBUG
        = true
#    else
        = false
#    endif
        ;

      connectivity = dealii::internal::p4est::functions<2>::connectivity_new(
        (set_vertex_info == true ? this->n_vertices() : 0),
        this->n_cells(0),
        this->n_vertices(),
        num_vtt);

      set_vertex_and_cell_info(*this,
                               vertex_touch_count,
                               vertex_to_cell,
                               coarse_cell_to_p4est_tree_permutation,
                               set_vertex_info,
                               connectivity);

      Assert(p4est_connectivity_is_valid(connectivity) == 1,
             ExcInternalError());

      // now create a forest out of the connectivity data structure
      parallel_forest = dealii::internal::p4est::functions<2>::new_forest(
        this->mpi_communicator,
        connectivity,
        /* minimum initial number of quadrants per tree */ 0,
        /* minimum level of upfront refinement */ 0,
        /* use uniform upfront refinement */ 1,
        /* user_data_size = */ 0,
        /* user_data_constructor = */ nullptr,
        /* user_pointer */ this);
    }



    template <>
    void
    Triangulation<3, 3>::copy_new_triangulation_to_p4est(
      std::integral_constant<int, 3>)
    {
      const int dim = 3, spacedim = 3;
      Assert(this->n_cells(0) > 0, ExcInternalError());
      Assert(this->n_levels() == 1, ExcInternalError());

      // data structures that counts how many cells touch each vertex
      // (vertex_touch_count), and which cells touch a given vertex (together
      // with the local numbering of that vertex within the cells that touch
      // it)
      std::vector<unsigned int> vertex_touch_count;
      std::vector<std::list<
        std::pair<Triangulation<3>::active_cell_iterator, unsigned int>>>
        vertex_to_cell;
      get_vertex_to_cell_mappings(*this, vertex_touch_count, vertex_to_cell);
      const dealii::internal::p4est::types<2>::locidx num_vtt =
        std::accumulate(vertex_touch_count.begin(),
                        vertex_touch_count.end(),
                        0u);

      std::vector<unsigned int> edge_touch_count;
      std::vector<std::list<
        std::pair<Triangulation<3>::active_cell_iterator, unsigned int>>>
        edge_to_cell;
      get_edge_to_cell_mappings(*this, edge_touch_count, edge_to_cell);
      const dealii::internal::p4est::types<2>::locidx num_ett =
        std::accumulate(edge_touch_count.begin(), edge_touch_count.end(), 0u);

      // now create a connectivity object with the right sizes for all arrays
      const bool set_vertex_info
#    ifdef DEBUG
        = true
#    else
        = false
#    endif
        ;

      connectivity = dealii::internal::p4est::functions<3>::connectivity_new(
        (set_vertex_info == true ? this->n_vertices() : 0),
        this->n_cells(0),
        this->n_active_lines(),
        num_ett,
        this->n_vertices(),
        num_vtt);

      set_vertex_and_cell_info(*this,
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

      const unsigned int deal_to_p4est_line_index[12] = {
        4, 5, 0, 1, 6, 7, 2, 3, 8, 9, 10, 11};

      for (Triangulation<dim, spacedim>::active_cell_iterator cell =
             this->begin_active();
           cell != this->end();
           ++cell)
        {
          const unsigned int index =
            coarse_cell_to_p4est_tree_permutation[cell->index()];
          for (unsigned int e = 0; e < GeometryInfo<3>::lines_per_cell; ++e)
            connectivity->tree_to_edge[index * GeometryInfo<3>::lines_per_cell +
                                       deal_to_p4est_line_index[e]] =
              cell->line(e)->index();
        }

      // now also set edge-to-tree
      // information
      connectivity->ett_offset[0] = 0;
      std::partial_sum(edge_touch_count.begin(),
                       edge_touch_count.end(),
                       &connectivity->ett_offset[1]);

      Assert(connectivity->ett_offset[this->n_active_lines()] == num_ett,
             ExcInternalError());

      for (unsigned int v = 0; v < this->n_active_lines(); ++v)
        {
          Assert(edge_to_cell[v].size() == edge_touch_count[v],
                 ExcInternalError());

          std::list<
            std::pair<Triangulation<dim, spacedim>::active_cell_iterator,
                      unsigned int>>::const_iterator p =
            edge_to_cell[v].begin();
          for (unsigned int c = 0; c < edge_touch_count[v]; ++c, ++p)
            {
              connectivity->edge_to_tree[connectivity->ett_offset[v] + c] =
                coarse_cell_to_p4est_tree_permutation[p->first->index()];
              connectivity->edge_to_edge[connectivity->ett_offset[v] + c] =
                deal_to_p4est_line_index[p->second];
            }
        }

      Assert(p8est_connectivity_is_valid(connectivity) == 1,
             ExcInternalError());

      // now create a forest out of the connectivity data structure
      parallel_forest = dealii::internal::p4est::functions<3>::new_forest(
        this->mpi_communicator,
        connectivity,
        /* minimum initial number of quadrants per tree */ 0,
        /* minimum level of upfront refinement */ 0,
        /* use uniform upfront refinement */ 1,
        /* user_data_size = */ 0,
        /* user_data_constructor = */ nullptr,
        /* user_pointer */ this);
    }
#  endif



    namespace
    {
      // ensures the 2:1 mesh balance for periodic boundary conditions in the
      // artificial cell layer (the active cells are taken care of by p4est)
      template <int dim, int spacedim>
      bool
      enforce_mesh_balance_over_periodic_boundaries(
        Triangulation<dim, spacedim> &tria)
      {
        if (tria.get_periodic_face_map().size() == 0)
          return false;

        std::vector<bool> flags_before[2];
        tria.save_coarsen_flags(flags_before[0]);
        tria.save_refine_flags(flags_before[1]);

        std::vector<unsigned int> topological_vertex_numbering(
          tria.n_vertices());
        for (unsigned int i = 0; i < topological_vertex_numbering.size(); ++i)
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
        using cell_iterator =
          typename Triangulation<dim, spacedim>::cell_iterator;
        typename std::map<std::pair<cell_iterator, unsigned int>,
                          std::pair<std::pair<cell_iterator, unsigned int>,
                                    std::bitset<3>>>::const_iterator it;
        for (it = tria.get_periodic_face_map().begin();
             it != tria.get_periodic_face_map().end();
             ++it)
          {
            const cell_iterator &cell_1           = it->first.first;
            const unsigned int   face_no_1        = it->first.second;
            const cell_iterator &cell_2           = it->second.first.first;
            const unsigned int   face_no_2        = it->second.first.second;
            const std::bitset<3> face_orientation = it->second.second;

            if (cell_1->level() == cell_2->level())
              {
                for (unsigned int v = 0;
                     v < GeometryInfo<dim - 1>::vertices_per_cell;
                     ++v)
                  {
                    // take possible non-standard orientation of face on cell[0]
                    // into account
                    const unsigned int vface0 =
                      GeometryInfo<dim>::standard_to_real_face_vertex(
                        v,
                        face_orientation[0],
                        face_orientation[1],
                        face_orientation[2]);
                    const unsigned int vi0 =
                      topological_vertex_numbering[cell_1->face(face_no_1)
                                                     ->vertex_index(vface0)];
                    const unsigned int vi1 =
                      topological_vertex_numbering[cell_2->face(face_no_2)
                                                     ->vertex_index(v)];
                    const unsigned int min_index = std::min(vi0, vi1);
                    topological_vertex_numbering[cell_1->face(face_no_1)
                                                   ->vertex_index(vface0)] =
                      topological_vertex_numbering[cell_2->face(face_no_2)
                                                     ->vertex_index(v)] =
                        min_index;
                  }
              }
          }

        // There must not be any chains!
        for (unsigned int i = 0; i < topological_vertex_numbering.size(); ++i)
          {
            const unsigned int j = topological_vertex_numbering[i];
            if (j != i)
              Assert(topological_vertex_numbering[j] == j, ExcInternalError());
          }


        // this code is replicated from grid/tria.cc but using an indirection
        // for periodic boundary conditions
        bool             continue_iterating = true;
        std::vector<int> vertex_level(tria.n_vertices());
        while (continue_iterating)
          {
            // store highest level one of the cells adjacent to a vertex
            // belongs to
            std::fill(vertex_level.begin(), vertex_level.end(), 0);
            typename Triangulation<dim, spacedim>::active_cell_iterator
              cell = tria.begin_active(),
              endc = tria.end();
            for (; cell != endc; ++cell)
              {
                if (cell->refine_flag_set())
                  for (const unsigned int vertex :
                       GeometryInfo<dim>::vertex_indices())
                    vertex_level[topological_vertex_numbering
                                   [cell->vertex_index(vertex)]] =
                      std::max(vertex_level[topological_vertex_numbering
                                              [cell->vertex_index(vertex)]],
                               cell->level() + 1);
                else if (!cell->coarsen_flag_set())
                  for (const unsigned int vertex :
                       GeometryInfo<dim>::vertex_indices())
                    vertex_level[topological_vertex_numbering
                                   [cell->vertex_index(vertex)]] =
                      std::max(vertex_level[topological_vertex_numbering
                                              [cell->vertex_index(vertex)]],
                               cell->level());
                else
                  {
                    // if coarsen flag is set then tentatively assume
                    // that the cell will be coarsened. this isn't
                    // always true (the coarsen flag could be removed
                    // again) and so we may make an error here. we try
                    // to correct this by iterating over the entire
                    // process until we are converged
                    Assert(cell->coarsen_flag_set(), ExcInternalError());
                    for (const unsigned int vertex :
                         GeometryInfo<dim>::vertex_indices())
                      vertex_level[topological_vertex_numbering
                                     [cell->vertex_index(vertex)]] =
                        std::max(vertex_level[topological_vertex_numbering
                                                [cell->vertex_index(vertex)]],
                                 cell->level() - 1);
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
            for (cell = tria.last_active(); cell != endc; --cell)
              if (cell->refine_flag_set() == false)
                {
                  for (const unsigned int vertex :
                       GeometryInfo<dim>::vertex_indices())
                    if (vertex_level[topological_vertex_numbering
                                       [cell->vertex_index(vertex)]] >=
                        cell->level() + 1)
                      {
                        // remove coarsen flag...
                        cell->clear_coarsen_flag();

                        // ...and if necessary also refine the current
                        // cell, at the same time updating the level
                        // information about vertices
                        if (vertex_level[topological_vertex_numbering
                                           [cell->vertex_index(vertex)]] >
                            cell->level() + 1)
                          {
                            cell->set_refine_flag();
                            continue_iterating = true;

                            for (const unsigned int v :
                                 GeometryInfo<dim>::vertex_indices())
                              vertex_level[topological_vertex_numbering
                                             [cell->vertex_index(v)]] =
                                std::max(
                                  vertex_level[topological_vertex_numbering
                                                 [cell->vertex_index(v)]],
                                  cell->level() + 1);
                          }

                        // continue and see whether we may, for example,
                        // go into the inner 'if' above based on a
                        // different vertex
                      }
                }

            // clear coarsen flag if not all children were marked
            for (const auto &cell : tria.cell_iterators())
              {
                // nothing to do if we are already on the finest level
                if (cell->is_active())
                  continue;

                const unsigned int n_children       = cell->n_children();
                unsigned int       flagged_children = 0;
                for (unsigned int child = 0; child < n_children; ++child)
                  if (cell->child(child)->is_active() &&
                      cell->child(child)->coarsen_flag_set())
                    ++flagged_children;

                // if not all children were flagged for coarsening, remove
                // coarsen flags
                if (flagged_children < n_children)
                  for (unsigned int child = 0; child < n_children; ++child)
                    if (cell->child(child)->is_active())
                      cell->child(child)->clear_coarsen_flag();
              }
          }
        std::vector<bool> flags_after[2];
        tria.save_coarsen_flags(flags_after[0]);
        tria.save_refine_flags(flags_after[1]);
        return ((flags_before[0] != flags_after[0]) ||
                (flags_before[1] != flags_after[1]));
      }
    } // namespace



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::prepare_coarsening_and_refinement()
    {
      std::vector<bool> flags_before[2];
      this->save_coarsen_flags(flags_before[0]);
      this->save_refine_flags(flags_before[1]);

      bool         mesh_changed = false;
      unsigned int loop_counter = 0;
      do
        {
          this->dealii::Triangulation<dim, spacedim>::
            prepare_coarsening_and_refinement();
          this->update_periodic_face_map();
          // enforce 2:1 mesh balance over periodic boundaries
          mesh_changed = enforce_mesh_balance_over_periodic_boundaries(*this);

          // We can't be sure that we won't run into a situation where we can
          // not reconcile mesh smoothing and balancing of periodic faces. As we
          // don't know what else to do, at least abort with an error message.
          ++loop_counter;
          AssertThrow(
            loop_counter < 32,
            ExcMessage(
              "Infinite loop in "
              "parallel::distributed::Triangulation::prepare_coarsening_and_refinement() "
              "for periodic boundaries detected. Aborting."));
        }
      while (mesh_changed);

      // check if any of the refinement flags were changed during this
      // function and return that value
      std::vector<bool> flags_after[2];
      this->save_coarsen_flags(flags_after[0]);
      this->save_refine_flags(flags_after[1]);
      return ((flags_before[0] != flags_after[0]) ||
              (flags_before[1] != flags_after[1]));
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::copy_local_forest_to_triangulation()
    {
      // Disable mesh smoothing for recreating the deal.II triangulation,
      // otherwise we might not be able to reproduce the p4est mesh
      // exactly. We restore the original smoothing at the end of this
      // function. Note that the smoothing flag is used in the normal
      // refinement process.
      typename Triangulation<dim, spacedim>::MeshSmoothing save_smooth =
        this->smooth_grid;

      // We will refine manually to match the p4est further down, which
      // obeys a level difference of 2 at each vertex (see the balance call
      // to p4est). We can disable this here so we store fewer artificial
      // cells (in some cases).
      // For geometric multigrid it turns out that
      // we will miss level cells at shared vertices if we ignore this.
      // See tests/mpi/mg_06. In particular, the flag is still necessary
      // even though we force it for the original smooth_grid in the
      // constructor.
      if (settings & construct_multigrid_hierarchy)
        this->smooth_grid =
          dealii::Triangulation<dim,
                                spacedim>::limit_level_difference_at_vertices;
      else
        this->smooth_grid = dealii::Triangulation<dim, spacedim>::none;

      bool mesh_changed = false;

      // Remove all deal.II refinements. Note that we could skip this and
      // start from our current state, because the algorithm later coarsens as
      // necessary. This has the advantage of being faster when large parts
      // of the local partition changes (likely) and gives a deterministic
      // ordering of the cells (useful for snapshot/resume).
      // TODO: is there a more efficient way to do this?
      if (settings & mesh_reconstruction_after_repartitioning)
        while (this->n_levels() > 1)
          {
            // Instead of marking all active cells, we slice off the finest
            // level, one level at a time. This takes the same number of
            // iterations but solves an issue where not all cells on a
            // periodic boundary are indeed coarsened and we run into an
            // irrelevant Assert() in update_periodic_face_map().
            for (const auto &cell :
                 this->active_cell_iterators_on_level(this->n_levels() - 1))
              {
                cell->set_coarsen_flag();
              }
            try
              {
                dealii::Triangulation<dim, spacedim>::
                  execute_coarsening_and_refinement();
              }
            catch (
              const typename Triangulation<dim, spacedim>::DistortedCellList &)
              {
                // the underlying triangulation should not be checking for
                // distorted cells
                Assert(false, ExcInternalError());
              }
          }


      // query p4est for the ghost cells
      if (parallel_ghost != nullptr)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy(
            parallel_ghost);
          parallel_ghost = nullptr;
        }
      parallel_ghost = dealii::internal::p4est::functions<dim>::ghost_new(
        parallel_forest,
        (dim == 2 ? typename dealii::internal::p4est::types<dim>::balance_type(
                      P4EST_CONNECT_CORNER) :
                    typename dealii::internal::p4est::types<dim>::balance_type(
                      P8EST_CONNECT_CORNER)));

      Assert(parallel_ghost, ExcInternalError());


      // set all cells to artificial. we will later set it to the correct
      // subdomain in match_tree_recursively
      for (const auto &cell : this->cell_iterators_on_level(0))
        cell->recursively_set_subdomain_id(numbers::artificial_subdomain_id);

      do
        {
          for (const auto &cell : this->cell_iterators_on_level(0))
            {
              // if this processor stores no part of the forest that comes out
              // of this coarse grid cell, then we need to delete all children
              // of this cell (the coarse grid cell remains)
              if (tree_exists_locally<dim, spacedim>(
                    parallel_forest,
                    coarse_cell_to_p4est_tree_permutation[cell->index()]) ==
                  false)
                {
                  delete_all_children<dim, spacedim>(cell);
                  if (cell->is_active())
                    cell->set_subdomain_id(numbers::artificial_subdomain_id);
                }

              else
                {
                  // this processor stores at least a part of the tree that
                  // comes out of this cell.

                  typename dealii::internal::p4est::types<dim>::quadrant
                                                                      p4est_coarse_cell;
                  typename dealii::internal::p4est::types<dim>::tree *tree =
                    init_tree(cell->index());

                  dealii::internal::p4est::init_coarse_quadrant<dim>(
                    p4est_coarse_cell);

                  match_tree_recursively<dim, spacedim>(*tree,
                                                        cell,
                                                        p4est_coarse_cell,
                                                        *parallel_forest,
                                                        this->my_subdomain);
                }
            }

          // check mesh for ghost cells, refine as necessary. iterate over
          // every ghostquadrant, find corresponding deal coarsecell and
          // recurse.
          typename dealii::internal::p4est::types<dim>::quadrant *quadr;
          types::subdomain_id                                  ghost_owner = 0;
          typename dealii::internal::p4est::types<dim>::topidx ghost_tree  = 0;

          for (unsigned int g_idx = 0;
               g_idx < parallel_ghost->ghosts.elem_count;
               ++g_idx)
            {
              while (g_idx >= static_cast<unsigned int>(
                                parallel_ghost->proc_offsets[ghost_owner + 1]))
                ++ghost_owner;
              while (g_idx >= static_cast<unsigned int>(
                                parallel_ghost->tree_offsets[ghost_tree + 1]))
                ++ghost_tree;

              quadr = static_cast<
                typename dealii::internal::p4est::types<dim>::quadrant *>(
                sc_array_index(&parallel_ghost->ghosts, g_idx));

              unsigned int coarse_cell_index =
                p4est_tree_to_coarse_cell_permutation[ghost_tree];

              match_quadrant<dim, spacedim>(this,
                                            coarse_cell_index,
                                            *quadr,
                                            ghost_owner);
            }

          // fix all the flags to make sure we have a consistent mesh
          this->prepare_coarsening_and_refinement();

          // see if any flags are still set
          mesh_changed =
            std::any_of(this->begin_active(),
                        active_cell_iterator{this->end()},
                        [](const CellAccessor<dim, spacedim> &cell) {
                          return cell.refine_flag_set() ||
                                 cell.coarsen_flag_set();
                        });

          // actually do the refinement to change the local mesh by
          // calling the base class refinement function directly
          try
            {
              dealii::Triangulation<dim, spacedim>::
                execute_coarsening_and_refinement();
            }
          catch (
            const typename Triangulation<dim, spacedim>::DistortedCellList &)
            {
              // the underlying triangulation should not be checking for
              // distorted cells
              Assert(false, ExcInternalError());
            }
        }
      while (mesh_changed);

#  ifdef DEBUG
      // check if correct number of ghosts is created
      unsigned int num_ghosts = 0;

      for (const auto &cell : this->active_cell_iterators())
        {
          if (cell->subdomain_id() != this->my_subdomain &&
              cell->subdomain_id() != numbers::artificial_subdomain_id)
            ++num_ghosts;
        }

      Assert(num_ghosts == parallel_ghost->ghosts.elem_count,
             ExcInternalError());
#  endif



      // fill level_subdomain_ids for geometric multigrid
      // the level ownership of a cell is defined as the owner if the cell is
      // active or as the owner of child(0) we need this information for all our
      // ancestors and the same-level neighbors of our own cells (=level ghosts)
      if (settings & construct_multigrid_hierarchy)
        {
          // step 1: We set our own ids all the way down and all the others to
          // -1. Note that we do not fill other cells we could figure out the
          // same way, because we might accidentally set an id for a cell that
          // is not a ghost cell.
          for (unsigned int lvl = this->n_levels(); lvl > 0;)
            {
              --lvl;
              for (const auto &cell : this->cell_iterators_on_level(lvl))
                {
                  if ((cell->is_active() &&
                       cell->subdomain_id() ==
                         this->locally_owned_subdomain()) ||
                      (cell->has_children() &&
                       cell->child(0)->level_subdomain_id() ==
                         this->locally_owned_subdomain()))
                    cell->set_level_subdomain_id(
                      this->locally_owned_subdomain());
                  else
                    {
                      // not our cell
                      cell->set_level_subdomain_id(
                        numbers::artificial_subdomain_id);
                    }
                }
            }

          // step 2: make sure all the neighbors to our level_cells exist. Need
          // to look up in p4est...
          std::vector<std::vector<bool>> marked_vertices(this->n_levels());
          for (unsigned int lvl = 0; lvl < this->n_levels(); ++lvl)
            marked_vertices[lvl] = mark_locally_active_vertices_on_level(lvl);

          for (const auto &cell : this->cell_iterators_on_level(0))
            {
              typename dealii::internal::p4est::types<dim>::quadrant
                                 p4est_coarse_cell;
              const unsigned int tree_index =
                coarse_cell_to_p4est_tree_permutation[cell->index()];
              typename dealii::internal::p4est::types<dim>::tree *tree =
                init_tree(cell->index());

              dealii::internal::p4est::init_coarse_quadrant<dim>(
                p4est_coarse_cell);

              determine_level_subdomain_id_recursively<dim, spacedim>(
                *tree,
                tree_index,
                cell,
                p4est_coarse_cell,
                *parallel_forest,
                this->my_subdomain,
                marked_vertices);
            }

          // step 3: make sure we have the parent of our level cells
          for (unsigned int lvl = this->n_levels(); lvl > 0;)
            {
              --lvl;
              for (const auto &cell : this->cell_iterators_on_level(lvl))
                {
                  if (cell->has_children())
                    for (unsigned int c = 0;
                         c < GeometryInfo<dim>::max_children_per_cell;
                         ++c)
                      {
                        if (cell->child(c)->level_subdomain_id() ==
                            this->locally_owned_subdomain())
                          {
                            // at least one of the children belongs to us, so
                            // make sure we set the level subdomain id
                            const types::subdomain_id mark =
                              cell->child(0)->level_subdomain_id();
                            Assert(mark != numbers::artificial_subdomain_id,
                                   ExcInternalError()); // we should know the
                                                        // child(0)
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

        if (Utilities::MPI::n_mpi_processes(this->mpi_communicator) == 1)
          {
            Assert(static_cast<unsigned int>(
                     parallel_forest->local_num_quadrants) == total_local_cells,
                   ExcInternalError());
          }
        else
          {
            Assert(static_cast<unsigned int>(
                     parallel_forest->local_num_quadrants) <= total_local_cells,
                   ExcInternalError());
          }

        // count the number of owned, active cells and compare with p4est.
        unsigned int n_owned = 0;
        for (const auto &cell : this->active_cell_iterators())
          {
            if (cell->subdomain_id() == this->my_subdomain)
              ++n_owned;
          }

        Assert(static_cast<unsigned int>(
                 parallel_forest->local_num_quadrants) == n_owned,
               ExcInternalError());
      }

      this->smooth_grid = save_smooth;

      // finally, after syncing the parallel_forest with the triangulation,
      // also update the cell_relations, which will be used for
      // repartitioning, further refinement/coarsening, and unpacking
      // of stored or transferred data.
      update_cell_relations();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_coarsening_and_refinement()
    {
      // do not allow anisotropic refinement
#  ifdef DEBUG
      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_locally_owned() && cell->refine_flag_set())
          Assert(cell->refine_flag_set() ==
                   RefinementPossibilities<dim>::isotropic_refinement,
                 ExcMessage(
                   "This class does not support anisotropic refinement"));
#  endif


      // safety check: p4est has an upper limit on the level of a cell
      if (this->n_levels() ==
          dealii::internal::p4est::functions<dim>::max_level)
        {
          for (typename Triangulation<dim, spacedim>::active_cell_iterator
                 cell = this->begin_active(
                   dealii::internal::p4est::functions<dim>::max_level - 1);
               cell !=
               this->end(dealii::internal::p4est::functions<dim>::max_level -
                         1);
               ++cell)
            {
              AssertThrow(
                !(cell->refine_flag_set()),
                ExcMessage(
                  "Fatal Error: maximum refinement level of p4est reached."));
            }
        }

      this->prepare_coarsening_and_refinement();

      // signal that refinement is going to happen
      this->signals.pre_distributed_refinement();

      // now do the work we're supposed to do when we are in charge
      // make sure all flags are cleared on cells we don't own, since nothing
      // good can come of that if they are still around
      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_ghost() || cell->is_artificial())
          {
            cell->clear_refine_flag();
            cell->clear_coarsen_flag();
          }


      // count how many cells will be refined and coarsened, and allocate that
      // much memory
      RefineAndCoarsenList<dim, spacedim> refine_and_coarsen_list(
        *this, p4est_tree_to_coarse_cell_permutation, this->my_subdomain);

      // copy refine and coarsen flags into p4est and execute the refinement
      // and coarsening. this uses the refine_and_coarsen_list just built,
      // which is communicated to the callback functions through
      // p4est's user_pointer object
      Assert(parallel_forest->user_pointer == this, ExcInternalError());
      parallel_forest->user_pointer = &refine_and_coarsen_list;

      if (parallel_ghost != nullptr)
        {
          dealii::internal::p4est::functions<dim>::ghost_destroy(
            parallel_ghost);
          parallel_ghost = nullptr;
        }
      dealii::internal::p4est::functions<dim>::refine(
        parallel_forest,
        /* refine_recursive */ false,
        &RefineAndCoarsenList<dim, spacedim>::refine_callback,
        /*init_callback=*/nullptr);
      dealii::internal::p4est::functions<dim>::coarsen(
        parallel_forest,
        /* coarsen_recursive */ false,
        &RefineAndCoarsenList<dim, spacedim>::coarsen_callback,
        /*init_callback=*/nullptr);

      // make sure all cells in the lists have been consumed
      Assert(refine_and_coarsen_list.pointers_are_at_end(), ExcInternalError());

      // reset the pointer
      parallel_forest->user_pointer = this;

      // enforce 2:1 hanging node condition
      dealii::internal::p4est::functions<dim>::balance(
        parallel_forest,
        /* face and corner balance */
        (dim == 2 ? typename dealii::internal::p4est::types<dim>::balance_type(
                      P4EST_CONNECT_FULL) :
                    typename dealii::internal::p4est::types<dim>::balance_type(
                      P8EST_CONNECT_FULL)),
        /*init_callback=*/nullptr);

      // since refinement and/or coarsening on the parallel forest
      // has happened, we need to update the quadrant cell relations
      update_cell_relations();

      // signals that parallel_forest has been refined and cell relations have
      // been updated
      this->signals.post_p4est_refinement();

      // before repartitioning the mesh, save a copy of the current positions of
      // quadrants
      // only if data needs to be transferred later
      std::vector<typename dealii::internal::p4est::types<dim>::gloidx>
        previous_global_first_quadrant;

      if (this->cell_attached_data.n_attached_data_sets > 0)
        {
          previous_global_first_quadrant.resize(parallel_forest->mpisize + 1);
          std::memcpy(previous_global_first_quadrant.data(),
                      parallel_forest->global_first_quadrant,
                      sizeof(
                        typename dealii::internal::p4est::types<dim>::gloidx) *
                        (parallel_forest->mpisize + 1));
        }

      if (!(settings & no_automatic_repartitioning))
        {
          // partition the new mesh between all processors. If cell weights have
          // not been given balance the number of cells.
          if (this->signals.cell_weight.num_slots() == 0)
            dealii::internal::p4est::functions<dim>::partition(
              parallel_forest,
              /* prepare coarsening */ 1,
              /* weight_callback */ nullptr);
          else
            {
              // get cell weights for a weighted repartitioning.
              const std::vector<unsigned int> cell_weights = get_cell_weights();

              PartitionWeights<dim, spacedim> partition_weights(cell_weights);

              // attach (temporarily) a pointer to the cell weights through
              // p4est's user_pointer object
              Assert(parallel_forest->user_pointer == this, ExcInternalError());
              parallel_forest->user_pointer = &partition_weights;

              dealii::internal::p4est::functions<dim>::partition(
                parallel_forest,
                /* prepare coarsening */ 1,
                /* weight_callback */
                &PartitionWeights<dim, spacedim>::cell_weight);

              // release data
              dealii::internal::p4est::functions<dim>::reset_data(
                parallel_forest, 0, nullptr, nullptr);
              // reset the user pointer to its previous state
              parallel_forest->user_pointer = this;
            }
        }

      // pack data before triangulation gets updated
      if (this->cell_attached_data.n_attached_data_sets > 0)
        {
          this->data_transfer.pack_data(
            this->local_cell_relations,
            this->cell_attached_data.pack_callbacks_fixed,
            this->cell_attached_data.pack_callbacks_variable);
        }

      // finally copy back from local part of tree to deal.II
      // triangulation. before doing so, make sure there are no refine or
      // coarsen flags pending
      for (const auto &cell : this->active_cell_iterators())
        {
          cell->clear_refine_flag();
          cell->clear_coarsen_flag();
        }

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      // transfer data after triangulation got updated
      if (this->cell_attached_data.n_attached_data_sets > 0)
        {
          this->execute_transfer(parallel_forest,
                                 previous_global_first_quadrant.data());

          // also update the CellStatus information on the new mesh
          this->data_transfer.unpack_cell_status(this->local_cell_relations);
        }

#  ifdef DEBUG
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
          for (unsigned int lvl = 0; lvl < this->n_global_levels(); ++lvl)
            {
              std::vector<bool> active_verts =
                this->mark_locally_active_vertices_on_level(lvl);

              const unsigned int maybe_coarser_lvl =
                (lvl > 0) ? (lvl - 1) : lvl;
              typename Triangulation<dim, spacedim>::cell_iterator
                cell = this->begin(maybe_coarser_lvl),
                endc = this->end(lvl);
              for (; cell != endc; ++cell)
                if (cell->level() == static_cast<int>(lvl) || cell->is_active())
                  {
                    const bool is_level_artificial =
                      (cell->level_subdomain_id() ==
                       numbers::artificial_subdomain_id);
                    bool need_to_know = false;
                    for (const unsigned int vertex :
                         GeometryInfo<dim>::vertex_indices())
                      if (active_verts[cell->vertex_index(vertex)])
                        {
                          need_to_know = true;
                          break;
                        }

                    Assert(
                      !need_to_know || !is_level_artificial,
                      ExcMessage(
                        "Internal error: the owner of cell" +
                        cell->id().to_string() +
                        " is unknown even though it is needed for geometric multigrid."));
                  }
            }
        }
#  endif

      this->update_periodic_face_map();
      this->update_number_cache();

      // signal that refinement is finished
      this->signals.post_distributed_refinement();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::repartition()
    {
#  ifdef DEBUG
      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_locally_owned())
          Assert(
            !cell->refine_flag_set() && !cell->coarsen_flag_set(),
            ExcMessage(
              "Error: There shouldn't be any cells flagged for coarsening/refinement when calling repartition()."));
#  endif

      // signal that repartitioning is going to happen
      this->signals.pre_distributed_repartition();

      // before repartitioning the mesh, save a copy of the current positions of
      // quadrants
      // only if data needs to be transferred later
      std::vector<typename dealii::internal::p4est::types<dim>::gloidx>
        previous_global_first_quadrant;

      if (this->cell_attached_data.n_attached_data_sets > 0)
        {
          previous_global_first_quadrant.resize(parallel_forest->mpisize + 1);
          std::memcpy(previous_global_first_quadrant.data(),
                      parallel_forest->global_first_quadrant,
                      sizeof(
                        typename dealii::internal::p4est::types<dim>::gloidx) *
                        (parallel_forest->mpisize + 1));
        }

      if (this->signals.cell_weight.num_slots() == 0)
        {
          // no cell weights given -- call p4est's 'partition' without a
          // callback for cell weights
          dealii::internal::p4est::functions<dim>::partition(
            parallel_forest,
            /* prepare coarsening */ 1,
            /* weight_callback */ nullptr);
        }
      else
        {
          // get cell weights for a weighted repartitioning.
          const std::vector<unsigned int> cell_weights = get_cell_weights();

          PartitionWeights<dim, spacedim> partition_weights(cell_weights);

          // attach (temporarily) a pointer to the cell weights through p4est's
          // user_pointer object
          Assert(parallel_forest->user_pointer == this, ExcInternalError());
          parallel_forest->user_pointer = &partition_weights;

          dealii::internal::p4est::functions<dim>::partition(
            parallel_forest,
            /* prepare coarsening */ 1,
            /* weight_callback */
            &PartitionWeights<dim, spacedim>::cell_weight);

          // reset the user pointer to its previous state
          parallel_forest->user_pointer = this;
        }

      // pack data before triangulation gets updated
      if (this->cell_attached_data.n_attached_data_sets > 0)
        {
          this->data_transfer.pack_data(
            this->local_cell_relations,
            this->cell_attached_data.pack_callbacks_fixed,
            this->cell_attached_data.pack_callbacks_variable);
        }

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      // transfer data after triangulation got updated
      if (this->cell_attached_data.n_attached_data_sets > 0)
        {
          this->execute_transfer(parallel_forest,
                                 previous_global_first_quadrant.data());
        }

      this->update_periodic_face_map();

      // update how many cells, edges, etc, we store locally
      this->update_number_cache();

      // signal that repartitioning is finished
      this->signals.post_distributed_repartition();
    }



    template <int dim, int spacedim>
    const std::vector<types::global_dof_index> &
    Triangulation<dim, spacedim>::get_p4est_tree_to_coarse_cell_permutation()
      const
    {
      return p4est_tree_to_coarse_cell_permutation;
    }



    template <int dim, int spacedim>
    const std::vector<types::global_dof_index> &
    Triangulation<dim, spacedim>::get_coarse_cell_to_p4est_tree_permutation()
      const
    {
      return coarse_cell_to_p4est_tree_permutation;
    }



    template <int dim, int spacedim>
    std::vector<bool>
    Triangulation<dim, spacedim>::mark_locally_active_vertices_on_level(
      const int level) const
    {
      Assert(dim > 1, ExcNotImplemented());

      std::vector<bool> marked_vertices(this->n_vertices(), false);
      for (const auto &cell : this->cell_iterators_on_level(level))
        if (cell->level_subdomain_id() == this->locally_owned_subdomain())
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            marked_vertices[cell->vertex_index(v)] = true;

      /**
       * ensure that if one of the two vertices on a periodic face is marked
       * as active (i.e., belonging to an owned level cell), also the other
       * one is active
       */
      // When a connectivity in the code below is detected, the assignment
      // 'marked_vertices[v1] = marked_vertices[v2] = true' makes sure that
      // the information about the periodicity propagates back to vertices on
      // cells that are not owned locally. However, in the worst case we want
      // to connect to a vertex that is 'dim' hops away from the locally owned
      // cell. Depending on the order of the periodic face map, we might
      // connect to that point by chance or miss it. However, after looping
      // through all the periodic directions (which are at most as many as
      // the number of space dimensions) we can be sure that all connections
      // to vertices have been created.
      for (unsigned int repetition = 0; repetition < dim; ++repetition)
        for (const auto &it : this->get_periodic_face_map())
          {
            const cell_iterator & cell_1           = it.first.first;
            const unsigned int    face_no_1        = it.first.second;
            const cell_iterator & cell_2           = it.second.first.first;
            const unsigned int    face_no_2        = it.second.first.second;
            const std::bitset<3> &face_orientation = it.second.second;

            if (cell_1->level() == level && cell_2->level() == level)
              {
                for (unsigned int v = 0;
                     v < GeometryInfo<dim - 1>::vertices_per_cell;
                     ++v)
                  {
                    // take possible non-standard orientation of faces into
                    // account
                    const unsigned int vface0 =
                      GeometryInfo<dim>::standard_to_real_face_vertex(
                        v,
                        face_orientation[0],
                        face_orientation[1],
                        face_orientation[2]);
                    if (marked_vertices[cell_1->face(face_no_1)->vertex_index(
                          vface0)] ||
                        marked_vertices[cell_2->face(face_no_2)->vertex_index(
                          v)])
                      marked_vertices[cell_1->face(face_no_1)->vertex_index(
                        vface0)] =
                        marked_vertices[cell_2->face(face_no_2)->vertex_index(
                          v)] = true;
                  }
              }
          }

      return marked_vertices;
    }



    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim, spacedim>::coarse_cell_id_to_coarse_cell_index(
      const types::coarse_cell_id coarse_cell_id) const
    {
      return p4est_tree_to_coarse_cell_permutation[coarse_cell_id];
    }



    template <int dim, int spacedim>
    types::coarse_cell_id
    Triangulation<dim, spacedim>::coarse_cell_index_to_coarse_cell_id(
      const unsigned int coarse_cell_index) const
    {
      return coarse_cell_to_p4est_tree_permutation[coarse_cell_index];
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::add_periodicity(
      const std::vector<dealii::GridTools::PeriodicFacePair<cell_iterator>>
        &periodicity_vector)
    {
      Assert(triangulation_has_content == true,
             ExcMessage("The triangulation is empty!"));
      Assert(this->n_levels() == 1,
             ExcMessage("The triangulation is refined!"));

      // call the base class for storing the periodicity information; we must
      // do this before going to p4est and rebuilding the triangulation to get
      // the level subdomain ids correct in the multigrid case
      dealii::Triangulation<dim, spacedim>::add_periodicity(periodicity_vector);

      for (const auto &face_pair : periodicity_vector)
        {
          const cell_iterator first_cell  = face_pair.cell[0];
          const cell_iterator second_cell = face_pair.cell[1];
          const unsigned int  face_left   = face_pair.face_idx[0];
          const unsigned int  face_right  = face_pair.face_idx[1];

          // respective cells of the matching faces in p4est
          const unsigned int tree_left =
            coarse_cell_to_p4est_tree_permutation[first_cell->index()];
          const unsigned int tree_right =
            coarse_cell_to_p4est_tree_permutation[second_cell->index()];

          // p4est wants to know which corner the first corner on
          // the face with the lower id is mapped to on the face with
          // with the higher id. For d==2 there are only two possibilities
          // that are determined by it->orientation[1].
          // For d==3 we have to use GridTools::OrientationLookupTable.
          // The result is given below.

          unsigned int p4est_orientation = 0;
          if (dim == 2)
            p4est_orientation = face_pair.orientation[1];
          else
            {
              const unsigned int  face_idx_list[] = {face_left, face_right};
              const cell_iterator cell_list[]     = {first_cell, second_cell};
              unsigned int        lower_idx, higher_idx;
              if (face_left <= face_right)
                {
                  higher_idx = 1;
                  lower_idx  = 0;
                }
              else
                {
                  higher_idx = 0;
                  lower_idx  = 1;
                }

              // get the cell index of the first index on the face with the
              // lower id
              unsigned int first_p4est_idx_on_cell =
                p8est_face_corners[face_idx_list[lower_idx]][0];
              unsigned int first_dealii_idx_on_face =
                numbers::invalid_unsigned_int;
              for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face;
                   ++i)
                {
                  const unsigned int first_dealii_idx_on_cell =
                    GeometryInfo<dim>::face_to_cell_vertices(
                      face_idx_list[lower_idx],
                      i,
                      cell_list[lower_idx]->face_orientation(
                        face_idx_list[lower_idx]),
                      cell_list[lower_idx]->face_flip(face_idx_list[lower_idx]),
                      cell_list[lower_idx]->face_rotation(
                        face_idx_list[lower_idx]));
                  if (first_p4est_idx_on_cell == first_dealii_idx_on_cell)
                    {
                      first_dealii_idx_on_face = i;
                      break;
                    }
                }
              Assert(first_dealii_idx_on_face != numbers::invalid_unsigned_int,
                     ExcInternalError());
              // Now map dealii_idx_on_face according to the orientation
              constexpr unsigned int left_to_right[8][4] = {{0, 2, 1, 3},
                                                            {0, 1, 2, 3},
                                                            {3, 1, 2, 0},
                                                            {3, 2, 1, 0},
                                                            {2, 3, 0, 1},
                                                            {1, 3, 0, 2},
                                                            {1, 0, 3, 2},
                                                            {2, 0, 3, 1}};
              constexpr unsigned int right_to_left[8][4] = {{0, 2, 1, 3},
                                                            {0, 1, 2, 3},
                                                            {3, 1, 2, 0},
                                                            {3, 2, 1, 0},
                                                            {2, 3, 0, 1},
                                                            {2, 0, 3, 1},
                                                            {1, 0, 3, 2},
                                                            {1, 3, 0, 2}};
              const unsigned int     second_dealii_idx_on_face =
                lower_idx == 0 ? left_to_right[face_pair.orientation.to_ulong()]
                                              [first_dealii_idx_on_face] :
                                 right_to_left[face_pair.orientation.to_ulong()]
                                              [first_dealii_idx_on_face];
              const unsigned int second_dealii_idx_on_cell =
                GeometryInfo<dim>::face_to_cell_vertices(
                  face_idx_list[higher_idx],
                  second_dealii_idx_on_face,
                  cell_list[higher_idx]->face_orientation(
                    face_idx_list[higher_idx]),
                  cell_list[higher_idx]->face_flip(face_idx_list[higher_idx]),
                  cell_list[higher_idx]->face_rotation(
                    face_idx_list[higher_idx]));
              // map back to p4est
              const unsigned int second_p4est_idx_on_face =
                p8est_corner_face_corners[second_dealii_idx_on_cell]
                                         [face_idx_list[higher_idx]];
              p4est_orientation = second_p4est_idx_on_face;
            }

          dealii::internal::p4est::functions<dim>::connectivity_join_faces(
            connectivity,
            tree_left,
            tree_right,
            face_left,
            face_right,
            p4est_orientation);
        }


      Assert(dealii::internal::p4est::functions<dim>::connectivity_is_valid(
               connectivity) == 1,
             ExcInternalError());

      // now create a forest out of the connectivity data structure
      dealii::internal::p4est::functions<dim>::destroy(parallel_forest);
      parallel_forest = dealii::internal::p4est::functions<dim>::new_forest(
        this->mpi_communicator,
        connectivity,
        /* minimum initial number of quadrants per tree */ 0,
        /* minimum level of upfront refinement */ 0,
        /* use uniform upfront refinement */ 1,
        /* user_data_size = */ 0,
        /* user_data_constructor = */ nullptr,
        /* user_pointer */ this);

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      // The range of ghost_owners might have changed so update that information
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    std::size_t
    Triangulation<dim, spacedim>::memory_consumption() const
    {
      std::size_t mem =
        this->dealii::parallel::TriangulationBase<dim, spacedim>::
          memory_consumption() +
        MemoryConsumption::memory_consumption(triangulation_has_content) +
        MemoryConsumption::memory_consumption(connectivity) +
        MemoryConsumption::memory_consumption(parallel_forest) +
        MemoryConsumption::memory_consumption(
          this->cell_attached_data.n_attached_data_sets) +
        // MemoryConsumption::memory_consumption(cell_attached_data.pack_callbacks_fixed)
        // +
        // MemoryConsumption::memory_consumption(cell_attached_data.pack_callbacks_variable)
        // +
        // TODO[TH]: how?
        MemoryConsumption::memory_consumption(
          coarse_cell_to_p4est_tree_permutation) +
        MemoryConsumption::memory_consumption(
          p4est_tree_to_coarse_cell_permutation) +
        memory_consumption_p4est();

      return mem;
    }



    template <int dim, int spacedim>
    std::size_t
    Triangulation<dim, spacedim>::memory_consumption_p4est() const
    {
      return dealii::internal::p4est::functions<dim>::forest_memory_used(
               parallel_forest) +
             dealii::internal::p4est::functions<dim>::connectivity_memory_used(
               connectivity);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &other_tria)
    {
      Assert(
        (dynamic_cast<
          const dealii::parallel::distributed::Triangulation<dim, spacedim> *>(
          &other_tria)) ||
          (other_tria.n_global_levels() == 1),
        ExcNotImplemented());

      dealii::parallel::distributed::Triangulation<dim, spacedim>::clear();

      try
        {
          dealii::parallel::TriangulationBase<dim, spacedim>::
            copy_triangulation(other_tria);
        }
      catch (
        const typename dealii::Triangulation<dim, spacedim>::DistortedCellList
          &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      if (const dealii::parallel::distributed::Triangulation<dim, spacedim>
            *other_distributed =
              dynamic_cast<const dealii::parallel::distributed::
                             Triangulation<dim, spacedim> *>(&other_tria))
        {
          // copy parallel distributed specifics
          settings = other_distributed->settings;
          triangulation_has_content =
            other_distributed->triangulation_has_content;
          coarse_cell_to_p4est_tree_permutation =
            other_distributed->coarse_cell_to_p4est_tree_permutation;
          p4est_tree_to_coarse_cell_permutation =
            other_distributed->p4est_tree_to_coarse_cell_permutation;
          this->cell_attached_data = other_distributed->cell_attached_data;
          this->data_transfer      = other_distributed->data_transfer;

          // create deep copy of connectivity graph
          typename dealii::internal::p4est::types<dim>::connectivity
            *temp_connectivity = const_cast<
              typename dealii::internal::p4est::types<dim>::connectivity *>(
              other_distributed->connectivity);
          connectivity =
            dealii::internal::p4est::copy_connectivity<dim>(temp_connectivity);

          // create deep copy of parallel forest
          typename dealii::internal::p4est::types<dim>::forest *temp_forest =
            const_cast<typename dealii::internal::p4est::types<dim>::forest *>(
              other_distributed->parallel_forest);
          parallel_forest =
            dealii::internal::p4est::functions<dim>::copy_forest(temp_forest,
                                                                 false);
          parallel_forest->connectivity = connectivity;
          parallel_forest->user_pointer = this;
        }
      else
        {
          triangulation_has_content = true;
          setup_coarse_cell_to_p4est_tree_permutation();
          copy_new_triangulation_to_p4est(std::integral_constant<int, dim>());
        }

      try
        {
          copy_local_forest_to_triangulation();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }

      this->update_periodic_face_map();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::update_cell_relations()
    {
      // reorganize memory for local_cell_relations
      this->local_cell_relations.resize(parallel_forest->local_num_quadrants);
      this->local_cell_relations.shrink_to_fit();

      // recurse over p4est
      for (const auto &cell : this->cell_iterators_on_level(0))
        {
          // skip coarse cells that are not ours
          if (tree_exists_locally<dim, spacedim>(
                parallel_forest,
                coarse_cell_to_p4est_tree_permutation[cell->index()]) == false)
            continue;

          // initialize auxiliary top level p4est quadrant
          typename dealii::internal::p4est::types<dim>::quadrant
            p4est_coarse_cell;
          dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

          // determine tree to start recursion on
          typename dealii::internal::p4est::types<dim>::tree *tree =
            init_tree(cell->index());

          update_cell_relations_recursively<dim, spacedim>(
            this->local_cell_relations, *tree, cell, p4est_coarse_cell);
        }
    }



    template <int dim, int spacedim>
    std::vector<unsigned int>
    Triangulation<dim, spacedim>::get_cell_weights() const
    {
      // check if local_cell_relations have been previously gathered
      // correctly
      Assert(this->local_cell_relations.size() ==
               static_cast<unsigned int>(parallel_forest->local_num_quadrants),
             ExcInternalError());

      // Allocate the space for the weights. In fact we do not know yet, how
      // many cells we own after the refinement (only p4est knows that
      // at this point). We simply reserve n_active_cells space and if many
      // more cells are refined than coarsened than additional reallocation
      // will be done inside get_cell_weights_recursively.
      std::vector<unsigned int> weights;
      weights.reserve(this->n_active_cells());

      // Iterate over p4est and Triangulation relations
      // to find refined/coarsened/kept
      // cells. Then append cell_weight.
      // Note that we need to follow the p4est ordering
      // instead of the deal.II ordering to get the cell_weights
      // in the same order p4est will encounter them during repartitioning.
      for (const auto &cell_rel : this->local_cell_relations)
        {
          const auto &cell_it     = cell_rel.first;
          const auto &cell_status = cell_rel.second;

          switch (cell_status)
            {
              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_PERSIST:
                weights.push_back(1000);
                weights.back() += this->signals.cell_weight(
                  cell_it,
                  parallel::distributed::Triangulation<dim,
                                                       spacedim>::CELL_PERSIST);
                break;

              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_REFINE:
              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_INVALID:
                {
                  // calculate weight of parent cell
                  unsigned int parent_weight = 1000;
                  parent_weight += this->signals.cell_weight(
                    cell_it,
                    parallel::distributed::Triangulation<dim, spacedim>::
                      CELL_REFINE);
                  // assign the weight of the parent cell equally to all
                  // children
                  weights.push_back(parent_weight);
                  break;
                }

              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_COARSEN:
                weights.push_back(1000);
                weights.back() += this->signals.cell_weight(
                  cell_it,
                  parallel::distributed::Triangulation<dim,
                                                       spacedim>::CELL_COARSEN);
                break;

              default:
                Assert(false, ExcInternalError());
                break;
            }
        }

      return weights;
    }



    template <int spacedim>
    Triangulation<1, spacedim>::Triangulation(
      const MPI_Comm &mpi_communicator,
      const typename dealii::Triangulation<1, spacedim>::MeshSmoothing
        smooth_grid,
      const Settings /*settings*/)
      : dealii::parallel::DistributedTriangulationBase<1, spacedim>(
          mpi_communicator,
          smooth_grid,
          false)
    {
      Assert(false, ExcNotImplemented());
    }


    template <int spacedim>
    Triangulation<1, spacedim>::~Triangulation()
    {
      AssertNothrow(false, ExcNotImplemented());
    }



    template <int spacedim>
    const std::vector<types::global_dof_index> &
    Triangulation<1, spacedim>::get_p4est_tree_to_coarse_cell_permutation()
      const
    {
      static std::vector<types::global_dof_index> a;
      return a;
    }



    template <int spacedim>
    std::map<unsigned int, std::set<dealii::types::subdomain_id>>
    Triangulation<1, spacedim>::compute_level_vertices_with_ghost_neighbors(
      const unsigned int /*level*/) const
    {
      Assert(false, ExcNotImplemented());

      return std::map<unsigned int, std::set<dealii::types::subdomain_id>>();
    }



    template <int spacedim>
    std::vector<bool>
    Triangulation<1, spacedim>::mark_locally_active_vertices_on_level(
      const unsigned int) const
    {
      Assert(false, ExcNotImplemented());
      return std::vector<bool>();
    }



    template <int spacedim>
    unsigned int
    Triangulation<1, spacedim>::coarse_cell_id_to_coarse_cell_index(
      const types::coarse_cell_id) const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <int spacedim>
    types::coarse_cell_id
    Triangulation<1, spacedim>::coarse_cell_index_to_coarse_cell_id(
      const unsigned int) const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <int spacedim>
    void
    Triangulation<1, spacedim>::load(const std::string &, const bool)
    {
      Assert(false, ExcNotImplemented());
    }



    template <int spacedim>
    void
    Triangulation<1, spacedim>::save(const std::string &) const
    {
      Assert(false, ExcNotImplemented());
    }



    template <int spacedim>
    bool
    Triangulation<1, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      Assert(false, ExcNotImplemented());
      return false;
    }



    template <int spacedim>
    void
    Triangulation<1, spacedim>::update_cell_relations()
    {
      Assert(false, ExcNotImplemented());
    }

  } // namespace distributed
} // namespace parallel


#endif // DEAL_II_WITH_P4EST



namespace internal
{
  namespace parallel
  {
    namespace distributed
    {
      template <int dim, int spacedim>
      TemporarilyMatchRefineFlags<dim, spacedim>::TemporarilyMatchRefineFlags(
        Triangulation<dim, spacedim> &tria)
        : distributed_tria(
            dynamic_cast<
              dealii::parallel::distributed::Triangulation<dim, spacedim> *>(
              &tria))
      {
#ifdef DEAL_II_WITH_P4EST
        if (distributed_tria != nullptr)
          {
            // Save the current set of refinement flags, and adjust the
            // refinement flags to be consistent with the p4est oracle.
            distributed_tria->save_coarsen_flags(saved_coarsen_flags);
            distributed_tria->save_refine_flags(saved_refine_flags);

            for (const auto &pair : distributed_tria->local_cell_relations)
              {
                const auto &cell   = pair.first;
                const auto &status = pair.second;

                switch (status)
                  {
                    case dealii::Triangulation<dim, spacedim>::CELL_PERSIST:
                      // cell remains unchanged
                      cell->clear_refine_flag();
                      cell->clear_coarsen_flag();
                      break;

                    case dealii::Triangulation<dim, spacedim>::CELL_REFINE:
                      // cell will be refined
                      cell->clear_coarsen_flag();
                      cell->set_refine_flag();
                      break;

                    case dealii::Triangulation<dim, spacedim>::CELL_COARSEN:
                      // children of this cell will be coarsened
                      for (const auto &child : cell->child_iterators())
                        {
                          child->clear_refine_flag();
                          child->set_coarsen_flag();
                        }
                      break;

                    case dealii::Triangulation<dim, spacedim>::CELL_INVALID:
                      // do nothing as cell does not exist yet
                      break;

                    default:
                      Assert(false, ExcInternalError());
                      break;
                  }
              }
          }
#endif
      }



      template <int dim, int spacedim>
      TemporarilyMatchRefineFlags<dim, spacedim>::~TemporarilyMatchRefineFlags()
      {
#ifdef DEAL_II_WITH_P4EST
        if (distributed_tria)
          {
            // Undo the refinement flags modification.
            distributed_tria->load_coarsen_flags(saved_coarsen_flags);
            distributed_tria->load_refine_flags(saved_refine_flags);
          }
#else
        // pretend that this destructor does something to silence clang-tidy
        (void)distributed_tria;
#endif
      }
    } // namespace distributed
  }   // namespace parallel
} // namespace internal



/*-------------- Explicit Instantiations -------------------------------*/
#include "tria.inst"


DEAL_II_NAMESPACE_CLOSE
