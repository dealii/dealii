// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/mapping_collection.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  unsigned int find_closest_vertex(const MeshType<dim, spacedim> &mesh,
                                   const Point<spacedim>         &p,
                                   const std::vector<bool> &marked_vertices)
  {
    // first get the underlying
    // triangulation from the
    // mesh and determine vertices
    // and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const std::vector<Point<spacedim>> &vertices = tria.get_vertices();

    Assert(tria.get_vertices().size() == marked_vertices.size() ||
             marked_vertices.empty(),
           ExcDimensionMismatch(tria.get_vertices().size(),
                                marked_vertices.size()));

    // If p is an element of marked_vertices,
    // and q is that of used_Vertices,
    // the vector marked_vertices does NOT
    // contain unused vertices if p implies q.
    // I.e., if p is true q must be true
    // (if p is false, q could be false or true).
    // p implies q logic is encapsulated in ~p|q.
    Assert(
      marked_vertices.empty() ||
        std::equal(marked_vertices.begin(),
                   marked_vertices.end(),
                   tria.get_used_vertices().begin(),
                   [](bool p, bool q) { return !p || q; }),
      ExcMessage(
        "marked_vertices should be a subset of used vertices in the triangulation "
        "but marked_vertices contains one or more vertices that are not used vertices!"));

    // In addition, if a vector bools
    // is specified (marked_vertices)
    // marking all the vertices which
    // could be the potentially closest
    // vertex to the point, use it instead
    // of used vertices
    const std::vector<bool> &used =
      (marked_vertices.empty()) ? tria.get_used_vertices() : marked_vertices;

    // At the beginning, the first
    // used vertex is the closest one
    std::vector<bool>::const_iterator first =
      std::find(used.begin(), used.end(), true);

    // Assert that at least one vertex
    // is actually used
    Assert(first != used.end(), ExcInternalError());

    unsigned int best_vertex = std::distance(used.begin(), first);
    double       best_dist   = (p - vertices[best_vertex]).norm_square();

    // For all remaining vertices, test
    // whether they are any closer
    for (unsigned int j = best_vertex + 1; j < vertices.size(); ++j)
      if (used[j])
        {
          double dist = (p - vertices[j]).norm_square();
          if (dist < best_dist)
            {
              best_vertex = j;
              best_dist   = dist;
            }
        }

    return best_vertex;
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  unsigned int find_closest_vertex(const Mapping<dim, spacedim>  &mapping,
                                   const MeshType<dim, spacedim> &mesh,
                                   const Point<spacedim>         &p,
                                   const std::vector<bool> &marked_vertices)
  {
    // Take a shortcut in the simple case.
    if (mapping.preserves_vertex_locations() == true)
      return find_closest_vertex(mesh, p, marked_vertices);

    // first get the underlying
    // triangulation from the
    // mesh and determine vertices
    // and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    auto vertices = extract_used_vertices(tria, mapping);

    Assert(tria.get_vertices().size() == marked_vertices.size() ||
             marked_vertices.empty(),
           ExcDimensionMismatch(tria.get_vertices().size(),
                                marked_vertices.size()));

    // If p is an element of marked_vertices,
    // and q is that of used_Vertices,
    // the vector marked_vertices does NOT
    // contain unused vertices if p implies q.
    // I.e., if p is true q must be true
    // (if p is false, q could be false or true).
    // p implies q logic is encapsulated in ~p|q.
    Assert(
      marked_vertices.empty() ||
        std::equal(marked_vertices.begin(),
                   marked_vertices.end(),
                   tria.get_used_vertices().begin(),
                   [](bool p, bool q) { return !p || q; }),
      ExcMessage(
        "marked_vertices should be a subset of used vertices in the triangulation "
        "but marked_vertices contains one or more vertices that are not used vertices!"));

    // Remove from the map unwanted elements.
    if (marked_vertices.size())
      for (auto it = vertices.begin(); it != vertices.end();)
        {
          if (marked_vertices[it->first] == false)
            {
              vertices.erase(it++);
            }
          else
            {
              ++it;
            }
        }

    return find_closest_vertex(vertices, p);
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES((concepts::is_triangulation_or_dof_handler<MeshType>))
#ifndef _MSC_VER
  std::vector<typename MeshType::active_cell_iterator>
#else
  std::vector<
    typename dealii::internal::ActiveCellIterator<MeshType::dimension,
                                                  MeshType::space_dimension,
                                                  MeshType>::type>
#endif
    find_cells_adjacent_to_vertex(const MeshType    &mesh,
                                  const unsigned int vertex)
  {
    const int dim      = MeshType::dimension;
    const int spacedim = MeshType::space_dimension;

    // make sure that the given vertex is
    // an active vertex of the underlying
    // triangulation
    AssertIndexRange(vertex, mesh.get_triangulation().n_vertices());
    Assert(mesh.get_triangulation().get_used_vertices()[vertex],
           ExcVertexNotUsed(vertex));

    // use a set instead of a vector
    // to ensure that cells are inserted only
    // once
    std::set<typename dealii::internal::
               ActiveCellIterator<dim, spacedim, MeshType>::type>
      adjacent_cells;

    typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType>::type
      cell = mesh.begin_active(),
      endc = mesh.end();

    // go through all active cells and look if the vertex is part of that cell
    //
    // in 1d, this is all we need to care about. in 2d/3d we also need to worry
    // that the vertex might be a hanging node on a face or edge of a cell; in
    // this case, we would want to add those cells as well on whose faces the
    // vertex is located but for which it is not a vertex itself.
    //
    // getting this right is a lot simpler in 2d than in 3d. in 2d, a hanging
    // node can only be in the middle of a face and we can query the neighboring
    // cell from the current cell. on the other hand, in 3d a hanging node
    // vertex can also be on an edge but there can be many other cells on
    // this edge and we can not access them from the cell we are currently
    // on.
    //
    // so, in the 3d case, if we run the algorithm as in 2d, we catch all
    // those cells for which the vertex we seek is on a *subface*, but we
    // miss the case of cells for which the vertex we seek is on a
    // sub-edge for which there is no corresponding sub-face (because the
    // immediate neighbor behind this face is not refined), see for example
    // the bits/find_cells_adjacent_to_vertex_6 testcase. thus, if we
    // haven't yet found the vertex for the current cell we also need to
    // look at the mid-points of edges
    //
    // as a final note, deciding whether a neighbor is actually coarser is
    // simple in the case of isotropic refinement (we just need to look at
    // the level of the current and the neighboring cell). however, this
    // isn't so simple if we have used anisotropic refinement since then
    // the level of a cell is not indicative of whether it is coarser or
    // not than the current cell. ultimately, we want to add all cells on
    // which the vertex is, independent of whether they are coarser or
    // finer and so in the 2d case below we simply add *any* *active* neighbor.
    // in the worst case, we add cells multiple times to the adjacent_cells
    // list, but std::set throws out those cells already entered
    for (; cell != endc; ++cell)
      {
        for (const unsigned int v : cell->vertex_indices())
          if (cell->vertex_index(v) == vertex)
            {
              // OK, we found a cell that contains
              // the given vertex. We add it
              // to the list.
              adjacent_cells.insert(cell);

              // as explained above, in 2+d we need to check whether
              // this vertex is on a face behind which there is a
              // (possibly) coarser neighbor. if this is the case,
              // then we need to also add this neighbor
              if (dim >= 2)
                {
                  const auto reference_cell = cell->reference_cell();
                  for (const auto face :
                       reference_cell.faces_for_given_vertex(v))
                    if (!cell->at_boundary(face) &&
                        cell->neighbor(face)->is_active())
                      {
                        // there is a (possibly) coarser cell behind a
                        // face to which the vertex belongs. the
                        // vertex we are looking at is then either a
                        // vertex of that coarser neighbor, or it is a
                        // hanging node on one of the faces of that
                        // cell. in either case, it is adjacent to the
                        // vertex, so add it to the list as well (if
                        // the cell was already in the list then the
                        // std::set makes sure that we get it only
                        // once)
                        adjacent_cells.insert(cell->neighbor(face));
                      }
                }

              // in any case, we have found a cell, so go to the next cell
              goto next_cell;
            }

        // in 3d also loop over the edges
        if (dim >= 3)
          {
            for (unsigned int e = 0; e < cell->n_lines(); ++e)
              if (cell->line(e)->has_children())
                // the only place where this vertex could have been
                // hiding is on the mid-edge point of the edge we
                // are looking at
                if (cell->line(e)->child(0)->vertex_index(1) == vertex)
                  {
                    adjacent_cells.insert(cell);

                    // jump out of this tangle of nested loops
                    goto next_cell;
                  }
          }

        // in more than 3d we would probably have to do the same as
        // above also for even lower-dimensional objects
        Assert(dim <= 3, ExcNotImplemented());

        // move on to the next cell if we have found the
        // vertex on the current one
      next_cell:;
      }

    // if this was an active vertex then there needs to have been
    // at least one cell to which it is adjacent!
    Assert(adjacent_cells.size() > 0, ExcInternalError());

    // return the result as a vector, rather than the set we built above
    return std::vector<typename dealii::internal::
                         ActiveCellIterator<dim, spacedim, MeshType>::type>(
      adjacent_cells.begin(), adjacent_cells.end());
  }



  namespace
  {
    template <int dim, template <int, int> class MeshType, int spacedim>
    DEAL_II_CXX20_REQUIRES(
      (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
    void find_active_cell_around_point_internal(
      const MeshType<dim, spacedim> &mesh,
#ifndef _MSC_VER
      std::set<typename MeshType<dim, spacedim>::active_cell_iterator>
        &searched_cells,
      std::set<typename MeshType<dim, spacedim>::active_cell_iterator>
        &adjacent_cells)
#else
      std::set<
        typename dealii::internal::
          ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type>
        &searched_cells,
      std::set<
        typename dealii::internal::
          ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type>
        &adjacent_cells)
#endif
    {
#ifndef _MSC_VER
      using cell_iterator =
        typename MeshType<dim, spacedim>::active_cell_iterator;
#else
      using cell_iterator = typename dealii::internal::
        ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type;
#endif

      // update the searched cells
      searched_cells.insert(adjacent_cells.begin(), adjacent_cells.end());
      // now we to collect all neighbors
      // of the cells in adjacent_cells we
      // have not yet searched.
      std::set<cell_iterator> adjacent_cells_new;

      for (const auto &cell : adjacent_cells)
        {
          std::vector<cell_iterator> active_neighbors;
          get_active_neighbors<MeshType<dim, spacedim>>(cell, active_neighbors);
          for (unsigned int i = 0; i < active_neighbors.size(); ++i)
            if (searched_cells.find(active_neighbors[i]) ==
                searched_cells.end())
              adjacent_cells_new.insert(active_neighbors[i]);
        }
      adjacent_cells.clear();
      adjacent_cells.insert(adjacent_cells_new.begin(),
                            adjacent_cells_new.end());
      if (adjacent_cells.empty())
        {
          // we haven't found any other cell that would be a
          // neighbor of a previously found cell, but we know
          // that we haven't checked all cells yet. that means
          // that the domain is disconnected. in that case,
          // choose the first previously untouched cell we
          // can find
          cell_iterator it = mesh.begin_active();
          for (; it != mesh.end(); ++it)
            if (searched_cells.find(it) == searched_cells.end())
              {
                adjacent_cells.insert(it);
                break;
              }
        }
    }
  } // namespace



  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  typename MeshType<dim, spacedim>::active_cell_iterator
#else
  typename dealii::internal::
    ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type
#endif
    find_active_cell_around_point(const MeshType<dim, spacedim> &mesh,
                                  const Point<spacedim>         &p,
                                  const std::vector<bool> &marked_vertices,
                                  const double             tolerance)
  {
    return find_active_cell_around_point<dim, MeshType, spacedim>(
             get_default_linear_mapping(mesh.get_triangulation()),
             mesh,
             p,
             marked_vertices,
             tolerance)
      .first;
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#endif
    find_active_cell_around_point(const Mapping<dim, spacedim>  &mapping,
                                  const MeshType<dim, spacedim> &mesh,
                                  const Point<spacedim>         &p,
                                  const std::vector<bool> &marked_vertices,
                                  const double             tolerance)
  {
    using active_cell_iterator = typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type;

    // The best distance is set to the
    // maximum allowable distance from
    // the unit cell; we assume a
    // max. deviation of the given tolerance
    double                                      best_distance = tolerance;
    int                                         best_level    = -1;
    std::pair<active_cell_iterator, Point<dim>> best_cell;

    // Initialize best_cell.first to the end iterator
    best_cell.first = mesh.end();

    // Find closest vertex and determine
    // all adjacent cells
    std::vector<active_cell_iterator> adjacent_cells_tmp =
      find_cells_adjacent_to_vertex(
        mesh, find_closest_vertex(mapping, mesh, p, marked_vertices));

    // Make sure that we have found
    // at least one cell adjacent to vertex.
    Assert(adjacent_cells_tmp.size() > 0, ExcInternalError());

    // Copy all the cells into a std::set
    std::set<active_cell_iterator> adjacent_cells(adjacent_cells_tmp.begin(),
                                                  adjacent_cells_tmp.end());
    std::set<active_cell_iterator> searched_cells;

    // Determine the maximal number of cells
    // in the grid.
    // As long as we have not found
    // the cell and have not searched
    // every cell in the triangulation,
    // we keep on looking.
    const auto   n_active_cells = mesh.get_triangulation().n_active_cells();
    bool         found          = false;
    unsigned int cells_searched = 0;
    while (!found && cells_searched < n_active_cells)
      {
        for (const auto &cell : adjacent_cells)
          {
            if (cell->is_artificial() == false)
              {
                // marked_vertices are used to filter cell candidates
                if (marked_vertices.size() > 0)
                  {
                    bool any_vertex_marked = false;
                    for (const auto &v : cell->vertex_indices())
                      {
                        if (marked_vertices[cell->vertex_index(v)])
                          {
                            any_vertex_marked = true;
                            break;
                          }
                      }
                    if (!any_vertex_marked)
                      continue;
                  }

                try
                  {
                    const Point<dim> p_cell =
                      mapping.transform_real_to_unit_cell(cell, p);

                    // calculate the Euclidean norm of
                    // the distance vector to the unit cell.
                    const double dist =
                      cell->reference_cell().closest_point(p_cell).distance(
                        p_cell);

                    // We compare if the point is inside the
                    // unit cell (or at least not too far
                    // outside). If it is, it is also checked
                    // that the cell has a more refined state
                    if ((dist < best_distance) ||
                        ((dist == best_distance) &&
                         (cell->level() > best_level)))
                      {
                        found         = true;
                        best_distance = dist;
                        best_level    = cell->level();
                        best_cell     = std::make_pair(cell, p_cell);
                      }
                  }
                catch (
                  typename MappingQ<dim, spacedim>::ExcTransformationFailed &)
                  {
                    // ok, the transformation
                    // failed presumably
                    // because the point we
                    // are looking for lies
                    // outside the current
                    // cell. this means that
                    // the current cell can't
                    // be the cell around the
                    // point, so just ignore
                    // this cell and move on
                    // to the next
                  }
              }
          }

        // update the number of cells searched
        cells_searched += adjacent_cells.size();

        // if we have not found the cell in
        // question and have not yet searched every
        // cell, we expand our search to
        // all the not already searched neighbors of
        // the cells in adjacent_cells. This is
        // what find_active_cell_around_point_internal
        // is for.
        if (!found && cells_searched < n_active_cells)
          {
            find_active_cell_around_point_internal<dim, MeshType, spacedim>(
              mesh, searched_cells, adjacent_cells);
          }
      }

    return best_cell;
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::vector<std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                        Point<dim>>>
#else
  std::vector<std::pair<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
    Point<dim>>>
#endif
    find_all_active_cells_around_point(const Mapping<dim, spacedim>  &mapping,
                                       const MeshType<dim, spacedim> &mesh,
                                       const Point<spacedim>         &p,
                                       const double                   tolerance,
                                       const std::vector<bool> &marked_vertices)
  {
    const auto cell_and_point = find_active_cell_around_point(
      mapping, mesh, p, marked_vertices, tolerance);

    if (cell_and_point.first == mesh.end())
      return {};

    return find_all_active_cells_around_point(
      mapping, mesh, p, tolerance, cell_and_point);
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::vector<std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                        Point<dim>>>
#else
  std::vector<std::pair<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
    Point<dim>>>
#endif
    find_all_active_cells_around_point(
      const Mapping<dim, spacedim>  &mapping,
      const MeshType<dim, spacedim> &mesh,
      const Point<spacedim>         &p,
      const double                   tolerance,
      const std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                      Point<dim>>   &first_cell,
      const std::vector<
        std::set<typename MeshType<dim, spacedim>::active_cell_iterator>>
        *vertex_to_cells)
  {
    std::vector<
      std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                Point<dim>>>
      cells_and_points;

    // insert the fist cell and point into the vector
    cells_and_points.push_back(first_cell);

    const Point<dim> unit_point = cells_and_points.front().second;
    const auto       my_cell    = cells_and_points.front().first;

    std::vector<typename MeshType<dim, spacedim>::active_cell_iterator>
      cells_to_add;

    if (my_cell->reference_cell().is_hyper_cube())
      {
        // check if the given point is on the surface of the unit cell. If yes,
        // need to find all neighbors

        Tensor<1, dim> distance_to_center;
        unsigned int   n_dirs_at_threshold     = 0;
        unsigned int   last_point_at_threshold = numbers::invalid_unsigned_int;
        for (unsigned int d = 0; d < dim; ++d)
          {
            distance_to_center[d] = std::abs(unit_point[d] - 0.5);
            if (distance_to_center[d] > 0.5 - tolerance)
              {
                ++n_dirs_at_threshold;
                last_point_at_threshold = d;
              }
          }

        // point is within face -> only need neighbor
        if (n_dirs_at_threshold == 1)
          {
            unsigned int neighbor_index =
              2 * last_point_at_threshold +
              (unit_point[last_point_at_threshold] > 0.5 ? 1 : 0);
            if (!my_cell->at_boundary(neighbor_index))
              {
                const auto neighbor_cell = my_cell->neighbor(neighbor_index);

                if (neighbor_cell->is_active())
                  cells_to_add.push_back(neighbor_cell);
                else
                  for (const auto &child_cell :
                       neighbor_cell->child_iterators())
                    {
                      if (child_cell->is_active())
                        cells_to_add.push_back(child_cell);
                    }
              }
          }
        // corner point -> use all neighbors
        else if (n_dirs_at_threshold == dim)
          {
            unsigned int local_vertex_index = 0;
            for (unsigned int d = 0; d < dim; ++d)
              local_vertex_index += (unit_point[d] > 0.5 ? 1 : 0) << d;

            const auto fu = [&](const auto &tentative_cells) {
              for (const auto &cell : tentative_cells)
                if (cell != my_cell)
                  cells_to_add.push_back(cell);
            };

            const auto vertex_index = my_cell->vertex_index(local_vertex_index);

            if (vertex_to_cells != nullptr)
              fu((*vertex_to_cells)[vertex_index]);
            else
              fu(find_cells_adjacent_to_vertex(mesh, vertex_index));
          }
        // point on line in 3d: We cannot simply take the intersection between
        // the two vertices of cells because of hanging nodes. So instead we
        // list the vertices around both points and then select the
        // appropriate cells according to the result of read_to_unit_cell
        // below.
        else if (n_dirs_at_threshold == 2)
          {
            std::pair<unsigned int, unsigned int> vertex_indices[3];
            unsigned int                          count_vertex_indices = 0;
            unsigned int free_direction = numbers::invalid_unsigned_int;
            for (unsigned int d = 0; d < dim; ++d)
              {
                if (distance_to_center[d] > 0.5 - tolerance)
                  {
                    vertex_indices[count_vertex_indices].first = d;
                    vertex_indices[count_vertex_indices].second =
                      unit_point[d] > 0.5 ? 1 : 0;
                    ++count_vertex_indices;
                  }
                else
                  free_direction = d;
              }

            AssertDimension(count_vertex_indices, 2);
            Assert(free_direction != numbers::invalid_unsigned_int,
                   ExcInternalError());

            const unsigned int first_vertex =
              (vertex_indices[0].second << vertex_indices[0].first) +
              (vertex_indices[1].second << vertex_indices[1].first);
            for (unsigned int d = 0; d < 2; ++d)
              {
                const auto fu = [&](const auto &tentative_cells) {
                  for (const auto &cell : tentative_cells)
                    {
                      bool cell_not_yet_present = true;
                      for (const auto &other_cell : cells_to_add)
                        if (cell == other_cell)
                          {
                            cell_not_yet_present = false;
                            break;
                          }
                      if (cell_not_yet_present)
                        cells_to_add.push_back(cell);
                    }
                };

                const auto vertex_index =
                  my_cell->vertex_index(first_vertex + (d << free_direction));

                if (vertex_to_cells != nullptr)
                  fu((*vertex_to_cells)[vertex_index]);
                else
                  fu(find_cells_adjacent_to_vertex(mesh, vertex_index));
              }
          }
      }
    else
      {
        // Note: The non-hypercube path takes a very naive approach and
        // checks all possible neighbors. This can be made faster by 1)
        // checking if the point is in the inner cell and 2) identifying
        // the right lines/vertices so that the number of potential
        // neighbors is reduced.

        for (const auto v : my_cell->vertex_indices())
          {
            const auto fu = [&](const auto &tentative_cells) {
              for (const auto &cell : tentative_cells)
                {
                  bool cell_not_yet_present = true;
                  for (const auto &other_cell : cells_to_add)
                    if (cell == other_cell)
                      {
                        cell_not_yet_present = false;
                        break;
                      }
                  if (cell_not_yet_present)
                    cells_to_add.push_back(cell);
                }
            };

            const auto vertex_index = my_cell->vertex_index(v);

            if (vertex_to_cells != nullptr)
              fu((*vertex_to_cells)[vertex_index]);
            else
              fu(find_cells_adjacent_to_vertex(mesh, vertex_index));
          }
      }

    for (const auto &cell : cells_to_add)
      {
        if (cell != my_cell)
          try
            {
              const Point<dim> p_unit =
                mapping.transform_real_to_unit_cell(cell, p);
              if (cell->reference_cell().contains_point(p_unit, tolerance))
                cells_and_points.emplace_back(cell, p_unit);
            }
          catch (typename Mapping<dim>::ExcTransformationFailed &)
            {}
      }

    std::sort(
      cells_and_points.begin(),
      cells_and_points.end(),
      [](const std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                         Point<dim>> &a,
         const std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                         Point<dim>> &b) { return a.first < b.first; });

    return cells_and_points;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::active_cell_iterator> compute_active_cell_halo_layer(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::active_cell_iterator &)>
        &predicate)
  {
    std::vector<typename MeshType::active_cell_iterator> active_halo_layer;
    std::vector<bool> locally_active_vertices_on_subdomain(
      mesh.get_triangulation().n_vertices(), false);

    std::map<unsigned int, std::vector<unsigned int>> coinciding_vertex_groups;
    std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;
    GridTools::collect_coinciding_vertices(mesh.get_triangulation(),
                                           coinciding_vertex_groups,
                                           vertex_to_coinciding_vertex_group);

    // Find the cells for which the predicate is true
    // These are the cells around which we wish to construct
    // the halo layer
    for (const auto &cell : mesh.active_cell_iterators())
      if (predicate(cell)) // True predicate --> Part of subdomain
        for (const auto v : cell->vertex_indices())
          {
            locally_active_vertices_on_subdomain[cell->vertex_index(v)] = true;
            for (const auto vv : coinciding_vertex_groups
                   [vertex_to_coinciding_vertex_group[cell->vertex_index(v)]])
              locally_active_vertices_on_subdomain[vv] = true;
          }

    // Find the cells that do not conform to the predicate
    // but share a vertex with the selected subdomain
    // These comprise the halo layer
    for (const auto &cell : mesh.active_cell_iterators())
      if (!predicate(cell)) // False predicate --> Potential halo cell
        for (const auto v : cell->vertex_indices())
          if (locally_active_vertices_on_subdomain[cell->vertex_index(v)] ==
              true)
            {
              active_halo_layer.push_back(cell);
              break;
            }

    return active_halo_layer;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::cell_iterator> compute_cell_halo_layer_on_level(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::cell_iterator &)>
                        &predicate,
      const unsigned int level)
  {
    std::vector<typename MeshType::cell_iterator> level_halo_layer;
    std::vector<bool> locally_active_vertices_on_level_subdomain(
      mesh.get_triangulation().n_vertices(), false);

    // Find the cells for which the predicate is true
    // These are the cells around which we wish to construct
    // the halo layer
    for (typename MeshType::cell_iterator cell = mesh.begin(level);
         cell != mesh.end(level);
         ++cell)
      if (predicate(cell)) // True predicate --> Part of subdomain
        for (const unsigned int v : cell->vertex_indices())
          locally_active_vertices_on_level_subdomain[cell->vertex_index(v)] =
            true;

    // Find the cells that do not conform to the predicate
    // but share a vertex with the selected subdomain on that level
    // These comprise the halo layer
    for (typename MeshType::cell_iterator cell = mesh.begin(level);
         cell != mesh.end(level);
         ++cell)
      if (!predicate(cell)) // False predicate --> Potential halo cell
        for (const unsigned int v : cell->vertex_indices())
          if (locally_active_vertices_on_level_subdomain[cell->vertex_index(
                v)] == true)
            {
              level_halo_layer.push_back(cell);
              break;
            }

    return level_halo_layer;
  }


  namespace
  {
    template <typename MeshType>
    DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
    bool contains_locally_owned_cells(
      const std::vector<typename MeshType::active_cell_iterator> &cells)
    {
      for (typename std::vector<
             typename MeshType::active_cell_iterator>::const_iterator it =
             cells.begin();
           it != cells.end();
           ++it)
        {
          if ((*it)->is_locally_owned())
            return true;
        }
      return false;
    }

    template <typename MeshType>
    DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
    bool contains_artificial_cells(
      const std::vector<typename MeshType::active_cell_iterator> &cells)
    {
      for (typename std::vector<
             typename MeshType::active_cell_iterator>::const_iterator it =
             cells.begin();
           it != cells.end();
           ++it)
        {
          if ((*it)->is_artificial())
            return true;
        }
      return false;
    }
  } // namespace



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<
    typename MeshType::
      active_cell_iterator> compute_ghost_cell_halo_layer(const MeshType &mesh)
  {
    std::function<bool(const typename MeshType::active_cell_iterator &)>
      predicate = IteratorFilters::LocallyOwnedCell();

    const std::vector<typename MeshType::active_cell_iterator>
      active_halo_layer = compute_active_cell_halo_layer(mesh, predicate);

    // Check that we never return locally owned or artificial cells
    // What is left should only be the ghost cells
    Assert(contains_locally_owned_cells<MeshType>(active_halo_layer) == false,
           ExcMessage("Halo layer contains locally owned cells"));
    Assert(contains_artificial_cells<MeshType>(active_halo_layer) == false,
           ExcMessage("Halo layer contains artificial cells"));

    return active_halo_layer;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::active_cell_iterator> compute_active_cell_layer_within_distance(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::active_cell_iterator &)>
                  &predicate,
      const double layer_thickness)
  {
    std::vector<typename MeshType::active_cell_iterator>
      subdomain_boundary_cells, active_cell_layer_within_distance;
    std::vector<bool> vertices_outside_subdomain(
      mesh.get_triangulation().n_vertices(), false);

    const unsigned int spacedim = MeshType::space_dimension;

    unsigned int n_non_predicate_cells = 0; // Number of non predicate cells

    // Find the layer of cells for which predicate is true and that
    // are on the boundary with other cells. These are
    // subdomain boundary cells.

    // Find the cells for which the predicate is false
    // These are the cells which are around the predicate subdomain
    for (const auto &cell : mesh.active_cell_iterators())
      if (!predicate(cell)) // Negation of predicate --> Not Part of subdomain
        {
          for (const unsigned int v : cell->vertex_indices())
            vertices_outside_subdomain[cell->vertex_index(v)] = true;
          ++n_non_predicate_cells;
        }

    // If all the active cells conform to the predicate
    // or if none of the active cells conform to the predicate
    // there is no active cell layer around the predicate
    // subdomain (within any distance)
    if (n_non_predicate_cells == 0 ||
        n_non_predicate_cells == mesh.get_triangulation().n_active_cells())
      return std::vector<typename MeshType::active_cell_iterator>();

    // Find the cells that conform to the predicate
    // but share a vertex with the cell not in the predicate subdomain
    for (const auto &cell : mesh.active_cell_iterators())
      if (predicate(cell)) // True predicate --> Potential boundary cell of the
                           // subdomain
        for (const unsigned int v : cell->vertex_indices())
          if (vertices_outside_subdomain[cell->vertex_index(v)] == true)
            {
              subdomain_boundary_cells.push_back(cell);
              break; // No need to go through remaining vertices
            }

    // To cheaply filter out some cells located far away from the predicate
    // subdomain, get the bounding box of the predicate subdomain.
    std::pair<Point<spacedim>, Point<spacedim>> bounding_box =
      compute_bounding_box(mesh, predicate);

    // DOUBLE_EPSILON to compare really close double values
    const double DOUBLE_EPSILON = 100. * std::numeric_limits<double>::epsilon();

    // Add layer_thickness to the bounding box
    for (unsigned int d = 0; d < spacedim; ++d)
      {
        bounding_box.first[d] -= (layer_thickness + DOUBLE_EPSILON);  // minp
        bounding_box.second[d] += (layer_thickness + DOUBLE_EPSILON); // maxp
      }

    std::vector<Point<spacedim>>
      subdomain_boundary_cells_centers; // cache all the subdomain boundary
                                        // cells centers here
    std::vector<double>
      subdomain_boundary_cells_radii; // cache all the subdomain boundary cells
                                      // radii
    subdomain_boundary_cells_centers.reserve(subdomain_boundary_cells.size());
    subdomain_boundary_cells_radii.reserve(subdomain_boundary_cells.size());
    // compute cell radius for each boundary cell of the predicate subdomain
    for (typename std::vector<typename MeshType::active_cell_iterator>::
           const_iterator subdomain_boundary_cell_iterator =
             subdomain_boundary_cells.begin();
         subdomain_boundary_cell_iterator != subdomain_boundary_cells.end();
         ++subdomain_boundary_cell_iterator)
      {
        const std::pair<Point<spacedim>, double>
          &subdomain_boundary_cell_enclosing_ball =
            (*subdomain_boundary_cell_iterator)->enclosing_ball();

        subdomain_boundary_cells_centers.push_back(
          subdomain_boundary_cell_enclosing_ball.first);
        subdomain_boundary_cells_radii.push_back(
          subdomain_boundary_cell_enclosing_ball.second);
      }
    AssertThrow(subdomain_boundary_cells_radii.size() ==
                  subdomain_boundary_cells_centers.size(),
                ExcInternalError());

    // Find the cells that are within layer_thickness of predicate subdomain
    // boundary distance but are inside the extended bounding box. Most cells
    // might be outside the extended bounding box, so we could skip them. Those
    // cells that are inside the extended bounding box but are not part of the
    // predicate subdomain are possible candidates to be within the distance to
    // the boundary cells of the predicate subdomain.
    for (const auto &cell : mesh.active_cell_iterators())
      {
        // Ignore all the cells that are in the predicate subdomain
        if (predicate(cell))
          continue;

        const std::pair<Point<spacedim>, double> &cell_enclosing_ball =
          cell->enclosing_ball();

        const Point<spacedim> cell_enclosing_ball_center =
          cell_enclosing_ball.first;
        const double cell_enclosing_ball_radius = cell_enclosing_ball.second;

        bool cell_inside = true; // reset for each cell

        for (unsigned int d = 0; d < spacedim; ++d)
          cell_inside &=
            (cell_enclosing_ball_center[d] + cell_enclosing_ball_radius >
             bounding_box.first[d]) &&
            (cell_enclosing_ball_center[d] - cell_enclosing_ball_radius <
             bounding_box.second[d]);
        // cell_inside is true if its enclosing ball intersects the extended
        // bounding box

        // Ignore all the cells that are outside the extended bounding box
        if (cell_inside)
          for (unsigned int i = 0; i < subdomain_boundary_cells_radii.size();
               ++i)
            if (cell_enclosing_ball_center.distance_square(
                  subdomain_boundary_cells_centers[i]) <
                Utilities::fixed_power<2>(cell_enclosing_ball_radius +
                                          subdomain_boundary_cells_radii[i] +
                                          layer_thickness + DOUBLE_EPSILON))
              {
                active_cell_layer_within_distance.push_back(cell);
                break; // Exit the loop checking all the remaining subdomain
                       // boundary cells
              }
      }
    return active_cell_layer_within_distance;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<
    typename MeshType::
      active_cell_iterator> compute_ghost_cell_layer_within_distance(const MeshType
                                                                       &mesh,
                                                                     const double
                                                                       layer_thickness)
  {
    IteratorFilters::LocallyOwnedCell locally_owned_cell_predicate;
    std::function<bool(const typename MeshType::active_cell_iterator &)>
      predicate(locally_owned_cell_predicate);

    const std::vector<typename MeshType::active_cell_iterator>
      ghost_cell_layer_within_distance =
        compute_active_cell_layer_within_distance(mesh,
                                                  predicate,
                                                  layer_thickness);

    // Check that we never return locally owned or artificial cells
    // What is left should only be the ghost cells
    Assert(
      contains_locally_owned_cells<MeshType>(
        ghost_cell_layer_within_distance) == false,
      ExcMessage(
        "Ghost cells within layer_thickness contains locally owned cells."));
    Assert(
      contains_artificial_cells<MeshType>(ghost_cell_layer_within_distance) ==
        false,
      ExcMessage(
        "Ghost cells within layer_thickness contains artificial cells. "
        "The function compute_ghost_cell_layer_within_distance "
        "is probably called while using parallel::distributed::Triangulation. "
        "In such case please refer to the description of this function."));

    return ghost_cell_layer_within_distance;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::pair<
    Point<MeshType::space_dimension>,
    Point<MeshType::
            space_dimension>> compute_bounding_box(const MeshType &mesh,
                                                   const std::function<bool(
                                                     const typename MeshType::
                                                       active_cell_iterator &)>
                                                     &predicate)
  {
    std::vector<bool> locally_active_vertices_on_subdomain(
      mesh.get_triangulation().n_vertices(), false);

    const unsigned int spacedim = MeshType::space_dimension;

    // Two extreme points can define the bounding box
    // around the active cells that conform to the given predicate.
    Point<MeshType::space_dimension> maxp, minp;

    // initialize minp and maxp with the first predicate cell center
    for (const auto &cell : mesh.active_cell_iterators())
      if (predicate(cell))
        {
          minp = cell->center();
          maxp = cell->center();
          break;
        }

    // Run through all the cells to check if it belongs to predicate domain,
    // if it belongs to the predicate domain, extend the bounding box.
    for (const auto &cell : mesh.active_cell_iterators())
      if (predicate(cell)) // True predicate --> Part of subdomain
        for (const unsigned int v : cell->vertex_indices())
          if (locally_active_vertices_on_subdomain[cell->vertex_index(v)] ==
              false)
            {
              locally_active_vertices_on_subdomain[cell->vertex_index(v)] =
                true;
              for (unsigned int d = 0; d < spacedim; ++d)
                {
                  minp[d] = std::min(minp[d], cell->vertex(v)[d]);
                  maxp[d] = std::max(maxp[d], cell->vertex(v)[d]);
                }
            }

    return std::make_pair(minp, maxp);
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::list<std::pair<
    typename MeshType::cell_iterator,
    typename MeshType::cell_iterator>> get_finest_common_cells(const MeshType
                                                                 &mesh_1,
                                                               const MeshType
                                                                 &mesh_2)
  {
    Assert(have_same_coarse_mesh(mesh_1, mesh_2),
           ExcMessage("The two meshes must be represent triangulations that "
                      "have the same coarse meshes"));
    // We will allow the output to contain ghost cells when we have shared
    // Triangulations (i.e., so that each processor will get exactly the same
    // list of cell pairs), but not when we have two distributed
    // Triangulations (so that all active cells are partitioned by processor).
    // Non-parallel Triangulations have no ghost or artificial cells, so they
    // work the same way as shared Triangulations here.
    bool remove_ghost_cells = false;
#ifdef DEAL_II_WITH_MPI
    {
      constexpr int dim      = MeshType::dimension;
      constexpr int spacedim = MeshType::space_dimension;
      if (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>
                         *>(&mesh_1.get_triangulation()) != nullptr ||
          dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>
                         *>(&mesh_2.get_triangulation()) != nullptr)
        {
          Assert(&mesh_1.get_triangulation() == &mesh_2.get_triangulation(),
                 ExcMessage("This function can only be used with meshes "
                            "corresponding to distributed Triangulations when "
                            "both Triangulations are equal."));
          remove_ghost_cells = true;
        }
    }
#endif

    // the algorithm goes as follows: first, we fill a list with pairs of
    // iterators common to the two meshes on the coarsest level. then we
    // traverse the list; each time, we find a pair of iterators for which
    // both correspond to non-active cells, we delete this item and push the
    // pairs of iterators to their children to the back. if these again both
    // correspond to non-active cells, we will get to the later on for further
    // consideration
    using CellList = std::list<std::pair<typename MeshType::cell_iterator,
                                         typename MeshType::cell_iterator>>;
    CellList cell_list;

    // first push the coarse level cells
    typename MeshType::cell_iterator cell_1 = mesh_1.begin(0),
                                     cell_2 = mesh_2.begin(0);
    for (; cell_1 != mesh_1.end(0); ++cell_1, ++cell_2)
      cell_list.emplace_back(cell_1, cell_2);

    // then traverse list as described above
    typename CellList::iterator cell_pair = cell_list.begin();
    while (cell_pair != cell_list.end())
      {
        // if both cells in this pair have children, then erase this element
        // and push their children instead
        if (cell_pair->first->has_children() &&
            cell_pair->second->has_children())
          {
            Assert(cell_pair->first->refinement_case() ==
                     cell_pair->second->refinement_case(),
                   ExcNotImplemented());
            for (unsigned int c = 0; c < cell_pair->first->n_children(); ++c)
              cell_list.emplace_back(cell_pair->first->child(c),
                                     cell_pair->second->child(c));

            // erasing an iterator keeps other iterators valid, so already
            // advance the present iterator by one and then delete the element
            // we've visited before
            const auto previous_cell_pair = cell_pair;
            ++cell_pair;
            cell_list.erase(previous_cell_pair);
          }
        else
          {
            // at least one cell is active
            if (remove_ghost_cells &&
                ((cell_pair->first->is_active() &&
                  !cell_pair->first->is_locally_owned()) ||
                 (cell_pair->second->is_active() &&
                  !cell_pair->second->is_locally_owned())))
              {
                // we only exclude ghost cells for distributed Triangulations
                const auto previous_cell_pair = cell_pair;
                ++cell_pair;
                cell_list.erase(previous_cell_pair);
              }
            else
              ++cell_pair;
          }
      }

    // just to make sure everything is ok, validate that all pairs have at
    // least one active iterator or have different refinement_cases
    for (cell_pair = cell_list.begin(); cell_pair != cell_list.end();
         ++cell_pair)
      Assert(cell_pair->first->is_active() || cell_pair->second->is_active() ||
               (cell_pair->first->refinement_case() !=
                cell_pair->second->refinement_case()),
             ExcInternalError());

    return cell_list;
  }



  template <int dim, int spacedim>
  bool
  have_same_coarse_mesh(const Triangulation<dim, spacedim> &mesh_1,
                        const Triangulation<dim, spacedim> &mesh_2)
  {
    // make sure the two meshes have
    // the same number of coarse cells
    if (mesh_1.n_cells(0) != mesh_2.n_cells(0))
      return false;

    // if so, also make sure they have
    // the same vertices on the cells
    // of the coarse mesh
    typename Triangulation<dim, spacedim>::cell_iterator cell_1 =
                                                           mesh_1.begin(0),
                                                         cell_2 =
                                                           mesh_2.begin(0),
                                                         endc = mesh_1.end(0);
    for (; cell_1 != endc; ++cell_1, ++cell_2)
      {
        if (cell_1->n_vertices() != cell_2->n_vertices())
          return false;
        for (const unsigned int v : cell_1->vertex_indices())
          if (cell_1->vertex(v) != cell_2->vertex(v))
            return false;
      }

    // if we've gotten through all
    // this, then the meshes really
    // seem to have a common coarse
    // mesh
    return true;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  bool have_same_coarse_mesh(const MeshType &mesh_1, const MeshType &mesh_2)
  {
    return have_same_coarse_mesh(mesh_1.get_triangulation(),
                                 mesh_2.get_triangulation());
  }



  template <int dim, int spacedim>
  std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
            Point<dim>>
  find_active_cell_around_point(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &mesh,
    const Point<spacedim>                      &p,
    const double                                tolerance)
  {
    Assert((mapping.size() == 1) ||
             (mapping.size() == mesh.get_fe_collection().size()),
           ExcMessage("Mapping collection needs to have either size 1 "
                      "or size equal to the number of elements in "
                      "the FECollection."));

    using cell_iterator =
      typename DoFHandler<dim, spacedim>::active_cell_iterator;

    std::pair<cell_iterator, Point<dim>> best_cell;
    // If we have only one element in the MappingCollection,
    // we use find_active_cell_around_point using only one
    // mapping.
    if (mapping.size() == 1)
      {
        const std::vector<bool> marked_vertices = {};
        best_cell                               = find_active_cell_around_point(
          mapping[0], mesh, p, marked_vertices, tolerance);
      }
    else
      {
        // The best distance is set to the
        // maximum allowable distance from
        // the unit cell
        double best_distance = tolerance;
        int    best_level    = -1;


        // Find closest vertex and determine
        // all adjacent cells
        unsigned int vertex = find_closest_vertex(mesh, p);

        std::vector<cell_iterator> adjacent_cells_tmp =
          find_cells_adjacent_to_vertex(mesh, vertex);

        // Make sure that we have found
        // at least one cell adjacent to vertex.
        Assert(adjacent_cells_tmp.size() > 0, ExcInternalError());

        // Copy all the cells into a std::set
        std::set<cell_iterator> adjacent_cells(adjacent_cells_tmp.begin(),
                                               adjacent_cells_tmp.end());
        std::set<cell_iterator> searched_cells;

        // Determine the maximal number of cells
        // in the grid.
        // As long as we have not found
        // the cell and have not searched
        // every cell in the triangulation,
        // we keep on looking.
        const auto   n_cells        = mesh.get_triangulation().n_cells();
        bool         found          = false;
        unsigned int cells_searched = 0;
        while (!found && cells_searched < n_cells)
          {
            for (const auto &cell : adjacent_cells)
              {
                try
                  {
                    const Point<dim> p_cell =
                      mapping[cell->active_fe_index()]
                        .transform_real_to_unit_cell(cell, p);


                    // calculate the Euclidean norm of
                    // the distance vector to the unit cell.
                    const double dist =
                      cell->reference_cell().closest_point(p_cell).distance(
                        p_cell);

                    // We compare if the point is inside the
                    // unit cell (or at least not too far
                    // outside). If it is, it is also checked
                    // that the cell has a more refined state
                    if (dist < best_distance ||
                        (dist == best_distance && cell->level() > best_level))
                      {
                        found         = true;
                        best_distance = dist;
                        best_level    = cell->level();
                        best_cell     = std::make_pair(cell, p_cell);
                      }
                  }
                catch (
                  typename MappingQ<dim, spacedim>::ExcTransformationFailed &)
                  {
                    // ok, the transformation
                    // failed presumably
                    // because the point we
                    // are looking for lies
                    // outside the current
                    // cell. this means that
                    // the current cell can't
                    // be the cell around the
                    // point, so just ignore
                    // this cell and move on
                    // to the next
                  }
              }
            // update the number of cells searched
            cells_searched += adjacent_cells.size();
            // if we have not found the cell in
            // question and have not yet searched every
            // cell, we expand our search to
            // all the not already searched neighbors of
            // the cells in adjacent_cells.
            if (!found && cells_searched < n_cells)
              {
                find_active_cell_around_point_internal<dim,
                                                       DoFHandler,
                                                       spacedim>(
                  mesh, searched_cells, adjacent_cells);
              }
          }
      }

    return best_cell;
  }


  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<typename MeshType::active_cell_iterator> get_patch_around_cell(
    const typename MeshType::active_cell_iterator &cell)
  {
    Assert(cell->is_locally_owned(),
           ExcMessage("This function only makes sense if the cell for "
                      "which you are asking for a patch, is locally "
                      "owned."));

    std::vector<typename MeshType::active_cell_iterator> patch;
    patch.push_back(cell);
    for (const unsigned int face_number : cell->face_indices())
      if (cell->face(face_number)->at_boundary() == false)
        {
          if (cell->neighbor(face_number)->has_children() == false)
            patch.push_back(cell->neighbor(face_number));
          else
            // the neighbor is refined. in 2d/3d, we can simply ask for the
            // children of the neighbor because they can not be further refined
            // and, consequently, the children is active
            if (MeshType::dimension > 1)
              {
                for (unsigned int subface = 0;
                     subface < cell->face(face_number)->n_children();
                     ++subface)
                  patch.push_back(
                    cell->neighbor_child_on_subface(face_number, subface));
              }
            else
              {
                // in 1d, we need to work a bit harder: iterate until we find
                // the child by going from cell to child to child etc
                typename MeshType::cell_iterator neighbor =
                  cell->neighbor(face_number);
                while (neighbor->has_children())
                  neighbor = neighbor->child(1 - face_number);

                Assert(neighbor->neighbor(1 - face_number) == cell,
                       ExcInternalError());
                patch.push_back(neighbor);
              }
        }
    return patch;
  }



  template <class Container>
  std::vector<typename Container::cell_iterator>
  get_cells_at_coarsest_common_level(
    const std::vector<typename Container::active_cell_iterator> &patch)
  {
    Assert(patch.size() > 0,
           ExcMessage(
             "Vector containing patch cells should not be an empty vector!"));
    // In order to extract the set of cells with the coarsest common level from
    // the give vector of cells: First it finds the number associated with the
    // minimum level of refinement, namely "min_level"
    int min_level = patch[0]->level();

    for (unsigned int i = 0; i < patch.size(); ++i)
      min_level = std::min(min_level, patch[i]->level());
    std::set<typename Container::cell_iterator> uniform_cells;
    typename std::vector<
      typename Container::active_cell_iterator>::const_iterator patch_cell;
    // it loops through all cells of the input vector
    for (patch_cell = patch.begin(); patch_cell != patch.end(); ++patch_cell)
      {
        // If the refinement level of each cell i the loop be equal to the
        // min_level, so that that cell inserted into the set of uniform_cells,
        // as the set of cells with the coarsest common refinement level
        if ((*patch_cell)->level() == min_level)
          uniform_cells.insert(*patch_cell);
        else
          // If not, it asks for the parent of the cell, until it finds the
          // parent cell with the refinement level equal to the min_level and
          // inserts that parent cell into the set of uniform_cells, as the
          // set of cells with the coarsest common refinement level.
          {
            typename Container::cell_iterator parent = *patch_cell;

            while (parent->level() > min_level)
              parent = parent->parent();
            uniform_cells.insert(parent);
          }
      }

    return std::vector<typename Container::cell_iterator>(uniform_cells.begin(),
                                                          uniform_cells.end());
  }



  template <class Container>
  void
  build_triangulation_from_patch(
    const std::vector<typename Container::active_cell_iterator> &patch,
    Triangulation<Container::dimension, Container::space_dimension>
      &local_triangulation,
    std::map<
      typename Triangulation<Container::dimension,
                             Container::space_dimension>::active_cell_iterator,
      typename Container::active_cell_iterator> &patch_to_global_tria_map)

  {
    const std::vector<typename Container::cell_iterator> uniform_cells =
      get_cells_at_coarsest_common_level<Container>(patch);
    // First it creates triangulation from the vector of "uniform_cells"
    local_triangulation.clear();
    std::vector<Point<Container::space_dimension>> vertices;
    const unsigned int n_uniform_cells = uniform_cells.size();
    std::vector<CellData<Container::dimension>> cells(n_uniform_cells);
    unsigned int                                k = 0; // for enumerating cells
    unsigned int i = 0; // for enumerating vertices
    typename std::vector<typename Container::cell_iterator>::const_iterator
      uniform_cell;
    for (uniform_cell = uniform_cells.begin();
         uniform_cell != uniform_cells.end();
         ++uniform_cell)
      {
        for (const unsigned int v : (*uniform_cell)->vertex_indices())
          {
            Point<Container::space_dimension> position =
              (*uniform_cell)->vertex(v);
            bool repeat_vertex = false;

            for (unsigned int m = 0; m < i; ++m)
              {
                if (position == vertices[m])
                  {
                    repeat_vertex        = true;
                    cells[k].vertices[v] = m;
                    break;
                  }
              }
            if (repeat_vertex == false)
              {
                vertices.push_back(position);
                cells[k].vertices[v] = i;
                i                    = i + 1;
              }

          } // for vertices_per_cell
        k = k + 1;
      }
    local_triangulation.create_triangulation(vertices, cells, SubCellData());
    Assert(local_triangulation.n_active_cells() == uniform_cells.size(),
           ExcInternalError());
    local_triangulation.clear_user_flags();
    unsigned int index = 0;
    // Create a map between cells of class DoFHandler into class Triangulation
    std::map<typename Triangulation<Container::dimension,
                                    Container::space_dimension>::cell_iterator,
             typename Container::cell_iterator>
      patch_to_global_tria_map_tmp;
    for (typename Triangulation<Container::dimension,
                                Container::space_dimension>::cell_iterator
           coarse_cell = local_triangulation.begin();
         coarse_cell != local_triangulation.end();
         ++coarse_cell, ++index)
      {
        patch_to_global_tria_map_tmp.insert(
          std::make_pair(coarse_cell, uniform_cells[index]));
        // To ensure that the cells with the same coordinates (here, we compare
        // their centers) are mapped into each other.

        Assert(coarse_cell->center().distance(uniform_cells[index]->center()) <=
                 1e-15 * coarse_cell->diameter(),
               ExcInternalError());
      }
    bool refinement_necessary;
    // In this loop we start to do refinement on the above coarse triangulation
    // to reach to the same level of refinement as the patch cells are really on
    do
      {
        refinement_necessary = false;
        for (const auto &active_tria_cell :
             local_triangulation.active_cell_iterators())
          {
            if (patch_to_global_tria_map_tmp[active_tria_cell]->has_children())
              {
                active_tria_cell->set_refine_flag();
                refinement_necessary = true;
              }
            else
              for (unsigned int i = 0; i < patch.size(); ++i)
                {
                  // Even though vertices may not be exactly the same, the
                  // appropriate cells will match since == for TriAccessors
                  // checks only cell level and index.
                  if (patch_to_global_tria_map_tmp[active_tria_cell] ==
                      patch[i])
                    {
                      // adjust the cell vertices of the local_triangulation to
                      // match cell vertices of the global triangulation
                      for (const unsigned int v :
                           active_tria_cell->vertex_indices())
                        active_tria_cell->vertex(v) = patch[i]->vertex(v);

                      Assert(active_tria_cell->center().distance(
                               patch_to_global_tria_map_tmp[active_tria_cell]
                                 ->center()) <=
                               1e-15 * active_tria_cell->diameter(),
                             ExcInternalError());

                      active_tria_cell->set_user_flag();
                      break;
                    }
                }
          }

        if (refinement_necessary)
          {
            local_triangulation.execute_coarsening_and_refinement();

            for (typename Triangulation<
                   Container::dimension,
                   Container::space_dimension>::cell_iterator cell =
                   local_triangulation.begin();
                 cell != local_triangulation.end();
                 ++cell)
              {
                if (patch_to_global_tria_map_tmp.find(cell) !=
                    patch_to_global_tria_map_tmp.end())
                  {
                    if (cell->has_children())
                      {
                        // Note: Since the cell got children, then it should not
                        // be in the map anymore children may be added into the
                        // map, instead

                        // these children may not yet be in the map
                        for (unsigned int c = 0; c < cell->n_children(); ++c)
                          {
                            if (patch_to_global_tria_map_tmp.find(cell->child(
                                  c)) == patch_to_global_tria_map_tmp.end())
                              {
                                patch_to_global_tria_map_tmp.insert(
                                  std::make_pair(
                                    cell->child(c),
                                    patch_to_global_tria_map_tmp[cell]->child(
                                      c)));

                                // One might be tempted to assert that the cell
                                // being added here has the same center as the
                                // equivalent cell in the global triangulation,
                                // but it may not be the case.  For
                                // triangulations that have been perturbed or
                                // smoothed, the cell indices and levels may be
                                // the same, but the vertex locations may not.
                                // We adjust the vertices of the cells that have
                                // no children (ie the active cells) to be
                                // consistent with the global triangulation
                                // later on and add assertions at that time
                                // to guarantee the cells in the
                                // local_triangulation are physically at the
                                // same locations of the cells in the patch of
                                // the global triangulation.
                              }
                          }
                        // The parent cell whose children were added
                        // into the map should be deleted from the map
                        patch_to_global_tria_map_tmp.erase(cell);
                      }
                  }
              }
          }
      }
    while (refinement_necessary);


    // Last assertion check to make sure we have the right cells and centers
    // in the map, and hence the correct vertices of the triangulation
    for (typename Triangulation<Container::dimension,
                                Container::space_dimension>::cell_iterator
           cell = local_triangulation.begin();
         cell != local_triangulation.end();
         ++cell)
      {
        if (cell->user_flag_set())
          {
            Assert(patch_to_global_tria_map_tmp.find(cell) !=
                     patch_to_global_tria_map_tmp.end(),
                   ExcInternalError());

            Assert(cell->center().distance(
                     patch_to_global_tria_map_tmp[cell]->center()) <=
                     1e-15 * cell->diameter(),
                   ExcInternalError());
          }
      }


    typename std::map<
      typename Triangulation<Container::dimension,
                             Container::space_dimension>::cell_iterator,
      typename Container::cell_iterator>::iterator
      map_tmp_it  = patch_to_global_tria_map_tmp.begin(),
      map_tmp_end = patch_to_global_tria_map_tmp.end();
    // Now we just need to take the temporary map of pairs of type cell_iterator
    // "patch_to_global_tria_map_tmp" making pair of active_cell_iterators so
    // that filling out the final map "patch_to_global_tria_map"
    for (; map_tmp_it != map_tmp_end; ++map_tmp_it)
      patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;
  }



  template <int dim, int spacedim>
  std::map<
    types::global_dof_index,
    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>>
  get_dof_to_support_patch_map(DoFHandler<dim, spacedim> &dof_handler)
  {
    // This is the map from global_dof_index to
    // a set of cells on patch.  We first map into
    // a set because it is very likely that we
    // will attempt to add a cell more than once
    // to a particular patch and we want to preserve
    // uniqueness of cell iterators. std::set does this
    // automatically for us.  Later after it is all
    // constructed, we will copy to a map of vectors
    // since that is the preferred output for other
    // functions.
    std::map<types::global_dof_index,
             std::set<typename DoFHandler<dim, spacedim>::active_cell_iterator>>
      dof_to_set_of_cells_map;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> local_face_dof_indices;
    std::vector<types::global_dof_index> local_line_dof_indices;

    // a place to save the dof_handler user flags and restore them later
    // to maintain const of dof_handler.
    std::vector<bool> user_flags;


    // in 3d, we need pointers from active lines to the
    // active parent lines, so we construct it as needed.
    std::map<typename DoFHandler<dim, spacedim>::active_line_iterator,
             typename DoFHandler<dim, spacedim>::line_iterator>
      lines_to_parent_lines_map;
    if (dim == 3)
      {
        // save user flags as they will be modified and then later restored
        dof_handler.get_triangulation().save_user_flags(user_flags);
        const_cast<dealii::Triangulation<dim, spacedim> &>(
          dof_handler.get_triangulation())
          .clear_user_flags();


        typename DoFHandler<dim, spacedim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
        for (; cell != endc; ++cell)
          {
            // We only want lines that are locally_relevant
            // although it doesn't hurt to have lines that
            // are children of ghost cells since there are
            // few and we don't have to use them.
            if (cell->is_artificial() == false)
              {
                for (unsigned int l = 0; l < cell->n_lines(); ++l)
                  if (cell->line(l)->has_children())
                    for (unsigned int c = 0; c < cell->line(l)->n_children();
                         ++c)
                      {
                        lines_to_parent_lines_map[cell->line(l)->child(c)] =
                          cell->line(l);
                        // set flags to know that child
                        // line has an active parent.
                        cell->line(l)->child(c)->set_user_flag();
                      }
              }
          }
      }


    // We loop through all cells and add cell to the
    // map for the dofs that it immediately touches
    // and then account for all the other dofs of
    // which it is a part, mainly the ones that must
    // be added on account of adaptivity hanging node
    // constraints.
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        // Need to loop through all cells that could
        // be in the patch of dofs on locally_owned
        // cells including ghost cells
        if (cell->is_artificial() == false)
          {
            const unsigned int n_dofs_per_cell =
              cell->get_fe().n_dofs_per_cell();
            local_dof_indices.resize(n_dofs_per_cell);

            // Take care of adding cell pointer to each
            // dofs that exists on cell.
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              dof_to_set_of_cells_map[local_dof_indices[i]].insert(cell);

            // In the case of the adjacent cell (over
            // faces or edges) being more refined, we
            // want to add all of the children to the
            // patch since the support function at that
            // dof could be non-zero along that entire
            // face (or line).

            // Take care of dofs on neighbor faces
            for (const unsigned int f : cell->face_indices())
              {
                if (cell->face(f)->has_children())
                  {
                    for (unsigned int c = 0; c < cell->face(f)->n_children();
                         ++c)
                      {
                        //  Add cell to dofs of all subfaces
                        //
                        //   *-------------------*----------*---------*
                        //   |                   | add cell |         |
                        //   |                   |<- to dofs|         |
                        //   |                   |of subface|         |
                        //   |        cell       *----------*---------*
                        //   |                   | add cell |         |
                        //   |                   |<- to dofs|         |
                        //   |                   |of subface|         |
                        //   *-------------------*----------*---------*
                        //
                        Assert(cell->face(f)->child(c)->has_children() == false,
                               ExcInternalError());

                        const unsigned int n_dofs_per_face =
                          cell->get_fe().n_dofs_per_face(f, c);
                        local_face_dof_indices.resize(n_dofs_per_face);

                        cell->face(f)->child(c)->get_dof_indices(
                          local_face_dof_indices);
                        for (unsigned int i = 0; i < n_dofs_per_face; ++i)
                          dof_to_set_of_cells_map[local_face_dof_indices[i]]
                            .insert(cell);
                      }
                  }
                else if ((cell->face(f)->at_boundary() == false) &&
                         (cell->neighbor_is_coarser(f)))
                  {
                    // Add cell to dofs of parent face and all
                    // child faces of parent face
                    //
                    //   *-------------------*----------*---------*
                    //   |                   |          |         |
                    //   |                   |   cell   |         |
                    //   |      add cell     |          |         |
                    //   |      to dofs   -> *----------*---------*
                    //   |      of parent    | add cell |         |
                    //   |       face        |<- to dofs|         |
                    //   |                   |of subface|         |
                    //   *-------------------*----------*---------*
                    //

                    // Add cell to all dofs of parent face
                    std::pair<unsigned int, unsigned int>
                      neighbor_face_no_subface_no =
                        cell->neighbor_of_coarser_neighbor(f);
                    unsigned int face_no = neighbor_face_no_subface_no.first;
                    unsigned int subface = neighbor_face_no_subface_no.second;

                    const unsigned int n_dofs_per_face =
                      cell->get_fe().n_dofs_per_face(face_no);
                    local_face_dof_indices.resize(n_dofs_per_face);

                    cell->neighbor(f)->face(face_no)->get_dof_indices(
                      local_face_dof_indices);
                    for (unsigned int i = 0; i < n_dofs_per_face; ++i)
                      dof_to_set_of_cells_map[local_face_dof_indices[i]].insert(
                        cell);

                    // Add cell to all dofs of children of
                    // parent face
                    for (unsigned int c = 0;
                         c < cell->neighbor(f)->face(face_no)->n_children();
                         ++c)
                      {
                        if (c != subface) // don't repeat work on dofs of
                                          // original cell
                          {
                            const unsigned int n_dofs_per_face =
                              cell->get_fe().n_dofs_per_face(face_no, c);
                            local_face_dof_indices.resize(n_dofs_per_face);

                            Assert(cell->neighbor(f)
                                       ->face(face_no)
                                       ->child(c)
                                       ->has_children() == false,
                                   ExcInternalError());
                            cell->neighbor(f)
                              ->face(face_no)
                              ->child(c)
                              ->get_dof_indices(local_face_dof_indices);
                            for (unsigned int i = 0; i < n_dofs_per_face; ++i)
                              dof_to_set_of_cells_map[local_face_dof_indices[i]]
                                .insert(cell);
                          }
                      }
                  }
              }


            // If 3d, take care of dofs on lines in the
            // same pattern as faces above. That is, if
            // a cell's line has children, distribute
            // cell to dofs of children of line,  and
            // if cell's line has an active parent, then
            // distribute cell to dofs on parent line
            // and dofs on all children of parent line.
            if (dim == 3)
              {
                for (unsigned int l = 0; l < cell->n_lines(); ++l)
                  {
                    if (cell->line(l)->has_children())
                      {
                        for (unsigned int c = 0;
                             c < cell->line(l)->n_children();
                             ++c)
                          {
                            Assert(cell->line(l)->child(c)->has_children() ==
                                     false,
                                   ExcInternalError());

                            // dofs_per_line returns number of dofs
                            // on line not including the vertices of the line.
                            const unsigned int n_dofs_per_line =
                              2 * cell->get_fe().n_dofs_per_vertex() +
                              cell->get_fe().n_dofs_per_line();
                            local_line_dof_indices.resize(n_dofs_per_line);

                            cell->line(l)->child(c)->get_dof_indices(
                              local_line_dof_indices);
                            for (unsigned int i = 0; i < n_dofs_per_line; ++i)
                              dof_to_set_of_cells_map[local_line_dof_indices[i]]
                                .insert(cell);
                          }
                      }
                    // user flag was set above to denote that
                    // an active parent line exists so add
                    // cell to dofs of parent and all it's
                    // children
                    else if (cell->line(l)->user_flag_set() == true)
                      {
                        typename DoFHandler<dim, spacedim>::line_iterator
                          parent_line =
                            lines_to_parent_lines_map[cell->line(l)];
                        Assert(parent_line->has_children(), ExcInternalError());

                        // dofs_per_line returns number of dofs
                        // on line not including the vertices of the line.
                        const unsigned int n_dofs_per_line =
                          2 * cell->get_fe().n_dofs_per_vertex() +
                          cell->get_fe().n_dofs_per_line();
                        local_line_dof_indices.resize(n_dofs_per_line);

                        parent_line->get_dof_indices(local_line_dof_indices);
                        for (unsigned int i = 0; i < n_dofs_per_line; ++i)
                          dof_to_set_of_cells_map[local_line_dof_indices[i]]
                            .insert(cell);

                        for (unsigned int c = 0; c < parent_line->n_children();
                             ++c)
                          {
                            Assert(parent_line->child(c)->has_children() ==
                                     false,
                                   ExcInternalError());

                            const unsigned int n_dofs_per_line =
                              2 * cell->get_fe().n_dofs_per_vertex() +
                              cell->get_fe().n_dofs_per_line();
                            local_line_dof_indices.resize(n_dofs_per_line);

                            parent_line->child(c)->get_dof_indices(
                              local_line_dof_indices);
                            for (unsigned int i = 0; i < n_dofs_per_line; ++i)
                              dof_to_set_of_cells_map[local_line_dof_indices[i]]
                                .insert(cell);
                          }
                      }
                  } // for lines l
              }     // if dim == 3
          }         // if cell->is_locally_owned()
      }             // for cells


    if (dim == 3)
      {
        // finally, restore user flags that were changed above
        // to when we constructed the pointers to parent of lines
        // Since dof_handler is const, we must leave it unchanged.
        const_cast<dealii::Triangulation<dim, spacedim> &>(
          dof_handler.get_triangulation())
          .load_user_flags(user_flags);
      }

    // Finally, we copy map of sets to
    // map of vectors using the std::vector::assign() function
    std::map<
      types::global_dof_index,
      std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>>
      dof_to_cell_patches;

    typename std::map<
      types::global_dof_index,
      std::set<typename DoFHandler<dim, spacedim>::active_cell_iterator>>::
      iterator it     = dof_to_set_of_cells_map.begin(),
               it_end = dof_to_set_of_cells_map.end();
    for (; it != it_end; ++it)
      dof_to_cell_patches[it->first].assign(it->second.begin(),
                                            it->second.end());

    return dof_to_cell_patches;
  }

  /*
   * Internally used in collect_periodic_faces
   */
  template <typename CellIterator>
  void
  match_periodic_face_pairs(
    std::set<std::pair<CellIterator, unsigned int>> &pairs1,
    std::set<std::pair<std_cxx20::type_identity_t<CellIterator>, unsigned int>>
                                                &pairs2,
    const unsigned int                           direction,
    std::vector<PeriodicFacePair<CellIterator>> &matched_pairs,
    const dealii::Tensor<1, CellIterator::AccessorType::space_dimension>
                             &offset,
    const FullMatrix<double> &matrix,
    const double              abs_tol = 1e-10)
  {
    static const int space_dim = CellIterator::AccessorType::space_dimension;
    AssertIndexRange(direction, space_dim);

    if constexpr (running_in_debug_mode())
      {
        {
          constexpr int dim      = CellIterator::AccessorType::dimension;
          constexpr int spacedim = CellIterator::AccessorType::space_dimension;
          // For parallel::fullydistributed::Triangulation there might be
          // unmatched faces on periodic boundaries on the coarse grid. As a
          // result this assert is not fulfilled (which is not a bug!). See also
          // the discussion in the method collect_periodic_faces.
          if (!(((pairs1.size() > 0) &&
                 (dynamic_cast<const parallel::fullydistributed::
                                 Triangulation<dim, spacedim> *>(
                    &pairs1.begin()->first->get_triangulation()) != nullptr)) ||
                ((pairs2.size() > 0) &&
                 (dynamic_cast<const parallel::fullydistributed::
                                 Triangulation<dim, spacedim> *>(
                    &pairs2.begin()->first->get_triangulation()) != nullptr))))
            Assert(pairs1.size() == pairs2.size(),
                   ExcMessage("Unmatched faces on periodic boundaries"));
        }
      }

    unsigned int n_matches = 0;

    // Match with a complexity of O(n^2). This could be improved...
    using PairIterator =
      typename std::set<std::pair<CellIterator, unsigned int>>::const_iterator;
    for (PairIterator it1 = pairs1.begin(); it1 != pairs1.end(); ++it1)
      {
        for (PairIterator it2 = pairs2.begin(); it2 != pairs2.end(); ++it2)
          {
            const CellIterator cell1     = it1->first;
            const CellIterator cell2     = it2->first;
            const unsigned int face_idx1 = it1->second;
            const unsigned int face_idx2 = it2->second;
            if (const std::optional<types::geometric_orientation> orientation =
                  GridTools::orthogonal_equality(cell1->face(face_idx1),
                                                 cell2->face(face_idx2),
                                                 direction,
                                                 offset,
                                                 matrix,
                                                 abs_tol))
              {
                // We have a match, so insert the matching pairs and
                // remove the matched cell in pairs2 to speed up the
                // matching:
                const PeriodicFacePair<CellIterator> matched_face = {
                  {cell1, cell2},
                  {face_idx1, face_idx2},
                  orientation.value(),
                  matrix};
                matched_pairs.push_back(matched_face);
                pairs2.erase(it2);
                ++n_matches;
                break;
              }
          }
      }

    // Assure that all faces are matched if
    // parallel::fullydistributed::Triangulation is not used. This is related to
    // the fact that the faces might not be successfully matched on the coarse
    // grid (not a bug!). See also the comment above and in the method
    // collect_periodic_faces.
    {
      constexpr int dim      = CellIterator::AccessorType::dimension;
      constexpr int spacedim = CellIterator::AccessorType::space_dimension;
      if (!(((pairs1.size() > 0) &&
             (dynamic_cast<const parallel::fullydistributed::
                             Triangulation<dim, spacedim> *>(
                &pairs1.begin()->first->get_triangulation()) != nullptr)) ||
            ((pairs2.size() > 0) &&
             (dynamic_cast<
                const parallel::fullydistributed::Triangulation<dim, spacedim>
                  *>(&pairs2.begin()->first->get_triangulation()) != nullptr))))
        AssertThrow(n_matches == pairs1.size() && pairs2.empty(),
                    ExcMessage("Unmatched faces on periodic boundaries"));
    }
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void collect_periodic_faces(
    const MeshType          &mesh,
    const types::boundary_id b_id,
    const unsigned int       direction,
    std::vector<PeriodicFacePair<typename MeshType::cell_iterator>>
                                               &matched_pairs,
    const Tensor<1, MeshType::space_dimension> &offset,
    const FullMatrix<double>                   &matrix,
    const double                                abs_tol)
  {
    static const int dim       = MeshType::dimension;
    static const int space_dim = MeshType::space_dimension;
    AssertIndexRange(direction, space_dim);
    Assert(dim == space_dim, ExcNotImplemented());

    // Loop over all cells on the highest level and collect all boundary
    // faces 2*direction and 2*direction*1:

    std::set<std::pair<typename MeshType::cell_iterator, unsigned int>> pairs1;
    std::set<std::pair<typename MeshType::cell_iterator, unsigned int>> pairs2;

    for (typename MeshType::cell_iterator cell = mesh.begin(0);
         cell != mesh.end(0);
         ++cell)
      {
        const typename MeshType::face_iterator face_1 =
          cell->face(2 * direction);
        const typename MeshType::face_iterator face_2 =
          cell->face(2 * direction + 1);

        if (face_1->at_boundary() && face_1->boundary_id() == b_id)
          {
            const std::pair<typename MeshType::cell_iterator, unsigned int>
              pair1 = std::make_pair(cell, 2 * direction);
            pairs1.insert(pair1);
          }

        if (face_2->at_boundary() && face_2->boundary_id() == b_id)
          {
            const std::pair<typename MeshType::cell_iterator, unsigned int>
              pair2 = std::make_pair(cell, 2 * direction + 1);
            pairs2.insert(pair2);
          }
      }

    Assert(pairs1.size() == pairs2.size(),
           ExcMessage("Unmatched faces on periodic boundaries"));

    Assert(pairs1.size() > 0,
           ExcMessage("No new periodic face pairs have been found. "
                      "Are you sure that you've selected the correct boundary "
                      "id's and that the coarsest level mesh is colorized?"));

    [[maybe_unused]] const unsigned int size_old = matched_pairs.size();

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(
      pairs1, pairs2, direction, matched_pairs, offset, matrix, abs_tol);

    if constexpr (running_in_debug_mode())
      {
        // check for standard orientation
        const unsigned int size_new = matched_pairs.size();
        for (unsigned int i = size_old; i < size_new; ++i)
          {
            Assert(matched_pairs[i].orientation ==
                     numbers::default_geometric_orientation,
                   ExcMessage(
                     "Found a face match with non standard orientation. "
                     "This function is only suitable for meshes with cells "
                     "in default orientation"));
          }
      }
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void collect_periodic_faces(
    const MeshType          &mesh,
    const types::boundary_id b_id1,
    const types::boundary_id b_id2,
    const unsigned int       direction,
    std::vector<PeriodicFacePair<typename MeshType::cell_iterator>>
                                               &matched_pairs,
    const Tensor<1, MeshType::space_dimension> &offset,
    const FullMatrix<double>                   &matrix,
    const double                                abs_tol)
  {
    static const int dim       = MeshType::dimension;
    static const int space_dim = MeshType::space_dimension;
    AssertIndexRange(direction, space_dim);

    // Loop over all cells on the highest level and collect all boundary
    // faces belonging to b_id1 and b_id2:

    std::set<std::pair<typename MeshType::cell_iterator, unsigned int>> pairs1;
    std::set<std::pair<typename MeshType::cell_iterator, unsigned int>> pairs2;

    for (typename MeshType::cell_iterator cell = mesh.begin(0);
         cell != mesh.end(0);
         ++cell)
      {
        for (const unsigned int i : cell->face_indices())
          {
            const typename MeshType::face_iterator face = cell->face(i);
            if (face->at_boundary() && face->boundary_id() == b_id1)
              {
                const std::pair<typename MeshType::cell_iterator, unsigned int>
                  pair1 = std::make_pair(cell, i);
                pairs1.insert(pair1);
              }

            if (face->at_boundary() && face->boundary_id() == b_id2)
              {
                const std::pair<typename MeshType::cell_iterator, unsigned int>
                  pair2 = std::make_pair(cell, i);
                pairs2.insert(pair2);
              }
          }
      }

    // Assure that all faces are matched on the coarse grid. This requirement
    // can only fulfilled if a process owns the complete coarse grid. This is
    // not the case for a parallel::fullydistributed::Triangulation, i.e., this
    // requirement has not to be met (consider faces on the outside of a
    // ghost cell that are periodic but for which the ghost neighbor doesn't
    // exist).
    if (!(((pairs1.size() > 0) &&
           (dynamic_cast<
              const parallel::fullydistributed::Triangulation<dim, space_dim>
                *>(&pairs1.begin()->first->get_triangulation()) != nullptr)) ||
          ((pairs2.size() > 0) &&
           (dynamic_cast<
              const parallel::fullydistributed::Triangulation<dim, space_dim>
                *>(&pairs2.begin()->first->get_triangulation()) != nullptr))))
      Assert(pairs1.size() == pairs2.size(),
             ExcMessage("Unmatched faces on periodic boundaries"));

    Assert(
      (pairs1.size() > 0 ||
       (dynamic_cast<
          const parallel::fullydistributed::Triangulation<dim, space_dim> *>(
          &mesh.begin()->get_triangulation()) != nullptr)),
      ExcMessage("No new periodic face pairs have been found. "
                 "Are you sure that you've selected the correct boundary "
                 "id's and that the coarsest level mesh is colorized?"));

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(
      pairs1, pairs2, direction, matched_pairs, offset, matrix, abs_tol);
  }



  /*
   * Internally used in orthogonal_equality
   *
   * An orthogonal equality test for points:
   *
   * point1 and point2 are considered equal, if
   *   matrix.point1 + offset - point2
   * is parallel to the unit vector in <direction>
   */
  template <int spacedim>
  inline bool
  orthogonal_equality(const Point<spacedim>     &point1,
                      const Point<spacedim>     &point2,
                      const unsigned int         direction,
                      const Tensor<1, spacedim> &offset,
                      const FullMatrix<double>  &matrix,
                      const double               abs_tol = 1e-10)
  {
    AssertIndexRange(direction, spacedim);

    Assert(matrix.m() == matrix.n(), ExcInternalError());

    Point<spacedim> distance;

    if (matrix.m() == spacedim)
      for (unsigned int i = 0; i < spacedim; ++i)
        for (unsigned int j = 0; j < spacedim; ++j)
          distance[i] += matrix(i, j) * point1[j];
    else
      distance = point1;

    distance += offset - point2;

    for (unsigned int i = 0; i < spacedim; ++i)
      {
        // Only compare coordinate-components != direction:
        if (i == direction)
          continue;

        if (std::abs(distance[i]) > abs_tol)
          return false;
      }

    return true;
  }



  template <typename FaceIterator>
  std::optional<types::geometric_orientation>
  orthogonal_equality(
    const FaceIterator                                           &face1,
    const FaceIterator                                           &face2,
    const unsigned int                                            direction,
    const Tensor<1, FaceIterator::AccessorType::space_dimension> &offset,
    const FullMatrix<double>                                     &matrix,
    const double                                                  abs_tol)
  {
    Assert(matrix.m() == matrix.n(),
           ExcMessage("The supplied matrix must be a square matrix"));
    Assert(face1->reference_cell() == face2->reference_cell(),
           ExcMessage(
             "The faces to be matched must have the same reference cell."));

    // Do a full matching of the face vertices:
    AssertDimension(face1->n_vertices(), face2->n_vertices());

    std::vector<unsigned int> face1_vertices(face1->n_vertices(),
                                             numbers::invalid_unsigned_int),
      face2_vertices(face2->n_vertices(), numbers::invalid_unsigned_int);

    std::set<unsigned int> face2_vertices_set;
    for (unsigned int i = 0; i < face1->n_vertices(); ++i)
      face2_vertices_set.insert(i);

    for (unsigned int i = 0; i < face1->n_vertices(); ++i)
      {
        for (auto it = face2_vertices_set.begin();
             it != face2_vertices_set.end();
             ++it)
          {
            if (orthogonal_equality(face1->vertex(i),
                                    face2->vertex(*it),
                                    direction,
                                    offset,
                                    matrix,
                                    abs_tol))
              {
                face1_vertices[i] = *it;
                face2_vertices[i] = i;
                face2_vertices_set.erase(it);
                break; // jump out of the innermost loop
              }
          }
      }

    if (face2_vertices_set.empty())
      {
        // Just to be sure, did we fill both arrays with sensible data?
        Assert(face1_vertices.end() ==
                 std::find(face1_vertices.begin(),
                           face1_vertices.begin() + face1->n_vertices(),
                           numbers::invalid_unsigned_int),
               ExcInternalError());
        Assert(face2_vertices.end() ==
                 std::find(face2_vertices.begin(),
                           face2_vertices.begin() + face1->n_vertices(),
                           numbers::invalid_unsigned_int),
               ExcInternalError());

        const auto reference_cell = face1->reference_cell();
        // We want the relative orientation of face1 with respect to face2 so
        // the order is flipped here:
        return std::make_optional(reference_cell.get_combined_orientation(
          make_array_view(face2_vertices.cbegin(),
                          face2_vertices.cbegin() + face2->n_vertices()),
          make_array_view(face1_vertices.cbegin(),
                          face1_vertices.cbegin() + face1->n_vertices())));
      }
    else
      return std::nullopt;
  }
} // namespace GridTools


#include "grid/grid_tools_dof_handlers.inst"


DEAL_II_NAMESPACE_CLOSE
