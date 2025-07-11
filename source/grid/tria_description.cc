// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

DEAL_II_NAMESPACE_OPEN


template <int structdim>
CellData<structdim>::CellData(const unsigned int n_vertices)
  : vertices(n_vertices, numbers::invalid_unsigned_int)
  , material_id(0)
  , manifold_id(numbers::flat_manifold_id)
{}



template <int structdim>
bool
CellData<structdim>::operator==(const CellData<structdim> &other) const
{
  if (vertices.size() != other.vertices.size())
    return false;

  for (unsigned int i = 0; i < vertices.size(); ++i)
    if (vertices[i] != other.vertices[i])
      return false;

  if (material_id != other.material_id)
    return false;

  if (boundary_id != other.boundary_id)
    return false;

  if (manifold_id != other.manifold_id)
    return false;

  return true;
}



bool
SubCellData::check_consistency(const unsigned int dim) const
{
  switch (dim)
    {
      case 1:
        return ((boundary_lines.empty()) && (boundary_quads.empty()));
      case 2:
        return (boundary_quads.empty());
    }
  return true;
}

namespace TriangulationDescription
{
  namespace Utilities
  {
    namespace
    {
      /**
       * A temporal class similar to TriangulationDescription::Description.
       */
      template <int dim, int spacedim>
      struct DescriptionTemp
      {
        /**
         * Serialization function for packing and unpacking the content of this
         * class.
         */
        template <class Archive>
        void
        serialize(Archive &ar, const unsigned int /*version*/)
        {
          ar &coarse_cells;
          ar &coarse_cell_vertices;
          ar &coarse_cell_index_to_coarse_cell_id;
          ar &cell_infos;
        }

        /**
         * Send the DescriptionTemp objects to the right process.
         */
        void
        collect(
          const std::vector<unsigned int> &future_owners_of_locally_owned_cells,
          const std::vector<DescriptionTemp<dim, spacedim>> &description_temp,
          const MPI_Comm                                     comm,
          const bool vertices_have_unique_ids)
        {
          // Use the some-to-some version of the consensus algorithm framework
          // whereby we send requests to other processes that then deal with
          // them but do not send anything back.
          //
          // Note that the input (description_temp) *may* contain an entry for
          // the current process. As documented, the consensus algorithm will
          // simply copy that into the output queue, i.e., it will call
          // process_request() on it as well, and the data will simply come
          // back out on the local process.
          const auto create_request = [&](const unsigned int other_rank) {
            const auto ptr =
              std::find(future_owners_of_locally_owned_cells.begin(),
                        future_owners_of_locally_owned_cells.end(),
                        other_rank);

            Assert(ptr != future_owners_of_locally_owned_cells.end(),
                   ExcInternalError());

            const auto other_rank_index =
              std::distance(future_owners_of_locally_owned_cells.begin(), ptr);

            return description_temp[other_rank_index];
          };

          const auto process_request =
            [&](const unsigned int,
                const DescriptionTemp<dim, spacedim> &request) -> void {
            this->merge(request, vertices_have_unique_ids);
          };

          dealii::Utilities::MPI::ConsensusAlgorithms::selector<
            DescriptionTemp<dim, spacedim>>(
            future_owners_of_locally_owned_cells,
            create_request,
            process_request,
            comm);
        }

        /**
         * Merge the given DescriptionTemp object into the current one. Note
         * that no duplicate cells, vertices, ... are removed at this stage.
         * This done by the reduce() function.
         */
        void
        merge(const DescriptionTemp<dim, spacedim> &other,
              const bool                            vertices_have_unique_ids)
        {
          this->cell_infos.resize(
            std::max(other.cell_infos.size(), this->cell_infos.size()));

          if (vertices_have_unique_ids == false) // need to compare points
            {
              // map point to local vertex index
              std::map<Point<spacedim>,
                       unsigned int,
                       FloatingPointComparator<double>>
                map_point_to_local_vertex_index(
                  FloatingPointComparator<double>(1e-10));

              // ... initialize map with existing points
              for (unsigned int i = 0; i < this->coarse_cell_vertices.size();
                   ++i)
                map_point_to_local_vertex_index[coarse_cell_vertices[i]
                                                  .second] = i;

              // map local vertex indices within other to the new local indices
              std::map<unsigned int, unsigned int>
                map_old_to_new_local_vertex_index;

              // 1) re-enumerate vertices in other and insert into maps
              unsigned int counter = coarse_cell_vertices.size();
              for (const auto &p : other.coarse_cell_vertices)
                if (map_point_to_local_vertex_index.find(p.second) ==
                    map_point_to_local_vertex_index.end())
                  {
                    this->coarse_cell_vertices.emplace_back(counter, p.second);
                    map_point_to_local_vertex_index[p.second] =
                      map_old_to_new_local_vertex_index[p.first] = counter++;
                  }
                else
                  map_old_to_new_local_vertex_index[p.first] =
                    map_point_to_local_vertex_index[p.second];

              // 2) re-enumerate vertices of cells
              auto other_coarse_cells_copy = other.coarse_cells;

              for (auto &cell : other_coarse_cells_copy)
                for (auto &v : cell.vertices)
                  v = map_old_to_new_local_vertex_index[v];

              this->coarse_cells.insert(this->coarse_cells.end(),
                                        other_coarse_cells_copy.begin(),
                                        other_coarse_cells_copy.end());
            }
          else
            {
              this->coarse_cells.insert(this->coarse_cells.end(),
                                        other.coarse_cells.begin(),
                                        other.coarse_cells.end());
              this->coarse_cell_vertices.insert(
                this->coarse_cell_vertices.end(),
                other.coarse_cell_vertices.begin(),
                other.coarse_cell_vertices.end());
            }

          this->coarse_cell_index_to_coarse_cell_id.insert(
            this->coarse_cell_index_to_coarse_cell_id.end(),
            other.coarse_cell_index_to_coarse_cell_id.begin(),
            other.coarse_cell_index_to_coarse_cell_id.end());

          for (unsigned int i = 0; i < this->cell_infos.size(); ++i)
            this->cell_infos[i].insert(this->cell_infos[i].end(),
                                       other.cell_infos[i].begin(),
                                       other.cell_infos[i].end());
        }

        /**
         * Remove all duplicate information obtained during merge().
         */
        void
        reduce()
        {
          // make coarse cells unique
          {
            std::vector<std::tuple<types::coarse_cell_id,
                                   dealii::CellData<dim>,
                                   unsigned int>>
              temp;

            temp.reserve(this->coarse_cells.size());
            for (unsigned int i = 0; i < this->coarse_cells.size(); ++i)
              temp.emplace_back(this->coarse_cell_index_to_coarse_cell_id[i],
                                this->coarse_cells[i],
                                i);

            std::sort(temp.begin(),
                      temp.end(),
                      [](const auto &a, const auto &b) {
                        return std::get<0>(a) < std::get<0>(b);
                      });
            temp.erase(std::unique(temp.begin(),
                                   temp.end(),
                                   [](const auto &a, const auto &b) {
                                     return std::get<0>(a) == std::get<0>(b);
                                   }),
                       temp.end());
            std::sort(temp.begin(),
                      temp.end(),
                      [](const auto &a, const auto &b) {
                        return std::get<2>(a) < std::get<2>(b);
                      });

            this->coarse_cell_index_to_coarse_cell_id.resize(temp.size());
            this->coarse_cells.resize(temp.size());

            for (unsigned int i = 0; i < temp.size(); ++i)
              {
                this->coarse_cell_index_to_coarse_cell_id[i] =
                  std::get<0>(temp[i]);
                this->coarse_cells[i] = std::get<1>(temp[i]);
              }
          }

          // make coarse cell vertices unique
          {
            std::sort(this->coarse_cell_vertices.begin(),
                      this->coarse_cell_vertices.end(),
                      [](const std::pair<unsigned int, Point<spacedim>> &a,
                         const std::pair<unsigned int, Point<spacedim>> &b) {
                        return a.first < b.first;
                      });
            this->coarse_cell_vertices.erase(
              std::unique(
                this->coarse_cell_vertices.begin(),
                this->coarse_cell_vertices.end(),
                [](const std::pair<unsigned int, Point<spacedim>> &a,
                   const std::pair<unsigned int, Point<spacedim>> &b) {
                  if (a.first == b.first)
                    {
                      Assert(a.second.distance(b.second) <=
                               1e-7 *
                                 std::max(a.second.norm(), b.second.norm()),
                             ExcMessage(
                               "In the process of merging the vertices of "
                               "the coarse meshes used on different processes, "
                               "there were two processes that used the same "
                               "vertex index for points that are not the same. "
                               "This suggests that you are using different "
                               "coarse meshes on different processes. This "
                               "should not happen."));
                      return true;
                    }
                  return false;
                }),
              this->coarse_cell_vertices.end());
          }

          // make cells unique
          for (unsigned int i = 0; i < this->cell_infos.size(); ++i)
            {
              if (this->cell_infos[i].empty())
                continue;

              std::sort(this->cell_infos[i].begin(),
                        this->cell_infos[i].end(),
                        [](const auto &a, const auto &b) {
                          return a.id < b.id;
                        });

              std::vector<CellData<dim>> temp;
              temp.push_back(this->cell_infos[i][0]);

              for (unsigned int j = 1; j < this->cell_infos[i].size(); ++j)
                if (temp.back().id == cell_infos[i][j].id)
                  {
                    temp.back().subdomain_id =
                      std::min(temp.back().subdomain_id,
                               this->cell_infos[i][j].subdomain_id);
                    temp.back().level_subdomain_id =
                      std::min(temp.back().level_subdomain_id,
                               this->cell_infos[i][j].level_subdomain_id);
                  }
                else
                  {
                    temp.push_back(this->cell_infos[i][j]);
                  }

              this->cell_infos[i] = temp;
            }
        }

        /**
         * Convert the content of this class to a
         * TriangulationDescription::Description, which can be used during
         * Triangulation::create_triangulation().
         */
        Description<dim, spacedim>
        convert(const MPI_Comm comm,
                const typename Triangulation<dim, spacedim>::MeshSmoothing
                                                         mesh_smoothing,
                const TriangulationDescription::Settings settings)
        {
          Description<dim, spacedim> description;

          // copy communicator
          description.comm = comm;

          description.settings = settings;

          // use mesh smoothing from base triangulation
          description.smoothing = mesh_smoothing;

          std::map<unsigned int, unsigned int> map;

          for (unsigned int i = 0; i < this->coarse_cell_vertices.size(); ++i)
            {
              description.coarse_cell_vertices.push_back(
                this->coarse_cell_vertices[i].second);
              map[this->coarse_cell_vertices[i].first] = i;
            }

          description.coarse_cells = this->coarse_cells;

          for (auto &cell : description.coarse_cells)
            for (unsigned int v = 0; v < cell.vertices.size(); ++v)
              cell.vertices[v] = map[cell.vertices[v]];

          description.coarse_cell_index_to_coarse_cell_id =
            this->coarse_cell_index_to_coarse_cell_id;
          description.cell_infos = this->cell_infos;

          return description;
        }

        std::vector<dealii::CellData<dim>> coarse_cells;

        std::vector<std::pair<unsigned int, Point<spacedim>>>
          coarse_cell_vertices;

        std::vector<types::coarse_cell_id> coarse_cell_index_to_coarse_cell_id;

        std::vector<std::vector<CellData<dim>>> cell_infos;
      };

      /**
       * Set the user_flag of a cell and of all its parent cells.
       */
      template <int dim, int spacedim>
      void
      mark_cell_and_its_parents(
        const TriaIterator<CellAccessor<dim, spacedim>> &cell,
        std::vector<std::vector<bool>>                  &cell_marked)
      {
        cell_marked[cell->level()][cell->index()] = true;
        if (cell->level() != 0)
          mark_cell_and_its_parents(cell->parent(), cell_marked);
      }

      /**
       * A helper function for the
       * TriangulationDescription::Utilities::create_description_from_triangulation()
       * function.
       */
      template <typename DescriptionType, int dim, int spacedim>
      DescriptionType
      create_description_for_rank(
        const dealii::Triangulation<dim, spacedim> &tria,
        const std::function<types::subdomain_id(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &)>
          &subdomain_id_function,
        const std::function<types::subdomain_id(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &)>
          &level_subdomain_id_function,
        const std::map<unsigned int, std::vector<unsigned int>>
          &coinciding_vertex_groups,
        const std::map<unsigned int, unsigned int>
                          &vertex_to_coinciding_vertex_group,
        const MPI_Comm     comm,
        const unsigned int my_rank,
        const TriangulationDescription::Settings settings)
      {
        static_assert(
          std::is_same_v<DescriptionType, Description<dim, spacedim>> ||
            std::is_same_v<DescriptionType, DescriptionTemp<dim, spacedim>>,
          "Wrong template type.");
        Assert(
          !(settings & TriangulationDescription::Settings::
                         construct_multigrid_hierarchy) ||
            (tria.get_mesh_smoothing() &
             Triangulation<dim, spacedim>::limit_level_difference_at_vertices),
          ExcMessage(
            "Source triangulation has to be set up with "
            "limit_level_difference_at_vertices if the construction of the "
            "multigrid hierarchy is requested!"));

        const bool construct_multigrid =
          ((settings & TriangulationDescription::Settings::
                         construct_multigrid_hierarchy) != 0u);

        DescriptionType construction_data;
        if constexpr (std::is_same_v<DescriptionType,
                                     Description<dim, spacedim>>)
          {
            construction_data.comm      = comm;
            construction_data.smoothing = tria.get_mesh_smoothing();
            construction_data.settings  = settings;
          }
        else
          (void)comm;

        // A helper function that marks the indices of all vertices belonging
        // to a cell (also taking into account their periodic breathren) in
        // the bit vector passed as second argument.
        const auto
          add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells =
            [&coinciding_vertex_groups, &vertex_to_coinciding_vertex_group](
              const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                &cell,
              std::vector<bool> &vertices_on_locally_owned_cells) {
              for (const unsigned int v : cell->vertex_indices())
                {
                  const auto global_vertex_index = cell->vertex_index(v);
                  vertices_on_locally_owned_cells[global_vertex_index] = true;

                  const auto coinciding_vertex_group =
                    vertex_to_coinciding_vertex_group.find(global_vertex_index);
                  if (coinciding_vertex_group !=
                      vertex_to_coinciding_vertex_group.end())
                    for (const auto &co_vertex : coinciding_vertex_groups.at(
                           coinciding_vertex_group->second))
                      vertices_on_locally_owned_cells[co_vertex] = true;
                }
            };

        const auto add_vertices =
          [&tria](const std::vector<bool> &vertices_locally_relevant,
                  DescriptionType         &construction_data) {
            if constexpr (std::is_same_v<DescriptionType,
                                         Description<dim, spacedim>>)
              {
                std::vector<unsigned int> vertices_locally_relevant_indices(
                  vertices_locally_relevant.size());

                // enumerate locally relevant vertices
                unsigned int vertex_counter = 0;
                for (unsigned int i = 0; i < vertices_locally_relevant.size();
                     ++i)
                  if (vertices_locally_relevant[i])
                    {
                      construction_data.coarse_cell_vertices.push_back(
                        tria.get_vertices()[i]);
                      vertices_locally_relevant_indices[i] = vertex_counter++;
                    }

                // correct vertices of cells (make them local)
                for (auto &cell : construction_data.coarse_cells)
                  for (unsigned int v = 0; v < cell.vertices.size(); ++v)
                    cell.vertices[v] =
                      vertices_locally_relevant_indices[cell.vertices[v]];
              }
            else
              {
                for (unsigned int i = 0; i < vertices_locally_relevant.size();
                     ++i)
                  if (vertices_locally_relevant[i])
                    construction_data.coarse_cell_vertices.emplace_back(
                      i, tria.get_vertices()[i]);
              }
          };


        // 1) loop over levels (from fine to coarse) and mark on each level
        //    the locally relevant cells
        std::vector<std::vector<bool>> cell_marked(tria.n_levels());
        for (unsigned int l = 0; l < tria.n_levels(); ++l)
          cell_marked[l].resize(tria.n_raw_cells(l));

        for (int level = tria.get_triangulation().n_global_levels() - 1;
             level >= 0;
             --level)
          {
            // collect vertices connected to a (on any level) locally owned
            // cell
            std::vector<bool> vertices_owned_by_locally_owned_cells_on_level(
              tria.n_vertices());
            for (const auto &cell : tria.cell_iterators_on_level(level))
              if (construct_multigrid &&
                  (level_subdomain_id_function(cell) == my_rank))
                add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                  cell, vertices_owned_by_locally_owned_cells_on_level);

            for (const auto &cell : tria.active_cell_iterators())
              if (subdomain_id_function(cell) == my_rank)
                add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                  cell, vertices_owned_by_locally_owned_cells_on_level);

            // helper function to determine if cell is locally relevant
            // (i.e. a cell which is connected to a vertex via a locally
            // owned cell)
            const auto is_locally_relevant_on_level = [&](const auto &cell) {
              for (const auto v : cell->vertex_indices())
                if (vertices_owned_by_locally_owned_cells_on_level
                      [cell->vertex_index(v)])
                  return true;
              return false;
            };

            // mark all locally relevant cells
            for (const auto &cell : tria.cell_iterators_on_level(level))
              if (is_locally_relevant_on_level(cell))
                mark_cell_and_its_parents(cell, cell_marked);
          }

        // 2) set_up coarse-grid triangulation
        {
          std::vector<bool> vertices_locally_relevant(tria.n_vertices(), false);

          // a) loop over all cells
          for (const auto &cell : tria.cell_iterators_on_level(0))
            {
              if (!cell_marked[cell->level()][cell->index()])
                continue;

              // extract cell definition (with old numbering of vertices)
              dealii::CellData<dim> cell_data(cell->n_vertices());
              cell_data.material_id = cell->material_id();
              cell_data.manifold_id = cell->manifold_id();
              for (const auto v : cell->vertex_indices())
                cell_data.vertices[v] = cell->vertex_index(v);
              construction_data.coarse_cells.push_back(cell_data);

              // save indices of each vertex of this cell
              for (const auto v : cell->vertex_indices())
                vertices_locally_relevant[cell->vertex_index(v)] = true;

              // save translation for corase grid: lid -> gid
              construction_data.coarse_cell_index_to_coarse_cell_id.push_back(
                cell->id().get_coarse_cell_id());
            }

          add_vertices(vertices_locally_relevant, construction_data);
        }


        // 3) collect info of each cell
        construction_data.cell_infos.resize(
          tria.get_triangulation().n_global_levels());

        // collect local vertices on active level
        std::vector<bool> vertices_owned_by_locally_owned_active_cells(
          tria.n_vertices());
        for (const auto &cell : tria.active_cell_iterators())
          if (subdomain_id_function(cell) == my_rank)
            add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
              cell, vertices_owned_by_locally_owned_active_cells);

        // helper function to determine if cell is locally relevant
        // on active level
        const auto is_locally_relevant_on_active_level = [&](const auto &cell) {
          if (cell->is_active())
            for (const auto v : cell->vertex_indices())
              if (vertices_owned_by_locally_owned_active_cells
                    [cell->vertex_index(v)])
                return true;
          return false;
        };

        for (unsigned int level = 0;
             level < tria.get_triangulation().n_global_levels();
             ++level)
          {
            // collect local vertices on level
            std::vector<bool> vertices_owned_by_locally_owned_cells_on_level(
              tria.n_vertices());
            for (const auto &cell : tria.cell_iterators_on_level(level))
              if ((construct_multigrid &&
                   (level_subdomain_id_function(cell) == my_rank)) ||
                  (cell->is_active() && subdomain_id_function(cell) == my_rank))
                add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                  cell, vertices_owned_by_locally_owned_cells_on_level);

            // helper function to determine if cell is locally relevant
            // on level
            const auto is_locally_relevant_on_level = [&](const auto &cell) {
              for (const auto v : cell->vertex_indices())
                if (vertices_owned_by_locally_owned_cells_on_level
                      [cell->vertex_index(v)])
                  return true;
              return false;
            };

            auto &level_cell_infos = construction_data.cell_infos[level];
            for (const auto &cell : tria.cell_iterators_on_level(level))
              {
                // check if cell is locally relevant
                if (!cell_marked[cell->level()][cell->index()])
                  continue;

                CellData<dim> cell_info;

                // save coarse-cell id
                cell_info.id = cell->id().template to_binary<dim>();

                // save boundary_ids of each face of this cell
                for (const auto f : cell->face_indices())
                  {
                    types::boundary_id boundary_ind =
                      cell->face(f)->boundary_id();
                    if (boundary_ind != numbers::internal_face_boundary_id)
                      cell_info.boundary_ids.emplace_back(f, boundary_ind);
                  }

                // save manifold id
                {
                  // ... of cell
                  cell_info.manifold_id = cell->manifold_id();

                  // ... of lines
                  if (dim >= 2)
                    for (const auto line : cell->line_indices())
                      cell_info.manifold_line_ids[line] =
                        cell->line(line)->manifold_id();

                  // ... of quads
                  if (dim == 3)
                    for (const auto f : cell->face_indices())
                      cell_info.manifold_quad_ids[f] =
                        cell->quad(f)->manifold_id();
                }

                // subdomain and level subdomain id
                cell_info.subdomain_id       = numbers::artificial_subdomain_id;
                cell_info.level_subdomain_id = numbers::artificial_subdomain_id;

                if (is_locally_relevant_on_active_level(cell))
                  {
                    cell_info.subdomain_id = subdomain_id_function(cell);

                    cell_info.level_subdomain_id =
                      level_subdomain_id_function(cell);
                  }
                else if (is_locally_relevant_on_level(cell))
                  {
                    cell_info.level_subdomain_id =
                      level_subdomain_id_function(cell);
                  }
                else
                  {
                    // cell is locally relevant but an artificial cell
                  }

                level_cell_infos.emplace_back(cell_info);
              }
          }

        return construction_data;
      }
    } // namespace


    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const dealii::Triangulation<dim, spacedim> &tria,
      const MPI_Comm                              comm,
      const TriangulationDescription::Settings    settings,
      const unsigned int                          my_rank_in)
    {
      if (const auto ptria =
            dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &tria))
        {
          Assert(comm == ptria->get_mpi_communicator(),
                 ExcMessage("MPI communicators do not match."));
          Assert(my_rank_in == numbers::invalid_unsigned_int ||
                   my_rank_in == dealii::Utilities::MPI::this_mpi_process(comm),
                 ExcMessage(
                   "For creation from a parallel::Triangulation, "
                   "my_rank has to equal the rank of the current process "
                   "in the given communicator."));
        }

      if constexpr (running_in_debug_mode())
        {
          // If we are dealing with a sequential triangulation, then someone
          // will have needed to set the subdomain_ids by hand. Make sure that
          // all ids we see are less than the number of processes we are
          // supposed to split the triangulation into.
          if (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
                &tria) == nullptr)
            {
              const unsigned int n_mpi_processes =
                dealii::Utilities::MPI::n_mpi_processes(comm);
              for (const auto &cell : tria.active_cell_iterators())
                Assert(cell->subdomain_id() < n_mpi_processes,
                       ExcMessage(
                         "You can't have a cell with subdomain_id of " +
                         std::to_string(cell->subdomain_id()) +
                         " when splitting the triangulation using an MPI "
                         " communicator with only " +
                         std::to_string(n_mpi_processes) + " processes."));
            }
        }

      // First, figure out for what rank we are supposed to build the
      // TriangulationDescription::Description object
      const unsigned int my_rank =
        (my_rank_in == numbers::invalid_unsigned_int ?
           dealii::Utilities::MPI::this_mpi_process(comm) :
           my_rank_in);

      const auto subdomain_id_function = [](const auto &cell) {
        return cell->subdomain_id();
      };

      const auto level_subdomain_id_function = [](const auto &cell) {
        return cell->level_subdomain_id();
      };

      std::map<unsigned int, std::vector<unsigned int>>
                                           coinciding_vertex_groups;
      std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;
      GridTools::collect_coinciding_vertices(tria,
                                             coinciding_vertex_groups,
                                             vertex_to_coinciding_vertex_group);

      return create_description_for_rank<Description<dim, spacedim>>(
        tria,
        subdomain_id_function,
        level_subdomain_id_function,
        coinciding_vertex_groups,
        vertex_to_coinciding_vertex_group,
        comm,
        my_rank,
        settings);
    }



    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation_in_groups(
      const std::function<void(dealii::Triangulation<dim, spacedim> &)>
                                                    &serial_grid_generator,
      const std::function<void(dealii::Triangulation<dim, spacedim> &,
                               const MPI_Comm,
                               const unsigned int)> &serial_grid_partitioner,
      const MPI_Comm                                 comm,
      const int                                      group_size,
      const typename Triangulation<dim, spacedim>::MeshSmoothing smoothing,
      const TriangulationDescription::Settings                   settings)
    {
#ifndef DEAL_II_WITH_MPI
      (void)serial_grid_generator;
      (void)serial_grid_partitioner;
      (void)comm;
      (void)group_size;
      (void)smoothing;
      (void)settings;

      return Description<dim, spacedim>();
#else
      const unsigned int my_rank =
        dealii::Utilities::MPI::this_mpi_process(comm);
      const unsigned int group_root = (my_rank / group_size) * group_size;

      const int mpi_tag =
        dealii::Utilities::MPI::internal::Tags::fully_distributed_create;

      // check if process is root of the group
      if (my_rank == group_root)
        {
          // Step 1: create serial triangulation
          dealii::Triangulation<dim, spacedim> tria(
            (settings & TriangulationDescription::Settings::
                          construct_multigrid_hierarchy) ?
              static_cast<
                typename dealii::Triangulation<dim, spacedim>::MeshSmoothing>(
                smoothing |
                Triangulation<dim,
                              spacedim>::limit_level_difference_at_vertices) :
              smoothing);
          serial_grid_generator(tria);

          // Step 2: partition active cells and ...
          serial_grid_partitioner(tria, comm, group_size);

          // ... cells on the levels if multigrid is required
          if (settings &
              TriangulationDescription::Settings::construct_multigrid_hierarchy)
            GridTools::partition_multigrid_levels(tria);

          const unsigned int end_group =
            std::min(group_root + group_size,
                     dealii::Utilities::MPI::n_mpi_processes(comm));

          // 3) create Description for the other processes in group; since
          // we expect that this function is called for huge meshes, one
          // Description is created at a time and send away; only once the
          // Description has been sent away, the next rank is processed.
          for (unsigned int other_rank = group_root + 1; other_rank < end_group;
               ++other_rank)
            {
              // 3a) create construction data for other ranks
              const auto construction_data =
                create_description_from_triangulation(tria,
                                                      comm,
                                                      settings,
                                                      other_rank);
              // 3b) pack
              std::vector<char> buffer;
              dealii::Utilities::pack(construction_data, buffer, false);

              // 3c) send TriangulationDescription::Description
              const auto ierr = MPI_Send(buffer.data(),
                                         buffer.size(),
                                         MPI_CHAR,
                                         other_rank,
                                         mpi_tag,
                                         comm);
              AssertThrowMPI(ierr);
            }

          // 4) create TriangulationDescription::Description for this process
          // (root of group)
          return create_description_from_triangulation(tria,
                                                       comm,
                                                       settings,
                                                       my_rank);
        }
      else
        {
          // 3a) recv packed TriangulationDescription::Description from
          // group-root process
          //     (counter-part of 3c of root process)
          MPI_Status status;
          auto       ierr = MPI_Probe(group_root, mpi_tag, comm, &status);
          AssertThrowMPI(ierr);

          int len;
          MPI_Get_count(&status, MPI_CHAR, &len);

          std::vector<char> buf(len);
          ierr = MPI_Recv(buf.data(),
                          len,
                          MPI_CHAR,
                          status.MPI_SOURCE,
                          mpi_tag,
                          comm,
                          &status);
          AssertThrowMPI(ierr);

          // 3b) unpack TriangulationDescription::Description (counter-part of
          // 3b of root process)
          auto construction_data =
            dealii::Utilities::template unpack<Description<dim, spacedim>>(
              buf, false);

          // WARNING: serialization cannot handle the MPI communicator
          // which is the reason why we have to set it here explicitly
          construction_data.comm = comm;

          return construction_data;
        }
#endif
    }



    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const Triangulation<dim, spacedim>               &tria,
      const LinearAlgebra::distributed::Vector<double> &partition,
      const TriangulationDescription::Settings          settings)
    {
      const bool construct_multigrid =
        (partition.size() > 0) &&
        (settings &
         TriangulationDescription::Settings::construct_multigrid_hierarchy);

      Assert(
        construct_multigrid == false ||
          (tria.get_mesh_smoothing() &
           Triangulation<dim, spacedim>::limit_level_difference_at_vertices),
        ExcMessage(
          "Source triangulation has to be set up with "
          "limit_level_difference_at_vertices if the construction of the "
          "multigrid hierarchy is requested!"));

      std::vector<LinearAlgebra::distributed::Vector<double>> partitions_mg;

      // If desired, also create a multigrid hierarchy. For this, we have to
      // build a hierarchy of partitions (one for each level of the
      // triangulation) in which each cell is assigned to the same process
      // as its first child (if not active) or to the same process that already
      // owns the cell (for an active level-cell).
      if (construct_multigrid)
        {
          const auto tria_parallel =
            dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &tria);
          Assert(tria_parallel, ExcNotImplemented());

          // Give the level partitioners the right size:
          partitions_mg.resize(tria.n_global_levels());
          for (unsigned int l = 0; l < tria.n_global_levels(); ++l)
            partitions_mg[l].reinit(
              tria_parallel->global_level_cell_index_partitioner(l).lock());

          // Make sure we know about all of the owners of the active cells,
          // whether locally owned or not. Then we traverse the triangulation
          // from the finest level to the coarsest level.
          //
          // On each level, traverse the cell. If the cell is not locally
          // owned, we don't care about it. If it is active, we copy the
          // owner process from the cell's non-level owner. Otherwise,
          // use the owner of the first cell.
          partition.update_ghost_values();
          for (int level = tria.n_global_levels() - 1; level >= 0; --level)
            {
              for (const auto &cell : tria.cell_iterators_on_level(level))
                {
                  if (cell->is_locally_owned_on_level())
                    {
                      if (cell->is_active())
                        partitions_mg[level][cell->global_level_cell_index()] =
                          partition[cell->global_active_cell_index()];
                      else
                        partitions_mg[level][cell->global_level_cell_index()] =
                          partitions_mg[level + 1]
                                       [cell->child(0)
                                          ->global_level_cell_index()];
                    }
                }

              // Having touched all of the locally owned cells on the
              // current level, exchange information with the other processes
              // about the cells that are ghosts so that on the next coarser
              // level we can access information about children again:
              partitions_mg[level].update_ghost_values();
            }
        }

      // Forward to the other function.
      return create_description_from_triangulation(tria,
                                                   partition,
                                                   partitions_mg,
                                                   settings);
    }



    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const Triangulation<dim, spacedim>               &tria,
      const LinearAlgebra::distributed::Vector<double> &partition,
      const std::vector<LinearAlgebra::distributed::Vector<double>>
                                              &partitions_mg,
      const TriangulationDescription::Settings settings_in)
    {
#ifdef DEAL_II_WITH_MPI
      if (tria.get_mpi_communicator() == MPI_COMM_NULL)
        AssertDimension(partition.locally_owned_size(), 0);
#endif

      if (partition.size() == 0)
        {
          AssertDimension(partitions_mg.size(), 0);
          return create_description_from_triangulation(
            tria, tria.get_mpi_communicator(), settings_in);
        }

      // Update partitioner ghost elements because we will later want
      // to ask also about the future owners of ghost cells.
      partition.update_ghost_values();
      for (const auto &partition : partitions_mg)
        partition.update_ghost_values();

      // 1) Determine process ids that appear on locally owned cells. Create
      //    a sorted vector by first creating a std::set and then copying
      //    the result. (Note that we get only locally *owned* cells in
      //    the output because we only loop over the locally *owned*
      //    entries of the partitioning vector, even though
      //    'partition.local_element(i)' could also return locally relevant
      //    elements if 'i' were to exceed the number of locally owned
      //    elements.)
      const std::vector<unsigned int> future_owners_of_locally_owned_cells =
        [&partition, &partitions_mg]() {
          std::set<unsigned int> relevant_process_set;

          const unsigned int n_mpi_ranks =
            dealii::Utilities::MPI::n_mpi_processes(
              partition.get_mpi_communicator());

          for (unsigned int i = 0; i < partition.locally_owned_size(); ++i)
            {
              Assert(static_cast<unsigned int>(partition.local_element(i)) ==
                       partition.local_element(i),
                     ExcMessage(
                       "The elements of a partition vector must be integers."));
              Assert(
                partition.local_element(i) < n_mpi_ranks,
                ExcMessage(
                  "The elements of a partition vector must be between zero "
                  "and the number of processes in the communicator "
                  "to be used for partitioning the triangulation."));
              relevant_process_set.insert(
                static_cast<unsigned int>(partition.local_element(i)));
            }

          for (const auto &partition : partitions_mg)
            for (unsigned int i = 0; i < partition.locally_owned_size(); ++i)
              {
                Assert(
                  static_cast<unsigned int>(partition.local_element(i)) ==
                    partition.local_element(i),
                  ExcMessage(
                    "The elements of a partition vector must be integers."));
                Assert(
                  partition.local_element(i) < n_mpi_ranks,
                  ExcMessage(
                    "The elements of a partition vector must be between zero "
                    "and the number of processes in the communicator "
                    "to be used for partitioning the triangulation."));
                relevant_process_set.insert(
                  static_cast<unsigned int>(partition.local_element(i)));
              }

          return std::vector<unsigned int>(relevant_process_set.begin(),
                                           relevant_process_set.end());
        }();

      const bool construct_multigrid = (partitions_mg.size() > 0);

      const TriangulationDescription::Settings settings =
        (construct_multigrid ?
           static_cast<TriangulationDescription::Settings>(
             settings_in | TriangulationDescription::Settings::
                             construct_multigrid_hierarchy) :
           settings_in);


      // Set up a function that returns the future owner rank for a cell.
      // Same then also for the level owner. These functions work for
      // locally owned and ghost cells.
      const auto cell_to_future_owner =
        [&partition](const auto &cell) -> types::subdomain_id {
        if ((cell->is_active() && (cell->is_artificial() == false)))
          return static_cast<types::subdomain_id>(
            partition[cell->global_active_cell_index()]);
        else
          return numbers::artificial_subdomain_id;
      };

      const auto mg_cell_to_future_owner =
        [&construct_multigrid,
         &partitions_mg](const auto &cell) -> types::subdomain_id {
        if (construct_multigrid && (cell->is_artificial_on_level() == false))
          return static_cast<types::subdomain_id>(
            partitions_mg[cell->level()][cell->global_level_cell_index()]);
        else
          return numbers::artificial_subdomain_id;
      };

      // Create a description (locally owned cell and a layer of ghost cells
      // and all their parents). We first create a description in the
      // 'temporary' format (using class DescriptionTemp), which we will
      // later convert to its final form.
      std::vector<DescriptionTemp<dim, spacedim>> descriptions_per_rank;
      descriptions_per_rank.reserve(
        future_owners_of_locally_owned_cells.size());

      std::map<unsigned int, std::vector<unsigned int>>
                                           coinciding_vertex_groups;
      std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;
      GridTools::collect_coinciding_vertices(tria,
                                             coinciding_vertex_groups,
                                             vertex_to_coinciding_vertex_group);

      for (const auto rank : future_owners_of_locally_owned_cells)
        descriptions_per_rank.emplace_back(
          create_description_for_rank<DescriptionTemp<dim, spacedim>>(
            tria,
            cell_to_future_owner,
            mg_cell_to_future_owner,
            coinciding_vertex_groups,
            vertex_to_coinciding_vertex_group,
            tria.get_mpi_communicator(),
            rank,
            settings));

      // Collect description from all processes that used to own locally-owned
      // active cells of this process in a single description
      DescriptionTemp<dim, spacedim> description_merged;
      description_merged.collect(
        future_owners_of_locally_owned_cells,
        descriptions_per_rank,
        partition.get_mpi_communicator(),
        dynamic_cast<
          const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
          &tria) == nullptr);

      // remove redundant entries
      description_merged.reduce();

      // convert to actual description
      return description_merged.convert(partition.get_mpi_communicator(),
                                        tria.get_mesh_smoothing(),
                                        settings);
    }

  } // namespace Utilities
} // namespace TriangulationDescription



/*-------------- Explicit Instantiations -------------------------------*/
#include "grid/tria_description.inst"


DEAL_II_NAMESPACE_CLOSE
