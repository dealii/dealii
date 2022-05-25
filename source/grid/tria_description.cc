// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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
        return ((boundary_lines.size() == 0) && (boundary_quads.size() == 0));
      case 2:
        return (boundary_quads.size() == 0);
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
          const std::vector<unsigned int> &                  relevant_processes,
          const std::vector<DescriptionTemp<dim, spacedim>> &description_temp,
          const MPI_Comm &                                   comm,
          const bool vertices_have_unique_ids)
        {
          const auto create_request = [&](const unsigned int other_rank) {
            const auto ptr = std::find(relevant_processes.begin(),
                                       relevant_processes.end(),
                                       other_rank);

            Assert(ptr != relevant_processes.end(), ExcInternalError());

            const auto other_rank_index =
              std::distance(relevant_processes.begin(), ptr);

            return description_temp[other_rank_index];
          };

          const auto process_request =
            [&](const unsigned int,
                const DescriptionTemp<dim, spacedim> &request) {
              this->merge(request, vertices_have_unique_ids);
            };

          dealii::Utilities::MPI::ConsensusAlgorithms::selector<
            DescriptionTemp<dim, spacedim>>(relevant_processes,
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
              // comparator of points
              const auto comp = [](const auto &a, const auto &b) {
                const double tolerance = 1e-10;

                for (unsigned int i = 0; i < spacedim; ++i)
                  if (std::abs(a[i] - b[i]) > tolerance)
                    return a[i] < b[i];

                return false;
              };

              // map point to local vertex index
              std::map<Point<spacedim>, unsigned int, decltype(comp)>
                map_point_to_local_vertex_index(comp);

              // ... initialize map with existing points
              for (unsigned int i = 0; i < this->coarse_cell_vertices.size();
                   ++i)
                map_point_to_local_vertex_index[coarse_cell_vertices[i]
                                                  .second] = i;

              // map local vertex indices within other to the new local indices
              std::map<unsigned int, unsigned int>
                map_old_to_new_local_vertex_index;

              // 1) renumerate vertices in other and insert into maps
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

              // 2) renumerate vertices of cells
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
                      [](const auto &a, const auto &b) {
                        return a.first < b.first;
                      });
            this->coarse_cell_vertices.erase(
              std::unique(this->coarse_cell_vertices.begin(),
                          this->coarse_cell_vertices.end(),
                          [](const auto &a, const auto &b) {
                            if (a.first == b.first)
                              {
                                Assert(a.second.distance(b.second) < 10e-8,
                                       ExcInternalError());
                                return true;
                              }
                            return false;
                          }),
              this->coarse_cell_vertices.end());
          }

          // make cells unique
          for (unsigned int i = 0; i < this->cell_infos.size(); ++i)
            {
              if (this->cell_infos[i].size() == 0)
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
      set_user_flag_and_of_its_parents(
        const TriaIterator<CellAccessor<dim, spacedim>> &cell)
      {
        cell->set_user_flag();
        if (cell->level() != 0)
          set_user_flag_and_of_its_parents(cell->parent());
      }

      /**
       * Helper calls for the
       * TriangulationDescription::Utilities::create_description_from_triangulation()
       * functions.
       */
      template <int dim, int spacedim>
      class CreateDescriptionFromTriangulationHelper
      {
      public:
        CreateDescriptionFromTriangulationHelper(
          const dealii::Triangulation<dim, spacedim> &tria,
          const std::function<types::subdomain_id(
            const typename dealii::Triangulation<dim, spacedim>::cell_iterator
              &)> &                                   subdomain_id_function,
          const std::function<types::subdomain_id(
            const typename dealii::Triangulation<dim, spacedim>::cell_iterator
              &)> &                                level_subdomain_id_function,
          const MPI_Comm &                         comm,
          const TriangulationDescription::Settings settings)
          : tria(tria)
          , subdomain_id_function(subdomain_id_function)
          , level_subdomain_id_function(level_subdomain_id_function)
          , comm(comm)
          , settings(settings)
          , construct_multigrid((settings & TriangulationDescription::Settings::
                                              construct_multigrid_hierarchy) !=
                                0u)
        {
          Assert(
            !(settings & TriangulationDescription::Settings::
                           construct_multigrid_hierarchy) ||
              (tria.get_mesh_smoothing() &
               Triangulation<dim,
                             spacedim>::limit_level_difference_at_vertices),
            ExcMessage(
              "Source triangulation has to be set up with "
              "limit_level_difference_at_vertices if the construction of the "
              "multigrid hierarchy is requested!"));

          GridTools::collect_coinciding_vertices(
            tria, coinciding_vertex_groups, vertex_to_coinciding_vertex_group);
        }

        template <typename T>
        T
        create_description_for_rank(const unsigned int my_rank) const
        {
          T construction_data;

          set_additional_data(construction_data);

          // helper function, which collects all vertices belonging to a cell
          // (also taking into account periodicity)
          const auto
            add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells =
              [this](const auto &       cell,
                     std::vector<bool> &vertices_owned_by_locally_owned_cells) {
                // add local vertices
                for (const auto v : cell->vertex_indices())
                  {
                    vertices_owned_by_locally_owned_cells[cell->vertex_index(
                      v)] = true;
                    const auto coinciding_vertex_group =
                      vertex_to_coinciding_vertex_group.find(
                        cell->vertex_index(v));
                    if (coinciding_vertex_group !=
                        vertex_to_coinciding_vertex_group.end())
                      for (const auto &co_vertex : coinciding_vertex_groups.at(
                             coinciding_vertex_group->second))
                        vertices_owned_by_locally_owned_cells[co_vertex] = true;
                  }
              };

          // 1) collect locally relevant cells (set user_flag)
          std::vector<bool> old_user_flags;
          tria.save_user_flags(old_user_flags);

          // 1a) clear user_flags
          const_cast<dealii::Triangulation<dim, spacedim> &>(tria)
            .clear_user_flags();

          // 1b) loop over levels (from fine to coarse) and mark on each level
          //     the locally relevant cells
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
                  set_user_flag_and_of_its_parents(cell);
            }

          // 2) set_up coarse-grid triangulation
          {
            std::vector<bool> vertices_locally_relevant(tria.n_vertices(),
                                                        false);

            // a) loop over all cells
            for (const auto &cell : tria.cell_iterators_on_level(0))
              {
                if (!cell->user_flag_set())
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

            add_vertices(construction_data, vertices_locally_relevant);
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
          const auto is_locally_relevant_on_active_level =
            [&](const auto &cell) {
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
                    (cell->is_active() &&
                     subdomain_id_function(cell) == my_rank))
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
                  if (!(cell->user_flag_set()))
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
                  cell_info.subdomain_id = numbers::artificial_subdomain_id;
                  cell_info.level_subdomain_id =
                    numbers::artificial_subdomain_id;

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

          const_cast<dealii::Triangulation<dim, spacedim> &>(tria)
            .load_user_flags(old_user_flags);

          return construction_data;
        }

      private:
        void
        set_additional_data(Description<dim, spacedim> &construction_data) const
        {
          construction_data.comm      = comm;
          construction_data.smoothing = tria.get_mesh_smoothing();
          construction_data.settings  = settings;
        }

        void
        set_additional_data(DescriptionTemp<dim, spacedim> &) const
        {
          // nothing to do
        }

        void
        add_vertices(Description<dim, spacedim> &construction_data,
                     const std::vector<bool> &vertices_locally_relevant) const
        {
          std::vector<unsigned int> vertices_locally_relevant_indices(
            vertices_locally_relevant.size());

          // enumerate locally relevant vertices
          unsigned int vertex_counter = 0;
          for (unsigned int i = 0; i < vertices_locally_relevant.size(); ++i)
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

        void
        add_vertices(DescriptionTemp<dim, spacedim> &construction_data,
                     const std::vector<bool> &vertices_locally_relevant) const
        {
          for (unsigned int i = 0; i < vertices_locally_relevant.size(); ++i)
            if (vertices_locally_relevant[i])
              construction_data.coarse_cell_vertices.emplace_back(
                i, tria.get_vertices()[i]);
        }


        const dealii::Triangulation<dim, spacedim> &tria;
        const std::function<types::subdomain_id(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &)>
          subdomain_id_function;
        const std::function<types::subdomain_id(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &)>
          level_subdomain_id_function;

        const MPI_Comm &                         comm;
        const TriangulationDescription::Settings settings;
        const bool                               construct_multigrid;

        std::map<unsigned int, std::vector<unsigned int>>
                                             coinciding_vertex_groups;
        std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;
      };

    } // namespace


    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const dealii::Triangulation<dim, spacedim> &tria,
      const MPI_Comm &                            comm,
      const TriangulationDescription::Settings    settings,
      const unsigned int                          my_rank_in)
    {
      if (const auto tria_pdt = dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> *>(&tria))
        Assert(comm == tria_pdt->get_communicator(),
               ExcMessage("MPI communicators do not match."));

      // First, figure out for what rank we are supposed to build the
      // TriangulationDescription::Description object
      unsigned int my_rank = my_rank_in;
      Assert(my_rank == numbers::invalid_unsigned_int ||
               my_rank < dealii::Utilities::MPI::n_mpi_processes(comm),
             ExcMessage("Rank has to be smaller than available processes."));

      if (auto tria_pdt = dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> *>(&tria))
        {
          Assert(my_rank == numbers::invalid_unsigned_int ||
                   my_rank == dealii::Utilities::MPI::this_mpi_process(comm),
                 ExcMessage(
                   "For creation from a parallel::distributed::Triangulation, "
                   "my_rank has to equal global rank."));

          my_rank = dealii::Utilities::MPI::this_mpi_process(comm);
        }
      else if (auto tria_serial =
                 dynamic_cast<const dealii::Triangulation<dim, spacedim> *>(
                   &tria))
        {
          if (my_rank == numbers::invalid_unsigned_int)
            my_rank = dealii::Utilities::MPI::this_mpi_process(comm);
        }
      else
        {
          Assert(false,
                 ExcMessage("This type of triangulation is not supported!"));
        }

      const auto subdomain_id_function = [](const auto &cell) {
        return cell->subdomain_id();
      };

      const auto level_subdomain_id_function = [](const auto &cell) {
        return cell->level_subdomain_id();
      };

      return CreateDescriptionFromTriangulationHelper<dim, spacedim>(
               tria,
               subdomain_id_function,
               level_subdomain_id_function,
               comm,
               settings)
        .template create_description_for_rank<Description<dim, spacedim>>(
          my_rank);
    }



    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation_in_groups(
      const std::function<void(dealii::Triangulation<dim, spacedim> &)>
        &                                            serial_grid_generator,
      const std::function<void(dealii::Triangulation<dim, spacedim> &,
                               const MPI_Comm &,
                               const unsigned int)> &serial_grid_partitioner,
      const MPI_Comm &                               comm,
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
      const Triangulation<dim, spacedim> &              tria,
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

      if (construct_multigrid) // perform first child policy
        {
          const auto tria_parallel =
            dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &tria);

          Assert(tria_parallel, ExcInternalError());

          partition.update_ghost_values();

          partitions_mg.resize(tria.n_global_levels());

          for (unsigned int l = 0; l < tria.n_global_levels(); ++l)
            partitions_mg[l].reinit(
              tria_parallel->global_level_cell_index_partitioner(l).lock());

          for (int level = tria.n_global_levels() - 1; level >= 0; --level)
            {
              for (const auto &cell : tria.cell_iterators_on_level(level))
                {
                  if (cell->is_locally_owned_on_level() == false)
                    continue;

                  if (cell->is_active())
                    partitions_mg[level][cell->global_level_cell_index()] =
                      partition[cell->global_active_cell_index()];
                  else
                    partitions_mg[level][cell->global_level_cell_index()] =
                      partitions_mg[level + 1]
                                   [cell->child(0)->global_level_cell_index()];
                }

              partitions_mg[level].update_ghost_values();
            }
        }

      return create_description_from_triangulation(tria,
                                                   partition,
                                                   partitions_mg,
                                                   settings);
    }



    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const Triangulation<dim, spacedim> &              tria,
      const LinearAlgebra::distributed::Vector<double> &partition,
      const std::vector<LinearAlgebra::distributed::Vector<double>>
        &                                      partitions_mg,
      const TriangulationDescription::Settings settings_in)
    {
#ifdef DEAL_II_WITH_MPI
      if (tria.get_communicator() == MPI_COMM_NULL)
        AssertDimension(partition.local_size(), 0);
#endif

      if (partition.size() == 0)
        {
          AssertDimension(partitions_mg.size(), 0);
          return create_description_from_triangulation(tria,
                                                       tria.get_communicator(),
                                                       settings_in);
        }

      partition.update_ghost_values();

      for (const auto &partition : partitions_mg)
        partition.update_ghost_values();

      // 1) determine processes owning locally owned cells
      const std::vector<unsigned int> relevant_processes = [&]() {
        std::set<unsigned int> relevant_processes;

        for (unsigned int i = 0; i < partition.local_size(); ++i)
          relevant_processes.insert(
            static_cast<unsigned int>(partition.local_element(i)));

        for (const auto &partition : partitions_mg)
          for (unsigned int i = 0; i < partition.local_size(); ++i)
            relevant_processes.insert(
              static_cast<unsigned int>(partition.local_element(i)));

        return std::vector<unsigned int>(relevant_processes.begin(),
                                         relevant_processes.end());
      }();

      const bool construct_multigrid = partitions_mg.size() > 0;

      TriangulationDescription::Settings settings = settings_in;

      if (construct_multigrid)
        settings = static_cast<TriangulationDescription::Settings>(
          settings |
          TriangulationDescription::Settings::construct_multigrid_hierarchy);

      const auto subdomain_id_function = [&partition](const auto &cell) {
        if ((cell->is_active() && (cell->is_artificial() == false)))
          return static_cast<unsigned int>(
            partition[cell->global_active_cell_index()]);
        else
          return numbers::artificial_subdomain_id;
      };

      const auto level_subdomain_id_function =
        [&construct_multigrid, &partitions_mg](const auto &cell) {
          if (construct_multigrid && (cell->is_artificial_on_level() == false))
            return static_cast<unsigned int>(
              partitions_mg[cell->level()][cell->global_level_cell_index()]);
          else
            return numbers::artificial_subdomain_id;
        };

      CreateDescriptionFromTriangulationHelper<dim, spacedim> helper(
        tria,
        subdomain_id_function,
        level_subdomain_id_function,
        tria.get_communicator(),
        settings);

      // create a description (locally owned cell and a layer of ghost cells
      // and all their parents)
      std::vector<DescriptionTemp<dim, spacedim>> descriptions_per_rank;
      descriptions_per_rank.reserve(relevant_processes.size());

      for (const auto rank : relevant_processes)
        descriptions_per_rank.emplace_back(
          helper.template create_description_for_rank<
            DescriptionTemp<dim, spacedim>>(rank));

      // collect description from all processes that used to own locally-owned
      // active cells of this process in a single description
      DescriptionTemp<dim, spacedim> description_merged;
      description_merged.collect(
        relevant_processes,
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
#include "tria_description.inst"


DEAL_II_NAMESPACE_CLOSE
