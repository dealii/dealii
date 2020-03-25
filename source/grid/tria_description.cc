// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

DEAL_II_NAMESPACE_OPEN

namespace TriangulationDescription
{
  namespace Utilities
  {
    namespace
    {
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
       * Convert the binary representation of a CellId to coarse-cell id as
       * if the finest level were the coarsest level ("level coarse-grid id").
       */
      template <int dim>
      types::coarse_cell_id
      convert_cell_id_binary_type_to_level_coarse_cell_id(
        const typename CellId::binary_type &binary_representation)
      {
        // exploiting the structure of CellId::binary_type
        // see also the documentation of CellId

        // actual coarse-grid id
        const unsigned int coarse_cell_id  = binary_representation[0];
        const unsigned int n_child_indices = binary_representation[1] >> 2;

        const unsigned int children_per_value =
          sizeof(CellId::binary_type::value_type) * 8 / dim;
        unsigned int child_level  = 0;
        unsigned int binary_entry = 2;

        // path to the get to the cell
        std::vector<unsigned int> cell_indices;
        while (child_level < n_child_indices)
          {
            Assert(binary_entry < binary_representation.size(),
                   ExcInternalError());

            for (unsigned int j = 0; j < children_per_value; ++j)
              {
                unsigned int cell_index =
                  (((binary_representation[binary_entry] >> (j * dim))) &
                   (GeometryInfo<dim>::max_children_per_cell - 1));
                cell_indices.push_back(cell_index);
                ++child_level;
                if (child_level == n_child_indices)
                  break;
              }
            ++binary_entry;
          }

        // compute new coarse-grid id: c_{i+1} = c_{i}*2^dim + q;
        types::coarse_cell_id level_coarse_cell_id = coarse_cell_id;
        for (auto i : cell_indices)
          level_coarse_cell_id =
            level_coarse_cell_id * GeometryInfo<dim>::max_children_per_cell + i;

        return level_coarse_cell_id;
      }
    } // namespace


    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const dealii::Triangulation<dim, spacedim> &tria,
      const MPI_Comm                              comm,
      const bool         construct_multilevel_hierarchy,
      const unsigned int my_rank_in)
    {
      if (auto tria_pdt = dynamic_cast<
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
          Assert(
            my_rank == numbers::invalid_unsigned_int ||
              my_rank == dealii::Utilities::MPI::this_mpi_process(comm),
            ExcMessage(
              "If parallel::distributed::Triangulation as source triangulation, my_rank has to equal global rank."));

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

      Description<dim, spacedim> construction_data;

      // store the communicator
      construction_data.comm = comm;


      std::map<unsigned int, std::vector<unsigned int>>
                                           coinciding_vertex_groups;
      std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;

      GridTools::collect_coinciding_vertices(tria,
                                             coinciding_vertex_groups,
                                             vertex_to_coinciding_vertex_group);

      // helper function, which collects all vertices belonging to a cell
      // (also taking into account periodicity)
      auto add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells =
        [coinciding_vertex_groups, vertex_to_coinciding_vertex_group](
          TriaIterator<CellAccessor<dim, spacedim>> &cell,
          std::vector<bool> &vertices_owned_by_locally_owned_cells) {
          // add local vertices
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            {
              vertices_owned_by_locally_owned_cells[cell->vertex_index(v)] =
                true;
              const auto coinciding_vertex_group =
                vertex_to_coinciding_vertex_group.find(cell->vertex_index(v));
              if (coinciding_vertex_group !=
                  vertex_to_coinciding_vertex_group.end())
                for (const auto &co_vertex : coinciding_vertex_groups.at(
                       coinciding_vertex_group->second))
                  vertices_owned_by_locally_owned_cells[co_vertex] = true;
            }
        };

      // check if multilevel hierarchy should be constructed
      if (construct_multilevel_hierarchy == false)
        {
          Assert(
            tria.has_hanging_nodes() == false,
            ExcMessage(
              "Hanging nodes are only supported if multilevel hierarchy is constructed!"));

          construction_data.settings =
            TriangulationDescription::Settings::default_setting;

          // 1) collect vertices of active locally owned cells
          std::vector<bool> vertices_owned_by_locally_owned_cells(
            tria.n_vertices());
          for (auto cell : tria.cell_iterators())
            if (cell->is_active() && cell->subdomain_id() == my_rank)
              add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                cell, vertices_owned_by_locally_owned_cells);

          // helper function to determine if cell is locally relevant
          // (i.e. a cell which is connected to a vertex via a locally owned
          // active cell)
          auto is_locally_relevant = [&](
                                       TriaIterator<CellAccessor<dim, spacedim>>
                                         &cell) {
            for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
              if (vertices_owned_by_locally_owned_cells[cell->vertex_index(v)])
                return true;
            return false;
          };

          // 2) process all local and ghost cells: set up needed data
          // structures and collect all locally relevant vertices
          // for second sweep
          std::map<unsigned int, unsigned int> vertices_locally_relevant;
          construction_data.cell_infos.resize(1);

          for (auto cell : tria.cell_iterators())
            if (cell->is_active() && is_locally_relevant(cell))
              {
                // to be filled
                CellData<dim> cell_info;

                // a) extract cell definition (with old numbering of
                // vertices)
                dealii::CellData<dim> cell_data;
                cell_data.material_id = cell->material_id();
                cell_data.manifold_id = cell->manifold_id();
                for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                  cell_data.vertices[v] = cell->vertex_index(v);
                construction_data.coarse_cells.push_back(cell_data);

                // b) save indices of each vertex of this cell
                for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                  vertices_locally_relevant[cell->vertex_index(v)] =
                    numbers::invalid_unsigned_int;

                // c) save boundary_ids of each face of this cell
                for (const unsigned int f : GeometryInfo<dim>::face_indices())
                  {
                    types::boundary_id boundary_ind =
                      cell->face(f)->boundary_id();
                    if (boundary_ind != numbers::internal_face_boundary_id)
                      cell_info.boundary_ids.emplace_back(f, boundary_ind);
                  }

                // d) compute a new coarse-cell id (by ignoring all parent
                // level)
                types::coarse_cell_id new_coarse_cell_id =
                  convert_cell_id_binary_type_to_level_coarse_cell_id<dim>(
                    cell->id().template to_binary<dim>());

                //    store the coarse cell id
                cell_info.id = CellId(new_coarse_cell_id, 0, nullptr)
                                 .template to_binary<dim>();

                //    save coarse_cell_index -> coarse_cell_id mapping
                construction_data.coarse_cell_index_to_coarse_cell_id.push_back(
                  new_coarse_cell_id);

                // e) store manifold id of cell
                cell_info.manifold_id = cell->manifold_id();

                // ... of lines
                if (dim >= 2)
                  for (unsigned int line = 0;
                       line < GeometryInfo<dim>::lines_per_cell;
                       ++line)
                    cell_info.manifold_line_ids[line] =
                      cell->line(line)->manifold_id();

                // ... of quads
                if (dim == 3)
                  for (unsigned int quad = 0;
                       quad < GeometryInfo<dim>::quads_per_cell;
                       ++quad)
                    cell_info.manifold_quad_ids[quad] =
                      cell->quad(quad)->manifold_id();

                // f) store subdomain_id
                cell_info.subdomain_id = cell->subdomain_id();

                // g) store invalid level_subdomain_id (since multilevel
                //    hierarchy is not constructed)
                cell_info.level_subdomain_id = numbers::invalid_subdomain_id;

                construction_data.cell_infos[0].emplace_back(cell_info);
              }

          // 3) enumerate locally relevant vertices
          unsigned int vertex_counter = 0;
          for (auto &vertex : vertices_locally_relevant)
            {
              construction_data.coarse_cell_vertices.push_back(
                tria.get_vertices()[vertex.first]);
              vertex.second = vertex_counter++;
            }

          // 4) correct vertices of cells (make them local)
          for (auto &cell : construction_data.coarse_cells)
            for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
              cell.vertices[v] = vertices_locally_relevant[cell.vertices[v]];
        }
      else
        {
          construction_data.settings =
            TriangulationDescription::Settings::construct_multigrid_hierarchy;

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
              for (auto cell : tria.cell_iterators_on_level(level))
                if (cell->level_subdomain_id() == my_rank ||
                    (cell->active() && cell->subdomain_id() == my_rank))
                  add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                    cell, vertices_owned_by_locally_owned_cells_on_level);

              for (auto cell : tria.active_cell_iterators())
                if (cell->subdomain_id() == my_rank)
                  add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                    cell, vertices_owned_by_locally_owned_cells_on_level);

              // helper function to determine if cell is locally relevant
              // (i.e. a cell which is connected to a vertex via a locally
              // owned cell)
              auto is_locally_relevant_on_level =
                [&](TriaIterator<CellAccessor<dim, spacedim>> &cell) {
                  for (const unsigned int v :
                       GeometryInfo<dim>::vertex_indices())
                    if (vertices_owned_by_locally_owned_cells_on_level
                          [cell->vertex_index(v)])
                      return true;
                  return false;
                };

              // mark all locally relevant cells
              for (auto cell : tria.cell_iterators_on_level(level))
                if (is_locally_relevant_on_level(cell))
                  set_user_flag_and_of_its_parents(cell);
            }

          // 2) set_up coarse-grid triangulation
          {
            std::map<unsigned int, unsigned int> vertices_locally_relevant;

            // a) loop over all cells
            for (auto cell : tria.cell_iterators_on_level(0))
              {
                if (!cell->user_flag_set())
                  continue;

                // extract cell definition (with old numbering of vertices)
                dealii::CellData<dim> cell_data;
                cell_data.material_id = cell->material_id();
                cell_data.manifold_id = cell->manifold_id();
                for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                  cell_data.vertices[v] = cell->vertex_index(v);
                construction_data.coarse_cells.push_back(cell_data);

                // save indices of each vertex of this cell
                for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                  vertices_locally_relevant[cell->vertex_index(v)] =
                    numbers::invalid_unsigned_int;

                // save translation for corase grid: lid -> gid
                construction_data.coarse_cell_index_to_coarse_cell_id.push_back(
                  cell->id().get_coarse_cell_id());
              }

            // b) enumerate locally relevant vertices
            unsigned int vertex_counter = 0;
            for (auto &vertex : vertices_locally_relevant)
              {
                construction_data.coarse_cell_vertices.push_back(
                  tria.get_vertices()[vertex.first]);
                vertex.second = vertex_counter++;
              }

            // c) correct vertices of cells (make them local)
            for (auto &cell : construction_data.coarse_cells)
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                cell.vertices[v] = vertices_locally_relevant[cell.vertices[v]];
          }


          // 3) collect info of each cell
          construction_data.cell_infos.resize(
            tria.get_triangulation().n_global_levels());

          // collect local vertices on active level
          std::vector<bool> vertices_owned_by_locally_owned_active_cells(
            tria.n_vertices());
          for (auto cell : tria.active_cell_iterators())
            if (cell->subdomain_id() == my_rank)
              add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                cell, vertices_owned_by_locally_owned_active_cells);

          // helper function to determine if cell is locally relevant
          // on active level
          auto is_locally_relevant_on_active_level =
            [&](TriaIterator<CellAccessor<dim, spacedim>> &cell) {
              if (cell->active())
                for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
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
              for (auto cell : tria.cell_iterators_on_level(level))
                if (cell->level_subdomain_id() == my_rank ||
                    (cell->active() && cell->subdomain_id() == my_rank))
                  add_vertices_of_cell_to_vertices_owned_by_locally_owned_cells(
                    cell, vertices_owned_by_locally_owned_cells_on_level);

              // helper function to determine if cell is locally relevant
              // on level
              auto is_locally_relevant_on_level =
                [&](TriaIterator<CellAccessor<dim, spacedim>> &cell) {
                  for (const unsigned int v :
                       GeometryInfo<dim>::vertex_indices())
                    if (vertices_owned_by_locally_owned_cells_on_level
                          [cell->vertex_index(v)])
                      return true;
                  return false;
                };

              auto &level_cell_infos = construction_data.cell_infos[level];
              for (auto cell : tria.cell_iterators_on_level(level))
                {
                  // check if cell is locally relevant
                  if (!(cell->user_flag_set()))
                    continue;

                  CellData<dim> cell_info;

                  // save coarse-cell id
                  cell_info.id = cell->id().template to_binary<dim>();

                  // save boundary_ids of each face of this cell
                  for (const unsigned int f : GeometryInfo<dim>::face_indices())
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
                      for (unsigned int line = 0;
                           line < GeometryInfo<dim>::lines_per_cell;
                           ++line)
                        cell_info.manifold_line_ids[line] =
                          cell->line(line)->manifold_id();

                    // ... of quads
                    if (dim == 3)
                      for (unsigned int quad = 0;
                           quad < GeometryInfo<dim>::quads_per_cell;
                           ++quad)
                        cell_info.manifold_quad_ids[quad] =
                          cell->quad(quad)->manifold_id();
                  }

                  // subdomain and level subdomain id
                  cell_info.subdomain_id = numbers::artificial_subdomain_id;
                  cell_info.level_subdomain_id =
                    numbers::artificial_subdomain_id;

                  if (is_locally_relevant_on_active_level(cell))
                    {
                      cell_info.level_subdomain_id = cell->level_subdomain_id();
                      cell_info.subdomain_id       = cell->subdomain_id();
                    }
                  else if (is_locally_relevant_on_level(cell))
                    {
                      cell_info.level_subdomain_id = cell->level_subdomain_id();
                    }
                  // else: cell is locally relevant but an artificial cell

                  level_cell_infos.emplace_back(cell_info);
                }
            }

          const_cast<dealii::Triangulation<dim, spacedim> &>(tria)
            .load_user_flags(old_user_flags);
        }

      return construction_data;
    }



    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation_in_groups(
      const std::function<void(dealii::Triangulation<dim, spacedim> &)>
        &                                            serial_grid_generator,
      const std::function<void(dealii::Triangulation<dim, spacedim> &,
                               const MPI_Comm,
                               const unsigned int)> &serial_grid_partitioner,
      const MPI_Comm                                 comm,
      const int                                      group_size,
      const bool construct_multilevel_hierarchy)
    {
#ifndef DEAL_II_WITH_MPI
      (void)serial_grid_generator;
      (void)serial_grid_partitioner;
      (void)comm;
      (void)group_size;
      (void)construct_multilevel_hierarchy;

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
            construct_multilevel_hierarchy ?
              dealii::Triangulation<dim, spacedim>::none :
              dealii::Triangulation<dim, spacedim>::
                limit_level_difference_at_vertices);
          serial_grid_generator(tria);

          // Step 2: partition active cells and ...
          serial_grid_partitioner(tria, comm, group_size);

          // ... cells on the levels if multigrid is required
          if (construct_multilevel_hierarchy)
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
                create_description_from_triangulation(
                  tria, comm, construct_multilevel_hierarchy, other_rank);
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
          return create_description_from_triangulation(
            tria, comm, construct_multilevel_hierarchy, my_rank);
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

  } // namespace Utilities
} // namespace TriangulationDescription



/*-------------- Explicit Instantiations -------------------------------*/
#include "tria_description.inst"


DEAL_II_NAMESPACE_CLOSE
