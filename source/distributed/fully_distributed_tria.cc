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


#include <deal.II/base/memory_consumption.h>

#include <deal.II/grid/grid_tools.h>
//#include <deal.II/base/std_cxx17/optional.h>

#include <deal.II/distributed/fully_distributed_tria.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
namespace GridGenerator
{
  template <int dim, int spacedim>
  void
  hyper_cube(Triangulation<dim, spacedim> &tria,
             const double                  left,
             const double                  right,
             const bool                    colorize);
} // namespace GridGenerator

// namespace GridTools
//{
// template<typename DataType , typename MeshType >
// void exchange_cell_data_to_ghosts	(	const MeshType & 	mesh,
// const std::function< std_cxx17::optional< DataType >(const typename
// MeshType::active_cell_iterator &)> & 	pack, const std::function< void(const
// typename MeshType::active_cell_iterator &, const DataType &)> & 	unpack )	;
//}

namespace parallel
{
  namespace fullydistributed
  {
    template <int dim, int spacedim>
    Triangulation<dim, spacedim>::Triangulation(MPI_Comm mpi_communicator)
      : parallel::DistributedTriangulationBase<dim, spacedim>(mpi_communicator)
      , settings(Settings::default_setting)
      , currently_processing_create_triangulation_for_internal_usage(false)
      , currently_processing_prepare_coarsening_and_refinement_for_internal_usage(
          false)
    {}



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::create_triangulation(
      const ConstructionData<dim, spacedim> &construction_data)
    {
      // check if the communicator of this parallel triangulation has been used
      // to construct the ConstructionData
      Assert(construction_data.comm == this->mpi_communicator,
             ExcMessage("MPI communicators do not match!"));

      // store internally the settings
      settings = construction_data.settings;

      // set the smoothing properties
      if (settings & construct_multigrid_hierarchy)
        this->set_mesh_smoothing(
          static_cast<
            typename dealii::Triangulation<dim, spacedim>::MeshSmoothing>(
            dealii::Triangulation<dim>::none |
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices));
      else
        this->set_mesh_smoothing(
          static_cast<
            typename dealii::Triangulation<dim, spacedim>::MeshSmoothing>(
            dealii::Triangulation<dim>::none));

      // clear internal data structures
      this->coarse_cell_id_to_coarse_cell_index_vector.clear();
      this->coarse_cell_index_to_coarse_cell_id_vector.clear();

      // check if no locally relevant coarse-grid cells have been provided
      if (construction_data.coarse_cell_vertices.empty())
        {
          // 1) create a dummy hypercube
          currently_processing_create_triangulation_for_internal_usage = true;
          GridGenerator::hyper_cube(*this, 0, 1, false);
          currently_processing_create_triangulation_for_internal_usage = false;

          // 2) mark cell as artificial
          auto cell = this->begin();
          cell->set_subdomain_id(dealii::numbers::artificial_subdomain_id);
          cell->set_level_subdomain_id(
            dealii::numbers::artificial_subdomain_id);

          // 3) set up dummy mapping between locally relevant coarse-grid cells
          //    and global cells
          this->coarse_cell_id_to_coarse_cell_index_vector.emplace_back(
            numbers::invalid_coarse_cell_id, 0);
          this->coarse_cell_index_to_coarse_cell_id_vector.emplace_back(
            numbers::invalid_coarse_cell_id);
        }
      else
        {
          // 1) store `coarse-cell index to coarse-cell id`-mapping
          this->coarse_cell_index_to_coarse_cell_id_vector =
            construction_data.coarse_cell_index_to_coarse_cell_id;

          // 2) set up `coarse-cell id to coarse-cell index`-mapping
          std::map<types::coarse_cell_id, unsigned int>
            coarse_cell_id_to_coarse_cell_index_vector;
          for (unsigned int i = 0;
               i < construction_data.coarse_cell_index_to_coarse_cell_id.size();
               ++i)
            coarse_cell_id_to_coarse_cell_index_vector
              [construction_data.coarse_cell_index_to_coarse_cell_id[i]] = i;

          for (auto i : coarse_cell_id_to_coarse_cell_index_vector)
            this->coarse_cell_id_to_coarse_cell_index_vector.emplace_back(i);

          // 3) create coarse grid
          dealii::parallel::Triangulation<dim, spacedim>::create_triangulation(
            construction_data.coarse_cell_vertices,
            construction_data.coarse_cells,
            SubCellData());

          Assert(this->n_cells() ==
                   this->coarse_cell_id_to_coarse_cell_index_vector.size(),
                 ExcInternalError());
          Assert(this->n_cells() ==
                   this->coarse_cell_index_to_coarse_cell_id_vector.size(),
                 ExcInternalError());

          // create a copy of cell_infos such that we can sort them
          auto cell_infos = construction_data.cell_infos;

          // sort cell_infos on each level separately
          for (auto &cell_info : cell_infos)
            std::sort(cell_info.begin(),
                      cell_info.end(),
                      [&](CellData<dim> a, CellData<dim> b) {
                        const CellId a_id(a.id);
                        const CellId b_id(b.id);

                        const auto a_coarse_cell_index =
                          this->coarse_cell_id_to_coarse_cell_index(
                            a_id.get_coarse_cell_id());
                        const auto b_coarse_cell_index =
                          this->coarse_cell_id_to_coarse_cell_index(
                            b_id.get_coarse_cell_id());

                        // according to their coarse-cell index and if that is
                        // same according to their cell id (the result is that
                        // cells on each level are sorted according to their
                        // index on that level - what we need in the following
                        // operations)
                        if (a_coarse_cell_index != b_coarse_cell_index)
                          return a_coarse_cell_index < b_coarse_cell_index;
                        else
                          return a_id < b_id;
                      });

          // 4) create all levels via a sequence of refinements
          for (unsigned int level = 0; level < cell_infos.size(); ++level)
            {
              // a) set manifold ids here (because new vertices have to be
              //    positioned correctly during each refinement step)
              {
                auto cell      = this->begin(level);
                auto cell_info = cell_infos[level].begin();
                for (; cell_info != cell_infos[level].end(); ++cell_info)
                  {
                    while (cell_info->id !=
                           cell->id().template to_binary<dim>())
                      ++cell;
                    if (spacedim == 3)
                      for (unsigned int quad = 0;
                           quad < GeometryInfo<spacedim>::quads_per_cell;
                           ++quad)
                        cell->quad(quad)->set_manifold_id(
                          cell_info->manifold_quad_ids[quad]);

                    if (spacedim >= 2)
                      for (unsigned int line = 0;
                           line < GeometryInfo<spacedim>::lines_per_cell;
                           ++line)
                        cell->line(line)->set_manifold_id(
                          cell_info->manifold_line_ids[line]);

                    cell->set_manifold_id(cell_info->manifold_id);
                  }
              }

              // b) perform refinement on all levels but on the finest
              if (level + 1 != cell_infos.size())
                {
                  // find cells that should have children and mark them for
                  // refinement
                  auto coarse_cell    = this->begin(level);
                  auto fine_cell_info = cell_infos[level + 1].begin();

                  // loop over all cells on the next level
                  for (; fine_cell_info != cell_infos[level + 1].end();
                       ++fine_cell_info)
                    {
                      // find the parent of that cell
                      while (!coarse_cell->id().is_parent_of(
                        CellId(fine_cell_info->id)))
                        ++coarse_cell;

                      // set parent for refinement
                      coarse_cell->set_refine_flag();
                    }

                  // execute refinement
                  currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
                    true;
                  dealii::Triangulation<dim, spacedim>::
                    execute_coarsening_and_refinement();
                  currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
                    false;
                }
            }

          // 4a) set all cells artificial (and set the actual
          //     (level_)subdomain_ids in the next step)
          for (auto cell = this->begin(); cell != this->end(); ++cell)
            {
              if (cell->active())
                cell->set_subdomain_id(
                  dealii::numbers::artificial_subdomain_id);

              cell->set_level_subdomain_id(
                dealii::numbers::artificial_subdomain_id);
            }

          // 4b) set actual (level_)subdomain_ids as well as boundary ids
          for (unsigned int level = 0; level < cell_infos.size(); ++level)
            {
              auto cell      = this->begin(level);
              auto cell_info = cell_infos[level].begin();
              for (; cell_info != cell_infos[level].end(); ++cell_info)
                {
                  // find cell that has the correct cell
                  while (cell_info->id != cell->id().template to_binary<dim>())
                    ++cell;

                  // subdomain id
                  if (cell->active())
                    cell->set_subdomain_id(cell_info->subdomain_id);

                  // level subdomain id
                  if (settings & construct_multigrid_hierarchy)
                    cell->set_level_subdomain_id(cell_info->level_subdomain_id);

                  // boundary ids
                  for (auto pair : cell_info->boundary_ids)
                    {
                      Assert(cell->at_boundary(pair.first),
                             ExcMessage("Cell face is not on the boundary!"));
                      cell->face(pair.first)->set_boundary_id(pair.second);
                    }
                }
            }
        }

      update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::create_triangulation(
      const std::vector<Point<spacedim>> &      vertices,
      const std::vector<dealii::CellData<dim>> &cells,
      const SubCellData &                       subcelldata)
    {
      AssertThrow(
        currently_processing_create_triangulation_for_internal_usage,
        ExcMessage(
          "Use the other create_triangulation() function to create triangulations of type parallel::fullydistributed::Triangulation.!"));

      dealii::Triangulation<dim, spacedim>::create_triangulation(vertices,
                                                                 cells,
                                                                 subcelldata);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::update_number_cache()
    {
      if (settings & construct_multigrid_hierarchy)
        parallel::Triangulation<dim, spacedim>::fill_level_ghost_owners();

      parallel::Triangulation<dim, spacedim>::update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_coarsening_and_refinement()
    {
      // 1) preparation
      {
        // reset the coarsening/refinement flags of non-local cells
        for (auto cell : this->cell_iterators())
          if (cell->active() && !cell->is_locally_owned())
            {
              cell->clear_refine_flag();
              cell->clear_coarsen_flag();

              // start assumption:
              if (cell->level() > 0)
                cell->set_coarsen_flag();
            }

        currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
          true;
        this->prepare_coarsening_and_refinement();
        currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
          false;


        // update refinement/coarsening flags of non-local cells via
        // ghost-cell exchange
        GridTools::exchange_cell_data_to_ghosts<std::pair<bool, bool>>(
          *this,
          [](const TriaIterator<CellAccessor<dim, spacedim>> &cell) {
            if (cell->refine_flag_set())
              {
                Assert(cell->refine_flag_set() ==
                         RefinementPossibilities<dim>::isotropic_refinement,
                       ExcMessage(
                         "This class does not support anisotropic refinement"));
                return std::pair<bool, bool>(true, false);
              }
            else if (!cell->coarsen_flag_set())
              return std::pair<bool, bool>(false, true);
            else
              return std::pair<bool, bool>(false, false);
          },
          [](const TriaIterator<CellAccessor<dim, spacedim>> &cell,
             const std::pair<bool, bool>                      refinement_flag) {
            if (refinement_flag.first == true)
              {
                cell->clear_coarsen_flag();
                cell->set_refine_flag();
              }
            else if (refinement_flag.second == true)
              cell->clear_coarsen_flag();
          });
      }

      // 1b) store current active cells
      std::set<typename CellId::binary_type> active_cells;
      for (auto cell : this->active_cell_iterators())
        active_cells.insert(cell->id().template to_binary<dim>());

      // 2) execute actual coarsening and/or refinement
      currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
        true;
      dealii::Triangulation<dim, spacedim>::execute_coarsening_and_refinement();
      currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
        false;

      // 3a) post-processing:
      for (auto cell : this->active_cell_iterators())
        if (cell->level() > 0 &&
            (active_cells.find(
               cell->parent()->id().template to_binary<dim>()) !=
             active_cells.end()))
          {
            // cell has been refined
            cell->set_subdomain_id(cell->parent()->level_subdomain_id());
            cell->set_level_subdomain_id(cell->parent()->level_subdomain_id());
          }
        else if (active_cells.find(cell->id().template to_binary<dim>()) !=
                 active_cells.end())
          {
            // nothing has happened with this cell: nothing to do (?)
            cell->set_subdomain_id(cell->level_subdomain_id());
          }
        else
          {
            // cell has become active, since it children have been coarsened
            // make level_subdomain_id to subdomain_id
            cell->set_subdomain_id(cell->level_subdomain_id());
          }

      // 3b) update artificial cells (keep only one ghost layer)
      {
        // collect vertices belonging to active local cells
        auto add_vertices_of_cell_to_vertices_owned_by_loclly_owned_cells =
          [](TriaIterator<CellAccessor<dim, spacedim>> &cell,
             std::vector<bool> &vertices_owned_by_loclly_owned_cells) {
            // add vertices belonging to a periodic neighbor
            for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; i++)
              if (cell->has_periodic_neighbor(i))
                {
                  const auto face_t = cell->face(i);
                  const auto face_n = cell->periodic_neighbor(i)->face(
                    cell->periodic_neighbor_face_no(i));
                  for (unsigned int j = 0;
                       j < GeometryInfo<dim>::vertices_per_face;
                       j++)
                    {
                      vertices_owned_by_loclly_owned_cells[face_t->vertex_index(
                        j)] = true;
                      vertices_owned_by_loclly_owned_cells[face_n->vertex_index(
                        j)] = true;
                    }
                }

            // add local vertices
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                 v++)
              vertices_owned_by_loclly_owned_cells[cell->vertex_index(v)] =
                true;
          };

        std::vector<bool> vertices_owned_by_loclly_owned_cells(
          this->n_vertices());
        for (auto cell : this->cell_iterators())
          if (cell->active() &&
              (cell->is_locally_owned() ||
               cell->level_subdomain_id() == this->locally_owned_subdomain()))
            add_vertices_of_cell_to_vertices_owned_by_loclly_owned_cells(
              cell, vertices_owned_by_loclly_owned_cells);

        auto is_locally_relevant = [&](auto &cell) {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               v++)
            if (vertices_owned_by_loclly_owned_cells[cell->vertex_index(v)])
              return true;
          return false;
        };

        // mark all active cells which are not connected via a vertex to
        // an active local cells as artificial
        for (auto cell : this->cell_iterators())
          if (cell->active() && !is_locally_relevant(cell))
            cell->set_subdomain_id(numbers::artificial_subdomain_id);
      }

      // 3c) update cache
      update_number_cache();
    } // namespace fullydistributed



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::prepare_coarsening_and_refinement()
    {
      Assert(
        currently_processing_prepare_coarsening_and_refinement_for_internal_usage,
        ExcMessage("No coarsening and refinement is supported!"));

      return dealii::Triangulation<dim, spacedim>::
        prepare_coarsening_and_refinement();
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::has_hanging_nodes() const
    {
      Assert(false, ExcNotImplemented());
      return false;
    }



    template <int dim, int spacedim>
    std::size_t
    Triangulation<dim, spacedim>::memory_consumption() const
    {
      const std::size_t mem =
        this->dealii::parallel::TriangulationBase<dim, spacedim>::
          memory_consumption() +
        MemoryConsumption::memory_consumption(
          coarse_cell_id_to_coarse_cell_index_vector) +
        MemoryConsumption::memory_consumption(
          coarse_cell_index_to_coarse_cell_id_vector);
      return mem;
    }



    template <int dim, int spacedim>
    std::map<unsigned int, std::set<dealii::types::subdomain_id>>
    Triangulation<dim, spacedim>::compute_vertices_with_ghost_neighbors() const
    {
      // 1) collect for each vertex on periodic faces all vertices it coincides
      //    with
      std::map<unsigned int, std::set<unsigned int>>
        vertex_to_coinciding_vertices;
      {
        // 1a) collect nodes coinciding due to periodicity
        std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex;
        for (auto &cell : this->active_cell_iterators())
          if (cell->is_locally_owned() || cell->is_ghost())
            for (auto i = 0u; i < GeometryInfo<dim>::faces_per_cell; ++i)
              if (cell->has_periodic_neighbor(i) &&
                  cell->periodic_neighbor(i)->active())
                {
                  auto face_t = cell->face(i);
                  auto face_n = cell->periodic_neighbor(i)->face(
                    cell->periodic_neighbor_face_no(i));
                  for (auto j = 0u; j < GeometryInfo<dim>::vertices_per_face;
                       ++j)
                    {
                      auto         v_t  = face_t->vertex_index(j);
                      auto         v_n  = face_n->vertex_index(j);
                      unsigned int temp = std::min(v_t, v_n);
                      {
                        auto it = vertex_to_coinciding_vertex.find(v_t);
                        if (it != vertex_to_coinciding_vertex.end())
                          temp = std::min(temp, it->second);
                      }
                      {
                        auto it = vertex_to_coinciding_vertex.find(v_n);
                        if (it != vertex_to_coinciding_vertex.end())
                          temp = std::min(temp, it->second);
                      }
                      vertex_to_coinciding_vertex[v_t] = temp;
                      vertex_to_coinciding_vertex[v_n] = temp;
                    }
                }

        // 1b) compress map: let vertices point to the coinciding vertex with
        //     the smallest index
        for (auto &p : vertex_to_coinciding_vertex)
          {
            if (p.first == p.second)
              continue;
            unsigned int temp = p.second;
            while (temp != vertex_to_coinciding_vertex[temp])
              temp = vertex_to_coinciding_vertex[temp];
            p.second = temp;
          }

#ifdef DEBUG
        // check if map is actually compressed
        for (auto p : vertex_to_coinciding_vertex)
          {
            if (p.first == p.second)
              continue;
            auto pp = vertex_to_coinciding_vertex.find(p.second);
            if (pp->first == pp->second)
              continue;
            AssertThrow(false, ExcMessage("Map has to be compressed!"));
          }
#endif

        // 1c) create a map: smallest index of coinciding index -> all
        // coinciding indices
        std::map<unsigned int, std::set<unsigned int>>
          smallest_coinciding_vertex_to_coinciding_vertices;
        for (auto p : vertex_to_coinciding_vertex)
          smallest_coinciding_vertex_to_coinciding_vertices[p.second] =
            std::set<unsigned int>();

        for (auto p : vertex_to_coinciding_vertex)
          smallest_coinciding_vertex_to_coinciding_vertices[p.second].insert(
            p.first);

        // 1d) create a map: vertex -> all coinciding indices
        for (auto &s : smallest_coinciding_vertex_to_coinciding_vertices)
          for (auto &ss : s.second)
            vertex_to_coinciding_vertices[ss] = s.second;
      }

      // 2) collect vertices belonging to local cells
      std::vector<bool> vertex_of_own_cell(this->n_vertices(), false);
      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_locally_owned())
          for (auto v = 0u; v < GeometryInfo<dim>::vertices_per_cell; ++v)
            vertex_of_own_cell[cell->vertex_index(v)] = true;

      // 3) for for each vertex belonging to a locally owned cell all ghost
      //    neighbors (including the periodic own)
      std::map<unsigned int, std::set<dealii::types::subdomain_id>> result;

      // loop over all active ghost cells
      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_ghost())
          {
            const auto owner = cell->subdomain_id();

            // loop over all its vertices
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                 ++v)
              {
                // set owner if vertex belongs to a local cell
                if (vertex_of_own_cell[cell->vertex_index(v)])
                  result[cell->vertex_index(v)].insert(owner);

                // mark also nodes coinciding due to periodicity
                auto coinciding_vertices =
                  vertex_to_coinciding_vertices.find(cell->vertex_index(v));
                if (coinciding_vertices != vertex_to_coinciding_vertices.end())
                  for (auto coinciding_vertex : coinciding_vertices->second)
                    if (vertex_of_own_cell[coinciding_vertex])
                      result[coinciding_vertex].insert(owner);
              }
          }

      return result;
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      return (settings & construct_multigrid_hierarchy);
    }



    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim, spacedim>::coarse_cell_id_to_coarse_cell_index(
      const types::coarse_cell_id coarse_cell_id) const
    {
      const auto coarse_cell_index = std::lower_bound(
        coarse_cell_id_to_coarse_cell_index_vector.begin(),
        coarse_cell_id_to_coarse_cell_index_vector.end(),
        coarse_cell_id,
        [](const std::pair<types::coarse_cell_id, unsigned int> &pair,
           const types::coarse_cell_id &val) { return pair.first < val; });
      Assert(coarse_cell_index !=
               coarse_cell_id_to_coarse_cell_index_vector.cend(),
             ExcMessage("Coarse cell index not found!"));
      return coarse_cell_index->second;
    }



    template <int dim, int spacedim>
    types::coarse_cell_id
    Triangulation<dim, spacedim>::coarse_cell_index_to_coarse_cell_id(
      const unsigned int coarse_cell_index) const
    {
      Assert(coarse_cell_index <
               coarse_cell_index_to_coarse_cell_id_vector.size(),
             ExcIndexRange(coarse_cell_index,
                           0,
                           coarse_cell_index_to_coarse_cell_id_vector.size()));

      const auto coarse_cell_id =
        coarse_cell_index_to_coarse_cell_id_vector[coarse_cell_index];
      Assert(coarse_cell_id != numbers::invalid_coarse_cell_id,
             ExcMessage("You are trying to access a dummy cell!"));
      return coarse_cell_id;
    }



  } // namespace fullydistributed
} // namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "fully_distributed_tria.inst"


DEAL_II_NAMESPACE_CLOSE
