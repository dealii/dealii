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
      parallel::Triangulation<dim, spacedim>::update_number_cache();

      if (settings & construct_multigrid_hierarchy)
        parallel::Triangulation<dim, spacedim>::fill_level_ghost_owners();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_coarsening_and_refinement()
    {
      Assert(false, ExcNotImplemented());
    }



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
      AssertThrow(coarse_cell_id != numbers::invalid_coarse_cell_id,
                  ExcMessage("You are trying to access a dummy cell!"));
      return coarse_cell_id;
    }



  } // namespace fullydistributed
} // namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "fully_distributed_tria.inst"


DEAL_II_NAMESPACE_CLOSE
