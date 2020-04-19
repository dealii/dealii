// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparsity_tools.h>


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_MPI
namespace parallel
{
  namespace shared
  {
    template <int dim, int spacedim>
    Triangulation<dim, spacedim>::Triangulation(
      MPI_Comm mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                     smooth_grid,
      const bool     allow_artificial_cells,
      const Settings settings)
      : dealii::parallel::TriangulationBase<dim, spacedim>(mpi_communicator,
                                                           smooth_grid,
                                                           false)
      , settings(settings)
      , allow_artificial_cells(allow_artificial_cells)
    {
      const auto partition_settings =
        (partition_zoltan | partition_metis | partition_zorder |
         partition_custom_signal) &
        settings;
      (void)partition_settings;
      Assert(partition_settings == partition_auto ||
               partition_settings == partition_metis ||
               partition_settings == partition_zoltan ||
               partition_settings == partition_zorder ||
               partition_settings == partition_custom_signal,
             ExcMessage("Settings must contain exactly one type of the active "
                        "cell partitioning scheme."));

      if (settings & construct_multigrid_hierarchy)
        Assert(allow_artificial_cells,
               ExcMessage("construct_multigrid_hierarchy requires "
                          "allow_artificial_cells to be set to true."));
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      return (settings & construct_multigrid_hierarchy);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::partition()
    {
#  ifdef DEBUG
      // Check that all meshes are the same (or at least have the same
      // total number of active cells):
      const unsigned int max_active_cells =
        Utilities::MPI::max(this->n_active_cells(), this->get_communicator());
      Assert(
        max_active_cells == this->n_active_cells(),
        ExcMessage(
          "A parallel::shared::Triangulation needs to be refined in the same "
          "way on all processors, but the participating processors don't "
          "agree on the number of active cells."));
#  endif

      auto partition_settings = (partition_zoltan | partition_metis |
                                 partition_zorder | partition_custom_signal) &
                                settings;
      if (partition_settings == partition_auto)
#  ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
        partition_settings = partition_zoltan;
#  elif defined DEAL_II_WITH_METIS
        partition_settings = partition_metis;
#  else
        partition_settings = partition_zorder;
#  endif

      if (partition_settings == partition_zoltan)
        {
#  ifndef DEAL_II_TRILINOS_WITH_ZOLTAN
          AssertThrow(false,
                      ExcMessage(
                        "Choosing 'partition_zoltan' requires the library "
                        "to be compiled with support for Zoltan! "
                        "Instead, you might use 'partition_auto' to select "
                        "a partitioning algorithm that is supported "
                        "by your current configuration."));
#  else
          GridTools::partition_triangulation(
            this->n_subdomains, *this, SparsityTools::Partitioner::zoltan);
#  endif
        }
      else if (partition_settings == partition_metis)
        {
#  ifndef DEAL_II_WITH_METIS
          AssertThrow(false,
                      ExcMessage(
                        "Choosing 'partition_metis' requires the library "
                        "to be compiled with support for METIS! "
                        "Instead, you might use 'partition_auto' to select "
                        "a partitioning algorithm that is supported "
                        "by your current configuration."));
#  else
          GridTools::partition_triangulation(this->n_subdomains,
                                             *this,
                                             SparsityTools::Partitioner::metis);
#  endif
        }
      else if (partition_settings == partition_zorder)
        {
          GridTools::partition_triangulation_zorder(this->n_subdomains, *this);
        }
      else if (partition_settings == partition_custom_signal)
        {
          // User partitions mesh manually
        }
      else
        {
          AssertThrow(false, ExcInternalError())
        }

      // do not partition multigrid levels if user is
      // defining a custom partition
      if ((settings & construct_multigrid_hierarchy) &&
          !(settings & partition_custom_signal))
        dealii::GridTools::partition_multigrid_levels(*this);

      true_subdomain_ids_of_cells.resize(this->n_active_cells());

      // loop over all cells and mark artificial:
      typename parallel::shared::Triangulation<dim,
                                               spacedim>::active_cell_iterator
        cell = this->begin_active(),
        endc = this->end();

      if (allow_artificial_cells)
        {
          // get active halo layer of (ghost) cells
          // parallel::shared::Triangulation<dim>::
          std::function<bool(
            const typename parallel::shared::Triangulation<dim, spacedim>::
              active_cell_iterator &)>
            predicate = IteratorFilters::SubdomainEqualTo(this->my_subdomain);

          const std::vector<typename parallel::shared::Triangulation<
            dim,
            spacedim>::active_cell_iterator>
            active_halo_layer_vector =
              dealii::GridTools::compute_active_cell_halo_layer(*this,
                                                                predicate);
          std::set<typename parallel::shared::Triangulation<dim, spacedim>::
                     active_cell_iterator>
            active_halo_layer(active_halo_layer_vector.begin(),
                              active_halo_layer_vector.end());

          for (unsigned int index = 0; cell != endc; cell++, index++)
            {
              // store original/true subdomain ids:
              true_subdomain_ids_of_cells[index] = cell->subdomain_id();

              if (cell->is_locally_owned() == false &&
                  active_halo_layer.find(cell) == active_halo_layer.end())
                cell->set_subdomain_id(numbers::artificial_subdomain_id);
            }

          // loop over all cells in multigrid hierarchy and mark artificial:
          if (settings & construct_multigrid_hierarchy)
            {
              true_level_subdomain_ids_of_cells.resize(this->n_levels());

              std::function<bool(
                const typename parallel::shared::Triangulation<dim, spacedim>::
                  cell_iterator &)>
                predicate = IteratorFilters::LocallyOwnedLevelCell();
              for (unsigned int lvl = 0; lvl < this->n_levels(); ++lvl)
                {
                  true_level_subdomain_ids_of_cells[lvl].resize(
                    this->n_cells(lvl));

                  const std::vector<typename parallel::shared::Triangulation<
                    dim,
                    spacedim>::cell_iterator>
                    level_halo_layer_vector =
                      dealii::GridTools::compute_cell_halo_layer_on_level(
                        *this, predicate, lvl);
                  std::set<typename parallel::shared::
                             Triangulation<dim, spacedim>::cell_iterator>
                    level_halo_layer(level_halo_layer_vector.begin(),
                                     level_halo_layer_vector.end());

                  typename parallel::shared::Triangulation<dim, spacedim>::
                    cell_iterator cell = this->begin(lvl),
                                  endc = this->end(lvl);
                  for (unsigned int index = 0; cell != endc; cell++, index++)
                    {
                      // Store true level subdomain IDs before setting
                      // artificial
                      true_level_subdomain_ids_of_cells[lvl][index] =
                        cell->level_subdomain_id();

                      // for active cells, we must have knowledge of level
                      // subdomain ids of all neighbors to our subdomain, not
                      // just neighbors on the same level. if the cells
                      // subdomain id was not set to artitficial above, we will
                      // also keep its level subdomain id since it is either
                      // owned by this processor or in the ghost layer of the
                      // active mesh.
                      if (cell->is_active() &&
                          cell->subdomain_id() !=
                            numbers::artificial_subdomain_id)
                        continue;

                      // we must have knowledge of our parent in the hierarchy
                      if (cell->has_children())
                        {
                          bool keep_cell = false;
                          for (unsigned int c = 0;
                               c < GeometryInfo<dim>::max_children_per_cell;
                               ++c)
                            if (cell->child(c)->level_subdomain_id() ==
                                this->my_subdomain)
                              {
                                keep_cell = true;
                                break;
                              }
                          if (keep_cell)
                            continue;
                        }

                      // we must have knowledge of our neighbors on the same
                      // level
                      if (!cell->is_locally_owned_on_level() &&
                          level_halo_layer.find(cell) != level_halo_layer.end())
                        continue;

                      // mark all other cells to artificial
                      cell->set_level_subdomain_id(
                        numbers::artificial_subdomain_id);
                    }
                }
            }
        }
      else
        {
          // just store true subdomain ids
          for (unsigned int index = 0; cell != endc; cell++, index++)
            true_subdomain_ids_of_cells[index] = cell->subdomain_id();
        }

#  ifdef DEBUG
      {
        // Assert that each cell is owned by a processor
        unsigned int n_my_cells = 0;
        typename parallel::shared::Triangulation<dim,
                                                 spacedim>::active_cell_iterator
          cell = this->begin_active(),
          endc = this->end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            n_my_cells += 1;

        const unsigned int total_cells =
          Utilities::MPI::sum(n_my_cells, this->get_communicator());
        Assert(total_cells == this->n_active_cells(),
               ExcMessage("Not all cells are assigned to a processor."));
      }

      // If running with multigrid, assert that each level
      // cell is owned by a processor
      if (settings & construct_multigrid_hierarchy)
        {
          unsigned int n_my_cells = 0;
          typename parallel::shared::Triangulation<dim, spacedim>::cell_iterator
            cell = this->begin(),
            endc = this->end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned_on_level())
              n_my_cells += 1;

          const unsigned int total_cells =
            Utilities::MPI::sum(n_my_cells, this->get_communicator());
          Assert(total_cells == this->n_cells(),
                 ExcMessage("Not all cells are assigned to a processor."));
        }
#  endif
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::with_artificial_cells() const
    {
      return allow_artificial_cells;
    }



    template <int dim, int spacedim>
    const std::vector<types::subdomain_id> &
    Triangulation<dim, spacedim>::get_true_subdomain_ids_of_cells() const
    {
      return true_subdomain_ids_of_cells;
    }



    template <int dim, int spacedim>
    const std::vector<types::subdomain_id> &
    Triangulation<dim, spacedim>::get_true_level_subdomain_ids_of_cells(
      const unsigned int level) const
    {
      Assert(level < true_level_subdomain_ids_of_cells.size(),
             ExcInternalError());
      Assert(true_level_subdomain_ids_of_cells[level].size() ==
               this->n_cells(level),
             ExcInternalError());
      return true_level_subdomain_ids_of_cells[level];
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_coarsening_and_refinement()
    {
      dealii::Triangulation<dim, spacedim>::execute_coarsening_and_refinement();
      partition();
      this->update_number_cache();
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
      partition();
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
    Triangulation<dim, spacedim>::copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &other_tria)
    {
      Assert(
        (dynamic_cast<
           const dealii::parallel::distributed::Triangulation<dim, spacedim> *>(
           &other_tria) == nullptr),
        ExcMessage(
          "Cannot use this function on parallel::distributed::Triangulation."));

      dealii::parallel::TriangulationBase<dim, spacedim>::copy_triangulation(
        other_tria);
      partition();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::update_number_cache()
    {
      parallel::TriangulationBase<dim, spacedim>::update_number_cache();

      if (settings & construct_multigrid_hierarchy)
        parallel::TriangulationBase<dim, spacedim>::fill_level_ghost_owners();
    }
  } // namespace shared
} // namespace parallel

#else

namespace parallel
{
  namespace shared
  {
    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::with_artificial_cells() const
    {
      Assert(false, ExcNotImplemented());
      return true;
    }

    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      return false;
    }

    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    Triangulation<dim, spacedim>::get_true_subdomain_ids_of_cells() const
    {
      Assert(false, ExcNotImplemented());
      return true_subdomain_ids_of_cells;
    }



    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    Triangulation<dim, spacedim>::get_true_level_subdomain_ids_of_cells(
      const unsigned int) const
    {
      Assert(false, ExcNotImplemented());
      return true_level_subdomain_ids_of_cells;
    }
  } // namespace shared
} // namespace parallel


#endif


/*-------------- Explicit Instantiations -------------------------------*/
#include "shared_tria.inst"

DEAL_II_NAMESPACE_CLOSE
