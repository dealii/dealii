// ---------------------------------------------------------------------
// $Id: tria.cc 32807 2014-04-22 15:01:57Z heister $
//
// Copyright (C) 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/distributed/shared_tria.h>




#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_MPI
namespace parallel
{
  namespace shared
  {

    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::Triangulation (MPI_Comm mpi_communicator,
                                                const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid,
                                                const bool allow_artificial_cells):
      dealii::parallel::Triangulation<dim,spacedim>(mpi_communicator,smooth_grid,false),
      allow_artificial_cells(allow_artificial_cells)
    {
    }

    template <int dim, int spacedim>
    void Triangulation<dim,spacedim>::partition()
    {
      dealii::GridTools::partition_triangulation (this->n_subdomains, *this);

      true_subdomain_ids_of_cells.resize(this->n_active_cells());

      // loop over all cells and mark artificial:
      typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator
      cell = this->begin_active(),
      endc = this->end();

      if (allow_artificial_cells)
        {
          // get halo layer of (ghost) cells
          // parallel::shared::Triangulation<dim>::
          std_cxx11::function<bool (const typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator &)> predicate
            = IteratorFilters::SubdomainEqualTo(this->my_subdomain);

          const std::vector<typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator>
          active_halo_layer_vector = GridTools::compute_active_cell_halo_layer (*this, predicate);
          std::set<typename parallel::shared::Triangulation<dim,spacedim>::active_cell_iterator>
          active_halo_layer(active_halo_layer_vector.begin(), active_halo_layer_vector.end());

          for (unsigned int index=0; cell != endc; cell++, index++)
            {
              // store original/true subdomain ids:
              true_subdomain_ids_of_cells[index] = cell->subdomain_id();

              if (cell->is_locally_owned() == false &&
                  active_halo_layer.find(cell) == active_halo_layer.end())
                cell->set_subdomain_id(numbers::artificial_subdomain_id);
            }
        }
      else
        {
          // just store true subdomain ids
          for (unsigned int index=0; cell != endc; cell++, index++)
            true_subdomain_ids_of_cells[index] = cell->subdomain_id();

        }
    }

    template <int dim, int spacedim>
    bool
    Triangulation<dim,spacedim>::with_artificial_cells() const
    {
      return allow_artificial_cells;
    }

    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    Triangulation<dim,spacedim>::get_true_subdomain_ids_of_cells() const
    {
      return true_subdomain_ids_of_cells;
    }

    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::~Triangulation ()
    {

    }

    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::execute_coarsening_and_refinement ()
    {
      dealii::Triangulation<dim,spacedim>::execute_coarsening_and_refinement ();
      partition();
      this->update_number_cache ();
    }

    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::create_triangulation (const std::vector< Point< spacedim > > &vertices,
                                                       const std::vector< CellData< dim > > &cells,
                                                       const SubCellData &subcelldata)
    {
      try
        {
          dealii::Triangulation<dim,spacedim>::
          create_triangulation (vertices, cells, subcelldata);
        }
      catch (const typename dealii::Triangulation<dim,spacedim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }
      partition();
      this->update_number_cache ();
    }

  }
}

#else

namespace parallel
{
  namespace shared
  {
    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::Triangulation ()
      :
      dealii::parallel::Triangulation<dim,spacedim>(MPI_COMM_SELF)
    {
      Assert (false, ExcNotImplemented());
    }

    template <int dim, int spacedim>
    bool
    Triangulation<dim,spacedim>::with_artificial_cells() const
    {
      Assert (false, ExcNotImplemented());
      return true;
    }

    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    Triangulation<dim,spacedim>::get_true_subdomain_ids_of_cells() const
    {
      Assert (false, ExcNotImplemented());
      return true_subdomain_ids_of_cells;
    }

  }
}


#endif


/*-------------- Explicit Instantiations -------------------------------*/
#include "shared_tria.inst"

DEAL_II_NAMESPACE_CLOSE
