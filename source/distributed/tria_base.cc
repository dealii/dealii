// ---------------------------------------------------------------------
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
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/tria_base.h>


#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{

  template <int dim, int spacedim>
  Triangulation<dim,spacedim>::Triangulation (MPI_Comm mpi_communicator,
                                              const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid,
                                              const bool check_for_distorted_cells)
    :
    dealii::Triangulation<dim,spacedim>(smooth_grid,check_for_distorted_cells),
    mpi_communicator (Utilities::MPI::
                      duplicate_communicator(mpi_communicator)),
    my_subdomain (Utilities::MPI::this_mpi_process (this->mpi_communicator)),
    n_subdomains(Utilities::MPI::n_mpi_processes(mpi_communicator))
  {
#ifndef DEAL_II_WITH_MPI
    Assert(false, ExcMessage("You compiled deal.II without MPI support, for "
                             "which parallel::Triangulation is not available."));
#endif
    number_cache.n_locally_owned_active_cells.resize (n_subdomains);
  }

  template <int dim, int spacedim>
  void
  Triangulation<dim,spacedim>::copy_triangulation (const dealii::Triangulation<dim, spacedim> &old_tria)
  {
#ifndef DEAL_II_WITH_MPI
    Assert(false, ExcNotImplemented());
#endif
    if (const dealii::parallel::Triangulation<dim,spacedim> *
        old_tria_x = dynamic_cast<const dealii::parallel::Triangulation<dim,spacedim> *>(&old_tria))
      {
        mpi_communicator = Utilities::MPI::duplicate_communicator (old_tria_x->get_communicator ());
      }
  }



  template <int dim, int spacedim>
  std::size_t
  Triangulation<dim,spacedim>::memory_consumption() const
  {
    std::size_t mem=
      this->dealii::Triangulation<dim,spacedim>::memory_consumption()
      + MemoryConsumption::memory_consumption(mpi_communicator)
      + MemoryConsumption::memory_consumption(my_subdomain)
      + MemoryConsumption::memory_consumption(number_cache.n_locally_owned_active_cells)
      + MemoryConsumption::memory_consumption(number_cache.n_global_active_cells)
      + MemoryConsumption::memory_consumption(number_cache.n_global_levels);
    return mem;

  }

  template <int dim, int spacedim>
  Triangulation<dim,spacedim>::~Triangulation ()
  {
#ifdef DEAL_II_WITH_MPI
    // get rid of the unique communicator used here again
    MPI_Comm_free (&this->mpi_communicator);
#endif
  }

  template <int dim, int spacedim>
  Triangulation<dim,spacedim>::NumberCache::NumberCache()
    :
    n_global_active_cells(0),
    n_global_levels(0)
  {}

  template <int dim, int spacedim>
  unsigned int
  Triangulation<dim,spacedim>::n_locally_owned_active_cells () const
  {
    return number_cache.n_locally_owned_active_cells[my_subdomain];
  }

  template <int dim, int spacedim>
  unsigned int
  Triangulation<dim,spacedim>::n_global_levels () const
  {
    return number_cache.n_global_levels;
  }

  template <int dim, int spacedim>
  types::global_dof_index
  Triangulation<dim,spacedim>::n_global_active_cells () const
  {
    return number_cache.n_global_active_cells;
  }

  template <int dim, int spacedim>
  const std::vector<unsigned int> &
  Triangulation<dim,spacedim>::n_locally_owned_active_cells_per_processor () const
  {
    return number_cache.n_locally_owned_active_cells;
  }

  template <int dim, int spacedim>
  MPI_Comm
  Triangulation<dim,spacedim>::get_communicator () const
  {
    return mpi_communicator;
  }

#ifdef DEAL_II_WITH_MPI
  template <int dim, int spacedim>
  void
  Triangulation<dim,spacedim>::update_number_cache ()
  {
    Assert (number_cache.n_locally_owned_active_cells.size()
            ==
            Utilities::MPI::n_mpi_processes (this->mpi_communicator),
            ExcInternalError());

    std::fill (number_cache.n_locally_owned_active_cells.begin(),
               number_cache.n_locally_owned_active_cells.end(),
               0);

    number_cache.ghost_owners.clear ();
    number_cache.level_ghost_owners.clear ();

    if (this->n_levels() == 0)
      {
        // Skip communication done below if we do not have any cells
        // (meaning the Triangulation is empty on all processors). This will
        // happen when called from the destructor of Triangulation, which
        // can get called during exception handling causing a hang in this
        // function.
        number_cache.n_global_active_cells = 0;
        number_cache.n_global_levels = 0;
        return;
      }


    {
      // find ghost owners
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end();
           ++cell)
        if (cell->is_ghost())
          number_cache.ghost_owners.insert(cell->subdomain_id());

      Assert(number_cache.ghost_owners.size() < Utilities::MPI::n_mpi_processes(mpi_communicator), ExcInternalError());
    }

    if (this->n_levels() > 0)
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell = this->begin_active();
           cell != this->end(); ++cell)
        if (cell->subdomain_id() == my_subdomain)
          ++number_cache.n_locally_owned_active_cells[my_subdomain];

    unsigned int send_value
      = number_cache.n_locally_owned_active_cells[my_subdomain];
    MPI_Allgather (&send_value,
                   1,
                   MPI_UNSIGNED,
                   &number_cache.n_locally_owned_active_cells[0],
                   1,
                   MPI_UNSIGNED,
                   this->mpi_communicator);

    number_cache.n_global_active_cells
      = std::accumulate (number_cache.n_locally_owned_active_cells.begin(),
                         number_cache.n_locally_owned_active_cells.end(),
                         /* ensure sum is computed with correct data type:*/
                         static_cast<types::global_dof_index>(0));
    number_cache.n_global_levels = Utilities::MPI::max(this->n_levels(), this->mpi_communicator);
  }
#else
  template <int dim, int spacedim>
  void
  Triangulation<dim,spacedim>::update_number_cache ()
  {
    Assert (false, ExcNotImplemented());
  }

#endif

  template <int dim, int spacedim>
  types::subdomain_id
  Triangulation<dim,spacedim>::locally_owned_subdomain () const
  {
    Assert (dim > 1, ExcNotImplemented());
    return my_subdomain;
  }

  template <int dim, int spacedim>
  const std::set<unsigned int> &
  Triangulation<dim,spacedim>::
  ghost_owners () const
  {
    return number_cache.ghost_owners;
  }

  template <int dim, int spacedim>
  const std::set<unsigned int> &
  Triangulation<dim,spacedim>::
  level_ghost_owners () const
  {
    return number_cache.level_ghost_owners;
  }

} // end namespace parallel




/*-------------- Explicit Instantiations -------------------------------*/
#include "tria_base.inst"

DEAL_II_NAMESPACE_CLOSE
