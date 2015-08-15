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
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/tria.h>


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
                                                const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid):
      dealii::parallel::Triangulation<dim,spacedim>(mpi_communicator,smooth_grid,false)
    {
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
      dealii::GridTools::partition_triangulation (this->n_subdomains, *this);
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
      dealii::GridTools::partition_triangulation (this->n_subdomains, *this);
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
    {
      Assert (false, ExcNotImplemented());
    }


    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::~Triangulation ()
    {
      Assert (false, ExcNotImplemented());
    }
  }
}


#endif


/*-------------- Explicit Instantiations -------------------------------*/
#include "shared_tria.inst"

DEAL_II_NAMESPACE_CLOSE
