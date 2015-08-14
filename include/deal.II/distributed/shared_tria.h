// ---------------------------------------------------------------------
// $Id: tria.h 32739 2014-04-08 16:39:47Z denis.davydov $
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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

#ifndef __deal2__distributed__shared_tria_h
#define __deal2__distributed__shared_tria_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/base/std_cxx1x/function.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <set>
#include <vector>
#include <list>
#include <utility>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif


DEAL_II_NAMESPACE_OPEN

template <int, int> class Triangulation;


namespace parallel
{

#ifdef DEAL_II_WITH_MPI


  namespace shared
  {

    /**
     * This is an extension of dealii::Triangulation class to automatically
     * partition triangulation when run with MPI.
     * Different from the parallel::distributed::Triangulation, the entire mesh
     * is stored on each processor. However, cells are labeled according to
     * the id of the processor which "owns" them. The partitioning is done
     * automatically inside the DoFHandler by calling Metis.
     * This enables distributing DoFs among processors and therefore splitting
     * matrices and vectors across processors.
     * The usage of this class is demonstrated in Step-18.
     *
     * @author Denis Davydov, 2015
     * @ingroup distributed
     *
     */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::parallel::Triangulation<dim,spacedim>
    {
    public:
      typedef typename dealii::Triangulation<dim,spacedim>::active_cell_iterator active_cell_iterator;
      typedef typename dealii::Triangulation<dim,spacedim>::cell_iterator        cell_iterator;

      /**
       * Constructor.
       */
      Triangulation (MPI_Comm mpi_communicator,
                     const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing =
                       (dealii::Triangulation<dim,spacedim>::none) );

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

      /**
       * Coarsen and refine the mesh according to refinement and
       * coarsening flags set.
       *
       * This step is equivalent to the dealii::Triangulation class
       * with an addition of calling dealii::GridTools::partition_triangulation() at the end.
       */
      virtual void execute_coarsening_and_refinement ();

      /**
        * Create a triangulation.
        *
        * This function also partitions triangulation based on the
        * MPI communicator provided to constructor.
        */
      virtual void create_triangulation (const std::vector< Point< spacedim > > &vertices,
                                         const std::vector< CellData< dim > > &cells,
                                         const SubCellData &subcelldata);

    };
  }
#else

  namespace shared
  {

    /**
     * Dummy class the compiler chooses for parallel shared
     * triangulations if we didn't actually configure deal.II with the
     * MPI library. The existence of this class allows us to refer
     * to parallel::shared::Triangulation objects throughout the
     * library even if it is disabled.
     *
     * Since the constructor of this class is private, no such objects
     * can actually be created if we don't have p4est available.
     */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::parallel::Triangulation<dim,spacedim>
    {
    private:
      /**
       * Constructor.
       */
      Triangulation ();
    public:

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

    };
  }


#endif
}

DEAL_II_NAMESPACE_CLOSE

#endif
