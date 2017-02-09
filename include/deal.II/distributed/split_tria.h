// ---------------------------------------------------------------------
// $Id: tria.h 32739 2014-04-08 16:39:47Z denis.davydov $
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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

#ifndef dealii__distributed__split_tria_h
#define dealii__distributed__split_tria_h


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

namespace parallel
{

#ifdef DEAL_II_WITH_MPI


  namespace split
  {

    /**
     * Split triangulation for parallel computations where the coarse mesh
     * is distributed. Work in progress.
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
       *
       * If @p allow_aritifical_cells is true, this class will behave similar
       * to parallel::distributed::Triangulation in that there will be locally
       * owned, ghost and artificial cells.
       *
       * Otherwise all non-locally owned cells are considered ghost.
       */
      Triangulation (MPI_Comm mpi_communicator,
                     const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing =
                       (dealii::Triangulation<dim,spacedim>::none));

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

      /**
       * Coarsen and refine the mesh according to refinement and coarsening
       * flags set.
       *
       * This step is equivalent to the dealii::Triangulation class with an
       * addition of calling dealii::GridTools::partition_triangulation() at
       * the end.
       */
      virtual void execute_coarsening_and_refinement ();

//      /**
//       * Create a triangulation.
//       *
//       * This function also partitions triangulation based on the MPI
//       * communicator provided to the constructor.
//       */
//      virtual void create_triangulation (const std::vector< Point< spacedim > > &vertices,
//                                         const std::vector< CellData< dim > > &cells,
//                                         const SubCellData &subcelldata);

      /**
       * Copy @p other_tria to this triangulation.
       *
       * This function also partitions triangulation based on the MPI
       * communicator provided to the constructor.
       *
       * @note This function can not be used with parallel::distributed::Triangulation,
       * since it only stores those cells that it owns, one layer of ghost cells around
       * the ones it locally owns, and a number of artificial cells.
       */
      virtual void copy_triangulation (const dealii::Triangulation<dim, spacedim> &other_tria);

      /**
       * Read the data of this object from a stream for the purpose of
       * serialization. Throw away the previous content.
       *
       * This function first does the same work as in dealii::Triangulation::load,
       * then partitions the triangulation based on the MPI communicator
       * provided to the constructor.
       */
      template <class Archive>
      void load (Archive &ar, const unsigned int version);



    private:
      /**
       * This function calls GridTools::partition_triangulation () and if
       * requested in the constructor of the class marks artificial cells.
       */
      void partition();

      /**
       * A vector containing subdomain IDs of cells obtained by partitioning
       * using METIS. In case allow_artificial_cells is false, this vector is
       * consistent with IDs stored in cell->subdomain_id() of the
       * triangulation class. When allow_artificial_cells is true, cells which
       * are artificial will have cell->subdomain_id() == numbers::artificial;
       *
       * The original parition information is stored to allow using sequential
       * DoF distribution and partitioning functions with semi-artificial
       * cells.
       */
      std::vector<types::subdomain_id> true_subdomain_ids_of_cells;
    };

    template <int dim, int spacedim>
    template <class Archive>
    void
    Triangulation<dim,spacedim>::load (Archive &ar,
                                       const unsigned int version)
    {
      dealii::Triangulation<dim,spacedim>::load (ar, version);
      partition();
      this->update_number_cache ();
    }
  }
#else

  TODO
  namespace shared
  {

    /**
     * Dummy class the compiler chooses for parallel shared triangulations if
     * we didn't actually configure deal.II with the MPI library. The
     * existence of this class allows us to refer to
     * parallel::shared::Triangulation objects throughout the library even if
     * it is disabled.
     *
     * Since the constructor of this class is private, no such objects can
     * actually be created if MPI is not available.
     */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::parallel::Triangulation<dim,spacedim>
    {
    public:

      /**
       * A dummy function to return empty vector.
       */
      const std::vector<types::subdomain_id> &get_true_subdomain_ids_of_cells() const;

      /**
       * A dummy function which always returns true.
       */
      bool with_artificial_cells() const;
    private:
      /**
       * Constructor.
       */
      Triangulation ();

      /**
       * A dummy vector.
       */
      std::vector<types::subdomain_id> true_subdomain_ids_of_cells;
    };
  }


#endif
}

DEAL_II_NAMESPACE_CLOSE

#endif
