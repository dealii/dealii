//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__distribute_grid_refinement_h
#define __deal2__distribute_grid_refinement_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/distributed/tria.h>

#include <vector>
#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
				     // forward declarations
    template <int dim, int spacedim> class Triangulation;


/**
 * Collection of functions controlling refinement and coarsening of
 * parallel::distributed::Triangulation objects. This namespace
 * provides similar functionality to the dealii::GridRefinement
 * namespace, except that it works for meshes that are parallel and
 * distributed.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 2009
 */
    namespace GridRefinement
    {
				       /**
					* Like
					* dealii::GridRefinement::refine_and_coarsen_fixed_number,
					* but for parallel distributed
					* triangulation.
					*
					* The vector of criteria needs to be a
					* vector of refinement criteria for
					* all cells active on the current
					* triangulation,
					* i.e. <code>tria.n_active_cells()</code>
					* (and not
					* <code>tria.n_locally_owned_active_cells()</code>). However,
					* the function will only look at the
					* indicators that correspond to those
					* cells that are actually locally
					* owned, and ignore the indicators for
					* all other cells. The function will
					* then coordinate among all processors
					* that store part of the triangulation
					* so that at the end @p
					* top_fraction_of_cells are refined,
					* where the fraction is enforced as a
					* fraction of
					* Triangulation::n_global_active_cells,
					* not
					* Triangulation::n_locally_active_cells
					* on each processor individually. In
					* other words, it may be that on some
					* processors, no cells are refined at
					* all.
					*
					* The same is true for the fraction of
					* cells that is coarsened.
					*/
      template <int dim, class Vector, int spacedim>
      void
      refine_and_coarsen_fixed_number (
	parallel::distributed::Triangulation<dim,spacedim> &tria,
	const Vector                &criteria,
	const double                top_fraction_of_cells,
	const double                bottom_fraction_of_cells);

				       /**
					* Like
					* dealii::GridRefinement::refine_and_coarsen_fixed_fraction,
					* but for parallel distributed
					* triangulation.
					*
					* The vector of criteria needs to be a
					* vector of refinement criteria for
					* all cells active on the current
					* triangulation,
					* <code>tria.n_active_cells()</code>
					* (and not
					* <code>tria.n_locally_owned_active_cells()</code>). However,
					* the function will only look at the
					* indicators that correspond to those
					* cells that are actually locally
					* owned, and ignore the indicators for
					* all other cells. The function will
					* then coordinate among all processors
					* that store part of the triangulation
					* so that at the end the smallest
					* fraction of
					* Triangulation::n_global_active_cells
					* (not
					* Triangulation::n_locally_active_cells
					* on each processor individually) is
					* refined that together make up a
					* total of @p top_fraction_of_error of
					* the total error. In other words, it
					* may be that on some processors, no
					* cells are refined at all.
					*
					* The same is true for the fraction of
					* cells that is coarsened.
					*/
      template <int dim, class Vector, int spacedim>
      void
      refine_and_coarsen_fixed_fraction (
	parallel::distributed::Triangulation<dim,spacedim> &tria,
	const Vector                &criteria,
	const double                top_fraction_of_error,
	const double                bottom_fraction_of_error);
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif //__deal2__distributed_grid_refinement_h
