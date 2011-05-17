//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__distributed_solution_transfer_h
#define __deal2__distributed_solution_transfer_h

#include <deal.II/base/config.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{

  namespace distributed
  {
/**
 * Transfers a discrete FE function (like a solution vector) by
 * interpolation while refining and/or
 * coarsening a distributed grid and
 * handles the necessary communication.
 *
 * <h3>Usage</h3>
 * @verbatim
 * SolutionTransfer<dim, Vector<double> > soltrans(dof_handler);
 *                                     // flag some cells for refinement
 *                                     // and coarsening, e.g.
 * GridRefinement::refine_and_coarsen_fixed_fraction(
 *   *tria, error_indicators, 0.3, 0.05);
 *                                     // prepare the triangulation,
 * tria->prepare_coarsening_and_refinement();
 *                                     // prepare the SolutionTransfer object
 *                                     // for coarsening and refinement and give
 *                                     // the solution vector that we intend to
 *                                     // interpolate later,
 * soltrans.prepare_for_coarsening_and_refinement(solution);
 *                                     // actually execute the refinement,
 * tria->execute_coarsening_and_refinement ();
 *                                     // redistribute dofs,
 * dof_handler->distribute_dofs (fe);
 *                                     // and interpolate the solution
 * Vector<double> interpolated_solution(dof_handler->n_dofs());
 * soltrans.interpolate(interpolated_solution);
 * @endverbatim
 * @ingroup distributed
 * @author Timo Heister, 2009
 */
    template<int dim, typename VECTOR, class DH=DoFHandler<dim> >
    class SolutionTransfer
    {
      public:
					 /**
					  * Constructor, takes the current DoFHandler
					  * as argument.
					  */
	SolutionTransfer(const DH &dof);
					 /**
					  * Destructor.
					  */
	~SolutionTransfer();

					 /**
					  * Prepares the @p SolutionTransfer for
					  * coarsening and refinement. It
					  * stores the dof indices of each cell and
					  * stores the dof values of the vectors in
					  * @p all_in in each cell that'll be coarsened.
					  * @p all_in includes all vectors
					  * that are to be interpolated
					  * onto the new (refined and/or
					  * coarsenend) grid.
					  */
	void prepare_for_coarsening_and_refinement (const std::vector<const VECTOR*> &all_in);

					 /**
					  * Same as previous function
					  * but for only one discrete function
					  * to be interpolated.
					  */
	void prepare_for_coarsening_and_refinement (const VECTOR &in);

					 /**
					  *
					  */
	void interpolate (std::vector<VECTOR*> &all_out);

					 /**
					  * Same as the previous function.
					  * It interpolates only one function.
					  * It assumes the vectors having the
					  * right sizes (i.e. <tt>in.size()==n_dofs_old</tt>,
					  * <tt>out.size()==n_dofs_refined</tt>)
					  *
					  * Multiple calling of this function is
					  * NOT allowed. Interpolating
					  * several functions can be performed
					  * in one step by using
					  * <tt>interpolate (all_in, all_out)</tt>
					  */
	void interpolate (VECTOR &out);


					 /**
					  * return the size in bytes that need
					  * to be stored per cell.
					  */
	unsigned int get_data_size() const;



      private:
					 /**
					  * Pointer to the degree of
					  * freedom handler to work
					  * with.
					  */
	SmartPointer<const DH,SolutionTransfer<dim,VECTOR,DH> > dof_handler;

					 /**
					  * A vector that stores
					  * pointers to all the
					  * vectors we are supposed to
					  * copy over from the old to
					  * the new mesh.
					  */
	std::vector<const VECTOR*> input_vectors;

					 /**
					  * The offset that the
					  * Triangulation has assigned
					  * to this object starting at
					  * which we are allowed to
					  * write.
					  */
	unsigned int offset;

					 /**
					  * A callback function used
					  * to pack the data on the
					  * current mesh into objects
					  * that can later be
					  * retrieved after
					  * refinement, coarsening and
					  * repartitioning.
					  */
	void pack_callback(const typename Triangulation<dim,dim>::cell_iterator &cell,
			   const typename Triangulation<dim,dim>::CellStatus status,
			   void* data);

					 /**
					  * A callback function used
					  * to unpack the data on the
					  * current mesh that has been
					  * packed up previously on
					  * the mesh before
					  * refinement, coarsening and
					  * repartitioning.
					  */
	void unpack_callback(const typename Triangulation<dim,dim>::cell_iterator &cell,
			     const typename Triangulation<dim,dim>::CellStatus status,
			     const void* data,
			     std::vector<VECTOR*> &all_out);

    };


  }
}



DEAL_II_NAMESPACE_CLOSE

#endif
