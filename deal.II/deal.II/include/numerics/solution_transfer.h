//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__solution_transfer_h
#define __deal2__solution_transfer_h


/*----------------------------   solutiontransfer.h     ----------------------*/


#include <base/config.h>
#include <base/smartpointer.h>
#include <lac/vector.h>
#include <base/exceptions.h>
#include <dofs/dof_handler.h>

#include <vector>


/**
 * Transfers a discrete FE function (like a solution vector) by interpolation
 * while refining and/or coarsening a grid. During interpolation the
 * vector is reinitialized to the new size and filled with the interpolated
 * values.
 * 
 * <h3>Usage</h3>
 *
 * As the interpolation while
 * coarsening is much more complicated to organize 
 * (see further documentation below) than interpolation while pure refinement,
 * @p SolutionTransfer offers two possible usages.
 * <ul>
 * <li> If the grid will only be purely refined
 * (i.e. not locally coarsened) then use @p SolutionTransfer as follows
 * @verbatim
 * SolutionTransfer<dim, double> soltrans(*dof_handler);
 *                                     // flag some cells for refinement, e.g.
 * GridRefinement::refine_and_coarsen_fixed_fraction(
 *   *tria, error_indicators, 0.3, 0);
 *                                     // prepare the triangulation
 *                                     // for refinement,
 * tria->prepare_coarsening_and_refinement();
 *                                     // tell the SolutionTransfer object
 *                                     // that we intent to do pure refinement,
 * soltrans.prepare_for_pure_refinement();
 *                                     // actually execute the refinement,
 * tria->execute_coarsening_and_refinement();
 *                                     // and redistribute dofs.
 * dof_handler->distribute_dofs (fe);
 * @endverbatim
 *
 * Then there are three different possibilities of how to proceed. Either
 * @verbatim
 *                                     // resize and interpolate a solution
 *                                     // vector `in-place'
 * soltrans.refine_interpolate(solution);
 * @endverbatim
 * or, when the old solution vector is still needed,
 * @verbatim
 *                                     // take a copy of the solution vector
 * Vector<double> solution_old(solution);
 *                                     // resize solution vector to the correct
 *                                     // size, as the @p refine_interpolate
 *                                     // function requires the vectors to be
 *                                     // of right sizes
 * solution.reinit(dof_handler->n_dofs());
 *                                     // and finally interpolate
 * soltrans.refine_interpolate(solution_old, solution);
 * @endverbatim
 *
 * Although the @p refine_interpolate functions are allowed to be
 * called multiple times, e.g. for interpolating several solution
 * vectors, there is following possibility of interpolating several
 * functions simultaneously.
 * @verbatim
 * vector<Vector<double> > solutions_old(n_vectors, Vector<double> (n));
 * ...
 * vector<Vector<double> > solutions(n_vectors, Vector<double> (n));
 * soltrans.refine_interpolate(solutions_old, solutions);
 * @endverbatim
 * <li> If the grid will be refined AND coarsened
 * then use @p SolutionTransfer as follows
 * @verbatim
 * SolutionTransfer<dim, double> soltrans(*dof_handler);
 *                                     // flag some cells for refinement
 *                                     // and coarsening, e.g.
 * GridRefinement::refine_and_coarsen_fixed_fraction(
 *   *tria, error_indicators, 0.3, 0.05);
 *                                     // prepare the triangulation,
 * tria->prepare_coarsening_and_refinement();
 *                                     // prepare the SolutionTransfer object
 *                                     // for coarsening and refinement and give
 *                                     // the solution vector that we intent to
 *                                     // interpolate later,
 * soltrans.prepare_for_coarsening_and_refinement(solution);
 *                                     // actually execute the refinement,
 * tria->execute_coarsening_and_refinement ();
 *                                     // redistribute dofs,
 * dof_handler->distribute_dofs (fe);
 *                                     // and interpolate the solution
 * Vector<double> interpolate_solution(dof_handler->n_dofs());
 * soltrans.interpolate(solution, interpolated_solution);
 * @endverbatim
 *
 * Multiple calls to the function 
 * <tt>interpolate (const Vector<number> &in, Vector<number> &out)</tt>
 * are NOT allowed. Interpolating several functions can be performed in one step
 * by using 
 * <tt>void interpolate (const vector<Vector<number> >&all_in, vector<Vector<number> >&all_out) const</tt>,
 * and using the respective @p prepare_for_coarsening_and_refinement function 
 * taking several vectors as input before actually refining and coarsening the
 * triangulation (see there).
 * </ul>
 *
 * For deleting all stored data in @p SolutionTransfer and reinitializing it
 * use the <tt>clear()</tt> function.
 *
 * Note that the @p user_pointer of some cells are used. 
 * Be sure that you don't need them otherwise.
 *
 * The template argument @p number denotes the data type of the vectors you want
 * to transfer.
 *
 *
 * <h3>Implementation</h3>
 *
 * <ul>
 * <li> Solution transfer while pure refinement. Assume that we have got a
 * solution vector on the current (original) grid.
 * Each entry of this vector belongs to one of the
 * DoFs of the discretisation. If we now refine the grid then the calling of
 * DoFHandler::distribute_dofs() will change at least some of the
 * DoF indices. Hence we need to store the DoF indices of all active cells
 * before the refinement. The @p user_pointer of each active cell
 * is used to point to the vector of these DoF indices of that cell, all other
 * @p user_pointers are cleared. All this is
 * done by prepare_for_pure_refinement().
 *
 * In the function <tt>refine_interpolate(in,out)</tt> and on each cell where the
 * @p user_pointer is set (i.e. the cells that were active in the original grid)
 * we can now access the local values of the solution vector @p in
 * on that cell by using the stored DoF indices. These local values are
 * interpolated and set into the vector @p out that is at the end the
 * discrete function @p in interpolated on the refined mesh.
 *
 * The <tt>refine_interpolate(in,out)</tt> function can be called multiply for
 * arbitrary many discrete functions (solution vectors) on the original grid. 
 *
 * <li> Solution transfer while coarsening and refinement. After 
 * calling Triangulation@p ::prepare_coarsening_and_refinement the
 * coarsen flags of either all or none of the children of a 
 * (father-)cell are set. While coarsening 
 * (Triangulation@p ::execute_coarsening_and_refinement)
 * the cells that are not needed any more will be deleted from the Triangulation.
 * 
 * For the interpolation from the (to be coarsenend) children to their father
 * the children cells are needed. Hence this interpolation
 * and the storing of the interpolated values of each of the discrete functions
 * that we want to interpolate needs to take place before these children cells
 * are coarsened (and deleted!!). Again the @p user_pointers of the cells are
 * set to point to these values (see below). 
 * Additionally the DoF indices of the cells
 * that will not be coarsened need to be stored according to the solution
 * transfer while pure refinement (cf there). All this is performed by
 * <tt>prepare_for_coarsening_and_refinement(all_in)</tt> where the 
 * <tt>vector<Vector<number> >vector all_in</tt> includes
 * all discrete functions to be interpolated onto the new grid.
 *
 * As we need two different kinds of pointers (<tt>vector<unsigned int> *</tt> for the Dof
 * indices and <tt>vector<Vector<number> > *</tt> for the interpolated DoF values)
 * we use the @p Pointerstruct that includes both of these pointers and
 * the @p user_pointer of each cell points to these @p Pointerstructs. 
 * On each cell only one of the two different pointers is used at one time 
 * hence we could use the
 * <tt>void * user_pointer</tt> as <tt>vector<unsigned int> *</tt> at one time and as 
 * <tt>vector<Vector<number> > *</tt> at the other but using this @p Pointerstruct
 * in between makes the use of these pointers more safe and gives better
 * possibility to expand their usage.
 * 
 * In <tt>interpolate(all_in, all_out)</tt> the refined cells are treated according
 * to the solution transfer while pure refinement. Additionally, on each
 * cell that is coarsened (hence previously was a father cell), 
 * the values of the discrete
 * functions in @p all_out are set to the stored local interpolated values 
 * that are accessible due to the 'vector<Vector<number> > *' pointer in 
 * @p Pointerstruct that is pointed to by the @p user_pointer of that cell.
 * It is clear that <tt>interpolate(all_in, all_out)</tt> only can be called with
 * the <tt>vector<Vector<number> > all_in</tt> that previously was the parameter
 * of the <tt>prepare_for_coarsening_and_refinement(all_in)</tt> function. Hence 
 * <tt>interpolate(all_in, all_out)</tt> can (in contrast to 
 * <tt>refine_interpolate(in, out)</tt>) only be called once.
 * </ul>
 *
 * @author Ralf Hartmann, 1999
 */
template<int dim, typename number>
class SolutionTransfer
{
  public:
    
				     /**
				      * Constructor, takes the current DoFHandler
				      * as argument.
				      */
    SolutionTransfer(const DoFHandler<dim> &dof);

    				     /**
				      * Destructor
				      */
    ~SolutionTransfer();
    
				     /**
				      * Reinit this class to the state that
				      * it has
				      * directly after calling the Constructor
				      */
    void clear();

				     /**
				      * Prepares the @p SolutionTransfer for
				      * pure refinement. It
				      * stores the dof indices of each cell.
				      * After calling this function 
				      * only calling the @p refine_interpolate
				      * functions is allowed.
				      */
    void prepare_for_pure_refinement();

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
    void prepare_for_coarsening_and_refinement (const std::vector<Vector<number> > &all_in);
    
				     /**
				      * Same as previous function
				      * but for only one discrete function
				      * to be interpolated.
				      */
    void prepare_for_coarsening_and_refinement (const Vector<number> &in);
		      
				     /**
				      * This function
				      * interpolates the discrete function @p in,
				      * which is a vector on the grid before the
				      * refinement, to the function @p out
				      * which then is a vector on the refined grid.
				      * It assumes the vectors having the
				      * right sizes (i.e. <tt>in.size()==n_dofs_old</tt>,
				      * <tt>out.size()==n_dofs_refined</tt>)
       				      *
				      * Calling this function is allowed only
				      * if @p prepare_for_pure_refinement is called
				      * and the refinement is
				      * executed before.
				      * Multiple calling of this function is
				      * allowed. e.g. for interpolating several
				      * functions.
				      */
    void refine_interpolate (const Vector<number> &in,
			     Vector<number> &out) const;
    
				     /**
				      * Same as <tt>interpolate(in,out)</tt>
				      * but it interpolates
				      * just 'in-place'. Therefore @p vec will be
				      * reinitialized to the new needed vector
				      * dimension.
				      */
    void refine_interpolate (Vector<number> &vec) const;
      
				     /**
				      * This function
				      * interpolates the discrete functions
				      * that are stored in @p all_out onto
				      * the refined and/or coarsenend grid.
				      * It assumes the vectors in @p all_in
				      * denote the same vectors
				      * as in @p all_in as parameter of
				      * <tt>prepare_for_refinement_and_coarsening(all_in)</tt>.
				      * However, there is no way of verifying
				      * this internally, so be careful here.
				      *
				      * Calling this function is
				      * allowed only if first
				      * Triangulation@p ::prepare_coarsening_and_refinement,
				      * second
				      * @p SolutionTransfer::prepare_for_coarsening_and_refinement,
				      * an then third
				      * Triangulation@p ::execute_coarsening_and_refinement
				      * are called before. Multiple
				      * calling of this function is
				      * NOT allowed. Interpolating
				      * several functions can be
				      * performed in one step.
				      *
				      * The number of output vectors
				      * is assumed to be the same as
				      * the number of input
				      * vectors. Also, the sizes of
				      * the output vectors are assumed
				      * to be of the right size
				      * (@p n_dofs_refined). Otherwise
				      * an assertion will be thrown.
				      */
    void interpolate (const std::vector<Vector<number> >&all_in,
		      std::vector<Vector<number> >      &all_out) const;
      
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
    void interpolate (const Vector<number> &in,
		      Vector<number>       &out) const;

    				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception
				      */
    DeclException0(ExcNotPrepared);

				     /**
				      * Exception
				      */
    DeclException0(ExcAlreadyPrepForRef);

				     /**
				      * Exception
				      */
    DeclException0(ExcAlreadyPrepForCoarseAndRef);
				     
				     /**
				      * Exception
				      */
    DeclException0(ExcTriaPrepCoarseningNotCalledBefore);
    
				     /**
				      * Exception
				      */
    DeclException0(ExcNoInVectorsGiven);

				     /**
				      * Exception
				      */
    DeclException0(ExcVectorsDifferFromInVectors);

				     /**
				      * Exception
				      */
    DeclException2(ExcWrongVectorSize,
		   int, int,
		   << "The size of the vector is " << arg1
		   << " although it should be " << arg2 << ".");
				  
  private:

				     /**
				      * Pointer to the degree of freedom handler
				      * to work with.
				      */
    SmartPointer<const DoFHandler<dim> > dof_handler;
    
				     /**
				      * Stores the number of DoFs before the
				      * refinement and/or coarsening.
				      */
    unsigned int n_dofs_old;

				     /**
				      * Declaration of
				      * @p PreparationState that
				      * denotes the three possible
				      * states of the
				      * @p SolutionTransfer: being
				      * prepared for 'pure
				      * refinement', prepared for
				      * 'coarsening and refinement' or
				      * not prepared.
				      */
    enum PreparationState {
	  none, pure_refinement, coarsening_and_refinement
    };

				     /**
				      * Definition of the respective variable.
				      */
    PreparationState prepared_for;


				     /**
				      * Is used for @p prepare_for_refining
				      * (of course also for
				      * @p repare_for_refining_and_coarsening)
				      * and stores all dof indices
				      * of the cells that'll be refined
				      */
    std::vector<std::vector<unsigned int> > indices_on_cell;

				     /**
				      * All cell data (the dof indices and
				      * the dof values)
				      * should be accessable from each cell.
				      * As each cell has got only one
				      * @p user_pointer, multiple pointers to the
				      * data need to be packetized in a structure.
				      * Note that in our case on each cell
				      * either the
				      * <tt>vector<unsigned int> indices</tt> (if the cell
				      * will be refined) or the
				      * <tt>vector<double> dof_values</tt> (if the
				      * children of this cell will be deleted)
				      * is needed, hence one @p user_pointer should
				      * be sufficient, but to allow some errorchecks
				      * and to preserve the user from making
				      * user errors the @p user_pointer will be
				      * 'multiplied' by this structure.
				      */
    struct Pointerstruct {
	unsigned int memory_consumption () const;
	
	std::vector<unsigned int>    *indices_ptr;
	std::vector<Vector<number> > *dof_values_ptr;
    };

				     /**
				      * Vector of all @p Pointerstructs (cf. there).
				      * It makes it
				      * easier to delete all these structs
				      * (without going over all <tt>cell->user_pointer</tt>)
				      * after they are not used any more, and
				      * collecting all these structures in a vector
				      * helps avoiding fraqmentation of the memory.
				      */
    std::vector<Pointerstruct> all_pointerstructs;

				     /**
				      * Is used for
				      * @p prepare_for_refining_and_coarsening
				      * The interpolated dof values
				      * of all cells that'll be coarsened
				      * will be stored in this vector.
				      */
    std::vector<std::vector<Vector<number> > > dof_values_on_cell;
};


/*----------------------------   solutiontransfer.h     ---------------------------*/

#endif
/*----------------------------   solutiontransfer.h     ---------------------------*/
