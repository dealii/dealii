/*----------------------------   solutiontransfer.h     ----------------------*/
/*      $Id$      */
/*           Ralf Hartmann, University of Heidelberg                          */
#ifndef __solutiontransfer_H
#define __solutiontransfer_H
/*----------------------------   solutiontransfer.h     ----------------------*/

#include <base/exceptions.h>
#include <vector.h>
template <typename number> class Vector;
template <int dim> class DoFHandler;

/**
 * Transfers a discrete FE function (like a solution vector) by interpolation
 * while refining and/or coarsening a grid. During interpolation the
 * vector is reinitialized to the new size and filled with the interpolated
 * values.
 * 
 * \subsection{Usage}
 *
 * As the interpolation while
 * coarsening is much more complicated to organize 
 * (see further documentation below) than interpolation while pure refinement,
 * #SolutionTransfer# offers two possible usages.
 * \begin{itemize}
 * \item If the grid will only be purely refined
 * (i.e. not locally coarsened) then use #SolutionTransfer# as follows
 * \begin{verbatim}
 * SolutionTransfer<dim, double> soltrans(*dof_handler);
 * soltrans.prepare_for_pure_refinement();
 *                                     // some refinement e.g.
 * tria->refine_and_coarsen_fixed_fraction(error_indicator, 0.3, 0);
 * tria->execute_coarsening_and_refinement();
 * dof_handler->distribute_dofs (fe);
 * soltrans.refine_interpolate(solution);
 *                                     // if necessary interpolate some
 *                                     // more functions
 * soltrans.refine_interpolate(sol2);
 * ...
 * \end{verbatim}
 * \item If the grid will be coarsenend and refined
 * then use #SolutionTransfer# as follows
 * \begin{verbatim}
 * SolutionTransfer<dim, double> soltrans(*dof_handler);
 *                                     // some refinement e.g.
 * tria->refine_and_coarsen_fixed_fraction(error_indicator, 0.3, 0);
 *                                     // very important:
 * tria->prepare_coarsening_and_refinement();
 * soltrans.prepare_for_coarsening_and_refinement(solution);
 * tria->execute_coarsening_and_refinement ();
 * dof_handler->distribute_dofs (fe);
 * soltrans.interpolate(solution);
 * \end{verbatim}
 
 * Multiple calling of this function #void interpolate (Vector<number> &out) const#
 * is NOT allowed. Interpolating several functions can be performed in one step
 * by using #void interpolate (vector<Vector<number> >&all_out) const#, and
 * using the respective #prepare_for_coarsening_and_refinement# function taking
 * several vectors as input before actually refining and coarsening the
 * triangulation (see there).
 * \end{itemize}
 *
 * For deleting all stored data in and reinitializing the
 * #SolutionTransfer# use the #clear()# function.
 *
 * Note that the #user_pointer# of some cells are used. Be sure that you don't need
 * them otherwise.
 *
 * The template argument #number# denotes the data type of the vectors you want
 * to transfer.
 *
 * @author Ralf Hartmann, 1999 */
template<int dim, typename number>
class SolutionTransfer
{
  public:	
				     /**
				      * Constructor, takes the current #DoFHandler#
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
				      * Prepares the #SolutionTransfer# for
				      * pure refinement. It
				      * stores the dof indices of each cell.
				      * After calling this function 
				      * only calling the #refine_interpolate#
				      * functions is allowed.
				      */
    void prepare_for_pure_refinement();

   				     /**
				      * Prepares the #SolutionTransfer# for
				      * coarsening and refinement. It
				      * stores the dof indices of each cell and
				      * stores the dof values of the vectors in
				      * #all_in# in each cell that'll be coarsened.
				      * #all_in# includes all vectors
				      * that are to be interpolated
				      * onto the new (refined and/or
				      * coarsenend) grid.
				      */
    void prepare_for_coarsening_and_refinement (const vector<Vector<number> > &all_in);
    
				     /**
				      * Same as previous function
				      * but for only one discrete function
				      * to interpolate.
				      */
    void prepare_for_coarsening_and_refinement (const Vector<number> &in);
		      
				     /**
				      * This function
				      * interpolates the discrete function #in#,
				      * which is a vector on the grid before the
				      * refinement, to the function #out#
				      * which then is a vector on the refined grid.
				      * It assumes the vectors having the
				      * right sizes (i.e. in.size()==n_dofs_old,
				      * out.size()==n_dofs_refined)
       				      *
				      * Calling this function is allowed only
				      * if #prepare_for_pure_refinement# is called
				      * and the refinement is
				      * executed before.
				      * Multiple calling of this function is
				      * allowed. e.g. for interpolating several
				      * functions.
				      */
    void refine_interpolate (const Vector<number> &in,
			     Vector<number> &out) const;
    
				     /**
				      * Same as #interpolate(Vector<number> in,
				      * Vector<number> out)# but it interpolates
				      * just 'in-place'. Therefore #vec# will be
				      * reinitialized to the new needed vector
				      * dimension.
				      */
    void refine_interpolate (Vector<number> &vec) const;
      
				     /**
				      * This function
				      * interpolates the discrete functions
				      * that are stored in #all_out# onto
				      * the refined and/or coarsenend grid.
				      * It assumes the vectors in #all_in#
				      * denote the same vectors
				      * as in #all_in# as parameter
				      * of #prepare_for_refinement_and_coarsening
				      * (vector<Vector<number> > &all_in)#.
				      * However, there is no way of verifying
				      * this internally, so be careful here.
				      *
				      * Calling this function is allowed only
				      * if first #Triangulation::prepare_coarsening_
				      * and_refinement#, second
				      * #SolutionTransfer::prepare_for_coarsening_
				      * and_refinement#, an then third 
				      * #Triangulation::execute_coarsening_
				      * and_refinement# are called before.
				      * Multiple calling of this function is
				      * NOT allowed. Interpolating
				      * several functions can be performed
				      * in one step.
				      */
    void interpolate (const vector<Vector<number> >&all_in,
		      vector<Vector<number> >      &all_out) const;
      
				     /**
				      * Same as the previous function.
				      * It interpolates only one function.
				      * 
				      * Multiple calling of this function is
				      * NOT allowed. Interpolating
				      * several functions can be performed
				      * in one step by using #void 
				      * interpolate (vector<Vector<number>
				      * >&all_out) const#
				      */
    void interpolate (const Vector<number> &in,
		      Vector<number>       &out) const;

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
    DeclException0(ExcAlreadyPrepForRefAndCoarse);
				     
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
		   << "although it should be " << arg2 << ".");
    
				     /**
				      * Exception
				      */
    DeclException0(ExcInternalError);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVectorSize,
		    int, int,
		    << "The data vector has the size " << arg1
		    << ", but " << arg2 << " was expected.");
				  
  private:

				     /**
				      * Pointer to the degree of freedom handler
				      * to work with.
				      */
    const DoFHandler<dim> *dof_handler;
    
				     /**
				      * Stores the number of DoFs before the
				      * refinement and/or coarsening.
				      */
    unsigned int n_dofs_old;

				     /**
				      * Denotes whether the #SolutionTransfer#
				      * is 'prepared for pure refinement'
				      * or not.
				      */
    bool prepared_for_pure_refinement;

				     /**
				      * Denotes whether the #SolutionTransfer#
				      * is 'prepared for coarsening and refinement'
				      * or not.
				      */
    bool prepared_for_coarsening_and_refinement;

				     /**
				      * Is used for #prepare_for_refining#
				      * (of course also for
				      * #repare_for_refining_and_coarsening#)
				      * and stores all dof indices
				      * of the cells that'll be refined
				      */
    vector<vector<int> > indices_on_cell;

				     /**
				      * All cell data (the dof indices and
				      * the dof values)
				      * should be accessable from each cell.
				      * As each cell has got only one
				      * #user_pointer#, multiple pointers to the
				      * data need to be packetized in a structure.
				      * Note that in our case on each cell
				      * either the
				      * #vector<int> indices# (if the cell
				      * will be refined) or the
				      * #vector<double> dof_values# (if the
				      * children of this cell will be deleted)
				      * is needed, hence one user_pointer should
				      * be sufficient, but to allow some errorchecks
				      * and to preserve the user from making
				      * user errors the #user_pointer# will be
				      * 'multiplied' by this structure.
				      */
    struct Pointerstruct {
	vector<int> *indices_ptr;
	vector<Vector<number> > *dof_values_ptr;
    };

				     /**
				      * Vector of all #Pointerstructs# (cf. there).
				      * It makes it
				      * easier to delete all these structs
				      * (without going over all #cell->user_pointer#)
				      * after they are not used any more, and
				      * collecting all these structures in a vector
				      * helps avoiding fraqmentation of the memory.
				      */
    vector<Pointerstruct> all_pointerstructs;

				     /**
				      * Is used for
				      * #prepare_for_refining_and_coarsening#
				      * The interpolated dof values
				      * of all cells that'll be coarsened
				      * will be stored in this vector.
				      */
    vector<vector<Vector<number> > > dof_values_on_cell;

				     /**
				      * After calling 
				      * #prepare_for_refinement_and_coarsening
				      * (vector<Vector<number> > &all_in)#
				      * this pointer points to the vector
				      * #all_in# for later comparison with
				      * the vector #all_out#
				      */
//    const vector<Vector<number> > * vecs_ptr;
};






/*----------------------------   solutiontransfer.h     ---------------------------*/
/* end of #ifndef __solutiontransfer_H */
#endif
/*----------------------------   solutiontransfer.h     ---------------------------*/
