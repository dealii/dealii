/*----------------------------   solutiontransfer.h     ----------------------*/
/*      $Id$      */
/*           Ralf Hartmann, University of Heidelberg                          */
#ifndef __solutiontransfer_H
#define __solutiontransfer_H
/*----------------------------   solutiontransfer.h     ----------------------*/


#include <lac/forward-declarations.h>
#include <basic/forward-declarations.h>

#include <base/exceptions.h>
#include <vector>




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
 * soltrans.refine_interpolate(solution, interpolated_solution);
 *                                     // if necessary interpolate some
 *                                     // more functions
 * soltrans.refine_interpolate(sol2, interpolated_sol2);
 * ...
 * \end{verbatim}
 * \item If the grid will be coarsenend and refined
 * then use #SolutionTransfer# as follows
 * \begin{verbatim}
 * SolutionTransfer<dim, double> soltrans(*dof_handler);
 *                                     // some refinement e.g.
 * tria->refine_and_coarsen_fixed_fraction(error_indicator, 0.3, 0.05);
 *                                     // very important:
 * tria->prepare_coarsening_and_refinement();
 * soltrans.prepare_for_coarsening_and_refinement(solution);
 * tria->execute_coarsening_and_refinement ();
 * dof_handler->distribute_dofs (fe);
 * soltrans.interpolate(solution, interpolated_solution);
 * \end{verbatim}
 
 * Multiple calling of the function 
 * #interpolate (const Vector<number> &in, Vector<number> &out)#
 * is NOT allowed. Interpolating several functions can be performed in one step
 * by using #void interpolate (const vector<Vector<number> >&all_in,
 * vector<Vector<number> >&all_out) const#, and
 * using the respective #prepare_for_coarsening_and_refinement# function taking
 * several vectors as input before actually refining and coarsening the
 * triangulation (see there).
 * \end{itemize}
 *
 * For deleting all stored data in #SolutionTransfer# and reinitializing it
 * use the #clear()# function.
 *
 * Note that the #user_pointer# of some cells are used. 
 * Be sure that you don't need them otherwise.
 *
 * The template argument #number# denotes the data type of the vectors you want
 * to transfer.
 *
 *
 * \subsection{Implementation}
 *
 * \begin{itemize}
 * \item Solution transfer while pure refinement. Assume that we have got a
 * solution vector on the current (original) grid.
 * Each entry of this vector belongs to one of the
 * DoFs of the discretisation. If we now refine the grid then the calling of
 * #DoFHandler::distribute_dofs(...)# will change at least some of the
 * DoF indices. Hence we need to store the DoF indices of all active cells
 * before the refinement. The #user_pointer# of each active cell
 * is used to point to the vector of these DoF indices of that cell, all other
 * #user_pointers# are cleared. All this is
 * done by #prepare_for_pure_refinement()#.
 *
 * In the function #refine_interpolate(in,out)# and on each cell where the
 * #user_pointer# is set (i.e. the cells that were active in the original grid)
 * we can now access the local values of the solution vector #in#
 * on that cell by using the stored DoF indices. These local values are
 * interpolated and set into the vector #out# that is at the end the
 * discrete function #in# interpolated on the refined mesh.
 *
 * The #refine_interpolate(in,out)# function can be called multiplely for
 * arbitrary many discrete functions (solution vectors) on the original grid. 
 *
 * \item Solution transfer while coarsening and refinement. After 
 * calling #Triangulation::prepare_coarsening_and_refinement# the
 * coarsen flags of either all or none of the children of a 
 * (father-)cell are set. While coarsening 
 * (#Triangulation::execute_coarsening_and_refinement#)
 * the cells that are not needed any more will be deleted from the #Triangulation#.
 * 
 * For the interpolation from the (to be coarsenend) children to their father
 * the children cells are needed. Hence this interpolation
 * and the storing of the interpolated values of each of the discrete functions
 * that we want to interpolate needs to take place before these children cells
 * are coarsened (and deleted!!). Again the #user_pointers# of the cells are
 * set to point to these values (see below). 
 * Additionally the DoF indices of the cells
 * that will not be coarsened need to be stored according to the solution
 * transfer while pure refinement (cf there). All this is performed by
 * #prepare_for_coarsening_and_refinement(all_in)# where the 
 * #vector<Vector<number> >vector all_in# includes
 * all discrete functions to be interpolated onto the new grid.
 *
 * As we need two different kinds of pointers (#vector<int> *# for the Dof
 * indices and #vector<Vector<number> > *# for the interpolated DoF values)
 * we use the #Pointerstruct# that includes both of these pointers and
 * the #user_pointer# of each cell points to these #Pointerstructs#. 
 * On each cell only one of the two different pointers is used at one time 
 * hence we could use the
 * #void * user_pointer# as #vector<int> *# at one time and as 
 * #vector<Vector<number> > *# at the other but using this #Pointerstruct#
 * in between makes the use of these pointers more safe and gives better
 * possibility to expand their usage.
 * 
 * In #interpolate(all_in, all_out)# the refined cells are treated according
 * to the solution transfer while pure refinement. Additionally, on each
 * cell that is coarsened (hence previously was a father cell), 
 * the values of the discrete
 * functions in #all_out# are set to the stored local interpolated values 
 * that are accessible due to the 'vector<Vector<number> > *' pointer in 
 * #Pointerstruct# that is pointed to by the #user_pointer# of that cell.
 * It is clear that #interpolate(all_in, all_out)# only can be called with
 * the #vector<Vector<number> > all_in# that previously was the parameter
 * of the #prepare_for_coarsening_and_refinement(all_in)# function. Hence 
 * #interpolate(all_in, all_out)# can (in contrast to 
 * #refine_interpolate(in, out)#) only be called once.
 * \end{itemize}
 *
 * @author Ralf Hartmann, 1999
 */
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
				      * to be interpolated.
				      */
    void prepare_for_coarsening_and_refinement (const Vector<number> &in);
		      
				     /**
				      * This function
				      * interpolates the discrete function #in#,
				      * which is a vector on the grid before the
				      * refinement, to the function #out#
				      * which then is a vector on the refined grid.
				      * It assumes the vectors having the
				      * right sizes (i.e. #in.size()==n_dofs_old#,
				      * #out.size()==n_dofs_refined#)
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
				      * Same as #interpolate(in,out)#
				      * but it interpolates
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
				      * as in #all_in# as parameter of
				      * #prepare_for_refinement_and_coarsening(all_in)#.
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
				      * in one step by using
				      * #interpolate (all_out)#
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
		   << "although it should be " << arg2 << ".");
    
				     /**
				      * Exception
				      */
    DeclException0(ExcInternalError);
				  
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
				      * #PreparationState# denotes the
				      * three possible states of the
				      * #SolutionTransfer#: being
				      * prepared for 'pure refinement',
				      * prepared for 'coarsening and
				      * refinement' or not prepared.
				      */
    enum PreparationState {
	  none, pure_refinement, coarsening_and_refinement
    };

				     /**
				      * ???
				      */
    PreparationState prepared_for;


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
				      * is needed, hence one #user_pointer# should
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
};






/*----------------------------   solutiontransfer.h     ---------------------------*/
/* end of #ifndef __solutiontransfer_H */
#endif
/*----------------------------   solutiontransfer.h     ---------------------------*/
