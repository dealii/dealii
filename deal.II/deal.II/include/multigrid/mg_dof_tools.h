/*----------------------------   mg_dof_tools.h     ---------------------------*/
/*      $Id$                 */
#ifndef __mg_dof_tools_H
#define __mg_dof_tools_H
/*----------------------------   mg_dof_tools.h     ---------------------------*/


#include <lac/forward_declarations.h>
#include <grid/forward_declarations.h>

//TODO: Consider incorporating these functions in DoFTools

/**
 * This is a collection of functions operating on, and manipulating
 * the numbers of degrees of freedom in a multilevel triangulation. It
 * is similar in purpose and function to the #DoFTools# class, but
 * operates on #MGDoFHandler# objects instead of #DoFHandler#
 * objects. See there and the documentation of the member functions
 * for more information.
 *
 * All member functions are static, so there is no need to create an
 * object of class #MGDoFTools#.
 *
 * @author Wolfgang Bangerth and others, 1999
 */
class MGDoFTools 
{
  public:
				     /**
				      * Write the sparsity structure
				      * of the matrix belonging to the
				      * specified #level#. The sparsity pattern
				      * is not compressed, so before 
				      * creating the actual matrix
				      * you have to compress the
				      * matrix yourself, using
				      * #SparseMatrixStruct::compress()#.
				      *
				      * There is no need to consider
				      * hanging nodes here, since only
				      * one level is considered.
				      */
    template <int dim>
    static void
    make_sparsity_pattern (const MGDoFHandler<dim> &dof_handler,
			   SparsityPattern         &sparsity,
			   const unsigned int       level);

				     /**
				      * Make a sparsity pattern including fluxes
				      * of discontinuous Galerkin methods.
				      * @see make_sparsity_pattern
				      * $see DoFTools
				      */
    template <int dim>
    static void
    make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof_handler,
				SparsityPattern         &sparsity,
				const unsigned int       level);


};



/*----------------------------   mg_dof_tools.h     ---------------------------*/
/* end of #ifndef __mg_dof_tools_H */
#endif
/*----------------------------   mg_dof_tools.h     ---------------------------*/
