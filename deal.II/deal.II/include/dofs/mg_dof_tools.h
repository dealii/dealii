/*----------------------------   mg_dof_tools.h     ---------------------------*/
/*      $Id$                 */
#ifndef __mg_dof_tools_H
#define __mg_dof_tools_H
/*----------------------------   mg_dof_tools.h     ---------------------------*/


#include <lac/forward_declarations.h>
#include <grid/forward_declarations.h>


/**
 * This is a collection of functions operating on, and manipulating
 * the numbers of degrees of freedom in a multilevel triangulation. It
 * is similar in purpose and function to the #DoFTools# class, but
 * operates on #MGDoFHandler# objects instead of #DoFHandler#
 * objects. See there and the documentation of the member functions
 * for more information.
 *
 * All member functions are static, so there is no need to create an
 * object of class #DoFTools#.
 *
 * @author Wolfgang Bangerth and others, 1999
 */
class MGDoFTools 
{
  public:
				     /**
				      * Write the sparsity structure
				      * of the matrix belonging to the
				      * specified #level# including
				      * constrained degrees of freedom
				      * into the matrix structure. The
				      * sparsity pattern does not
				      * include entries introduced by
				      * the elimination of constrained
				      * nodes.  The sparsity pattern
				      * is not compressed, since if
				      * you want to call
				      * #ConstraintMatrix::condense(1)#
				      * afterwards, new entries have
				      * to be added. However, if you
				      * don't want to call
				      * #ConstraintMatrix::condense(1)#,
				      * you have to compress the
				      * matrix yourself, using
				      * #SparseMatrixStruct::compress()#.
				      */
    template <int dim>
    static void
    make_sparsity_pattern (const MGDoFHandler<dim> &dof_handler,
			   const unsigned int       level,
			   SparseMatrixStruct      &sparsity);
};



/*----------------------------   mg_dof_tools.h     ---------------------------*/
/* end of #ifndef __mg_dof_tools_H */
#endif
/*----------------------------   mg_dof_tools.h     ---------------------------*/
