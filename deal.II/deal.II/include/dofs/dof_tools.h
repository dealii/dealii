/*----------------------------   dof_tools.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dof_tools_H
#define __dof_tools_H
/*----------------------------   dof_tools.h     ---------------------------*/
// Copyright Guido Kanschat, 1999


#include <vector>




/**
 * Operations on DoF-numbers.
 * This is a collection of functions manipulating the numbers of
 * degrees of freedom. The documentation of the member functions will
 * provide more information.
 *
 * All member functions are static, so there is no need to create an
 * object of class #DoFTools#.
 * @author Guido Kanschat, 1999
 */
class DoFTools
{
  public:
				     /**
				      * Extract the indices of the degrees
				      * of freedom belonging to certain
				      * components. The bit vector #select#
				      * defines, which components of an
				      * #FESystem# are to be extracted
				      * from the DoFHandler #dof#. The
				      * respective entries in #selected_dofs#
				      * are then flagged #true#, while all
				      * others are set to #false#.
				      *
				      * The size of #select# shall equal
				      * the number of components in the
				      * finite element used by #dof#.
				      */
    template<int dim>
    static void extract_dofs(const DoFHandler<dim> &dof,
			     const vector<bool>    &select,
			     vector<bool>          &selected_dofs);

				     /**
				      * Do the same thing as #extract_dofs#
				      * for one level of a multi-grid DoF
				      * numbering.
				      */
    template<int dim>
    static void extract_level_dofs(const unsigned int       level,
				   const MGDoFHandler<dim> &dof,
				   const vector<bool>      &select,
				   vector<bool>            &selected_dofs);
};




/*----------------------------   dof_tools.h     ---------------------------*/
/* end of #ifndef __dof_tools_H */
#endif
/*----------------------------   dof_tools.h     ---------------------------*/
