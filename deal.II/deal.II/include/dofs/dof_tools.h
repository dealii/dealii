// $Id$
// Copyright Guido Kanschat, 1999

#ifndef __deal_dof_tools_H
#define __deal_dof_tools_H

#include <bvector.h>

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
				      * Extract DoFs of components.
				      * The bit vector #select#
				      * defines, which components of an
				      * #FESystem# are to be extracted
				      * from the DoFHandler #dof#. The
				      * numbers of these dofs are then
				      * entered consecutively into
				      * #selected_dofs#.
				      *
				      * A prior ordering of dof-values
				      * is not destroyed by this process.
				      */
    template<int dim>
    static void extract_dofs(const DoFHandler<dim>& dof,
			     const bit_vector& select,
			     bit_vector& selected_dofs);
};


#endif
