/*----------------------------   dof_tools.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dof_tools_H
#define __dof_tools_H
/*----------------------------   dof_tools.h     ---------------------------*/
// Copyright Wolfgang Bangerth, Guido Kanschat, and others 1999


#include <lac/forward-declarations.h>
#include <vector>



// This is part of the original documentation of DoFHandler
// I tried to put the essential part into the documentation of
// the function, where it belongs. Still, I keep it for a while,
// to be able to recover missing information.
/*
 * \subsection{Setting up sparsity patterns}
 *
 * When assembling system matrices, the entries are usually of the form
 * $a_{ij} = a(\phi_i, \phi_j)$, where $a$ is a bilinear functional, often an
 * integral. When using sparse matrices, we therefore only need to reserve space
 * for those $a_{ij}$ only, which are nonzero, which is the same as to say that
 * that the basis functions $\phi_i$ and $\phi_j$ have a nonempty intersection of
 * their support. Since the support of basis functions is bound only on cells
 * on which they are located or to which they are adjacent, to determine the
 * sparsity pattern it is sufficient to loop over all cells and connect all
 * basis functions on each cell with all other basis functions on that cell.
 * There may be finite elements for which not all basis functions on a cell
 * connect with each other, but no use of this case is made since no examples
 * where this occurs are known to the author.
 *
 * When setting up sparsity patterns for matrices on the boundary, the same
 * procedure is done, except for the fact that the loop only goes over faces
 * on the boundary and the basis functions thereon. It is assumed that all
 * other basis functions on a cell adjacent to the boundary vanish at the
 * boundary itself, except for those which are located on the boundary.
 *
 */

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
				      * Locate non-zero entries of the system matrix.
				      *
				      * This function computes the possible
				      * positions of non-zero entries
				      * in the global system
				      * matrix. We assume that a
				      * finite element basis function
				      * is non-zero on a cell only if
				      * its dual degree of freedom is
				      * associated with the interior,
				      * a face, an edge or a vertex of
				      * this cell. As a result, the
				      * matrix entry between two basis
				      * functions can be non-zero only
				      * if they correspond to degrees
				      * of freedom of at least one
				      * common
				      * cell. Therefore,
				      * #make_sparsity_pattern# just
				      * loops over all cells and
				      * enters all couplings local to
				      * that cell.
				      *
				      * Since this process is purely
				      * local, the sparsity pattern
				      * does not provide for entries introduced by
				      * the elimination of hanging nodes.
				      * They have to be taken care of
				      * by a call to 
				      * #ConstraintMatrix::condense()#
				      * afterwards.
				      *
				      * Remember using
				      * #SparseMatrixStruct::compress()#
				      * after generating the pattern.
				      */
    template<int dim>
    static void make_sparsity_pattern (const DoFHandler<dim>& dof,
				       SparseMatrixStruct &);

				     /**
				      * Locate non-zero entries for
				      * mixed methods.
				      * This function does mostly the same as
				      * the other
				      * #make_sparsity_pattern#, but
				      * it is specialized for
				      * mixed finite elements and allows to
				      * specify which variables couple in which
				      * equation. For example, if wanted to solve
				      * the Stokes equations,
				      * \begin{verbatim}
				      *    -\Delta \vec u + \nabla p = 0,
				      *    \div u                    = 0
				      * \end{verbatim}
				      * in, two space dimensions, using
				      * stable Q2/Q1 mixed elements
				      * (using the #FESystem# class), then
				      * you don't want all degrees of freedom
				      * to couple in each equation. You rather
				      * may want to give the following pattern
				      * of couplings:
				      * \begin{verbatim}
				      *    1  0  1
				      *    0  1  1
				      *    1  1  0
				      * \end{verbatim}
				      * where "1" indicates that two variables
				      * (i.e. components of the #FESystem#)
				      * couple in the respective equation, and
				      * a "0" means no coupling, in which case
				      * it is not necessary to allocate space
				      * in the matrix structure. Obviously, the
				      * mask refers to components of the
				      * composed #FESystem#, rather than to the
				      * degrees of freedom contained in there.
				      *
				      * This function is designed to accept
				      * a mask, like the one shown above,
				      * through the #mask# parameter, which
				      * contains boolean values. It builds
				      * the matrix structure just like the
				      * previous function, but does not create
				      * elements if not specified by the mask.
				      */
    template<int dim>
    static void make_sparsity_pattern (const DoFHandler<dim>& dof,
				       const vector<vector<bool> > &mask,
				       SparseMatrixStruct          &);

    				     /**
				      * Write the sparsity structure of the
				      * matrix composed of the basis functions
				      * on the boundary into the
				      * matrix structure. The sparsity pattern
				      * does not include entries introduced by
				      * the elimination of constrained nodes.
				      * The sparsity
				      * pattern is not compressed, since if
				      * you want to call
				      * #ConstraintMatrix::condense(1)#
				      * afterwards, new entries have to be
				      * added. However, if you don't want to call
				      * #ConstraintMatrix::condense(1)#, you
				      * have to compress the matrix yourself,
				      * using #SparseMatrixStruct::compress()#.
				      *
				      * Since this function is obviously useless
				      * in one spatial dimension, it is not
				      * implemented.
				      */
    template<int dim>
    static void
    make_boundary_sparsity_pattern (const DoFHandler<dim>& dof,
				    const vector<int> &dof_to_boundary_mapping,
				    SparseMatrixStruct &); 

				     /**
				      * Write the sparsity structure of the
				      * matrix composed of the basis functions
				      * on the boundary into the
				      * matrix structure. In contrast to the
				      * previous function, only those parts
				      * of the boundary are considered of which
				      * the boundary indicator is listed in the
				      * set of numbers passed to this function.
				      *
				      * In fact, rather than a #set# of boundary
				      * indicators, a #map# needs to be passed,
				      * since most of the functions handling with
				      * boundary indicators take a mapping of
				      * boundary indicators and the respective
				      * boundary functions. The boundary function,
				      * however, is ignored in this function.
				      * If you have no functions at hand, but only
				      * the boundary indicators, set the function
				      * pointers to null pointers.
				      *
				      * Since this function is obviously useless
				      * in one spatial dimension, it is not
				      * implemented.
				      */
    template<int dim>
    static void
    make_boundary_sparsity_pattern (const DoFHandler<dim>& dof,
				    const typename DoFHandler<dim>::FunctionMap &boundary_indicators,
				    const vector<int>  &dof_to_boundary_mapping,
				    SparseMatrixStruct &sparsity); 

				     /**
				      * Generate sparsity pattern for fluxes.
				      * This is a replacement of the
				      * function
				      * #make_sparsity_pattern# for
				      * discontinuous methods. Since
				      * the fluxes include couplings
				      * between neighboring elements,
				      * the normal couplings and these
				      * extra matrix entries are
				      * considered.
				      */
    template<int dim>
    static void make_flux_sparsity_pattern (const DoFHandler<dim>&,
					    SparseMatrixStruct &);
    
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
