/*----------------------------   dof_levels.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dof_levels_H
#define __dof_levels_H
/*----------------------------   dof_levels.h     ---------------------------*/



#include <vector>


/**
 * Store the indices of the degrees of freedom which are located on the lines.
 * Declare it to have a template parameter, but do not actually declare
 * other types than those explicitely instantiated.
 */
template <int N>
class DoFLevel;





/**
 * Store the indices of the degrees of freedom which are located on the lines.
 *
 * \subsection{Information for all #DoFLevel# classes}
 *
 * The #DoFLevel<N># classes 
 * store the global indices of the degrees of freedom for each cell on a
 * certain level. The index or number of a degree of freedom is the zero-based
 * index of the according value in the solution vector and the row and column
 * index in the global matrix or the multigrid matrix for this level. These
 * indices refer to the unconstrained vectors and matrices, where we have not
 * taken account of the constraints introduced by hanging nodes. If more than
 * one value corresponds to one basis function, for example for vector equations
 * where the solution is vector valued and thus has several degrees of freedom
 * for each basis function, we nonetheless store only one index. This can then
 * be viewed as the index into a block vector, where each block contains the
 * different values according to a degree of freedom. It is left to the derived
 * classes, whether the values in a block are stored consecutively or distributed
 * (e.g. if the solution function is $u=(u_1, u_2)$, we could store the values
 * in the solution vector like
 * $\ldots, u_1^m, u_2^m, u_1^{m+1}, u_2^{m+1},\ldots$ with $m$ denoting the
 * $m$th basis function, or $\ldots, u_1^m, u_1^{m+1}, u_1^{m+2}, \ldots,
 * u_2^m, u_2^{m+1}, u_2^{m+2}, \ldots$, respectively). Likewise, the
 * constraint matrix returned by #DoFHandler::make_constraint_matrix ()# is then
 * to be understood as a block matrix.
 *
 * The storage format of the degrees of freedom indices (short: DoF indices) is
 * somewhat like a mirror of the data structures of the triangulation classes.
 * There is a hierarchy of #DoFLevel<dim># classes for the different dimensions
 * which have objects named #line_dofs#, #quad_dofs# and so on, in which the
 * indices of DoFs located on lines and quads, respectively, are stored. The
 * indices are stored levelwise. The layout in
 * these arrays is as follows: if for a selected finite element (use
 * #DoFHandler::distribute_dofs()# to select a finite element) the number of
 * DoFs on each line (without those in the vertices) is #N#, then the length
 * of the #line_dofs# array is #N# times the number of lines on this level. The
 * DoF indices for the #i#th line are at the positions #N*i...(N+1)*i-1#.
 *
 * The DoF indices for vertices are not stored this way, since they need
 * different treatment in multigrid environments. If no multigrid is used, the
 * indices are stored in the #vertex_dofs# array of the #DoFHandler# class.
 */
class DoFLevel<1> {
  public:
				     /**
				      * Store the global indices of the degrees
				      * of freedom. See \Ref{DoFLevel} for
				      * detailed information.
				      */
    vector<int> line_dofs;
};




/**
 * Store the indices of the degrees of freedom which are located on quads.
 * See \Ref{DoFLevel<1>} for more information.
 */
class DoFLevel<2> : public DoFLevel<1> {
  public:
				     /**
				      * Store the global indices of the degrees
				      * of freedom. See \Ref{DoFLevel} for
				      * detailed information.
				      */
    vector<int> quad_dofs;
};



/*----------------------------   dof_levels.h     ---------------------------*/
/* end of #ifndef __dof_levels_H */
#endif
/*----------------------------   dof_levels.h     ---------------------------*/
