//----------------------------  dof_tools.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools.h  ---------------------------
#ifndef __deal2__dof_tools_h
#define __deal2__dof_tools_h


// Copyright Wolfgang Bangerth, Guido Kanschat, and others 1999, 2000


#include <base/exceptions.h>
#include <lac/forward_declarations.h>
#include <grid/forward_declarations.h>
#include <vector>


/**
 * This is a collection of functions operating on, and manipulating
 * the numbers of degrees of freedom. The documentation of the member
 * functions will provide more information, but for functions that
 * exist in multiple versions, there are sections in this global
 * documentation stating some commonalities.
 *
 * All member functions are static, so there is no need to create an
 * object of class #DoFTools#.
 *
 * \subsection{Setting up sparsity patterns}
 *
 * When assembling system matrices, the entries are usually of the form
 * $a_{ij} = a(\phi_i, \phi_j)$, where $a$ is a bilinear functional, often an
 * integral. When using sparse matrices, we therefore only need to reserve space
 * for those $a_{ij}$ only, which are nonzero, which is the same as to say that
 * the basis functions $\phi_i$ and $\phi_j$ have a nonempty intersection of
 * their support. Since the support of basis functions is bound only on cells
 * on which they are located or to which they are adjacent, to
 * determine the sparsity pattern it is sufficient to loop over all
 * cells and connect all basis functions on each cell with all other
 * basis functions on that cell.  There may be finite elements for
 * which not all basis functions on a cell connect with each other,
 * but no use of this case is made since no examples where this occurs
 * are known to the author.
 *
 * When setting up sparsity patterns for matrices on the boundary, the same
 * procedure is done, except for the fact that the loop only goes over faces
 * on the boundary and the basis functions thereon. It is assumed that all
 * other basis functions on a cell adjacent to the boundary vanish at the
 * boundary itself, except for those which are located on the boundary.
 *
 * @author Wolfgang Bangerth and others, 1998, 1999, 2000
 */
class DoFTools
{
  public:
				     /**
				      * Locate non-zero entries of the
				      * system matrix.
				      *
				      * This function computes the
				      * possible positions of non-zero
				      * entries in the global system
				      * matrix. We assume that a
				      * certain finite element basis
				      * function is non-zero on a cell
				      * only if its degree of freedom
				      * is associated with the
				      * interior, a face, an edge or a
				      * vertex of this cell. As a
				      * result, the matrix entry
				      * between two basis functions
				      * can be non-zero only if they
				      * correspond to degrees of
				      * freedom of at least one common
				      * cell. Therefore,
				      * #make_sparsity_pattern# just
				      * loops over all cells and
				      * enters all couplings local to
				      * that cell. As the generation
				      * of the sparsity pattern is
				      * irrespective of the equation
				      * which is solved later on, the
				      * resulting sparsity pattern is
				      * symmetric.
				      *
				      * Since this process is purely
				      * local, the sparsity pattern
				      * does not provide for entries
				      * introduced by the elimination
				      * of hanging nodes.  They have
				      * to be taken care of by a call
				      * to
				      * #ConstraintMatrix::condense()#
				      * afterwards.
				      *
				      * Remember using
				      * #SparsityPattern::compress()#
				      * after generating the pattern.
				      *
				      * The actual type of the
				      * sparsity pattern may be
				      * #SparsityPattern#,
				      * #BlockSparsityPattern#, or any
				      * other class that satisfies
				      * similar requirements.
				      */
    template <int dim, class SparsityPattern>
    static void make_sparsity_pattern (const DoFHandler<dim> &dof,
				       SparsityPattern       &sparsity_pattern);

				     /**
				      * Locate non-zero entries for
				      * mixed methods.  This function
				      * does mostly the same as the
				      * other #make_sparsity_pattern#,
				      * but it is specialized for
				      * mixed finite elements and
				      * allows to specify which
				      * variables couple in which
				      * equation. For example, if
				      * wanted to solve the Stokes
				      * equations,
				      *
				      *
				      * \begin{verbatim}
				      * -\Delta \vec u + \nabla p = 0,
				      * \div u                    = 0
				      * \end{verbatim}
				      *
				      * in two space dimensions,
				      * using stable Q2/Q1 mixed
				      * elements (using the #FESystem#
				      * class), then you don't want
				      * all degrees of freedom to
				      * couple in each equation. You
				      * rather may want to give the
				      * following pattern of
				      * couplings:
				      *
				      * \begin{verbatim}
				      *   1 0 1
				      *   0 1 1
				      *   1 1 0
				      * \end{verbatim}
				      * where "1" indicates that two
				      * variables (i.e. components of
				      * the #FESystem#) couple in the
				      * respective equation, and a "0"
				      * means no coupling, in which
				      * case it is not necessary to
				      * allocate space in the matrix
				      * structure. Obviously, the mask
				      * refers to components of the
				      * composed #FESystem#, rather
				      * than to the degrees of freedom
				      * contained in there.
				      *
				      * This function is designed to
				      * accept a mask, like the one
				      * shown above, through the
				      * #mask# parameter, which
				      * contains boolean values. It
				      * builds the matrix structure
				      * just like the previous
				      * function, but does not create
				      * elements if not specified by
				      * the mask. If the mask is
				      * symmetric, then so will be the
				      * resulting sparsity pattern.
				      *
				      * The actual type of the
				      * sparsity pattern may be
				      * #SparsityPattern#,
				      * #BlockSparsityPattern#, or any
				      * other class that satisfies
				      * similar requirements.
				      */
    template<int dim, class SparsityPattern>
    static void make_sparsity_pattern (const DoFHandler<dim>       &dof,
				       const vector<vector<bool> > &mask,
				       SparsityPattern             &sparsity_pattern);

    				     /**
				      * Write the sparsity structure
				      * of the matrix composed of the
				      * basis functions on the
				      * boundary into the matrix
				      * structure. The sparsity
				      * pattern does not include
				      * entries introduced by the
				      * elimination of constrained
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
				      * #SparsityPattern::compress()#.
				      *
				      * Since this function is
				      * obviously useless in one
				      * spatial dimension, it is not
				      * implemented.
				      */
    template<int dim>
    static void
    make_boundary_sparsity_pattern (const DoFHandler<dim>      &dof,
				    const vector<unsigned int> &dof_to_boundary_mapping,
				    SparsityPattern            &sparsity_pattern); 

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
				      * In fact, rather than a #set#
				      * of boundary indicators, a
				      * #map# needs to be passed,
				      * since most of the functions
				      * handling with boundary
				      * indicators take a mapping of
				      * boundary indicators and the
				      * respective boundary
				      * functions. The boundary
				      * function, however, is ignored
				      * in this function.  If you have
				      * no functions at hand, but only
				      * the boundary indicators, set
				      * the function pointers to null
				      * pointers.
				      *
				      * Since this function is
				      * obviously useless in one
				      * spatial dimension, it is not
				      * implemented.
				      */
    template<int dim>
    static void
    make_boundary_sparsity_pattern (const DoFHandler<dim>& dof,
				    const typename DoFHandler<dim>::FunctionMap &boundary_indicators,
				    const vector<unsigned int>  &dof_to_boundary_mapping,
				    SparsityPattern    &sparsity); 

				     /**
				      * Generate sparsity pattern for
				      * fluxes, i.e. formulations of
				      * the discrete problem with
				      * discontinuous elements which
				      * couple across faces of cells.
				      * This is a replacement of the
				      * function
				      * #make_sparsity_pattern# for
				      * discontinuous methods. Since
				      * the fluxes include couplings
				      * between neighboring elements,
				      * the normal couplings and these
				      * extra matrix entries are
				      * considered.
				      *
				      * This function uses the user
				      * flags of the triangulation.
				      */
    template<int dim>
    static void
    make_flux_sparsity_pattern (const DoFHandler<dim> &dof_handler,
				SparsityPattern       &sparsity_pattern);
    
				     /**
				      * Make up the constraints which
				      * is result from the use of hanging
				      * nodes. The object into which these
				      * are inserted is later
				      * used to condensate the global
				      * system matrices and to prolong
				      * the solution vectors from the true
				      * degrees of freedom also to the
				      * constraint nodes.
				      *
				      * Since this method does not make sense in
				      * one dimension, the function returns
				      * immediately. The object is not cleared
				      * before use, so you should make sure
				      * it containts only constraints you still
				      * want; otherwise call the #clear#
				      * function.
				      *
				      * To condense a given sparsity pattern,
				      * use #ConstraintMatrix::condense#.
				      * Before doing so, you need to close
				      * the constraint object, which must be
				      * done after all constraints are entered.
				      * This function does not close the object
				      * since you may want to enter other
				      * constraints later on yourself.
				      */
    template <int dim>
    static void
    make_hanging_node_constraints (const DoFHandler<dim> &dof_handler,
				   ConstraintMatrix      &constraints);

				     /**
				      * Take a vector of values which live on
				      * cells (e.g. an error per cell) and
				      * distribute it to the dofs in such a
				      * way that a finite element field results,
				      * which can then be further processed,
				      * e.g. for output. You should note that
				      * the resulting field will not be
				      * continuous at hanging nodes. This can,
				      * however, easily be arranged by calling
				      * the appropraite #distribute# function
				      * of a #ConstraintMatrix# object created
				      * for this #DoFHandler# object.
				      *
				      * It is assumed that the number of
				      * elements in #cell_data# equals the
				      * number of active cells. The size of
				      * #dof_data# is adjusted to the right
				      * size.
				      *
				      * Note that the input vector may be
				      * a vector of any data type as long
				      * as it is convertible to #double#.
				      * The output vector, being a data
				      * vector on the grid, always consists
				      * of elements of type #double#.
				      *
				      * In case the finite element used by
				      * this DoFHandler consists of more than
				      * one component, you should give which
				      * component in the output vector should
				      * be used to store the finite element
				      * field in; the default is zero (no other
				      * value is allowed if the finite element
				      * consists only of one component). All
				      * other components of the vector remain
				      * untouched, i.e. their contents are
				      * not changed.
				      *
				      * It is assumed that the output vector
				      * #dof_data# already has the right size,
				      * i.e. #n_dofs()# elements.
				      */
    template <int dim, typename Number>
    static void
    distribute_cell_to_dof_vector (const DoFHandler<dim> &dof_handler,
				   const Vector<Number>  &cell_data,
				   Vector<double>        &dof_data,
				   const unsigned int     component = 0);

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
				      * The size of #component_select#
				      * shall equal the number of
				      * components in the finite
				      * element used by #dof#. The
				      * size of #selected_dofs# shall
				      * equal
				      * #dof_handler.n_dofs()#. Previous
				      * contents of this array or
				      * overwritten.
				      */
    template <int dim>
    static void
    extract_dofs (const DoFHandler<dim> &dof_handler,
		  const vector<bool>    &component_select,
		  vector<bool>          &selected_dofs);

				     /**
				      * Do the same thing as
				      * #extract_dofs# for one level
				      * of a multi-grid DoF numbering.
				      */
    template <int dim>
    static void
    extract_level_dofs (const unsigned int       level,
			const MGDoFHandler<dim> &dof,
			const vector<bool>      &select,
			vector<bool>            &selected_dofs);

				     /**
				      * Extract all degrees of freedom
				      * which are at the boundary and
				      * belong to specified components
				      * of the solution. The function
				      * returns its results in the
				      * last parameter which contains
				      * #true# is a degree of freedom
				      * is at the boundary and belongs
				      * to one of the selected
				      * components, and #false#
				      * otherwise.
				      *
				      * The size of #component_select#
				      * shall equal the number of
				      * components in the finite
				      * element used by #dof#. The
				      * size of #selected_dofs# shall
				      * equal
				      * #dof_handler.n_dofs()#. Previous
				      * contents of this array or
				      * overwritten.
				      */
    template <int dim>
    static void
    extract_boundary_dofs (const DoFHandler<dim> &dof_handler,
			   const vector<bool>    &component_select,
			   vector<bool>          &selected_dofs);

				     /**
				      * Select all dofs that will be
				      * constrained by interface
				      * constraints, i.e. all hanging
				      * nodes.
				      *
				      * The size of #selected_dofs#
				      * shall equal
				      * #dof_handler.n_dofs()#. Previous
				      * contents of this array or
				      * overwritten.
				      */
    template <int dim>
    static void
    extract_hanging_node_dofs (const DoFHandler<dim> &dof_handler,
			       vector<bool>          &selected_dofs);
    
				     /**
				      * This function can be used when
				      * different variables shall be
				      * discritized on different
				      * grids, where one grid is
				      * coarser than the other. This
				      * idea might seem nonsensical at
				      * first, but has reasonable
				      * applications in inverse
				      * (parameter estimation)
				      * problems, where there might
				      * not be enough information to
				      * recover the parameter on the
				      * same grid as the state
				      * variable; furthermore, the
				      * smoothness properties of state
				      * variable and parameter might
				      * not be too much related, so
				      * using different grids might be
				      * an alternative to using
				      * stronger regularization of the
				      * problem.
				      *
				      * The basic idea of this
				      * function is explained in the
				      * following. Let us, for
				      * convenience, denote by
				      * ``parameter grid'' the coarser
				      * of the two grids, and by
				      * ``state grid'' the finer of
				      * the two. We furthermore assume
				      * that the finer grid can be
				      * obtained by refinement of the
				      * coarser one, i.e. the fine
				      * grid is at least as much
				      * refined as the coarse grid at
				      * each point of the
				      * domain. Then, each shape
				      * function on the coarse grid
				      * can be represented as a linear
				      * combination of shape functions
				      * on the fine grid (assuming
				      * identical ansatz
				      * spaces). Thus, if we
				      * discretize as usual, using
				      * shape functions on the fine
				      * grid, we can consider the
				      * restriction that the parameter
				      * variable shall in fact be
				      * discretized by shape functions
				      * on the coarse grid as a
				      * constraint. These constraints
				      * are linear and happen to have
				      * the form managed by the
				      * ``ConstraintMatrix'' class.
				      *
				      * The construction of these
				      * constraints is done as
				      * follows: for each of the
				      * degrees of freedom (i.e. shape
				      * functions) on the coarse grid,
				      * we compute its representation
				      * on the fine grid, i.e. how the
				      * linear combination of shape
				      * functions on the fine grid
				      * looks like that resembles the
				      * shape function on the coarse
				      * grid. From this information,
				      * we can then compute the
				      * constraints which have to hold
				      * if a solution of a linear
				      * equation on the fine grid
				      * shall be representable on the
				      * coarse grid. The exact
				      * algorithm how these
				      * constraints can be computed is
				      * rather complicated and is best
				      * understood by reading the
				      * source code, which contains
				      * many comments.
				      *
				      * Before explaining the use of
				      * this function, we would like
				      * to state that the total number
				      * of degrees of freedom used for
				      * the discretization is not
				      * reduced by the use of this
				      * function, i.e. even though we
				      * discretize one variable on a
				      * coarser grid, the total number
				      * of degrees of freedom is that
				      * of the fine grid. This seems
				      * to be counter-productive,
				      * since it does not give us a
				      * benefit from using a coarser
				      * grid. The reason why it may be
				      * useful to choose this approach
				      * nonetheless is three-fold:
				      * first, as stated above, there
				      * might not be enough
				      * information to recover a
				      * parameter on a fine grid,
				      * i.e. we chose to discretize it
				      * on the coarse grid not to save
				      * DoFs, but for other
				      * reasons. Second, the
				      * ``ConstraintMatrix'' includes
				      * the constraints into the
				      * linear system of equations, by
				      * which constrained nodes become
				      * dummy nodes; we may therefore
				      * exclude them from the linear
				      * algebra, for example by
				      * sorting them to the back of
				      * the DoF numbers and simply
				      * calling the solver for the
				      * upper left block of the matrix
				      * which works on the
				      * non-constrained nodes only,
				      * thus actually realizing the
				      * savings in numerical effort
				      * from the reduced number of
				      * actual degrees of freedom. The
				      * third reason is that for some
				      * or other reason we have chosen
				      * to use two different grids, it
				      * may be actually quite
				      * difficult to write a function
				      * that assembles the system
				      * matrix for finite element
				      * spaces on different grids;
				      * using the approach of
				      * constraints as with this
				      * function allows to use
				      * standard techniques when
				      * discretizing on only one grid
				      * (the finer one) without having
				      * to take care of the fact that
				      * one or several of the variable
				      * actually belong to different
				      * grids.
				      *
				      * The use of this function is as
				      * follows: it accepts as
				      * parameters two DoF Handlers,
				      * the first of which refers to
				      * the coarse grid and the second
				      * of which is the fine grid. On
				      * both, a finite element is
				      * represented by the DoF handler
				      * objects, which will usually
				      * have several components, which
				      * may belong to different finite
				      * elements. The second and
				      * fourth parameter of this
				      * function therefore state which
				      * variable on the coarse grid
				      * shall be used to restrict the
				      * stated component on the fine
				      * grid. Of course, the finite
				      * elements used for the
				      * respective components on the
				      * two grids need to be the
				      * same. An example may clarify
				      * this: consider the parameter
				      * estimation mentioned briefly
				      * above; there, on the fine grid
				      * the whole discretization is
				      * done, thus the variables are
				      * ``u'', ``q'', and the Lagrange
				      * multiplier ``lambda'', which
				      * are discretized using
				      * continuous linear, piecewise
				      * constant discontinuous, and
				      * continuous linear elements,
				      * respectively. Only the
				      * parameter ``q'' shall be
				      * represented on the coarse
				      * grid, thus the DoFHandler
				      * object on the coarse grid
				      * represents only one variable,
				      * discretized using piecewise
				      * constant discontinuous
				      * elements. Then, the parameter
				      * denoting the component on the
				      * coarse grid would be zero (the
				      * only possible choice, since
				      * the variable on the coarse
				      * grid is scalar), and one on
				      * the fine grid (corresponding
				      * to the variable ``q''; zero
				      * would be ``u'', two would be
				      * ``lambda''). Furthermore, an
				      * object of type #IntergridMap#
				      * is needed; this could in
				      * principle be generated by the
				      * function itself from the two
				      * DoFHandler objects, but since
				      * it is probably available
				      * anyway in programs that use
				      * this function, we shall use it
				      * instead of re-generating
				      * it. Finally, the computed
				      * constraints are entered into a
				      * variable of type
				      * #ConstraintMatrix#; the
				      * constraints are added,
				      * i.e. previous contents which
				      * may have, for example, be
				      * obtained from hanging nodes,
				      * are not deleted, so that you
				      * only need one object of this
				      * type.
				      */
    template <int dim>
    static void
    compute_intergrid_constraints (const DoFHandler<dim>              &coarse_grid,
				   const unsigned int                  coarse_component,
				   const DoFHandler<dim>              &fine_grid,
				   const unsigned int                  fine_component,
				   const InterGridMap<DoFHandler,dim> &coarse_to_fine_grid_map,
				   ConstraintMatrix                   &constraints);

				   
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongSize,
		    int, int,
		    << "The dimension " << arg1 << " of the vector is wrong. "
		    << "It should be " << arg2);
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidComponent,
		    int, int,
		    << "The component you gave (" << arg1 << ") "
		    << "is invalid with respect to the number "
		    << "of components in the finite element "
		    << "(" << arg2 << ")");
				     /**
				      * Exception
				      */
    DeclException0 (ExcFiniteElementsDontMatch);
				     /**
				      * Exception
				      */
    DeclException0 (ExcGridNotCoarser);
				     /**
				      * Exception
				      */
    DeclException0 (ExcGridsDontMatch);
};


#endif
