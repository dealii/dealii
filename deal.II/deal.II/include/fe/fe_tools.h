//----------------------------  fe_tools.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools.h  ---------------------------
#ifndef __deal2__fe_tools_H
#define __deal2__fe_tools_H



template <typename number> class FullMatrix;
template <int dim> class FiniteElement;
template <int dim> class DoFHandler;
template <typename number> class Vector;
template <int dim> class FE_Q;

#include <base/config.h>
#include <base/exceptions.h>

/**
 * This class performs interpolations and extrapolations of discrete
 * functions of one @p{FiniteElement} @p{fe1} to another @p{FiniteElement}
 * @p{fe2}.
 *
 * It also provides the local interpolation matrices that interpolate
 * on each cell. Furthermore it provides the difference matrix
 * $id-I_h$ that is needed for evaluating $(id-I_h)z$ for e.g. the
 * dual solution $z$.
 *
 * @author Ralf Hartmann, 2000
 */
class FETools
{
  public:
				     /**
				      * Gives the interpolation matrix
				      * that interpolates a @p{fe1}-
				      * function to a @p{fe2}-function on
				      * each cell. The interpolation_matrix
				      * needs to be of size
				      * @p{(fe2.dofs_per_cell, fe1.dofs_per_cell)}.
				      *
				      * Note, that if the finite element
				      * space @p{fe1} is a subset of
				      * the finite element space
				      * @p{fe2} than the @p{interpolation_matrix}
				      * is an embedding matrix.
				      */
    template <int dim, typename number>
    static void get_interpolation_matrix(const FiniteElement<dim> &fe1,
					 const FiniteElement<dim> &fe2,
					 FullMatrix<number> &interpolation_matrix);
    
				     /**
				      * Gives the interpolation matrix
				      * that interpolates a @p{fe1}-
				      * function to a @p{fe2}-function, and
				      * interpolates this to a second
				      * @p{fe1}-function on
				      * each cell. The interpolation_matrix
				      * needs to be of size
				      * @p{(fe1.dofs_per_cell, fe1.dofs_per_cell)}.
				      *
				      * Note, that this function only
				      * makes sense if the finite element
				      * space due to @p{fe1} is not a subset of
				      * the finite element space due to
				      * @p{fe2}, as if it were a subset then
				      * the @p{interpolation_matrix} would be 
				      * only the unit matrix.
				      */
    template <int dim, typename number>
    static void get_back_interpolation_matrix(const FiniteElement<dim> &fe1,
					      const FiniteElement<dim> &fe2,
					      FullMatrix<number> &interpolation_matrix);

				     /**
				      * Gives the unit matrix minus the
				      * back interpolation matrix.
				      * The @p{difference_matrix}
				      * needs to be of size
				      * @p{(fe1.dofs_per_cell, fe1.dofs_per_cell)}.
				      *
				      * This function gives
				      * the matrix that transforms a
				      * @p{fe1} function $z$ to $z-I_hz$
				      * where $I_h$ denotes the interpolation
				      * operator from the @p{fe1} space to
				      * the @p{fe2} space. This matrix hence
				      * is useful to evaluate
				      * error-representations where $z$
				      * denotes the dual solution.
				      */
    template <int dim, typename number>
    static void get_interpolation_difference_matrix(
      const FiniteElement<dim> &fe1,
      const FiniteElement<dim> &fe2,
      FullMatrix<number> &difference_matrix);

				     /**
				      * Gives the interpolation of a the
				      * @p{dof1}-function @p{u1} to a
				      * @p{dof2}-function @p{u2}. @p{dof1} and
				      * @p{dof2} need to be @ref{DoFHandler}s
				      * based on the same triangulation.
				      *
				      * If the elements @p{fe1} and @p{fe2}
				      * are either both continuous or
				      * both discontinuous then this
				      * interpolation is the usual point
				      * interpolation. The same is true
				      * if @p{fe1} is a continuous and
				      * @p{fe2} is a discontinuous finite
				      * element. For the case that @p{fe1}
				      * is a discontinuous and @p{fe2} is
				      * a continuous finite element
				      * there is no point interpolation
				      * defined at the discontinuities.
				      * Therefore the meanvalue is taken
				      * at the DoF values on the
				      * discontinuities.
				      */
    template <int dim, typename number>
    static void interpolate(const DoFHandler<dim> &dof1,
			    const Vector<number> &u1,
			    const DoFHandler<dim> &dof2,
			    Vector<number> &u2);
    
				     /**
				      * Gives the interpolation of the @p{fe1}-
				      * function @p{u1} to a @p{fe2}-function, and
				      * interpolates this to a second
				      * @p{fe1}-function named @p{u1_interpolated}.
				      *
				      * Note, that this function only
				      * makes sense if the finite element
				      * space due to @p{fe1} is not a subset of
				      * the finite element space due to
				      * @p{fe2}, as if it were a subset then
				      * @p{u1_interpolated} would be equal to @p{u1}.
				      */
    template <int dim, typename number>
    static void back_interpolate(const DoFHandler<dim> &dof1,
				 const Vector<number> &u1,
				 const FiniteElement<dim> &fe2,
				 Vector<number> &u1_interpolated);

				     /**
				      * Gives $(Id-I_h)z2$ for a given
				      * @p{fe2}-function @p{z2}, where $I_h$
				      * is the interpolation from @p{fe2}
				      * to @p{fe1}. $(Id-I_h)z2$ is
				      * denoted by @p{z2_difference}.
				      */
    template <int dim, typename number>
    static void interpolation_difference(const DoFHandler<dim> &dof1,
					 const Vector<number> &z1,
					 const FiniteElement<dim> &fe2,
					 Vector<number> &z1_difference);    
    
				     /**
				      * Gives the patchwise
				      * extrapolation of a @p{dof1}
				      * function @p{z1} to a @p{dof2}
				      * function @p{z2}.  @p{dof1} and
				      * @p{dof2} need to be @ref{DoFHandler}
				      * based on the same triangulation.
				      *
				      * This function is interesting
				      * for e.g. extrapolating
				      * patchwise a piecewise linear
				      * solution to a piecewise
				      * quadratic solution.
				      *
				      * Note that the resulting field
				      * does not satisfy continuity
				      * requirements of the given
				      * finite elements. You have to
				      * apply
				      * @ref{ConstraintMatrix}@p{::distribute}
				      * with the hanging node
				      * constraints of the second DoF
				      * handler object to make the
				      * output continuous again.
				      */
    template <int dim, typename number>
    static void extrapolate(const DoFHandler<dim> &dof1,
			    const Vector<number>  &z1,
			    const DoFHandler<dim> &dof2,
			    Vector<number>        &z2);    

				     /**
				      * The numbering of the degrees
				      * of freedom in continous finite
				      * elements is hierarchic,
				      * i.e. in such a way that we
				      * first number the vertex dofs,
				      * in the order of the vertices
				      * as defined by the
				      * triangulation, then the line
				      * dofs in the order and
				      * respecting the direction of
				      * the lines, then the dofs on
				      * quads, etc. However, we could
				      * have, as well, numbered them
				      * in a lexicographic way,
				      * i.e. with indices first
				      * running in x-direction, then
				      * in y-direction and finally in
				      * z-direction. Discontinuous
				      * elements of class @ref{FE_DGQ}
				      * are numbered in this way, for
				      * example.
				      *
				      * This function constructs a
				      * table which lexicographic
				      * index each degree of freedom
				      * in the hierarchic numbering
				      * would have. It operates on the
				      * continuous finite element
				      * given as first argument, and
				      * outputs the lexicographic
				      * indices in the second.
				      *
				      * Note that since this function
				      * uses specifics of the
				      * continuous finite elements, it
				      * can only operate on objects of
				      * type @ref{FE_Q}.
				      *
				      * It is assumed that the size of
				      * the output argument already
				      * matches the correct size,
				      * which is equal to the number
				      * of degrees of freedom in the
				      * finite element.
				      */
    template <int dim>
    static void
    hierarchic_to_lexicographic_numbering (const FE_Q<dim>           &fe,
					   std::vector<unsigned int> &h2l);

				     /**
				      * This is the reverse function
				      * to the above one, generating
				      * the map from the lexicographic
				      * to the hierarchical
				      * numbering. All the remarks
				      * made about the above function
				      * are also valid here.
				      */
    template <int dim>
    static void
    lexicographic_to_hierarchic_numbering (const FE_Q<dim>           &fe,
					   std::vector<unsigned int> &l2h);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidFE);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcTriangulationMismatch);
  
				     /**
				      * Exception
				      */
    DeclException4 (ExcMatrixDimensionMismatch,
		    int, int, int, int,
		    << "This is a " << arg1 << "x" << arg2 << " matrix, "
		    << "but should be a " << arg1 << "x" << arg2 << " matrix.");
};



/*----------------------------   fe_tools.h     ---------------------------*/
/* end of #ifndef __deal2__fe_tools_H */
#endif
/*----------------------------   fe_tools.h     ---------------------------*/
