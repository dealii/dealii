//----------------------------  fe_tools.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools.h  ---------------------------
#ifndef __deal2__fe_tools_H
#define __deal2__fe_tools_H



#include <base/config.h>
#include <base/exceptions.h>

#include <vector>

template <typename number> class FullMatrix;
template <typename number> class Vector;
template <int dim> class FiniteElement;
template <int dim> class DoFHandler;
template <int dim> class hpDoFHandler;
template <int dim> class FE_Q;
class ConstraintMatrix;


#include <base/config.h>
#include <base/exceptions.h>
#include <vector>
#include <string>

/*!@addtogroup febase */
/*@{*/


/**
 * This class performs interpolations and extrapolations of discrete
 * functions of one @p FiniteElement @p fe1 to another @p FiniteElement
 * @p fe2.
 *
 * It also provides the local interpolation matrices that interpolate
 * on each cell. Furthermore it provides the difference matrix
 * $id-I_h$ that is needed for evaluating $(id-I_h)z$ for e.g. the
 * dual solution $z$.
 *
 * @author Ralf Hartmann, 2000; Wolfgang Bangerth, 2003, Guido Kanschat, 2000, 2004
 */
class FETools
{
  public:
				     /**
				      * Compute the vector required to
				      * renumber the dofs of a cell by
				      * component. Furthermore,
				      * compute the vector storing the
				      * start indices of each
				      * component in the local block
				      * vector.
				      *
				      * The second vector is organized
				      * such that there is a vector
				      * for each base element
				      * containing the start index for
				      * each component served by this
				      * base element.
				      *
				      * While the first vector is
				      * checked to have the correct
				      * size, the second one is
				      * reinitialized for convenience.
				      */
    template<int dim>
    static void compute_component_wise(
      const FiniteElement<dim>&                fe,
      std::vector<unsigned int>&               renumbering,
      std::vector<std::vector<unsigned int> >& start_indices);
    
				     /**
				      * Gives the interpolation matrix
				      * that interpolates a @p fe1-
				      * function to a @p fe2-function on
				      * each cell. The interpolation_matrix
				      * needs to be of size
				      * <tt>(fe2.dofs_per_cell, fe1.dofs_per_cell)</tt>.
				      *
				      * Note, that if the finite element
				      * space @p fe1 is a subset of
				      * the finite element space
				      * @p fe2 then the @p interpolation_matrix
				      * is an embedding matrix.
				      */
    template <int dim, typename number>
    static
    void
    get_interpolation_matrix(const FiniteElement<dim> &fe1,
                             const FiniteElement<dim> &fe2,
                             FullMatrix<number> &interpolation_matrix);
    
				     /**
				      * Gives the interpolation matrix
				      * that interpolates a @p fe1-
				      * function to a @p fe2-function, and
				      * interpolates this to a second
				      * @p fe1-function on
				      * each cell. The interpolation_matrix
				      * needs to be of size
				      * <tt>(fe1.dofs_per_cell, fe1.dofs_per_cell)</tt>.
				      *
				      * Note, that this function only
				      * makes sense if the finite element
				      * space due to @p fe1 is not a subset of
				      * the finite element space due to
				      * @p fe2, as if it were a subset then
				      * the @p interpolation_matrix would be 
				      * only the unit matrix.
				      */
    template <int dim, typename number>
    static
    void
    get_back_interpolation_matrix(const FiniteElement<dim> &fe1,
                                  const FiniteElement<dim> &fe2,
                                  FullMatrix<number> &interpolation_matrix);

				     /**
				      * Gives the unit matrix minus the
				      * back interpolation matrix.
				      * The @p difference_matrix
				      * needs to be of size
				      * <tt>(fe1.dofs_per_cell, fe1.dofs_per_cell)</tt>.
				      *
				      * This function gives
				      * the matrix that transforms a
				      * @p fe1 function $z$ to $z-I_hz$
				      * where $I_h$ denotes the interpolation
				      * operator from the @p fe1 space to
				      * the @p fe2 space. This matrix hence
				      * is useful to evaluate
				      * error-representations where $z$
				      * denotes the dual solution.
				      */
    template <int dim, typename number>
    static
    void
    get_interpolation_difference_matrix(const FiniteElement<dim> &fe1,
                                        const FiniteElement<dim> &fe2,
                                        FullMatrix<number> &difference_matrix);

				     /**
				      * Compute the local
				      * $L^2$-projection matrix from
				      * fe1 to fe2.
				      */
    template <int dim, typename number>
    static void get_projection_matrix(const FiniteElement<dim> &fe1,
				      const FiniteElement<dim> &fe2,
				      FullMatrix<number> &matrix);

    
				     /**
				      * Gives the interpolation of a the
				      * @p dof1-function @p u1 to a
				      * @p dof2-function @p u2. @p dof1 and
				      * @p dof2 need to be DoFHandlers
				      * based on the same triangulation.
				      *
				      * If the elements @p fe1 and @p fe2
				      * are either both continuous or
				      * both discontinuous then this
				      * interpolation is the usual point
				      * interpolation. The same is true
				      * if @p fe1 is a continuous and
				      * @p fe2 is a discontinuous finite
				      * element. For the case that @p fe1
				      * is a discontinuous and @p fe2 is
				      * a continuous finite element
				      * there is no point interpolation
				      * defined at the discontinuities.
				      * Therefore the meanvalue is taken
				      * at the DoF values on the
				      * discontinuities.
				      *
				      * Note that for continuous
				      * elements on grids with hanging
				      * nodes (i.e. locally refined
				      * grids) this function does not
				      * give the expected output.
				      * Indeed, the resulting output
				      * vector does not necessarily
				      * respect continuity
				      * requirements at hanging nodes:
				      * if, for example, you are
				      * interpolating a Q2 field to a
				      * Q1 field, then at hanging
				      * nodes the output field will
				      * have the function value of the
				      * input field, which however is
				      * not usually the mean value of
				      * the two adjacent nodes. It is
				      * thus not part of the Q1
				      * function space on the whole
				      * triangulation, although it is
				      * of course Q1 on each cell.
				      *
				      * For this case (continuous
				      * elements on grids with hanging
				      * nodes), please use the
				      * @p interpolate function with
				      * an additional
				      * @p ConstraintMatrix argument,
				      * see below, or make the field
				      * conforming yourself by calling
				      * the @p distribute function of
				      * your hanging node constraints
				      * object.
				      */
    template <int dim,
              template <int> class DH1,
              template <int> class DH2,
              class InVector, class OutVector>
    static
    void
    interpolate (const DH1<dim> &dof1,
                 const InVector &u1,
                 const DH2<dim> &dof2,
                 OutVector      &u2);
    
				     /**
				      * Gives the interpolation of a
				      * the @p dof1-function @p u1
				      * to a @p dof2-function
				      * @p u2. @p dof1 and @p dof2
				      * need to be DoFHandlers
				      * based on the same
				      * triangulation.
				      * @p constraints is a hanging
				      * node constraints object
				      * corresponding to
				      * @p dof2. This object is
				      * particular important when
				      * interpolating onto continuous
				      * elements on grids with hanging
				      * nodes (locally refined grids).
				      *
				      * If the elements @p fe1 and @p fe2
				      * are either both continuous or
				      * both discontinuous then this
				      * interpolation is the usual point
				      * interpolation. The same is true
				      * if @p fe1 is a continuous and
				      * @p fe2 is a discontinuous finite
				      * element. For the case that @p fe1
				      * is a discontinuous and @p fe2 is
				      * a continuous finite element
				      * there is no point interpolation
				      * defined at the discontinuities.
				      * Therefore the meanvalue is taken
				      * at the DoF values on the
				      * discontinuities.
				      */
    template <int dim, class InVector, class OutVector>
    static void interpolate (const DoFHandler<dim>&  dof1,
			     const InVector&         u1,
			     const DoFHandler<dim>&  dof2,
			     const ConstraintMatrix& constraints,
			     OutVector&              u2);


                                     /**
                                      * Same as last function, except
                                      * that one or both of the dof
                                      * handler objects might be of
                                      * type @p hpDoFHandler.
                                      */
    template <int dim,
              template <int> class DH1,
              template <int> class DH2,              
              class InVector, class OutVector>
    static void interpolate (const DH1<dim>         &dof1,
			     const InVector         &u1,
			     const DH2<dim>         &dof2,
			     const ConstraintMatrix &constraints,
			     OutVector&              u2);
    

				     /**
				      * Gives the interpolation of the
				      * @p fe1-function @p u1 to a
				      * @p fe2-function, and
				      * interpolates this to a second
				      * @p fe1-function named
				      * @p u1_interpolated.
				      *
				      * Note, that this function does
				      * not work on continuous
				      * elements at hanging nodes. For
				      * that case use the
				      * @p back_interpolate function,
				      * below, that takes an
				      * additional
				      * @p ConstraintMatrix object.
				      *
				      * Furthermore note, that for the
				      * specific case when the finite
				      * element space corresponding to
				      * @p fe1 is a subset of the
				      * finite element space
				      * corresponding to @p fe2, this
				      * function is simply an identity
				      * mapping.
				      */
    template <int dim, class InVector, class OutVector>
    static void back_interpolate (const DoFHandler<dim>    &dof1,
				  const InVector           &u1,
				  const FiniteElement<dim> &fe2,
				  OutVector                &u1_interpolated);

                                     /**
                                      * Same as last function, except
                                      * that the dof handler objects
                                      * might be of type
                                      * @p hpDoFHandler.
                                      */
    template <int dim,
              template <int> class DH,
              class InVector, class OutVector>
    static void back_interpolate (const DH<dim>            &dof1,
				  const InVector           &u1,
				  const FiniteElement<dim> &fe2,
				  OutVector                &u1_interpolated);

                                     /**
				      * Gives the interpolation of the
				      * @p dof1-function @p u1 to a
				      * @p dof2-function, and
				      * interpolates this to a second
				      * @p dof1-function named
				      * @p u1_interpolated.
				      * @p constraints1 and
				      * @p constraints2 are the
				      * hanging node constraints
				      * corresponding to @p dof1 and
				      * @p dof2, respectively. These
				      * objects are particular
				      * important when continuous
				      * elements on grids with hanging
				      * nodes (locally refined grids)
				      * are involved.
				      *
				      * Furthermore note, that for the
				      * specific case when the finite
				      * element space corresponding to
				      * @p dof1 is a subset of the
				      * finite element space
				      * corresponding to @p dof2, this
				      * function is simply an identity
				      * mapping.
				      */
    template <int dim, class InVector, class OutVector>
    static void back_interpolate (const DoFHandler<dim>&  dof1,
				  const ConstraintMatrix& constraints1,
				  const InVector&         u1,
				  const DoFHandler<dim>&  dof2,
				  const ConstraintMatrix& constraints2,
				  OutVector&              u1_interpolated);

				     /**
				      * Gives $(Id-I_h)z1$ for a given
				      * @p dof1-function @p z1, where $I_h$
				      * is the interpolation from @p fe1
				      * to @p fe2. $(Id-I_h)z1$ is
				      * denoted by @p z1_difference.
				      *
				      * Note, that this function does
				      * not work on continuous
				      * elements at hanging nodes. For
				      * that case use the
				      * @p interpolation_difference
				      * function, below, that takes an
				      * additional
				      * @p ConstraintMatrix object.
				      */
    template <int dim, class InVector, class OutVector>
    static void interpolation_difference(const DoFHandler<dim> &dof1,
					 const InVector &z1,
					 const FiniteElement<dim> &fe2,
					 OutVector &z1_difference);    
    
				     /**
				      * Gives $(Id-I_h)z1$ for a given
				      * @p dof1-function @p z1, where $I_h$
				      * is the interpolation from @p fe1
				      * to @p fe2. $(Id-I_h)z1$ is
				      * denoted by @p z1_difference.
				      * @p constraints1 and
				      * @p constraints2 are the
				      * hanging node constraints
				      * corresponding to @p dof1 and
				      * @p dof2, respectively. These
				      * objects are particular
				      * important when continuous
				      * elements on grids with hanging
				      * nodes (locally refined grids)
				      * are involved.
				      */
    template <int dim, class InVector, class OutVector>
    static void interpolation_difference(const DoFHandler<dim>&  dof1,
					 const ConstraintMatrix& constraints1,
					 const InVector&         z1,
					 const DoFHandler<dim>&  dof2,
					 const ConstraintMatrix& constraints2,
					 OutVector&              z1_difference);
    
				     /**
				      * $L^2$ projection for
				      * discontinuous
				      * elements. Operates the same
				      * direction as interpolate.
				      *
				      * The global projection can be
				      * computed by local matrices if
				      * the finite element spaces are
				      * discontinuous. With continuous
				      * elements, this is impossible,
				      * since a global mass matrix
				      * must be inverted.
				      */
    template <int dim, class InVector, class OutVector>
    static void project_dg (const DoFHandler<dim>& dof1,
			    const InVector&        u1,
			    const DoFHandler<dim>& dof2,
			    OutVector&             u2);
    
				     /**
				      * Gives the patchwise
				      * extrapolation of a @p dof1
				      * function @p z1 to a @p dof2
				      * function @p z2.  @p dof1 and
				      * @p dof2 need to be DoFHandler
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
				      * finite elements.
				      * 
				      * When you use continuous
				      * elements on grids with hanging
				      * nodes, please use the
				      * @p extrapolate function with
				      * an additional
				      * @p ConstraintMatrix argument,
				      * see below.
				      *
				      * Since this function operates
				      * on patches of cells, it is
				      * required that the underlying
				      * grid is refined at least once
				      * for every coarse grid cell. If
				      * this is not the case, an
				      * exception will be raised.
				      */
    template <int dim, class InVector, class OutVector>
    static void extrapolate (const DoFHandler<dim>& dof1,
			     const InVector&        z1,
			     const DoFHandler<dim>& dof2,
			     OutVector&             z2);    

				     /**
				      * Gives the patchwise
				      * extrapolation of a @p dof1
				      * function @p z1 to a @p dof2
				      * function @p z2.  @p dof1 and
				      * @p dof2 need to be DoFHandler
				      * based on the same triangulation.
				      * @p constraints is a hanging
				      * node constraints object
				      * corresponding to
				      * @p dof2. This object is
				      * particular important when
				      * interpolating onto continuous
				      * elements on grids with hanging
				      * nodes (locally refined grids).
				      *
				      * Otherwise, the same holds as
				      * for the other @p extrapolate
				      * function.
				      */
    template <int dim, class InVector, class OutVector>
    static void extrapolate (const DoFHandler<dim>&  dof1,
			     const InVector&         z1,
			     const DoFHandler<dim>&  dof2,
			     const ConstraintMatrix& constraints,
			     OutVector&              z2);    

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
				      * elements of class FE_DGQ()
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
				      * type FE_Q().
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
				      * Given a name in the form which
				      * is returned by the
				      * @p FiniteElement::get_name
				      * function, regenerate such a
				      * finite element.
				      *
				      * This function is useful to
				      * convert the name given in an
				      * input file to an actual finite
				      * element, without having to
				      * parse the name yourself.
				      *
				      * Note that the given name must
				      * match exactly what one would
				      * get from the finite element to
				      * be created, since otherwise
				      * the parsing would fail. If no
				      * finite element can be
				      * reconstructed from this
				      * string, an exception of type
				      * @p FETools::ExcInvalidFEName
				      * is thrown.
				      *
				      * There is one exception,
				      * however, where the names must
				      * not match exactly: while the
				      * finite elements write the
				      * space dimension in the form of
				      * a template argument after the
				      * name of the class (for example
				      * <tt>FE_Q<2></tt>, you can omit the
				      * dimension argument altogether,
				      * or replace it with the string
				      * <tt>@<dim@></tt>. The reason
				      * is that the dimension argument
				      * may be cumbersome if the name
				      * of a finite element is given
				      * in an input file that may be
				      * used to control operation of
				      * the program in different space
				      * dimensions. Running the
				      * program in another space
				      * dimension would then require
				      * changing the input file as
				      * well. With above exception,
				      * there is a canonical spelling
				      * that doesn't require this.
				      *
				      * The function returns a pointer
				      * to a newly create finite
				      * element. It is in the callers
				      * responsibility to destroy the
				      * object pointed to at an
				      * appropriate time.
				      */
    template <int dim>
    static
    FiniteElement<dim> *
    get_fe_from_name (const std::string &name);

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidFEName,
		    std::string,
		    << "Can't re-generate a finite element from the string <"
		    << arg1 << ">.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidFE);

                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcFEMustBePrimitive);
				     /**
				      * Exception
				      */
    DeclException0 (ExcTriangulationMismatch);

				     /**
				      * Exception
				      */
    DeclException1 (ExcHangingNodesNotAllowed,
		    int,
		    << "You are using continuous elements on a grid with "
		    << "hanging nodes but without providing hanging node "
		    << "constraints. Use the respective function with "
		    << "additional ConstraintMatrix argument(s), instead.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcGridNotRefinedAtLeastOnce);
				     /**
				      * Exception
				      */
    DeclException4 (ExcMatrixDimensionMismatch,
		    int, int, int, int,
		    << "This is a " << arg1 << "x" << arg2 << " matrix, "
		    << "but should be a " << arg3 << "x" << arg4 << " matrix.");

  private:
				     /**
				      * Return a finite element that
				      * is created using the
				      * characters of the input
				      * parameters. The second part of
				      * the return value indicates how
				      * many characters have been used
				      * up in the creation of the
				      * finite element, so that the
				      * calling site can continue
				      * parsing finite element lists
				      * (for example for
				      * FESystem objects) at the
				      * position after which the
				      * present element's name ends.
				      *
				      * If no finite element could be
				      * created from the string at the
				      * beginning of the given string,
				      * then an exception is thrown,
				      * just as for the
				      * get_fe_from_name()
				      * function.
				      */
    template <int dim>
    static
    std::pair<FiniteElement<dim> *, unsigned int>
    get_fe_from_name_aux (const std::string &name);    
};

/*@}*/

/*----------------------------   fe_tools.h     ---------------------------*/
/* end of #ifndef __deal2__fe_tools_H */
#endif
/*----------------------------   fe_tools.h     ---------------------------*/
