//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_tools_H
#define __deal2__fe_tools_H



#include <base/config.h>
#include <base/exceptions.h>
#include <base/geometry_info.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

template <typename number> class FullMatrix;
template <typename number> class Vector;
template <int dim> class Quadrature;
template <int dim, int spacedim> class FiniteElement;
template <int dim, int spacedim> class DoFHandler;
namespace hp
{
  template <int dim, int spacedim> class DoFHandler;
}
template <int dim> class FiniteElementData;
class ConstraintMatrix;


#include <base/config.h>
#include <base/exceptions.h>
#include <vector>
#include <string>

/*!@addtogroup feall */
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
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, Guido Kanschat;
 * 2000, 2003, 2004, 2005, 2006
 */
class FETools
{
  public:
				     /**
				      * A base class for factory
				      * objects creating finite
				      * elements of a given
				      * degree. Derived classes are
				      * called whenever one wants to
				      * have a transparent way to
				      * create a finite element
				      * object.
				      *
				      * This class is used in the
				      * FETools::get_fe_from_name()
				      * and FETools::add_fe_name()
				      * functions.
				      *
				      * @author Guido Kanschat, 2006
				      */
    template <int dim, int spacedim=dim>
    class FEFactoryBase
    {
      public:
					 /**
					  * Create a FiniteElement and
					  * return a pointer to it.
					  */
	virtual FiniteElement<dim,spacedim>*
	get (const unsigned int degree) const = 0;

                                         /**
					  * Virtual destructor doing
					  * nothing but making the
					  * compiler happy.
					  */
	virtual ~FEFactoryBase();
    };
    
				     /**
				      * A concrete class for factory
				      * objects creating finite
				      * elements of a given degree.
				      *
				      * The class's get() function
				      * generates a finite element
				      * object of the type given as
				      * template argument, and with
				      * the degree (however the finite
				      * element class wishes to
				      * interpret this number) given
				      * as argument to get().
				      * 
				      * @author Guido Kanschat, 2006
				      */
    template <class FE>
    class FEFactory : public FEFactoryBase<FE::dimension,FE::dimension>
    {
      public:
					 /**
					  * Create a FiniteElement and
					  * return a pointer to it.
					  */
	virtual FiniteElement<FE::dimension,FE::dimension>*
	get (const unsigned int degree) const;
    };
    
				     /**
				      * @warning In most cases, you
				      * will probably want to use
				      * compute_base_renumbering().
				      *
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
    template<int dim, int spacedim>
    static void compute_component_wise(
      const FiniteElement<dim,spacedim>&                fe,
      std::vector<unsigned int>&               renumbering,
      std::vector<std::vector<unsigned int> >& start_indices);

				     /**
				      * Compute the vector required to
				      * renumber the dofs of a cell by
				      * block. Furthermore,
				      * compute the vector storing the
				      * start indices of each local block
				      * vector.
				      *
				      * @todo Which way does this
				      * vector map the numbers?
				      */
    template<int dim, int spacedim>
    static void compute_block_renumbering (
      const FiniteElement<dim,spacedim>&  fe,
      std::vector<unsigned int>& renumbering,
      std::vector<unsigned int>& start_indices);
    
				     /**
				      * @name Generation of local matrices
				      * @{
				      */
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
    template <int dim, typename number, int spacedim>
    static
    void
    get_interpolation_matrix(const FiniteElement<dim,spacedim> &fe1,
                             const FiniteElement<dim,spacedim> &fe2,
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
    template <int dim, typename number, int spacedim>
    static
    void
    get_back_interpolation_matrix(const FiniteElement<dim,spacedim> &fe1,
                                  const FiniteElement<dim,spacedim> &fe2,
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
    template <int dim, typename number, int spacedim>
    static
    void
    get_interpolation_difference_matrix(const FiniteElement<dim,spacedim> &fe1,
                                        const FiniteElement<dim,spacedim> &fe2,
                                        FullMatrix<number> &difference_matrix);

				     /**
				      * Compute the local
				      * $L^2$-projection matrix from
				      * fe1 to fe2.
				      */
    template <int dim, typename number, int spacedim>
    static void get_projection_matrix(const FiniteElement<dim,spacedim> &fe1,
				      const FiniteElement<dim,spacedim> &fe2,
				      FullMatrix<number> &matrix);
    
				     /**
				      * Compute the matrix of nodal
				      * values of a finite element
				      * applied to all its shape
				      * functions.
				      *
				      * This function is supposed to
				      * help building finite elements
				      * from polynomial spaces and
				      * should be called inside the
				      * constructor of an
				      * element. Applied to a
				      * completely initialized finite
				      * element, the result should be
				      * the unit matrix by definition
				      * of the node values.
				      *
				      * Using this matrix allows the
				      * construction of the basis of
				      * shape functions in two steps.
				      * <ol>
				      *
				      * <li>Define the space of shape
				      * functions using an arbitrary
				      * basis <i>w<sub>j</sub></i> and
				      * compute the matrix <i>M</i> of
				      * node functionals
				      * <i>N<sub>i</sub></i> applied
				      * to these basis functions.
				      *
				      * <li>Compute the basis
				      * <i>v<sub>j</sub></i> of the
				      * finite element shape function
				      * space by applying
				      * <i>M<sup>-1</sup></i> to the
				      * basis <i>w<sub>j</sub></i>.
				      * </ol>
				      *
				      * @note The FiniteElement must
				      * provide generalized support
				      * points and and interpolation
				      * functions.
				      */
    template <int dim, int spacedim>
    static void compute_node_matrix(FullMatrix<double>& M,
				    const FiniteElement<dim,spacedim>& fe);
    
				     /**
				      * For all possible (isotropic
				      * and anisotropic) refinement
				      * cases compute the embedding
				      * matrices from a coarse cell to
				      * the child cells. Each column
				      * of the resulting matrices
				      * contains the representation of
				      * a coarse grid basis functon by
				      * the fine grid basis; the
				      * matrices are split such that
				      * there is one matrix for every
				      * child.
				      *
				      * This function computes the
				      * coarse grid function in a
				      * sufficiently large number of
				      * quadrature points and fits the
				      * fine grid functions using
				      * least squares
				      * approximation. Therefore, the
				      * use of this function is
				      * restricted to the case that
				      * the finite element spaces are
				      * actually nested.
				      *
				      * Note, that
				      * <code>matrices[refinement_case-1][child]</code>
				      * includes the embedding (or prolongation)
				      * matrix of child
				      * <code>child</code> for the
				      * RefinementCase
				      * <code>refinement_case</code>. Here,
				      * we use
				      * <code>refinement_case-1</code>
				      * instead of
				      * <code>refinement_case</code>
				      * as for
				      * RefinementCase::no_refinement(=0)
				      * there are no prolongation
				      * matrices available.
				      *
				      * Typically this function is
				      * called by the various
				      * implementations of
				      * FiniteElement classes in order
				      * to fill the respective
				      * FiniteElement::prolongation
				      * matrices.
				      *
				      * @param <tt>fe</tt> The finite element
				      * class for which we compute the
				      * embedding matrices.
				      *
				      * @param
				      * <tt>matrices</tt> A reference to
				      * <tt>RefinementCase<dim>::isotropic_refinement</tt>
				      * vectors of FullMatrix
				      * objects. Each vector
				      * corresponds to one
				      * RefinementCase @p
				      * refinement_case and is of the
				      * vector size
				      * <tt>GeometryInfo<dim>::n_children(refinement_case)</tt>. This
				      * is the format used in
				      * FiniteElement, where we want
				      * to use this function mostly.
				      *
				      * @param
				      * <tt>isotropic_only</tt>: set
				      * this <tt>true</tt> if you only
				      * want to compute matrices for
				      * isotropic refinement.
				      */
    template <int dim, typename number, int spacedim>
    static void compute_embedding_matrices(const FiniteElement<dim,spacedim> &fe,
					   std::vector<std::vector<FullMatrix<number> > >& matrices,
					   const bool isotropic_only = false);

				     /**
				      * Compute the embedding matrices
				      * on faces needed for constraint
				      * matrices.
				      *
				      * @param fe The finite element
				      * for which to compute these
				      * matrices.  @param matrices An
				      * array of
				      * <i>GeometryInfo<dim>::subfaces_per_face
				      * = 2<sup>dim-1</sup></i>
				      * FullMatrix objects,holding the
				      * embedding matrix for each
				      * subface.  @param face_coarse
				      * The number of the face on the
				      * coarse side of the face for
				      * which this is computed.
				      * @param face_fine The number of
				      * the face on the refined side
				      * of the face for which this is
				      * computed.
				      *
				      * @warning This function will be
				      * used in computing constraint
				      * matrices. It is not
				      * sufficiently tested yet.
				     */
    template <int dim, typename number, int spacedim>
    static void
    compute_face_embedding_matrices(const FiniteElement<dim,spacedim>& fe,
				    FullMatrix<number> (&matrices)[GeometryInfo<dim>::max_children_per_face],
				    const unsigned int face_coarse,
				    const unsigned int face_fine);

				     /**
				      * For all possible (isotropic
				      * and anisotropic) refinement
				      * cases compute the
				      * <i>L<sup>2</sup></i>-projection
				      * matrices from the children to
				      * a coarse cell.
				      *
				      * Note, that
				      * <code>matrices[refinement_case-1][child]</code>
				      * includes the projection (or restriction)
				      * matrix of child
				      * <code>child</code> for the
				      * RefinementCase
				      * <code>refinement_case</code>. Here,
				      * we use
				      * <code>refinement_case-1</code>
				      * instead of
				      * <code>refinement_case</code>
				      * as for
				      * RefinementCase::no_refinement(=0)
				      * there are no projection
				      * matrices available.
				      *
				      * Typically this function is
				      * called by the various
				      * implementations of
				      * FiniteElement classes in order
				      * to fill the respective
				      * FiniteElement::restriction
				      * matrices.
				      *
				      * @arg fe The finite element
				      * class for which we compute the
				      * projection matrices.  @arg
				      * matrices A reference to
				      * <tt>RefinementCase<dim>::isotropic_refinement</tt>
				      * vectors of FullMatrix
				      * objects. Each vector
				      * corresponds to one
				      * RefinementCase @p
				      * refinement_case and is of the
				      * vector size
				      * <tt>GeometryInfo<dim>::n_children(refinement_case)</tt>. This
				      * is the format used in
				      * FiniteElement, where we want
				      * to use this function mostly.
				      */
    template <int dim, typename number, int spacedim>
    static void compute_projection_matrices(const FiniteElement<dim,spacedim> &fe,
					    std::vector<std::vector<FullMatrix<number> > >& matrices);

//TODO:[WB] Replace this documentation by something comprehensible
    
				     /**
                                      * Projects scalar data defined in
                                      * quadrature points to a finite element
                                      * space on a single cell.
                                      *
                                      * What this function does is the
                                      * following: assume that there is scalar
                                      * data <tt>u<sub>q</sub>, 0 <= q <
                                      * Q:=quadrature.size()</tt>
                                      * defined at the quadrature points of a
                                      * cell, with the points defined by the
                                      * given <tt>rhs_quadrature</tt>
                                      * object. We may then want to ask for
                                      * that finite element function (on a
                                      * single cell) <tt>v<sub>h</sub></tt> in
                                      * the finite-dimensional space defined
                                      * by the given FE object that is the
                                      * projection of <tt>u</tt> in the
                                      * following sense:
                                      *
                                      * Usually, the projection
                                      * <tt>v<sub>h</sub></tt> is that
                                      * function that satisfies
                                      * <tt>(v<sub>h</sub>,w)=(u,w)</tt> for
                                      * all discrete test functions
                                      * <tt>w</tt>. In the present case, we
                                      * can't evaluate the right hand side,
                                      * since <tt>u</tt> is only defined in
                                      * the quadrature points given by
                                      * <tt>rhs_quadrature</tt>, so we replace
                                      * it by a quadrature
                                      * approximation. Likewise, the left hand
                                      * side is approximated using the
                                      * <tt>lhs_quadrature</tt> object; if
                                      * this quadrature object is chosen
                                      * appropriately, then the integration of
                                      * the left hand side can be done
                                      * exactly, without any
                                      * approximation. The use of different
                                      * quadrature objects is necessary if the
                                      * quadrature object for the right hand
                                      * side has too few quadrature points --
                                      * for example, if data <tt>q</tt> is
                                      * only defined at the cell center, then
                                      * the corresponding one-point quadrature
                                      * formula is obviously insufficient to
                                      * approximate the scalar product on the
                                      * left hand side by a definite form.
                                      *
                                      * After these quadrature approximations,
                                      * we end up with a nodal representation
                                      * <tt>V<sub>h</sub></tt> of
                                      * <tt>v<sub>h</sub></tt> that satisfies
                                      * the following system of linear
                                      * equations: <tt>M V<sub>h</sub> = Q
                                      * U</tt>, where
                                      * <tt>M<sub>ij</sub>=(phi_i,phi_j)</tt>
                                      * is the mass matrix approximated by
                                      * <tt>lhs_quadrature</tt>, and
                                      * <tt>Q</tt> is the matrix
                                      * <tt>Q<sub>iq</sub>=phi<sub>i</sub>(x<sub>q</sub>)
                                      * w<sub>q</sub></tt> where
                                      * <tt>w<sub>q</sub></tt> are quadrature
                                      * weights; <tt>U</tt> is the vector of
                                      * quadrature point data
                                      * <tt>u<sub>q</sub></tt>.
                                      *
                                      * In order to then get the nodal
                                      * representation <tt>V<sub>h</sub></tt>
                                      * of the projection of <tt>U</tt>, one
                                      * computes <tt>V<sub>h</sub> = X U,
                                      * X=M<sup>-1</sup> Q</tt>. The purpose
                                      * of this function is to compute the
                                      * matrix <tt>X</tt> and return it
                                      * through the last argument of this
                                      * function.
                                      *
                                      * Note that this function presently only
                                      * supports scalar data. An extension of
                                      * the mass matrix is of course trivial,
                                      * but one has to define the order of
                                      * data in the vector <tt>U</tt> if it
                                      * contains vector valued data in all
                                      * quadrature points.
                                      *
                                      * A use for this function is described
                                      * in the introduction to the @ref step_18 "step-18"
                                      * example program.
                                      *
                                      * The opposite of this function,
                                      * interpolation of a finite element
                                      * function onto quadrature points is
                                      * essentially what the
                                      * <tt>FEValues::get_function_values</tt>
                                      * functions do; to make things a little
                                      * simpler, the
                                      * <tt>FETools::compute_interpolation_to_quadrature_points_matrix</tt>
                                      * provides the matrix form of this.
				      *
				      * Note that this function works
				      * on a single cell, rather than
				      * an entire triangulation. In
				      * effect, it therefore doesn't
				      * matter if you use a continuous
				      * or discontinuous version of
				      * the finite element.
				      *
				      * It is worth noting that there
				      * are a few confusing cases of
				      * this function. The first one
				      * is that it really only makes
				      * sense to project onto a finite
				      * element that has at most as
				      * many degrees of freedom per
				      * cell as there are quadrature
				      * points; the projection of N
				      * quadrature point data into a
				      * space with M>N unknowns is
				      * well-defined, but often yields
				      * funny and non-intuitive
				      * results. Secondly, one would
				      * think that if the quadrature
				      * point data is defined in the
				      * support points of the finite
				      * element, i.e. the quadrature
				      * points of
				      * <tt>ths_quadrature</tt> equal
				      * <tt>fe.get_unit_support_points()</tt>,
				      * then the projection should be
				      * the identity, i.e. each degree
				      * of freedom of the finite
				      * element equals the value of
				      * the given data in the support
				      * point of the corresponding
				      * shape function. However, this
				      * is not generally the case:
				      * while the matrix <tt>Q</tt> in
				      * that case is the identity
				      * matrix, the mass matrix
				      * <tt>M</tt> is not equal to the
				      * identity matrix, except for
				      * the special case that the
				      * quadrature formula
				      * <tt>lhs_quadrature</tt> also
				      * has its quadrature points in
				      * the support points of the
				      * finite element.
                                      */
    template <int dim, int spacedim>
    static
    void
    compute_projection_from_quadrature_points_matrix (const FiniteElement<dim,spacedim> &fe,
                                                      const Quadrature<dim>    &lhs_quadrature,
                                                      const Quadrature<dim>    &rhs_quadrature,
                                                      FullMatrix<double>       &X);

                                     /**
                                      * Given a (scalar) local finite element
                                      * function, compute the matrix that maps
                                      * the vector of nodal values onto the
                                      * vector of values of this function at
                                      * quadrature points as given by the
                                      * second argument. In a sense, this
                                      * function does the opposite of the @p
                                      * compute_projection_from_quadrature_points_matrix
                                      * function.
                                      */
    template <int dim, int spacedim>
    static
    void
    compute_interpolation_to_quadrature_points_matrix (const FiniteElement<dim,spacedim> &fe,
                                                       const Quadrature<dim>    &quadrature,
                                                       FullMatrix<double>       &I_q);
    

				     //@}
				     /**
				      * @name Functions which should be in DoFTools
				      */
				     //@{
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
    template <int dim, int spacedim,
              template <int,int> class DH1,
              template <int,int> class DH2,
              class InVector, class OutVector>
    static
    void
    interpolate (const DH1<dim,spacedim> &dof1,
                 const InVector          &u1,
                 const DH2<dim,spacedim> &dof2,
                 OutVector               &u2);
    
				     /**
				      * Gives the interpolation of a
				      * the @p dof1-function @p u1 to
				      * a @p dof2-function @p u2. @p
				      * dof1 and @p dof2 need to be
				      * DoFHandlers (or
				      * hp::DoFHandlers) based on the
				      * same triangulation.  @p
				      * constraints is a hanging node
				      * constraints object
				      * corresponding to @p dof2. This
				      * object is particular important
				      * when interpolating onto
				      * continuous elements on grids
				      * with hanging nodes (locally
				      * refined grids).
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
    template <int dim, int spacedim, 
              template <int, int> class DH1,
              template <int, int> class DH2,              
              class InVector, class OutVector>
    static void interpolate (const DH1<dim,spacedim>  &dof1,
			     const InVector           &u1,
			     const DH2<dim,spacedim>  &dof2,
			     const ConstraintMatrix   &constraints,
			     OutVector&                u2);    

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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void back_interpolate (const DoFHandler<dim,spacedim>    &dof1,
				  const InVector           &u1,
				  const FiniteElement<dim,spacedim> &fe2,
				  OutVector                &u1_interpolated);

                                     /**
                                      * Same as last function, except
                                      * that the dof handler objects
                                      * might be of type
                                      * @p hp::DoFHandler.
                                      */
    template <int dim,
              template <int> class DH,
              class InVector, class OutVector, int spacedim>
    static void back_interpolate (const DH<dim>            &dof1,
				  const InVector           &u1,
				  const FiniteElement<dim,spacedim> &fe2,
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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void back_interpolate (const DoFHandler<dim,spacedim>&  dof1,
				  const ConstraintMatrix& constraints1,
				  const InVector&         u1,
				  const DoFHandler<dim,spacedim>&  dof2,
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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void interpolation_difference(const DoFHandler<dim,spacedim> &dof1,
					 const InVector &z1,
					 const FiniteElement<dim,spacedim> &fe2,
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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void interpolation_difference(const DoFHandler<dim,spacedim>&  dof1,
					 const ConstraintMatrix& constraints1,
					 const InVector&         z1,
					 const DoFHandler<dim,spacedim>&  dof2,
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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void project_dg (const DoFHandler<dim,spacedim>& dof1,
			    const InVector&        u1,
			    const DoFHandler<dim,spacedim>& dof2,
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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void extrapolate (const DoFHandler<dim,spacedim>& dof1,
			     const InVector&        z1,
			     const DoFHandler<dim,spacedim>& dof2,
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
    template <int dim, class InVector, class OutVector, int spacedim>
    static void extrapolate (const DoFHandler<dim,spacedim>&  dof1,
			     const InVector&         z1,
			     const DoFHandler<dim,spacedim>&  dof2,
			     const ConstraintMatrix& constraints,
			     OutVector&              z2);    
				     //@}
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
				      * can only operate on
				      * FiniteElementData<dim> objects
				      * inherent in FE_Q(). However,
				      * this function does not take a
				      * FE_Q object as it is also
				      * invoked by the FE_Q()
				      * constructor.
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
    hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe_data,
					   std::vector<unsigned int>    &h2l);

				     /**
				      * Like the previous function but
				      * instead of returning its
				      * result through the last
				      * argument return it as a value.
				      */
    template <int dim>
    static
    std::vector<unsigned int>
    hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe_data);    
    
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
    lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data,
					   std::vector<unsigned int>    &l2h);

				     /**
				      * Like the previous function but
				      * instead of returning its
				      * result through the last
				      * argument return it as a value.
				      */
    template <int dim>
    static
    std::vector<unsigned int>
    lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data);
    
				     /**
				      * Parse the name of a finite
				      * element and generate a finite
				      * element object accordingly.
				      *
				      * The name must be in the form which
				      * is returned by the
				      * @p FiniteElement::get_name
				      * function, where a few
				      * modifications are allowed:
				      *
				      * <ul><li> Dimension template
				      * parameters &lt;2&gt; etc. can
				      * be omitted. Alternatively, the
				      * explicit number can be
				      * replaced by <tt>dim</tt> or
				      * <tt>d</tt>. If a number is
				      * given, it <b>must</b> match
				      * the template parameter of this
				      * function.
				      *
				      * <li> The powers used for
				      * FESystem may either be numbers
				      * or can be
				      * replaced by <tt>dim</tt> or
				      * <tt>d</tt>.
				      * </ul>
				      *
				      * If no finite element can be
				      * reconstructed from this
				      * string, an exception of type
				      * @p FETools::ExcInvalidFEName
				      * is thrown.
				      *
				      * The function returns a pointer
				      * to a newly create finite
				      * element. It is in the caller's
				      * responsibility to destroy the
				      * object pointed to at an
				      * appropriate later time.
				      *
				      * Since the value of the template
				      * argument can't be deduced from the
				      * (string) argument given to this
				      * function, you have to explicitly
				      * specify it when you call this
				      * function.
				      *
				      * This function knows about all
				      * the standard elements defined
				      * in the library. However, it
				      * doesn't by default know about
				      * elements that you may have
				      * defined in your program. To
				      * make your own elements known
				      * to this function, use the
				      * add_fe_name() function.
				      * This function does not work 
				      * if one wants to get a codimension
				      * 1 finite element.
				      */
    template <int dim>
    static
    FiniteElement<dim, dim> *
    get_fe_from_name (const std::string &name);


				     /**
				      * Extend the list of finite
				      * elements that can be generated
				      * by get_fe_from_name() by the
				      * one given as @p name. If
				      * get_fe_from_name() is later
				      * called with this name, it will
				      * use the object given as second
				      * argument to create a finite
				      * element object.
				      *
				      * The format of the @p name
				      * parameter should include the
				      * name of a finite
				      * element. However, it is safe
				      * to use either the class name
				      * alone or to use the result of
				      * FiniteElement::get_name (which
				      * includes the space dimension
				      * as well as the polynomial
				      * degree), since everything
				      * after the first non-name
				      * character will be ignored.
				      *
				      * The FEFactory object should be
				      * an object newly created with
				      * <tt>new</tt>. FETools will
				      * take ownership of this object
				      * and delete it once it is not
				      * used anymore.
				      *
				      * In most cases, if you want
				      * objects of type
				      * <code>MyFE</code> be created
				      * whenever the name
				      * <code>my_fe</code> is given to
				      * get_fe_from_name, you will
				      * want the second argument to
				      * this function be of type
				      * FEFactory@<MyFE@>, but you can
				      * of course create your custom
				      * finite element factory class.
				      *
				      * This function takes over
				      * ownership of the object given
				      * as second argument, i.e. you
				      * should never attempt to
				      * destroy it later on. The
				      * object will be deleted at the
				      * end of the program's lifetime.
				      * 
				      * If the name of the element
				      * is already in use, an exception
				      * is thrown. Thus, functionality
				      * of get_fe_from_name() can only
				      * be added, not changed.
				      *
				      * @note This function
				      * manipulates a global table
				      * (one table for each space
				      * dimension). It is thread safe
				      * in the sense that every access
				      * to this table is secured by a
				      * lock. Nevertheless, since each
				      * name can be added only once,
				      * user code has to make sure
				      * that only one thread adds a
				      * new element.
				      *
				      * Note also that this table
				      * exists once for each space
				      * dimension. If you have a
				      * program that works with finite
				      * elements in different space
				      * dimensions (for example, @ref
				      * step_4 "step-4" does something
				      * like this), then you should
				      * call this function for each
				      * space dimension for which you
				      * want your finite element added
				      * to the map.
				      */
    template <int dim, int spacedim>
    static void add_fe_name (const std::string& name,
			     const FEFactoryBase<dim,spacedim>* factory);
    
				     /**
				      * The string used for
				      * get_fe_from_name() cannot be
				      * translated to a finite
				      * element.
				      *
				      * Either the string is badly
				      * formatted or you are using a
				      * custom element that must be
				      * added using add_fe_name()
				      * first.
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcInvalidFEName,
		    std::string,
		    << "Can't re-generate a finite element from the string '"
		    << arg1 << "'.");
    
				     /**
				      * The string used for
				      * get_fe_from_name() cannot be
				      * translated to a finite
				      * element.
				      *
				      * Dimension arguments in finite
				      * element names should be
				      * avoided. If they are there,
				      * the dimension should be
				      * <tt>dim</tt> or
				      * <tt>d</tt>. Here, you gave a
				      * numeric dimension argument,
				      * which does not match the
				      * template dimension of the
				      * finite element class.
				      *
				      * @ingroup Exceptions
				      */
    DeclException2 (ExcInvalidFEDimension,
		    char, int,
		    << "The dimension " << arg1
		    << " in the finite element string must match "
		    << "the space dimension "
		    << arg2 << ".");
    
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcInvalidFE);

                                     /**
                                      * The finite element must be
				      * @ref GlossPrimitive "primitive".
				      *
				      * @ingroup Exceptions
                                      */
    DeclException0 (ExcFENotPrimitive);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcTriangulationMismatch);

				     /**
				      * A continuous element is used
				      * on a mesh with hanging nodes,
				      * but the constraint matrices
				      * are missing.
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcHangingNodesNotAllowed,
		    int,
		    << "You are using continuous elements on a grid with "
		    << "hanging nodes but without providing hanging node "
		    << "constraints. Use the respective function with "
		    << "additional ConstraintMatrix argument(s), instead.");
				     /**
				      * You need at least two grid levels.
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcGridNotRefinedAtLeastOnce);
				     /**
				      * The dimensions of the matrix
				      * used did not match the
				      * expected dimensions.
				      *
				      * @ingroup Exceptions
				      */
    DeclException4 (ExcMatrixDimensionMismatch,
		    int, int, int, int,
		    << "This is a " << arg1 << "x" << arg2 << " matrix, "
		    << "but should be a " << arg3 << "x" << arg4 << " matrix.");

				     /**
				      * Exception thrown if an
				      * embedding matrix was computed
				      * inaccurately.
				      *
				      * @ingroup Exceptions
				      */
    DeclException1(ExcLeastSquaresError, double,
		   << "Least squares fit leaves a gap of " << arg1);
  private:
				     /**
				      * Return a finite element that
				      * is created using the beginning
				      * of <tt>name</tt> and eat away
				      * the part of <tt>name</tt>
				      * defining this element.
				      */
    template <int dim, int spacedim>
    static
    FiniteElement<dim,spacedim> *
    get_fe_from_name_aux (std::string &name);
    
};


template<class FE>
FiniteElement<FE::dimension, FE::dimension>*
FETools::FEFactory<FE>::get (const unsigned int degree) const
{
  return new FE(degree);
}



/*@}*/

DEAL_II_NAMESPACE_CLOSE

/*----------------------------   fe_tools.h     ---------------------------*/
/* end of #ifndef __deal2__fe_tools_H */
#endif
/*----------------------------   fe_tools.h     ---------------------------*/
