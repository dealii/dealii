//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_raviart_thomas_h
#define __deal2__fe_raviart_thomas_h

#include <base/config.h>
#include <base/table.h>
#include <base/polynomials_raviart_thomas.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <base/geometry_info.h>
#include <fe/fe.h>
#include <fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class MappingQ;


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of Raviart-Thomas (RT) elements, conforming with the
 * space H<sup>div</sup>. These elements generate vector fields with
 * normel components continuous between mesh cells.
 *
 * We follow the usual definition of the degree of RT elements, which
 * denotes the polynomial degree of the largest complete polynomial
 * subspace contained in the RT space. Then, approciamtion order of
 * the function itself is <i>degree+1</i>, as with usual polynomial
 * spaces.
 *
 * This class is not implemented for the codimension one case
 * (<tt>spacedim != dim</tt>).
 *
 * @todo Even if this element is implemented for two and three space
 * dimensions, the definition of the node values relies on
 * consistently oriented faces in 3D. Therefore, care should be taken
 * on complicated meshes.
 *
 * <h3>Interpolation</h3>
 *
 * The @ref GlossInterpolation "interpolation" operators associated
 * with the RT element are constructed such that interpolation and
 * computing the divergence are commuting operations. We require this
 * from interpolating arbitrary functions as well as the #restriction
 * matrices.  It can be achieved by two interpolation schemes, the
 * simplified one in FE_RaviartThomasNodal and the original one here:
 *
 * <h4>Node values on edges/faces</h4>
 *
 * On edges or faces, the @ref GlossNodes "node values" are the moments of
 * the normal component of the interpolated function with respect to
 * the traces of the RT polynomials. Since the normal trace of the RT
 * space of degree <i>k</i> on an edge/face is the space
 * <i>Q<sub>k</sub></i>, the moments are taken with respect to this
 * space.
 *
 * <h4>Interior node values</h4>
 *
 * Higher order RT spaces have interior nodes. These are moments taken
 * with respect to the gradient of functions in <i>Q<sub>k</sub></i>
 * on the cell (this space is the matching space for RT<sub>k</sub> in
 * a mixed formulation).
 *
 * <h4>Generalized support points</h4>
 *
 * The node values above rely on integrals, which will be computed by
 * quadrature rules themselves. The generalized support points are a
 * set of points such that this quadrature can be performed with
 * sufficient accuracy. The points needed are thode of
 * QGauss<sub>k+1</sub> on each face as well as QGauss<sub>k</sub> in
 * the interior of the cell (or none for RT<sub>0</sub>).
 *
 *
 * @author Guido Kanschat, 2005, based on previous Work by Wolfgang Bangerth
 */
template <int dim>
class FE_RaviartThomas
  :
  public FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim>
{
  public:
				     /**
				      * Constructor for the Raviart-Thomas
				      * element of degree @p p.
				      */
    FE_RaviartThomas (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_RaviartThomas<dim>(degree)</tt>, with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;


				     /**
				      * Number of base elements in a
				      * mixed discretization. Here,
				      * this is of course equal to
				      * one.
				      */
    virtual unsigned int n_base_elements () const;
    
				     /**
				      * Access to base element
				      * objects. Since this element is
				      * atomic, <tt>base_element(0)</tt> is
				      * @p this, and all other
				      * indices throw an error.
				      */
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;

                                     /**
                                      * Multiplicity of base element
                                      * @p index. Since this is an
                                      * atomic element,
                                      * <tt>element_multiplicity(0)</tt>
                                      * returns one, and all other
                                      * indices will throw an error.
                                      */
    virtual unsigned int element_multiplicity (const unsigned int index) const;
    
				     /**
				      * Check whether a shape function
				      * may be non-zero on a face.
				      *
				      * Right now, this is only
				      * implemented for RT0 in
				      * 1D. Otherwise, returns always
				      * @p true.
				      */
    virtual bool has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const;
    
    virtual void interpolate(std::vector<double>&                local_dofs,
			     const std::vector<double>& values) const;
    virtual void interpolate(std::vector<double>&                local_dofs,
			     const std::vector<Vector<double> >& values,
			     unsigned int offset = 0) const;
    virtual void interpolate(
      std::vector<double>& local_dofs,
      const VectorSlice<const std::vector<std::vector<double> > >& values) const;
    virtual unsigned int memory_consumption () const;
    virtual FiniteElement<dim> * clone() const;
    
  private:
				     /**
				      * The order of the
				      * Raviart-Thomas element. The
				      * lowest order elements are
				      * usually referred to as RT0,
				      * even though their shape
				      * functions are piecewise
				      * linears.
				      */  
    const unsigned int rt_order;

    				     /**
				      * Only for internal use. Its
				      * full name is
				      * @p get_dofs_per_object_vector
				      * function and it creates the
				      * @p dofs_per_object vector that is
				      * needed within the constructor to
				      * be passed to the constructor of
				      * @p FiniteElementData.
				      */
    static std::vector<unsigned int>
    get_dpo_vector (const unsigned int degree);

				     /**
				      * Initialize the @p
				      * generalized_support_points
				      * field of the FiniteElement
				      * class and fill the tables with
				      * interpolation weights
				      * (#boundary_weights and
				      * #interior_weights). Called
				      * from the constructor.
				      */
    void initialize_support_points (const unsigned int rt_degree);

				     /**
				      * Initialize the interpolation
				      * from functions on refined mesh
				      * cells onto the father
				      * cell. According to the
				      * philosophy of the
				      * Raviart-Thomas element, this
				      * restriction operator preserves
				      * the divergence of a function
				      * weakly.
				      */
    void initialize_restriction ();
    
				     /**
				      * Given a set of flags indicating
				      * what quantities are requested
				      * from a @p FEValues object,
				      * return which of these can be
				      * precomputed once and for
				      * all. Often, the values of
				      * shape function at quadrature
				      * points can be precomputed, for
				      * example, in which case the
				      * return value of this function
				      * would be the logical and of
				      * the input @p flags and
				      * @p update_values.
				      *
				      * For the present kind of finite
				      * element, this is exactly the
				      * case.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * This is the opposite to the
				      * above function: given a set of
				      * flags indicating what we want
				      * to know, return which of these
				      * need to be computed each time
				      * we visit a new cell.
				      *
				      * If for the computation of one
				      * quantity something else is
				      * also required (for example, we
				      * often need the covariant
				      * transformation when gradients
				      * need to be computed), include
				      * this in the result as well.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;
    
				     /**
				      * Fields of cell-independent data.
				      *
				      * For information about the
				      * general purpose of this class,
				      * see the documentation of the
				      * base class.
				      */
    class InternalData : public FiniteElement<dim>::InternalDataBase
    {
      public:
					 /**
					  * Array with shape function
					  * values in quadrature
					  * points. There is one row
					  * for each shape function,
					  * containing values for each
					  * quadrature point. Since
					  * the shape functions are
					  * vector-valued (with as
					  * many components as there
					  * are space dimensions), the
					  * value is a tensor.
					  *
					  * In this array, we store
					  * the values of the shape
					  * function in the quadrature
					  * points on the unit
					  * cell. The transformation
					  * to the real space cell is
					  * then simply done by
					  * multiplication with the
					  * Jacobian of the mapping.
					  */
	std::vector<std::vector<Tensor<1,dim> > > shape_values;

					 /**
					  * Array with shape function
					  * gradients in quadrature
					  * points. There is one
					  * row for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  *
					  * We store the gradients in
					  * the quadrature points on
					  * the unit cell. We then
					  * only have to apply the
					  * transformation (which is a
					  * matrix-vector
					  * multiplication) when
					  * visiting an actual cell.
					  */
	std::vector<std::vector<Tensor<2,dim> > > shape_gradients;
    };

				     /**
				      * These are the factors
				      * multiplied to a function in
				      * the
				      * #generalized_face_support_points
				      * when computing the
				      * integration. They are
				      * organized such that there is
				      * one row for each generalized
				      * face support point and one
				      * column for each degree of
				      * freedom on the face.
				      */
    Table<2, double> boundary_weights;
				     /**
				      * Precomputed factors for
				      * interpolation of interior
				      * degrees of freedom. The
				      * rationale for this Table is
				      * the same as for
				      * #boundary_weights. Only, this
				      * table has a third coordinate
				      * for the space direction of the
				      * component evaluated.
				      */
    Table<3, double> interior_weights;
    
				     /**
				      * Allow access from other
				      * dimensions.
				      */
    template <int dim1> friend class FE_RaviartThomas;
};



/**
 * The Raviart-Thomas elements with node functionals defined as point
 * values in Gauss points.
 *
 * <h3>Description of node values</h3>
 *
 * For this Raviart-Thomas element, the node values are not cell and
 * face moments with respect to certain polynomials, but the values in
 * quadrature points. Following the general scheme for numbering
 * degrees of freedom, the node values on edges are first, edge by
 * edge, according to the natural ordering of the edges of a cell. The
 * interior degrees of freedom are last.
 *
 * For an RT-element of degree <i>k</i>, we choose
 * <i>(k+1)<sup>d-1</sup></i> Gauss points on each face. These points
 * are ordered lexicographically with respect to the orientation of
 * the face. This way, the normal component which is in
 * <i>Q<sub>k</sub></i> is uniquely determined. Furthermore, since
 * this Gauss-formula is exact on <i>Q<sub>2k+1</sub></i>, these node
 * values correspond to the exact integration of the moments of the
 * RT-space.
 *
 * In the interior of the cells, the moments are with respect to an
 * anisotropic <i>Q<sub>k</sub></i> space, where the test functions
 * are one degree lower in the direction corresponding to the vector
 * component under consideration. This is emulated by using an
 * anisotropic Gauss formula for integration.
 *
 * @todo The current implementation is for Cartesian meshes
 * only. You must use MappingCartesian.
 * 
 * @todo Even if this element is implemented for two and three space
 * dimensions, the definition of the node values relies on
 * consistently oriented faces in 3D. Therefore, care should be taken
 * on complicated meshes.
 *
 * @note The degree stored in the member variable
 * FiniteElementData<dim>::degree is higher by one than the
 * constructor argument!
 * 
 * @author Guido Kanschat, 2005, Zhu Liang, 2008
 */
template <int dim>
class FE_RaviartThomasNodal
  :
  public FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim>
{
  public:
				     /**
				      * Constructor for the Raviart-Thomas
				      * element of degree @p p.
				      */
    FE_RaviartThomasNodal (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_RaviartThomasNodal<dim>(degree)</tt>, with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;

    virtual FiniteElement<dim>* clone () const;

    virtual void interpolate(std::vector<double>&                local_dofs,
			     const std::vector<double>& values) const;
    virtual void interpolate(std::vector<double>&                local_dofs,
			     const std::vector<Vector<double> >& values,
			     unsigned int offset = 0) const;    
    virtual void interpolate(
      std::vector<double>& local_dofs,
      const VectorSlice<const std::vector<std::vector<double> > >& values) const;

       
    virtual void get_face_interpolation_matrix (const FiniteElement<dim> &source,
						FullMatrix<double>       &matrix) const;

    virtual void get_subface_interpolation_matrix (const FiniteElement<dim> &source,
						   const unsigned int        subface,
						   FullMatrix<double>       &matrix) const;
    virtual bool hp_constraints_are_implemented () const;
    
    virtual std::vector<std::pair<unsigned int, unsigned int> >
    hp_vertex_dof_identities (const FiniteElement<dim> &fe_other) const;

    virtual std::vector<std::pair<unsigned int, unsigned int> >
    hp_line_dof_identities (const FiniteElement<dim> &fe_other) const;

    virtual std::vector<std::pair<unsigned int, unsigned int> >
    hp_quad_dof_identities (const FiniteElement<dim> &fe_other) const;
    
    virtual FiniteElementDomination::Domination
    compare_for_face_domination (const FiniteElement<dim> &fe_other) const;

    virtual UpdateFlags update_once (const UpdateFlags flags) const;
    virtual UpdateFlags update_each (const UpdateFlags flags) const;
  private:
				     /**
				      * Only for internal use. Its
				      * full name is
				      * @p get_dofs_per_object_vector
				      * function and it creates the
				      * @p dofs_per_object vector that is
				      * needed within the constructor to
				      * be passed to the constructor of
				      * @p FiniteElementData.
				      */
    static std::vector<unsigned int>
    get_dpo_vector (const unsigned int degree);

				     /**
				      * Compute the vector used for
				      * the
				      * @p restriction_is_additive
				      * field passed to the base
				      * class's constructor.
				      */
    static std::vector<bool>
    get_ria_vector (const unsigned int degree);
				     /**
				      * Initialize the
				      * FiniteElement<dim>::generalized_support_points
				      * and FiniteElement<dim>::generalized_face_support_points
				      * fields. Called from the
				      * constructor.
				      */
    void initialize_support_points (const unsigned int rt_degree);
};


/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
std::vector<unsigned int> FE_RaviartThomas<1>::get_dpo_vector (const unsigned int);
template <>
void
FE_RaviartThomas<1>::initialize_restriction();

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
