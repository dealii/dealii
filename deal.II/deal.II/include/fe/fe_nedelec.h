//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_nedelec_h
#define __deal2__fe_nedelec_h

#include <base/config.h>
#include <base/table.h>
#include <base/tensor.h>
#include <base/tensor_base.h>
#include <base/polynomials_nedelec.h>
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
 * Implementation of Nédélec elements, conforming with the
 * space H<sup>curl</sup>. These elements generate vector fields with
 * tangential components continuous between mesh cells.
 *
 * We follow the usual definition of the degree of Nédélec elements,
 * which denotes the polynomial degree of the lowest complete polynomial
 * subspace contained in the Nédélec space. Then, approximation order of
 * the function itself is <i>degree</i>. In this scheme, the lowest
 * order element would be created by the call FE_Nedelec<dim>(0).
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
 * with the Nédélec element are constructed such that interpolation and
 * computing the curl are commuting operations. We require this
 * from interpolating arbitrary functions as well as the #restriction
 * matrices.
 *
 * <h4>Node values</h4>
 *
 * The @ref GlossNodes "node values" on edges are the moments of the
 * tangential component of the interpolated function with respect to
 * the traces of the Nédélec polynomials. Higher-order Nédélec spaces
 * also have face and interior nodes.
 *
 * <h4>Generalized support points</h4>
 *
 * The node values above rely on integrals, which will be computed by
 * quadrature rules themselves. The generalized support points are a
 * set of points such that this quadrature can be performed with
 * sufficient accuracy. The points needed are thode of
 * QGauss<sub>k+1</sub> on each edge and QGauss<sub>k+2</sub> on each face and in
 * the interior of the cell (or none for N<sub>1</sub>).
 *
 *
 * @author Markus Bürg, 2009
 */
template <int dim>
class FE_Nedelec : public FE_PolyTensor<PolynomialsNedelec<dim>, dim> {
   public:
				     /**
				      * Constructor for the Nédélec
				      * element of degree @p p.
				      */
      FE_Nedelec (const unsigned int p);

				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_Nedelec<dim>(degree)</tt>, with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;


				     /**
				      * Check whether a shape function
				      * may be non-zero on a face.
				      */
    virtual bool has_support_on_face (const unsigned int shape_index, const unsigned int face_index) const;

				     /**
				      * Return whether this element implements its
				      * hanging node constraints in the new way, which
				      * has to be used to make elements "hp compatible".
				      *
				      * For the <tt>FE_Nedelec</tt> class the result is
				      * always true (independent of the degree of the
				      * element), as it implements the complete set of
				      * functions necessary for hp capability.
				      */
    virtual bool hp_constraints_are_implemented () const;

				     /**
				      * If, on a vertex, several finite elements are
				      * active, the hp code first assigns the degrees
				      * of freedom of each of these FEs different global
				      * indices. It then calls this function to find
				      * out which of them should get identical values,
				      * and consequently can receive the same global DoF
				      * index. This function therefore returns a list
				      * of identities between DoFs of the present finite
				      * element object with the DoFs of fe_other, which
				      * is a reference to a finite element object representing
				      * one of the other finite elements active on this
				      * particular vertex. The function computes which
				      * of the degrees of freedom of the two finite
				      * element objects are equivalent, and returns a
				      * list of pairs of global dof indices in identities.
				      * The first index of each pair denotes one of the
				      * vertex dofs of the present element, whereas the
				      * second is the corresponding index of the other
				      * finite element.
				      */
    virtual std::vector<std::pair<unsigned int, unsigned int> > hp_vertex_dof_identities (const FiniteElement<dim>& fe_other) const;

				     /**
				      * Same as hp_vertex_dof_indices(), except that
				      * the function treats degrees of freedom on lines.
				      */
    virtual std::vector<std::pair<unsigned int, unsigned int> > hp_line_dof_identities (const FiniteElement<dim>& fe_other) const;

				     /**
				      * Same as hp_vertex_dof_indices(), except that
				      * the function treats degrees of freedom on lines.
				      */
    virtual std::vector<std::pair<unsigned int, unsigned int> > hp_quad_dof_identities (const FiniteElement<dim>& fe_other) const;

				     /**
				      * Return the matrix interpolating from a face of one
				      * element to the face of the neighboring element. The
				      * size of the matrix is then <tt>source.dofs_per_face</tt>
				      * times <tt>this->dofs_per_face</tt>.
				      *
				      * Derived elements will have to implement this function.
				      * They may only provide interpolation matrices for certain
				      * source finite elements, for example those from the same
				      * family. If they don't implement interpolation from a given
				      * element, then they must throw an exception of type
				      * <tt>FiniteElement<dim>::ExcInterpolationNotImplemented</tt>.
				      */
    virtual void get_face_interpolation_matrix (const FiniteElement<dim>& source,
						FullMatrix<double>& matrix) const;

				     /**
				      * Return the matrix interpolating from a face of one element
				      * to the subface of the neighboring element. The size of
				      * the matrix is then <tt>source.dofs_per_face</tt> times
				      * <tt>this->dofs_per_face</tt>.
				      *
				      * Derived elements will have to implement this function.
				      * They may only provide interpolation matrices for certain
				      * source finite elements, for example those from the same
				      * family. If they don't implement interpolation from a given
				      * element, then they must throw an exception of type
				      * <tt>ExcInterpolationNotImplemented</tt>.
				      */
    virtual void get_subface_interpolation_matrix (const FiniteElement<dim>& source,
						   const unsigned int subface,
						   FullMatrix<double>& matrix) const;

    virtual void interpolate (std::vector<double>& local_dofs, const std::vector<double>& values) const;

    virtual void interpolate (std::vector<double>& local_dofs, const std::vector<Vector<double> >& values, unsigned int offset = 0) const;
    virtual void interpolate (std::vector<double>& local_dofs, const VectorSlice<const std::vector<std::vector<double> > >& values) const;
    virtual unsigned int memory_consumption () const;
    virtual FiniteElement<dim> * clone() const;
    
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
    static std::vector<unsigned int> get_dpo_vector (const unsigned int degree);

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
    void initialize_support_points (const unsigned int degree);

				     /**
				      * Initialize the interpolation
				      * from functions on refined mesh
				      * cells onto the father
				      * cell. According to the
				      * philosophy of the
				      * Nédélec element, this
				      * restriction operator preserves
				      * the curl of a function
				      * weakly.
				      */
    void initialize_restriction ();

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
	std::vector<std::vector<Tensor<1, dim> > > shape_values;

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
	std::vector<std::vector<Tensor<2, dim> > > shape_gradients;
    };

    const unsigned int deg;
				     /**
				      * These are the factors
				      * multiplied to a function in
				      * the
				      * #generalized_face_support_points
				      * when computing the
				      * integration.
				      */
    Table<2, double> boundary_weights;
    
				     /**
				      * Allow access from other
				      * dimensions.
				      */
    template <int dim1> friend class FE_Nedelec;
};

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
std::vector<unsigned int> FE_Nedelec<1>::get_dpo_vector (const unsigned int);
template <>
void
FE_Nedelec<1>::initialize_restriction();

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
