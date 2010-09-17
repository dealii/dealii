//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_dg_vector_h
#define __deal2__fe_dg_vector_h

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


/**
 * DG elements based on vector valued polynomials.
 *
 * These elements use vector valued polynomial spaces as they have
 * been introduced for H<sup>div</sup> and H<sup>curl</sup> conforming
 * finite elements, but do not use the usual continuity of these
 * elements. Thus, they are suitable for DG and hybrid formulations
 * involving these function spaces.
 *
 * The template argument <tt>POLY</tt> refers to a vector valued
 * polynomial space like PolynomialsRaviartThomas or
 * PolynomialsNedelec. Note that the dimension of the polynomial space
 * and the argument <tt>dim</tt> must coincide.
 *
 * @ingroup fe
 * @author Guido Kanschat
 * @date 2010
 */
template <class POLY, int dim, int spacedim=dim>
class FE_DGVector
  :
  public FE_PolyTensor<POLY, dim, spacedim>
{
  public:
				     /**
				      * Constructor for the vector
				      * element of degree @p p.
				      */
    FE_DGVector (const unsigned int p, MappingType m);
  public:
    
    FiniteElement<dim, spacedim>* clone() const;
    
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
				      * Check whether a shape function
				      * may be non-zero on a face.
				      *
				      * Returns always
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
				      * Initialize the @p
				      * generalized_support_points
				      * field of the FiniteElement
				      * class and fill the tables with
				      * #interior_weights. Called
				      * from the constructor.
				      */
    void initialize_support_points (const unsigned int degree);

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
    Table<3, double> interior_weights;    
};


/* -------------- declaration of explicit specializations ------------- */

DEAL_II_NAMESPACE_CLOSE

#endif
