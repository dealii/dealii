//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_bdm_h
#define __deal2__fe_bdm_h

#include <base/config.h>
#include <base/table.h>
#include <base/polynomials_bdm.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <base/geometry_info.h>
#include <fe/fe.h>
#include <fe/fe_poly_tensor.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * The Brezzi-Douglas-Marini element.
 *
 * <h3>Degrees of freedom</h3>
 *
 * @todo This is for 2D only.
 *
 * @todo Transformation works only for uniform, Cartesian meshes
 *
 * The matching pressure psace for FE_BDM of order <i>k</i> is the
 * element FE_DGP of order <i>k</i>.
 *
 * The BDM element of order @p p has <i>p+1</i> degrees of freedom on
 * each face. These are implemented as the function values in the
 * <i>p+1</i> Gauss points on each face.
 *
 * Additionally, for order greater or equal 2, we have additional
 * <i>p(p-1)</i>, the number of vector valued polynomials in
 * <i>P<sub>p</sub></i>, interior degrees of freedom. These are the
 * vector function values in the first <i>p(p-1)/2</i> of the
 * <i>p<sup>2</sup></i> Gauss points in the cell.
 */
template <int dim>
class FE_BDM
  :
  public FE_PolyTensor<PolynomialsBDM<dim>, dim>
{
  public:
				     /**
				      * Constructor for the BDM
				      * element of degree @p p.
				      */
    FE_BDM (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_BDM<dim>(degree)</tt>, with
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
				     /**
				      * The values in the interior
				      * support points of the
				      * polynomials needed as test
				      * functions. The outer vector is
				      * indexed by quadrature points,
				      * the inner by the test
				      * function.
				      */
    std::vector<std::vector<double> > test_values;
};

DEAL_II_NAMESPACE_CLOSE

#endif
