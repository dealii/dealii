// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__fe_dg_vector_h
#define dealii__fe_dg_vector_h

#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * DG elements based on vector valued polynomials.
 *
 * These elements use vector valued polynomial spaces as they have been
 * introduced for H<sup>div</sup> and H<sup>curl</sup> conforming finite
 * elements, but do not use the usual continuity of these elements. Thus, they
 * are suitable for DG and hybrid formulations involving these function
 * spaces.
 *
 * The template argument <tt>PolynomialType</tt> refers to a vector valued
 * polynomial space like PolynomialsRaviartThomas or PolynomialsNedelec. Note
 * that the dimension of the polynomial space and the argument <tt>dim</tt>
 * must coincide.
 *
 * @ingroup febase
 * @author Guido Kanschat
 * @date 2010
 */
template <class PolynomialType, int dim, int spacedim=dim>
class FE_DGVector
  :
  public FE_PolyTensor<PolynomialType, dim, spacedim>
{
public:
  /**
   * Constructor for the vector element of degree @p p.
   */
  FE_DGVector (const unsigned int p, MappingType m);
public:

  FiniteElement<dim, spacedim> *clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_RaviartThomas<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string get_name () const;


  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   *
   * For this element, we always return @p true.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  virtual void interpolate(std::vector<double>                &local_dofs,
                           const std::vector<double> &values) const;
  virtual void interpolate(std::vector<double>                &local_dofs,
                           const std::vector<Vector<double> > &values,
                           unsigned int offset = 0) const;
  virtual void interpolate(
    std::vector<double> &local_dofs,
    const VectorSlice<const std::vector<std::vector<double> > > &values) const;
  virtual std::size_t memory_consumption () const;

private:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector (const unsigned int degree);

  /**
   * Initialize the @p generalized_support_points field of the FiniteElement
   * class and fill the tables with @p interior_weights. Called from the
   * constructor.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  void initialize_support_points (const unsigned int degree);

  /**
   * Initialize the interpolation from functions on refined mesh cells onto
   * the father cell. According to the philosophy of the Raviart-Thomas
   * element, this restriction operator preserves the divergence of a function
   * weakly.
   */
  void initialize_restriction ();

  /**
   * Fields of cell-independent data.
   *
   * For information about the general purpose of this class, see the
   * documentation of the base class.
   */
  class InternalData : public FiniteElement<dim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function values in quadrature points. There is one row
     * for each shape function, containing values for each quadrature point.
     * Since the shape functions are vector-valued (with as many components as
     * there are space dimensions), the value is a tensor.
     *
     * In this array, we store the values of the shape function in the
     * quadrature points on the unit cell. The transformation to the real
     * space cell is then simply done by multiplication with the Jacobian of
     * the mapping.
     */
    std::vector<std::vector<Tensor<1,dim> > > shape_values;

    /**
     * Array with shape function gradients in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     *
     * We store the gradients in the quadrature points on the unit cell. We
     * then only have to apply the transformation (which is a matrix-vector
     * multiplication) when visiting an actual cell.
     */
    std::vector<std::vector<Tensor<2,dim> > > shape_gradients;
  };
  Table<3, double> interior_weights;
};



/**
 * A vector-valued DG element based on the polynomials space of FE_Nedelec.
 *
 * @ingroup fe
 * @author Guido Kanschat
 * @date 2011
 */
template <int dim, int spacedim=dim>
class FE_DGNedelec : public FE_DGVector<PolynomialsNedelec<dim>, dim, spacedim>
{
public:
  /**
   * Constructor for the discontinuous N&eacute;d&eacute;lec element of degree
   * @p p.
   */
  FE_DGNedelec (const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGNedelec<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string get_name () const;
};



/**
 * A vector-valued DG element based on the polynomials space of
 * FE_RaviartThomas.
 *
 * @ingroup fe
 * @author Guido Kanschat
 * @date 2011
 */
template <int dim, int spacedim=dim>
class FE_DGRaviartThomas : public FE_DGVector<PolynomialsRaviartThomas<dim>, dim, spacedim>
{
public:
  /**
   * Constructor for the Raviart-Thomas element of degree @p p.
   */
  FE_DGRaviartThomas (const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGRaviartThomas<dim>(degree)</tt>, with @p dim and @p
   * degree replaced by appropriate values.
   */
  virtual std::string get_name () const;
};



/**
 * A vector-valued DG element based on the polynomials space of FE_BDM.
 *
 * @ingroup fe
 * @author Guido Kanschat
 * @date 2011
 */
template <int dim, int spacedim=dim>
class FE_DGBDM : public FE_DGVector<PolynomialsBDM<dim>, dim, spacedim>
{
public:
  /**
   * Constructor for the discontinuous BDM element of degree @p p.
   */
  FE_DGBDM (const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGBDM<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string get_name () const;
};


DEAL_II_NAMESPACE_CLOSE

#endif
