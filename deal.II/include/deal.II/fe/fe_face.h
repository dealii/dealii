// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#ifndef __deal2__fe_face_h
#define __deal2__fe_face_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/fe/fe_poly_face.h>

DEAL_II_NAMESPACE_OPEN


/**
 * A finite element, which is a tensor product polynomial on each face and
 * undefined in the interior of the cells. The basis functions on the faces
 * are Lagrange polynomials based on the support points of the
 * (dim-1)-dimensional Gauss--Lobatto quadrature rule. For element degree one
 * and two, the polynomials hence correspond to the usual Lagrange polynomials
 * on equidistant points.
 *
 * This finite element is the trace space of FE_RaviartThomas on the faces and
 * serves in hybridized methods, e.g. in combination with the FE_DGQ
 * element. Its use is demonstrated in the step-51 tutorial program.
 *
 * @note Since these are only finite elements on faces, only
 * FEFaceValues and FESubfaceValues will be able to extract reasonable
 * values from any face polynomial. In order to make the use of
 * FESystem simpler, FEValues objects will not fail using this finite
 * element space, but all shape function values extracted will equal
 * to zero.
 *
 * @ingroup fe
 * @author Guido Kanschat, Martin Kronbichler
 * @date 2009, 2011, 2013
 */
template <int dim, int spacedim=dim>
class FE_FaceQ : public FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Lagrange
   * polynomials in each coordinate direction.
   */
  FE_FaceQ(unsigned int p);

  virtual FiniteElement<dim,spacedim> *clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * Check for non-zero values on a face.
   *
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero values on the face @p face_index.
   *
   * Implementation of the interface in FiniteElement
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};




/**
 * A finite element, which is a Legendre on each face (i.e., FE_DGP)
 * and undefined in the interior of the cells. The basis functions on
 * the faces are from Polynomials::Legendre.
 *
 * This element is used in a hybridized method together with the FE_DGP
 * element for the interior degrees of freedom.
 *
 * @note Since these are only finite elements on faces, only
 * FEFaceValues and FESubfaceValues will be able to extract reasonable
 * values from any face polynomial. In order to make the use of
 * FESystem simpler, FEValues objects will not fail using this finite
 * element space, but all shape function values extracted will equal
 * to zero.
 *
 * @ingroup fe
 * @author Martin Kronbichler
 * @date 2013
 */
template <int dim, int spacedim=dim>
class FE_FaceP : public FE_PolyFace<PolynomialSpace<dim-1>, dim, spacedim>
{
public:
  /**
   * Constructor for complete basis of polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Legendre
   * polynomials in each coordinate direction.
   */
  FE_FaceP(unsigned int p);

  /**
   * Clone method.
   */
  virtual FiniteElement<dim,spacedim> *clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceP<dim>(degree)</tt> , with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * Check for non-zero values on a face.
   *
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero values on the face @p face_index.
   *
   * Implementation of the interface in FiniteElement
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};


DEAL_II_NAMESPACE_CLOSE

#endif
