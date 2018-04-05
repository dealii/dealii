// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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

#ifndef dealii_fe_trace_h
#define dealii_fe_trace_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/fe/fe_poly_face.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_face.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A finite element, which is the trace of FE_Q elements, that is a tensor
 * product of polynomials on the faces, undefined in the interior of the cells
 * and continuous. The basis functions on the faces are formed by a tensor
 * product of 1D Lagrange polynomials with equidistant points up to degree 2
 * and Gauss-Lobatto points starting from degree 3.
 *
 * This finite element is the trace space of FE_Q on the faces.
 *
 * @note Since these are only finite elements on faces, only FEFaceValues and
 * FESubfaceValues will be able to extract reasonable values from any face
 * polynomial. In order to make the use of FESystem simpler, FEValues objects
 * will not fail using this finite element space, but all shape function
 * values extracted will equal to zero.
 */

template <int dim, int spacedim=dim>
class FE_TraceQ : public FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Legendre
   * polynomials in each coordinate direction.
   */
  FE_TraceQ(unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  virtual
  std::unique_ptr<FiniteElement<dim,spacedim> >
  clone() const;

  /**
   * Implementation of the corresponding function in the FiniteElement
   * class.  Since the current element is interpolatory, the nodal
   * values are exactly the support point values. Furthermore, since
   * the current element is scalar, the support point values need to
   * be vectors of length 1.
   */
  virtual
  void
  convert_generalized_support_point_values_to_dof_values (const std::vector<Vector<double> > &support_point_values,
                                                          std::vector<double>                &nodal_values) const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Return a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return the matrix interpolating from a face of one element to the face
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
   * Return the matrix interpolating from a face of one element to the face
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
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementDomination::Domination
   * and in particular the
   * @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

private:
  /**
   * Store a copy of FE_Q for delegating the hp-constraints functionality.
   */
  FE_Q<dim, spacedim> fe_q;

  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};



/**
 * FE_TraceQ in 1D, i.e., with degrees of freedom on the element vertices.
 */
template <int spacedim>
class FE_TraceQ<1,spacedim> : public FE_FaceQ<1,spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_TraceQ (const unsigned int p);

  /**
   * Return the name of the element
   */
  std::string get_name() const;
};


DEAL_II_NAMESPACE_CLOSE

#endif
