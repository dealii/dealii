// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii__fe_raviart_thomas_bubbles_h
#define dealii__fe_raviart_thomas_bubbles_h

#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/polynomials_rt_bubbles.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * Implementation of curl-enhanced Raviart-Thomas (RT_Bubbles) elements,
 * conforming with H<sup>div</sup> space. The node functionals are defined as
 * point values in Gauss-Lobatto points. These elements generate vector fields
 * with normal components continuous between mesh cells. The main use of these
 * is in providing the localization of interactions of degrees of freedom
 * around the nodes when an appropriate quadrature rule is used
 *
 * The elements are defined through enrichment of classical Raviart-Thomas
 * elements with extra curls, so that the H<sup>div</sup> conformity is
 * preserved, and the total number of degrees of freedom of RT_Bubbles(k)
 * becomes equal to the number of DoFs of <i>dim * Q<sub>k</sub></i>.
 * The lowest order element is FE_RT_bubbles(1). The matching pressure space
 * for FE_RT_Bubbles(k) is FE_DGQ(k-1).
 *
 * For this enhanced Raviart-Thomas element, the node values are not cell
 * and face moments with respect to certain polynomials, but the values in
 * Gauss-Lobatto quadrature points. The node values on edges are first,
 * edge by edge, according to the natural ordering of the edges of a cell.
 * The interior degrees of freedom are last.
 *
 * For an RT_Bubbles-element of degree <i>k</i>, we choose <i>(k+1)<sup>d-1</sup></i>
 * Gauss-Lobatto points on each face. These points are ordered lexicographically
 * with respect to the orientation of the face. In the interior of the cells,
 * the moments are computed using an anisotropic Gauss-Lobatto formula for
 * integration.
 *
 * @todo The current implementation is for Cartesian meshes only. You must use
 * MappingCartesian.
 *
 * @todo Even if this element is implemented for two and three space
 * dimensions, the definition of the node values relies on consistently
 * oriented faces in 3D. Therefore, care should be taken on complicated
 * meshes.
 * @todo Implement restriction matrices
 *
 * @author Eldar Khattatov, Ilona Ambartsumyan, 2017
 */
template <int dim>
class FE_RT_Bubbles
  :
  public FE_PolyTensor<PolynomialsRT_Bubbles<dim>, dim>
{
public:
  /**
   * Constructor for the RT_Bubbles element of degree @p p.
   */
  FE_RT_Bubbles (const unsigned int p);

  /**
   * Returns a string that uniquely identifies a finite element. This class
   * returns <tt>FE_RT_Bubbles<dim>(degree)</tt>, with @p dim and @p
   * degree replaced by appropriate values.
   */
  virtual std::string get_name () const;

  virtual std::unique_ptr<FiniteElement<dim,dim>> clone () const;

  // documentation inherited from the base class
  virtual
  void
  convert_generalized_support_point_values_to_nodal_values (const std::vector<Vector<double> > &support_point_values,
                                                            std::vector<double>                &nodal_values) const;

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
   * Compute the vector used for the @p restriction_is_additive field passed
   * to the base class's constructor.
   */
  static std::vector<bool>
  get_ria_vector (const unsigned int degree);

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   *
   * Right now, this is only implemented for RT0 in 1D. Otherwise, returns
   * always @p true.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;
  /**
   * Initialize the FiniteElement<dim>::generalized_support_points and
   * FiniteElement<dim>::generalized_face_support_points fields. Called from
   * the constructor.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  void initialize_support_points (const unsigned int rt_degree);
};


/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
