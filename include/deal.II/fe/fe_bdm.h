// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_bdm_h
#define dealii_fe_bdm_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * The Brezzi-Douglas-Marini element.
 *
 * <h3>Degrees of freedom</h3>
 *
 * @todo The 3D version exhibits some numerical instabilities, in particular
 * for higher order
 *
 * @todo Restriction matrices are missing.
 *
 * The matching pressure space for FE_BDM of order <i>k</i> is the element
 * FE_DGP of order <i>k-1</i>.
 *
 * The BDM element of order @p p has <i>p+1</i> degrees of freedom on each
 * face. These are implemented as the function values in the <i>p+1</i> Gauss
 * points on each face.
 *
 * Additionally, for order greater or equal 2, we have additional
 * <i>p(p-1)</i>, the number of vector valued polynomials in
 * <i>P<sub>p</sub></i>, interior degrees of freedom. These are the vector
 * function values in the first <i>p(p-1)/2</i> of the <i>p<sup>2</sup></i>
 * Gauss points in the cell.
 *
 * @ingroup fe
 */
template <int dim>
class FE_BDM : public FE_PolyTensor<dim>
{
public:
  /**
   * Constructor for the BDM element of degree @p p.
   */
  FE_BDM(const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_BDM<dim>(degree)</tt>, with @p dim and @p degree replaced
   * by appropriate values.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

private:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * Compute the vector used for the @p restriction_is_additive field passed
   * to the base class's constructor.
   */
  static std::vector<bool>
  get_ria_vector(const unsigned int degree);
  /**
   * Initialize the FiniteElement<dim>::generalized_support_points and
   * FiniteElement<dim>::generalized_face_support_points fields. Called from
   * the constructor. See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  void
  initialize_support_points(const unsigned int bdm_degree);
  /**
   * The values in the face support points of the polynomials needed as test
   * functions. The outer vector is indexed by quadrature points, the inner by
   * the test function. The test function space is PolynomialsP<dim-1>.
   */
  std::vector<std::vector<double>> test_values_face;
  /**
   * The values in the interior support points of the polynomials needed as
   * test functions. The outer vector is indexed by quadrature points, the
   * inner by the test function. The test function space is PolynomialsP<dim>.
   */
  std::vector<std::vector<double>> test_values_cell;
};

DEAL_II_NAMESPACE_CLOSE

#endif
