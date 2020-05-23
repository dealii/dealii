// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2019 by the deal.II authors
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

#ifndef dealii_fe_raviart_thomas_h
#define dealii_fe_raviart_thomas_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of Raviart-Thomas (RT) elements. The Raviart-Thomas space
 * is designed to solve problems in which the solution only lives in the
 * space
 * $H^\text{div}=\{ {\mathbf u} \in L_2: \text{div}\, {\mathbf u} \in L_2\}$,
 * rather than in the more commonly used space
 * $H^1=\{ u \in L_2: \nabla u \in L_2\}$. In other words, the solution must
 * be a vector field whose divergence is square integrable, but for which the
 * gradient may not be square integrable. The typical application for this
 * space (and these elements) is to the mixed formulation of the Laplace
 * equation and related situations, see for example step-20. The defining
 * characteristic of functions in $H^\text{div}$ is that they are in
 * general discontinuous -- but that if you draw a line in 2d (or a
 * surface in 3d), then the <i>normal</i> component of the vector
 * field must be continuous across the line (or surface) even though
 * the tangential component may not be. As a consequence, the
 * Raviart-Thomas element is constructed in such a way that (i) it is
 * @ref vector_valued "vector-valued", (ii) the shape functions are
 * discontinuous, but (iii) the normal component of the vector field
 * represented by each shape function is continuous across the faces
 * of cells.
 *
 * Other properties of the Raviart-Thomas element are that (i) it is
 * @ref GlossPrimitive "not a primitive element"; (ii) the shape functions
 * are defined so that certain integrals over the faces are either zero
 * or one, rather than the common case of certain point values being
 * either zero or one. (There is, however, the FE_RaviartThomasNodal
 * element that uses point values.)
 *
 * We follow the commonly used -- though confusing -- definition of the "degree"
 * of RT elements. Specifically, the "degree" of the element denotes
 * the polynomial degree of the <i>largest complete polynomial subspace</i>
 * contained in the finite element space, even if the space may contain shape
 * functions of higher polynomial degree. The lowest order element is
 * consequently FE_RaviartThomas(0), i.e., the Raviart-Thomas element "of
 * degree zero", even though the functions of this space are in general
 * polynomials of degree one in each variable. This choice of "degree"
 * implies that the approximation order of the function itself is
 * <i>degree+1</i>, as with usual polynomial spaces. The numbering so chosen
 * implies the sequence
 * @f[
 *   Q_{k+1}
 *   \stackrel{\text{grad}}{\rightarrow}
 *   \text{Nedelec}_k
 *   \stackrel{\text{curl}}{\rightarrow}
 *   \text{RaviartThomas}_k
 *   \stackrel{\text{div}}{\rightarrow}
 *   DGQ_{k}
 * @f]
 *
 * This class is not implemented for the codimension one case (<tt>spacedim !=
 * dim</tt>).
 *
 *
 * <h3>Interpolation</h3>
 *
 * The
 * @ref GlossInterpolation "interpolation"
 * operators associated with the RT element are constructed such that
 * interpolation and computing the divergence are commuting operations. We
 * require this from interpolating arbitrary functions as well as the
 * #restriction matrices.  It can be achieved by two interpolation schemes,
 * the simplified one in FE_RaviartThomasNodal and the original one here:
 *
 * <h4>Node values on edges/faces</h4>
 *
 * On edges or faces, the
 * @ref GlossNodes "node values"
 * are the moments of the normal component of the interpolated function with
 * respect to the traces of the RT polynomials. Since the normal trace of the
 * RT space of degree <i>k</i> on an edge/face is the space
 * <i>Q<sub>k</sub></i>, the moments are taken with respect to this space.
 *
 * <h4>Interior node values</h4>
 *
 * Higher order RT spaces have interior nodes. These are moments taken with
 * respect to the gradient of functions in <i>Q<sub>k</sub></i> on the cell
 * (this space is the matching space for RT<sub>k</sub> in a mixed
 * formulation).
 *
 * <h4>Generalized support points</h4>
 *
 * The node values above rely on integrals, which will be computed by
 * quadrature rules themselves. The generalized support points are a set of
 * points such that this quadrature can be performed with sufficient accuracy.
 * The points needed are those of QGauss<sub>k+1</sub> on each face as well as
 * QGauss<sub>k+1</sub> in the interior of the cell (or none for
 * RT<sub>0</sub>).
 *
 *
 * @author Guido Kanschat, 2005, based on previous work by Wolfgang Bangerth.
 */
template <int dim>
class FE_RaviartThomas : public FE_PolyTensor<dim>
{
public:
  /**
   * Constructor for the Raviart-Thomas element of degree @p p.
   */
  FE_RaviartThomas(const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_RaviartThomas<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  // documentation inherited from the base class
  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   *
   * Right now, this is only implemented for RT0 in 1D. Otherwise, returns
   * always @p true.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * Return a list of constant modes of the element. This method is currently
   * not correctly implemented because it returns ones for all components.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  virtual std::size_t
  memory_consumption() const override;

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
   * Initialize the @p generalized_support_points field of the FiniteElement
   * class and fill the tables with interpolation weights (#boundary_weights
   * and #interior_weights). Called from the constructor.
   */
  void
  initialize_support_points(const unsigned int rt_degree);

  /**
   * Initialize the interpolation from functions on refined mesh cells onto
   * the father cell. According to the philosophy of the Raviart-Thomas
   * element, this restriction operator preserves the divergence of a function
   * weakly.
   */
  void
  initialize_restriction();

  /**
   * These are the factors multiplied to a function in the
   * #generalized_face_support_points when computing the integration. They are
   * organized such that there is one row for each generalized face support
   * point and one column for each degree of freedom on the face.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  Table<2, double> boundary_weights;

  /**
   * Precomputed factors for interpolation of interior degrees of freedom. The
   * rationale for this Table is the same as for #boundary_weights. Only, this
   * table has a third coordinate for the space direction of the component
   * evaluated.
   */
  Table<3, double> interior_weights;

  // Allow access from other dimensions.
  template <int dim1>
  friend class FE_RaviartThomas;
};



/**
 * The Raviart-Thomas elements with node functionals defined as point values
 * in Gauss points.
 *
 * <h3>Description of node values</h3>
 *
 * For this Raviart-Thomas element, the node values are not cell and face
 * moments with respect to certain polynomials, but the values in quadrature
 * points. Following the general scheme for numbering degrees of freedom, the
 * node values on edges are first, edge by edge, according to the natural
 * ordering of the edges of a cell. The interior degrees of freedom are last.
 *
 * For an RT-element of degree <i>k</i>, we choose <i>(k+1)<sup>d-1</sup></i>
 * Gauss points on each face. These points are ordered lexicographically with
 * respect to the orientation of the face. This way, the normal component
 * which is in <i>Q<sub>k</sub></i> is uniquely determined. Furthermore, since
 * this Gauss-formula is exact on <i>Q<sub>2k+1</sub></i>, these node values
 * correspond to the exact integration of the moments of the RT-space.
 *
 * In the interior of the cells, the moments are with respect to an
 * anisotropic <i>Q<sub>k</sub></i> space, where the test functions are one
 * degree lower in the direction corresponding to the vector component under
 * consideration. This is emulated by using an anisotropic Gauss formula for
 * integration.
 *
 * @todo The current implementation is for Cartesian meshes only. You must use
 * MappingCartesian.
 *
 * @todo Even if this element is implemented for two and three space
 * dimensions, the definition of the node values relies on consistently
 * oriented faces in 3D. Therefore, care should be taken on complicated
 * meshes.
 *
 * @note The degree stored in the member variable
 * FiniteElementData<dim>::degree is higher by one than the constructor
 * argument!
 *
 * @author Guido Kanschat, 2005, Zhu Liang, 2008
 */
template <int dim>
class FE_RaviartThomasNodal : public FE_PolyTensor<dim>
{
public:
  /**
   * Constructor for the Raviart-Thomas element of degree @p p.
   */
  FE_RaviartThomasNodal(const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_RaviartThomasNodal<dim>(degree)</tt>, with @p dim and @p
   * degree replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  // documentation inherited from the base class
  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim> &source,
                                FullMatrix<double> &matrix) const override;

  virtual void
  get_subface_interpolation_matrix(const FiniteElement<dim> &source,
                                   const unsigned int        subface,
                                   FullMatrix<double> &matrix) const override;
  virtual bool
  hp_constraints_are_implemented() const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim> &fe_other,
                         const unsigned int codim = 0) const override final;

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
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   *
   * Right now, this is only implemented for RT0 in 1D. Otherwise, returns
   * always @p true.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;
  /**
   * Initialize the FiniteElement<dim>::generalized_support_points and
   * FiniteElement<dim>::generalized_face_support_points fields. Called from
   * the constructor.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  void
  initialize_support_points(const unsigned int rt_degree);
};


/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
FE_RaviartThomas<1>::initialize_restriction();

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
