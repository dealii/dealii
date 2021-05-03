// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2019 by the deal.II authors
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

#ifndef dealii_fe_nedelec_h
#define dealii_fe_nedelec_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup fe */
/*@{*/

/**
 * @warning Several aspects of the implementation are experimental. For the
 * moment, it is safe to use the element on globally refined meshes with
 * consistent orientation of faces. See the todo entries below for more
 * detailed caveats.
 *
 * Implementation of N&eacute;d&eacute;lec elements. The N&eacute;d&eacute;lec
 * space is designed to solve problems in which the solution only lives in the
 * space
 * $H^\text{curl}=\{ {\mathbf u} \in L_2: \text{curl}\, {\mathbf u} \in L_2\}$,
 * rather than in the more commonly used space
 * $H^1=\{ u \in L_2: \nabla u \in L_2\}$. In other words, the solution must
 * be a vector field whose curl is square integrable, but for which the
 * gradient may not be square integrable. The typical application for this
 * space (and these elements) is to the Maxwell equations and corresponding
 * simplifications, such as the reduced version of the Maxwell equation
 * that only involves the electric field $\mathbf E$ which has to satisfy
 * the equation $\text{curl}\, \text{curl}\, {\mathbf E} = 0$ in the
 * time independent case when no currents are present, or the equation
 * $\text{curl}\,\text{curl}\,{\mathbf A} = 4\pi{\mathbf j}$ that the
 * magnetic vector potential $\mathbf A$ has to satisfy in the
 * time independent case.
 *
 * The defining
 * characteristic of functions in $H^\text{curl}$ is that they are in
 * general discontinuous -- but that if you draw a line in 2d (or a
 * surface in 3d), then the <i>tangential</i> component(s) of the vector
 * field must be continuous across the line (or surface) even though
 * the normal component may not be. As a consequence, the
 * N&eacute;d&eacute;lec element is constructed in such a way that (i) it is
 * @ref vector_valued "vector-valued", (ii) the shape functions are
 * discontinuous, but (iii) the tangential component(s) of the vector field
 * represented by each shape function are continuous across the faces
 * of cells.
 *
 * Other properties of the N&eacute;d&eacute;lec element are that (i) it is
 * @ref GlossPrimitive "not a primitive element"; (ii) the shape functions
 * are defined so that certain integrals over the faces are either zero
 * or one, rather than the common case of certain point values being
 * either zero or one.
 *
 * We follow the commonly used -- though confusing -- definition of the "degree"
 * of N&eacute;d&eacute;lec elements. Specifically, the "degree" of the element
 * denotes the polynomial degree of the <i>largest complete polynomial
 * subspace</i> contained in the finite element space, even if the space may
 * contain shape functions of higher polynomial degree. The lowest order element
 * is consequently FE_Nedelec(0), i.e., the Raviart-Thomas element "of degree
 * zero", even though the functions of this space are in general polynomials of
 * degree one in each variable. This choice of "degree" implies that the
 * approximation order of the function itself is <i>degree+1</i>, as with usual
 * polynomial spaces. The numbering so chosen implies the sequence
 * @f[
 *   Q_{k+1}
 *   \stackrel{\text{grad}}{\rightarrow}
 *   \text{Nedelec}_k
 *   \stackrel{\text{curl}}{\rightarrow}
 *   \text{RaviartThomas}_k
 *   \stackrel{\text{div}}{\rightarrow}
 *   DGQ_{k}
 * @f]
 * Note that this follows the convention of Brezzi and Raviart,
 * though not the one used in the original paper by N&eacute;d&eacute;lec.
 *
 * This class is not implemented for the codimension one case (<tt>spacedim !=
 * dim</tt>).
 *
 * @todo Even if this element is implemented for two and three space
 * dimensions, the definition of the node values relies on consistently
 * oriented faces in 3D. Therefore, care should be taken on complicated
 * meshes.
 *
 *
 * <h3>Interpolation</h3>
 *
 * The
 * @ref GlossInterpolation "interpolation"
 * operators associated with the N&eacute;d&eacute;lec element are constructed
 * such that interpolation and computing the curl are commuting operations on
 * rectangular mesh cells. We require this from interpolating arbitrary
 * functions as well as the #restriction matrices.
 *
 * <h4>Node values</h4>
 *
 * The
 * @ref GlossNodes "node values"
 * for an element of degree <i>k</i> on the reference cell are:
 * <ol>
 * <li> On edges: the moments of the tangential component with respect to
 * polynomials of degree <i>k</i>.
 * <li> On faces: the moments of the tangential components with respect to
 * <tt>dim</tt>-1 dimensional FE_Nedelec polynomials of degree <i>k</i>-1.
 * <li> In cells: the moments with respect to gradients of polynomials in FE_Q
 * of degree <i>k</i>.
 * </ol>
 *
 * <h4>Generalized support points</h4>
 *
 * The node values above rely on integrals, which will be computed by
 * quadrature rules themselves. The generalized support points are a set of
 * points such that this quadrature can be performed with sufficient accuracy.
 * The points needed are those of QGauss<sub>k+1</sub> on each edge and
 * QGauss<sub>k+2</sub> on each face and in the interior of the cell (or none
 * for N<sub>1</sub>).
 */
template <int dim>
class FE_Nedelec : public FE_PolyTensor<dim>
{
public:
  /**
   * Constructor for the Nedelec element of given @p order. The maximal
   * polynomial degree of the shape functions is `order+1` (in each variable;
   * the total polynomial degree may be higher). If `order = 0`, the element is
   * linear and has degrees of freedom only on the edges. If `order >=1` the
   * element has degrees of freedom on the edges, faces and volume. For example
   * the 3D version of FE_Nedelec has 12 degrees of freedom for `order = 0`
   * and 54 for `degree = 1`. It is important to have enough quadrature points
   * in order to perform the quadrature with sufficient accuracy.
   * For example [QGauss<dim>(order + 2)](@ref QGauss) can be used for the
   * quadrature formula, where `order` is the order of FE_Nedelec.
   */
  FE_Nedelec(const unsigned int order);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Nedelec<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;


  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp-compatible".
   *
   * For the <tt>FE_Nedelec</tt> class the result is always true (independent
   * of the degree of the element), as it implements the complete set of
   * functions necessary for hp-capability.
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * If, on a vertex, several finite elements are active, the hp-code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of
   * n_dofs_per_vertex() of the two finite elements. The first index of each
   * pair denotes one of the vertex dofs of the present element, whereas the
   * second is the corresponding index of the other finite element.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                         const unsigned int        face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face of
   * the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * <tt>FiniteElement<dim>::ExcInterpolationNotImplemented</tt>.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim> &source,
                                FullMatrix<double> &      matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the subface
   * of the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * <tt>ExcInterpolationNotImplemented</tt>.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim> &source,
    const unsigned int        subface,
    FullMatrix<double> &      matrix,
    const unsigned int        face_no = 0) const override;

  /**
   * Projection from a fine grid space onto a coarse grid space. If this
   * projection operator is associated with a matrix @p P, then the
   * restriction of this matrix @p P_i to a single child cell is returned
   * here.
   *
   * The matrix @p P is the concatenation or the sum of the cell matrices @p
   * P_i, depending on the #restriction_is_additive_flags. This distinguishes
   * interpolation (concatenation) and projection with respect to scalar
   * products (summation).
   *
   * Row and column indices are related to coarse grid and fine grid spaces,
   * respectively, consistent with the definition of the associated operator.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Embedding matrix between grids.
   *
   * The identity operator from a coarse grid space into a fine grid space is
   * associated with a matrix @p P. The restriction of this matrix @p P_i to a
   * single child cell is returned here.
   *
   * The matrix @p P is the concatenation, not the sum of the cell matrices @p
   * P_i. That is, if the same non-zero entry <tt>j,k</tt> exists in two
   * different child matrices @p P_i, the value should be the same in both
   * matrices and it is copied into the matrix @p P only once.
   *
   * Row and column indices are related to fine grid and coarse grid spaces,
   * respectively, consistent with the definition of the associated operator.
   *
   * These matrices are used by routines assembling the prolongation matrix
   * for multi-level methods.  Upon assembling the transfer matrix between
   * cells using this matrix array, zero elements in the prolongation matrix
   * are discarded and will not fill up the transfer matrix.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * Return a list of constant modes of the element.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  virtual std::size_t
  memory_consumption() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

private:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   *
   * If the optional argument <tt>dg</tt> is true, the vector returned will
   * have all degrees of freedom assigned to the cell, none on the faces and
   * edges.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree, bool dg = false);

  /**
   * Initialize the @p generalized_support_points field of the FiniteElement
   * class and fill the tables with interpolation weights (#boundary_weights
   * and interior_weights). Called from the constructor.
   */
  void
  initialize_support_points(const unsigned int order);

  /**
   * Initialize the interpolation from functions on refined mesh cells onto
   * the father cell. According to the philosophy of the Nédélec element,
   * this restriction operator preserves the curl of a function weakly.
   */
  void
  initialize_restriction();

  /**
   * These are the factors multiplied to a function in the
   * #generalized_face_support_points when computing the integration.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  Table<2, double> boundary_weights;

  /**
   * Mutex for protecting initialization of restriction and embedding matrix.
   */
  mutable Threads::Mutex mutex;

  /**
   * Initialize the permutation pattern and the pattern of sign change.
   *
   * @note This function is not fully filled with the correct implementation
   * yet. It needs to be consistently implemented in a future release to work
   * on meshes that contain cells with flipped faces.
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();

  // Allow access from other dimensions.
  template <int dim1>
  friend class FE_Nedelec;
};

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
FE_Nedelec<1>::initialize_restriction();

#endif // DOXYGEN

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
