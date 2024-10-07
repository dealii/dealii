// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_nedelec_h
#define dealii_fe_nedelec_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup fe
 * @{
 */

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
 * @ref vector_valued "vector-valued",
 * (ii) the shape functions are
 * discontinuous, but (iii) the tangential component(s) of the vector field
 * represented by each shape function are continuous across the faces
 * of cells.
 *
 * Other properties of the N&eacute;d&eacute;lec element are that (i) it is
 * @ref GlossPrimitive "not a primitive element"
 * ; (ii) the shape functions
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
 * oriented faces in 3d. Therefore, care should be taken on complicated
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
 * The points needed are those of QGauss(k+1) on each edge and
 * QGauss(k+2) on each face and in the interior of the cell (or none
 * for FE_Nedelec(0)).
 *
 * <h3> Depictions of shape functions </h3>
 *
 * The following subsections depict the shape functions defined by this class on
 * the unit cell. The figures below illustrate the direction and magnitude of
 * these shape functions.
 *
 * <h4>FE_Nedelec(0)</h4>
 *
 * For the lowest order N&eacute;d&eacute;lec element, we have a single shape
 * function associated with each edge (i.e., the tangential component of each
 * shape function is non-zero on only one edge).
 *
 * In 2d, these shape functions look as follows: <table> <tr> <td
 * align="center">
 * @image html fe_nedelec_shape_function_0_00.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_0_01.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(0) element, shape function 0
 * </td>
 *
 * <td align="center"> FE_Nedelec(0) element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_0_02.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_0_03.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(0) element, shape function 2
 * </td>
 *
 * <td align="center"> FE_Nedelec(0) element, shape function 3 </td> </tr>
 * </table>
 *
 * <h4>FE_Nedelec(1)</h4>
 *
 * For higher order N&eacute;d&eacute;lec cells, we have shape functions
 * associated with the edges, faces, and the volume.
 *
 * In 2d, for example, with FE_Nedelec(1), we have 2 shape functions associated
 * with each edge, and 4 shape functions associated with the cell, which
 * correspond to the shape functions with no non-zero tangential components on
 * the boundary of the cell.
 *
 * These shape functions look
 * as follows: <table> <tr> <td align="center">
 * @image html fe_nedelec_shape_function_1_00.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_1_01.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(1) element, shape function 0
 * </td>
 *
 * <td align="center"> FE_Nedelec(1) element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_1_02.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_1_03.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(1) element, shape function 2
 * </td>
 *
 * <td align="center"> FE_Nedelec(1) element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_1_04.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_1_05.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(1) element, shape function 4
 * </td>
 *
 * <td align="center"> FE_Nedelec(1) element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_1_06.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_1_07.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(1) element, shape function 6
 * </td>
 *
 * <td align="center"> FE_Nedelec(1) element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_1_08.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_1_09.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(1) element, shape function 8
 * </td>
 *
 * <td align="center"> FE_Nedelec(1) element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_1_10.png
 * </td>
 * <td align="center">
 * @image html fe_nedelec_shape_function_1_11.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(1) element, shape function 10
 * </td>
 *
 * <td align="center"> FE_Nedelec(1) element, shape function 11 </td> </table>
 *
 * <h4>FE_Nedelec(2)</h4>
 *
 * For higher order N&eacute;d&eacute;lec cells, we have shape functions
 * associated with the edges, faces, and the volume.
 *
 * In 2d, with FE_Nedelec(2), we have 3 shape functions associated with each
 * edge, and 12 shape functions associated with the cell.
 *
 * These shape functions look
 * as follows: <table> <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_00.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_01.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 0
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_02.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_03.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 2
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_04.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_05.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 4
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_06.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_07.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 6
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_08.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_09.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 8
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_10.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_11.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 10
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_12.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_13.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 12
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_14.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_15.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 14
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 15 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_16.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_17.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 16
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 17 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_18.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_19.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 18
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 19 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_20.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_21.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 20
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 21 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html fe_nedelec_shape_function_2_22.png
 * </td>
 *
 * <td align="center">
 * @image html fe_nedelec_shape_function_2_23.png
 * </td> </tr> <tr> <td align="center"> FE_Nedelec(2) element, shape function 22
 * </td>
 *
 * <td align="center"> FE_Nedelec(2) element, shape function 23 </td> </table>
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
   * the 3d version of FE_Nedelec has 12 degrees of freedom for `order = 0`
   * and 54 for `degree = 1`. It is important to have enough quadrature points
   * in order to perform the quadrature with sufficient accuracy.
   * For example
   * [QGauss<dim>(order + 2)](@ref QGauss)
   * can be used for the
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
                                FullMatrix<double>       &matrix,
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
    FullMatrix<double>       &matrix,
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
    std::vector<double>               &nodal_values) const override;

  /**
   * Return a list of constant modes of the element.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  virtual std::size_t
  memory_consumption() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  /**
   * For a finite element of degree larger than @p sub_degree, we return a
   * vector which maps the numbering on an FE of degree @p sub_degree into the
   * numbering on this element.
   *
   * Note that for the Nedelec element, by @p sub_degree,
   * we refer to the maximal polynomial degree (in any coordinate direction) as
   * opposed to the Nedelec degree.
   */
  std::vector<unsigned int>
  get_embedding_dofs(const unsigned int sub_degree) const;

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
   * Mutex variables used for protecting the initialization of restriction
   * and embedding matrices.
   */
  mutable Threads::Mutex restriction_matrix_mutex;
  mutable Threads::Mutex prolongation_matrix_mutex;

  /**
   * Initialize the permutation pattern and the pattern of sign change.
   *
   * @note Currently this function is implemented for the finite elements of
   * the order k < 4.
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();

  // Allow access from other dimensions.
  template <int dim1>
  friend class FE_Nedelec;

  // Information related to face orientation ----------------------------------
  //
  // The order of the Nedelec elements equals the tensor degree minus one,
  // k = n - 1. In the three-dimensional space the Nedelec elements of the
  // lowermost order, k = 0, have only 12 line (edge) dofs. The Nedelec elements
  // of the higher orders, k > 0, have 3*(k+1)*(k+2)^2 dofs in total if dim=3.
  // The dofs in a cell are distributed between lines (edges), quads (faces),
  // and the hex (the interior of the cell) as the following:
  //
  // 12*(k+1) line dofs; (k+1) dofs per line.
  // 2*6*k*(k+1) quad dofs; 2*k*(k+1) dofs per quad.
  // 3*(k+2)^2*(k+1) hex dofs;
  //
  // The dofs are indexed in the following order: first all line dofs,
  // then all quad dofs, and then all hex dofs.
  //
  // The line dofs need only sign adjustments. No permutation of line
  // dofs is needed. The line dofs are treated by
  // internal::FE_PolyTensor::get_dof_sign_change_nedelec(...)
  // in fe_poly_tensor.cc.
  //
  // The hex dofs need no adjustments: they are not shared between
  // neighbouring mesh cells.
  //
  // The two-dimensional Nedelec finite elements share no quad dofs between
  // neighbouring mesh cells. The zero-order three-dimensional Nedelec finite
  // elements have no quad dofs. Consequently, here we treat only quad dofs of
  // the three-dimensional Nedelec finite elements of the higher orders, k>0.
  // The questions how the curl looks like in the higher-dimensional spaces and
  // what does it mean to be curl-conforming if dim>3 we leave unanswered.
  //
  // In FE_Nedelec<dim>::initialize_quad_dof_index_permutation_and_sign_change()
  // we need to change some entries in the following two vectors of tables:
  // adjust_quad_dof_index_for_face_orientation_table
  // and
  // adjust_quad_dof_sign_for_face_orientation_table.
  // These tables specify the permutations and sign adjustments of the quad dofs
  // only. The tables are already filled with zeros meaning no permutations or
  // sign change are required. We need to change some entries of the tables such
  // that the shape functions that correspond to the quad dofs and are shared
  // between neighbouring cells have consistent orientations.
  //
  // The swap tables below describe the dof permutations and sign changes that
  // need to be done. The function
  // FE_Nedelec<dim>::initialize_quad_dof_index_permutation_and_sign_change()
  // simply reads the information in the swap tables below and puts it into
  // tables
  // adjust_quad_dof_index_for_face_orientation_table
  // and
  // adjust_quad_dof_sign_for_face_orientation_table.
  // A good question is: why don't we put the information into the tables of
  // deal.II right away? The answer is the following. The information on the
  // necessary dof permutations and sign changes is derived by plotting the
  // shape functions and observing them on faces of different orientations.
  // It is convenient to put the observations first in the format of the swap
  // tables below and then convert the swap tables into the format used by
  // deal.II.
  //
  // The dofs on a quad are indexed as the following:
  //
  // | x0, x1, x2, x3, ..., xk  | y0, y1, y2, y3 ..., yk  |
  // |                          |                         |
  // |<------ k*(k+1) --------->|<------ k*(k+1) -------->|
  // |                                                    |
  // |<------------------- 2*k*(k+1) -------------------->|
  //
  // Only one type of dof permutation is needed: swap between two dofs; one
  // dof being xi, another yj. That is, if x4 is replaced with y7,
  // then y7 must be replaced with x4. Such swaps can be ordered as
  // illustrated by the following example:
  //
  //                *
  //        y0, y9, y1, y0, ..., yk
  //      ---------------------------                     (swap)
  //        x0, x1, x2, x3, ..., xk
  //        *
  //
  // An x-dof below the line is swapped with the corresponding y-dof above the
  // line. A dof marked by the asterisk must change its sign before the swap.
  //
  // The x-dofs are assumed to have the normal order. There is no need to
  // encode it. Therefore, the swap tables need to encode the following
  // information: indices of the y-dofs, the sign change of the x-dofs, and
  // sign change of the y-dofs. The swap above is encoded as the following:
  //
  // swap = { 0, 9, 1, 0, ...., yk,  // indices of the y-dofs
  //          1, 0, 0, 0, ...., 0,   // sign change of the x-dofs,
  //          0, 0, 1, 0, ...., 0};  // sign change of the y-dofs.
  //
  // If no swap is needed, -1 is placed instead of the y-dof index.
  //
  // Such swaps are assembled into the swap table:
  //
  // swap_table = {swap_0, swap_1, ... swap_7};
  //
  // Each swap table contains eight swaps - one swap for each possible quad
  // orientation. The deal.II encodes the orientation of a quad using
  // three boolean parameters:
  // face_orientation - true if face is in standard orientation
  // and false otherwise;
  // face_rotation - rotation by 90 deg counterclockwise if true;
  // face_flip - rotation by 180 deg counterclockwise if true.
  // See the documentation of GeometryInfo<dim>.
  //
  // The combined face orientation is computes as
  // orientation_no = face_flip*4 + face_rotation*2 + face_orientation*1;
  // See tria_orientation.h.
  //
  // The parameter orientation_no (0...7) indexes the swaps in a swap table.
  //
  // Nedelec elements of order k have their own swap table, swap_table_k.
  // Recall, the swap_table_0 is empty as the Nedelec finite elements of the
  // lowermost order have no quad dofs.

  const std::vector<std::vector<std::vector<int>>> swap_table_1 = {
    // 0   1
    {{0, 1}, // 0
     {0, 0},
     {0, 0}},
    {{-1, -1}, // 1
     {0, 0},
     {0, 0}},
    {{-1, -1}, // 2
     {0, 0},
     {1, 0}},
    {{0, 1}, // 3
     {1, 0},
     {0, 0}},
    {{0, 1}, // 4
     {1, 0},
     {1, 0}},
    {{-1, -1}, // 5
     {1, 0},
     {1, 0}},
    {{-1, -1}, // 6
     {1, 0},
     {0, 0}},
    {{0, 1}, // 7
     {0, 0},
     {1, 0}}};

  const std::vector<std::vector<std::vector<int>>> swap_table_2 = {
    // 0   1   2   3   4   5
    {{0, 3, 1, 4, 2, 5}, // 0
     {0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 0}},
    {{-1, -1, -1, -1, -1, -1}, // 1
     {0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 0}},
    {{-1, -1, -1, -1, -1, -1}, // 2
     {0, 1, 0, 1, 0, 1},
     {1, 0, 1, 1, 0, 1}},
    {{0, 3, 1, 4, 2, 5}, // 3
     {1, 1, 0, 0, 1, 1},
     {0, 0, 0, 1, 1, 1}},
    {{0, 3, 1, 4, 2, 5}, // 4
     {1, 0, 0, 1, 1, 0},
     {1, 0, 1, 0, 1, 0}},
    {{-1, -1, -1, -1, -1, -1}, // 5
     {1, 0, 0, 1, 1, 0},
     {1, 0, 1, 0, 1, 0}},
    {{-1, -1, -1, -1, -1, -1}, // 6
     {1, 1, 0, 0, 1, 1},
     {0, 0, 0, 1, 1, 1}},
    {{0, 3, 1, 4, 2, 5}, // 7
     {0, 1, 0, 1, 0, 1},
     {1, 0, 1, 1, 0, 1}}};

  const std::vector<std::vector<std::vector<int>>> swap_table_3 = {
    // 0   1   2   3   4   5   6   7   8   9  10  11
    {{0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11}, // 0
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 2
     {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0},
     {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0}},
    {{0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11}, // 3
     {1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0},
     {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0}},
    {{0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11}, // 4
     {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
     {1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 5
     {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
     {1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 6
     {1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0},
     {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0}},
    {{0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11}, // 7
     {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0},
     {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0}}};

  const std::vector<std::vector<std::vector<int>>> swap_table_4 = {
    // Swap sign_X and sign_Y rows if k=4. Why?...
    // 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    // 19
    {{0,  5,  10, 15, 1,  6,  11, 16, 2,  7,
      12, 17, 3,  8,  13, 18, 4,  9,  14, 19}, // 0
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 2
     {1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1},
     {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}},
    {{0,  5,  10, 15, 1,  6,  11, 16, 2,  7,
      12, 17, 3,  8,  13, 18, 4,  9,  14, 19}, // 3
     {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1}},
    {{0,  5,  10, 15, 1,  6,  11, 16, 2,  7,
      12, 17, 3,  8,  13, 18, 4,  9,  14, 19}, // 4
     {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
     {1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 5
     {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
     {1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0}},
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 6
     {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1}},
    {{0,  5,  10, 15, 1,  6,  11, 16, 2,  7,
      12, 17, 3,  8,  13, 18, 4,  9,  14, 19}, // 7
     {1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1},
     {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}}};

  // Finally, the swap tables are organized into a vector.

  const std::vector<std::vector<std::vector<std::vector<int>>>>
    swap_tables_vector = {
      {{}}, // There are no quad dofs in FE_Nedelec of the lowermost order.
      {swap_table_1},
      {swap_table_2},
      {swap_table_3},
      {swap_table_4}};

  // TODO: understand why the lines sign_X and sign_Y must be swapped in the
  // case k=4. After that, formulate the tendencies in the swap tables as
  // closed-form expressions.
};

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
FE_Nedelec<1>::initialize_restriction();

#endif // DOXYGEN

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
