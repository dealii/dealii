// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_nedelec_sz_h
#define dealii_fe_nedelec_sz_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/polynomials_integrated_legendre_sz.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup fe
 * @{
 */

/**
 * This class represents an implementation of the
 * H<sup>curl</sup>-conforming N&eacute;d&eacute;lec element described in the
 * PhD thesis @cite Zaglmayr2006. It its used in the same context as
 * described at the top of the description for the FE_Nedelec class.
 *
 * This element overcomes the sign conflict issues present in
 * traditional N&eacute;d&eacute;lec elements that arise from the edge
 * and face parameterizations used in the basis functions. Therefore,
 * this element should provide consistent results for general
 * quadrilateral and hexahedral elements for which the relative
 * orientations of edges and faces as seen from all adjacent cells are
 * often difficult to establish.
 *
 * The way this element addresses the sign conflict problem is to
 * assign local edges and faces a globally defined orientation. The
 * local edge orientation is always chosen such that the first vertex
 * defining the edge is the one that has the highest global vertex
 * numbering, with the second edge vertex being that which has the
 * lowest global vertex numbering.
 *
 * Similarly, the face orientation is always chosen such that the first
 * vertex is chosen to be that with the highest global vertex numbering of the
 * four vertices making up the face. The third vertex is then chosen to be that
 * which is geometrically opposite the first vertex, and the second and fourth
 * vertices are decided such that the second has a higher global vertex
 * numbering than the fourth.
 *
 * To overcome the sign conflict in the case of non-conforming meshes, hanging
 * edges (and faces in 3D) must be addressed. This element chooses the
 * orientation of hanging edges (and faces) based on the parent cell. Moreover,
 * the associated constraint matrix that restricts the additional DoFs
 * originating from the hanging edges and faces must be adapted to account for
 * the orientation of the edges and faces.
 *
 * Further details on this element, including some benchmarking, can be found
 * in the paper @cite Kynch2017. For details on the implementation of the
 * hanging node constraints, see the paper @cite Kinnewig2023.
 */
template <int dim, int spacedim = dim>
class FE_NedelecSZ : public FiniteElement<dim, dim>
{
public:
  static_assert(dim == spacedim,
                "FE_NedelecSZ is only implemented for dim==spacedim!");

  /**
   * Constructor for the NedelecSZ element of given @p order. The maximal
   * polynomial degree of the shape functions is `order+1` (in each variable;
   * the total polynomial degree may be higher). If `order = 0`, the element is
   * linear and has degrees of freedom only on the edges. If `order >= 1` the
   * element has degrees of freedom on the edges, faces and volume. For example
   * the 3d version of FE_NedelecSZ has 12 degrees of freedom for `order = 0`
   * and 54 for `degree = 1`. It is important to have enough quadrature points
   * in order to perform the quadrature with sufficient accuracy.
   * For example
   * [QGauss<dim>(order + 2)](@ref QGauss)
   * can be used for the
   * quadrature formula, where `order` is the order of FE_NedelecSZ.
   */
  FE_NedelecSZ(const unsigned int order);

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  /**
   * Return the value of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   *
   * Since this element is vector-valued, this function will
   * throw an exception.
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the value of the <tt>component</tt>th vector component of the
   * <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim>  &p,
                        const unsigned int component) const override;

  /**
   * Return the gradient of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   *
   * Since this element is vector-valued, this function will
   * throw an exception.
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the gradient of the <tt>component</tt>th vector component of the
   * <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim>  &p,
                       const unsigned int component) const override;

  /**
   * Return the tensor of second derivatives of the <tt>i</tt>th shape
   * function at point <tt>p</tt> on the unit cell. See the FiniteElement base
   * class for more information about the semantics of this function.
   *
   * Since this element is vector-valued, this function will
   * throw an exception.
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the second derivative of the <tt>component</tt>th vector component
   * of the <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim>  &p,
                            const unsigned int component) const override;

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

protected:
  /**
   * The mapping kind to be used to map shape functions from the reference
   * cell to the mesh cell.
   */
  MappingKind mapping_kind;

  /**
   * Compute the value and the derivatives of the Nedelec functions at
   * the points given in <tt>p_list</tt>.
   */
  void
  evaluate(const std::vector<Point<dim>> &p_list,
           const UpdateFlags              update_flags,
           std::unique_ptr<
             typename dealii::FiniteElement<dim, spacedim>::InternalDataBase>
             &data_ptr) const;

  virtual std::unique_ptr<
    typename dealii::FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim>        &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * Compute information about the shape functions on the cell denoted by the
   * first argument. Note that this function must recompute the cell-dependent
   * degrees of freedom, and so is not thread-safe at this time.
   */
  virtual void
  fill_fe_values(
    const typename Triangulation<dim, dim>::cell_iterator &cell,
    const CellSimilarity::Similarity                       cell_similarity,
    const Quadrature<dim>                                 &quadrature,
    const Mapping<dim, dim>                               &mapping,
    const typename Mapping<dim, dim>::InternalDataBase    &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, dim>
                                                             &mapping_data,
    const typename FiniteElement<dim, dim>::InternalDataBase &fedata,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
      &data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  /**
   * Compute information about the shape functions on the cell and face denoted
   * by the first two arguments. Note that this function must recompute the
   * cell-dependent degrees of freedom, and so is not thread-safe at this time.
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, dim>::cell_iterator &cell,
    const unsigned int                                     face_no,
    const hp::QCollection<dim - 1>                        &quadrature,
    const Mapping<dim, dim>                               &mapping,
    const typename Mapping<dim, dim>::InternalDataBase    &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, dim>
                                                             &mapping_data,
    const typename FiniteElement<dim, dim>::InternalDataBase &fedata,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
      &data) const override;

  /**
   * Not implemented.
   */
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, dim>::cell_iterator &cell,
    const unsigned int                                     face_no,
    const unsigned int                                     sub_no,
    const Quadrature<dim - 1>                             &quadrature,
    const Mapping<dim, dim>                               &mapping,
    const typename Mapping<dim, dim>::InternalDataBase    &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, dim>
                                                             &mapping_data,
    const typename FiniteElement<dim, dim>::InternalDataBase &fedata,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
      &data) const override;

  /**
   * Derived Internal data which is used to store cell-independent data.
   * Note that due to the nature of this element, a number of useful
   * pre-computed quantities are stored for the computation of cell-dependent
   * shape functions.
   *
   * The main quantities which are stored are associated with edge and face
   * parameterizations. These are:
   * <ul>
   * <li> $\lambda_{i}$ - trilinear function, equal to one at the $i$-th vertex
   * and zero at all other vertices.</li>
   * <li> $\sigma_{i}$ - linear functional associated with the $i$-th vertex.</li>
   * </ul>
   *
   * The definitions of these functionals, as well as the edge and face
   * parameterizations and edge and face extension parameters, can be found on
   * page 82 of Zaglmayr's thesis. The details of the definition of the
   * globally-defined edge and face orientations can be found on page 67.
   */
  class InternalData : public FiniteElement<dim, dim>::InternalDataBase
  {
  public:
    /**
     * Storage for shape functions on the reference element. We only pre-compute
     * cell-based DoFs, as the edge- and face-based DoFs depend on the cell.
     *
     * Due to the cell-dependent DoFs, this variable is declared mutable.
     */
    mutable std::vector<std::vector<Tensor<1, dim>>> shape_values;

    /**
     * Storage for shape function gradients on the reference element. We only
     * pre-compute cell-based DoFs, as the edge- and face-based DoFs depend on
     * the cell.
     *
     * Due to the cell-dependent DoFs, this variable is declared mutable.
     */
    mutable std::vector<std::vector<DerivativeForm<1, dim, dim>>> shape_grads;

    /**
     * Storage for shape function hessians on the reference element. We only
     * pre-compute cell-based DoFs, as the edge- and face-based DoFs depend on
     * the cell.
     *
     * Due to the cell-dependent DoFs, this variable is declared mutable.
     */
    mutable std::vector<std::vector<DerivativeForm<2, dim, dim>>>
      shape_hessians;

    /**
     * Storage for all possible edge parameterization between vertices. These
     * are required in the computation of edge- and face-based DoFs, which are
     * cell-dependent.
     *
     * The edge parameterization of an edge, E, starting at vertex i and ending
     * at vertex $j$ is given by $\sigma_{E} = \sigma_{i} - \sigma{j}$.
     *
     * sigma_imj_values[q][i][j] stores the value of the edge parameterization
     * connected by vertices $i$ and $j$ at the q-th quadrature point.
     *
     * Note that not all of the $i$ and $j$ combinations result in valid edges
     * on the hexahedral cell, but they are computed in this fashion for use
     * with non-standard edge and face orientations.
     */
    std::vector<std::vector<std::vector<double>>> sigma_imj_values;

    /**
     * Storage for gradients of all possible edge parameterizations between
     * vertices. These are required in the computation of edge- and face-based
     * DoFs, which are cell-dependent. Note that the components of the gradient
     * are constant.
     *
     * The edge parameterization of an edge, $E$, starting at vertex $i$ and
     * ending at vertex $j$ is given by $\sigma_{E} = \sigma_{i} - \sigma{j}$.
     *
     * sigma_imj_grads[i][j][d] stores the gradient of the edge parameterization
     * connected by vertices $i$ and $j$ in component $d$.
     *
     * Note that the gradient of the edge parameterization is constant on an
     * edge, so we do not need to store it at every quadrature point.
     */
    std::vector<std::vector<std::vector<double>>> sigma_imj_grads;

    /**
     * Storage for values of edge parameterizations at quadrature points. These
     * are stored for the 12 edges such that the global vertex numbering would
     * follow the order defined by the "standard" deal.II cell.
     *
     * edge_sigma_values[m][q] stores the edge parameterization value at the
     * q-th quadrature point on edge m.
     *
     * These values change with the orientation of the edges of a physical cell
     * and so must take the "sign" into account when used for computation.
     */
    std::vector<std::vector<double>> edge_sigma_values;

    /**
     * Storage for gradients of edge parameterization at quadrature points.
     * These are stored for the 12 edges such that the global vertex numbering
     * would follow the order defined by the "standard" deal.II cell.
     *
     * edge_sigma_grads[m][d] stores the gradient of the edge parameterization
     * for component d on edge m.
     *
     * These values change with the orientation of the edges of a physical cell
     * and so must take the "sign" into account when used for computation.
     */
    std::vector<std::vector<double>> edge_sigma_grads;

    /**
     * Storage for edge extension parameters at quadrature points. These are
     * stored for the 12 edges such that the global vertex numbering would
     * follow the order defined by the "standard" deal.II cell.
     *
     * The edge extension parameter of an edge, $E$, starting at vertex $i$ and
     * ending at vertex $j$ is given by $\lambda_{E} = \lambda_{i} +
     * \lambda_{j}$.
     *
     * Note that under this definition, the values of $\lambda_{E}$ do not
     * change with the orientation of the edge.
     *
     * edge_lambda_values[m][q] stores the edge extension parameter value at
     * the $q$-th quadrature point on edge $m$.
     */
    std::vector<std::vector<double>> edge_lambda_values;

    /**
     * Storage for gradients of edge extension parameters in 2d. In this case
     * they are constant. These are stored for the 12 edges such that the global
     * vertex numbering* would follow the order defined by the "standard"
     * deal.II cell.
     *
     * edge_lambda_grads_2d[m][d] stores the gradient of the edge extension
     * parameter for component $d$ on edge $m$.
     */
    std::vector<std::vector<double>> edge_lambda_grads_2d;

    /**
     * Storage for gradients of edge extension parameters in 3d. In this case
     * they are non-constant. These are stored for the 12 edges such that the
     * global vertex numbering* would follow the order defined by the
     * "standard" deal.II cell.
     *
     * edge_lambda_grads_3d[m][q][d] stores the gradient of the edge extension
     * parameter for component $d$ at the $q$-th quadrature point on edge m.
     */
    std::vector<std::vector<std::vector<double>>> edge_lambda_grads_3d;

    /**
     * Storage for 2nd derivatives of edge extension parameters in 3d, which are
     * constant across the cell. These are stored for the 12 edges such that the
     * global vertex numbering* would follow the order defined by the
     * "standard" deal.II cell.
     *
     * edge_lambda_gradgrads_3d[m][d1][d2] stores the 2nd derivatives of the
     * edge extension parameters with respect to components d1 and d2 on edge
     * $m$.
     */
    std::vector<std::vector<std::vector<double>>> edge_lambda_gradgrads_3d;

    /**
     * Storage for the face extension parameters. These are stored for the 6
     * faces such that the global vertex numbering would follow the order
     * defined by the "standard" deal.II cell.
     *
     * The face extension parameter of a face, F, defined by the vertices
     * v1, v2, v3, v4 is given by
     * $\lambda_{F} = \lambda_{v1} + \lambda_{v2} + \lambda_{v3} +
     * \lambda_{v4}$.
     *
     * Note that under this definition, the values of $\lambda_{F}$ do not
     * change with the orientation of the face.
     *
     * face_lambda_values[m][q] stores the face extension parameter value at
     * the $q$-th quadrature point on face $m$.
     */
    std::vector<std::vector<double>> face_lambda_values;

    /**
     * Storage for gradients of face extension parameters. These are stored for
     * the 6 faces such that the global vertex numbering would follow the order
     * defined by the "standard" deal.II cell.
     *
     * face_lambda_grads[m][d] stores the gradient of the face extension
     * parameters for component $d$ on face $m$.
     */
    std::vector<std::vector<double>> face_lambda_grads;
  };

private:
  /**
   * Internal function to return a vector of "dofs per object"
   * where the components of the returned vector refer to:
   * 0 = vertex
   * 1 = edge
   * 2 = face (which is a cell in 2d)
   * 3 = cell
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * Internal storage for all required integrated Legendre polynomials.
   */
  std::vector<Polynomials::Polynomial<double>> IntegratedLegendrePolynomials;

  /**
   * Internal function to populate the internal array of integrated Legendre
   * polynomials.
   */
  void
  create_polynomials(const unsigned int degree);

  /**
   * Returns the number of DoFs in the basis set.
   */
  unsigned int
  compute_num_dofs(const unsigned int degree) const;

  /**
   * Populates cell-dependent edge-based shape functions on the given
   * InternalData object.
   */
  void
  fill_edge_values(const typename Triangulation<dim, dim>::cell_iterator &cell,
                   const Quadrature<dim> &quadrature,
                   const InternalData    &fedata) const;

  /**
   * Populates the cell-dependent face-based shape functions on the given
   * InternalData object.
   */
  void
  fill_face_values(const typename Triangulation<dim, dim>::cell_iterator &cell,
                   const Quadrature<dim> &quadrature,
                   const InternalData    &fedata) const;

  /**
   * Mutex variables used for protecting the initialization of restriction
   * and embedding matrices.
   */
  mutable Threads::Mutex prolongation_matrix_mutex;
};



/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
