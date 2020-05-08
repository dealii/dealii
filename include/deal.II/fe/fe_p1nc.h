// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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

#ifndef dealii_fe_p1nc_h
#define dealii_fe_p1nc_h

#include <deal.II/base/config.h>

#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of the scalar version of the P1 nonconforming finite
 * element, a piecewise linear element on quadrilaterals in 2D.
 * This implementation is only for 2D cells in a 2D space (i.e., codimension 0).
 *
 * Unlike the usual continuous, $H^1$ conforming finite elements,
 * the P1 nonconforming element does not enforce continuity across edges.
 * However, it requires the continuity in an integral sense:
 * any function in the space should have the same integral values
 * on two sides of the common edge shared by two adjacent elements.
 *
 * Thus, each function in the nonconforming element space can be
 * discontinuous, and consequently not included in $H^1_0$, just like
 * the basis functions in Discontinuous Galerkin (DG) finite element
 * spaces. On the other hand, basis functions in DG spaces are
 * completely discontinuous across edges without any relation between
 * the values from both sides.  This is a reason why usual weak
 * formulations for DG schemes contain additional penalty terms for
 * jump across edges to control discontinuity.  However, nonconforming
 * elements usually do not need additional terms in their weak
 * formulations because their integrals along edges are the same from
 * both sides, i.e., there is <i>some level</i> of continuity.
 *
 * <h3>Dice Rule</h3>
 * Since any function in the P1 nonconforming space is piecewise linear on each
 * element, the function value at the midpoint of each edge is same as the mean
 * value on the edge. Thus the continuity of the integral value across each edge
 * is equivalent to the continuity of the midpoint value of each edge in this
 * case.
 *
 * Thus for the P1 nonconforming element, the function values at midpoints on
 * edges of a cell are important. The first attempt to define (local) degrees of
 * freedom (DoFs) on a quadrilateral is by using midpoint values of a function.
 *
 * However, these 4 functionals are not linearly independent
 * because a linear function on 2D is uniquely determined by only 3 independent
 * values. A simple observation reads that any linear function on a
 * quadrilateral should satisfy the 'dice rule': the sum of two function values
 * at the midpoints of the edge pair on opposite sides of a cell is equal to the
 * sum of those at the midpoints of the other edge pair. This is called the
 * 'dice rule' because the number of points on opposite sides of a dice always
 * adds up to the same number as well (in the case of dice, to seven).
 *
 * In formulas, the dice rule is written as $\phi(m_0) + \phi(m_1) = \phi(m_2) +
 * \phi(m_3)$ for all $\phi$ in the function space where $m_j$ is the midpoint
 * of the edge $e_j$. Here, we assume the standard numbering convention for
 * edges used in deal.II and described in class GeometryInfo.
 *
 * Conversely if 4 values at midpoints satisfying the dice rule are given,
 * then there always exists the unique linear function which coincides with 4
 * midpoints values.
 *
 * Due to the dice rule, three values at any three midpoints can determine
 * the last value at the last midpoint.
 * It means that the number of independent local functionals on a cell is 3,
 * and this is also the dimension of the linear polynomial space on a cell in
 * 2D.
 *
 * <h3>Shape functions</h3>
 * Before introducing the degrees of freedom, we present 4 local shape functions
 * on a cell. Due to the dice rule, we need a special construction for shape
 * functions. Although the following 4 shape functions are not linearly
 * independent within a cell, they are helpful to define the global basis
 * functions which are linearly independent on the whole domain. Again, we
 * assume the standard numbering for vertices used in deal.II.
 *
 *  @verbatim
 *  2---------|---------3
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 *  -                   -
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 *  0---------|---------1
 *  @endverbatim
 *
 * For each vertex $v_j$ of given cell, there are two edges of which $v_j$ is
 * one of end points. Consider a linear function such that it has value 0.5 at
 * the midpoints of two adjacent edges, and 0.0 at the two midpoints of the
 * other edges. Note that the set of these values satisfies the dice rule which
 * is described above. We denote such a function associated with vertex $v_j$ by
 * $\phi_j$. Then the set of 4 shape functions is a partition of unity on a
 * cell: $\sum_{j=0}^{3} \phi_j = 1$. (This is easy to see: at each edge
 * midpoint, the sum of the four function adds up to one because two functions
 * have value 0.5 and the other value 0.0. Because the function is globally
 * linear, the only function that can have value 1 at four points must also be
 * globally equal to one.)
 *
 * The following figures represent $\phi_j$ for $j=0,\cdots,3$ with their
 * midpoint values:
 *
 * <ul>
 * <li> shape function $\phi_0$:
 *  @verbatim
 *  +--------0.0--------+
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 * 0.5                 0.0
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 *  +--------0.5--------+
 *  @endverbatim
 *
 * <li> shape function $\phi_1$:
 *  @verbatim
 *  +--------0.0--------+
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 * 0.0                 0.5
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 *  +--------0.5--------+
 *  @endverbatim
 *
 * <li> shape function $\phi_2$:
 *  @verbatim
 *  +--------0.5--------+
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 * 0.5                 0.0
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 *  +--------0.0--------+
 *  @endverbatim
 *
 * <li> shape function $\phi_3$:
 *  @verbatim
 *  +--------0.5--------+
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 * 0.0                 0.5
 *  |                   |
 *  |                   |
 *  |                   |
 *  |                   |
 *  +--------0.0--------+
 *  @endverbatim
 *
 * </ul>
 *
 * The local DoFs are defined by the coefficients of the shape functions
 * associated with vertices, respectively. Although these 4 local DoFs are not
 * linearly independent within a single cell, this definition is a good start
 * point for the definition of the global DoFs.
 *
 * We want to emphasize that the shape functions are constructed on each cell,
 * not on the reference cell only. Usual finite elements are defined based on a
 * 'parametric' concept. It means that a function space for a finite element is
 * defined on one reference cell, and it is transformed into each cell via a
 * mapping from the reference cell. However the P1 nonconforming element does
 * not follow such concept. It defines a function space with linear shape
 * functions on each cell without any help of a function space on the reference
 * cell. In other words, the element is defined in real space, not via a mapping
 * from a reference cell. In this, it is similar to the FE_DGPNonparametric
 * element.
 *
 * Thus this implementation does not have to compute shape values on the
 * reference cell. Rather, the shape values are computed by construction of the
 * shape functions on each cell independently.
 *
 * <h3>Degrees of freedom</h3>
 * We next have to consider the <i>global</i> basis functions for the element
 * because the system of equations which we ultimately have to solve is for a
 * global system, not local. The global basis functions associated with a node
 * are defined by a cell-wise composition of local shape functions associated
 * with the node on each element.
 *
 * There is a theoretical result about the linear independence of the global
 * basis functions depending on the type of the boundary condition we consider.
 *
 * When homogeneous Dirichlet boundary conditions are given,
 * the global basis functions associated with interior nodes are linearly
 * independent. Then, the number of DoFs is equal to the number of interior
 * nodes, and consequently the same as the number of DoFs for the standard
 * bilinear $Q_1$ finite element.
 *
 * When Neumann boundary conditions are given,
 * the global basis functions associated with all nodes (including boundary
 * nodes) are actually not linearly independent. There exists one redundancy.
 * Thus in this case, the number of DoFs is equal to the number of all nodes
 * minus 1. This is, again as for the regular $Q_1$ element.
 *
 * <h3>Unit support points</h3>
 * For a smooth function, we construct a piecewise linear function which belongs
 * to the element space by using its nodal values as DoF values.
 *
 * Note that for the P1 nonconforming element, two nodal values of a smooth
 * function and its interpolant do not coincide in general, in contrast to
 * ordinary Lagrange finite elements. Of course, it is meaningless to refer
 * 'nodal value' because the element space has nonconformity. But it is also
 * true even though the single global basis function associated with a node is
 * considered the unique 'nodal value' at the node. For instance, consider the
 * basis function associated with a node. Consider two lines representing the
 * level sets for value 0.5 and 0, respectively, by connecting two midpoints.
 * Then we cut the quad into two sub-triangles by the diagonal which is placed
 * along those two lines. It gives another level set for value 0.25 which
 * coincides with the cutting diagonal. Therefore these three level sets are all
 * parallel in the quad and it gives the value 0.75 at the base node, not
 * value 1. This is true whether the quadrilateral is a rectangle,
 * parallelogram, or any other shape.
 *
 * <h3>Reference</h3>
 * The original paper for the P1 nonconforming element  by Park and Sheen
 * is accessible at https://doi.org/10.1137/S0036142902404923 ,
 * see @cite park2003p .
 *
 * @author Jaeryun Yim, 2015, 2016.
 */
class FE_P1NC : public FiniteElement<2, 2>
{
public:
  /**
   * Constructor for the P1 nonconforming element.
   * It is only for 2D and codimension = 0.
   */
  FE_P1NC();

  virtual std::string
  get_name() const override;

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags flags) const override;

  virtual std::unique_ptr<FiniteElement<2, 2>>
  clone() const override;

  /**
   * Destructor.
   */
  virtual ~FE_P1NC() override = default;



private:
  /**
   * Return the vector consists of the numbers of degrees of freedom per
   * objects.
   */
  static std::vector<unsigned int>
  get_dpo_vector();

  /**
   * Return the coefficients of 4 local linear shape functions $\phi_j(x,y) = a
   * x + b y + c$ on given cell. For each local shape function, the array
   * consists of three coefficients is in order of a,b and c.
   */
  static std::array<std::array<double, 3>, 4>
  get_linear_shape_coefficients(const Triangulation<2, 2>::cell_iterator &cell);

  /**
   * Do the work which is needed before cellwise data computation.
   * Since the shape functions are constructed independently on each cell,
   * the data on the reference cell is not necessary.
   * It returns an empty variable type of @ InternalDataBase and updates @
   * update_flags, and computes trivially zero Hessian for each cell if it is
   * needed.
   */
  virtual std::unique_ptr<FiniteElement<2, 2>::InternalDataBase>
  get_data(
    const UpdateFlags update_flags,
    const Mapping<2, 2> &,
    const Quadrature<2> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  virtual std::unique_ptr<FiniteElement<2, 2>::InternalDataBase>
  get_face_data(
    const UpdateFlags update_flags,
    const Mapping<2, 2> &,
    const Quadrature<1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  virtual std::unique_ptr<FiniteElement<2, 2>::InternalDataBase>
  get_subface_data(
    const UpdateFlags update_flags,
    const Mapping<2, 2> &,
    const Quadrature<1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * Compute the data on the current cell.
   */
  virtual void
  fill_fe_values(
    const Triangulation<2, 2>::cell_iterator &cell,
    const CellSimilarity::Similarity          cell_similarity,
    const Quadrature<2> &                     quadrature,
    const Mapping<2, 2> &                     mapping,
    const Mapping<2, 2>::InternalDataBase &   mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<2, 2>
      &                                          mapping_data,
    const FiniteElement<2, 2>::InternalDataBase &fe_internal,
    internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * Compute the data on the face of the current cell.
   */
  virtual void
  fill_fe_face_values(
    const Triangulation<2, 2>::cell_iterator &cell,
    const unsigned int                        face_no,
    const Quadrature<1> &                     quadrature,
    const Mapping<2, 2> &                     mapping,
    const Mapping<2, 2>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<2, 2>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * Compute the data on the subface of the current cell.
   */
  virtual void
  fill_fe_subface_values(
    const Triangulation<2, 2>::cell_iterator &cell,
    const unsigned int                        face_no,
    const unsigned int                        sub_no,
    const Quadrature<1> &                     quadrature,
    const Mapping<2, 2> &                     mapping,
    const Mapping<2, 2>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<2, 2>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * Create the constraints matrix for hanging edges.
   */
  void
  initialize_constraints();
};



/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
