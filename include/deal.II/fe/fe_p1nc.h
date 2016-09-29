// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2016 by the deal.II authors
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

#ifndef dealii__fe_p1nc_h
#define dealii__fe_p1nc_h

#include <deal.II/base/config.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation for the scalar version of the P1 nonconforming finite
 * element, a piecewise linear element on quadrilaterals in 2D.
 * This implementation is only for 2D and codimension = 0.
 *
 * Unlike any continuous conforming finite element which belongs to $H^1_0$,
 * the P1 nonconforming element does not enforce the continuity across edges.
 * But it requires the continuity just in integral sense:
 * any function in the space should have the same integral values
 * on two sides of the common edge shared by two adjacent elements.
 *
 * Thus each function in the nonconforming element space can be discontinuous,
 * not included in $H^1_0$, as functions in Discontinuous Galerkin (DG) finite
 * element spaces.
 * Although any function in DG space also has nonconformity,
 * it is completely discontinuous across edges without any relation.
 * This is a reason why usual weak formulations for DG schemes contain
 * additional penalty terms for jump across edges to control discontinuity.
 * However nonconforming elements usually do not need additional terms
 * in their weak formulations due to the continuity in integral on edges.
 *
 * <h3>Dice Rule</h3>
 * Since any function in the P1 nonconforming space is piecewise linear on each element,
 * the function value at the midpoint of each edge is same as the mean value on the edge.
 * Thus the continuity of the integral value across each edge is equivalent to
 * the continuity of the midpoint value of each edge in this case.
 *
 * Thus for the P1 nonconforming element, the function values at midpoints on edges of a cell are important.
 * The first attempt to define (local) degrees of freedom (DOFs) on a quadrilateral
 * is by using midpoint values of a function as usual nonconforming finite elements.
 *
 * However, these 4 functionals are not linearly independent
 * because a linear function on 2D is uniquely determined by only 3 independent values.
 * A simple observation reads that any linear function on a quadrilateral should satisfies the 'dice rule':
 * the sum of two function values at two midpoints of the edge pair on opposite
 * position is equal to the sum of those of the other edge pair.
 * This is called the 'dice rule' because the number of points on opposite sides of a dice always
 * adds up to the same number as well (in the case of dice, to seven).
 *
 * In formulas, the dice rule is written as $\phi(m_0) + \phi(m_1) = \phi(m_2) + \phi(m_3)$
 * for all $\phi$ in the function space where $m_j$ is the midpoint of the edge $e_j$.
 * Here, we assume the standard numbering convention for edges used in deal.II
 * and described for class GeometryInfo.
 *
 * Conversely if 4 values at midpoints satisfying the dice rule are just given,
 * then there always exists the unique linear function which coincides with 4 midpoints values.
 *
 * Due to the dice rule, three values at any three midpoints can determine
 * the last value at the last midpoint.
 * It means that the number of independent local functionals on a cell is 3,
 * and it is same as the dimension of the linear polynomial space on a cell in 2D.
 *
 * <h3>Shape functions</h3>
 * Before introducing the DOFs, we present 4 local shape functions on a cell.
 * Due to the dice rule, we need a special construction for shape functions.
 * Although the following 4 shape functions are not linearly independent within a cell,
 * they are helpful to define the global basis functions which are linearly independent on whole domain.
 * Again, we assume the standard numbering for vertices used in deal.II.
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
 * For each vertex $v_j$ of given cell, there are two edges of which $v_j$ is one of end points.
 * Consider a linear function such that 0.5 value at two midpoints of such edges,
 * and 0.0 at two midpoints of other edges.
 * Note that the set of these values satisfies the dice rule which is described above.
 * We denote such a function associated with vertex $v_j$ by $\phi_j$.
 * Then the set of 4 shape functions is a partition of unity on a cell: $\sum_{j=0}^{3} \phi_j = 1$.
 *
 * The following figures represent $\phi_j$, $j=0,\cdots,3$ with its values at midpoints.
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
 * The local DOFs are defined by the coefficients of the shape functions associated vertices, respectively.
 * Although these 4 local DOFs are not linearly independent within a single cell as well,
 * this definition is a good start point for the definition of the global DOFs.
 *
 * We want to emphasize that the shape functions are constructed on each cell, not on the reference cell only.
 * Usual finite elements are defined based on a 'parametric' concept.
 * That means, a function space for a finite element is defined on one reference cell, and it is transfomed
 * into each cell via a mapping from the reference cell.
 * However the P1 nonconforming element does not follow such concept. It defines a function space with
 * linear shape functions on each cell without any help of a function space on the reference cell.
 * In other words, the element is defined in real space, not via a mapping from a reference cell.
 *
 * Thus this implementation does not have to compute shape values on the reference cell.
 * Rather than, the shape values are computed by construction of the shape functions
 * on each cell independently.
 *
 * <h3>DOFs</h3>
 * We have to consider the basis function for the element space in global domain
 * because the system of equations which we have to solve at last is for a global system, not local.
 * The global basis function associated with a node is defined by a cell-wise composition of
 * local shape functions associated with the node on each element.
 * And we define a global DOF associated with a node by a coefficient of the basis function associated with that node.
 *
 * There is a theoretical result about the linear independency of the global basis functions
 * depending on the type of the boundary condition we consider.
 *
 * When the homogeneous Dirichlet boundary condition is given,
 * the global basis functions associated with interior nodes are linearly independent.
 * And the number of DOFs is equal to the number of interior nodes,
 * same as the number of DOFs for the standard bilinear finite element @p Q_1.
 *
 * When the Neumann boundary condition is given,
 * the global basis functions associated with all nodes (including boundary nodes)
 * are actually not linearly independent. There exists 1 redundancy.
 * Thus in this case, the number of DOFs is equal to the number of all nodes minus 1.
 *
 * <h3>Unit support points</h3>
 * For a smooth function, we construct a piecewise linear function which belongs to the element space by
 * using its nodal values as DOF values.
 * This interpolation is implemented by using appropriate @p unit_support_points.
 *
 * Note that two nodal values of a smooth function and its interpolant do not coincide in general,
 * contrast with ordinary Lagrange finite elements.
 * Of course, it is meaningless to refer 'nodal value' because the element space has nonconformity.
 * But it is also true even though the single global basis function associated a node is considered
 * with a valid expression 'nodal value'.
 * For instance, consider the basis function associated with a node.
 * Consider two lines representing the level sets for value 0.5 and 0, respectively, by connecting two midpoints.
 * Then we cut the quad into two sub-triangles by the diagonal which is placed along those two lines.
 * It gives another level set for value 0.25 which coincides with the cutting diagonal.
 * Therefore these three level sets are all parallel and it gives the value 0.75 at the base node, not value 1.
 * Even though a general quad is given, this is also true.
 *
 * <h3>References</h3>
 * You can find the paper about the P1NC element
 * Park & Sheen (2003). P1-nonconforming quadrilateral finite element methods for second-order elliptic problems.
 * SIAM Journal on Numerical Analysis, 41(2), 624-640,
 * available at http://epubs.siam.org/doi/abs/10.1137/S0036142902404923.
 *
 */

class FE_P1NC : public FiniteElement<2,2>
{

public:
  /**
   * Constructor for the P1 nonconforming element.
   * It is only for 2D and codimension = 0.
   */
  FE_P1NC() ;

  virtual std::string get_name () const ;

  virtual UpdateFlags requires_update_flags (const UpdateFlags flags) const ;

  virtual FiniteElement<2,2> *clone () const ;

  /**
   * Destructor.
   */
  virtual ~FE_P1NC ();



private:

  /**
   * Return the vector consists of the numbers of degrees of freedom per objects.
   */
  static std::vector<unsigned int> get_dpo_vector ();

  /**
   * Return the coefficients of 4 local linear shape functions $\phi_j(x,y) = a x + b y + c$ on given cell.
   * For each local shape function, the array consists of three coefficients is in order of a,b and c.
   */
  static std_cxx11::array<std_cxx11::array<double,3>,4>
  get_linear_shape_coefficients (const Triangulation<2,2>::cell_iterator &cell);

  /**
   * Do the work which is needed before cellwise data computation.
   * Since the shape functions are constructed independently on each cell,
   * the data on the reference cell is not necessary.
   * It returns an empty variable type of @ InternalDataBase and updates @ update_flags.
   */
  virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &,
            const Quadrature<2> &,
            dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const ;

  /**
   * Compute the data on the current cell.
   */
  virtual void
  fill_fe_values (const Triangulation<2,2>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                   ,
                  const Quadrature<2>                               &quadrature,
                  const Mapping<2,2>                                &mapping,
                  const Mapping<2,2>::InternalDataBase &,
                  const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                  const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                  internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;

  /**
   * Compute the data on the face of the current cell.
   */
  virtual void
  fill_fe_face_values (const Triangulation<2,2>::cell_iterator           &cell,
                       const unsigned int                                                   face_no,
                       const Quadrature<1>                                             &quadrature,
                       const Mapping<2,2>                                         &mapping,
                       const Mapping<2,2>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                       const InternalDataBase                                              &fe_internal,
                       dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;

  /**
   * Compute the data on the subface of the current cell.
   */
  virtual void
  fill_fe_subface_values (const Triangulation<2,2>::cell_iterator           &cell,
                          const unsigned int                                                   face_no,
                          const unsigned int                                                   sub_no,
                          const Quadrature<1>                                             &quadrature,
                          const Mapping<2,2>                                         &mapping,
                          const Mapping<2,2>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                          const InternalDataBase                                              &fe_internal,
                          dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;

  /**
   * Create the constraints matrix for hanging edges.
   */
  void initialize_constraints () ;

  class InternalData : public FiniteElement<2,2>::InternalDataBase
  {
  public:
    Table<2,Tensor<2,2> > shape_hessians;
  };

};




/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
