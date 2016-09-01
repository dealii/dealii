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
 *
 * Unlike any continuous conforming finite element belongs to H^1_0,
 * the P1 nonconforming element does not enforce the continuity across edges.
 * But it requires the continuity just in integral sense:
 * any function in the space should have the same integral values on two sides of the common edge shared by two adjacent elements.

 * Thus each function in the nonconforming element space can be discontinuous, not included in H^1_0, as functions in Discontinuous Galerkin (DG) finite element spaces.
 * Although any function in DG element space also has nonconformity, it is completely discontinuous across edges without any relation.
 * It is a reason why usual weak formulations for DG schemes contain additional penalty terms for jump across edges to control discontinuity.
 * However nonconforming elements usually do not need additional terms in their weak formulations due to the continuity in integral on edges.

 * <h3>DOFs and Dice Rule</h3>
 * Since any function in the P1 nonconforming space is piecewise linear on each element,
 * the function value at the mipoint of each edge is same to the mean value on the edge.
 * Thus the continuity of the integral value across each edge is equivalent to
 * the continuity of the midpoint value of each edge in this case.
 *
 * The degrees of freedom (DOFs) on a quadrilateral are defined by midpoint values on edges.
 * But these four DOFs are not independent in fact.
 * The simple observation reads that any linear function on a quadrilateral
 * satisfies 'dice rule': the sum of two function values at two midpoints of the edge pair on opposite
 * position is equal to the sum of those of the another edge pair.
 *
 * \phi(m_0) + \phi(m_1) = \phi(m_2) + \phi(m_3).
 *
 * Conversely if 4 values at midpoints satisfying the dice rule are just given,
 * then there always exists the unique linear function which coincides with 4 midpoints values.
 *
 * Due to the dice rule, three values at any three midpoints can determine the last value at the last midpoint.
 * It means that the genuine number of (independent) DOFs on a quad is 3,
 * and it is the same number to the dimension of the linear polynomial space in 2D.


 * <h3>Shape functions</h3>

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

 * For each vertex v_j of given quad, there are two edges of which v_j is one of end points.
 * Consider a linear function such that 0.5 value at two midpoints of such edges,
 * and 0.0 at two midpoints of other edges.
 * Note that the set of these values satisfies the dice rule which is described above.
 * We denote such a function assoicated with vertex v_j by \phi_j.

 * The canonical (local) basis functions are given as any three shape functions of
 * the following four linear functions:

 * <ul>
 * <li> shape function \phi_0:
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

 * <li> shape function \phi_1:
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

 * <li> shape function \phi_2:
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

 * <li> shape function \phi_3:
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

 * </ul>

 * Note that above shape functions are constructed on each cell, not on the reference cell only.
 * @p get_linear_shape computes the coefficients for shape functions when @p fill_fe_values is called on each cell.

 * The (global) basis function associated with a node is defined by the composition of
 * (local) basis functions associated with the node on each element.
 * When a problem with homogeneous Dirichlet boundary condition is dealt,
 * the total number of DOFs is equal to the number of interior nodes, as the standard bilinear finite element @p Q_1.

 * <h3>Unit support points</h3>
 * Contrast with ordinary Lagrange finite elements, DOF value with respect to the P1 nonconforming element at given node does not coincide with the function value at that node.
 * For instance, the (global) basis function associated with a node has 0.75 at that node, not 1.0.
 * Thus we need an interpolation operator which maps any smooth function into a function with proper DOF values in the P1 element space.
 * One natural interpolant associated with given smooth function is the linear function whose midpoint value at each edge is defined by
 * the average of two values at endpoints of the edge.
 * It provides appropriate weights used in @p unit_support_points.

 * <h3>References</h3>
 * You can find the paper about the P1NC element at
 * http://epubs.siam.org/doi/abs/10.1137/S0036142902404923.

 **/

class FE_P1NC : public FiniteElement<2,2>
{

public:
  /**
   * Constructor for nonparametric version of P1 nonconforming element.
   */
  FE_P1NC() ;

  /**
   * Return the name of the class for the element.
   */
  virtual std::string get_name () const ;

  /**
   * Return the update flags which are needed.
   */
  virtual UpdateFlags     requires_update_flags (const UpdateFlags flags) const ;

  /**
   * Copy constructor.
   */
  virtual FiniteElement<2,2> *clone () const ;

  /**
   * Destructor.
   */
  virtual ~FE_P1NC ();



private:

  static
  std::vector<ComponentMask>
  get_nonzero_component();


  /**
   * Return the vector consists of the numbers of degrees of freedom per objects.
   */
  static
  std::vector<unsigned int>
  get_dpo_vector ();


  /**
   * Compute the linear shape functions phi(x,y) = ax + by + c
   * such that each midpoint value on two connecting edges is a half,
   * and two other midpoint values are all zero.
   */
  static
  void
  get_linear_shape (const Triangulation<2,2>::cell_iterator &cell,
                    std::vector<double> &a,
                    std::vector<double> &b,
                    std::vector<double> &c);



  /**
   * Do the work which is needed before cellwise data computation.
   * Since the basis functions are constructed independently on each cell,
   * the data on the reference cell is not necessary.
   */
  virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &,
            const Quadrature<2> &,
            dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const ;




  /**
   * Compute the data on the current cell.
   */
  virtual
  void
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
  virtual
  void
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
  virtual
  void
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

};




/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
