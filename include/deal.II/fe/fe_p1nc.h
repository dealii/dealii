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
 * element, a piecewise linear finite element on quadrilaterals in 2D.
 *
 * Unlike any continuous conforming finite element,
 * it does not have the continuity across edges.
 * But it requires the continuity in weak sense:
 * a function in the space should have the same integral values on two sides of the common edge shared by two adjacent elements.
 *
 * Since any function in the space is piecewise linear on each element,
 * the continuity of the integral value across the edge is equivalent to
 * the continuity of the value at the midpoint of the edge.
 *
 * The degrees of freedom on a quadrilateral are given by midpoint values on edges.
 * However these four dofs in 2D are not independent in fact.
 * The simple observation reads that any linear function on a quadrilateral
 * satisfies 'the dice rule': the sum of two function values at two midpoints of the edge pair on opposite
 * position is equal to the sum of those on the another edge pair.
 *
 * \phi(m_0) + \phi(m_1) = \phi(m_2) + \phi(m_3).
 *
 * Due to the dice rule, three values at any three midpoints determine the last value at the last midpoint.
 * It means that the genuine number of independent dofs on a quad is 3,
 * and it is the same number to the dimension of the linear polynomial space in 2D.


 * Shape functions

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

 * For each vertex v_j, there are two edges of which v_j is one of the end points.
 * Consider the linear function such that one half at two midpoints of such edges,
 * and zero at two midpoints of other edges.
 * Note that this situation satisfies the dice rule which is described above.
 * We denote such a function by \phi_j.

 * The canonical (local) basis functions are given as any three shape functions of
 * the following four linear functions:

 * shape function \phi_0

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

 * shape function \phi_1

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

 * shape function \phi_2

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

 * shape function \phi_3

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


 * Note that this shape functions are constructed on each cell, not on the reference cell only.
 * get_linear_shape computes the coefficients for shape functions when fill_fe_values is called on each cell. 

 * The (global) basis function associated with a node is defined by the composition of
 * (local) basis functions associated with the node on each element.
 * In case of the problem with homogeneous Dirichlet boundary condition, 
 * the number of DOFs is equal to the number of interior nodes.

 * (TODO: unit_support_points)

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




/** @}*/

DEAL_II_NAMESPACE_CLOSE

#endif
