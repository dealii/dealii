// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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

#ifndef dealii_fe_dgp_nonparametric_h
#define dealii_fe_dgp_nonparametric_h

#include <deal.II/base/config.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup fe */
/*@{*/

/**
 * Discontinuous finite elements evaluated at the mapped quadrature points.
 *
 * Warning: this class does not work properly, yet. Don't use it!
 *
 * This finite element implements complete polynomial spaces, that is,
 * $d$-dimensional polynomials of order $k$.
 *
 * The polynomials are not mapped. Therefore, they are constant, linear,
 * quadratic, etc. on any grid cell.
 *
 * Since the polynomials are evaluated at the quadrature points of the actual
 * grid cell, no grid transfer and interpolation matrices are available.
 *
 * The purpose of this class is experimental, therefore the implementation
 * will remain incomplete.
 *
 * Besides, this class is not implemented for the codimension one case
 * (<tt>spacedim != dim</tt>).
 *
 *
 * <h3>Visualization of shape functions</h3> In 2d, the shape functions of
 * this element look as follows.
 *
 * <h4>$P_0$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_0$ element,
 * shape function 0 </td>
 *
 * <td align="center"></tr> </table>
 *
 * <h4>$P_1$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_1$ element, shape function 0 </td>
 *
 * <td align="center"> $P_1$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_1$ element,
 * shape function 2 </td>
 *
 * <td align="center"></td> </tr> </table>
 *
 *
 * <h4>$P_2$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 0 </td>
 *
 * <td align="center"> $P_2$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 2 </td>
 *
 * <td align="center"> $P_2$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 4 </td>
 *
 * <td align="center"> $P_2$ element, shape function 5 </td> </tr> </table>
 *
 *
 * <h4>$P_3$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 0 </td>
 *
 * <td align="center"> $P_3$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 2 </td>
 *
 * <td align="center"> $P_3$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 4 </td>
 *
 * <td align="center"> $P_3$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 6 </td>
 *
 * <td align="center"> $P_3$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 8 </td>
 *
 * <td align="center"> $P_3$ element, shape function 9 </td> </tr> </table>
 *
 *
 * <h4>$P_4$ element</h4> <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 0 </td>
 *
 * <td align="center"> $P_4$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 2 </td>
 *
 * <td align="center"> $P_4$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 4 </td>
 *
 * <td align="center"> $P_4$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 6 </td>
 *
 * <td align="center"> $P_4$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 8 </td>
 *
 * <td align="center"> $P_4$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0010.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0011.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 10 </td>
 *
 * <td align="center"> $P_4$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0012.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0013.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 12 </td>
 *
 * <td align="center"> $P_4$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0014.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_4$ element,
 * shape function 14 </td>
 *
 * <td align="center"></td> </tr> </table>
 *
 *
 * <h3> Implementation details </h3>
 *
 * This element does not have an InternalData class, unlike all other
 * elements, because the InternalData classes are used to store things that
 * can be computed once and reused multiple times (such as the values of shape
 * functions at quadrature points on the reference cell). However, because the
 * element is not mapped, this element has nothing that could be computed on
 * the reference cell -- everything needs to be computed on the real cell --
 * and consequently there is nothing we'd like to store in such an object. We
 * can thus simply use the members already provided by
 * FiniteElement::InternalDataBase without adding anything in a derived class
 * in this class.
 *
 * @author Guido Kanschat, 2002
 */
template <int dim, int spacedim=dim>
class FE_DGPNonparametric : public FiniteElement<dim,spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p k.
   */
  FE_DGPNonparametric (const unsigned int k);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGPNonparametric<dim>(degree)</tt>, with @p dim and @p
   * degree replaced by appropriate values.
   */
  virtual std::string get_name () const;

  virtual
  std::unique_ptr<FiniteElement<dim,spacedim> >
  clone() const;

  // for documentation, see the FiniteElement base class
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const;

  /**
   * This function is intended to return the value of a shape function at a
   * point on the reference cell. However, since the current element does not
   * implement shape functions by mapping from a reference cell, no shape
   * functions exist on the reference cell.
   *
   * Consequently, as discussed in the corresponding function in the base
   * class, FiniteElement::shape_value(), this function throws an exception of
   * type FiniteElement::ExcUnitShapeValuesDoNotExist.
   */
  virtual double shape_value (const unsigned int i,
                              const Point<dim> &p) const;

  /**
   * This function is intended to return the value of a shape function at a
   * point on the reference cell. However, since the current element does not
   * implement shape functions by mapping from a reference cell, no shape
   * functions exist on the reference cell.
   *
   * Consequently, as discussed in the corresponding function in the base
   * class, FiniteElement::shape_value_component(), this function throws an
   * exception of type FiniteElement::ExcUnitShapeValuesDoNotExist.
   */
  virtual double shape_value_component (const unsigned int i,
                                        const Point<dim> &p,
                                        const unsigned int component) const;

  /**
   * This function is intended to return the gradient of a shape function at a
   * point on the reference cell. However, since the current element does not
   * implement shape functions by mapping from a reference cell, no shape
   * functions exist on the reference cell.
   *
   * Consequently, as discussed in the corresponding function in the base
   * class, FiniteElement::shape_grad(), this function throws an exception of
   * type FiniteElement::ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<1,dim> shape_grad (const unsigned int  i,
                                    const Point<dim>   &p) const;

  /**
   * This function is intended to return the gradient of a shape function at a
   * point on the reference cell. However, since the current element does not
   * implement shape functions by mapping from a reference cell, no shape
   * functions exist on the reference cell.
   *
   * Consequently, as discussed in the corresponding function in the base
   * class, FiniteElement::shape_grad_component(), this function throws an
   * exception of type FiniteElement::ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
                                              const Point<dim> &p,
                                              const unsigned int component) const;

  /**
   * This function is intended to return the Hessian of a shape function at a
   * point on the reference cell. However, since the current element does not
   * implement shape functions by mapping from a reference cell, no shape
   * functions exist on the reference cell.
   *
   * Consequently, as discussed in the corresponding function in the base
   * class, FiniteElement::shape_grad_grad(), this function throws an
   * exception of type FiniteElement::ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
                                         const Point<dim> &p) const;

  /**
   * This function is intended to return the Hessian of a shape function at a
   * point on the reference cell. However, since the current element does not
   * implement shape functions by mapping from a reference cell, no shape
   * functions exist on the reference cell.
   *
   * Consequently, as discussed in the corresponding function in the base
   * class, FiniteElement::shape_grad_grad_component(), this function throws
   * an exception of type FiniteElement::ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
                                                   const Point<dim> &p,
                                                   const unsigned int component) const;

  /**
   * Return the polynomial degree of this finite element, i.e. the value
   * passed to the constructor.
   */
  unsigned int get_degree () const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * If, on a vertex, several finite elements are active, the hp code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of dofs_per_vertex
   * of the two finite elements. The first index of each pair denotes one of
   * the vertex dofs of the present element, whereas the second is the
   * corresponding index of the other finite element.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   *
   * For the FE_DGPNonparametric class the result is always true (independent
   * of the degree of the element), as it has no hanging nodes (being a
   * discontinuous element).
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementDomination::Domination
   * and in particular the
   * @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * @}
   */

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   *
   * This function is made virtual, since finite element objects are usually
   * accessed through pointers to their base class, rather than the class
   * itself.
   */
  virtual std::size_t memory_consumption () const;

protected:

  /**
   * Prepare internal data structures and fill in values independent of the
   * cell.
   */
  virtual
  std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>
  get_data (const UpdateFlags                                                    update_flags,
            const Mapping<dim,spacedim>                                         &mapping,
            const Quadrature<dim>                                               &quadrature,
            dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
  void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                                     cell_similarity,
                  const Quadrature<dim>                                               &quadrature,
                  const Mapping<dim,spacedim>                                         &mapping,
                  const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                  const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
  void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                       const unsigned int                                                   face_no,
                       const Quadrature<dim-1>                                             &quadrature,
                       const Mapping<dim,spacedim>                                         &mapping,
                       const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                       const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                       dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
  void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                          const unsigned int                                                   face_no,
                          const unsigned int                                                   sub_no,
                          const Quadrature<dim-1>                                             &quadrature,
                          const Mapping<dim,spacedim>                                         &mapping,
                          const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                          const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                          dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

private:

  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static
  std::vector<unsigned int>
  get_dpo_vector (const unsigned int degree);

  /**
   * Pointer to an object representing the polynomial space used here.
   */
  const PolynomialSpace<dim> polynomial_space;

  /**
   * Allow access from other dimensions.
   */
  template <int, int> friend class FE_DGPNonparametric;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
