// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_face_h
#define dealii_fe_face_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_poly_face.h>

DEAL_II_NAMESPACE_OPEN


/**
 * A finite element that is a tensor product polynomial on each face and
 * undefined in the interior of the cells. The basis functions on the faces
 * are Lagrange polynomials based on the support points of the
 * (dim-1)-dimensional Gauss--Lobatto quadrature rule. For element degree one
 * and two, the polynomials hence correspond to the usual Lagrange polynomials
 * on equidistant points.
 *
 * Although the name does not give it away, the element is discontinuous at
 * locations where faces of cells meet. In particular, this finite element is
 * the trace space of FE_RaviartThomas on the faces and serves in hybridized
 * methods, e.g. in combination with the FE_DGQ element. Its use is
 * demonstrated in the step-51 tutorial program.
 *
 * @note Since this element is defined only on faces, only FEFaceValues and
 * FESubfaceValues will provide useful information. On the other hand, if you
 * use this element with FEValues for cell integration, then the values
 * and derivatives of shape functions will have invalid values and will not
 * likely produce anything useful. In order to make the use of this element
 * as part of an FESystem simpler, using a (cell) FEValues object will not fail
 * outright, but those components of shape functions of the combined element
 * that correspond to FE_FaceQ will have the invalid values mentioned above.
 *
 * @ingroup fe
 */
template <int dim, int spacedim = dim>
class FE_FaceQ
  : public FE_PolyFace<TensorProductPolynomials<dim - 1>, dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Lagrange
   * polynomials in each coordinate direction.
   */
  FE_FaceQ(const unsigned int p);

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  /**
   * Implementation of the corresponding function in the FiniteElement
   * class.  Since the current element is interpolatory, the nodal
   * values are exactly the support point values. Furthermore, since
   * the current element is scalar, the support point values need to
   * be vectors of length 1.
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double>                 &matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double>                 &matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * @name Functions to support hp
   * @{
   */

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
   *
   * The set of such constraints is non-empty only for dim==1.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * The set of such constraints is non-empty only for dim==2.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * The set of such constraints is non-empty only for dim==3.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp-compatible".
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * @}
   */

  /**
   * Return a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * Specialization of FE_FaceQ for 1d. In that case, the finite element only
 * consists of one degree of freedom in each of the two faces (= vertices) of
 * a cell, irrespective of the degree. However, this element still accepts a
 * degree in its constructor and also returns that degree. This way,
 * dimension-independent programming with trace elements is also possible in
 * 1d (even though there is no computational benefit at all from it in 1d).
 *
 * @ingroup fe
 */
template <int spacedim>
class FE_FaceQ<1, spacedim> : public FiniteElement<1, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_FaceQ(const unsigned int p);

  virtual std::unique_ptr<FiniteElement<1, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<1, spacedim> &source,
                                FullMatrix<double>               &matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<1, spacedim> &source,
    const unsigned int                subface,
    FullMatrix<double>               &matrix,
    const unsigned int                face_no = 0) const override;

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
   */
  virtual bool
  hp_constraints_are_implemented() const override;

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
   *
   * The set of such constraints is non-empty only for dim==1.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<1, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * The set of such constraints is non-empty only for dim==2.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<1, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * The set of such constraints is non-empty only for dim==3.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<1, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * Return a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

protected:
  /*
   * NOTE: The following functions have their definitions inlined into the class
   * declaration because we otherwise run into a compiler error with MS Visual
   * Studio.
   */


  virtual std::unique_ptr<typename FiniteElement<1, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags /*update_flags*/,
    const Mapping<1, spacedim> & /*mapping*/,
    const Quadrature<1> & /*quadrature*/,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      & /*output_data*/) const override
  {
    return std::make_unique<
      typename FiniteElement<1, spacedim>::InternalDataBase>();
  }

  using FiniteElement<1, spacedim>::get_face_data;

  std::unique_ptr<typename FiniteElement<1, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags update_flags,
    const Mapping<1, spacedim> & /*mapping*/,
    const hp::QCollection<0> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      & /*output_data*/) const override
  {
    AssertDimension(quadrature.size(), 1);

    // generate a new data object and initialize some fields
    auto data_ptr =
      std::make_unique<typename FiniteElement<1, spacedim>::InternalDataBase>();
    data_ptr->update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature[0].size();
    AssertDimension(n_q_points, 1);
    (void)n_q_points;

    // No derivatives of this element are implemented.
    if (data_ptr->update_each & update_gradients ||
        data_ptr->update_each & update_hessians)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

    return data_ptr;
  }

  std::unique_ptr<typename FiniteElement<1, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags           update_flags,
    const Mapping<1, spacedim> &mapping,
    const Quadrature<0>        &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override
  {
    return get_face_data(update_flags,
                         mapping,
                         hp::QCollection<0>(quadrature),
                         output_data);
  }

  virtual void
  fill_fe_values(
    const typename Triangulation<1, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                          cell_similarity,
    const Quadrature<1>                                      &quadrature,
    const Mapping<1, spacedim>                               &mapping,
    const typename Mapping<1, spacedim>::InternalDataBase    &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<1, spacedim>
                                                                &mapping_data,
    const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<1, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<1, spacedim>::cell_iterator &cell,
    const unsigned int                                        face_no,
    const hp::QCollection<0>                                 &quadrature,
    const Mapping<1, spacedim>                               &mapping,
    const typename Mapping<1, spacedim>::InternalDataBase    &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<1, spacedim>
                                                                &mapping_data,
    const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<1, spacedim>::cell_iterator &cell,
    const unsigned int                                        face_no,
    const unsigned int                                        sub_no,
    const Quadrature<0>                                      &quadrature,
    const Mapping<1, spacedim>                               &mapping,
    const typename Mapping<1, spacedim>::InternalDataBase    &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<1, spacedim>
                                                                &mapping_data,
    const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * A finite element that is a Legendre element of complete polynomials on
 * each face (i.e., it is the face equivalent of what FE_DGP is on cells) and
 * undefined in the interior of the cells. The basis functions on the faces
 * are from Polynomials::Legendre.
 *
 * Although the name does not give it away, the element is discontinuous at
 * locations where faces of cells meet. The element serves in hybridized
 * methods, e.g. in combination with the FE_DGP element. An example of
 * hybridizes methods can be found in the step-51 tutorial program.
 *
 * @note Since this element is defined only on faces, only FEFaceValues and
 * FESubfaceValues will provide useful information. On the other hand, if you
 * use this element with FEValues for cell integration, then the values
 * and derivatives of shape functions will have invalid values and will not
 * likely produce anything useful. In order to make the use of this element
 * as part of an FESystem simpler, using a (cell) FEValues object will not fail
 * outright, but those components of shape functions of the combined element
 * that correspond to FE_FaceP will have the invalid values mentioned above.
 *
 * @ingroup fe
 */
template <int dim, int spacedim = dim>
class FE_FaceP : public FE_PolyFace<PolynomialSpace<dim - 1>, dim, spacedim>
{
public:
  /**
   * Constructor for complete basis of polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Legendre
   * polynomials in each coordinate direction.
   */
  FE_FaceP(unsigned int p);

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceP<dim>(degree)</tt> , with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double>                 &matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double>                 &matrix,
    const unsigned int                  face_no = 0) const override;

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
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * Return a list of constant modes of the element. For this element, the
   * first entry on each face is true, all other are false (as the constant
   * function is represented by the first base function of Legendre
   * polynomials).
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * FE_FaceP in 1d, i.e., with degrees of freedom on the element vertices.
 * See the documentation of the general template for more information.
 */
template <int spacedim>
class FE_FaceP<1, spacedim> : public FE_FaceQ<1, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_FaceP(const unsigned int p);

  /**
   * Return the name of the element
   */
  std::string
  get_name() const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif
